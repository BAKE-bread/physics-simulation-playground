import { clamp, SpatialHash, Random } from '../core/math.js';
import { Solver } from '../core/simulation.js';
import { LbmSolver } from './lbm.js';

const toWorldX = (i, n, width) => -width / 2 + (i / (n - 1)) * width;
const toWorldY = (j, n, height) => -height / 2 + (j / (n - 1)) * height;

/**
 * Lazily-created offscreen paint target. Fields are rasterised into an
 * ImageData at grid resolution and upsampled by drawImage with smoothing,
 * which removes the blocky per-cell look of the old fillRect renderer.
 * Guarded so headless numeric tests (no DOM) can still step the solvers.
 */
function ensurePainter(owner, n) {
  if (typeof document === 'undefined') return null;
  if (!owner.painter) {
    const canvas = document.createElement('canvas');
    owner.painter = { canvas, ctx: canvas.getContext('2d'), image: null, pixels: null, size: 0 };
  }
  const painter = owner.painter;
  if (painter.size !== n) {
    painter.size = n;
    canvasSize(painter.canvas, n, n);
    painter.image = painter.ctx.createImageData(n, n);
    painter.pixels = new Uint32Array(painter.image.data.buffer);
  }
  return painter;
}

function canvasSize(canvas, w, h) {
  canvas.width = w;
  canvas.height = h;
}

/** Pack RGBA into the little-endian Uint32 ImageData layout. */
const pack = (r, g, b, a) => ((a << 24) | (b << 16) | (g << 8) | r) >>> 0;

class GridSmoke {
  constructor(width, height, n = 110) {
    this.width = width;
    this.height = height;
    this.iterations = 34;
    this.time = 0;
    this.painter = null;
    this.setResolution(n);
  }

  setResolution(n) {
    if (this.n === n && this.u) return;
    this.n = n;
    this.dx = this.width / (n - 2);
    this.size = n * n;
    this.u = new Float32Array(this.size);
    this.v = new Float32Array(this.size);
    this.uOld = new Float32Array(this.size);
    this.vOld = new Float32Array(this.size);
    this.dye = new Float32Array(this.size);
    this.dyeOld = new Float32Array(this.size);
    this.pressure = new Float32Array(this.size);
    this.pressureNext = new Float32Array(this.size);
    this.divergence = new Float32Array(this.size);
    this.sources = [];
    this.painter = null;
  }

  index(i, j) { return i + j * this.n; }

  clear() {
    for (const field of [this.u, this.v, this.uOld, this.vOld, this.dye, this.dyeOld, this.pressure, this.pressureNext, this.divergence]) field.fill(0);
    this.sources = [];
    this.time = 0;
  }

  clampGrid(value) { return clamp(value, 0, this.n - 1.001); }

  sample(field, x, y) {
    const gx = this.clampGrid(x);
    const gy = this.clampGrid(y);
    const x0 = Math.floor(gx);
    const y0 = Math.floor(gy);
    const x1 = Math.min(x0 + 1, this.n - 1);
    const y1 = Math.min(y0 + 1, this.n - 1);
    const tx = gx - x0;
    const ty = gy - y0;
    const a = field[this.index(x0, y0)] * (1 - tx) + field[this.index(x1, y0)] * tx;
    const b = field[this.index(x0, y1)] * (1 - tx) + field[this.index(x1, y1)] * tx;
    return a * (1 - ty) + b * ty;
  }

  enforceBoundary(field) {
    const n = this.n;
    for (let i = 0; i < n; i += 1) {
      field[this.index(i, 0)] = 0;
      field[this.index(i, n - 1)] = 0;
      field[this.index(0, i)] = 0;
      field[this.index(n - 1, i)] = 0;
    }
  }

  inject(x, y, amount = 1, vx = 0, vy = 0) {
    this.sources.push({ x, y, amount, vx, vy });
  }

  /** Source radius scales with the grid so the physical plume size is stable. */
  sourceRadius() {
    return Math.max(2.2, this.n * 0.062);
  }

  applySources() {
    const radius = this.sourceRadius();
    const falloff = 2.4 / (radius * radius * 0.32);
    for (const source of this.sources) {
      const cx = ((source.x + this.width / 2) / this.width) * (this.n - 1);
      const cy = ((source.y + this.height / 2) / this.height) * (this.n - 1);
      for (let j = Math.max(1, Math.floor(cy - radius)); j <= Math.min(this.n - 2, Math.ceil(cy + radius)); j += 1) {
        for (let i = Math.max(1, Math.floor(cx - radius)); i <= Math.min(this.n - 2, Math.ceil(cx + radius)); i += 1) {
          const d2 = (i - cx) ** 2 + (j - cy) ** 2;
          if (d2 > radius * radius) continue;
          const weight = Math.exp(-d2 * falloff);
          const index = this.index(i, j);
          this.dye[index] += source.amount * weight;
          this.u[index] += source.vx * weight;
          this.v[index] += source.vy * weight;
        }
      }
    }
    this.sources = [];
  }

  advect(destination, source, velocityU, velocityV, dt, dissipation = 1) {
    const scale = dt / this.dx;
    for (let j = 1; j < this.n - 1; j += 1) {
      for (let i = 1; i < this.n - 1; i += 1) {
        const index = this.index(i, j);
        destination[index] = this.sample(source, i - scale * velocityU[index], j - scale * velocityV[index]) * dissipation;
      }
    }
    this.enforceBoundary(destination);
  }

  project() {
    const dx = this.dx;
    for (let j = 1; j < this.n - 1; j += 1) {
      for (let i = 1; i < this.n - 1; i += 1) {
        const index = this.index(i, j);
        this.divergence[index] = (this.u[this.index(i + 1, j)] - this.u[this.index(i - 1, j)]
          + this.v[this.index(i, j + 1)] - this.v[this.index(i, j - 1)]) / (2 * dx);
      }
    }
    this.pressure.fill(0);
    for (let iteration = 0; iteration < this.iterations; iteration += 1) {
      for (let j = 1; j < this.n - 1; j += 1) {
        for (let i = 1; i < this.n - 1; i += 1) {
          const index = this.index(i, j);
          this.pressureNext[index] = (this.pressure[this.index(i - 1, j)] + this.pressure[this.index(i + 1, j)]
            + this.pressure[this.index(i, j - 1)] + this.pressure[this.index(i, j + 1)]
            - this.divergence[index] * dx * dx) * 0.25;
        }
      }
      [this.pressure, this.pressureNext] = [this.pressureNext, this.pressure];
      this.enforceBoundary(this.pressure);
    }
    for (let j = 1; j < this.n - 1; j += 1) {
      for (let i = 1; i < this.n - 1; i += 1) {
        const index = this.index(i, j);
        this.u[index] -= (this.pressure[this.index(i + 1, j)] - this.pressure[this.index(i - 1, j)]) / (2 * dx);
        this.v[index] -= (this.pressure[this.index(i, j + 1)] - this.pressure[this.index(i, j - 1)]) / (2 * dx);
      }
    }
    this.enforceBoundary(this.u);
    this.enforceBoundary(this.v);
  }

  step(dt) {
    this.time += dt;
    // A gentle continuous source makes the incompressible projection visible before interaction.
    const y = -1.2 + Math.sin(this.time * 1.25) * 0.55;
    this.inject(-4.95, y, 0.19, 2.15, 0.24 * Math.cos(this.time * 1.25));
    this.applySources();
    this.uOld.set(this.u); this.vOld.set(this.v);
    this.advect(this.u, this.uOld, this.uOld, this.vOld, dt, 0.999);
    this.advect(this.v, this.vOld, this.uOld, this.vOld, dt, 0.999);
    this.project();
    this.dyeOld.set(this.dye);
    this.advect(this.dye, this.dyeOld, this.u, this.v, dt, 0.996);
  }

  divergenceNorm() {
    let total = 0;
    let count = 0;
    for (let j = 1; j < this.n - 1; j += 1) {
      for (let i = 1; i < this.n - 1; i += 1) {
        const div = (this.u[this.index(i + 1, j)] - this.u[this.index(i - 1, j)]
          + this.v[this.index(i, j + 1)] - this.v[this.index(i, j - 1)]) / (2 * this.dx);
        total += div * div;
        count += 1;
      }
    }
    return Math.sqrt(total / Math.max(count, 1));
  }

  render(ctx) {
    const painter = ensurePainter(this, this.n);
    if (painter) {
      const n = this.n;
      const pixels = painter.pixels;
      const dye = this.dye;
      for (let index = 0; index < this.size; index += 1) {
        const d = dye[index];
        if (d < 0.012) { pixels[index] = 0; continue; }
        const t = clamp(d / 1.25, 0, 1);
        const glow = t * t;
        const r = Math.round(14 + 96 * t + 66 * glow);
        const g = Math.round(52 + 158 * t + 38 * glow);
        const b = Math.round(84 + 150 * t + 20 * glow);
        const a = Math.round(clamp(d * 1.5, 0, 1) * 235 + 20);
        pixels[index] = pack(Math.min(r, 255), Math.min(g, 255), Math.min(b, 255), a);
      }
      painter.ctx.putImageData(painter.image, 0, 0);
      ctx.save();
      ctx.imageSmoothingEnabled = true;
      ctx.drawImage(painter.canvas, -this.width / 2, -this.height / 2, this.width, this.height);
      ctx.restore();
    }
    ctx.strokeStyle = 'rgba(180, 241, 255, 0.3)';
    ctx.lineWidth = 0.012;
    const stride = Math.max(3, Math.round(this.n / 22));
    for (let j = 2; j < this.n - 2; j += stride) {
      for (let i = 2; i < this.n - 2; i += stride) {
        const index = this.index(i, j);
        const x = toWorldX(i, this.n, this.width);
        const y = toWorldY(j, this.n, this.height);
        ctx.beginPath();
        ctx.moveTo(x, y);
        ctx.lineTo(x + this.u[index] * 0.18, y + this.v[index] * 0.18);
        ctx.stroke();
      }
    }
  }

  diagnostics() {
    let mass = 0;
    for (const density of this.dye) mass += density;
    return {
      entities: this.n * this.n,
      contacts: this.iterations,
      energy: mass * this.dx * this.dx,
      residual: this.divergenceNorm(),
      method: `半拉格朗日平流 + 压力投影 · ${this.n}² 网格`,
    };
  }
}

class ParticleWater {
  constructor(width, height) {
    this.width = width;
    this.height = height;
    this.radius = 0.24;
    this.hash = new SpatialHash(this.radius);
    this.particles = [];
    this.gravity = -9.81;
    // These are calibrated in the unnormalised 2D compact-kernel units below.
    this.restDensity = 0.1;
    this.stiffness = 8000;
    this.nearStiffness = 30000;
    this.viscosity = 0.075;
    this.substeps = 2;
    this.lastPointer = null;
    this.random = new Random(0xfeed);
    this.reset();
  }

  reset() {
    this.particles = [];
    const columns = 25;
    const rows = 14;
    const separation = 0.205;
    // The block origin wanders a little between resets/shuffles so repeated
    // runs are not pixel-identical.
    const originX = -4.9 + this.random.range(-0.35, 1.35);
    const originY = -2.9 + this.random.range(0, 0.55);
    for (let y = 0; y < rows; y += 1) {
      for (let x = 0; x < columns; x += 1) {
        this.particles.push({
          x: originX + x * separation, y: originY + y * separation,
          px: 0, py: 0, vx: this.random.range(-0.35, 0.35), vy: this.random.range(-0.1, 0.25),
          density: 0, nearDensity: 0, hitFloor: false,
        });
      }
    }
    this.lastPointer = null;
  }

  /** Shuffle button: reseed, then rebuild the block at a fresh position. */
  randomize() {
    this.random = new Random((Math.random() * 2 ** 32) >>> 0);
    this.reset();
  }

  buildHash() {
    this.hash.clear();
    this.particles.forEach((p, index) => this.hash.insert(index, p.x, p.y));
  }

  bounds(p) {
    const left = -this.width / 2 + this.radius * 0.55;
    const right = this.width / 2 - this.radius * 0.55;
    const floor = -this.height / 2 + this.radius * 0.55;
    p.hitFloor = false;
    if (p.x < left) p.x = left;
    if (p.x > right) p.x = right;
    if (p.y < floor) { p.y = floor; p.hitFloor = true; }
  }

  step(dt) {
    const h = dt / this.substeps;
    for (let substep = 0; substep < this.substeps; substep += 1) {
      for (const p of this.particles) {
        p.px = p.x; p.py = p.y;
        p.vy += this.gravity * h;
        p.x += p.vx * h;
        p.y += p.vy * h;
      }
      this.buildHash();
      for (const p of this.particles) { p.density = 0; p.nearDensity = 0; }
      for (let i = 0; i < this.particles.length; i += 1) {
        const a = this.particles[i];
        this.hash.forEachNearby(a.x, a.y, this.radius, (j) => {
          if (j <= i) return;
          const b = this.particles[j];
          const distance = Math.hypot(b.x - a.x, b.y - a.y);
          if (distance >= this.radius) return;
          const q = 1 - distance / this.radius;
          const q2 = q * q;
          a.density += q2; b.density += q2;
          a.nearDensity += q2 * q; b.nearDensity += q2 * q;
        });
      }
      for (let i = 0; i < this.particles.length; i += 1) {
        const a = this.particles[i];
        const pressureA = this.stiffness * Math.max(0, a.density - this.restDensity);
        const nearPressureA = this.nearStiffness * a.nearDensity;
        this.hash.forEachNearby(a.x, a.y, this.radius, (j) => {
          if (j <= i) return;
          const b = this.particles[j];
          const dx = b.x - a.x;
          const dy = b.y - a.y;
          const distance = Math.hypot(dx, dy);
          if (distance >= this.radius || distance < 1e-7) return;
          const q = 1 - distance / this.radius;
          const pressureB = this.stiffness * Math.max(0, b.density - this.restDensity);
          const nearPressureB = this.nearStiffness * b.nearDensity;
          const displacement = (h * h) * ((pressureA + pressureB) * q + (nearPressureA + nearPressureB) * q * q) * 0.5;
          const nx = dx / distance;
          const ny = dy / distance;
          a.x -= nx * displacement * 0.5; a.y -= ny * displacement * 0.5;
          b.x += nx * displacement * 0.5; b.y += ny * displacement * 0.5;
        });
      }
      for (const p of this.particles) this.bounds(p);
      for (const p of this.particles) {
        p.vx = (p.x - p.px) / h;
        p.vy = (p.y - p.py) / h;
        if (p.hitFloor && p.vy < 0) { p.vy *= -0.17; p.vx *= 0.84; }
      }
      this.buildHash();
      for (let i = 0; i < this.particles.length; i += 1) {
        const a = this.particles[i];
        this.hash.forEachNearby(a.x, a.y, this.radius, (j) => {
          if (j <= i) return;
          const b = this.particles[j];
          const distance = Math.hypot(b.x - a.x, b.y - a.y);
          if (distance >= this.radius) return;
          const q = 1 - distance / this.radius;
          const dvx = (b.vx - a.vx) * this.viscosity * q * h;
          const dvy = (b.vy - a.vy) * this.viscosity * q * h;
          a.vx += dvx; a.vy += dvy; b.vx -= dvx; b.vy -= dvy;
        });
      }
    }
  }

  pointer(kind, x, y) {
    if (kind === 'down') this.lastPointer = { x, y };
    if (kind === 'move' && this.lastPointer) {
      const vx = clamp((x - this.lastPointer.x) * 28, -12, 12);
      const vy = clamp((y - this.lastPointer.y) * 28, -12, 12);
      for (const p of this.particles) {
        const dx = p.x - x;
        const dy = p.y - y;
        const d2 = dx * dx + dy * dy;
        if (d2 < 0.75 ** 2) {
          const weight = 1 - Math.sqrt(d2) / 0.75;
          p.vx += vx * weight * 0.16;
          p.vy += vy * weight * 0.16;
        }
      }
      this.lastPointer = { x, y };
    }
    if (kind === 'up' || kind === 'leave') this.lastPointer = null;
  }

  render(ctx) {
    for (const p of this.particles) {
      const color = clamp(p.density / (this.restDensity * 1.35), 0, 1);
      ctx.fillStyle = `rgba(${Math.round(57 + 36 * color)}, ${Math.round(165 + 75 * color)}, 242, 0.78)`;
      ctx.beginPath(); ctx.arc(p.x, p.y, this.radius * 0.47, 0, Math.PI * 2); ctx.fill();
    }
  }

  diagnostics() {
    let density = 0;
    let energy = 0;
    for (const p of this.particles) { density += p.density; energy += 0.5 * (p.vx * p.vx + p.vy * p.vy); }
    return { entities: this.particles.length, contacts: this.hash.cells.size, energy, residual: Math.abs(density / this.particles.length - this.restDensity) / this.restDensity, method: 'SPH 核密度 + 双密度松弛' };
  }
}

class ShallowWater {
  constructor(width, height, n = 124) {
    this.width = width;
    this.height = height;
    this.depth = 0.58;
    this.gravity = 9.81;
    this.painter = null;
    this.setResolution(n);
  }

  setResolution(n) {
    if (this.n === n && this.eta) return;
    this.n = n;
    this.dx = this.width / (n - 1);
    this.size = n * n;
    this.eta = new Float32Array(this.size);
    this.etaNext = new Float32Array(this.size);
    this.u = new Float32Array(this.size);
    this.v = new Float32Array(this.size);
    this.uNext = new Float32Array(this.size);
    this.vNext = new Float32Array(this.size);
    this.pending = [];
    this.painter = null;
  }

  index(i, j) { return i + j * this.n; }

  clear() { this.eta.fill(0); this.etaNext.fill(0); this.u.fill(0); this.v.fill(0); this.pending = []; }

  disturb(x, y, amplitude, vx = 0, vy = 0) { this.pending.push({ x, y, amplitude, vx, vy }); }

  applyDisturbances() {
    const radius = Math.max(3, this.n * 0.055);
    const falloff = 2.2 / (radius * radius * 0.45);
    for (const source of this.pending) {
      const gx = ((source.x + this.width / 2) / this.width) * (this.n - 1);
      const gy = ((source.y + this.height / 2) / this.height) * (this.n - 1);
      for (let j = Math.max(1, Math.floor(gy - radius)); j <= Math.min(this.n - 2, Math.ceil(gy + radius)); j += 1) {
        for (let i = Math.max(1, Math.floor(gx - radius)); i <= Math.min(this.n - 2, Math.ceil(gx + radius)); i += 1) {
          const weight = Math.exp(-((i - gx) ** 2 + (j - gy) ** 2) * falloff);
          const index = this.index(i, j);
          this.eta[index] += source.amplitude * weight;
          this.u[index] += source.vx * weight;
          this.v[index] += source.vy * weight;
        }
      }
    }
    this.pending = [];
  }

  step(dt) {
    this.applyDisturbances();
    const waveCfl = Math.sqrt(this.gravity * this.depth) * dt / this.dx;
    const substeps = Math.max(1, Math.ceil(waveCfl / 0.42));
    const h = dt / substeps;
    for (let sub = 0; sub < substeps; sub += 1) {
      for (let j = 1; j < this.n - 1; j += 1) {
        for (let i = 1; i < this.n - 1; i += 1) {
          const index = this.index(i, j);
          this.uNext[index] = (this.u[index] - h * this.gravity * (this.eta[this.index(i + 1, j)] - this.eta[this.index(i - 1, j)]) / (2 * this.dx)) * 0.998;
          this.vNext[index] = (this.v[index] - h * this.gravity * (this.eta[this.index(i, j + 1)] - this.eta[this.index(i, j - 1)]) / (2 * this.dx)) * 0.998;
        }
      }
      for (let j = 1; j < this.n - 1; j += 1) {
        for (let i = 1; i < this.n - 1; i += 1) {
          const index = this.index(i, j);
          const divergence = (this.uNext[this.index(i + 1, j)] - this.uNext[this.index(i - 1, j)]
            + this.vNext[this.index(i, j + 1)] - this.vNext[this.index(i, j - 1)]) / (2 * this.dx);
          this.etaNext[index] = (this.eta[index] - h * this.depth * divergence) * 0.9995;
        }
      }
      for (let i = 0; i < this.n; i += 1) {
        const low = this.index(i, 0); const high = this.index(i, this.n - 1);
        this.etaNext[low] = 0; this.etaNext[high] = 0; this.uNext[low] = 0; this.uNext[high] = 0; this.vNext[low] = 0; this.vNext[high] = 0;
        const left = this.index(0, i); const right = this.index(this.n - 1, i);
        this.etaNext[left] = 0; this.etaNext[right] = 0; this.uNext[left] = 0; this.uNext[right] = 0; this.vNext[left] = 0; this.vNext[right] = 0;
      }
      [this.eta, this.etaNext] = [this.etaNext, this.eta];
      [this.u, this.uNext] = [this.uNext, this.u];
      [this.v, this.vNext] = [this.vNext, this.v];
    }
  }

  pointer(kind, x, y) {
    if (kind === 'down') this.last = { x, y };
    if (kind === 'move' && this.last) {
      this.disturb(x, y, 0.035, (x - this.last.x) * 4, (y - this.last.y) * 4);
      this.last = { x, y };
    }
    if (kind === 'up' || kind === 'leave') this.last = null;
  }

  render(ctx) {
    const painter = ensurePainter(this, this.n);
    if (painter) {
      const n = this.n;
      const pixels = painter.pixels;
      for (let j = 0; j < n; j += 1) {
        for (let i = 0; i < n; i += 1) {
          const index = this.index(i, j);
          const eta = this.eta[index];
          const accent = clamp((eta + 0.16) / 0.32, 0, 1);
          // Specular crest highlight from the local surface slope.
          const grad = Math.abs(this.eta[this.index(Math.min(i + 1, n - 1), j)] - this.eta[this.index(Math.max(i - 1, 0), j)]);
          const crest = clamp(grad * 26, 0, 1) * accent;
          const r = Math.round(14 + accent * 34 + crest * 120);
          const g = Math.round(62 + accent * 108 + crest * 90);
          const b = Math.round(102 + accent * 128 + crest * 30);
          pixels[index] = pack(Math.min(r, 255), Math.min(g, 255), Math.min(b, 255), 255);
        }
      }
      painter.ctx.putImageData(painter.image, 0, 0);
      ctx.save();
      ctx.imageSmoothingEnabled = true;
      ctx.drawImage(painter.canvas, -this.width / 2, -this.height / 2, this.width, this.height);
      ctx.restore();
    }
    ctx.strokeStyle = 'rgba(185, 246, 255, 0.24)'; ctx.lineWidth = 0.012;
    const stride = Math.max(4, Math.round(this.n / 15));
    for (let j = 4; j < this.n - 4; j += stride) {
      for (let i = 4; i < this.n - 4; i += stride) {
        const index = this.index(i, j); const x = toWorldX(i, this.n, this.width); const y = toWorldY(j, this.n, this.height);
        ctx.beginPath(); ctx.moveTo(x, y); ctx.lineTo(x + this.u[index] * 0.12, y + this.v[index] * 0.12); ctx.stroke();
      }
    }
  }

  diagnostics() {
    let mass = 0; let energy = 0;
    for (let i = 0; i < this.size; i += 1) {
      mass += this.eta[i];
      energy += 0.5 * this.gravity * this.eta[i] ** 2 + 0.5 * this.depth * (this.u[i] ** 2 + this.v[i] ** 2);
    }
    return {
      entities: this.size,
      contacts: Math.ceil(Math.sqrt(this.gravity * this.depth) / this.dx),
      energy,
      residual: Math.abs(mass) / this.size,
      method: `线性浅水方程 · 有限差分 · ${this.n}² 网格`,
    };
  }
}

/** A dispatcher deliberately keeps each discretisation independent and swappable. */
export class FluidSolver extends Solver {
  constructor({ width = 12, height = 7 } = {}) {
    super();
    this.width = width;
    this.height = height;
    this.mode = 'smoke';
    this.smoke = new GridSmoke(width, height);
    this.water = new ParticleWater(width, height);
    this.shallow = new ShallowWater(width, height);
    this.lbm = new LbmSolver(width, height);
  }

  current() {
    if (this.mode === 'sph') return this.water;
    if (this.mode === 'shallow') return this.shallow;
    if (this.mode === 'lbm') return this.lbm;
    return this.smoke;
  }

  reset(mode = this.mode) { this.mode = mode; this.current().clear?.(); this.current().reset?.(); }
  step(dt) { this.current().step(dt); }
  render(ctx) { this.current().render(ctx); }
  pointer(kind, x, y) {
    if (this.mode === 'smoke') {
      if (kind === 'down') this.last = { x, y };
      if (kind === 'move' && this.last) {
        this.smoke.inject(x, y, 0.42, (x - this.last.x) * 20, (y - this.last.y) * 20);
        this.last = { x, y };
      }
      if (kind === 'up' || kind === 'leave') this.last = null;
    } else this.current().pointer(kind, x, y);
  }
  diagnostics() { return this.current().diagnostics(); }
}

export { GridSmoke, ParticleWater, ShallowWater };
