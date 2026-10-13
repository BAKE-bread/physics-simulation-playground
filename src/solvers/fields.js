import { clamp, Random } from '../core/math.js';
import { Solver } from '../core/simulation.js';
import { FlockSolver } from './flock.js';

/** Colour packer for the little-endian Uint32 ImageData layout. */
const pack = (r, g, b, a) => ((a << 24) | (b << 16) | (g << 8) | r) >>> 0;

function ensurePainter(owner, n) {
  if (typeof document === 'undefined') return null;
  if (!owner.painter) {
    const canvas = document.createElement('canvas');
    owner.painter = { canvas, ctx: canvas.getContext('2d'), image: null, pixels: null, size: 0 };
  }
  const painter = owner.painter;
  if (painter.size !== n) {
    painter.size = n;
    painter.canvas.width = n;
    painter.canvas.height = n;
    painter.image = painter.ctx.createImageData(n, n);
    painter.pixels = new Uint32Array(painter.image.data.buffer);
  }
  return painter;
}

/**
 * Transient heat conduction on an unstructured triangular mesh of a plate
 * with a circular hole — the textbook FEM showcase. Linear P1 elements,
 * lumped mass, backward Euler time integration solved with Jacobi sweeps.
 * Outer boundary is held at T = 0, the hole is adiabatic (natural BC).
 */
class HeatPlateFem {
  constructor(width, height) {
    this.width = width;
    this.height = height;
    this.plateW = 9.2;
    this.plateH = 4.6;
    this.hole = { x: 1.1, y: 0.15, r: 0.8 };
    this.heater = { x: -3.35, y: 1.25, r: 0.42 };
    this.kappa = 0.42;
    this.level = 0;
    this.setQuality(2);
  }

  setQuality(level) {
    if (level === this.level) return;
    this.level = level;
    const h = { 1: 0.46, 2: 0.34, 3: 0.26 }[level] ?? 0.34;
    this.meshH = h;
    this.buildMesh(h);
    this.reset();
  }

  buildMesh(h) {
    // Preserve the current temperature field across re-meshing (hole drags).
    const oldNodes = this.nodes;
    const oldT = this.temperature;
    const oldGrid = this.gridIndex ? { index: this.gridIndex, nx: this.gridNx, ny: this.gridNy } : null;

    const nx = Math.round(this.plateW / h) + 1;
    const ny = Math.round(this.plateH / h) + 1;
    const xs = [];
    const ys = [];
    const gridIndex = new Int32Array(nx * ny).fill(-1);
    const nodes = [];
    for (let j = 0; j < ny; j += 1) {
      for (let i = 0; i < nx; i += 1) {
        let x = -this.plateW / 2 + i * (this.plateW / (nx - 1));
        let y = -this.plateH / 2 + j * (this.plateH / (ny - 1));
        const d = Math.hypot(x - this.hole.x, y - this.hole.y);
        if (d < this.hole.r - h * 0.3) continue;
        // Snap the annulus of near-hole nodes onto the exact circle so the
        // adiabatic boundary is a smooth polygon instead of jagged steps.
        if (d < this.hole.r + h * 1.15 && d > 1e-9) {
          x = this.hole.x + ((x - this.hole.x) / d) * this.hole.r;
          y = this.hole.y + ((y - this.hole.y) / d) * this.hole.r;
        }
        gridIndex[i + j * nx] = nodes.length;
        nodes.push({
          x, y, gi: i, gj: j,
          dirichlet: i === 0 || j === 0 || i === nx - 1 || j === ny - 1,
        });
        xs.push(x); ys.push(y);
      }
    }
    const elements = [];
    for (let j = 0; j < ny - 1; j += 1) {
      for (let i = 0; i < nx - 1; i += 1) {
        const n00 = gridIndex[i + j * nx];
        const n10 = gridIndex[i + 1 + j * nx];
        const n01 = gridIndex[i + (j + 1) * nx];
        const n11 = gridIndex[i + 1 + (j + 1) * nx];
        if (n00 >= 0 && n10 >= 0 && n11 >= 0) elements.push([n00, n10, n11]);
        if (n00 >= 0 && n11 >= 0 && n01 >= 0) elements.push([n00, n11, n01]);
      }
    }
    // P1 stiffness assembly into per-node adjacency lists.
    const n = nodes.length;
    const adj = Array.from({ length: n }, () => []);
    const mass = new Float64Array(n);
    for (const [a, b, c] of elements) {
      const ax = xs[a]; const ay = ys[a];
      const bx = xs[b]; const by = ys[b];
      const cx = xs[c]; const cy = ys[c];
      const area2 = Math.abs((bx - ax) * (cy - ay) - (cx - ax) * (by - ay));
      if (area2 < 1e-12) continue;
      const area = area2 / 2;
      const grads = [
        [(by - cy) / area2, (cx - bx) / area2],
        [(cy - ay) / area2, (ax - cx) / area2],
        [(ay - by) / area2, (bx - ax) / area2],
      ];
      const tri = [a, b, c];
      for (let p = 0; p < 3; p += 1) {
        mass[tri[p]] += area / 3;
        for (let q = 0; q < 3; q += 1) {
          const k = this.kappa * area * (grads[p][0] * grads[q][0] + grads[p][1] * grads[q][1]);
          adj[tri[p]].push([tri[q], k]);
        }
      }
    }
    // Merge duplicate neighbour entries into flat CSR-ish arrays.
    this.nodeCount = n;
    this.nodes = nodes;
    this.elements = elements;
    this.mass = mass;
    this.nbrIndex = new Array(n);
    this.nbrValue = new Array(n);
    for (let i = 0; i < n; i += 1) {
      const merged = new Map();
      for (const [j, k] of adj[i]) merged.set(j, (merged.get(j) ?? 0) + k);
      this.nbrIndex[i] = Int32Array.from(merged.keys());
      this.nbrValue[i] = Float64Array.from(merged.values());
    }
    this.temperature = new Float64Array(n);
    this.tempNext = new Float64Array(n);
    this.mOverDt = null;
    this.mOverDtDt = 0;
    this.gridIndex = gridIndex;
    this.gridNx = nx;
    this.gridNy = ny;
    // Carry the old field onto the new mesh by matching grid coordinates.
    if (oldNodes && oldT && oldGrid) {
      for (let i = 0; i < n; i += 1) {
        const node = nodes[i];
        const gi = clamp(node.gi, 0, oldGrid.nx - 1);
        const gj = clamp(node.gj, 0, oldGrid.ny - 1);
        const oldIndex = oldGrid.index[gi + gj * oldGrid.nx];
        if (oldIndex >= 0) this.temperature[i] = oldT[oldIndex];
      }
    }
  }

  moveHole(x, y) {
    const margin = this.hole.r + 0.35;
    this.hole.x = clamp(x, -this.plateW / 2 + margin, this.plateW / 2 - margin);
    this.hole.y = clamp(y, -this.plateH / 2 + margin, this.plateH / 2 - margin);
    this.buildMesh(this.meshH);
  }

  reset() {
    this.temperature.fill(0);
    for (let i = 0; i < this.nodeCount; i += 1) {
      const node = this.nodes[i];
      if (node.dirichlet) continue;
      const d2 = (node.x + 2.6) ** 2 + (node.y + 0.7) ** 2;
      this.temperature[i] = 2.1 * Math.exp(-d2 * 1.4);
    }
    this.lastResidual = 0;
  }

  heatAt(x, y, amount, radius = 0.55) {
    for (let i = 0; i < this.nodeCount; i += 1) {
      const node = this.nodes[i];
      if (node.dirichlet) continue;
      const d = Math.hypot(node.x - x, node.y - y);
      if (d < radius) {
        const weight = 1 - d / radius;
        this.temperature[i] = Math.min(this.temperature[i] + amount * weight, 4.5);
      }
    }
  }

  step(dt) {
    // A small steady heater keeps the plate alive without interaction.
    this.heatAt(this.heater.x, this.heater.y, 14 * dt, this.heater.r);
    const { temperature: T, tempNext: Tn, mass, nodeCount, nodes } = this;
    if (!this.mOverDt || this.mOverDtDt !== dt) {
      this.mOverDt = new Float64Array(nodeCount);
      for (let i = 0; i < nodeCount; i += 1) this.mOverDt[i] = mass[i] / dt;
      this.mOverDtDt = dt;
    }
    const mOverDt = this.mOverDt;
    Tn.set(T);
    let residual = 0;
    for (let sweep = 0; sweep < 26; sweep += 1) {
      residual = 0;
      for (let i = 0; i < nodeCount; i += 1) {
        if (nodes[i].dirichlet) { Tn[i] = 0; continue; }
        const idx = this.nbrIndex[i];
        const val = this.nbrValue[i];
        let sum = mOverDt[i] * T[i];
        let diag = mOverDt[i];
        for (let k = 0; k < idx.length; k += 1) {
          const j = idx[k];
          if (j === i) diag += val[k];
          else sum -= val[k] * Tn[j];
        }
        const next = sum / diag;
        residual = Math.max(residual, Math.abs(next - Tn[i]));
        Tn[i] = next;
      }
    }
    [this.temperature, this.tempNext] = [Tn, T];
    this.lastResidual = residual;
  }

  pointer(kind, x, y) {
    if (kind === 'down') {
      this.dragHole = Math.hypot(x - this.hole.x, y - this.hole.y) < this.hole.r * 1.25;
      if (this.dragHole) this.moveHole(x, y);
      else this.heatAt(x, y, 1.6);
    } else if (kind === 'move') {
      if (this.dragHole) this.moveHole(x, y);
      else this.heatAt(x, y, 0.55);
    } else if (kind === 'up' || kind === 'leave') {
      this.dragHole = false;
    }
  }

  render(ctx) {
    const xs = this.nodes;
    for (const [a, b, c] of this.elements) {
      const t = clamp((this.temperature[a] + this.temperature[b] + this.temperature[c]) / 4.8, 0, 1);
      const hot = t * t;
      ctx.fillStyle = `rgb(${Math.round(16 + 239 * hot + 30 * t)}, ${Math.round(40 + 110 * t + 60 * hot)}, ${Math.round(62 + 30 * t - 40 * hot)})`;
      ctx.beginPath();
      ctx.moveTo(xs[a].x, xs[a].y);
      ctx.lineTo(xs[b].x, xs[b].y);
      ctx.lineTo(xs[c].x, xs[c].y);
      ctx.closePath();
      ctx.fill();
    }
    // Mesh wireframe — the FEM discretisation is part of the show.
    ctx.strokeStyle = 'rgba(190, 226, 240, 0.09)';
    ctx.lineWidth = 0.008;
    ctx.beginPath();
    for (const [a, b, c] of this.elements) {
      ctx.moveTo(xs[a].x, xs[a].y);
      ctx.lineTo(xs[b].x, xs[b].y);
      ctx.lineTo(xs[c].x, xs[c].y);
      ctx.closePath();
    }
    ctx.stroke();
    ctx.strokeStyle = 'rgba(233, 246, 252, 0.5)';
    ctx.lineWidth = 0.02;
    ctx.beginPath();
    ctx.arc(this.hole.x, this.hole.y, this.hole.r, 0, Math.PI * 2);
    ctx.stroke();
  }

  diagnostics() {
    let heat = 0;
    for (let i = 0; i < this.nodeCount; i += 1) heat += this.mass[i] * this.temperature[i];
    return {
      entities: this.nodeCount,
      contacts: this.elements.length,
      energy: heat,
      residual: this.lastResidual,
      method: `P1 三角形 FEM · 后向欧拉 + Jacobi · ${this.nodeCount} 节点`,
    };
  }
}

/**
 * 2D wave equation on a stretched membrane: classic leapfrog FDM with a
 * CFL-stable explicit scheme and light damping. Pointer plucks the surface.
 */
class WaveMembrane {
  constructor(width, height) {
    this.width = width;
    this.height = height;
    this.speed = 2.6;
    this.damping = 0.011;
    this.random = new Random(0x9a7e);
    this.level = 0;
    this.setQuality(2);
  }

  setQuality(level) {
    if (level === this.level) return;
    this.level = level;
    const n = { 1: 96, 2: 132, 3: 176 }[level] ?? 132;
    this.n = n;
    this.dx = this.width / (n - 1);
    this.size = n * n;
    this.hCur = new Float32Array(this.size);
    this.hPrev = new Float32Array(this.size);
    this.hNext = new Float32Array(this.size);
    this.painter = null;
    this.reset();
  }

  index(i, j) { return i + j * this.n; }

  reset() {
    this.hCur.fill(0);
    this.hPrev.fill(0);
    this.time = 0;
    this.dripTimer = 1.2;
    this.pluck(this.random.range(-3, 3), this.random.range(-2, 2), 0.65);
  }

  pluck(x, y, amplitude) {
    const gx = ((x + this.width / 2) / this.width) * (this.n - 1);
    const gy = ((y + this.height / 2) / this.height) * (this.n - 1);
    const radius = Math.max(3, this.n * 0.055);
    for (let j = Math.max(1, Math.floor(gy - radius)); j <= Math.min(this.n - 2, Math.ceil(gy + radius)); j += 1) {
      for (let i = Math.max(1, Math.floor(gx - radius)); i <= Math.min(this.n - 2, Math.ceil(gx + radius)); i += 1) {
        const d2 = (i - gx) ** 2 + (j - gy) ** 2;
        if (d2 > radius * radius) continue;
        this.hCur[this.index(i, j)] += amplitude * Math.exp(-d2 * 3.2 / (radius * radius));
      }
    }
  }

  step(dt) {
    this.time += dt;
    this.dripTimer -= dt;
    if (this.dripTimer <= 0) {
      this.dripTimer = 2.4;
      this.pluck(this.random.range(-4.5, 4.5), this.random.range(-2.6, 2.6), 0.5);
    }
    const r2 = (this.speed * dt / this.dx) ** 2;
    const n = this.n;
    for (let j = 1; j < n - 1; j += 1) {
      for (let i = 1; i < n - 1; i += 1) {
        const index = this.index(i, j);
        const lap = this.hCur[index - 1] + this.hCur[index + 1] + this.hCur[index - n] + this.hCur[index + n] - 4 * this.hCur[index];
        this.hNext[index] = (2 * this.hCur[index] - this.hPrev[index] + r2 * lap) - this.damping * (this.hCur[index] - this.hPrev[index]);
      }
    }
    [this.hPrev, this.hCur, this.hNext] = [this.hCur, this.hNext, this.hPrev];
  }

  pointer(kind, x, y) {
    if (kind === 'down') this.pluck(x, y, 0.75);
    else if (kind === 'move') this.pluck(x, y, 0.22);
  }

  render(ctx) {
    const painter = ensurePainter(this, this.n);
    if (!painter) return;
    const n = this.n;
    const pixels = painter.pixels;
    for (let j = 0; j < n; j += 1) {
      for (let i = 0; i < n; i += 1) {
        const index = this.index(i, j);
        const h = this.hCur[index];
        // Cheap Lambert shading from the local slope for a relief look.
        const sx = this.hCur[this.index(Math.min(i + 1, n - 1), j)] - this.hCur[this.index(Math.max(i - 1, 0), j)];
        const sy = this.hCur[this.index(i, Math.min(j + 1, n - 1))] - this.hCur[this.index(i, Math.max(j - 1, 0))];
        const light = clamp(0.66 - sx * 3.4 - sy * 2.6, 0.18, 1.5);
        const t = clamp(h * 1.1 + 0.5, 0, 1);
        const crest = clamp(h * 1.7, 0, 1);
        const r = Math.round((20 + 34 * t) * light + crest * 46);
        const g = Math.round((64 + 108 * t) * light + crest * 92);
        const b = Math.round((112 + 132 * t) * light + crest * 96);
        pixels[index] = pack(Math.min(r, 255), Math.min(g, 255), Math.min(b, 255), 255);
      }
    }
    painter.ctx.putImageData(painter.image, 0, 0);
    ctx.save();
    ctx.imageSmoothingEnabled = true;
    ctx.drawImage(painter.canvas, -this.width / 2, -this.height / 2, this.width, this.height);
    ctx.restore();
  }

  diagnostics() {
    let energy = 0;
    let peak = 0;
    for (let i = 0; i < this.size; i += 1) {
      energy += 0.5 * this.hCur[i] ** 2;
      peak = Math.max(peak, Math.abs(this.hCur[i]));
    }
    return {
      entities: this.size,
      contacts: Math.round(((this.speed * (1 / 120)) / this.dx) ** 2 * 100) / 100,
      energy,
      residual: peak,
      method: `波动方程 · 蛙跳有限差分 · ${this.n}² 网格`,
    };
  }
}

/**
 * Gray–Scott reaction–diffusion: two coupled PDEs whose Turing patterns
 * grow coral-like structures. Pointer paints fresh reagent V.
 */
class GrayScott {
  constructor(width, height) {
    this.width = width;
    this.height = height;
    this.Du = 0.16;
    this.Dv = 0.08;
    this.feed = 0.0545;
    this.kill = 0.062;
    this.random = new Random(0x6a11);
    this.level = 0;
    this.setQuality(2);
  }

  setQuality(level) {
    if (level === this.level) return;
    this.level = level;
    const n = { 1: 104, 2: 140, 3: 180 }[level] ?? 140;
    this.n = n;
    this.size = n * n;
    this.stepsPerDt = { 1: 6, 2: 8, 3: 10 }[level] ?? 8;
    this.u = new Float32Array(this.size);
    this.v = new Float32Array(this.size);
    this.uNext = new Float32Array(this.size);
    this.vNext = new Float32Array(this.size);
    this.painter = null;
    this.reset();
  }

  index(i, j) { return i + j * this.n; }

  reset() {
    this.u.fill(1);
    this.v.fill(0);
    for (let s = 0; s < 7; s += 1) {
      this.seed(this.random.range(-4.2, 4.2), this.random.range(-2.4, 2.4), this.random.range(0.28, 0.6));
    }
  }

  seed(x, y, radiusWorld = 0.4) {
    const gx = ((x + this.width / 2) / this.width) * (this.n - 1);
    const gy = ((y + this.height / 2) / this.height) * (this.n - 1);
    const radius = Math.max(2.5, (radiusWorld / this.width) * this.n);
    for (let j = Math.max(1, Math.floor(gy - radius)); j <= Math.min(this.n - 2, Math.ceil(gy + radius)); j += 1) {
      for (let i = Math.max(1, Math.floor(gx - radius)); i <= Math.min(this.n - 2, Math.ceil(gx + radius)); i += 1) {
        const d2 = (i - gx) ** 2 + (j - gy) ** 2;
        if (d2 > radius * radius) continue;
        const weight = 1 - Math.sqrt(d2) / radius;
        const index = this.index(i, j);
        this.v[index] = Math.min(this.v[index] + 0.5 * weight, 1);
        this.u[index] = Math.max(this.u[index] - 0.4 * weight, 0);
      }
    }
  }

  iterate() {
    const n = this.n;
    const { u, v, uNext, vNext, Du, Dv, feed, kill } = this;
    for (let j = 0; j < n; j += 1) {
      const jm = Math.max(j - 1, 0) * n;
      const jp = Math.min(j + 1, n - 1) * n;
      const j0 = j * n;
      for (let i = 0; i < n; i += 1) {
        const index = j0 + i;
        const im = Math.max(i - 1, 0);
        const ip = Math.min(i + 1, n - 1);
        const cu = u[index];
        const cv = v[index];
        const lapU = u[j0 + im] + u[j0 + ip] + u[jm + i] + u[jp + i] - 4 * cu;
        const lapV = v[j0 + im] + v[j0 + ip] + v[jm + i] + v[jp + i] - 4 * cv;
        const reaction = cu * cv * cv;
        uNext[index] = cu + Du * lapU - reaction + feed * (1 - cu);
        vNext[index] = cv + Dv * lapV + reaction - (feed + kill) * cv;
      }
    }
    [this.u, this.uNext] = [this.uNext, this.u];
    [this.v, this.vNext] = [this.vNext, this.v];
  }

  step() {
    for (let s = 0; s < this.stepsPerDt; s += 1) this.iterate();
  }

  pointer(kind, x, y) {
    if (kind === 'down' || kind === 'move') this.seed(x, y, 0.42);
  }

  render(ctx) {
    const painter = ensurePainter(this, this.n);
    if (!painter) return;
    const pixels = painter.pixels;
    for (let index = 0; index < this.size; index += 1) {
      const t = clamp(this.v[index] * 3.1, 0, 1);
      const ember = t * t;
      const r = Math.round(10 + 40 * t + 200 * ember);
      const g = Math.round(26 + 110 * t + 96 * ember);
      const b = Math.round(48 + 160 * t - 20 * ember);
      pixels[index] = pack(Math.min(r, 255), Math.min(g, 255), Math.max(Math.min(b, 255), 0), 255);
    }
    painter.ctx.putImageData(painter.image, 0, 0);
    ctx.save();
    ctx.imageSmoothingEnabled = true;
    ctx.drawImage(painter.canvas, -this.width / 2, -this.height / 2, this.width, this.height);
    ctx.restore();
  }

  diagnostics() {
    let total = 0;
    let peak = 0;
    for (let index = 0; index < this.size; index += 1) {
      total += this.v[index];
      peak = Math.max(peak, this.v[index]);
    }
    return {
      entities: this.size,
      contacts: this.stepsPerDt,
      energy: total / this.size,
      residual: peak,
      method: `Gray–Scott 反应扩散 · 显式 Euler · ${this.n}² 网格`,
    };
  }
}

/** Field-family dispatcher mirroring the fluid solver's structure. */
export class PdeLabSolver extends Solver {
  constructor({ width = 12, height = 7 } = {}) {
    super();
    this.mode = 'heat';
    this.heat = new HeatPlateFem(width, height);
    this.wave = new WaveMembrane(width, height);
    this.grayscott = new GrayScott(width, height);
    this.flock = new FlockSolver({ width, height });
  }

  current() {
    if (this.mode === 'wave') return this.wave;
    if (this.mode === 'grayscott') return this.grayscott;
    if (this.mode === 'flock') return this.flock;
    return this.heat;
  }

  setQuality(level) {
    this.heat.setQuality(level);
    this.wave.setQuality(level);
    this.grayscott.setQuality(level);
    this.flock.setQuality(level);
  }

  randomize() {
    const current = this.current();
    if (typeof current.randomize === 'function') current.randomize();
    else current.reset();
  }

  reset(mode = this.mode) { this.mode = mode; this.current().reset(); }
  step(dt) { this.current().step(dt); }
  render(ctx) { this.current().render(ctx); }
  pointer(kind, x, y) { this.current().pointer(kind, x, y); }
  diagnostics() { return this.current().diagnostics(); }
}

export { HeatPlateFem, WaveMembrane, GrayScott };
