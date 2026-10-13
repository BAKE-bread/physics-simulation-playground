import { clamp } from '../core/math.js';
import { Solver } from '../core/simulation.js';

// D2Q9 lattice constants.
const EX = [0, 1, 0, -1, 0, 1, -1, -1, 1];
const EY = [0, 0, 1, 0, -1, 1, 1, -1, -1];
const OPP = [0, 3, 4, 1, 2, 7, 8, 5, 6];
const W = [4 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 36, 1 / 36, 1 / 36, 1 / 36];

/** Colour packer for the little-endian Uint32 ImageData layout. */
const pack = (r, g, b, a) => ((a << 24) | (b << 16) | (g << 8) | r) >>> 0;

/**
 * Lattice-Boltzmann (D2Q9, BGK) channel flow past a circular cylinder.
 * The classic gateway CFD demo: a von Kármán vortex street develops behind
 * the obstacle at moderate Reynolds number. Combined pull-scheme streaming
 * and collision keeps the inner loop allocation-free.
 */
export class LbmSolver extends Solver {
  constructor(width = 12, height = 7) {
    super();
    this.width = width;
    this.height = height;
    this.u0 = 0.085;          // lattice inlet speed; keep < ~0.1 for low Mach
    this.reynolds = 160;
    this.obstacle = { x: -2.6, y: 0 };
    this.dragging = false;
    this.lastPointer = null;
    this.painter = null;
    this.level = 0;
    this.setQuality(2);
  }

  /** Quality level 1..3 controls lattice size and steps per fixed dt. */
  setQuality(level) {
    if (level === this.level) return;
    this.level = level;
    const presets = { 1: [150, 88, 3], 2: [190, 110, 4], 3: [230, 134, 5] };
    const [nx, ny, steps] = presets[level] ?? presets[2];
    this.nx = nx;
    this.ny = ny;
    this.stepsPerDt = steps;
    this.radiusLattice = ny / 10;
    const nu = (this.u0 * this.radiusLattice) / this.reynolds;
    this.tau = 0.5 + 3 * nu;
    const size = nx * ny;
    this.size = size;
    this.f1 = new Float32Array(9 * size);
    this.f2 = new Float32Array(9 * size);
    this.ux = new Float32Array(size);
    this.uy = new Float32Array(size);
    this.rho = new Float32Array(size);
    this.barrier = new Uint8Array(size);
    this.brushW = new Float32Array(size);
    this.brushIndices = [];
    this.brushForce = { x: 0, y: 0 };
    this.painter = null;
    this.reset();
  }

  reset() {
    // Uniform inlet flow as initial condition.
    for (let i = 0; i < this.size; i += 1) this.setEquilibrium(this.f1, i, 1, this.u0, 0);
    this.buildBarrier();
  }

  setEquilibrium(f, i, rho, ux, uy) {
    const size = this.size;
    const usq = ux * ux + uy * uy;
    for (let q = 0; q < 9; q += 1) {
      const eu = EX[q] * ux + EY[q] * uy;
      f[q * size + i] = W[q] * rho * (1 + 3 * eu + 4.5 * eu * eu - 1.5 * usq);
    }
  }

  buildBarrier() {
    const cx = ((this.obstacle.x + this.width / 2) / this.width) * (this.nx - 1);
    const cy = ((this.obstacle.y + this.height / 2) / this.height) * (this.ny - 1);
    const r = this.radiusLattice;
    let count = 0;
    for (let y = 0; y < this.ny; y += 1) {
      for (let x = 0; x < this.nx; x += 1) {
        const i = x + y * this.nx;
        const inside = (x - cx) ** 2 + (y - cy) ** 2 <= r * r;
        const was = this.barrier[i] === 1;
        if (inside) { this.barrier[i] = 1; count += 1; }
        else if (was) {
          // A freshly-freed cell gets a sane local state instead of stale f.
          this.barrier[i] = 0;
          this.setEquilibrium(this.f1, i, 1, this.u0, 0);
        }
      }
    }
    this.barrierCount = count;
  }

  moveObstacle(x, y) {
    this.obstacle.x = clamp(x, -this.width / 2 + 1, this.width / 2 - 1);
    this.obstacle.y = clamp(y, -this.height / 2 + 0.9, this.height / 2 - 0.9);
    this.buildBarrier();
  }

  latticeStep() {
    const { nx, ny, size, f1, f2, barrier, tau, ux, uy, rho, u0 } = this;
    const invTau = 1 / tau;
    for (let y = 0; y < ny; y += 1) {
      for (let x = 0; x < nx; x += 1) {
        const i = x + y * nx;
        if (barrier[i]) continue;
        if (x === 0) {
          this.setEquilibrium(f2, i, 1, u0, 0);
          ux[i] = u0; uy[i] = 0; rho[i] = 1;
          continue;
        }
        if (x === nx - 1) {
          for (let q = 0; q < 9; q += 1) f2[q * size + i] = f1[q * size + i - 1];
          continue;
        }
        // Pull streaming with bounce-back at walls and the barrier surface.
        let f0; let fA; let fB; let fC; let fD; let fE; let fF; let fG; let fH;
        f0 = f1[i];
        fA = pull(f1, 1, x - 1, y, i, nx, ny, size, barrier);
        fB = pull(f1, 2, x, y - 1, i, nx, ny, size, barrier);
        fC = pull(f1, 3, x + 1, y, i, nx, ny, size, barrier);
        fD = pull(f1, 4, x, y + 1, i, nx, ny, size, barrier);
        fE = pull(f1, 5, x - 1, y - 1, i, nx, ny, size, barrier);
        fF = pull(f1, 6, x + 1, y - 1, i, nx, ny, size, barrier);
        fG = pull(f1, 7, x + 1, y + 1, i, nx, ny, size, barrier);
        fH = pull(f1, 8, x - 1, y + 1, i, nx, ny, size, barrier);
        const r = f0 + fA + fB + fC + fD + fE + fF + fG + fH;
        let vx = (fA - fC + fE - fF - fG + fH) / r;
        let vy = (fB - fD + fE + fF - fG - fH) / r;
        // Pointer brush adds a gentle body force instead of resetting local
        // populations, so existing vortices are perturbed, never wiped out.
        const bw = this.brushW[i];
        if (bw !== 0) {
          vx += this.brushForce.x * bw;
          vy += this.brushForce.y * bw;
          const speed = Math.hypot(vx, vy);
          // Keep the local Mach number well below the BGK stability ceiling.
          if (speed > 0.16) { vx = (vx / speed) * 0.16; vy = (vy / speed) * 0.16; }
        }
        rho[i] = r; ux[i] = vx; uy[i] = vy;
        const usq = vx * vx + vy * vy;
        f2[i] = relax(f0, W[0] * r * (1 - 1.5 * usq), invTau);
        f2[size + i] = relax(fA, W[1] * r * (1 + 3 * vx + 4.5 * vx * vx - 1.5 * usq), invTau);
        f2[2 * size + i] = relax(fB, W[2] * r * (1 + 3 * vy + 4.5 * vy * vy - 1.5 * usq), invTau);
        f2[3 * size + i] = relax(fC, W[3] * r * (1 - 3 * vx + 4.5 * vx * vx - 1.5 * usq), invTau);
        f2[4 * size + i] = relax(fD, W[4] * r * (1 - 3 * vy + 4.5 * vy * vy - 1.5 * usq), invTau);
        f2[5 * size + i] = relax(fE, W[5] * r * (1 + 3 * (vx + vy) + 4.5 * (vx + vy) ** 2 - 1.5 * usq), invTau);
        f2[6 * size + i] = relax(fF, W[6] * r * (1 + 3 * (-vx + vy) + 4.5 * (-vx + vy) ** 2 - 1.5 * usq), invTau);
        f2[7 * size + i] = relax(fG, W[7] * r * (1 - 3 * (vx + vy) + 4.5 * (vx + vy) ** 2 - 1.5 * usq), invTau);
        f2[8 * size + i] = relax(fH, W[8] * r * (1 + 3 * (vx - vy) + 4.5 * (vx - vy) ** 2 - 1.5 * usq), invTau);
      }
    }
    [this.f1, this.f2] = [this.f2, this.f1];
  }

  /** Soft force brush: nudges velocities inside a compact stencil. */
  updateBrush(wx, wy, dvx, dvy) {
    this.clearBrush();
    const cx = ((wx + this.width / 2) / this.width) * (this.nx - 1);
    const cy = ((wy + this.height / 2) / this.height) * (this.ny - 1);
    const radius = this.ny * 0.06;
    const scaleU = (this.nx - 1) / this.width; // world → lattice velocity
    this.brushForce.x = clamp(dvx * scaleU * 0.045, -0.016, 0.016);
    this.brushForce.y = clamp(dvy * scaleU * 0.045, -0.016, 0.016);
    for (let y = Math.max(1, Math.floor(cy - radius)); y <= Math.min(this.ny - 2, Math.ceil(cy + radius)); y += 1) {
      for (let x = Math.max(1, Math.floor(cx - radius)); x <= Math.min(this.nx - 2, Math.ceil(cx + radius)); x += 1) {
        const i = x + y * this.nx;
        if (this.barrier[i]) continue;
        const d2 = (x - cx) ** 2 + (y - cy) ** 2;
        if (d2 > radius * radius) continue;
        this.brushW[i] = (1 - Math.sqrt(d2) / radius) * 0.9;
        this.brushIndices.push(i);
      }
    }
  }

  clearBrush() {
    for (const i of this.brushIndices) this.brushW[i] = 0;
    this.brushIndices = [];
    this.brushForce.x = 0;
    this.brushForce.y = 0;
  }

  step() {
    for (let s = 0; s < this.stepsPerDt; s += 1) this.latticeStep();
  }

  pointer(kind, x, y) {
    const radiusWorld = (this.radiusLattice / (this.ny - 1)) * this.height;
    if (kind === 'down') {
      const near = Math.hypot(x - this.obstacle.x, y - this.obstacle.y) < radiusWorld * 1.5;
      this.dragging = near;
      if (near) this.moveObstacle(x, y);
      this.lastPointer = { x, y };
    } else if (kind === 'move' && this.lastPointer) {
      if (this.dragging) this.moveObstacle(x, y);
      else this.updateBrush(x, y, x - this.lastPointer.x, y - this.lastPointer.y);
      this.lastPointer = { x, y };
    } else if (kind === 'up' || kind === 'leave') {
      this.dragging = false;
      this.lastPointer = null;
      this.clearBrush();
    }
  }

  render(ctx) {
    const painter = ensurePainter(this);
    if (painter) {
      const { nx, ny, ux, uy, barrier } = this;
      const pixels = painter.pixels;
      for (let y = 0; y < ny; y += 1) {
        for (let x = 0; x < nx; x += 1) {
          const i = x + y * nx;
          if (barrier[i]) { pixels[i] = pack(16, 34, 46, 255); continue; }
          // Vorticity (curl) colouring makes the vortex street pop.
          const xm = Math.max(x - 1, 1); const xp = Math.min(x + 1, nx - 2);
          const ym = Math.max(y - 1, 1); const yp = Math.min(y + 1, ny - 2);
          const curl = (uy[xp + y * nx] - uy[xm + y * nx]) - (ux[x + yp * nx] - ux[x + ym * nx]);
          const s = clamp(curl * 26, -1, 1);
          const a = Math.abs(s);
          const wake = clamp((ux[i] / this.u0 - 1) * 1.4, -1, 1);
          let r = 10 + 12 * Math.max(wake, 0);
          let g = 30 + 10 * Math.max(wake, 0);
          let b = 46 + 14 * Math.max(wake, 0);
          if (s > 0) { r += 30 * a; g += 190 * a; b += 220 * a; }
          else { r += 225 * a; g += 120 * a; b += 42 * a; }
          pixels[i] = pack(Math.min(Math.round(r), 255), Math.min(Math.round(g), 255), Math.min(Math.round(b), 255), 255);
        }
      }
      painter.ctx.putImageData(painter.image, 0, 0);
      ctx.save();
      ctx.imageSmoothingEnabled = true;
      ctx.drawImage(painter.canvas, -this.width / 2, -this.height / 2, this.width, this.height);
      ctx.restore();
    }
    // Obstacle outline with a soft glow, so the draggable target reads clearly.
    const radiusWorld = (this.radiusLattice / (this.ny - 1)) * this.height;
    ctx.save();
    ctx.strokeStyle = 'rgba(233, 246, 252, 0.85)';
    ctx.fillStyle = 'rgba(13, 32, 44, 0.55)';
    ctx.lineWidth = 0.028;
    ctx.shadowColor = 'rgba(116, 215, 238, 0.8)';
    ctx.shadowBlur = 14;
    ctx.beginPath();
    ctx.arc(this.obstacle.x, this.obstacle.y, radiusWorld, 0, Math.PI * 2);
    ctx.fill();
    ctx.stroke();
    ctx.restore();
  }

  diagnostics() {
    let energy = 0;
    let densityError = 0;
    let count = 0;
    for (let i = 0; i < this.size; i += 1) {
      if (this.barrier[i]) continue;
      energy += 0.5 * (this.ux[i] ** 2 + this.uy[i] ** 2);
      densityError += Math.abs(this.rho[i] - 1);
      count += 1;
    }
    return {
      entities: this.size,
      contacts: this.barrierCount,
      energy: count ? energy / count : 0,
      residual: count ? densityError / count : 0,
      method: `LBM D2Q9 · BGK · Re≈${this.reynolds} · ${this.nx}×${this.ny} 格点`,
    };
  }
}

function pull(f, q, x, y, i, nx, ny, size, barrier) {
  if (y < 0 || y >= ny || barrier[x + y * nx]) return f[OPP[q] * size + i];
  return f[q * size + x + y * nx];
}

function relax(value, equilibrium, invTau) {
  return value + invTau * (equilibrium - value);
}

function ensurePainter(owner) {
  if (typeof document === 'undefined') return null;
  if (!owner.painter) {
    const canvas = document.createElement('canvas');
    owner.painter = { canvas, ctx: canvas.getContext('2d'), image: null, pixels: null, nx: 0, ny: 0 };
  }
  const painter = owner.painter;
  if (painter.nx !== owner.nx || painter.ny !== owner.ny) {
    painter.nx = owner.nx;
    painter.ny = owner.ny;
    painter.canvas.width = owner.nx;
    painter.canvas.height = owner.ny;
    painter.image = painter.ctx.createImageData(owner.nx, owner.ny);
    painter.pixels = new Uint32Array(painter.image.data.buffer);
  }
  return painter;
}
