import { clamp } from '../core/math.js';
import { Solver } from '../core/simulation.js';

function makeNode(x, y, mass = 1, pinned = false) {
  return { x, y, restX: x, restY: y, vx: 0, vy: 0, mass, invMass: pinned ? 0 : 1 / mass, baseInvMass: pinned ? 0 : 1 / mass, fx: 0, fy: 0 };
}

function makeTriangle(i0, i1, i2, nodes) {
  const a = nodes[i0]; const b = nodes[i1]; const c = nodes[i2];
  const d00 = b.restX - a.restX;
  const d01 = c.restX - a.restX;
  const d10 = b.restY - a.restY;
  const d11 = c.restY - a.restY;
  const determinant = d00 * d11 - d01 * d10;
  const area = Math.abs(determinant) * 0.5;
  return {
    i0, i1, i2, area,
    inv00: d11 / determinant, inv01: -d01 / determinant,
    inv10: -d10 / determinant, inv11: d00 / determinant,
    strain: 0,
  };
}

/**
 * Corotated linear triangular FEM, integrated with small explicit symplectic
 * substeps. It intentionally complements XPBD: the UI can compare an energy
 * based element force model with a constraint projection model.
 */
export class FemSolver extends Solver {
  constructor({ width = 12, height = 7 } = {}) {
    super();
    this.width = width;
    this.height = height;
    this.gravity = -9.81;
    this.young = 1500;
    this.poisson = 0.28;
    this.damping = 0.998;
    this.substeps = 5;
    this.nodes = [];
    this.triangles = [];
    this.drag = null;
    this.obstacle = { x: 0.7, y: -2.48, radius: 0.72 };
    this.maxStrain = 0;
    this.reset();
  }

  get material() {
    const mu = this.young / (2 * (1 + this.poisson));
    const lambda = (this.young * this.poisson) / ((1 + this.poisson) * (1 - 2 * this.poisson));
    return { mu, lambda };
  }

  reset() {
    this.nodes = [];
    this.triangles = [];
    this.drag = null;
    const columns = 17;
    const rows = 8;
    const dx = 0.37;
    const dy = 0.37;
    const originX = -4.95;
    const originY = 1.65;
    const index = (x, y) => y * columns + x;
    for (let y = 0; y < rows; y += 1) {
      for (let x = 0; x < columns; x += 1) this.nodes.push(makeNode(originX + x * dx, originY - y * dy, 0.12, x === 0));
    }
    for (let y = 0; y < rows - 1; y += 1) {
      for (let x = 0; x < columns - 1; x += 1) {
        const a = index(x, y); const b = index(x + 1, y); const c = index(x, y + 1); const d = index(x + 1, y + 1);
        this.triangles.push(makeTriangle(a, b, c, this.nodes), makeTriangle(b, d, c, this.nodes));
      }
    }
    // Assign mass after assembling elements so mass follows the rest area.
    for (const node of this.nodes) node.mass = 0;
    for (const triangle of this.triangles) {
      for (const indexNode of [triangle.i0, triangle.i1, triangle.i2]) this.nodes[indexNode].mass += triangle.area * 2.1 / 3;
    }
    for (const node of this.nodes) {
      node.invMass = node.baseInvMass === 0 ? 0 : 1 / node.mass;
      node.baseInvMass = node.invMass;
    }
  }

  addTriangleForce(triangle, mu, lambda) {
    const a = this.nodes[triangle.i0]; const b = this.nodes[triangle.i1]; const c = this.nodes[triangle.i2];
    const ds00 = b.x - a.x; const ds01 = c.x - a.x;
    const ds10 = b.y - a.y; const ds11 = c.y - a.y;
    const f00 = ds00 * triangle.inv00 + ds01 * triangle.inv10;
    const f01 = ds00 * triangle.inv01 + ds01 * triangle.inv11;
    const f10 = ds10 * triangle.inv00 + ds11 * triangle.inv10;
    const f11 = ds10 * triangle.inv01 + ds11 * triangle.inv11;
    // Closest 2D rotation of F, avoiding a full SVD for this teaching-scale solver.
    const angle = Math.atan2(f10 - f01, f00 + f11);
    const cos = Math.cos(angle); const sin = Math.sin(angle);
    const trace = (cos * f00 + sin * f10 - 1) + (-sin * f01 + cos * f11 - 1);
    const p00 = 2 * mu * (f00 - cos) + lambda * trace * cos;
    const p01 = 2 * mu * (f01 + sin) - lambda * trace * sin;
    const p10 = 2 * mu * (f10 - sin) + lambda * trace * sin;
    const p11 = 2 * mu * (f11 - cos) + lambda * trace * cos;
    const h00 = -triangle.area * (p00 * triangle.inv00 + p01 * triangle.inv01);
    const h01 = -triangle.area * (p00 * triangle.inv10 + p01 * triangle.inv11);
    const h10 = -triangle.area * (p10 * triangle.inv00 + p11 * triangle.inv01);
    const h11 = -triangle.area * (p10 * triangle.inv10 + p11 * triangle.inv11);
    b.fx += h00; b.fy += h10;
    c.fx += h01; c.fy += h11;
    a.fx -= h00 + h01; a.fy -= h10 + h11;
    const detF = f00 * f11 - f01 * f10;
    triangle.strain = Math.min(2, Math.sqrt((f00 - cos) ** 2 + (f01 + sin) ** 2 + (f10 - sin) ** 2 + (f11 - cos) ** 2));
    return Math.abs(detF - 1);
  }

  collide(node) {
    const margin = 0.035;
    const floor = -this.height / 2 + margin;
    const left = -this.width / 2 + margin;
    const right = this.width / 2 - margin;
    if (node.y < floor) { node.y = floor; if (node.vy < 0) node.vy *= -0.13; node.vx *= 0.91; }
    if (node.x < left) { node.x = left; if (node.vx < 0) node.vx *= -0.12; }
    if (node.x > right) { node.x = right; if (node.vx > 0) node.vx *= -0.12; }
    const dx = node.x - this.obstacle.x; const dy = node.y - this.obstacle.y;
    const distance = Math.hypot(dx, dy);
    if (distance < this.obstacle.radius + margin) {
      const nx = distance > 1e-8 ? dx / distance : 1; const ny = distance > 1e-8 ? dy / distance : 0;
      node.x = this.obstacle.x + nx * (this.obstacle.radius + margin);
      node.y = this.obstacle.y + ny * (this.obstacle.radius + margin);
      const normalVelocity = node.vx * nx + node.vy * ny;
      if (normalVelocity < 0) { node.vx -= normalVelocity * 1.08 * nx; node.vy -= normalVelocity * 1.08 * ny; }
    }
  }

  step(dt) {
    const { mu, lambda } = this.material;
    const h = dt / this.substeps;
    this.maxStrain = 0;
    for (let substep = 0; substep < this.substeps; substep += 1) {
      for (const node of this.nodes) { node.fx = 0; node.fy = node.mass * this.gravity; }
      for (const triangle of this.triangles) this.maxStrain = Math.max(this.maxStrain, this.addTriangleForce(triangle, mu, lambda));
      for (const node of this.nodes) {
        if (node.invMass === 0) continue;
        node.vx = clamp((node.vx + node.fx * node.invMass * h) * this.damping, -35, 35);
        node.vy = clamp((node.vy + node.fy * node.invMass * h) * this.damping, -35, 35);
        node.x += node.vx * h; node.y += node.vy * h;
        this.collide(node);
      }
      if (this.drag) {
        const node = this.nodes[this.drag.index];
        node.x = this.drag.x; node.y = this.drag.y; node.vx = 0; node.vy = 0;
      }
    }
  }

  pointer(kind, x, y) {
    if (kind === 'down') {
      let selected = -1; let distanceBest = 0.4;
      for (let i = 0; i < this.nodes.length; i += 1) {
        const node = this.nodes[i]; const distance = Math.hypot(node.x - x, node.y - y);
        if (distance < distanceBest) { selected = i; distanceBest = distance; }
      }
      if (selected >= 0) {
        const node = this.nodes[selected]; this.drag = { index: selected, x, y, restore: node.invMass }; node.invMass = 0;
      }
    } else if (kind === 'move' && this.drag) {
      this.drag.x = x; this.drag.y = y;
    } else if ((kind === 'up' || kind === 'leave') && this.drag) {
      const node = this.nodes[this.drag.index]; node.invMass = this.drag.restore; this.drag = null;
    }
  }

  render(ctx) {
    for (const triangle of this.triangles) {
      const a = this.nodes[triangle.i0]; const b = this.nodes[triangle.i1]; const c = this.nodes[triangle.i2];
      const color = clamp(triangle.strain * 8, 0, 1);
      ctx.fillStyle = `rgba(${Math.round(69 + 182 * color)}, ${Math.round(190 - 86 * color)}, ${Math.round(211 - 136 * color)}, 0.72)`;
      ctx.beginPath(); ctx.moveTo(a.x, a.y); ctx.lineTo(b.x, b.y); ctx.lineTo(c.x, c.y); ctx.closePath(); ctx.fill();
      ctx.strokeStyle = 'rgba(214, 247, 251, 0.17)'; ctx.lineWidth = 0.008; ctx.stroke();
    }
    ctx.fillStyle = '#203f54'; ctx.beginPath(); ctx.arc(this.obstacle.x, this.obstacle.y, this.obstacle.radius, 0, Math.PI * 2); ctx.fill();
    ctx.fillStyle = '#e6fbff';
    for (const node of this.nodes) { if (node.invMass === 0) { ctx.beginPath(); ctx.arc(node.x, node.y, 0.04, 0, Math.PI * 2); ctx.fill(); } }
  }

  diagnostics() {
    let energy = 0;
    for (const node of this.nodes) energy += 0.5 * node.mass * (node.vx ** 2 + node.vy ** 2);
    return { entities: this.nodes.length, contacts: this.triangles.length, energy, residual: this.maxStrain, method: '共旋线性三角 FEM · 半隐式积分' };
  }
}
