import { clamp, length, Random } from '../core/math.js';
import { Solver } from '../core/simulation.js';

function particle(x, y, mass = 1, pinned = false) {
  return { x, y, px: x, py: y, vx: 0, vy: 0, invMass: pinned ? 0 : 1 / mass, baseInvMass: pinned ? 0 : 1 / mass, pinned };
}

function distanceConstraint(i, j, rest, compliance = 2e-6, kind = 'stretch') {
  return { i, j, rest, compliance, lambda: 0, kind };
}

function polygonArea(points, indices) {
  let total = 0;
  for (let k = 0; k < indices.length; k += 1) {
    const a = points[indices[k]];
    const b = points[indices[(k + 1) % indices.length]];
    total += a.x * b.y - a.y * b.x;
  }
  return total * 0.5;
}

/** XPBD distance projection; compliance has physical time-step scaling. */
function solveDistance(points, constraint, dt) {
  const a = points[constraint.i];
  const b = points[constraint.j];
  const dx = b.x - a.x;
  const dy = b.y - a.y;
  const distance = length(dx, dy);
  if (distance < 1e-8) return 0;
  const nx = dx / distance;
  const ny = dy / distance;
  const c = distance - constraint.rest;
  const alpha = constraint.compliance / (dt * dt);
  const denominator = a.invMass + b.invMass + alpha;
  if (denominator < 1e-10) return Math.abs(c);
  const deltaLambda = (-c - alpha * constraint.lambda) / denominator;
  constraint.lambda += deltaLambda;
  a.x -= a.invMass * deltaLambda * nx;
  a.y -= a.invMass * deltaLambda * ny;
  b.x += b.invMass * deltaLambda * nx;
  b.y += b.invMass * deltaLambda * ny;
  return Math.abs(c);
}

function solveArea(points, constraint, dt) {
  const indices = constraint.indices;
  const gradients = [];
  let denominator = 0;
  for (let k = 0; k < indices.length; k += 1) {
    const prev = points[indices[(k - 1 + indices.length) % indices.length]];
    const next = points[indices[(k + 1) % indices.length]];
    const current = points[indices[k]];
    const gx = 0.5 * (next.y - prev.y);
    const gy = 0.5 * (prev.x - next.x);
    gradients.push({ gx, gy });
    denominator += current.invMass * (gx * gx + gy * gy);
  }
  const c = polygonArea(points, indices) - constraint.rest;
  const alpha = constraint.compliance / (dt * dt);
  denominator += alpha;
  if (denominator < 1e-10) return Math.abs(c);
  const deltaLambda = (-c - alpha * constraint.lambda) / denominator;
  constraint.lambda += deltaLambda;
  for (let k = 0; k < indices.length; k += 1) {
    const p = points[indices[k]];
    p.x += p.invMass * deltaLambda * gradients[k].gx;
    p.y += p.invMass * deltaLambda * gradients[k].gy;
  }
  return Math.abs(c);
}

export class SoftSolver extends Solver {
  constructor({ width = 12, height = 7 } = {}) {
    super();
    this.width = width;
    this.height = height;
    this.gravity = -9.81;
    this.iterations = 10;
    this.substeps = 2;
    this.compliance = 2e-6;
    this.points = [];
    this.constraints = [];
    this.areaConstraint = null;
    this.triangles = [];
    this.mode = 'cloth';
    this.drag = null;
    this.random = new Random(0x51f7);
    this.maxError = 0;
    this.reset();
  }

  reset(mode = this.mode) {
    this.mode = mode;
    this.points = [];
    this.constraints = [];
    this.triangles = [];
    this.areaConstraint = null;
    this.drag = null;
    if (mode === 'balloon') this.buildBalloon();
    else this.buildCloth();
  }

  buildCloth() {
    const cols = 17;
    const rows = 11;
    const spacing = 0.43;
    const originX = -3.45;
    const originY = 2.55;
    const at = (x, y) => y * cols + x;
    for (let y = 0; y < rows; y += 1) {
      for (let x = 0; x < cols; x += 1) {
        this.points.push(particle(originX + x * spacing, originY - y * spacing, 0.55, y === 0 && (x === 0 || x === cols - 1)));
      }
    }
    const add = (a, b, rest, compliance, kind) => this.constraints.push(distanceConstraint(a, b, rest, compliance, kind));
    for (let y = 0; y < rows; y += 1) {
      for (let x = 0; x < cols; x += 1) {
        if (x + 1 < cols) add(at(x, y), at(x + 1, y), spacing, this.compliance, 'stretch');
        if (y + 1 < rows) add(at(x, y), at(x, y + 1), spacing, this.compliance, 'stretch');
        if (x + 1 < cols && y + 1 < rows) {
          add(at(x, y), at(x + 1, y + 1), spacing * Math.SQRT2, this.compliance * 1.5, 'shear');
          add(at(x + 1, y), at(x, y + 1), spacing * Math.SQRT2, this.compliance * 1.5, 'shear');
          this.triangles.push([at(x, y), at(x + 1, y), at(x, y + 1)], [at(x + 1, y), at(x + 1, y + 1), at(x, y + 1)]);
        }
        if (x + 2 < cols) add(at(x, y), at(x + 2, y), spacing * 2, this.compliance * 8, 'bend');
        if (y + 2 < rows) add(at(x, y), at(x, y + 2), spacing * 2, this.compliance * 8, 'bend');
      }
    }
    // Start below the hanging mesh; otherwise an initial collision projection
    // would tear apart the unrelaxed cloth before gravity has acted.
    this.obstacle = { x: 1.1, y: -2.55, radius: 0.78 };
  }

  buildBalloon() {
    const count = 24;
    const center = { x: -2.2, y: 1.65 };
    const radius = 1.12;
    for (let i = 0; i < count; i += 1) {
      const angle = (i / count) * Math.PI * 2;
      this.points.push(particle(center.x + Math.cos(angle) * radius, center.y + Math.sin(angle) * radius, 0.45));
    }
    const edge = 2 * radius * Math.sin(Math.PI / count);
    const bend = 2 * radius * Math.sin((2 * Math.PI) / count);
    for (let i = 0; i < count; i += 1) {
      this.constraints.push(distanceConstraint(i, (i + 1) % count, edge, this.compliance, 'stretch'));
      this.constraints.push(distanceConstraint(i, (i + 2) % count, bend, this.compliance * 5, 'bend'));
    }
    this.areaConstraint = { indices: Array.from({ length: count }, (_, index) => index), rest: Math.PI * radius * radius, compliance: 8e-8, lambda: 0 };
    this.obstacle = { x: 1.2, y: -1.65, radius: 0.95 };
    this.points[5].vx = 2.5;
  }

  projectBounds(p) {
    const margin = 0.08;
    const left = -this.width / 2 + margin;
    const right = this.width / 2 - margin;
    const floor = -this.height / 2 + margin;
    if (p.x < left) p.x = left;
    if (p.x > right) p.x = right;
    if (p.y < floor) p.y = floor;
    const dx = p.x - this.obstacle.x;
    const dy = p.y - this.obstacle.y;
    const distance = Math.hypot(dx, dy);
    const target = this.obstacle.radius + margin;
    if (distance < target) {
      const nx = distance > 1e-7 ? dx / distance : 1;
      const ny = distance > 1e-7 ? dy / distance : 0;
      p.x = this.obstacle.x + nx * target;
      p.y = this.obstacle.y + ny * target;
    }
  }

  step(dt) {
    const h = dt / this.substeps;
    this.maxError = 0;
    for (let substep = 0; substep < this.substeps; substep += 1) {
      for (const constraint of this.constraints) constraint.lambda = 0;
      if (this.areaConstraint) this.areaConstraint.lambda = 0;
      for (const p of this.points) {
        if (p.invMass === 0) continue;
        p.px = p.x;
        p.py = p.y;
        p.vy += this.gravity * h;
        p.vx *= 0.999;
        p.vy *= 0.999;
        p.x += p.vx * h;
        p.y += p.vy * h;
      }
      if (this.drag) {
        const p = this.points[this.drag.index];
        p.x = this.drag.x;
        p.y = this.drag.y;
      }
      for (let iteration = 0; iteration < this.iterations; iteration += 1) {
        for (const constraint of this.constraints) this.maxError = Math.max(this.maxError, solveDistance(this.points, constraint, h));
        if (this.areaConstraint) this.maxError = Math.max(this.maxError, solveArea(this.points, this.areaConstraint, h));
        for (const p of this.points) this.projectBounds(p);
        if (this.drag) {
          const p = this.points[this.drag.index];
          p.x = this.drag.x;
          p.y = this.drag.y;
        }
      }
      for (const p of this.points) {
        if (p.invMass === 0) continue;
        p.vx = (p.x - p.px) / h;
        p.vy = (p.y - p.py) / h;
      }
    }
    // Report the post-projection constraint residual, not the larger value
    // observed while a Gauss--Seidel sweep was still converging.
    this.maxError = 0;
    for (const constraint of this.constraints) {
      const a = this.points[constraint.i];
      const b = this.points[constraint.j];
      this.maxError = Math.max(this.maxError, Math.abs(Math.hypot(b.x - a.x, b.y - a.y) - constraint.rest));
    }
  }

  pointer(kind, x, y) {
    if (kind === 'down') {
      // Grabbing inside the obstacle ball moves the ball itself.
      const obstacleDistance = Math.hypot(this.obstacle.x - x, this.obstacle.y - y);
      if (obstacleDistance < this.obstacle.radius * 0.95) {
        this.obstacleDrag = { offsetX: this.obstacle.x - x, offsetY: this.obstacle.y - y };
        return;
      }
      let index = -1;
      let best = 0.55;
      for (let i = 0; i < this.points.length; i += 1) {
        const p = this.points[i];
        const distance = Math.hypot(p.x - x, p.y - y);
        if (distance < best) { best = distance; index = i; }
      }
      if (index >= 0) {
        const p = this.points[index];
        this.drag = { index, x, y, restore: p.invMass };
        p.invMass = 0;
      }
    } else if (kind === 'move') {
      if (this.obstacleDrag) {
        this.obstacle.x = clamp(x + this.obstacleDrag.offsetX, -this.width / 2 + this.obstacle.radius, this.width / 2 - this.obstacle.radius);
        this.obstacle.y = clamp(y + this.obstacleDrag.offsetY, -this.height / 2 + this.obstacle.radius, this.height / 2 - this.obstacle.radius);
        return;
      }
      if (this.drag) {
        this.drag.x = x;
        this.drag.y = y;
      }
    } else if (kind === 'up' || kind === 'leave') {
      if (this.drag) {
        const p = this.points[this.drag.index];
        p.invMass = this.drag.restore;
        this.drag = null;
      }
      this.obstacleDrag = null;
    }
  }

  render(ctx) {
    if (this.mode === 'cloth') {
      ctx.fillStyle = 'rgba(79, 196, 209, 0.13)';
      for (const triangle of this.triangles) {
        const a = this.points[triangle[0]];
        const b = this.points[triangle[1]];
        const c = this.points[triangle[2]];
        ctx.beginPath(); ctx.moveTo(a.x, a.y); ctx.lineTo(b.x, b.y); ctx.lineTo(c.x, c.y); ctx.closePath(); ctx.fill();
      }
    } else {
      ctx.fillStyle = 'rgba(79, 196, 209, 0.56)';
      ctx.beginPath();
      this.points.forEach((p, index) => index === 0 ? ctx.moveTo(p.x, p.y) : ctx.lineTo(p.x, p.y));
      ctx.closePath();
      ctx.fill();
    }
    ctx.fillStyle = '#203f54';
    ctx.beginPath();
    ctx.arc(this.obstacle.x, this.obstacle.y, this.obstacle.radius, 0, Math.PI * 2);
    ctx.fill();
    // Dashed rim signals the obstacle is draggable.
    ctx.strokeStyle = 'rgba(186, 246, 255, 0.42)';
    ctx.lineWidth = 0.02;
    ctx.setLineDash([0.1, 0.1]);
    ctx.stroke();
    ctx.setLineDash([]);
    ctx.strokeStyle = 'rgba(186, 246, 255, 0.75)';
    ctx.lineWidth = 0.018;
    for (const constraint of this.constraints) {
      if (constraint.kind === 'bend' && this.mode === 'cloth') continue;
      const a = this.points[constraint.i];
      const b = this.points[constraint.j];
      const strain = clamp(Math.abs(Math.hypot(b.x - a.x, b.y - a.y) / constraint.rest - 1) * 12, 0, 1);
      ctx.strokeStyle = `rgba(${Math.round(116 + strain * 130)}, ${Math.round(217 - strain * 115)}, ${Math.round(238 - strain * 150)}, 0.82)`;
      ctx.beginPath(); ctx.moveTo(a.x, a.y); ctx.lineTo(b.x, b.y); ctx.stroke();
    }
    ctx.fillStyle = '#d9f8ff';
    for (const p of this.points) {
      ctx.beginPath(); ctx.arc(p.x, p.y, p.invMass === 0 ? 0.08 : 0.042, 0, Math.PI * 2); ctx.fill();
    }
  }

  diagnostics() {
    let kinetic = 0;
    for (const p of this.points) {
      if (p.invMass > 0) kinetic += 0.5 * (1 / p.invMass) * (p.vx ** 2 + p.vy ** 2);
    }
    const areaError = this.areaConstraint
      ? Math.abs(polygonArea(this.points, this.areaConstraint.indices) - this.areaConstraint.rest) / this.areaConstraint.rest
      : null;
    return {
      entities: this.points.length,
      contacts: this.constraints.length + (this.areaConstraint ? 1 : 0),
      energy: kinetic,
      residual: areaError ?? this.maxError,
      method: this.areaConstraint ? 'XPBD 距离 + 面积守恒' : 'XPBD 拉伸/剪切/弯曲网格',
    };
  }
}

export { polygonArea, solveDistance, solveArea };
