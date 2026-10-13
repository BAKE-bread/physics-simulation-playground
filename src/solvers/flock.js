import { clamp, SpatialHash, Random } from '../core/math.js';
import { Solver } from '../core/simulation.js';

/**
 * Reynolds boids: separation / alignment / cohesion steering with a uniform
 * spatial hash for neighbour lookups. Emergent flocking is a nice contrast
 * to the continuum models, and the pointer acts as a moving attractor.
 */
export class FlockSolver extends Solver {
  constructor({ width = 12, height = 7 } = {}) {
    super();
    this.width = width;
    this.height = height;
    this.perception = 0.62;
    this.hash = new SpatialHash(this.perception);
    this.maxSpeed = 3.4;
    this.maxForce = 9.5;
    this.random = new Random(0xb01d5);
    this.attractor = null;
    this.level = 0;
    this.setQuality(2);
  }

  setQuality(level) {
    if (level === this.level) return;
    this.level = level;
    this.count = { 1: 110, 2: 160, 3: 220 }[level] ?? 160;
    this.reset();
  }

  reset() {
    this.boids = [];
    for (let i = 0; i < this.count; i += 1) {
      const angle = this.random.range(0, Math.PI * 2);
      this.boids.push({
        x: this.random.range(-this.width / 2 + 1, this.width / 2 - 1),
        y: this.random.range(-this.height / 2 + 0.8, this.height / 2 - 0.8),
        vx: Math.cos(angle) * this.random.range(1.2, 2.6),
        vy: Math.sin(angle) * this.random.range(1.2, 2.6),
      });
    }
  }

  /** Shuffle button: relaunch with a fresh deterministic-seed sequence. */
  randomize() {
    this.random = new Random((Math.random() * 2 ** 32) >>> 0);
    this.reset();
  }

  steerTowards(boid, tx, ty, weight, fx, fy) {
    const dx = tx - boid.x;
    const dy = ty - boid.y;
    const d = Math.hypot(dx, dy);
    if (d < 1e-6) return [fx, fy];
    return [fx + (dx / d) * weight, fy + (dy / d) * weight];
  }

  step(dt) {
    const { boids, perception, hash, maxSpeed, maxForce } = this;
    hash.clear();
    boids.forEach((boid, index) => hash.insert(index, boid.x, boid.y));
    const margin = 1.0;
    for (let i = 0; i < boids.length; i += 1) {
      const boid = boids[i];
      let sepX = 0; let sepY = 0; let sepCount = 0;
      let aliX = 0; let aliY = 0; let count = 0;
      let cohX = 0; let cohY = 0;
      hash.forEachNearby(boid.x, boid.y, perception, (j) => {
        if (j === i) return;
        const other = boids[j];
        const dx = boid.x - other.x;
        const dy = boid.y - other.y;
        const d = Math.hypot(dx, dy);
        if (d >= perception || d < 1e-6) return;
        count += 1;
        aliX += other.vx; aliY += other.vy;
        cohX += other.x; cohY += other.y;
        if (d < perception * 0.55) {
          sepX += dx / (d * d);
          sepY += dy / (d * d);
          sepCount += 1;
        }
      });
      let fx = 0; let fy = 0;
      if (sepCount > 0) {
        const d = Math.hypot(sepX, sepY) || 1;
        fx += (sepX / d) * 1.7;
        fy += (sepY / d) * 1.7;
      }
      if (count > 0) {
        aliX /= count; aliY /= count;
        const d = Math.hypot(aliX, aliY) || 1;
        fx += (aliX / d) * 1.0;
        fy += (aliY / d) * 1.0;
        cohX /= count; cohY /= count;
        [fx, fy] = this.steerTowards(boid, cohX, cohY, 0.9, fx, fy);
      }
      // Soft walls steer the flock back into view instead of hard clamping.
      if (boid.x < -this.width / 2 + margin) fx += 2.4 * (1 - (boid.x + this.width / 2) / margin + 1);
      if (boid.x > this.width / 2 - margin) fx -= 2.4 * (1 - (this.width / 2 - boid.x) / margin + 1);
      if (boid.y < -this.height / 2 + margin) fy += 2.4 * (1 - (boid.y + this.height / 2) / margin + 1);
      if (boid.y > this.height / 2 - margin) fy -= 2.4 * (1 - (this.height / 2 - boid.y) / margin + 1);
      if (this.attractor) {
        [fx, fy] = this.steerTowards(boid, this.attractor.x, this.attractor.y, 1.9, fx, fy);
      }
      const forceMag = Math.hypot(fx, fy);
      if (forceMag > maxForce) { fx = (fx / forceMag) * maxForce; fy = (fy / forceMag) * maxForce; }
      boid.vx += fx * dt * 3;
      boid.vy += fy * dt * 3;
      const speed = Math.hypot(boid.vx, boid.vy);
      if (speed > maxSpeed) { boid.vx = (boid.vx / speed) * maxSpeed; boid.vy = (boid.vy / speed) * maxSpeed; }
      if (speed < 1.1 && speed > 1e-6) { boid.vx = (boid.vx / speed) * 1.1; boid.vy = (boid.vy / speed) * 1.1; }
      boid.x = clamp(boid.x + boid.vx * dt, -this.width / 2 + 0.12, this.width / 2 - 0.12);
      boid.y = clamp(boid.y + boid.vy * dt, -this.height / 2 + 0.12, this.height / 2 - 0.12);
    }
  }

  pointer(kind, x, y) {
    if (kind === 'down' || kind === 'move') this.attractor = { x, y };
    else this.attractor = null;
  }

  render(ctx) {
    for (const boid of this.boids) {
      const speed = Math.hypot(boid.vx, boid.vy);
      const t = clamp(speed / this.maxSpeed, 0, 1);
      const angle = Math.atan2(boid.vy, boid.vx);
      ctx.save();
      ctx.translate(boid.x, boid.y);
      ctx.rotate(angle);
      ctx.fillStyle = `rgba(${Math.round(90 + 130 * t)}, ${Math.round(200 + 40 * t)}, ${Math.round(240 - 60 * t)}, 0.9)`;
      ctx.beginPath();
      ctx.moveTo(0.15, 0);
      ctx.lineTo(-0.11, 0.075);
      ctx.lineTo(-0.06, 0);
      ctx.lineTo(-0.11, -0.075);
      ctx.closePath();
      ctx.fill();
      ctx.restore();
    }
    if (this.attractor) {
      ctx.strokeStyle = 'rgba(255, 189, 89, 0.5)';
      ctx.lineWidth = 0.015;
      ctx.beginPath();
      ctx.arc(this.attractor.x, this.attractor.y, 0.22, 0, Math.PI * 2);
      ctx.stroke();
    }
  }

  diagnostics() {
    let energy = 0;
    let orderX = 0;
    let orderY = 0;
    for (const boid of this.boids) {
      energy += 0.5 * (boid.vx ** 2 + boid.vy ** 2);
      const speed = Math.hypot(boid.vx, boid.vy) || 1;
      orderX += boid.vx / speed;
      orderY += boid.vy / speed;
    }
    const n = Math.max(this.boids.length, 1);
    const polarisation = Math.hypot(orderX, orderY) / n;
    return {
      entities: this.boids.length,
      contacts: this.hash.cells.size,
      energy: energy / n,
      residual: 1 - polarisation,
      method: `Boids 三规则转向 · 空间哈希 · ${this.boids.length} 个体`,
    };
  }
}
