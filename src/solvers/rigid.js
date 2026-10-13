import { clamp, cross, dot, inverseRotate, length, normalize, rotate, Random } from '../core/math.js';
import { Solver } from '../core/simulation.js';

const TAU = Math.PI * 2;

function makeBody({ type, x, y, radius = 0.35, hx = 0.4, hy = 0.4, angle = 0, mass = 1, restitution = 0.16, friction = 0.52, color = '#ffbd59' }) {
  const isStatic = !Number.isFinite(mass) || mass <= 0;
  const inertia = type === 'circle'
    ? 0.5 * mass * radius * radius
    : (mass * ((hx * 2) ** 2 + (hy * 2) ** 2)) / 12;
  return {
    type, x, y, px: x, py: y, angle, vx: 0, vy: 0, omega: 0,
    radius, hx, hy, mass, invMass: isStatic ? 0 : 1 / mass,
    inertia, invInertia: isStatic ? 0 : 1 / inertia,
    restitution, friction, color, isStatic,
  };
}

function vertices(body) {
  const local = [[-body.hx, -body.hy], [body.hx, -body.hy], [body.hx, body.hy], [-body.hx, body.hy]];
  return local.map(([x, y]) => {
    const rotated = rotate(x, y, body.angle);
    return { x: body.x + rotated.x, y: body.y + rotated.y };
  });
}

function axesFor(body) {
  const x = rotate(1, 0, body.angle);
  const y = rotate(0, 1, body.angle);
  return [x, y];
}

function projectBox(body, axis) {
  const center = dot(body.x, body.y, axis.x, axis.y);
  const worldX = rotate(1, 0, body.angle);
  const worldY = rotate(0, 1, body.angle);
  const radius = body.hx * Math.abs(dot(worldX.x, worldX.y, axis.x, axis.y))
    + body.hy * Math.abs(dot(worldY.x, worldY.y, axis.x, axis.y));
  return { min: center - radius, max: center + radius };
}

function supportBox(body, nx, ny) {
  const local = inverseRotate(nx, ny, body.angle);
  const point = rotate(local.x >= 0 ? body.hx : -body.hx, local.y >= 0 ? body.hy : -body.hy, body.angle);
  return { x: body.x + point.x, y: body.y + point.y };
}

function circleCircle(a, b) {
  const dx = b.x - a.x;
  const dy = b.y - a.y;
  const distance = length(dx, dy);
  const target = a.radius + b.radius;
  if (distance >= target) return null;
  const normal = distance > 1e-7 ? { x: dx / distance, y: dy / distance } : { x: 1, y: 0 };
  return {
    a, b, nx: normal.x, ny: normal.y, penetration: target - distance,
    px: a.x + normal.x * (a.radius - (target - distance) * 0.5),
    py: a.y + normal.y * (a.radius - (target - distance) * 0.5),
  };
}

/** Returns a normal pointing from the circle to the box. */
function circleBox(circle, box) {
  const local = inverseRotate(circle.x - box.x, circle.y - box.y, box.angle);
  const closestX = clamp(local.x, -box.hx, box.hx);
  const closestY = clamp(local.y, -box.hy, box.hy);
  let dx = closestX - local.x;
  let dy = closestY - local.y;
  let distance = length(dx, dy);
  let penetration;

  if (distance < 1e-7) {
    const toX = box.hx - Math.abs(local.x);
    const toY = box.hy - Math.abs(local.y);
    if (toX < toY) {
      dx = local.x >= 0 ? 1 : -1;
      dy = 0;
      distance = 0;
      penetration = circle.radius + toX;
    } else {
      dx = 0;
      dy = local.y >= 0 ? 1 : -1;
      distance = 0;
      penetration = circle.radius + toY;
    }
  } else {
    penetration = circle.radius - distance;
  }
  if (penetration <= 0) return null;

  const inverseDistance = distance > 1e-7 ? 1 / distance : 1;
  const normal = rotate(dx * inverseDistance, dy * inverseDistance, box.angle);
  const pointLocal = { x: closestX, y: closestY };
  const point = rotate(pointLocal.x, pointLocal.y, box.angle);
  return { a: circle, b: box, nx: normal.x, ny: normal.y, penetration, px: box.x + point.x, py: box.y + point.y };
}

function boxBox(a, b) {
  let minOverlap = Infinity;
  let bestAxis = null;
  for (const axis of [...axesFor(a), ...axesFor(b)]) {
    const pa = projectBox(a, axis);
    const pb = projectBox(b, axis);
    const overlap = Math.min(pa.max, pb.max) - Math.max(pa.min, pb.min);
    if (overlap <= 0) return null;
    if (overlap < minOverlap) {
      minOverlap = overlap;
      const direction = dot(b.x - a.x, b.y - a.y, axis.x, axis.y) >= 0 ? 1 : -1;
      bestAxis = { x: axis.x * direction, y: axis.y * direction };
    }
  }
  const pa = supportBox(a, bestAxis.x, bestAxis.y);
  const pb = supportBox(b, -bestAxis.x, -bestAxis.y);
  return {
    a, b, nx: bestAxis.x, ny: bestAxis.y, penetration: minOverlap,
    px: (pa.x + pb.x) * 0.5, py: (pa.y + pb.y) * 0.5,
  };
}

function contactFor(a, b) {
  if (a.type === 'circle' && b.type === 'circle') return circleCircle(a, b);
  if (a.type === 'circle' && b.type === 'box') return circleBox(a, b);
  if (a.type === 'box' && b.type === 'circle') {
    const contact = circleBox(b, a);
    if (!contact) return null;
    return { ...contact, a, b, nx: -contact.nx, ny: -contact.ny };
  }
  return boxBox(a, b);
}

function velocityAt(body, rx, ry) {
  return { x: body.vx - body.omega * ry, y: body.vy + body.omega * rx };
}

function applyImpulse(body, ix, iy, rx, ry, sign) {
  if (body.invMass === 0) return;
  body.vx += sign * ix * body.invMass;
  body.vy += sign * iy * body.invMass;
  body.omega += sign * cross(rx, ry, ix, iy) * body.invInertia;
}

export class RigidSolver extends Solver {
  constructor({ width = 12, height = 7 } = {}) {
    super();
    this.width = width;
    this.height = height;
    this.gravity = -9.81;
    this.iterations = 8;
    this.bodies = [];
    this.contacts = [];
    this.random = new Random(0xdecafbad);
    this.grab = null;
    this.preset = 'stack';
    this.reset();
  }

  addCircle(x, y, radius = 0.28, speed = 0) {
    const body = makeBody({ type: 'circle', x, y, radius, mass: Math.PI * radius * radius, color: '#83d9ff', restitution: 0.32 });
    body.vx = speed;
    this.bodies.push(body);
    return body;
  }

  addBox(x, y, hx = 0.35, hy = 0.35, angle = 0, mass = null) {
    const body = makeBody({ type: 'box', x, y, hx, hy, angle, mass: mass ?? 4 * hx * hy, color: '#ffbd59' });
    this.bodies.push(body);
    return body;
  }

  addStaticBox(x, y, hx, hy, angle = 0) {
    const body = makeBody({ type: 'box', x, y, hx, hy, angle, mass: Infinity, color: '#547590', friction: 0.65 });
    this.bodies.push(body);
    return body;
  }

  clearDynamic() {
    this.bodies = this.bodies.filter((body) => body.isStatic);
  }

  makeBounds() {
    const t = 0.25;
    const bounds = [
      this.addStaticBox(0, -this.height / 2 - t * 0.5, this.width / 2 + t, t * 0.5),
      this.addStaticBox(-this.width / 2 - t * 0.5, 0, t * 0.5, this.height / 2 + t),
      this.addStaticBox(this.width / 2 + t * 0.5, 0, t * 0.5, this.height / 2 + t),
    ];
    for (const body of bounds) body.isBound = true;
  }

  reset(preset = this.preset, freshSeed = false) {
    this.preset = preset;
    const baseSeed = freshSeed ? (Math.random() * 2 ** 32) >>> 0 : (preset === 'collision' ? 0x2222 : 0x1111);
    this.random = new Random(baseSeed);
    this.bodies = [];
    this.makeBounds();
    if (preset === 'collision') {
      // Shuffle also re-poses the ramp slightly; the canonical reset keeps it fixed.
      const rampX = freshSeed ? this.random.range(-0.45, 0.45) : 0;
      const rampAngle = freshSeed ? this.random.range(-0.07, 0.07) : 0;
      this.addStaticBox(rampX, -2.25 + (freshSeed ? this.random.range(-0.12, 0.12) : 0), 3.3, 0.16, -0.14 + rampAngle);
      const a = this.addBox(-4.4, 1.2 + this.random.range(-0.5, 0.5), 0.52, 0.34, 0.25 + this.random.range(-0.3, 0.3), 1.2);
      a.vx = 8.5 + this.random.range(-1.5, 1.5);
      a.omega = -3 + this.random.range(-2, 2);
      const b = this.addCircle(2.7 + this.random.range(-0.6, 0.6), 1.0 + this.random.range(-0.4, 0.4), 0.43, -1.4);
      b.vx = -1.4 + this.random.range(-0.8, 0.8);
      for (let i = 0; i < 5; i += 1) this.addCircle(-0.8 + i * 0.46 + this.random.range(-0.12, 0.12), -0.8 + i * 0.1, 0.2 + (i % 2) * 0.04);
    } else {
      const rampX = freshSeed ? this.random.range(-0.55, 0.55) : 0;
      const rampY = freshSeed ? this.random.range(-0.14, 0.14) : 0;
      const rampAngle = freshSeed ? this.random.range(-0.09, 0.09) : 0;
      this.addStaticBox(1.6 + rampX, -2.2 + rampY, 2.0, 0.14, 0.13 + rampAngle);
      for (let row = 0; row < 6; row += 1) {
        for (let col = 0; col < 5 - Math.floor(row / 2); col += 1) {
          const jitterX = this.random.range(-0.09, 0.09);
          const jitterAngle = this.random.range(-0.14, 0.14);
          const body = row % 2 === 0
            ? this.addBox(-1.2 + col * 0.76 + jitterX, -1.7 + row * 0.72, 0.31, 0.31, (col - row) * 0.05 + jitterAngle)
            : this.addCircle(-1.0 + col * 0.78 + jitterX, -1.65 + row * 0.72, 0.28);
          body.vx = this.random.range(-0.15, 0.15);
        }
      }
    }
    this.contacts = [];
    this.grab = null;
  }

  /** Shuffle button: same preset, new randomized initial layout. */
  randomize() {
    this.reset(this.preset, true);
  }

  buildContacts() {
    const contacts = [];
    let maxPenetration = 0;
    for (let i = 0; i < this.bodies.length; i += 1) {
      const a = this.bodies[i];
      for (let j = i + 1; j < this.bodies.length; j += 1) {
        const b = this.bodies[j];
        if (a.invMass === 0 && b.invMass === 0) continue;
        const contact = contactFor(a, b);
        if (contact) {
          contacts.push(contact);
          maxPenetration = Math.max(maxPenetration, contact.penetration);
        }
      }
    }
    this.maxPenetration = maxPenetration;
    return contacts;
  }

  positionalCorrection(contact) {
    const { a, b, nx, ny, penetration } = contact;
    const denominator = a.invMass + b.invMass;
    if (denominator === 0) return;
    const correction = Math.max(penetration - 0.004, 0) * 0.72 / denominator;
    a.x -= nx * correction * a.invMass;
    a.y -= ny * correction * a.invMass;
    b.x += nx * correction * b.invMass;
    b.y += ny * correction * b.invMass;
  }

  solveVelocity(contact) {
    const { a, b, nx, ny, px, py } = contact;
    const rax = px - a.x;
    const ray = py - a.y;
    const rbx = px - b.x;
    const rby = py - b.y;
    const va = velocityAt(a, rax, ray);
    const vb = velocityAt(b, rbx, rby);
    const rvx = vb.x - va.x;
    const rvy = vb.y - va.y;
    const alongNormal = dot(rvx, rvy, nx, ny);
    if (alongNormal > 0) return;

    const raCrossN = cross(rax, ray, nx, ny);
    const rbCrossN = cross(rbx, rby, nx, ny);
    const denominator = a.invMass + b.invMass + raCrossN ** 2 * a.invInertia + rbCrossN ** 2 * b.invInertia;
    if (denominator < 1e-9) return;
    const restitution = Math.min(a.restitution, b.restitution);
    const impulse = -(1 + restitution) * alongNormal / denominator;
    const ix = impulse * nx;
    const iy = impulse * ny;
    applyImpulse(a, ix, iy, rax, ray, -1);
    applyImpulse(b, ix, iy, rbx, rby, 1);

    const postA = velocityAt(a, rax, ray);
    const postB = velocityAt(b, rbx, rby);
    const tangentRaw = { x: postB.x - postA.x - dot(postB.x - postA.x, postB.y - postA.y, nx, ny) * nx,
      y: postB.y - postA.y - dot(postB.x - postA.x, postB.y - postA.y, nx, ny) * ny };
    const tangent = normalize(tangentRaw.x, tangentRaw.y);
    const raCrossT = cross(rax, ray, tangent.x, tangent.y);
    const rbCrossT = cross(rbx, rby, tangent.x, tangent.y);
    const tangentDenominator = a.invMass + b.invMass + raCrossT ** 2 * a.invInertia + rbCrossT ** 2 * b.invInertia;
    if (tangentDenominator < 1e-9) return;
    const jt = -dot(postB.x - postA.x, postB.y - postA.y, tangent.x, tangent.y) / tangentDenominator;
    const friction = Math.sqrt(a.friction * b.friction);
    const capped = clamp(jt, -impulse * friction, impulse * friction);
    applyImpulse(a, capped * tangent.x, capped * tangent.y, rax, ray, -1);
    applyImpulse(b, capped * tangent.x, capped * tangent.y, rbx, rby, 1);
  }

  /**
   * A narrow-domain fallback guards against discrete tunnelling. It is kept
   * separate from normal SAT contacts so the ordinary contact solver remains
   * responsible for stacks and slopes, while a single missed contact cannot
   * eject an object from an interactive lesson.
   */
  safetyContain(body) {
    if (body.isStatic) return;
    const radius = body.type === 'circle' ? body.radius : Math.hypot(body.hx, body.hy);
    const left = -this.width / 2 + 0.15 + radius;
    const right = this.width / 2 - 0.15 - radius;
    const floor = -this.height / 2 + 0.15 + radius;
    if (body.x < left) { body.x = left; if (body.vx < 0) body.vx *= -body.restitution; }
    if (body.x > right) { body.x = right; if (body.vx > 0) body.vx *= -body.restitution; }
    if (body.y < floor) {
      body.y = floor;
      if (body.vy < 0) body.vy *= -body.restitution;
      body.vx *= 0.82;
    }
    body.vx = clamp(body.vx, -35, 35);
    body.vy = clamp(body.vy, -35, 35);
    body.omega = clamp(body.omega, -45, 45);
  }

  step(dt) {
    for (const body of this.bodies) {
      if (body.isStatic) continue;
      body.px = body.x;
      body.py = body.y;
      body.vy += this.gravity * dt;
      body.vx *= 0.9995;
      body.vy *= 0.9995;
      body.omega *= 0.998;
      body.x += body.vx * dt;
      body.y += body.vy * dt;
      body.angle = (body.angle + body.omega * dt) % TAU;
    }

    if (this.grab && !this.grab.kinematic) {
      const { body, x, y, offsetX, offsetY } = this.grab;
      const targetX = x + offsetX;
      const targetY = y + offsetY;
      body.vx += (targetX - body.x) * 65 * dt;
      body.vy += (targetY - body.y) * 65 * dt;
    }

    this.contacts = this.buildContacts();
    for (let pass = 0; pass < 2; pass += 1) {
      for (const contact of this.contacts) this.positionalCorrection(contact);
    }
    for (let pass = 0; pass < this.iterations; pass += 1) {
      for (const contact of this.contacts) this.solveVelocity(contact);
    }
    for (const body of this.bodies) this.safetyContain(body);
  }

  pointer(kind, x, y) {
    if (kind === 'down') {
      let best = null;
      let bestDistance = Infinity;
      for (const body of this.bodies) {
        if (body.isStatic) continue;
        const distance = Math.hypot(body.x - x, body.y - y);
        const capture = body.type === 'circle' ? body.radius * 1.4 : Math.max(body.hx, body.hy) * 1.5;
        if (distance < capture && distance < bestDistance) {
          best = body;
          bestDistance = distance;
        }
      }
      if (best) {
        this.grab = { body: best, x, y, offsetX: best.x - x, offsetY: best.y - y, kinematic: false };
        return;
      }
      // Otherwise allow dragging static obstacles (ramps), but never the walls.
      for (const body of this.bodies) {
        if (!body.isStatic || body.isBound) continue;
        const distance = Math.hypot(body.x - x, body.y - y);
        if (distance < Math.max(body.hx, body.hy) * 1.6) {
          this.grab = { body, x, y, offsetX: body.x - x, offsetY: body.y - y, kinematic: true };
          return;
        }
      }
    } else if (kind === 'move' && this.grab) {
      this.grab.x = x;
      this.grab.y = y;
      if (this.grab.kinematic) {
        const body = this.grab.body;
        body.x = clamp(x + this.grab.offsetX, -this.width / 2 + body.hx * 0.4, this.width / 2 - body.hx * 0.4);
        body.y = clamp(y + this.grab.offsetY, -this.height / 2 + 0.35, this.height / 2 - 0.35);
      }
    } else if (kind === 'up' || kind === 'leave') {
      this.grab = null;
    }
  }

  render(ctx) {
    for (const body of this.bodies) {
      ctx.save();
      ctx.translate(body.x, body.y);
      ctx.rotate(body.angle);
      ctx.fillStyle = body.isStatic ? '#1a354b' : body.color;
      ctx.strokeStyle = body.isStatic ? '#41637c' : '#11283b';
      ctx.lineWidth = 0.024;
      if (body.type === 'circle') {
        ctx.beginPath();
        ctx.arc(0, 0, body.radius, 0, TAU);
        ctx.fill();
        ctx.stroke();
        ctx.beginPath();
        ctx.moveTo(0, 0);
        ctx.lineTo(body.radius * 0.78, 0);
        ctx.stroke();
      } else {
        ctx.fillRect(-body.hx, -body.hy, body.hx * 2, body.hy * 2);
        ctx.strokeRect(-body.hx, -body.hy, body.hx * 2, body.hy * 2);
        if (body.isStatic && !body.isBound) {
          // Dashed rim signals the ramp is draggable.
          ctx.strokeStyle = 'rgba(186, 246, 255, 0.45)';
          ctx.setLineDash([0.09, 0.09]);
          ctx.strokeRect(-body.hx, -body.hy, body.hx * 2, body.hy * 2);
          ctx.setLineDash([]);
        }
      }
      ctx.restore();
    }
    if (this.grab) {
      ctx.fillStyle = '#ed5f64';
      ctx.beginPath();
      ctx.arc(this.grab.x, this.grab.y, 0.08, 0, TAU);
      ctx.fill();
    }
  }

  diagnostics() {
    let kinetic = 0;
    for (const body of this.bodies) {
      if (!body.isStatic) kinetic += 0.5 * body.mass * (body.vx ** 2 + body.vy ** 2) + 0.5 * body.inertia * body.omega ** 2;
    }
    return {
      entities: this.bodies.filter((body) => !body.isStatic).length,
      contacts: this.contacts.length,
      energy: kinetic,
      residual: this.maxPenetration ?? 0,
      method: '顺序冲量 · SAT/GJK-style 接触',
    };
  }
}

export { makeBody, contactFor, vertices };
