/** Small allocation-conscious 2D helpers shared by every solver. */
export const EPSILON = 1e-8;

export const clamp = (value, min, max) => Math.max(min, Math.min(max, value));
export const lerp = (a, b, t) => a + (b - a) * t;
export const dot = (ax, ay, bx, by) => ax * bx + ay * by;
export const cross = (ax, ay, bx, by) => ax * by - ay * bx;
export const length = (x, y) => Math.hypot(x, y);
export const lengthSquared = (x, y) => x * x + y * y;
export const perpendicular = (x, y) => ({ x: -y, y: x });

export function normalize(x, y) {
  const magnitude = Math.hypot(x, y);
  return magnitude > EPSILON ? { x: x / magnitude, y: y / magnitude } : { x: 0, y: 0 };
}

export function rotate(x, y, angle) {
  const c = Math.cos(angle);
  const s = Math.sin(angle);
  return { x: c * x - s * y, y: s * x + c * y };
}

export function inverseRotate(x, y, angle) {
  const c = Math.cos(angle);
  const s = Math.sin(angle);
  return { x: c * x + s * y, y: -s * x + c * y };
}

export function finite(value, fallback = 0) {
  return Number.isFinite(value) ? value : fallback;
}

/** Deterministic random generator makes reset/test behavior reproducible. */
export class Random {
  constructor(seed = 0x51f15e) {
    this.state = seed >>> 0;
  }

  next() {
    let value = this.state || 0x6d2b79f5;
    value ^= value << 13;
    value ^= value >>> 17;
    value ^= value << 5;
    this.state = value >>> 0;
    return this.state / 0x100000000;
  }

  range(min, max) {
    return min + (max - min) * this.next();
  }
}

/** Uniform grid with reusable buckets for O(n · k) local-neighbour searches. */
export class SpatialHash {
  constructor(cellSize) {
    this.cellSize = cellSize;
    this.cells = new Map();
  }

  clear() {
    this.cells.clear();
  }

  key(ix, iy) {
    return `${ix},${iy}`;
  }

  cell(x, y) {
    return [Math.floor(x / this.cellSize), Math.floor(y / this.cellSize)];
  }

  insert(index, x, y) {
    const [ix, iy] = this.cell(x, y);
    const key = this.key(ix, iy);
    const bucket = this.cells.get(key);
    if (bucket) bucket.push(index);
    else this.cells.set(key, [index]);
  }

  forEachNearby(x, y, radius, visit) {
    const [cx, cy] = this.cell(x, y);
    const reach = Math.ceil(radius / this.cellSize);
    for (let iy = cy - reach; iy <= cy + reach; iy += 1) {
      for (let ix = cx - reach; ix <= cx + reach; ix += 1) {
        const bucket = this.cells.get(this.key(ix, iy));
        if (!bucket) continue;
        for (const index of bucket) visit(index);
      }
    }
  }
}
