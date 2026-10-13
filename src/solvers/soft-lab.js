import { Solver } from '../core/simulation.js';
import { SoftSolver } from './soft.js';
import { FemSolver } from './fem.js';

/** Keeps XPBD and energy-based FEM selectable through the same soft-body tab. */
export class SoftLabSolver extends Solver {
  constructor(options) {
    super();
    this.xpbd = new SoftSolver(options);
    this.fem = new FemSolver(options);
    this.mode = 'cloth';
  }

  get active() { return this.mode === 'fem' ? this.fem : this.xpbd; }
  get gravity() { return this.xpbd.gravity; }
  set gravity(value) { this.xpbd.gravity = value; this.fem.gravity = value; }
  get iterations() { return this.xpbd.iterations; }
  set iterations(value) { this.xpbd.iterations = value; }
  get substeps() { return this.xpbd.substeps; }
  set substeps(value) { this.xpbd.substeps = value; this.fem.substeps = Math.max(5, value * 2); }
  reset(mode = this.mode) { this.mode = mode; this.active.reset(mode); }
  step(dt) { this.active.step(dt); }
  render(ctx) { this.active.render(ctx); }
  pointer(kind, x, y) { this.active.pointer(kind, x, y); }
  diagnostics() { return this.active.diagnostics(); }
}
