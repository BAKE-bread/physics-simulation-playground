import { clamp } from './math.js';

/**
 * Shared fixed-step clock. Rendering stays decoupled from numerical time so a
 * dropped animation frame never silently changes the physical timestep.
 */
export class FixedStepClock {
  constructor({ dt = 1 / 120, maxFrameTime = 0.1, maxSubsteps = 10 } = {}) {
    this.dt = dt;
    this.maxFrameTime = maxFrameTime;
    this.maxSubsteps = maxSubsteps;
    this.accumulator = 0;
  }

  reset() {
    this.accumulator = 0;
  }

  advance(frameSeconds, step) {
    this.accumulator += clamp(frameSeconds, 0, this.maxFrameTime);
    let steps = 0;
    while (this.accumulator >= this.dt && steps < this.maxSubsteps) {
      step(this.dt);
      this.accumulator -= this.dt;
      steps += 1;
    }
    if (steps === this.maxSubsteps) this.accumulator = 0;
    return steps;
  }
}

/** Interface convention used by the UI and by headless numeric tests. */
export class Solver {
  reset() {}
  step() {}
  render() {}
  pointer() {}
  diagnostics() { return {}; }
}
