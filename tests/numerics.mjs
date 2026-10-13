import assert from 'node:assert/strict';
import { RigidSolver } from '../src/solvers/rigid.js';
import { SoftSolver } from '../src/solvers/soft.js';
import { FemSolver } from '../src/solvers/fem.js';
import { FluidSolver } from '../src/solvers/fluid.js';
import { PdeLabSolver } from '../src/solvers/fields.js';
import { FlockSolver } from '../src/solvers/flock.js';

const finite = (value, label) => assert.ok(Number.isFinite(value), `${label} must be finite`);

function testRigid() {
  const solver = new RigidSolver();
  for (let step = 0; step < 720; step += 1) solver.step(1 / 120);
  const dynamic = solver.bodies.filter((body) => !body.isStatic);
  assert.ok(dynamic.length >= 10, 'rigid preset should retain a nontrivial body set');
  for (const body of dynamic) {
    finite(body.x, 'rigid x'); finite(body.y, 'rigid y'); finite(body.vx, 'rigid vx'); finite(body.omega, 'rigid omega');
    assert.ok(body.y >= -3.9, 'rigid body should remain near the bounded domain');
  }
  finite(solver.diagnostics().energy, 'rigid kinetic energy');
  // Shuffled layouts must be just as stable as the canonical one.
  solver.randomize();
  for (let step = 0; step < 360; step += 1) solver.step(1 / 120);
  for (const body of solver.bodies) {
    if (body.isStatic) continue;
    finite(body.x, 'rigid shuffled x'); finite(body.y, 'rigid shuffled y');
    assert.ok(body.y >= -3.9, 'shuffled bodies should stay in the domain');
  }
  // Dragging the static ramp is kinematic and must not corrupt the stack.
  const ramp = solver.bodies.find((body) => body.isStatic && !body.isBound);
  assert.ok(ramp, 'a draggable ramp should exist');
  solver.pointer('down', ramp.x, ramp.y);
  solver.pointer('move', ramp.x - 0.8, ramp.y + 0.3);
  solver.pointer('up', 0, 0);
  for (let step = 0; step < 120; step += 1) solver.step(1 / 120);
  finite(solver.diagnostics().energy, 'rigid energy after ramp drag');
}

function testSoft() {
  const solver = new SoftSolver();
  for (let step = 0; step < 480; step += 1) solver.step(1 / 120);
  for (const point of solver.points) {
    finite(point.x, 'soft x'); finite(point.y, 'soft y'); finite(point.vx, 'soft vx');
  }
  assert.ok(solver.diagnostics().residual < 0.25, 'cloth constraints should remain bounded');
  solver.reset('balloon');
  for (let step = 0; step < 480; step += 1) solver.step(1 / 120);
  assert.ok(solver.diagnostics().residual < 0.08, 'area-preserving balloon should retain its area');

  const fem = new FemSolver();
  for (let step = 0; step < 720; step += 1) fem.step(1 / 120);
  for (const node of fem.nodes) {
    finite(node.x, 'FEM x'); finite(node.y, 'FEM y'); finite(node.vx, 'FEM vx');
  }
  assert.ok(fem.diagnostics().residual < 0.5, 'corotated FEM should retain a bounded area distortion at default load');
}

function testFlock() {
  const flock = new FlockSolver({ width: 12, height: 7 });
  for (let step = 0; step < 360; step += 1) flock.step(1 / 120);
  for (const boid of flock.boids) {
    finite(boid.x, 'flock x'); finite(boid.y, 'flock y'); finite(boid.vx, 'flock vx');
    assert.ok(Math.abs(boid.x) <= 6 && Math.abs(boid.y) <= 3.5, 'boids should stay inside the world');
  }
  const diagnostic = flock.diagnostics();
  finite(diagnostic.energy, 'flock energy');
  assert.ok(diagnostic.residual >= 0 && diagnostic.residual <= 1, 'flock polarisation residual should lie in [0, 1]');
}

function testFluids() {
  const solver = new FluidSolver();
  const plan = {
    smoke: { steps: 180, limit: 2 },
    sph: { steps: 180, limit: 0.5 },
    shallow: { steps: 180, limit: 1 },
    lbm: { steps: 90, limit: 0.2 },
  };
  for (const [mode, { steps, limit }] of Object.entries(plan)) {
    solver.reset(mode);
    for (let step = 0; step < steps; step += 1) solver.step(1 / 120);
    const diagnostic = solver.diagnostics();
    finite(diagnostic.energy, `${mode} energy`);
    finite(diagnostic.residual, `${mode} residual`);
    assert.ok(diagnostic.residual < limit, `${mode} residual should be bounded`);
  }
  // Grid re-resolution must not destabilise the projection solver.
  solver.reset('smoke');
  solver.smoke.setResolution(144);
  for (let step = 0; step < 90; step += 1) solver.step(1 / 120);
  finite(solver.diagnostics().residual, 'smoke high-res residual');
  // The pointer brush must perturb the LBM wake without destabilising it.
  solver.reset('lbm');
  solver.pointer('down', 1.5, 0.5);
  for (let step = 0; step < 60; step += 1) {
    solver.pointer('move', 1.5 + Math.sin(step * 0.1) * 0.8, 0.5);
    solver.step(1 / 120);
  }
  solver.pointer('up', 0, 0);
  for (let step = 0; step < 60; step += 1) solver.step(1 / 120);
  const brushed = solver.diagnostics();
  finite(brushed.residual, 'lbm brushed residual');
  assert.ok(brushed.residual < 0.3, 'lbm brush should keep the flow stable');
  // SPH shuffle rebuilds the block at a fresh position without instability.
  solver.reset('sph');
  solver.water.randomize();
  for (let step = 0; step < 180; step += 1) solver.step(1 / 120);
  for (const p of solver.water.particles) {
    finite(p.x, 'sph shuffled x'); finite(p.y, 'sph shuffled y');
  }
  assert.ok(solver.diagnostics().residual < 0.5, 'sph shuffled density should be bounded');
}

function testFields() {
  const lab = new PdeLabSolver({ width: 12, height: 7 });

  lab.reset('heat');
  for (let step = 0; step < 240; step += 1) lab.step(1 / 120);
  for (const t of lab.heat.temperature) {
    finite(t, 'heat temperature');
    assert.ok(t >= 0 && t <= 4.6, 'plate temperature should stay within the physical clamp');
  }
  finite(lab.diagnostics().energy, 'heat energy');
  // Dragging the hole re-meshes; the carried-over field must stay stable.
  lab.pointer('down', lab.heat.hole.x, lab.heat.hole.y);
  lab.pointer('move', lab.heat.hole.x - 1.2, lab.heat.hole.y + 0.5);
  lab.pointer('up', 0, 0);
  for (let step = 0; step < 120; step += 1) lab.step(1 / 120);
  for (const t of lab.heat.temperature) finite(t, 'heat temperature after hole drag');
  finite(lab.diagnostics().residual, 'heat residual after hole drag');

  lab.reset('wave');
  for (let step = 0; step < 360; step += 1) lab.step(1 / 120);
  const wave = lab.diagnostics();
  finite(wave.energy, 'wave energy');
  assert.ok(wave.residual < 5, 'membrane displacement should remain bounded under the explicit scheme');

  lab.reset('grayscott');
  for (let step = 0; step < 40; step += 1) lab.step(1 / 120);
  for (const value of lab.grayscott.v) {
    finite(value, 'gray-scott v');
    assert.ok(value >= -1e-4 && value <= 1.2, 'reagent concentration should remain bounded');
  }

  lab.reset('flock');
  for (let step = 0; step < 240; step += 1) lab.step(1 / 120);
  for (const boid of lab.flock.boids) finite(boid.x, 'field-family flock x');
}

testRigid();
testSoft();
testFlock();
testFluids();
testFields();
console.log('Numerical smoke tests passed.');
