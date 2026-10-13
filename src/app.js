import { FixedStepClock } from './core/simulation.js';
import { RigidSolver } from './solvers/rigid.js';
import { SoftLabSolver } from './solvers/soft-lab.js';
import { FluidSolver } from './solvers/fluid.js';
import { PdeLabSolver } from './solvers/fields.js';

const WORLD = { width: 12, height: 7 };
const canvas = document.querySelector('#simCanvas');
const ctx = canvas.getContext('2d');
const scenarioSelect = document.querySelector('#scenarioSelect');
const scenarioHint = document.querySelector('#scenarioHint');
const scenarioTitle = document.querySelector('#scenarioTitle');
const familyLabel = document.querySelector('#familyLabel');
const interactionHint = document.querySelector('#interactionHint');
const legendText = document.querySelector('#legendText');
const legendText2 = document.querySelector('#legendText2');
const playButton = document.querySelector('#playButton');
const stepButton = document.querySelector('#stepButton');
const resetButton = document.querySelector('#resetButton');
const shuffleButton = document.querySelector('#shuffleButton');
const speedInput = document.querySelector('#speedInput');
const gravityInput = document.querySelector('#gravityInput');
const qualityInput = document.querySelector('#qualityInput');
const adaptiveInput = document.querySelector('#adaptiveInput');

const solvers = {
  rigid: new RigidSolver(WORLD),
  soft: new SoftLabSolver(WORLD),
  fluid: new FluidSolver(WORLD),
  field: new PdeLabSolver(WORLD),
};

const scenarios = {
  rigid: [
    { id: 'stack', label: '堆叠与摩擦接触', tag: 'RIGID BODY', hint: '多体堆叠会持续生成接触流形；拖动物体观察冲量、转动与静摩擦。', interaction: '拖动物体或斜坡（虚线框）；“打乱”换个初始布局。', legend: '动态刚体', legend2: '可拖斜坡', shuffle: true },
    { id: 'collision', label: '斜坡冲量碰撞', tag: 'RIGID BODY', hint: '旋转矩形与圆盘在斜坡上发生非中心碰撞，突出角动量耦合。', interaction: '拖动物体或斜坡；“打乱”改变入射状态。', legend: '冲量刚体', legend2: '可拖斜坡', shuffle: true },
  ],
  soft: [
    { id: 'cloth', label: '悬挂布料与障碍', tag: 'COMPLIANT XPBD', hint: '拉伸、剪切和弯曲约束组成网格；颜色随局部应变增强。', interaction: '拖住任一节点，或拖动障碍球（虚线圆）改变缠绕。', legend: 'XPBD 网格', legend2: '可拖障碍球' },
    { id: 'balloon', label: '保面积软体水球', tag: 'COMPLIANT XPBD', hint: '闭环边/弯曲约束与全局面积约束共同形成可压缩程度可控的软体。', interaction: '拖动边界节点，或推动障碍球挤压水球。', legend: '面积约束软体', legend2: '可拖障碍球' },
    { id: 'fem', label: '共旋 FEM 悬臂梁', tag: 'TRIANGULAR FEM', hint: '三角形单元通过共旋弹性能量传递内力；可与 XPBD 的位置投影直接对照。', interaction: '拖动网格节点，观察应变能与局部应变着色。', legend: '三角 FEM 网格', legend2: '应变着色' },
  ],
  fluid: [
    { id: 'smoke', label: '不可压缩烟流', tag: 'EULERIAN CFD', hint: '规则网格采用半拉格朗日平流和压力泊松投影；分辨率随质量档位提升。', interaction: '在画布中拖动，注入染料与动量。', legend: '染料浓度场', legend2: '速度矢量' },
    { id: 'lbm', label: '绕流圆柱 · 卡门涡街', tag: 'LATTICE BOLTZMANN', hint: 'D2Q9 格子 Boltzmann 在中等雷诺数下自发形成涡街；涡量着色让脱落清晰可见。', interaction: '拖动圆柱改变绕流；在流场中按住滑动可轻柔搅动尾迹。', legend: '涡量（红顺/蓝逆）', legend2: '可拖圆柱' },
    { id: 'sph', label: 'SPH 粒子水槽', tag: 'LAGRANGIAN SPH', hint: '空间哈希邻域、核密度和双密度压力松弛在实时性与可塑性间折中。', interaction: '搅动水粒子；“打乱”整体重排水块初态。', legend: 'SPH 粒子', legend2: '密度着色', shuffle: true },
    { id: 'shallow', label: '浅水波与扰动', tag: 'SHALLOW WATER', hint: '线性浅水方程在规则网格上按 CFL 条件自动细分步长，适合长波过程。', interaction: '在水面拖动，制造波高与流速扰动。', legend: '水位高度场', legend2: '流速矢量' },
  ],
  field: [
    { id: 'heat', label: '含孔板热传导 FEM', tag: 'HEAT PDE · FEM', hint: 'P1 三角单元组装刚度矩阵，后向欧拉 + Jacobi 推进瞬态热方程；孔洞为绝热自然边界。', interaction: '拖动涂抹热源；按住圆孔内部可整体移动孔位。', legend: '温度场', legend2: '可移孔洞' },
    { id: 'wave', label: '薄膜波动方程', tag: 'WAVE PDE · FDM', hint: '蛙跳格式显式推进双曲波动方程，CFL 数显示在诊断栏；斜率光照突出波形。', interaction: '点击或拖动，拨动膜面激起涟漪。', legend: '位移高度场', legend2: '斜率光照' },
    { id: 'grayscott', label: 'Gray–Scott 反应扩散', tag: 'TURING PATTERNS', hint: '两个耦合 PDE 自发长出珊瑚状图灵斑图，是抛物型系统自组织的经典范例。', interaction: '涂抹画布播撒反应物；“打乱”重新播种。', legend: '反应物 V 浓度', legend2: '自组织斑图', shuffle: true },
    { id: 'flock', label: 'Boids 群体运动', tag: 'AGENT FIELD', hint: '分离、对齐、聚集三条局部规则涌现出群体有序；极化度残差显示有序程度。', interaction: '按住画布放置引导点；“打乱”重组群体。', legend: '个体（按速度着色）', legend2: '引导点', shuffle: true },
  ],
};

const familyOrder = ['rigid', 'soft', 'fluid', 'field'];

let family = 'rigid';
let active = solvers[family];
let running = true;
let targetQuality = Number(qualityInput.value);
let appliedQuality = targetQuality;
let viewport = { width: 1, height: 1, scale: 1, dpr: 1 };
let fps = 60;
let fpsFrames = 0;
let fpsElapsed = 0;
let lastFrame = performance.now();
let lastSteps = 0;
const clock = new FixedStepClock({ dt: 1 / 120, maxSubsteps: 12 });

function currentScenario() {
  return scenarios[family].find((scenario) => scenario.id === scenarioSelect.value) ?? scenarios[family][0];
}

function syncScenarioUI() {
  const scenario = currentScenario();
  scenarioTitle.textContent = scenario.label;
  scenarioHint.textContent = scenario.hint;
  familyLabel.textContent = scenario.tag;
  interactionHint.textContent = scenario.interaction;
  legendText.textContent = scenario.legend;
  legendText2.textContent = scenario.legend2;
  // Hide the shuffle button where it would just duplicate 重置.
  shuffleButton.hidden = !scenario.shuffle;
}

function populateScenarios(nextFamily, preferredId) {
  scenarioSelect.textContent = '';
  for (const scenario of scenarios[nextFamily]) {
    const option = document.createElement('option');
    option.value = scenario.id;
    option.textContent = scenario.label;
    scenarioSelect.append(option);
  }
  scenarioSelect.value = scenarios[nextFamily].some((scenario) => scenario.id === preferredId) ? preferredId : scenarios[nextFamily][0].id;
}

function applyGravity() {
  const gravity = -Number(gravityInput.value);
  solvers.rigid.gravity = gravity;
  solvers.soft.gravity = gravity;
  solvers.fluid.water.gravity = gravity;
  solvers.fluid.shallow.gravity = -gravity;
  document.querySelector('#gravityValue').textContent = Math.abs(gravity).toFixed(1);
}

function qualityText(level, auto) {
  const names = ['', '低', '中', '高'];
  return auto && level !== targetQuality ? `${names[level]}（自适应）` : names[level];
}

/**
 * Iteration budgets adapt freely, but grid/lattice resolutions only change
 * on explicit quality edits or resets — re-gridding silently mid-run would
 * wipe the field the user is watching.
 */
function applyQuality(level = appliedQuality, withResolution = false) {
  appliedQuality = level;
  solvers.rigid.iterations = 3 + level * 4;
  solvers.soft.iterations = 4 + level * 5;
  solvers.soft.substeps = level >= 3 ? 3 : 2;
  solvers.fluid.smoke.iterations = 10 + level * 12;
  solvers.fluid.water.substeps = 1 + level;
  clock.maxSubsteps = 6 + level * 3;
  if (withResolution) {
    solvers.fluid.smoke.setResolution({ 1: 72, 2: 110, 3: 144 }[level]);
    solvers.fluid.shallow.setResolution({ 1: 88, 2: 124, 3: 168 }[level]);
    solvers.fluid.lbm.setQuality(level);
    solvers.field.setQuality(level);
  }
  document.querySelector('#qualityValue').textContent = qualityText(level, adaptiveInput.checked);
}

function resetActive() {
  active.reset(scenarioSelect.value);
  clock.reset();
  applyGravity();
  applyQuality(appliedQuality, true);
}

function selectFamily(nextFamily) {
  family = nextFamily;
  active = solvers[family];
  populateScenarios(family);
  for (const button of document.querySelectorAll('.family-tab')) {
    const selected = button.dataset.family === family;
    button.classList.toggle('is-active', selected);
    button.setAttribute('aria-selected', String(selected));
  }
  syncScenarioUI();
  resetActive();
}

function setRunning(next) {
  running = next;
  playButton.textContent = running ? '暂停' : '继续';
  playButton.setAttribute('aria-pressed', String(!running));
  document.querySelector('#runStatus').textContent = running ? '运行中' : '已暂停';
  document.querySelector('.status-dot').style.background = running ? '#66d699' : '#ffbd59';
  clock.reset();
}

function stepOnce() {
  if (running) setRunning(false);
  active.step(clock.dt);
  lastSteps = 1;
  render();
  updateDiagnostics();
}

function resizeCanvas() {
  const rect = canvas.getBoundingClientRect();
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  const width = Math.max(1, Math.round(rect.width * dpr));
  const height = Math.max(1, Math.round(rect.height * dpr));
  if (canvas.width !== width || canvas.height !== height) {
    canvas.width = width;
    canvas.height = height;
  }
  viewport = { width: rect.width, height: rect.height, dpr, scale: Math.min(rect.width / WORLD.width, rect.height / WORLD.height) };
}

function drawWorldFrame() {
  ctx.save();
  ctx.strokeStyle = 'rgba(109, 177, 200, 0.11)';
  ctx.lineWidth = 0.008;
  for (let x = -5; x <= 5; x += 1) {
    ctx.beginPath(); ctx.moveTo(x, -WORLD.height / 2); ctx.lineTo(x, WORLD.height / 2); ctx.stroke();
  }
  for (let y = -3; y <= 3; y += 1) {
    ctx.beginPath(); ctx.moveTo(-WORLD.width / 2, y); ctx.lineTo(WORLD.width / 2, y); ctx.stroke();
  }
  ctx.strokeStyle = 'rgba(117, 202, 225, 0.42)';
  ctx.lineWidth = 0.024;
  ctx.strokeRect(-WORLD.width / 2, -WORLD.height / 2, WORLD.width, WORLD.height);
  ctx.restore();
}

function render() {
  const { width, height, dpr, scale } = viewport;
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, width, height);
  const backdrop = ctx.createLinearGradient(0, 0, 0, height);
  backdrop.addColorStop(0, '#0b202c');
  backdrop.addColorStop(1, '#07131d');
  ctx.fillStyle = backdrop;
  ctx.fillRect(0, 0, width, height);
  ctx.setTransform(dpr * scale, 0, 0, -dpr * scale, dpr * width / 2, dpr * height / 2);
  drawWorldFrame();
  active.render(ctx);
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
}

function formatValue(value, places = 2) {
  if (!Number.isFinite(value)) return '—';
  if (Math.abs(value) >= 10000 || (Math.abs(value) > 0 && Math.abs(value) < 0.001)) return value.toExponential(2);
  return value.toFixed(places);
}

function updateDiagnostics() {
  const diagnostic = active.diagnostics();
  document.querySelector('#entitiesMetric').textContent = diagnostic.entities;
  document.querySelector('#contactsMetric').textContent = diagnostic.contacts;
  document.querySelector('#energyMetric').textContent = formatValue(diagnostic.energy);
  document.querySelector('#residualMetric').textContent = formatValue(diagnostic.residual, 4);
  document.querySelector('#stepsMetric').textContent = lastSteps;
  document.querySelector('#fpsMetric').textContent = `${Math.round(fps)}`;
  document.querySelector('#fpsReadout').textContent = `${Math.round(fps)} FPS`;
  document.querySelector('#methodMetric').textContent = diagnostic.method;
}

function adaptQuality() {
  if (!adaptiveInput.checked) {
    if (appliedQuality !== targetQuality) applyQuality(targetQuality);
    return;
  }
  if (fps < 42 && appliedQuality > 1) applyQuality(appliedQuality - 1);
  if (fps > 58 && appliedQuality < targetQuality) applyQuality(appliedQuality + 1);
}

function tick(now) {
  const elapsed = Math.min((now - lastFrame) / 1000, 0.1);
  lastFrame = now;
  if (running) lastSteps = clock.advance(elapsed * Number(speedInput.value), (dt) => active.step(dt));
  else lastSteps = 0;
  fpsFrames += 1;
  fpsElapsed += elapsed;
  if (fpsElapsed >= 0.5) {
    fps = fpsFrames / fpsElapsed;
    fpsFrames = 0;
    fpsElapsed = 0;
    adaptQuality();
  }
  render();
  updateDiagnostics();
  requestAnimationFrame(tick);
}

function toWorld(event) {
  const rect = canvas.getBoundingClientRect();
  return {
    x: (event.clientX - rect.left - rect.width / 2) / viewport.scale,
    y: -(event.clientY - rect.top - rect.height / 2) / viewport.scale,
  };
}

for (const tab of document.querySelectorAll('.family-tab')) {
  tab.addEventListener('click', () => selectFamily(tab.dataset.family));
}
scenarioSelect.addEventListener('change', () => { syncScenarioUI(); resetActive(); });
playButton.addEventListener('click', () => setRunning(!running));
stepButton.addEventListener('click', stepOnce);
resetButton.addEventListener('click', resetActive);
shuffleButton.addEventListener('click', () => {
  if (typeof active.randomize === 'function') active.randomize();
  else active.reset(scenarioSelect.value);
  clock.reset();
});
speedInput.addEventListener('input', () => { document.querySelector('#speedValue').textContent = `${Number(speedInput.value).toFixed(2)}×`; });
gravityInput.addEventListener('input', applyGravity);
qualityInput.addEventListener('input', () => { targetQuality = Number(qualityInput.value); applyQuality(targetQuality, true); });
adaptiveInput.addEventListener('change', () => applyQuality(adaptiveInput.checked ? appliedQuality : targetQuality));

// Keyboard shortcuts: Space 暂停/继续, S 单步, R 重置, 1-4 切换家族.
window.addEventListener('keydown', (event) => {
  if (event.target instanceof HTMLInputElement || event.target instanceof HTMLSelectElement) return;
  if (event.code === 'Space') { event.preventDefault(); setRunning(!running); }
  else if (event.key === 's' || event.key === 'S') stepOnce();
  else if (event.key === 'r' || event.key === 'R') resetActive();
  else if (/^[1-4]$/.test(event.key)) selectFamily(familyOrder[Number(event.key) - 1]);
});

canvas.addEventListener('pointerdown', (event) => {
  canvas.setPointerCapture(event.pointerId);
  const point = toWorld(event);
  active.pointer('down', point.x, point.y);
});
canvas.addEventListener('pointermove', (event) => {
  if (!(event.buttons & 1)) return;
  const point = toWorld(event);
  active.pointer('move', point.x, point.y);
});
canvas.addEventListener('pointerup', (event) => {
  const point = toWorld(event);
  active.pointer('up', point.x, point.y);
});
canvas.addEventListener('pointercancel', () => active.pointer('leave', 0, 0));
canvas.addEventListener('lostpointercapture', () => active.pointer('leave', 0, 0));
new ResizeObserver(resizeCanvas).observe(canvas);

populateScenarios(family, 'stack');
syncScenarioUI();
applyGravity();
applyQuality(appliedQuality, true);
resizeCanvas();
requestAnimationFrame(tick);
