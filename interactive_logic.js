/**
 * Interactive Algorithm Visualization Engine
 * PathOpt: Frontier Reduction Algorithm for Shortest Paths
 */

// Graph Examples
const EXAMPLES = {
  small: {
    name: "Small Graph (n=8)",
    vertices: [
      { id: 's', x: 100, y: 250, label: 's' },
      { id: 'A', x: 220, y: 150, label: 'A' },
      { id: 'B', x: 220, y: 350, label: 'B' },
      { id: 'C', x: 340, y: 250, label: 'C' },
      { id: 'D', x: 460, y: 150, label: 'D' },
      { id: 'E', x: 460, y: 350, label: 'E' },
      { id: 'F', x: 580, y: 250, label: 'F' },
      { id: 'G', x: 700, y: 250, label: 'G' }
    ],
    edges: [
      { from: 's', to: 'A', weight: 1 },
      { from: 's', to: 'B', weight: 3 },
      { from: 's', to: 'C', weight: 5 },
      { from: 'A', to: 'D', weight: 2 },
      { from: 'A', to: 'E', weight: 4 },
      { from: 'B', to: 'D', weight: 1 },
      { from: 'B', to: 'F', weight: 2 },
      { from: 'C', to: 'E', weight: 1 },
      { from: 'C', to: 'G', weight: 3 },
      { from: 'D', to: 'G', weight: 1 },
      { from: 'E', to: 'F', weight: 1 },
      { from: 'E', to: 'G', weight: 2 },
      { from: 'F', to: 'G', weight: 1 }
    ],
    k: 2,
    t: 4,
    source: 's'
  },
  medium: {
    name: "Medium Graph (n=16)",
    vertices: [
      { id: 's', x: 80, y: 80, label: 's' },
      { id: '1', x: 230, y: 80, label: '1' },
      { id: '2', x: 380, y: 80, label: '2' },
      { id: '3', x: 530, y: 80, label: '3' },
      { id: '4', x: 80, y: 200, label: '4' },
      { id: '5', x: 230, y: 200, label: '5' },
      { id: '6', x: 380, y: 200, label: '6' },
      { id: '7', x: 530, y: 200, label: '7' },
      { id: '8', x: 80, y: 320, label: '8' },
      { id: '9', x: 230, y: 320, label: '9' },
      { id: '10', x: 380, y: 320, label: '10' },
      { id: '11', x: 530, y: 320, label: '11' },
      { id: '12', x: 80, y: 440, label: '12' },
      { id: '13', x: 230, y: 440, label: '13' },
      { id: '14', x: 380, y: 440, label: '14' },
      { id: '15', x: 530, y: 440, label: '15' }
    ],
    edges: [
      { from: 's', to: '1', weight: 2 }, { from: 's', to: '4', weight: 5 },
      { from: '1', to: '2', weight: 1 }, { from: '1', to: '5', weight: 3 },
      { from: '2', to: '3', weight: 2 }, { from: '2', to: '6', weight: 4 },
      { from: '3', to: '7', weight: 1 }, { from: '4', to: '5', weight: 2 },
      { from: '4', to: '8', weight: 3 }, { from: '5', to: '6', weight: 1 },
      { from: '5', to: '9', weight: 2 }, { from: '6', to: '7', weight: 1 },
      { from: '6', to: '10', weight: 3 }, { from: '7', to: '11', weight: 2 },
      { from: '8', to: '9', weight: 1 }, { from: '8', to: '12', weight: 4 },
      { from: '9', to: '10', weight: 2 }, { from: '9', to: '13', weight: 1 },
      { from: '10', to: '11', weight: 1 }, { from: '10', to: '14', weight: 2 },
      { from: '11', to: '15', weight: 1 }, { from: '12', to: '13', weight: 2 },
      { from: '13', to: '14', weight: 1 }, { from: '14', to: '15', weight: 1 }
    ],
    k: 2,
    t: 5,
    source: 's'
  },
  sparse: {
    name: "Sparse Network (n=12)",
    vertices: [
      { id: 's', x: 100, y: 250, label: 's' },
      { id: '1', x: 250, y: 100, label: '1' },
      { id: '2', x: 250, y: 400, label: '2' },
      { id: '3', x: 350, y: 180, label: '3' },
      { id: '4', x: 350, y: 320, label: '4' },
      { id: '5', x: 450, y: 130, label: '5' },
      { id: '6', x: 450, y: 370, label: '6' },
      { id: '7', x: 550, y: 100, label: '7' },
      { id: '8', x: 550, y: 250, label: '8' },
      { id: '9', x: 550, y: 400, label: '9' },
      { id: '10', x: 650, y: 180, label: '10' },
      { id: '11', x: 650, y: 320, label: '11' }
    ],
    edges: [
      { from: 's', to: '1', weight: 2 }, { from: 's', to: '2', weight: 4 },
      { from: '1', to: '3', weight: 3 }, { from: '1', to: '5', weight: 5 },
      { from: '2', to: '4', weight: 2 }, { from: '2', to: '6', weight: 3 },
      { from: '3', to: '5', weight: 1 }, { from: '3', to: '8', weight: 4 },
      { from: '4', to: '6', weight: 1 }, { from: '4', to: '8', weight: 3 },
      { from: '5', to: '7', weight: 2 }, { from: '5', to: '10', weight: 3 },
      { from: '6', to: '9', weight: 2 }, { from: '6', to: '11', weight: 4 },
      { from: '7', to: '10', weight: 1 }, { from: '8', to: '10', weight: 2 },
      { from: '8', to: '11', weight: 2 }, { from: '9', to: '11', weight: 1 }
    ],
    k: 2,
    t: 4,
    source: 's'
  }
};

// Algorithm State
class AlgorithmState {
  constructor(example) {
    this.example = example;
    this.n = example.vertices.length;
    this.k = example.k;
    this.t = example.t;
    
    // Distance and predecessor arrays
    this.dist = {};
    this.pred = {};
    example.vertices.forEach(v => {
      this.dist[v.id] = v.id === example.source ? 0 : Infinity;
      this.pred[v.id] = null;
    });
    
    // Algorithm state
    this.phase = 'Not started';
    this.round = 0;
    this.frontier = new Set([example.source]);
    this.settled = new Set();
    this.pivots = new Set();
    this.processedVertices = new Set();
    
    // Metrics
    this.relaxations = 0;
    this.comparisons = 0;
    
    // Execution log
    this.log = [];
  }
  
  addLog(message, highlight = false) {
    this.log.push({ message, highlight, step: this.log.length + 1 });
  }
}

// Global state
let currentExample = 'small';
let state = null;
let canvas, ctx;
let animationSpeed = 5;
let isRunning = false;
let autoRunInterval = null;

// Initialize on page load
document.addEventListener('DOMContentLoaded', () => {
  canvas = document.getElementById('graph-canvas');
  if (!canvas) {
    console.error('Canvas element not found!');
    return;
  }
  
  ctx = canvas.getContext('2d');
  
  setupCanvas();
  setupControls();
  setupExampleSelector();
  
  // Load initial example
  loadExample('small');
});

// Canvas Setup
function setupCanvas() {
  // Set canvas size
  const container = canvas.parentElement;
  canvas.width = container.clientWidth - 40; // Account for padding
  canvas.height = 500;
  
  // Handle resize
  window.addEventListener('resize', () => {
    const container = canvas.parentElement;
    canvas.width = container.clientWidth - 40;
    canvas.height = 500;
    if (state) drawGraph();
  });
}

// Control Setup
function setupControls() {
  document.getElementById('btn-init').addEventListener('click', initializeAlgorithm);
  document.getElementById('btn-reset').addEventListener('click', resetAlgorithm);
  document.getElementById('btn-step').addEventListener('click', stepForward);
  document.getElementById('btn-run').addEventListener('click', toggleAutoRun);
  
  const speedSlider = document.getElementById('speed-slider');
  speedSlider.addEventListener('input', (e) => {
    animationSpeed = parseInt(e.target.value);
    document.getElementById('speed-value').textContent = animationSpeed;
  });
}

// Example Selector Setup
function setupExampleSelector() {
  const cards = document.querySelectorAll('.example-card');
  cards.forEach(card => {
    card.addEventListener('click', () => {
      const exampleId = card.getAttribute('data-example');
      if (exampleId === 'custom') {
        alert('Custom graph builder coming soon! Try the 3 pre-built examples.');
        return;
      }
      
      cards.forEach(c => c.classList.remove('active'));
      card.classList.add('active');
      loadExample(exampleId);
    });
  });
}

// Load Example
function loadExample(exampleId) {
  currentExample = exampleId;
  const example = EXAMPLES[exampleId];
  
  if (!example) return;
  
  // Stop any running animation
  if (autoRunInterval) {
    clearInterval(autoRunInterval);
    autoRunInterval = null;
    isRunning = false;
    document.getElementById('btn-run').textContent = 'Auto Run';
  }
  
  state = new AlgorithmState(example);
  state.addLog(`Loaded ${example.name}: n=${state.n}, m=${example.edges.length}`, true);
  state.addLog(`Click "Initialize" to start the algorithm`);
  
  drawGraph();
  updateUI();
  
  // Disable step buttons until initialized
  document.getElementById('btn-step').disabled = true;
  document.getElementById('btn-run').disabled = true;
}

// Initialize Algorithm
function initializeAlgorithm() {
  if (!state) return;
  
  // Reset if already run
  if (state.phase === 'Complete') {
    loadExample(currentExample);
    state = new AlgorithmState(EXAMPLES[currentExample]);
  }
  
  state.phase = 'Initialized';
  state.addLog(`✓ Algorithm initialized on ${state.example.name}`, true);
  state.addLog(`Parameters: n=${state.n}, m=${state.example.edges.length}, k=${state.k}, t=${state.t}`);
  state.addLog(`Source vertex: ${state.example.source} with dist[${state.example.source}]=0`);
  state.addLog(`Ready to start FindPivots phase...`);
  
  document.getElementById('btn-step').disabled = false;
  document.getElementById('btn-run').disabled = false;
  
  updateUI();
  drawGraph();
}

// Reset Algorithm
function resetAlgorithm() {
  if (autoRunInterval) {
    clearInterval(autoRunInterval);
    autoRunInterval = null;
    isRunning = false;
    document.getElementById('btn-run').textContent = 'Auto Run';
  }
  
  loadExample(currentExample);
}

// Step Forward
function stepForward() {
  if (!state || state.phase === 'Not started') return;
  
  if (state.phase === 'Complete') {
    state.addLog('⚠ Algorithm already complete! Click Reset to run again.');
    return;
  }
  
  if (state.phase === 'Initialized') {
    startFindPivots();
  } else if (state.phase === 'FindPivots' && state.round < state.k) {
    executeRelaxationRound();
  } else if (state.phase === 'FindPivots' && state.round >= state.k) {
    identifyPivots();
  }
  
  updateUI();
  drawGraph();
}

// Toggle Auto Run
function toggleAutoRun() {
  if (state.phase === 'Complete') {
    alert('Algorithm complete! Click Reset to run again.');
    return;
  }
  
  if (isRunning) {
    clearInterval(autoRunInterval);
    autoRunInterval = null;
    isRunning = false;
    document.getElementById('btn-run').textContent = 'Auto Run';
  } else {
    isRunning = true;
    document.getElementById('btn-run').textContent = 'Pause';
    const delay = 1000 / animationSpeed;
    
    autoRunInterval = setInterval(() => {
      stepForward();
      if (state.phase === 'Complete') {
        clearInterval(autoRunInterval);
        autoRunInterval = null;
        isRunning = false;
        document.getElementById('btn-run').textContent = 'Auto Run';
      }
    }, delay);
  }
}

// Algorithm Steps

function startFindPivots() {
  state.phase = 'FindPivots';
  state.round = 0;
  state.addLog(`━━━ Phase 1: FindPivots ━━━`, true);
  state.addLog(`Starting k=${state.k} rounds of edge relaxations`);
  state.addLog(`Each round relaxes edges from newly settled vertices`);
}

function executeRelaxationRound() {
  state.round++;
  state.addLog(`━━━ Round ${state.round}/${state.k} ━━━`, true);
  
  // Determine which vertices to process this round
  let toProcess = [];
  if (state.round === 1) {
    // First round: process from source
    toProcess = [state.example.source];
    state.processedVertices.add(state.example.source);
  } else {
    // Subsequent rounds: process newly settled vertices from last round
    toProcess = Array.from(state.settled).filter(v => !state.processedVertices.has(v));
    toProcess.forEach(v => state.processedVertices.add(v));
  }
  
  if (toProcess.length === 0) {
    state.addLog(`No new vertices to process this round`);
    return;
  }
  
  state.addLog(`Processing ${toProcess.length} vertex(vertices): {${toProcess.join(', ')}}`);
  
  // Relax all outgoing edges from these vertices
  let relaxed = 0;
  let updated = 0;
  
  toProcess.forEach(u => {
    const edges = state.example.edges.filter(e => e.from === u);
    edges.forEach(edge => {
      state.relaxations++;
      state.comparisons++;
      
      const newDist = state.dist[edge.from] + edge.weight;
      if (newDist < state.dist[edge.to]) {
        const oldDist = state.dist[edge.to];
        state.dist[edge.to] = newDist;
        state.pred[edge.to] = edge.from;
        state.settled.add(edge.to);
        relaxed++;
        updated++;
        
        state.addLog(
          `  Relaxed ${edge.from}→${edge.to} (weight ${edge.weight}): dist[${edge.to}] = ${oldDist === Infinity ? '∞' : oldDist} → ${newDist}`
        );
      } else {
        relaxed++;
      }
    });
  });
  
  state.addLog(`Round ${state.round} complete: ${relaxed} edges relaxed, ${updated} distances updated`);
  state.addLog(`Settled set |W| = ${state.settled.size} vertices`);
}

function identifyPivots() {
  state.addLog(`━━━ Identifying Pivots ━━━`, true);
  state.addLog(`Building shortest-path forest rooted at frontier vertices...`);
  
  // Build shortest-path forest
  const forest = new Map();
  state.example.vertices.forEach(v => forest.set(v.id, []));
  
  state.settled.forEach(v => {
    if (state.pred[v] && state.frontier.has(state.pred[v])) {
      forest.get(state.pred[v]).push(v);
    }
  });
  
  // Compute tree sizes recursively
  const treeSizes = new Map();
  
  function computeTreeSize(v) {
    if (treeSizes.has(v)) return treeSizes.get(v);
    
    let size = 1;
    const children = forest.get(v) || [];
    children.forEach(child => {
      size += computeTreeSize(child);
    });
    treeSizes.set(v, size);
    return size;
  }
  
  // Compute tree sizes for all frontier vertices
  state.frontier.forEach(v => {
    const size = computeTreeSize(v);
    state.addLog(`Tree(${v}): ${size} vertices`);
    
    if (size >= state.k) {
      state.pivots.add(v);
      state.addLog(`  → ${v} is a PIVOT (tree size ${size} ≥ k=${state.k})`, true);
    }
  });
  
  if (state.pivots.size === 0) {
    state.addLog(`No pivots found (all tree sizes < k=${state.k})`);
    state.addLog(`Frontier remains unchanged: |S| = ${state.frontier.size}`);
  } else {
    state.addLog(`✓ Pivots identified: {${Array.from(state.pivots).join(', ')}}`, true);
    state.addLog(`Frontier reduction: ${state.frontier.size} vertices → ${state.pivots.size} pivots`);
  }
  
  state.phase = 'Complete';
  finalizeAlgorithm();
}

function finalizeAlgorithm() {
  state.addLog(`━━━ Algorithm Complete ━━━`, true);
  state.addLog(`Total edge relaxations: ${state.relaxations}`);
  state.addLog(`Total comparisons: ${state.comparisons}`);
  state.addLog(`Vertices settled: ${state.settled.size}/${state.n}`);
  state.addLog(`Pivots found: ${state.pivots.size}`);
  state.addLog(` `);
  state.addLog(`✓ PathOpt frontier reduction successful!`, true);
}

// Drawing Functions

function drawGraph() {
  if (!ctx || !state) return;
  
  // Clear canvas
  ctx.fillStyle = 'rgba(15, 23, 42, 0.9)';
  ctx.fillRect(0, 0, canvas.width, canvas.height);
  
  // Calculate scale based on canvas size
  const maxX = Math.max(...state.example.vertices.map(v => v.x));
  const maxY = Math.max(...state.example.vertices.map(v => v.y));
  const scaleX = (canvas.width - 80) / maxX;
  const scaleY = (canvas.height - 80) / maxY;
  const scale = Math.min(scaleX, scaleY, 1);
  const offsetX = 40;
  const offsetY = 40;
  
  // Draw edges first (behind vertices)
  state.example.edges.forEach(edge => {
    const from = state.example.vertices.find(v => v.id === edge.from);
    const to = state.example.vertices.find(v => v.id === edge.to);
    
    if (!from || !to) return;
    
    const fromX = from.x * scale + offsetX;
    const fromY = from.y * scale + offsetY;
    const toX = to.x * scale + offsetX;
    const toY = to.y * scale + offsetY;
    
    drawEdge(fromX, fromY, toX, toY, edge.weight, state.pred[edge.to] === edge.from);
  });
  
  // Draw vertices on top
  state.example.vertices.forEach(vertex => {
    const x = vertex.x * scale + offsetX;
    const y = vertex.y * scale + offsetY;
    drawVertex(vertex, x, y);
  });
}

function drawEdge(fromX, fromY, toX, toY, weight, isInTree) {
  // Edge line
  const gradient = ctx.createLinearGradient(fromX, fromY, toX, toY);
  if (isInTree) {
    gradient.addColorStop(0, 'rgba(16, 185, 129, 0.6)');
    gradient.addColorStop(1, 'rgba(6, 182, 212, 0.6)');
    ctx.lineWidth = 3;
  } else {
    gradient.addColorStop(0, 'rgba(139, 92, 246, 0.2)');
    gradient.addColorStop(1, 'rgba(6, 182, 212, 0.2)');
    ctx.lineWidth = 2;
  }
  
  ctx.strokeStyle = gradient;
  
  // Calculate arrow endpoint (stop before vertex circle)
  const angle = Math.atan2(toY - fromY, toX - fromX);
  const radius = 22;
  const endX = toX - Math.cos(angle) * radius;
  const endY = toY - Math.sin(angle) * radius;
  const startX = fromX + Math.cos(angle) * radius;
  const startY = fromY + Math.sin(angle) * radius;
  
  // Draw line
  ctx.beginPath();
  ctx.moveTo(startX, startY);
  ctx.lineTo(endX, endY);
  ctx.stroke();
  
  // Draw arrowhead
  const arrowSize = 8;
  ctx.fillStyle = isInTree ? 'rgba(16, 185, 129, 0.8)' : 'rgba(139, 92, 246, 0.4)';
  ctx.beginPath();
  ctx.moveTo(endX, endY);
  ctx.lineTo(
    endX - arrowSize * Math.cos(angle - Math.PI / 6),
    endY - arrowSize * Math.sin(angle - Math.PI / 6)
  );
  ctx.lineTo(
    endX - arrowSize * Math.cos(angle + Math.PI / 6),
    endY - arrowSize * Math.sin(angle + Math.PI / 6)
  );
  ctx.closePath();
  ctx.fill();
  
  // Draw weight label
  const midX = (fromX + toX) / 2;
  const midY = (fromY + toY) / 2;
  const offsetDist = 15;
  const labelX = midX + offsetDist * Math.sin(angle);
  const labelY = midY - offsetDist * Math.cos(angle);
  
  ctx.fillStyle = 'rgba(0, 0, 0, 0.6)';
  ctx.beginPath();
  ctx.arc(labelX, labelY, 12, 0, 2 * Math.PI);
  ctx.fill();
  
  ctx.fillStyle = '#fff';
  ctx.font = 'bold 11px monospace';
  ctx.textAlign = 'center';
  ctx.textBaseline = 'middle';
  ctx.fillText(weight, labelX, labelY);
}

function drawVertex(vertex, x, y) {
  const id = vertex.id;
  const isSource = id === state.example.source;
  const isFrontier = state.frontier.has(id);
  const isPivot = state.pivots.has(id);
  const isSettled = state.settled.has(id);
  
  // Determine color
  let fillColor, strokeColor;
  if (isSource) {
    fillColor = 'rgba(245, 158, 11, 0.9)'; // Orange for source
    strokeColor = '#f59e0b';
  } else if (isPivot) {
    fillColor = 'rgba(6, 182, 212, 0.9)'; // Cyan for pivots
    strokeColor = '#06b6d4';
  } else if (isFrontier) {
    fillColor = 'rgba(139, 92, 246, 0.8)'; // Purple for frontier
    strokeColor = '#8b5cf6';
  } else if (isSettled) {
    fillColor = 'rgba(16, 185, 129, 0.7)'; // Green for settled
    strokeColor = '#10b981';
  } else {
    fillColor = 'rgba(30, 30, 60, 0.8)'; // Dark for unvisited
    strokeColor = 'rgba(139, 92, 246, 0.3)';
  }
  
  // Draw circle
  ctx.fillStyle = fillColor;
  ctx.strokeStyle = strokeColor;
  ctx.lineWidth = 3;
  ctx.beginPath();
  ctx.arc(x, y, 20, 0, 2 * Math.PI);
  ctx.fill();
  ctx.stroke();
  
  // Draw label
  ctx.fillStyle = '#fff';
  ctx.font = 'bold 14px sans-serif';
  ctx.textAlign = 'center';
  ctx.textBaseline = 'middle';
  ctx.fillText(vertex.label, x, y);
  
  // Draw distance below vertex
  if (state.dist[id] !== Infinity) {
    ctx.fillStyle = 'rgba(0, 0, 0, 0.7)';
    ctx.fillRect(x - 18, y + 24, 36, 18);
    
    ctx.fillStyle = '#06b6d4';
    ctx.font = 'bold 11px monospace';
    ctx.fillText(`d=${state.dist[id]}`, x, y + 33);
  }
}

// UI Update Functions

function updateUI() {
  if (!state) return;
  
  // Update state display
  document.getElementById('state-phase').textContent = state.phase;
  document.getElementById('state-round').textContent = state.phase === 'FindPivots' ? `${state.round}/${state.k}` : '—';
  document.getElementById('state-frontier').textContent = state.frontier.size;
  document.getElementById('state-settled').textContent = state.settled.size;
  document.getElementById('state-k').textContent = state.k;
  document.getElementById('state-t').textContent = state.t;
  
  // Update metrics
  document.getElementById('metric-relaxations').textContent = state.relaxations;
  document.getElementById('metric-comparisons').textContent = state.comparisons;
  document.getElementById('metric-pivots').textContent = state.pivots.size;
  document.getElementById('metric-settled').textContent = state.settled.size;
  
  // Update execution log
  updateExecutionLog();
  
  // Update distance table
  updateDistanceTable();
}

function updateExecutionLog() {
  const logContainer = document.getElementById('calculation-log');
  const logContent = logContainer.querySelector('h3').nextElementSibling;
  
  // Remove all log entries except the title
  while (logContainer.children.length > 1) {
    logContainer.removeChild(logContainer.lastChild);
  }
  
  // Add new log entries (show last 15 to avoid overflow)
  const recentLogs = state.log.slice(-15);
  recentLogs.forEach(entry => {
    const div = document.createElement('div');
    div.className = 'log-entry' + (entry.highlight ? ' highlight' : '');
    
    const stepSpan = document.createElement('span');
    stepSpan.className = 'step-num';
    stepSpan.textContent = `[${entry.step}]`;
    
    div.appendChild(stepSpan);
    div.appendChild(document.createTextNode(entry.message));
    logContainer.appendChild(div);
  });
  
  // Auto-scroll to bottom
  logContainer.scrollTop = logContainer.scrollHeight;
}

function updateDistanceTable() {
  const tableDiv = document.getElementById('distance-table');
  
  if (state.phase === 'Not started') {
    tableDiv.innerHTML = '<p style="color: var(--color-secondary); font-size: 0.9rem;">Initialize to see distances</p>';
    return;
  }
  
  // Create distance table
  let html = '<div style="font-family: monospace; font-size: 0.9rem;">';
  
  state.example.vertices.forEach(v => {
    const dist = state.dist[v.id];
    const distStr = dist === Infinity ? '∞' : dist;
    const pred = state.pred[v.id] || '—';
    
    let color = 'var(--color-secondary)';
    if (state.pivots.has(v.id)) color = 'var(--color-cyan)';
    else if (state.frontier.has(v.id)) color = 'var(--color-purple)';
    else if (state.settled.has(v.id)) color = '#10b981';
    
    html += `<div style="display: flex; justify-content: space-between; padding: 0.4rem 0.5rem; border-bottom: 1px solid rgba(139, 92, 246, 0.1); color: ${color};">`;
    html += `<span><strong>${v.label}:</strong></span>`;
    html += `<span>d=${distStr}, π=${pred}</span>`;
    html += `</div>`;
  });
  
  html += '</div>';
  tableDiv.innerHTML = html;
}
