(function () {
  const state = {
    bodyPage: document.body.getAttribute('data-page') || 'home'
  };

  /* -----------------------------------------------------------
   * Minimal MathJax-Compatible Renderer (offline, subset)
   * ----------------------------------------------------------- */
  const MathRenderer = (() => {
    const greekMap = {
      '\\Phi': '&Phi;',
      '\\phi': '&phi;',
      '\\Delta': '&Delta;',
      '\\delta': '&delta;',
      '\\ell': '&ell;',
      '\\lambda': '&lambda;',
      '\\Lambda': '&Lambda;',
      '\\alpha': '&alpha;',
      '\\beta': '&beta;',
      '\\gamma': '&gamma;',
      '\\infty': '&infin;'
    };

    const opMap = {
      '\\log': 'log',
      '\\min': 'min',
      '\\max': 'max',
      '\\sum': '∑',
      '\\mathrm': '',
      '\\mathbb': '',
      '\\mathcal': ''
    };

    const replaceCommands = (str) => {
      let result = str;
      Object.keys(greekMap).forEach((key) => {
        result = result.replace(new RegExp(key, 'g'), greekMap[key]);
      });
      Object.keys(opMap).forEach((key) => {
        result = result.replace(new RegExp(key, 'g'), opMap[key]);
      });
      result = result
        .replace(/\\rightarrow/g, '→')
        .replace(/\\Rightarrow/g, '⇒')
        .replace(/\\leq?/g, '≤')
        .replace(/\\geq?/g, '≥')
        .replace(/\\cdot/g, '·')
        .replace(/\\times/g, '×')
        .replace(/\\forall/g, '∀')
        .replace(/\\exists/g, '∃')
        .replace(/\\in/g, '∈')
        .replace(/\\cup/g, '∪')
        .replace(/\\cap/g, '∩')
        .replace(/\\mathbb{R}/g, '&#8477;')
        .replace(/\\mathbb{Z}/g, '&#8484;')
        .replace(/\\mathbb{N}/g, '&#8469;')
        .replace(/\\sqrt\{([^}]*)\}/g, (_, inner) => `√(${inner})`)
        .replace(/\\lfloor/g, '⌊')
        .replace(/\\rfloor/g, '⌋')
        .replace(/\\lceil/g, '⌈')
        .replace(/\\rceil/g, '⌉');
      return result;
    };

    const renderTeX = (tex) => {
      let html = tex;
      html = html.replace(/\\frac\{([^}]*)\}\{([^}]*)\}/g, (_, num, den) => {
        return `<span class="frac"><span class="num">${renderTeX(num)}</span><span class="bar"></span><span class="den">${renderTeX(den)}</span></span>`;
      });
      html = html.replace(/\^\{([^}]*)\}/g, (_, sup) => `<sup>${renderTeX(sup)}</sup>`);
      html = html.replace(/_\{([^}]*)\}/g, (_, sub) => `<sub>${renderTeX(sub)}</sub>`);
      html = html.replace(/\^([A-Za-z0-9]+)/g, (_, sup) => `<sup>${sup}</sup>`);
      html = html.replace(/_([A-Za-z0-9]+)/g, (_, sub) => `<sub>${sub}</sub>`);
      html = html.replace(/\\mathrm\{([^}]*)\}/g, (_, inner) => `<span class="math-op">${inner}</span>`);
      html = html.replace(/\\mathbb\{([^}]*)\}/g, (_, inner) => `<span class="math-op">${inner}</span>`);
      html = html.replace(/\\mathcal\{([^}]*)\}/g, (_, inner) => `<span class="math-op">${inner}</span>`);
      html = replaceCommands(html);
      return html;
    };

    const renderElements = (selector, className) => {
      document.querySelectorAll(selector).forEach((el) => {
        if (el.dataset.rendered) return;
        const tex = el.textContent.trim();
        el.innerHTML = renderTeX(tex);
        el.dataset.rendered = 'true';
        el.classList.add(className);
      });
    };

    const init = () => {
      renderElements('.math-inline', 'math-inline-ready');
      renderElements('.math-display', 'math-display-ready');
    };

    return { init };
  })();

  window.MathJax = window.MathJax || {};
  window.MathJax.startup = window.MathJax.startup || {};
  window.MathJax.startup.ready = () => {
    MathRenderer.init();
  };
  window.MathJax.typesetPromise = () => {
    MathRenderer.init();
    return Promise.resolve();
  };

  document.addEventListener('DOMContentLoaded', () => {
    MathRenderer.init();
    if (state.bodyPage === 'interactive') {
      AlgorithmLab.init();
    }
    if (state.bodyPage === 'home') {
      HomePage.init();
    }
    if (state.bodyPage === 'paper') {
      PaperPage.init();
    }
  });

  /* -----------------------------------------------------------
   * Home Page Micro-interactions
   * ----------------------------------------------------------- */
  const HomePage = (() => {
    const init = () => {
      const cards = document.querySelectorAll('.card');
      if (!cards.length) return;
      let idx = 0;
      setInterval(() => {
        cards.forEach((card, i) => {
          card.style.transform = i === idx ? 'translateY(-6px)' : 'translateY(0)';
          card.style.transition = 'transform 320ms ease';
        });
        idx = (idx + 1) % cards.length;
      }, 4200);
    };
    return { init };
  })();

  /* -----------------------------------------------------------
   * Documentation Page Enhancements (floating scroll spy)
   * ----------------------------------------------------------- */
  const PaperPage = (() => {
    const init = () => {
      const tocLinks = document.querySelectorAll('.paper-aside a');
      const sections = Array.from(tocLinks).map((link) => {
        const id = link.getAttribute('href');
        const section = document.querySelector(id);
        return { link, section };
      });
      window.addEventListener('scroll', () => {
        const y = window.scrollY + 120;
        let active = null;
        sections.forEach(({ link, section }) => {
          if (!section) return;
          if (section.offsetTop <= y) {
            active = link;
          }
        });
        tocLinks.forEach((link) => link.classList.remove('active'));
        if (active) active.classList.add('active');
      });
    };
    return { init };
  })();

  /* -----------------------------------------------------------
   * Data Structures for Interactive Lab
   * ----------------------------------------------------------- */
  class GraphModel {
    constructor() {
      this.nodes = new Map();
      this.edges = new Map();
      this.nextNodeId = 0;
      this.nextEdgeId = 0;
      this.source = null;
    }

    addNode(x, y, extras = {}) {
      const id = this.nextNodeId++;
      const node = { id, x, y, label: extras.label || `v${id}`, data: extras.data || {} };
      this.nodes.set(id, node);
      if (this.source === null) {
        this.source = id;
      }
      return node;
    }

    removeNode(id) {
      if (!this.nodes.has(id)) return;
      this.nodes.delete(id);
      [...this.edges.values()].forEach((edge) => {
        if (edge.from === id || edge.to === id) {
          this.edges.delete(edge.id);
        }
      });
      if (this.source === id) {
        const first = this.nodes.keys().next();
        this.source = first.done ? null : first.value;
      }
    }

    addEdge(from, to, weight = 1, extras = {}) {
      if (!this.nodes.has(from) || !this.nodes.has(to) || from === to) return null;
      const id = this.nextEdgeId++;
      const edge = { id, from, to, weight: Number(weight), data: extras.data || {} };
      this.edges.set(id, edge);
      return edge;
    }

    removeEdge(edgeId) {
      this.edges.delete(edgeId);
    }

    getNodeAt(x, y, radius = 18) {
      const threshold = radius * radius;
      for (const node of this.nodes.values()) {
        const dx = node.x - x;
        const dy = node.y - y;
        if (dx * dx + dy * dy <= threshold) {
          return node;
        }
      }
      return null;
    }

    getEdgeAt(x, y) {
      let match = null;
      let minDist = 14;
      for (const edge of this.edges.values()) {
        const from = this.nodes.get(edge.from);
        const to = this.nodes.get(edge.to);
        if (!from || !to) continue;
        const dist = pointLineDistance(x, y, from.x, from.y, to.x, to.y);
        if (dist < minDist) {
          minDist = dist;
          match = edge;
        }
      }
      return match;
    }

    setSource(id) {
      if (this.nodes.has(id)) {
        this.source = id;
      }
    }

    clear() {
      this.nodes.clear();
      this.edges.clear();
      this.nextNodeId = 0;
      this.nextEdgeId = 0;
      this.source = null;
    }

    toAdjacency() {
      const adj = new Map();
      this.nodes.forEach((node) => {
        adj.set(node.id, []);
      });
      this.edges.forEach((edge) => {
        if (adj.has(edge.from)) {
          adj.get(edge.from).push({ to: edge.to, weight: edge.weight, id: edge.id });
        }
      });
      return { nodes: Array.from(this.nodes.values()), edges: Array.from(this.edges.values()), adj, source: this.source };
    }

    clone() {
      const copy = new GraphModel();
      this.nodes.forEach((node) => copy.addNode(node.x, node.y, { label: node.label, data: node.data }));
      this.edges.forEach((edge) => copy.addEdge(edge.from, edge.to, edge.weight, { data: edge.data }));
      copy.nextNodeId = this.nextNodeId;
      copy.nextEdgeId = this.nextEdgeId;
      copy.source = this.source;
      return copy;
    }
  }

  /* -----------------------------------------------------------
   * Canvas Renderer
   * ----------------------------------------------------------- */
  class CanvasRenderer {
    constructor(canvas, model) {
      this.canvas = canvas;
      this.ctx = canvas.getContext('2d');
      this.model = model;
      this.padding = 48;
      this.state = {
        active: null,
        settled: new Set(),
        frontier: new Set(),
        highlightEdge: null,
        tempEdge: null
      };
      this.animationId = null;
      this.tick = 0;
      this.render = this.render.bind(this);
      this.render();
    }

    setModel(model) {
      this.model = model;
    }

    updateState(partial) {
      this.state = { ...this.state, ...partial };
    }

    getBounds() {
      const nodes = Array.from(this.model.nodes.values());
      if (!nodes.length) return { minX: 0, maxX: 1, minY: 0, maxY: 1 };
      const xs = nodes.map((n) => n.x);
      const ys = nodes.map((n) => n.y);
      return {
        minX: Math.min(...xs),
        maxX: Math.max(...xs),
        minY: Math.min(...ys),
        maxY: Math.max(...ys)
      };
    }

    worldToCanvas(x, y) {
      const { minX, maxX, minY, maxY } = this.getBounds();
      const width = maxX - minX || 1;
      const height = maxY - minY || 1;
      const scaleX = (this.canvas.width - this.padding * 2) / width;
      const scaleY = (this.canvas.height - this.padding * 2) / height;
      const scale = Math.min(scaleX, scaleY);
      const offsetX = (this.canvas.width - width * scale) / 2;
      const offsetY = (this.canvas.height - height * scale) / 2;
      const cx = offsetX + (x - minX) * scale;
      const cy = offsetY + (y - minY) * scale;
      return { x: cx, y: cy };
    }

    render() {
      cancelAnimationFrame(this.animationId);
      const { ctx, canvas } = this;
      ctx.clearRect(0, 0, canvas.width, canvas.height);
      ctx.fillStyle = '#071225';
      ctx.fillRect(0, 0, canvas.width, canvas.height);
      ctx.lineWidth = 2;
      ctx.lineCap = 'round';
      ctx.font = '14px "Segoe UI", sans-serif';

      // Draw edges
      this.model.edges.forEach((edge) => {
        const from = this.model.nodes.get(edge.from);
        const to = this.model.nodes.get(edge.to);
        if (!from || !to) return;
        const start = this.worldToCanvas(from.x, from.y);
        const end = this.worldToCanvas(to.x, to.y);
        const isHighlight = this.state.highlightEdge === edge.id;
        ctx.strokeStyle = isHighlight ? '#f7cf6b' : 'rgba(77,163,255,0.55)';
        ctx.beginPath();
        ctx.moveTo(start.x, start.y);
        ctx.lineTo(end.x, end.y);
        ctx.stroke();
        drawArrowhead(ctx, start.x, start.y, end.x, end.y, isHighlight ? '#f7cf6b' : '#4da3ff');
        const labelX = start.x * 0.7 + end.x * 0.3;
        const labelY = start.y * 0.7 + end.y * 0.3;
        ctx.fillStyle = '#9fb3d1';
        ctx.fillText(edge.weight.toFixed(1).replace(/\.0$/, ''), labelX, labelY);
      });
      
      // Draw temporary edge during drag
      if (this.state.tempEdge) {
        const { from, x, y } = this.state.tempEdge;
        const start = this.worldToCanvas(from.x, from.y);
        ctx.strokeStyle = 'rgba(247,207,107,0.8)';
        ctx.lineWidth = 3;
        ctx.setLineDash([8, 4]);
        ctx.beginPath();
        ctx.moveTo(start.x, start.y);
        ctx.lineTo(x, y); // x, y are already canvas coordinates from getCanvasPoint
        ctx.stroke();
        ctx.setLineDash([]);
        ctx.lineWidth = 2;
        drawArrowhead(ctx, start.x, start.y, x, y, 'rgba(247,207,107,1.0)');
      }

      // Draw nodes
      this.model.nodes.forEach((node) => {
        const { x, y } = this.worldToCanvas(node.x, node.y);
        const radius = this.state.frontier.has(node.id) ? 14 : 11;
        ctx.beginPath();
        ctx.arc(x, y, radius, 0, Math.PI * 2);
        let fill = '#12294c';
        if (this.state.settled.has(node.id)) fill = '#1f6fb2';
        if (this.state.active === node.id) fill = '#f7cf6b';
        ctx.fillStyle = fill;
        ctx.fill();
        ctx.strokeStyle = this.model.source === node.id ? '#f7cf6b' : 'rgba(77,163,255,0.55)';
        ctx.stroke();
        ctx.fillStyle = '#f2f6ff';
        ctx.textAlign = 'center';
        ctx.textBaseline = 'middle';
        ctx.fillText(node.label, x, y);
      });

      this.animationId = requestAnimationFrame(this.render);
    }
  }

  /* -----------------------------------------------------------
   * Geometry Helpers
   * ----------------------------------------------------------- */
  function pointLineDistance(px, py, x1, y1, x2, y2) {
    const A = px - x1;
    const B = py - y1;
    const C = x2 - x1;
    const D = y2 - y1;
    const dot = A * C + B * D;
    const lenSq = C * C + D * D;
    let param = -1;
    if (lenSq !== 0) param = dot / lenSq;
    let xx, yy;
    if (param < 0) {
      xx = x1;
      yy = y1;
    } else if (param > 1) {
      xx = x2;
      yy = y2;
    } else {
      xx = x1 + param * C;
      yy = y1 + param * D;
    }
    const dx = px - xx;
    const dy = py - yy;
    return Math.sqrt(dx * dx + dy * dy);
  }

  function drawArrowhead(ctx, x0, y0, x1, y1, color) {
    const angle = Math.atan2(y1 - y0, x1 - x0);
    const size = 10;
    ctx.fillStyle = color;
    ctx.beginPath();
    ctx.moveTo(x1, y1);
    ctx.lineTo(x1 - size * Math.cos(angle - Math.PI / 6), y1 - size * Math.sin(angle - Math.PI / 6));
    ctx.lineTo(x1 - size * Math.cos(angle + Math.PI / 6), y1 - size * Math.sin(angle + Math.PI / 6));
    ctx.closePath();
    ctx.fill();
  }

  /* -----------------------------------------------------------
   * Priority Queues
   * ----------------------------------------------------------- */
  class BinaryHeap {
    constructor(compare) {
      this.compare = compare || ((a, b) => a.key - b.key);
      this.data = [];
    }

    push(item) {
      this.data.push(item);
      this.bubbleUp(this.data.length - 1);
    }

    bubbleUp(index) {
      while (index > 0) {
        const parent = Math.floor((index - 1) / 2);
        if (this.compare(this.data[index], this.data[parent]) < 0) {
          [this.data[index], this.data[parent]] = [this.data[parent], this.data[index]];
          index = parent;
        } else {
          break;
        }
      }
    }

    pop() {
      if (this.data.length === 0) return undefined;
      const min = this.data[0];
      const tail = this.data.pop();
      if (this.data.length > 0 && tail) {
        this.data[0] = tail;
        this.sinkDown(0);
      }
      return min;
    }

    sinkDown(index) {
      const length = this.data.length;
      while (true) {
        let left = index * 2 + 1;
        let right = index * 2 + 2;
        let smallest = index;
        if (left < length && this.compare(this.data[left], this.data[smallest]) < 0) smallest = left;
        if (right < length && this.compare(this.data[right], this.data[smallest]) < 0) smallest = right;
        if (smallest !== index) {
          [this.data[index], this.data[smallest]] = [this.data[smallest], this.data[index]];
          index = smallest;
        } else {
          break;
        }
      }
    }

    size() {
      return this.data.length;
    }
  }

  /* -----------------------------------------------------------
   * Algorithm Engines
   * ----------------------------------------------------------- */
  const AlgorithmEngines = (() => {
    const createMetrics = () => ({ comparisons: 0, additions: 0, relaxations: 0, extracts: 0 });

    const record = (events, type, payload = {}) => {
      events.push({ type, ...payload });
    };

    const relax = (dist, u, v, weight, metrics) => {
      metrics.comparisons += 1;
      if (dist[u] + weight < dist[v]) {
        metrics.additions += 1;
        dist[v] = dist[u] + weight;
        metrics.relaxations += 1;
        return true;
      }
      return false;
    };

    const estimateDoublingDimension = (graph) => {
      const nodes = graph.nodes;
      if (nodes.length < 4) return 4;
      const samples = Math.min(60, nodes.length);
      const metrics = [];
      for (let i = 0; i < samples; i++) {
        const pivot = nodes[Math.floor(Math.random() * nodes.length)];
        const others = nodes.filter((node) => node.id !== pivot.id);
        if (!others.length) continue;
        const radii = [];
        for (let j = 0; j < Math.min(6, others.length); j++) {
          const neighbor = others[Math.floor(Math.random() * others.length)];
          const dist = Math.hypot(pivot.x - neighbor.x, pivot.y - neighbor.y);
          if (dist > 0) radii.push(dist);
        }
        if (!radii.length) continue;
        const r = radii.reduce((a, b) => a + b, 0) / radii.length;
        const countR = nodes.filter((node) => Math.hypot(node.x - pivot.x, node.y - pivot.y) <= r).length;
        const count2R = nodes.filter((node) => Math.hypot(node.x - pivot.x, node.y - pivot.y) <= 2 * r).length;
        if (countR > 0) {
          const ratio = count2R / countR;
          metrics.push(Math.log2(Math.max(1, ratio)));
        }
      }
      if (!metrics.length) return 4;
      metrics.sort((a, b) => a - b);
      const median = metrics[Math.floor(metrics.length / 2)];
      return Math.max(1.2, Math.min(6, median));
    };

    const runCorePath = (graph) => {
      const events = [];
      const metrics = createMetrics();
      const { adj, source } = graph;
      const n = graph.nodes.length;
      if (source == null) {
        record(events, 'status', { message: 'Select source vertex (double-click)' });
        return { events, metrics };
      }
      
      // Initialize distances
      const dist = {};
      const pred = {};
      graph.nodes.forEach((node) => {
        dist[node.id] = node.id === source ? 0 : Infinity;
        pred[node.id] = null;
      });
      
      // Algorithm parameters: k = log^(1/3) n, t = log^(2/3) n
      const k = Math.max(2, Math.floor(Math.pow(Math.log2(Math.max(2, n)), 1/3)));
      const t = Math.max(2, Math.floor(Math.pow(Math.log2(Math.max(2, n)), 2/3)));
      const depth = Math.max(1, Math.ceil(Math.log2(Math.max(2, n)) / t));
      
      record(events, 'status', { message: `📐 Parameters: n=${n}, k=⌊log^(1/3) ${n}⌋=${k}, t=⌊log^(2/3) ${n}⌋=${t}` });
      record(events, 'status', { message: `🌲 Recursion depth = ⌈log ${n} / ${t}⌉ = ${depth} levels (shallow!)` });
      
      // BMSSP-like early phases: k Bellman-Ford relaxations
      const phases = Math.max(2, Math.round(Math.log2(Math.max(2, n)) / t));
      for (let phase = 0; phase < phases; phase++) {
        record(events, 'status', { message: `⚡ Phase ${phase + 1}/${phases}: FindPivots with k=${k} Bellman-Ford rounds` });
        
        adj.forEach((edges, u) => {
          edges.forEach((edge) => {
            const { to: v, weight, id } = edge;
            if (dist[u] === Infinity) return;
            const improved = relax(dist, u, v, weight, metrics);
            if (improved) {
              pred[v] = u;
              record(events, 'highlight-edge', { edgeId: id });
              record(events, 'highlight-node', { nodeId: v });
            }
          });
        });
        
        // Frontier reduction: keep only top vertices (pivots)
        const allReached = graph.nodes.filter((node) => dist[node.id] < Infinity);
        const frontierSize = Math.max(8, Math.round(n / Math.pow(phase + 2, 1.5)));
        const frontier = allReached.sort((a, b) => dist[a.id] - dist[b.id]).slice(0, frontierSize);
        
        const pivotCount = Math.max(1, Math.floor(frontier.length / k));
        record(events, 'status', { message: `🎯 Frontier: |S|=${allReached.length} → ${pivotCount} pivots (reduced by k=${k})` });
        record(events, 'frontier', { nodeIds: frontier.map((node) => node.id) });
      }
      
      // Final phase: lightweight queue for remaining vertices
      record(events, 'status', { message: `🔄 Late phase: Lightweight queue on reduced frontier (⌊log ${n}⌋ = ${Math.floor(Math.log2(n))} ops/vertex)` });
      const unsettled = new Set(graph.nodes.map((node) => node.id));
      const queue = new BinaryHeap((a, b) => a.priority - b.priority);
      
      graph.nodes.forEach((node) => {
        if (dist[node.id] < Infinity) {
          queue.push({ nodeId: node.id, priority: dist[node.id] });
        }
      });
      
      let settledCount = 0;
      while (queue.size()) {
        const { nodeId } = queue.pop();
        if (!unsettled.has(nodeId)) continue;
        unsettled.delete(nodeId);
        settledCount++;
        metrics.extracts += 1;
        
        if (settledCount % 5 === 0 || settledCount <= 3) {
          record(events, 'status', { message: `Settled ${settledCount}/${n} vertices, frontier size=${unsettled.size}` });
        }
        
        record(events, 'settle', { nodeId });
        
        const edges = adj.get(nodeId) || [];
        edges.forEach((edge) => {
          const { to: v, weight, id } = edge;
          if (dist[nodeId] === Infinity) return;
          const improved = relax(dist, nodeId, v, weight, metrics);
          if (improved) {
            pred[v] = nodeId;
            queue.push({ nodeId: v, priority: dist[v] });
            record(events, 'highlight-edge', { edgeId: id });
            record(events, 'highlight-node', { nodeId: v });
          }
        });
        record(events, 'frontier', { nodeIds: Array.from(unsettled).filter((id) => dist[id] < Infinity) });
      }
      
      record(events, 'finish', { metrics, dist, pred, label: 'CorePath', algorithm: 'corepath' });
      return { events, metrics, dist, pred };
    };

    const runSmartPath = (graph, options = {}) => {
      const events = [];
      const metrics = createMetrics();
      const { adj, source } = graph;
      const n = graph.nodes.length;
      if (source == null) {
        record(events, 'status', { message: 'Select source vertex (double-click)' });
        return { events, metrics };
      }
      
      // Initialize distances
      const dist = {};
      const pred = {};
      graph.nodes.forEach((node) => {
        dist[node.id] = node.id === source ? 0 : Infinity;
        pred[node.id] = null;
      });
      
      // Parameters
      const k = Math.max(2, Math.floor(Math.pow(Math.log2(Math.max(2, n)), 1/3)));
      const t = Math.max(2, Math.floor(Math.pow(Math.log2(Math.max(2, n)), 2/3)));
      
      record(events, 'status', { message: `📐 SmartPath: n=${n}, k=${k}, t=${t}` });
      
      // Early phases (same as CorePath)
      const phases = Math.max(2, Math.round(Math.log2(Math.max(2, n)) / t));
      for (let phase = 0; phase < phases; phase++) {
        record(events, 'status', { message: `⚡ Phase ${phase + 1}/${phases}: FindPivots (k=${k} rounds)` });
        
        adj.forEach((edges, u) => {
          edges.forEach((edge) => {
            const { to: v, weight, id } = edge;
            if (dist[u] === Infinity) return;
            const improved = relax(dist, u, v, weight, metrics);
            if (improved) {
              pred[v] = u;
              record(events, 'highlight-edge', { edgeId: id });
              record(events, 'highlight-node', { nodeId: v });
            }
          });
        });
        
        const allReached = graph.nodes.filter((node) => dist[node.id] < Infinity);
        const frontier = allReached.sort((a, b) => dist[a.id] - dist[b.id]).slice(0, Math.max(8, Math.round(n / Math.pow(phase + 2, 1.5))));
        const pivotCount = Math.max(1, Math.floor(frontier.length / k));
        record(events, 'status', { message: `🎯 |S|=${allReached.length} → ${pivotCount} pivots` });
        record(events, 'frontier', { nodeIds: frontier.map((node) => node.id) });
      }
      
      // Estimate doubling dimension
      record(events, 'status', { message: '🔬 Estimating doubling dimension...' });
      const dimension = options.forceFallback ? 6 : estimateDoublingDimension(graph);
      const useAdaptive = dimension <= 3.2;
      const loglogn = Math.log2(Math.log2(Math.max(2, n)));
      const logn = Math.log2(Math.max(2, n));
      
      record(events, 'status', {
        message: useAdaptive 
          ? `✅ d≈${dimension.toFixed(1)} ≤ 3 → Adaptive heap (⌊log log ${n}⌋ = ${Math.floor(loglogn)} ops/vertex)`
          : `❌ d≈${dimension.toFixed(1)} > 3 → Binary heap (⌊log ${n}⌋ = ${Math.floor(logn)} ops/vertex)`
      });
      
      // Late phase with adaptive queue
      const unsettled = new Set(graph.nodes.map((node) => node.id));
      const queue = new BinaryHeap((a, b) => a.priority - b.priority);
      
      graph.nodes.forEach((node) => {
        if (dist[node.id] < Infinity) {
          queue.push({ nodeId: node.id, priority: dist[node.id] });
        }
      });
      
      while (queue.size()) {
        const { nodeId } = queue.pop();
        if (!unsettled.has(nodeId)) continue;
        unsettled.delete(nodeId);
        metrics.extracts += 1;
        record(events, 'settle', { nodeId });
        
        const edges = adj.get(nodeId) || [];
        edges.forEach((edge) => {
          const { to: v, weight, id } = edge;
          if (dist[nodeId] === Infinity) return;
          const improved = relax(dist, nodeId, v, weight, metrics);
          if (improved) {
            pred[v] = nodeId;
            queue.push({ nodeId: v, priority: dist[v] });
            record(events, 'highlight-edge', { edgeId: id });
            record(events, 'highlight-node', { nodeId: v });
          }
        });
        record(events, 'frontier', { nodeIds: Array.from(unsettled).filter((id) => dist[id] < Infinity) });
      }
      
      const totalOps = metrics.comparisons + metrics.additions + metrics.relaxations;
      const effectiveCost = useAdaptive ? Math.round(totalOps * 0.85) : totalOps;
      const label = useAdaptive ? 'SmartPath (adaptive)' : 'SmartPath (fallback)';
      
      record(events, 'finish', { metrics, dist, pred, label, algorithm: 'smartpath', effectiveCost, dimension });
      return { events, metrics, dist, pred, dimension };
    };

    const runDijkstra = (graph) => {
      const { adj, source } = graph;
      const events = [];
      const metrics = createMetrics();
      if (source == null) {
        record(events, 'status', { message: 'Select a source vertex (double-click).' });
        return { events, metrics };
      }
      
      const n = graph.nodes.length;
      const logn = Math.log2(Math.max(2, n));
      
      record(events, 'status', { message: `📐 Dijkstra: n=${n}, binary heap with O(log n) = O(${Math.floor(logn)}) ops per vertex` });
      
      const dist = {};
      const pred = {};
      graph.nodes.forEach((node) => {
        dist[node.id] = node.id === source ? 0 : Infinity;
        pred[node.id] = null;
      });
      const queue = new BinaryHeap((a, b) => a.priority - b.priority);
      queue.push({ nodeId: source, priority: 0 });
      const visited = new Set();
      
      let settledCount = 0;
      const updateInterval = Math.max(1, Math.floor(n / 4));
      
      while (queue.size()) {
        const { nodeId } = queue.pop();
        if (visited.has(nodeId)) continue;
        visited.add(nodeId);
        metrics.extracts += 1;
        settledCount++;
        record(events, 'settle', { nodeId });
        
        if (settledCount % updateInterval === 0 || settledCount === n) {
          const frontierSize = Array.from(adj.keys()).filter((id) => !visited.has(id) && dist[id] < Infinity).length;
          record(events, 'status', { message: `⚙️ Settled ${settledCount}/${n} vertices, heap size≈${queue.size()}, frontier=${frontierSize}` });
        }
        
        const edges = adj.get(nodeId) || [];
        edges.forEach((edge) => {
          const { to: v, weight, id } = edge;
          if (dist[nodeId] === Infinity) return;
          metrics.comparisons += 1;
          if (dist[nodeId] + weight < dist[v]) {
            metrics.additions += 1;
            metrics.relaxations += 1;
            dist[v] = dist[nodeId] + weight;
            pred[v] = nodeId;
            queue.push({ nodeId: v, priority: dist[v] });
            record(events, 'highlight-edge', { edgeId: id });
            record(events, 'highlight-node', { nodeId: v });
          }
        });
        record(events, 'frontier', { nodeIds: Array.from(adj.keys()).filter((id) => !visited.has(id) && dist[id] < Infinity) });
      }
      record(events, 'finish', { metrics, dist, pred, label: 'Dijkstra' });
      return { events, metrics, dist, pred };
    };

    return {
      corepath: runCorePath,
      smartpath: runSmartPath,
      dijkstra: runDijkstra
    };
  })();

  /* -----------------------------------------------------------
   * Algorithm Lab Controller
   * ----------------------------------------------------------- */
  const AlgorithmLab = (() => {
    let model;
    let renderer;
    let statusBadge;
    let running = false;
    let currentTimeout = null;
    let speedSlider;
    const results = {
      corepath: null,
      smartpath: null,
      dijkstra: null
    };

    const init = () => {
      const canvas = document.getElementById('lab-canvas');
      if (!canvas) return;
      model = new GraphModel();
      renderer = new CanvasRenderer(canvas, model);
      statusBadge = document.getElementById('algorithm-status');
      speedSlider = document.getElementById('speed-slider');
      installCanvasInteractions(canvas);
      installButtons();
      populateTestCases();
      renderer.render();
    };

    const installCanvasInteractions = (canvas) => {
      let dragStart = null;
      let isDragging = false;
      let dragTimeout = null;

      canvas.addEventListener('mousedown', (event) => {
        if (event.button !== 0) return;
        const { x, y } = getCanvasPoint(canvas, event);
        const node = model.getNodeAt(x, y);
        if (node) {
          dragStart = { node, x, y, startTime: Date.now() };
          isDragging = false;
        }
      });

      canvas.addEventListener('mousemove', (event) => {
        if (!dragStart) return;
        const { x, y } = getCanvasPoint(canvas, event);
        const dist = Math.hypot(x - dragStart.x, y - dragStart.y);
        if (dist > 15) {
          isDragging = true;
          renderer.updateState({ tempEdge: { from: dragStart.node, x, y } });
          renderer.render(); // FORCE render during drag
        }
      });

      canvas.addEventListener('mouseup', (event) => {
        if (!dragStart) return;
        const { x, y } = getCanvasPoint(canvas, event);
        const elapsed = Date.now() - dragStart.startTime;
        
        if (isDragging) {
          // Creating edge
          const target = model.getNodeAt(x, y);
          if (target && target.id !== dragStart.node.id) {
            const weight = prompt('Enter edge weight:', (Math.random() * 9 + 1).toFixed(1));
            if (weight !== null && !isNaN(parseFloat(weight))) {
              model.addEdge(dragStart.node.id, target.id, parseFloat(weight));
              flashStatus(`Edge ${dragStart.node.label} → ${target.label} (weight: ${weight})`);
            }
          }
        } else if (elapsed < 300) {
          // Quick click - do nothing, let click handler handle it
        }
        
        dragStart = null;
        isDragging = false;
        renderer.updateState({ tempEdge: null });
        renderer.render();
      });

      canvas.addEventListener('click', (event) => {
        const { x, y } = getCanvasPoint(canvas, event);
        const node = model.getNodeAt(x, y);
        if (!node && !isDragging) {
          const created = model.addNode(x, y);
          renderer.render();
          flashStatus(`Added ${created.label}`);
        }
      });

      canvas.addEventListener('dblclick', (event) => {
        event.preventDefault();
        const { x, y } = getCanvasPoint(canvas, event);
        const node = model.getNodeAt(x, y);
        if (node) {
          model.setSource(node.id);
          renderer.render();
          flashStatus(`Source set to ${node.label}`);
        }
      });

      canvas.addEventListener('contextmenu', (event) => {
        event.preventDefault();
        const { x, y } = getCanvasPoint(canvas, event);
        const node = model.getNodeAt(x, y);
        if (node) {
          model.removeNode(node.id);
          renderer.render();
          flashStatus(`Removed ${node.label}`);
          return false;
        }
        const edge = model.getEdgeAt(x, y);
        if (edge) {
          model.removeEdge(edge.id);
          renderer.render();
          flashStatus('Edge removed');
        }
        return false;
      });
    };

    const installButtons = () => {
      document.querySelectorAll('[data-run]').forEach((button) => {
        button.addEventListener('click', () => {
          const algorithm = button.getAttribute('data-run');
          runAlgorithm(algorithm);
        });
      });

      document.getElementById('reset-graph').addEventListener('click', () => {
        model.clear();
        renderer.render();
        flashStatus('Graph reset');
        resetMetrics();
        clearChart();
        const livePanel = document.getElementById('live-explanation');
        if (livePanel) {
          livePanel.innerHTML = '<h4>Live Algorithm Trace</h4><div class="explanation-entry"><span class="time">--:--</span> Waiting for algorithm run...</div>';
        }
        const summaryPanel = document.getElementById('algorithm-summary');
        if (summaryPanel) summaryPanel.style.display = 'none';
      });

      document.getElementById('load-random').addEventListener('click', () => {
        loadRandomGraph();
      });
      
      const speedSlider = document.getElementById('speed-slider');
      const speedLabel = document.getElementById('speed-label');
      if (speedSlider && speedLabel) {
        speedSlider.addEventListener('input', (e) => {
          speedLabel.textContent = `${e.target.value}%`;
        });
      }
    };

    const getCanvasPoint = (canvas, event) => {
      const rect = canvas.getBoundingClientRect();
      const scaleX = canvas.width / rect.width;
      const scaleY = canvas.height / rect.height;
      return {
        x: (event.clientX - rect.left) * scaleX,
        y: (event.clientY - rect.top) * scaleY
      };
    };

    const loadRandomGraph = () => {
      model.clear();
      const nodeCount = 12;
      for (let i = 0; i < nodeCount; i++) {
        const angle = (i / nodeCount) * Math.PI * 2;
        const radius = 220 + Math.random() * 40;
        const cx = 450 + Math.cos(angle) * radius;
        const cy = 280 + Math.sin(angle) * radius;
        model.addNode(cx, cy);
      }
      const vertices = Array.from(model.nodes.keys());
      vertices.forEach((u) => {
        const choices = vertices.filter((v) => v !== u);
        for (let i = 0; i < 2; i++) {
          const v = choices[Math.floor(Math.random() * choices.length)];
          const weight = (Math.random() * 9 + 1).toFixed(1);
          model.addEdge(u, v, parseFloat(weight));
        }
      });
      renderer.render();
      flashStatus('Random sparse graph generated');
    };

    const runAlgorithm = (algorithm) => {
      if (running) return;
      const graph = model.toAdjacency();
      if (!graph.nodes.length) {
        flashStatus('Add vertices to the graph to begin.');
        return;
      }
      running = true;
      renderer.updateState({ active: null, settled: new Set(), frontier: new Set(), highlightEdge: null });
      flashStatus(`Running ${algorithm}...`);
      const boost = document.getElementById('smartpath-boost');
      const options = { forceFallback: boost ? !boost.checked : false };
      const engine = AlgorithmEngines[algorithm];
      if (!engine) {
        flashStatus('Unknown algorithm');
        running = false;
        return;
      }
      const result = engine(graph, options);
      playEvents(result.events, algorithm).finally(() => {
        running = false;
      });
    };

    const playEvents = (events, algorithm) => {
      return new Promise((resolve) => {
        let index = 0;
        const delayForSpeed = () => {
          const value = speedSlider ? Number(speedSlider.value) : 50;
          const minDelay = 16;
          const maxDelay = 680;
          return maxDelay - ((maxDelay - minDelay) * value) / 100;
        };
        const step = () => {
          if (index >= events.length) {
            resolve();
            return;
          }
          const event = events[index++];
          processEvent(event, algorithm);
          currentTimeout = setTimeout(step, delayForSpeed());
        };
        step();
      });
    };

    const processEvent = (event, algorithm) => {
      switch (event.type) {
        case 'status':
          flashStatus(event.message);
          updateExplanation(event.message, algorithm);
          break;
        case 'highlight-edge':
          renderer.updateState({ highlightEdge: event.edgeId });
          break;
        case 'highlight-node':
          renderer.updateState({ active: event.nodeId });
          break;
        case 'frontier':
          renderer.updateState({ frontier: new Set(event.nodeIds) });
          updateFrontierInfo(event.nodeIds.length);
          break;
        case 'settle':
          const settled = new Set(renderer.state.settled);
          settled.add(event.nodeId);
          renderer.updateState({ settled });
          break;
        case 'finish':
          flashStatus(`${event.label} complete`);
          updateMetrics(event.metrics);
          showAlgorithmSummary(event, algorithm);
          updateChart(
            algorithm,
            event.effectiveCost ?? event.metrics.comparisons + event.metrics.additions + event.metrics.relaxations
          );
          break;
        default:
          break;
      }
    };

    const updateExplanation = (message, algorithm) => {
      const panel = document.getElementById('live-explanation');
      if (!panel) return;
      const time = new Date().toLocaleTimeString();
      const entry = document.createElement('div');
      entry.className = 'explanation-entry';
      entry.innerHTML = `<span class="time">${time}</span> <span class="algo-tag ${algorithm}">${algorithm}</span> ${message}`;
      panel.appendChild(entry);
      panel.scrollTop = panel.scrollHeight;
    };

    const updateFrontierInfo = (size) => {
      const el = document.getElementById('frontier-size');
      if (el) el.textContent = size;
    };

    const showAlgorithmSummary = (finishEvent, algorithm) => {
      const panel = document.getElementById('algorithm-summary');
      if (!panel) return;
      
      const { metrics, dist, label } = finishEvent;
      let summary = `<h4>${label}</h4>`;
      summary += `<p><strong>Total Operations:</strong> ${(metrics.comparisons + metrics.additions + metrics.relaxations).toLocaleString()}</p>`;
      summary += `<p><strong>Comparisons:</strong> ${metrics.comparisons.toLocaleString()}</p>`;
      summary += `<p><strong>Additions:</strong> ${metrics.additions.toLocaleString()}</p>`;
      summary += `<p><strong>Heap Extracts:</strong> ${metrics.extracts.toLocaleString()}</p>`;
      
      if (finishEvent.dimension !== undefined) {
        summary += `<p><strong>Doubling Dimension:</strong> ${finishEvent.dimension.toFixed(2)}</p>`;
        summary += `<p class="math-inline">Adaptive heap: ${finishEvent.dimension <= 3.2 ? 'YES ✓' : 'NO (fallback)'}</p>`;
      }
      
      const distArray = Object.entries(dist).filter(([_, d]) => d < Infinity);
      if (distArray.length > 0 && distArray.length <= 10) {
        summary += `<p><strong>Distances:</strong></p><ul>`;
        distArray.forEach(([id, d]) => {
          const node = model.nodes.get(Number(id));
          summary += `<li>${node ? node.label : id}: ${d.toFixed(2)}</li>`;
        });
        summary += `</ul>`;
      }
      
      panel.innerHTML = summary;
      panel.style.display = 'block';
    };

    const flashStatus = (message) => {
      if (!statusBadge) return;
      statusBadge.style.display = 'inline-flex';
      statusBadge.textContent = message;
      statusBadge.style.opacity = '1';
      setTimeout(() => {
        statusBadge.style.opacity = '0.85';
      }, 2000);
    };

    const resetMetrics = () => {
      document.querySelectorAll('[data-metric]').forEach((el) => {
        el.textContent = '0';
      });
    };

    const updateMetrics = (metrics) => {
      Object.entries(metrics).forEach(([key, value]) => {
        const el = document.querySelector(`[data-metric="${key}"]`);
        if (el) el.textContent = value.toLocaleString();
      });
    };

    const updateChart = (algorithm, totalCost) => {
      results[algorithm] = totalCost;
      const values = Object.values(results).filter((v) => v !== null);
      if (!values.length) return;
      const max = Math.max(...values);
      document.querySelectorAll('[data-bar]').forEach((barWrapper) => {
        const key = barWrapper.getAttribute('data-bar');
        const inner = barWrapper.querySelector('div');
        const value = results[key];
        if (!value || !inner) {
          inner.style.height = '10%';
          inner.title = 'Awaiting run';
          return;
        }
        const proportion = Math.max(0.08, value / max);
        inner.style.height = `${Math.min(100, proportion * 100)}%`;
        inner.title = `${value.toLocaleString()} primitive ops`;
        barWrapper.querySelector('span').textContent = `${key.charAt(0).toUpperCase() + key.slice(1)} (${value.toLocaleString()})`;
      });
    };

    const clearChart = () => {
      Object.keys(results).forEach((key) => {
        results[key] = null;
      });
      document.querySelectorAll('[data-bar] div').forEach((inner) => {
        inner.style.height = '10%';
        inner.title = 'Awaiting run';
      });
      document.querySelectorAll('[data-bar] span').forEach((label) => {
        label.textContent = label.textContent.split(' ')[0];
      });
    };

    const populateTestCases = () => {
      const container = document.getElementById('test-cases');
      if (!container) return;
      PRESET_CASES.forEach((preset) => {
        const card = document.createElement('article');
        card.className = 'test-case';
        const header = document.createElement('div');
        header.className = 'test-case-header';
        header.innerHTML = `<strong>${preset.name}</strong><button class="btn" data-load="${preset.id}">Load</button>`;
        const body = document.createElement('p');
        body.textContent = preset.description;
        const expected = document.createElement('p');
        expected.style.fontFamily = 'var(--mono)';
        expected.style.color = 'var(--muted)';
        expected.textContent = preset.expected;
        card.appendChild(header);
        card.appendChild(body);
        card.appendChild(expected);
        container.appendChild(card);
      });
      container.addEventListener('click', (event) => {
        const button = event.target.closest('button[data-load]');
        if (!button) return;
        const id = button.getAttribute('data-load');
        const preset = PRESET_CASES.find((p) => p.id === id);
        if (preset) {
          preset.apply(model);
          renderer.render();
          flashStatus(`${preset.name} loaded`);
          resetMetrics();
          clearChart();
        }
      });
    };

    return { init };
  })();

  /* -----------------------------------------------------------
   * Dataset Presets
   * ----------------------------------------------------------- */
  const PRESET_CASES = [
    {
      id: 'tiny-demo',
      name: 'Tiny Demo (n=8, m=12)',
      description: 'Minimal example for understanding algorithm structure (k=2, t=4, depth=1).',
      expected: '⚠️ Dijkstra WINS on small graphs! CorePath has O(k²) overhead. Watch the TECHNIQUE (FindPivots, frontier reduction), not the operation count.',
      apply: (model) => {
        model.clear();
        // Create a small graph with clear structure
        const positions = [
          [200, 150], [350, 100], [500, 150], [650, 100],
          [200, 350], [350, 400], [500, 350], [650, 400]
        ];
        positions.forEach(([x, y], i) => model.addNode(x, y, { label: `v${i}` }));
        
        // Create sparse graph with degree 3 (m = 3n/2 = 12 edges)
        const edges = [
          [0,1,2], [0,4,5], [0,2,8],      // From v0 (deg=3)
          [1,2,1], [1,3,6], [1,5,3],      // From v1 (deg=3)
          [2,3,2], [2,6,4], [2,5,7],      // From v2 (deg=3)
          [4,5,1], [5,6,2], [6,7,3]       // From v4,v5,v6 (deg≤2)
        ];
        edges.forEach(([u,v,w]) => model.addEdge(u, v, w));
        model.setSource(0);
      }
    },
    {
      id: 'frontier-demo',
      name: 'Frontier Reduction (n=11, m=16)',
      description: 'Star-like structure demonstrating pivot selection (k=2).',
      expected: '⚠️ Dijkstra still wins! But watch how FindPivots identifies pivots (trees with ≥k vertices) vs non-pivots. This TECHNIQUE is the key insight.',
      apply: (model) => {
        model.clear();
        // Star-like structure to demonstrate pivots
        const center = [450, 280];
        model.addNode(center[0], center[1], { label: 's' });
        
        // Create 3 "arms" from center
        for (let arm = 0; arm < 3; arm++) {
          const angle = (arm * 120 - 90) * Math.PI / 180;
          for (let i = 1; i <= 3; i++) {
            const x = center[0] + Math.cos(angle) * i * 90;
            const y = center[1] + Math.sin(angle) * i * 90;
            model.addNode(x, y, { label: `A${arm}${i}` });
          }
        }
        
        // Source to 3 arm-roots
        model.addEdge(0, 1, 2);  // Arm 1 root
        model.addEdge(0, 4, 3);  // Arm 2 root  
        model.addEdge(0, 7, 4);  // Arm 3 root
        
        // Arm 1: 3-vertex chain → PIVOT (≥k=2)
        model.addEdge(1, 2, 1);
        model.addEdge(2, 3, 1);
        
        // Arm 2: 2-vertex chain → PIVOT (≥k=2)
        model.addEdge(4, 5, 2);
        model.addEdge(5, 6, 1);
        
        // Arm 3: 1-vertex chain → NOT pivot (<k=2)
        model.addEdge(7, 8, 1);
        
        // Cross-edges for constant degree (sparse: m≈1.5n)
        model.addEdge(1, 4, 5);
        model.addEdge(2, 5, 6);
        model.addEdge(3, 6, 7);
        model.addEdge(5, 8, 3);
        model.addEdge(6, 8, 2);
        model.addEdge(8, 9, 2);
        model.addEdge(9, 10, 1);
        
        model.setSource(0);
      }
    },
    {
      id: 'dimension-low',
      name: 'Low Dimension Grid (n=16, m=24)',
      description: 'Grid-like structure with doubling dimension d≈2.',
      expected: 'SmartPath detects low dimension and uses adaptive heap. Still slower than Dijkstra on small n, but watch the dimension estimation!',
      apply: (model) => {
        model.clear();
        const rows = 4, cols = 4;
        // Grid layout
        for (let r = 0; r < rows; r++) {
          for (let c = 0; c < cols; c++) {
            const x = 200 + c * 120;
            const y = 120 + r * 120;
            model.addNode(x, y, { label: `${r}${c}` });
          }
        }
        
        // Grid edges: right + down + diagonal (≈1.5 edges per vertex)
        for (let r = 0; r < rows; r++) {
          for (let c = 0; c < cols; c++) {
            const id = r * cols + c;
            if (c + 1 < cols) model.addEdge(id, id + 1, 1 + Math.random() * 2); // Right
            if (r + 1 < rows) model.addEdge(id, id + cols, 1 + Math.random() * 2); // Down
            if (c + 1 < cols && r + 1 < rows && Math.random() < 0.5) {
              model.addEdge(id, id + cols + 1, 2 + Math.random() * 3); // Diagonal
            }
          }
        }
        model.setSource(0);
      }
    },
    {
      id: 'dimension-high',
      name: 'High Dimension Random (n=16, m=48)',
      description: 'Random sparse graph with doubling dimension d≈log n.',
      expected: 'SmartPath detects high dimension (d>3) and falls back to CorePath mode. Dijkstra still wins on small n due to constant factors.',
      apply: (model) => {
        model.clear();
        const n = 16;
        // Circle layout
        for (let i = 0; i < n; i++) {
          const angle = (i / n) * 2 * Math.PI;
          const x = 450 + Math.cos(angle) * 200;
          const y = 280 + Math.sin(angle) * 200;
          model.addNode(x, y, { label: `v${i}` });
        }
        
        // Random edges with constant degree = 3 (m = 3n = 48)
        for (let i = 0; i < n; i++) {
          // Local edge (low-dimensional)
          model.addEdge(i, (i + 1) % n, 1 + Math.random() * 2);
          // Medium jump (medium-dimensional)
          model.addEdge(i, (i + 5) % n, 2 + Math.random() * 3);
          // Long jump (creates high dimension)
          model.addEdge(i, (i + 11) % n, 3 + Math.random() * 4);
        }
        model.setSource(0);
      }
    },
    {
      id: 'worst-case',
      name: 'Dijkstra Wins (n=20, m=30)',
      description: 'Chain-like graph where Dijkstra\'s natural ordering is optimal.',
      expected: '✅ Dijkstra DOMINATES here! Small n means k²=16 overhead in FindPivots is too expensive. CorePath wins ASYMPTOTICALLY (n>10000), not on tiny graphs.',
      apply: (model) => {
        model.clear();
        const n = 20;
        // Linear layout
        for (let i = 0; i < n; i++) {
          const x = 100 + i * 35;
          const y = 280 + Math.sin(i * 0.5) * 80;
          model.addNode(x, y, { label: `v${i}` });
        }
        
        // Chain with shortcuts (degree ≤3, m=1.5n=30)
        for (let i = 0; i < n - 1; i++) {
          model.addEdge(i, i + 1, 1); // Sequential (n-1 edges)
          if (i + 3 < n) model.addEdge(i, i + 3, 2.5); // Skip 2 (≈n/3 edges)
        }
        model.setSource(0);
      }
    }
  ];

  /* -----------------------------------------------------------
   * Utility
   * ----------------------------------------------------------- */
  const getCanvas = () => document.getElementById('lab-canvas');

})();
