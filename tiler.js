/* Aperiodic Tiler — Cut & Project (with Thomas’s fixes)
 * numeric.js required (loaded in index.html)
 * All logic is client-side.
 */

// ---------- Utilities ----------
function dot(a, b) { let s = 0; for (let i = 0; i < a.length; i++) s += a[i] * b[i]; return s; }
function sub(a, b) { return a.map((v, i) => v - b[i]); }
function scale(a, s) { return a.map(v => v * s); }
function norm(a) { return Math.sqrt(dot(a, a)); }
function inv(M) { return numeric.inv(M); }
function mulMV(M, v) { return numeric.dot(M, v); }

function combinations(n, k) {
  const res = [];
  if (k === 0) return [[]];
  if (k > n) return res;
  const comb = Array.from({ length: k }, (_, i) => i);
  const next = () => {
    let i = k - 1;
    while (i >= 0 && comb[i] === n - k + i) i--;
    if (i < 0) return false;
    comb[i]++;
    for (let j = i + 1; j < k; j++) comb[j] = comb[j - 1] + 1;
    return true;
  };
  res.push(comb.slice());
  while (next()) res.push(comb.slice());
  return res;
}

function cartesianRange(d, k) {
  const vals = [];
  const range = [];
  for (let i = -k; i <= k; i++) range.push(i);
  (function rec(level, acc) {
    if (level === d) { vals.push(acc.slice()); return; }
    for (const v of range) { acc.push(v); rec(level + 1, acc); acc.pop(); }
  })(0, []);
  return vals;
}

// Gram–Schmidt on ROWS (returns orthonormal rows)
function gramSchmidtRows(E) {
  const d = E.length, n = E[0].length;
  const U = Array.from({ length: d }, () => Array(n).fill(0));
  for (let i = 0; i < d; i++) {
    let v = E[i].slice();
    for (let j = 0; j < i; j++) {
      const proj = dot(v, U[j]);
      v = sub(v, scale(U[j], proj));
    }
    const nv = norm(v);
    if (nv < 1e-12) throw new Error('Gram–Schmidt failed: dependent rows.');
    U[i] = scale(v, 1 / nv);
  }
  return U;
}

// ---------- Core logic (ports of SageMath) ----------

// generators_to_grid(E) — E: d×n -> G: n×d
function generators_to_grid(E) {
  const d = E.length, n = E[0].length;
  const U = gramSchmidtRows(E); // d×n orthonormal
  const G = Array.from({ length: n }, (_, i) =>
    Array.from({ length: d }, (_, j) => U[j][i])
  );
  return G;
}

// dual(G, S, k)
function dual(G, S = null, k = 10) {
  const n = G.length;
  const d = G[0].length;
  if (!S) S = Array.from({ length: n }, () => Math.random());
  if (S.length !== n) throw new Error('Shift S must have length n.');

  const tiles = [];
  const combos = combinations(n, d);
  const indices = cartesianRange(d, k);

  for (const t of combos) {
    const Grows = t.map(idx => G[idx]); // d×d
    let M;
    try { M = inv(Grows); } catch { continue; }
    for (const ind of indices) {
      const rhs = ind.map((v, j) => S[t[j]] + v);
      const intersec = mulMV(M, rhs); // length d
      const pos = Array(n).fill(0).map((_, j) => {
        const val = dot(intersec, G[j]) - S[j];
        return Math.ceil(val - 1e-12);
      });
      if (pos.reduce((m, v) => Math.max(m, Math.abs(v)), 0) > k) continue;
      for (let h = 0; h < d; h++) pos[t[h]] = ind[h];
      tiles.push({ t: t.slice(), pos });
    }
  }
  return tiles;
}

function orthogonal_projection(E) {
  const d = E.length;
  if (d !== 2) throw new Error('orthogonal_projection expects d=2 for visualization.');
  return gramSchmidtRows(E);
}

// grayscale color like Sage's couleur(i,j,n)
function couleur(i, j, n) {
  const z = (n * i - (i * (i + 1)) / 2 + (j - i) - 1) / (n * (n - 1) / 2);
  const v = Math.floor(Math.max(0, Math.min(255, Math.floor(z * 156 + 100))));
  const hex = v.toString(16).padStart(2, '0');
  return `#${hex}${hex}${hex}`;
}

// draw_tiling — SVG
function draw_tiling(T, A, container) {
  const width = container.clientWidth || 800;
  const height = container.clientHeight || 600;
  const padding = 12;

  const Ax = A[0], Ay = A[1];
  const polys = [];
  let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;

  for (const tile of T) {
    const [i, j] = tile.t;
    const base = tile.pos.slice();
    const corners = [[0,0],[0,1],[1,1],[1,0]];
    const pts = [];
    for (const [ai, aj] of corners) {
      const q = base.slice();
      q[i] += ai; q[j] += aj;
      const x = dot(Ax, q), y = dot(Ay, q);
      pts.push([x, y]);
      if (x < minX) minX = x;
      if (x > maxX) maxX = x;
      if (y < minY) minY = y;
      if (y > maxY) maxY = y;
    }
    polys.push({ points: pts, color: couleur(i, j, Ax.length) });
  }

  if (polys.length === 0 || !isFinite(minX)) {
    container.innerHTML = `<div class="hint">No tiles to display. Try increasing <b>k</b>.</div>`;
    return;
  }

  const spanX = Math.max(1e-6, maxX - minX);
  const spanY = Math.max(1e-6, maxY - minY);
  const scale = 0.92 * Math.min((width - 2 * padding) / spanX, (height - 2 * padding) / spanY);

  const mapPoint = ([x, y]) => {
    const X = padding + (x - minX) * scale;
    const Y = height - padding - (y - minY) * scale;
    return `${X.toFixed(2)},${Y.toFixed(2)}`;
  };

  let svg = `<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}">`;
  svg += `<rect width="100%" height="100%" fill="white"/>`;
  for (const poly of polys) {
    const d = poly.points.map(mapPoint).join(' ');
    svg += `<polygon points="${d}" fill="${poly.color}" fill-opacity="0.95" stroke="none"/>`;
  }
  svg += `</svg>`;
  container.innerHTML = svg;
}

// ---------- Slope & Shift generation (Thomas’s rules) ----------

// For odd n: use n-th roots (angles 2πk/n).
// For even n: use 2n-th roots and take the first n (angles 2πk/(2n), k=0..n-1).
function generateSlopeE(n, d) {
  const m = (n % 2 === 0) ? 2 * n : n;  // modulus for roots
  const rows = Array.from({ length: d }, () => Array(n).fill(0));
  for (let k = 0; k < n; k++) {
    const th = (2 * Math.PI * k) / m;
    if (d >= 1) rows[0][k] = Math.cos(th);
    if (d >= 2) rows[1][k] = Math.sin(th);
    if (d >= 3) rows[2][k] = Math.cos(2 * th); // simple extra harmonic
  }
  return rows;
}

// Shift: never 0. Symmetric shift is 1/n (odd) or 1/(2n) (even). Otherwise random.
function makeShift(n, mode) {
  if (mode === 'symmetric') {
    const s = (n % 2 === 1) ? (1 / n) : (1 / (2 * n));
    return Array(n).fill(s);
  }
  // random
  return Array.from({ length: n }, () => Math.random());
}

// ---------- Wiring ----------
const els = {
  dDim: () => document.getElementById('dDim'),
  nCols: () => document.getElementById('nCols'),
  kBound: () => document.getElementById('kBound'),
  dLabel: () => document.getElementById('dLabel'),
  nLabel: () => document.getElementById('nLabel'),
  kLabel: () => document.getElementById('kLabel'),
  seedShift: () => document.getElementById('seedShift'),
  btnGenerate: () => document.getElementById('btnGenerate'),
  status: () => document.getElementById('status'),
  viz: () => document.getElementById('viz'),
  meta: () => document.getElementById('metaInfo'),
};

function setStatus(text) { els.status().textContent = text; }

function generate() {
  try {
    const d = parseInt(els.dDim().value, 10);
    const n = parseInt(els.nCols().value, 10);
    const k = parseInt(els.kBound().value, 10);
    const shiftMode = els.seedShift().value;

    els.dLabel().textContent = String(d);
    els.nLabel().textContent = String(n);
    els.kLabel().textContent = String(k);

    if (n < d) throw new Error('Require n ≥ d.');
    if (d !== 2) {
      els.viz().innerHTML = `<div class="hint">This demo visualizes only <b>d = 2</b>. Set d=2 and click “Generate”.</div>`;
      els.meta().textContent = `d=${d}, n=${n}, k=${k}`;
      setStatus('Ready (set d=2 to visualize).');
      return;
    }

    setStatus('Building slope E (roots rule)…');
    const E = generateSlopeE(n, d); // d×n, even n uses 2n-th roots

    setStatus('Computing multigrid G…');
    const G = generators_to_grid(E); // n×d

    setStatus('Generating dual tiling…');
    const S = makeShift(n, shiftMode); // never zero; symmetric or random
    const T = dual(G, S, k);

    setStatus('Building projection…');
    const A = orthogonal_projection(E); // 2×n

    setStatus('Drawing…');
    draw_tiling(T, A, els.viz());

    els.meta().textContent = `d=2, n=${n} (${n % 2 === 0 ? '2n-th roots' : 'n-th roots'}), tiles=${T.length}, k=${k}, shift=${shiftMode}`;
    setStatus('Done.');
  } catch (err) {
    console.error(err);
    els.viz().innerHTML = '';
    els.meta().textContent = '';
    setStatus(`Error: ${err.message || err}`);
  }
}

window.addEventListener('DOMContentLoaded', () => {
  const updateLabels = () => {
    els.dLabel().textContent = els.dDim().value;
    els.nLabel().textContent = els.nCols().value;
    els.kLabel().textContent = els.kBound().value;
  };

  ['input', 'change'].forEach(ev => {
    els.dDim().addEventListener(ev, updateLabels);
    els.nCols().addEventListener(ev, updateLabels);
    els.kBound().addEventListener(ev, updateLabels);
  });

  els.btnGenerate().addEventListener('click', () => {
    setStatus('Calculating…');
    setTimeout(generate, 20);
  });

  // Auto-run once with defaults
  setTimeout(() => els.btnGenerate().click(), 60);
});
