/* SDC3b results website.
 *
 * Reads data/index.js and data/<case>.js written by website/build_data.py.
 * They are loaded as <script> tags, not fetch(), so the page also works from file://.
 * 2-D marginals arrive as base64 uint8 maps of s = -log10(1 - enclosed mass),
 * so a credible region of probability p is the set {s < -log10(1 - p)}.
 */

const RELEASE_URL = 'https://github.com/abonaldi/SDC3b/releases/latest/download';
const LEVELS = [
  { id: 'lvl-99', p: 0.9973, alpha: 0.06, width: 0.8, dash: 'dot' },
  { id: 'lvl-95', p: 0.9545, alpha: 0.12, width: 1.2, dash: 'solid' },
  { id: 'lvl-68', p: 0.6827, alpha: 0.28, width: 1.8, dash: 'solid' },
];
const AXIS_LABELS = ['x<sub>HI,1</sub>', 'x<sub>HI,2</sub>', 'x<sub>HI,3</sub>'];
const BAND_LABELS = ['151–166 MHz', '166–181 MHz', '181–196 MHz'];
const TABS = ['corner', 'intervals', 'metrics', 'data'];

const state = {
  index: null,
  cases: {},        // case -> data
  selected: {},     // case -> Set of team names
  levelCache: new Map(),
  caseName: 'PS1',
  tab: 'corner',
};

const $ = id => document.getElementById(id);

/* ---------- helpers ---------- */

function cssVar(name) {
  return getComputedStyle(document.documentElement).getPropertyValue(name).trim();
}

function rgba(hex, a) {
  const n = parseInt(hex.slice(1), 16);
  return `rgba(${(n >> 16) & 255},${(n >> 8) & 255},${n & 255},${a})`;
}

function fmt(v, digits = 3) {
  if (v === null || v === undefined) return '–';
  if (Math.abs(v) >= 1000) return '≫1';
  return v.toFixed(digits);
}

function levelGrid(caseName, team, key, n, sMax) {
  const id = `${caseName}|${team.name}|${key}`;
  if (state.levelCache.has(id)) return state.levelCache.get(id);
  const bytes = Uint8Array.from(atob(team.marg2d[key]), c => c.charCodeAt(0));
  const z = [];
  for (let r = 0; r < n; r++) {
    const row = new Array(n);
    for (let c = 0; c < n; c++) row[c] = bytes[r * n + c] * sMax / 255;
    z.push(row);
  }
  state.levelCache.set(id, z);
  return z;
}

function plotTheme() {
  return {
    fg: cssVar('--fg'), muted: cssVar('--muted'), grid: cssVar('--grid'),
    axis: cssVar('--axis'), truth: cssVar('--truth'), range: cssVar('--range'),
    surface: cssVar('--surface'), font: cssVar('--font'),
  };
}

function axisStyle(t, narrow = false) {
  // On narrow screens the panels are small; keep tick labels away from panel edges.
  const ticks = narrow ? { tickvals: [0.25, 0.5, 0.75] } : { nticks: 5 };
  return {
    showgrid: true, gridcolor: t.grid, zeroline: false,
    showline: true, linecolor: t.axis, mirror: false,
    ticks: 'outside', tickcolor: t.axis, ticklen: 4,
    tickfont: { size: narrow ? 10 : 11, color: t.muted }, fixedrange: true, ...ticks,
  };
}

const isNarrow = () => window.innerWidth < 600;

const PLOT_CONFIG = {
  responsive: true, displaylogo: false,
  modeBarButtonsToRemove: ['select2d', 'lasso2d', 'zoom2d', 'pan2d', 'zoomIn2d', 'zoomOut2d', 'autoScale2d', 'resetScale2d'],
  toImageButtonOptions: { format: 'png', scale: 2 },
};

/* ---------- data loading ---------- */

function loadData(key) {
  window.SDC3B_DATA = window.SDC3B_DATA || {};
  if (window.SDC3B_DATA[key]) return Promise.resolve(window.SDC3B_DATA[key]);
  return new Promise((resolve, reject) => {
    const el = document.createElement('script');
    el.src = `data/${key}.js`;
    el.onload = () => window.SDC3B_DATA[key] ? resolve(window.SDC3B_DATA[key]) : reject(new Error(`${el.src} is empty`));
    el.onerror = () => reject(new Error(`Could not load ${el.src}`));
    document.head.appendChild(el);
  });
}

async function loadCase(caseName) {
  if (!state.cases[caseName]) {
    const data = await loadData(caseName);
    state.cases[caseName] = data;
    const teams = data.teams;
    const initial = teams.slice(0, teams.length <= 6 ? teams.length : 5).map(t => t.name);
    state.selected[caseName] = new Set(initial);
  }
  return state.cases[caseName];
}

/* ---------- routing ---------- */

function readHash() {
  const [c, t] = location.hash.replace('#', '').split('/');
  const cases = state.index.cases.map(x => x.case);
  if (cases.includes(c)) state.caseName = c;
  if (TABS.includes(t)) state.tab = t;
}

function writeHash() {
  history.replaceState(null, '', `#${state.caseName}/${state.tab}`);
}

/* ---------- static parts ---------- */

function buildChrome() {
  const sw = $('case-switch');
  state.index.cases.forEach(({ case: c }) => {
    const b = document.createElement('button');
    b.textContent = c;
    b.setAttribute('role', 'tab');
    b.dataset.case = c;
    b.addEventListener('click', () => selectCase(c));
    sw.appendChild(b);
  });

  document.querySelectorAll('.tabs button').forEach(b =>
    b.addEventListener('click', () => selectTab(b.dataset.tab)));

  $('generated').textContent = state.index.generated;

  const dl = $('download-list');
  const rows = $('dataset-table');
  state.index.cases.forEach(({ case: c, info, simulation, model, n_teams }) => {
    const li = document.createElement('li');
    const a = document.createElement('a');
    a.href = `${RELEASE_URL}/SDC3b_posteriors_${c}.zip`;
    a.textContent = `SDC3b_posteriors_${c}.zip`;
    const note = document.createElement('span');
    note.className = 'muted';
    note.textContent = ` · ${n_teams} posteriors`;
    li.append(a, note);
    dl.appendChild(li);

    const tr = document.createElement('tr');
    [c, info, `${simulation}, ${model.replace('_', ' ')}`].forEach(v => {
      const td = document.createElement('td');
      td.textContent = v;
      tr.appendChild(td);
    });
    rows.appendChild(tr);
  });

  $('team-search').addEventListener('input', filterTeams);
  $('sel-all').addEventListener('click', () => setSelection(current().teams.map(t => t.name)));
  $('sel-none').addEventListener('click', () => setSelection([]));
  ['lvl-68', 'lvl-95', 'lvl-99', 'opt-fill', 'opt-truth', 'opt-range'].forEach(id =>
    $(id).addEventListener('change', scheduleRender));
  $('csv-btn').addEventListener('click', downloadCSV);

  matchMedia('(prefers-color-scheme: dark)').addEventListener('change', scheduleRender);
  matchMedia('(max-width: 599px)').addEventListener('change', scheduleRender);
}

function current() { return state.cases[state.caseName]; }

async function selectCase(c) {
  state.caseName = c;
  writeHash();
  document.querySelectorAll('#case-switch button').forEach(b =>
    b.setAttribute('aria-selected', String(b.dataset.case === c)));
  const data = await loadCase(c);
  if (state.caseName !== c) return;
  $('case-info').textContent = `${data.info} ${data.teams.length} submissions.`;
  buildTeamList(data);
  scheduleRender();
}

function selectTab(tab) {
  state.tab = tab;
  writeHash();
  document.querySelectorAll('.tabs button').forEach(b =>
    b.setAttribute('aria-selected', String(b.dataset.tab === tab)));
  TABS.forEach(t => { $(`panel-${t}`).hidden = t !== tab; });
  scheduleRender();
}

/* ---------- team list ---------- */

function buildTeamList(data) {
  const list = $('team-list');
  list.innerHTML = '';
  const sel = state.selected[state.caseName];
  data.teams.forEach(team => {
    const li = document.createElement('li');
    li.dataset.name = team.name.toLowerCase();
    const label = document.createElement('label');
    label.style.setProperty('--swatch', team.colour);
    const input = document.createElement('input');
    input.type = 'checkbox';
    input.checked = sel.has(team.name);
    input.addEventListener('change', () => {
      input.checked ? sel.add(team.name) : sel.delete(team.name);
      updateCount();
      scheduleRender();
    });
    const sw = document.createElement('span');
    sw.className = 'swatch';
    const name = document.createElement('span');
    name.className = 'name';
    name.textContent = team.name;
    label.append(input, sw, name);
    li.appendChild(label);
    list.appendChild(li);
  });
  filterTeams();
  updateCount();
}

function filterTeams() {
  const term = $('team-search').value.trim().toLowerCase();
  document.querySelectorAll('#team-list li').forEach(li => {
    li.hidden = term !== '' && !li.dataset.name.includes(term);
  });
}

function setSelection(names) {
  state.selected[state.caseName] = new Set(names);
  buildTeamList(current());
  scheduleRender();
}

function updateCount() {
  $('team-count').textContent = state.selected[state.caseName].size;
}

/* ---------- rendering ---------- */

let renderQueued = false;
function scheduleRender() {
  if (renderQueued) return;
  renderQueued = true;
  requestAnimationFrame(() => {
    renderQueued = false;
    const data = current();
    if (!data) return;
    if (state.tab === 'corner') renderCorner(data);
    if (state.tab === 'intervals') renderIntervals(data);
    if (state.tab === 'metrics') renderMetrics(data);
  });
}

function renderCorner(data) {
  const t = plotTheme();
  const sel = state.selected[state.caseName];
  const teams = data.teams.filter(tm => sel.has(tm.name));
  $('corner-empty').hidden = teams.length > 0;

  const coords = data.coords;
  const n = coords.length;
  const truth = data.truth;
  const fill = $('opt-fill').checked;
  const showTruth = $('opt-truth').checked;
  const showRange = $('opt-range').checked;
  const levels = LEVELS.filter(l => $(l.id).checked);

  // Panel (row, col) -> axis suffix; diagonal panels hold 1-D marginals.
  const ax = { '0,0': '', '1,0': '2', '1,1': '3', '2,0': '4', '2,1': '5', '2,2': '6' };
  const traces = [];
  const shapes = [];

  teams.forEach(team => {
    for (let p = 0; p < 3; p++) {
      const s = ax[`${p},${p}`];
      traces.push({
        x: coords, y: team.marg1d[p], xaxis: `x${s}`, yaxis: `y${s}`,
        type: 'scatter', mode: 'lines', name: team.name,
        line: { color: team.colour, width: 2 },
        fill: fill ? 'tozeroy' : 'none', fillcolor: rgba(team.colour, 0.08),
        hovertemplate: `${team.name}<br>x = %{x:.3f}<extra></extra>`,
      });
    }
    [[0, 1], [0, 2], [1, 2]].forEach(([i, j]) => {
      const s = ax[`${j},${i}`];
      const z = levelGrid(state.caseName, team, `${i}_${j}`, n, data.s_max);
      levels.forEach(l => {
        traces.push({
          x: coords, y: coords, z, xaxis: `x${s}`, yaxis: `y${s}`,
          type: 'contour', showscale: false, hoverinfo: 'skip',
          // Plotly shades the region that fails the constraint, so '>' fills {s < value}.
          contours: { type: 'constraint', operation: '>', value: -Math.log10(1 - l.p) },
          line: { color: team.colour, width: l.width, dash: l.dash, smoothing: 1 },
          fillcolor: fill ? rgba(team.colour, l.alpha) : 'rgba(0,0,0,0)',
        });
      });
    });
  });

  for (let p = 0; p < 3; p++) {
    const s = ax[`${p},${p}`];
    const [lo, hi] = truth.ranges[p];
    if (showRange) shapes.push({ type: 'rect', xref: `x${s}`, yref: `y${s} domain`, x0: lo, x1: hi, y0: 0, y1: 1, fillcolor: t.range, line: { width: 0 }, layer: 'below' });
    if (showTruth) shapes.push({ type: 'line', xref: `x${s}`, yref: `y${s} domain`, x0: truth.centre[p], x1: truth.centre[p], y0: 0, y1: 1, line: { color: t.truth, width: 1.2, dash: 'dash' } });
  }
  [[0, 1], [0, 2], [1, 2]].forEach(([i, j]) => {
    const s = ax[`${j},${i}`];
    if (showRange) shapes.push({ type: 'rect', xref: `x${s}`, yref: `y${s}`, x0: truth.ranges[i][0], x1: truth.ranges[i][1], y0: truth.ranges[j][0], y1: truth.ranges[j][1], fillcolor: t.range, line: { width: 0 }, layer: 'below' });
    if (showTruth) traces.push({
      x: [truth.centre[i]], y: [truth.centre[j]], xaxis: `x${s}`, yaxis: `y${s}`,
      type: 'scatter', mode: 'markers', hovertemplate: 'True value<br>(%{x:.3f}, %{y:.3f})<extra></extra>',
      marker: { symbol: 'star', size: 12, color: t.truth, line: { color: t.surface, width: 1 } },
    });
  });

  const narrow = isNarrow();
  const left = 0.08, bottom = 0.07, gap = narrow ? 0.035 : 0.025;
  const cell = (1 - left - 2 * gap) / 3;
  const cellY = (1 - bottom - 2 * gap) / 3;
  const xd = c => [left + c * (cell + gap), left + c * (cell + gap) + cell];
  const yd = r => [1 - r * (cellY + gap) - cellY, 1 - r * (cellY + gap)];
  const base = axisStyle(t, narrow);
  const layout = {
    paper_bgcolor: t.surface, plot_bgcolor: t.surface,
    font: { family: t.font, color: t.fg, size: 12 },
    showlegend: false, margin: narrow ? { t: 6, r: 6, b: 40, l: 40 } : { t: 10, r: 10, b: 50, l: 60 },
    hoverlabel: { font: { family: t.font } },
    shapes,
  };
  Object.entries(ax).forEach(([rc, s]) => {
    const [r, c] = rc.split(',').map(Number);
    const diag = r === c;
    layout[`xaxis${s}`] = {
      ...base, domain: xd(c), range: [0, 1], anchor: `y${s}`,
      showticklabels: r === 2, title: r === 2 ? { text: AXIS_LABELS[c], font: { size: 14 } } : undefined,
    };
    layout[`yaxis${s}`] = {
      ...base, domain: yd(r), anchor: `x${s}`,
      range: diag ? [0, 1.08] : [0, 1],
      showticklabels: !diag && c === 0, showgrid: !diag,
      title: !diag && c === 0 ? { text: AXIS_LABELS[r], font: { size: 14 } } : undefined,
    };
  });

  Plotly.react('corner-plot', traces, layout, { ...PLOT_CONFIG, toImageButtonOptions: { ...PLOT_CONFIG.toImageButtonOptions, filename: `SDC3b_${state.caseName}_corner` } });
}

function renderIntervals(data) {
  const t = plotTheme();
  const teams = data.teams;
  const nT = teams.length;
  const traces = [];
  const shapes = [];
  const gap = 0.03;
  const w = (1 - 2 * gap) / 3;

  for (let p = 0; p < 3; p++) {
    const s = p === 0 ? '' : String(p + 1);
    teams.forEach((team, k) => {
      const st = team.stats;
      const y = nT - 1 - k;
      const hover = `${team.name}<br>median %{x:.3f}<br>68%: ${fmt(st.ci68[p][0])}–${fmt(st.ci68[p][1])}<br>95%: ${fmt(st.ci95[p][0])}–${fmt(st.ci95[p][1])}<extra></extra>`;
      traces.push(
        { x: st.ci95[p], y: [y, y], xaxis: `x${s}`, yaxis: 'y', type: 'scatter', mode: 'lines',
          line: { color: team.colour, width: 2 }, hoverinfo: 'skip' },
        { x: st.ci68[p], y: [y, y], xaxis: `x${s}`, yaxis: 'y', type: 'scatter', mode: 'lines',
          line: { color: rgba(team.colour, 0.55), width: 11 }, hoverinfo: 'skip' },
        { x: [st.median[p]], y: [y], xaxis: `x${s}`, yaxis: 'y', type: 'scatter', mode: 'markers',
          marker: { color: team.colour, size: 9, symbol: 'diamond', line: { color: t.fg, width: 0.8 } },
          hovertemplate: hover },
      );
    });
    const [lo, hi] = data.truth.ranges[p];
    shapes.push(
      { type: 'rect', xref: `x${s}`, yref: 'paper', x0: lo, x1: hi, y0: 0, y1: 1, fillcolor: t.range, line: { width: 0 }, layer: 'below' },
      { type: 'line', xref: `x${s}`, yref: 'paper', x0: data.truth.centre[p], x1: data.truth.centre[p], y0: 0, y1: 1, line: { color: t.truth, width: 1.2, dash: 'dash' } },
    );
  }

  const narrow = isNarrow();
  const base = { ...axisStyle(t, narrow), ...(narrow ? { tickvals: [0, 0.5, 1], tickangle: 0 } : {}) };
  const layout = {
    paper_bgcolor: t.surface, plot_bgcolor: t.surface,
    font: { family: t.font, color: t.fg, size: 12 },
    showlegend: false, height: Math.max(220, nT * 28 + 110),
    margin: { t: 30, r: 10, b: 50, l: 10 }, shapes,
    hoverlabel: { font: { family: t.font } },
    yaxis: {
      ...base, showgrid: false, showline: false, ticks: '',
      tickvals: teams.map((_, k) => nT - 1 - k), ticktext: teams.map(tm => tm.name),
      tickfont: { size: 12, color: t.fg }, range: [-0.7, nT - 0.3], automargin: true,
    },
    annotations: [],
  };
  for (let p = 0; p < 3; p++) {
    const s = p === 0 ? '' : String(p + 1);
    layout[`xaxis${s}`] = {
      ...base, domain: [p * (w + gap), p * (w + gap) + w], range: [0, 1], anchor: 'y',
      title: { text: AXIS_LABELS[p], font: { size: 14 } },
    };
    if (!narrow) layout.annotations.push({
      xref: `x${s} domain`, yref: 'paper', x: 0.5, y: 1, yanchor: 'bottom', yshift: 6,
      text: BAND_LABELS[p], showarrow: false, font: { size: 12, color: t.muted },
    });
  }

  Plotly.react('interval-plot', traces, layout, { ...PLOT_CONFIG, toImageButtonOptions: { ...PLOT_CONFIG.toImageButtonOptions, filename: `SDC3b_${state.caseName}_intervals` } });
}

// Teams stay in paper order; the table is deliberately not sortable.
const COLUMNS = ['Team', 'Z<sub>1</sub>', 'Z<sub>2</sub>', 'Z<sub>3</sub>', 'D<sub>M</sub>', 'S<sub>SDC3b</sub>'];

function renderMetrics(data) {
  const table = $('metrics-table');
  const thead = document.createElement('thead');
  const hr = document.createElement('tr');
  COLUMNS.forEach(label => {
    const th = document.createElement('th');
    th.innerHTML = label;
    th.scope = 'col';
    hr.appendChild(th);
  });
  thead.appendChild(hr);

  const tbody = document.createElement('tbody');
  data.teams.forEach(tm => {
    const tr = document.createElement('tr');
    const name = document.createElement('td');
    const dot = document.createElement('span');
    dot.className = 'dot';
    dot.style.background = tm.colour;
    name.append(dot, tm.name);
    tr.appendChild(name);
    tm.stats.z.forEach(v => tr.appendChild(cell(fmt(v))));
    tr.appendChild(cell(fmt(tm.stats.mahalanobis)));
    tr.appendChild(cell(fmt(tm.stats.score, 4)));
    tbody.appendChild(tr);
  });

  table.replaceChildren(thead, tbody);
}

function cell(text) {
  const td = document.createElement('td');
  td.textContent = text;
  return td;
}

function downloadCSV() {
  const data = current();
  const head = ['team', 'Z1', 'Z2', 'Z3', 'D_M', 'S_SDC3b',
    'median1', 'median2', 'median3', 'ci68_1_lo', 'ci68_1_hi', 'ci68_2_lo', 'ci68_2_hi', 'ci68_3_lo', 'ci68_3_hi'];
  const lines = [head.join(',')];
  data.teams.forEach(tm => {
    const s = tm.stats;
    const row = [`"${tm.name}"`, ...s.z, s.mahalanobis, s.score, ...s.median, ...s.ci68.flat()];
    lines.push(row.map(v => v ?? '').join(','));
  });
  const blob = new Blob([lines.join('\n') + '\n'], { type: 'text/csv' });
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = `SDC3b_${data.case}_metrics.csv`;
  a.click();
  URL.revokeObjectURL(a.href);
}

/* ---------- start ---------- */

(async function main() {
  try {
    state.index = await loadData('index');
  } catch (err) {
    $('case-info').textContent = 'Could not load the results data (website/data/). Run website/build_data.py to create it.';
    return;
  }
  readHash();
  buildChrome();
  selectTab(state.tab);
  await selectCase(state.caseName);
})();
