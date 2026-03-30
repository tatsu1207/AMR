/* AMR Predictor - Frontend Logic */

const $ = (sel) => document.querySelector(sel);
const $$ = (sel) => document.querySelectorAll(sel);

let selectedFiles = [];
let currentJobs = []; // [{id, filename, status, results}]
let pollTimer = null;
let currentResults = null;
let currentSort = { key: null, asc: true };

// ── File Selection ──────────────────────────────────────────────────
const dropZone = $('#drop-zone');
const fileInput = $('#file-input');
const analyzeBtn = $('#analyze-btn');

dropZone.addEventListener('click', () => fileInput.click());
dropZone.addEventListener('dragover', (e) => { e.preventDefault(); dropZone.classList.add('drag-over'); });
dropZone.addEventListener('dragleave', () => dropZone.classList.remove('drag-over'));
dropZone.addEventListener('drop', (e) => {
    e.preventDefault();
    dropZone.classList.remove('drag-over');
    if (e.dataTransfer.files.length > 0) selectFiles(e.dataTransfer.files);
});

// Prevent browser from opening files dropped outside the drop zone
document.body.addEventListener('dragover', (e) => { if (!dropZone.contains(e.target)) e.preventDefault(); });
document.body.addEventListener('drop', (e) => { if (!dropZone.contains(e.target)) e.preventDefault(); });
fileInput.addEventListener('change', () => { if (fileInput.files.length > 0) selectFiles(fileInput.files); });

$('#clear-file').addEventListener('click', clearFile);
$('#analyze-btn').addEventListener('click', startAnalysis);

function selectFiles(fileList) {
    const validExts = ['.fasta', '.fa', '.fna', '.fasta.gz'];
    const valid = [];
    for (const file of fileList) {
        const name = file.name.toLowerCase();
        if (!validExts.some(ext => name.endsWith(ext))) {
            alert(`Skipped "${file.name}": not a FASTA file.`);
            continue;
        }
        if (file.size > 100 * 1024 * 1024) {
            alert(`Skipped "${file.name}": exceeds 100 MB.`);
            continue;
        }
        valid.push(file);
    }
    if (valid.length === 0) return;
    selectedFiles = valid;
    const totalSize = valid.reduce((s, f) => s + f.size, 0);
    $('#selected-filecount').textContent = valid.length === 1
        ? valid[0].name
        : `${valid.length} files selected`;
    $('#selected-filesize').textContent = formatSize(totalSize);
    const filelist = $('#selected-filelist');
    filelist.innerHTML = valid.map(f =>
        `<li>${f.name} <span class="filesize">(${formatSize(f.size)})</span></li>`
    ).join('');
    filelist.classList.toggle('hidden', valid.length <= 1);
    $('#file-selected').classList.remove('hidden');
    dropZone.classList.add('hidden');
    analyzeBtn.disabled = false;
}

function clearFile() {
    selectedFiles = [];
    fileInput.value = '';
    $('#file-selected').classList.add('hidden');
    dropZone.classList.remove('hidden');
    analyzeBtn.disabled = true;
}

function formatSize(bytes) {
    if (bytes < 1024) return bytes + ' B';
    if (bytes < 1024 * 1024) return (bytes / 1024).toFixed(1) + ' KB';
    return (bytes / (1024 * 1024)).toFixed(1) + ' MB';
}

// ── Upload & Analysis ───────────────────────────────────────────────
async function startAnalysis() {
    if (selectedFiles.length === 0) return;

    analyzeBtn.disabled = true;
    analyzeBtn.textContent = 'Uploading...';
    currentJobs = [];

    try {
        for (const file of selectedFiles) {
            const formData = new FormData();
            formData.append('file', file);
            const resp = await fetch('/api/upload', { method: 'POST', body: formData });
            const data = await resp.json();
            if (!resp.ok) {
                showError(data.detail || `Upload failed for ${file.name}`);
                return;
            }
            currentJobs.push({ id: data.job_id, filename: file.name, status: null, results: null });
        }

        showSection('progress');
        renderJobList();
        startPolling();
    } catch (err) {
        showError('Network error: ' + err.message);
    }
}

function renderJobList() {
    // Show per-job status when multiple files
    if (currentJobs.length <= 1) return;
    let container = $('#job-list');
    if (!container) {
        container = document.createElement('div');
        container.id = 'job-list';
        container.className = 'job-list';
        $('#progress-section').appendChild(container);
    }
    container.innerHTML = currentJobs.map(j => {
        const st = j.status;
        const label = st ? st.stage_label : 'Queued';
        const pct = st ? st.progress : 0;
        const cls = st && st.stage === 'complete' ? 'job-done' : st && st.stage === 'failed' ? 'job-failed' : '';
        return `<div class="job-item ${cls}">
            <span class="job-filename">${j.filename}</span>
            <span class="job-status">${label} (${pct}%)</span>
        </div>`;
    }).join('');
}

// ── Progress Polling ────────────────────────────────────────────────
function startPolling() {
    pollTimer = setInterval(async () => {
        try {
            let allDone = true;
            let anyFailed = false;
            let failMsg = '';

            for (const job of currentJobs) {
                if (job.status && (job.status.stage === 'complete' || job.status.stage === 'failed')) continue;
                const resp = await fetch(`/api/job/${job.id}`);
                job.status = await resp.json();
                if (job.status.stage === 'failed') { anyFailed = true; failMsg = job.status.error; }
                if (job.status.stage !== 'complete' && job.status.stage !== 'failed') allDone = false;
            }

            // Update progress bar with aggregate
            if (currentJobs.length === 1) {
                updateProgress(currentJobs[0].status);
            } else {
                const avgPct = Math.round(currentJobs.reduce((s, j) => s + (j.status ? j.status.progress : 0), 0) / currentJobs.length);
                const done = currentJobs.filter(j => j.status && j.status.stage === 'complete').length;
                $('#progress-fill').style.width = avgPct + '%';
                $('#progress-label').textContent = `${done} / ${currentJobs.length} samples complete`;
                $('#progress-detail').textContent = '';
                renderJobList();
            }

            if (allDone) {
                clearInterval(pollTimer);
                if (anyFailed && currentJobs.length === 1) {
                    showError(failMsg || 'Analysis failed');
                } else {
                    await loadResults();
                }
            }
        } catch (err) {
            console.error('Poll error:', err);
        }
    }, 1000);
}

function updateProgress(status) {
    $('#progress-fill').style.width = status.progress + '%';
    $('#progress-label').textContent = status.stage_label;

    let detail = '';
    if (status.species) detail = `Species: ${status.species}`;
    if (status.mlst_st) detail += ` | ST${status.mlst_st}`;
    $('#progress-detail').textContent = detail;

    // Update step indicators
    const stages = ['identifying_species', 'detecting_genes', 'analyzing_expression', 'detecting_mutations', 'running_predictions'];
    const currentIdx = stages.indexOf(status.stage);
    $$('.step').forEach((el, i) => {
        el.classList.remove('active', 'done');
        if (i < currentIdx) el.classList.add('done');
        else if (i === currentIdx) el.classList.add('active');
    });
}

// ── Results ─────────────────────────────────────────────────────────
// Per-strain state for filters/sorting
let strainStates = {}; // {idx: {results, sort: {key, asc}}}

async function loadResults() {
    try {
        for (const job of currentJobs) {
            if (job.status && job.status.stage === 'complete') {
                const resp = await fetch(`/api/results/${job.id}`);
                job.results = await resp.json();
            }
        }
        renderAllResults();
        showSection('results');
    } catch (err) {
        showError('Failed to load results: ' + err.message);
    }
}

function renderAllResults() {
    const container = $('#results-content');
    container.innerHTML = '';
    strainStates = {};

    const completedJobs = currentJobs.filter(j => j.results);
    const failedJobs = currentJobs.filter(j => j.status && j.status.stage === 'failed');

    // Build tabs if multiple samples, otherwise just show single result
    const showTabs = completedJobs.length > 1;

    // Header + actions
    container.insertAdjacentHTML('beforeend', `
        <div class="results-header">
            <h2>Prediction Results (${completedJobs.length} sample${completedJobs.length !== 1 ? 's' : ''})</h2>
            <div class="results-actions">
                ${showTabs ? '<button id="download-all-csv" class="btn btn-outline btn-sm">Download All CSV</button>' : ''}
                <button id="new-analysis" class="btn btn-outline btn-sm">New Analysis</button>
            </div>
        </div>
    `);

    // Failed jobs warning
    if (failedJobs.length > 0) {
        container.insertAdjacentHTML('beforeend', failedJobs.map(j =>
            `<div class="strain-failed"><strong>${j.filename}</strong>: ${j.status.error || 'Analysis failed'}</div>`
        ).join(''));
    }

    // Tab bar (only for multi-sample)
    if (showTabs) {
        const tabsHtml = `<div class="tab-bar">${completedJobs.map((j, i) => {
            const r = j.results;
            return `<button class="tab ${i === 0 ? 'tab-active' : ''}" data-tab="${i}">
                <span class="tab-name">${j.filename}</span>
                <span class="tab-info">${r.species_display} | ST${r.mlst_st}</span>
                <span class="tab-badges">
                    <span class="badge badge-resistant badge-sm">${r.n_resistant}R</span>
                    <span class="badge badge-susceptible badge-sm">${r.n_susceptible}S</span>
                </span>
            </button>`;
        }).join('')}</div>`;
        container.insertAdjacentHTML('beforeend', tabsHtml);
    }

    // Build all strain panels
    completedJobs.forEach((job, idx) => {
        strainStates[idx] = { results: job.results, sort: { key: null, asc: true } };
        const panelHtml = buildStrainPanel(job, idx, !showTabs || idx === 0);
        container.insertAdjacentHTML('beforeend', panelHtml);
        wireStrainEvents(idx);
    });

    // Wire tab clicks
    if (showTabs) {
        document.querySelectorAll('.tab').forEach(tab => {
            tab.onclick = () => switchTab(parseInt(tab.dataset.tab), completedJobs.length);
        });
        if ($('#download-all-csv')) {
            $('#download-all-csv').onclick = () => completedJobs.forEach(j => downloadCSV(j.results, j.filename));
        }
    }

    $('#new-analysis').onclick = resetUI;
}

function switchTab(activeIdx, total) {
    for (let i = 0; i < total; i++) {
        const panel = document.getElementById(`strain-panel-${i}`);
        const tab = document.querySelector(`.tab[data-tab="${i}"]`);
        if (i === activeIdx) {
            panel.classList.remove('hidden');
            tab.classList.add('tab-active');
        } else {
            panel.classList.add('hidden');
            tab.classList.remove('tab-active');
        }
    }
}

function buildStrainPanel(job, idx, visible) {
    const r = job.results;
    const genes = r.detected_genes || [];
    const classes = [...new Set(r.predictions.map(p => p.drug_class))].sort();
    const classOptions = classes.map(c => `<option value="${c}">${c}</option>`).join('');

    return `
        <div class="strain-panel ${visible ? '' : 'hidden'}" id="strain-panel-${idx}">
            <div class="results-summary">
                <div class="summary-card">
                    <span class="summary-label">Species</span>
                    <span class="summary-value species-name">${r.species_display}</span>
                </div>
                <div class="summary-card">
                    <span class="summary-label">MLST</span>
                    <span class="summary-value">ST${r.mlst_st}</span>
                </div>
                <div class="summary-card summary-resistant">
                    <span class="summary-label">Resistant</span>
                    <span class="summary-value">${r.n_resistant}</span>
                </div>
                <div class="summary-card summary-susceptible">
                    <span class="summary-label">Susceptible</span>
                    <span class="summary-value">${r.n_susceptible}</span>
                </div>
            </div>

            <div class="table-controls">
                <div class="filter-group">
                    <select data-filter="prediction" data-idx="${idx}">
                        <option value="all">All Predictions</option>
                        <option value="Resistant">Resistant Only</option>
                        <option value="Susceptible">Susceptible Only</option>
                    </select>
                    <select data-filter="drugclass" data-idx="${idx}">
                        <option value="all">All Drug Classes</option>
                        ${classOptions}
                    </select>
                    <select data-filter="confidence" data-idx="${idx}">
                        <option value="all">All Confidence</option>
                        <option value="High">High</option>
                        <option value="Moderate">Moderate</option>
                        <option value="Low">Low</option>
                    </select>
                    <input type="text" data-filter="search" data-idx="${idx}" placeholder="Search antibiotic..." />
                </div>
            </div>

            <div class="results-table-container">
                <table class="results-table">
                    <thead>
                        <tr>
                            <th class="sortable" data-sort="antibiotic" data-idx="${idx}">Antibiotic <span class="sort-icon"></span></th>
                            <th class="sortable" data-sort="drug_class" data-idx="${idx}">Drug Class <span class="sort-icon"></span></th>
                            <th class="sortable" data-sort="prediction" data-idx="${idx}">Prediction <span class="sort-icon"></span></th>
                            <th class="sortable" data-sort="probability" data-idx="${idx}">P(Resistant) <span class="sort-icon"></span></th>
                            <th class="sortable" data-sort="confidence" data-idx="${idx}">Confidence <span class="sort-icon"></span></th>
                            <th>Key Determinants</th>
                        </tr>
                    </thead>
                    <tbody id="results-tbody-${idx}"></tbody>
                </table>
            </div>

            <details class="genes-detail" ${genes.length > 0 ? 'open' : ''}>
                <summary>Detected Resistance Genes (${genes.length})</summary>
                <table class="genes-table">
                    <thead>
                        <tr>
                            <th>Gene</th>
                            <th>Drug Class</th>
                            <th>Identity %</th>
                            <th>Coverage %</th>
                            <th>Plasmid</th>
                        </tr>
                    </thead>
                    <tbody id="genes-tbody-${idx}"></tbody>
                </table>
            </details>

            <div style="text-align:right;margin-top:8px">
                <button class="btn btn-outline btn-sm" onclick="downloadCSV(strainStates[${idx}].results, '${job.filename.replace(/'/g, "\\'")}')">Download CSV</button>
            </div>
        </div>
    `;
}

function wireStrainEvents(idx) {
    applyStrainFilters(idx);
    renderStrainGenes(idx);

    document.querySelectorAll(`[data-filter][data-idx="${idx}"]`).forEach(el => {
        const evt = el.tagName === 'INPUT' ? 'input' : 'change';
        el.addEventListener(evt, () => applyStrainFilters(idx));
    });

    document.querySelectorAll(`.sortable[data-idx="${idx}"]`).forEach(th => {
        th.onclick = () => {
            const key = th.dataset.sort;
            const state = strainStates[idx];
            if (state.sort.key === key) {
                state.sort.asc = !state.sort.asc;
            } else {
                state.sort.key = key;
                state.sort.asc = true;
            }
            updateStrainSortIcons(idx);
            applyStrainFilters(idx);
        };
    });
}

function applyStrainFilters(idx) {
    const state = strainStates[idx];
    const getVal = (name) => {
        const el = document.querySelector(`[data-filter="${name}"][data-idx="${idx}"]`);
        return el ? el.value : 'all';
    };

    const predFilter = getVal('prediction');
    const classFilter = getVal('drugclass');
    const confFilter = getVal('confidence');
    const search = (getVal('search') || '').toLowerCase();

    let filtered = state.results.predictions.filter(p => {
        if (predFilter !== 'all' && p.prediction !== predFilter) return false;
        if (classFilter !== 'all' && p.drug_class !== classFilter) return false;
        if (confFilter !== 'all' && p.confidence !== confFilter) return false;
        if (search && !p.antibiotic.replace(/_/g, ' ').toLowerCase().includes(search)) return false;
        return true;
    });

    if (state.sort.key) {
        const key = state.sort.key;
        const dir = state.sort.asc ? 1 : -1;
        filtered.sort((a, b) => {
            let va = a[key], vb = b[key];
            if (key === 'probability') return (va - vb) * dir;
            if (key === 'confidence') {
                const order = { High: 3, Moderate: 2, Low: 1 };
                return ((order[va] || 0) - (order[vb] || 0)) * dir;
            }
            return String(va).localeCompare(String(vb)) * dir;
        });
    } else {
        filtered.sort((a, b) => {
            if (a.prediction !== b.prediction) return a.prediction === 'Resistant' ? -1 : 1;
            return a.drug_class.localeCompare(b.drug_class) || a.antibiotic.localeCompare(b.antibiotic);
        });
    }

    const tbody = document.getElementById(`results-tbody-${idx}`);
    tbody.innerHTML = '';
    for (const pred of filtered) {
        const tr = document.createElement('tr');
        const isR = pred.prediction === 'Resistant';
        const probColor = isR ? 'var(--resistant)' : 'var(--susceptible)';
        const confClass = `confidence-${pred.confidence.toLowerCase()}`;
        const determinants = [...pred.key_genes, ...pred.key_mutations].join(', ') || '-';
        tr.innerHTML = `
            <td><strong>${formatAntibiotic(pred.antibiotic)}</strong></td>
            <td style="font-size:12px;color:var(--gray-500)">${pred.drug_class}</td>
            <td><span class="badge ${isR ? 'badge-resistant' : 'badge-susceptible'}">${pred.prediction}</span></td>
            <td>
                <span class="prob-bar"><span class="prob-bar-fill" style="width:${pred.probability*100}%;background:${probColor}"></span></span>
                ${(pred.probability*100).toFixed(0)}%
            </td>
            <td><span class="${confClass}">${pred.confidence}</span></td>
            <td class="determinants">${determinants}</td>
        `;
        tbody.appendChild(tr);
    }
    if (filtered.length === 0) {
        const tr = document.createElement('tr');
        tr.innerHTML = '<td colspan="6" style="text-align:center;color:var(--gray-500);padding:24px">No results match the current filters</td>';
        tbody.appendChild(tr);
    }
}

function updateStrainSortIcons(idx) {
    const state = strainStates[idx];
    document.querySelectorAll(`.sortable[data-idx="${idx}"]`).forEach(th => {
        const icon = th.querySelector('.sort-icon');
        if (th.dataset.sort === state.sort.key) {
            icon.textContent = state.sort.asc ? ' \u25B2' : ' \u25BC';
        } else {
            icon.textContent = '';
        }
    });
}

function renderStrainGenes(idx) {
    const genes = strainStates[idx].results.detected_genes || [];
    const gtbody = document.getElementById(`genes-tbody-${idx}`);
    gtbody.innerHTML = '';
    for (const g of genes) {
        const tr = document.createElement('tr');
        tr.innerHTML = `
            <td><strong>${g.gene}</strong></td>
            <td>${g.drug_class}</td>
            <td>${g.identity.toFixed(1)}</td>
            <td>${g.coverage.toFixed(1)}</td>
            <td>${g.on_plasmid ? 'Yes' : '-'}</td>
        `;
        gtbody.appendChild(tr);
    }
}

function formatAntibiotic(name) {
    return name.replace(/_/g, ' ').replace(/\b\w/g, c => c.toUpperCase());
}

function downloadCSV(results, filename) {
    let csv = 'Antibiotic,Drug Class,Prediction,Probability,Confidence,Key Determinants\n';
    for (const p of results.predictions) {
        const det = [...p.key_genes, ...p.key_mutations].join('; ');
        csv += `"${formatAntibiotic(p.antibiotic)}","${p.drug_class}","${p.prediction}",${p.probability.toFixed(3)},"${p.confidence}","${det}"\n`;
    }
    const blob = new Blob([csv], { type: 'text/csv' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    const base = filename ? filename.replace(/\.[^.]+$/, '') : `${results.species}_ST${results.mlst_st}`;
    a.download = `amr_prediction_${base}.csv`;
    a.click();
    URL.revokeObjectURL(url);
}

// ── Section Management ──────────────────────────────────────────────
function showSection(name) {
    ['upload', 'progress', 'results', 'error'].forEach(s => {
        $(`#${s}-section`).classList.toggle('hidden', s !== name);
    });
}

function showError(msg) {
    $('#error-message').textContent = msg;
    showSection('error');
    $('#retry-btn').onclick = resetUI;
}

function resetUI() {
    clearFile();
    analyzeBtn.textContent = 'Analyze Genome';
    analyzeBtn.disabled = true;
    currentJobs = [];
    strainStates = {};
    if (pollTimer) clearInterval(pollTimer);
    const jobList = $('#job-list');
    if (jobList) jobList.remove();
    const content = $('#results-content');
    if (content) content.innerHTML = '';
    showSection('upload');
}
