"""Figures 1-3 for the BD-PRS x plasma protein manuscript (design: basic model on all proteins = primary analysis).

Input directory (one per frozen analysis version):
  PROVENANCE.txt  first line "FROZEN <analysis_id> <code_commit>", "DRAFT ...", "SYNTHETIC ..." or "MOCK <description>"
                  (written by run_all.R; MOCK = hand-made layout inputs from thesis values)
  primary.csv     protein, estimate, ci_low, ci_high, p, q, n        basic model, every protein analysed (one row each)
  models.csv      model, protein, beta, ci_low, ci_high, p, n         all models (MODELS) for the proteins that met the
                                                                      primary Bonferroni threshold
  flow.json       {"participants": [[label, value], ...], "proteins": [[label, n, note], ...]}

  python3 -I make_figures.py <input_dir> <out_dir>

MOCK inputs always produce watermarked figures and skip the strict checks that thesis-derived values cannot pass;
DRAFT and SYNTHETIC inputs are validated and watermarked; FROZEN inputs are validated fail-fast (schema, finite values, P in (0, 1], CI around the estimate, unique keys,
complete model grid, flow counts matching the results). Static print figures, legible in greyscale: identity is
carried by marker shape and fill, never by colour alone.
"""
import csv, json, math, os, sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

MM = 1 / 25.4
INK, INK2, MUTED, GRID, ACCENT = '#0b0b0b', '#52514e', '#b9b8b2', '#e4e3df', '#2a78d6'
MODELS = ['Basic', 'Extended 1', 'Extended 2', 'Extended 3', 'Sensitivity']
plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Liberation Sans', 'Arial', 'DejaVu Sans'],
    'font.size': 7, 'axes.labelsize': 7, 'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5,
    'axes.edgecolor': INK2, 'axes.linewidth': 0.6, 'xtick.color': INK2, 'ytick.color': INK2,
    'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'axes.labelcolor': INK,
    'pdf.fonttype': 42, 'svg.fonttype': 'none',
})


class InputError(Exception):
    pass


def num(x):
    return float(x) if x not in ('', None) else math.nan


def read_csv(path, cols):
    with open(path, encoding='utf-8', newline='') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    missing = set(cols) - set(reader.fieldnames or [])
    if missing:
        raise InputError(f'{os.path.basename(path)}: missing columns {sorted(missing)}')
    return rows


def is_selected(r, alpha):
    """Primary selection as decided in R (column 'bonferroni'); older inputs fall back to P < alpha."""
    if r.get('bonferroni', '') != '':
        return r['bonferroni'].upper() == 'TRUE'
    return num(r['p']) < alpha


def is_fdr(r, q_thr):
    if r.get('fdr', '') != '':
        return r['fdr'].upper() == 'TRUE'
    return num(r['q']) < q_thr


def nlp(r):
    """-log10 P: the value computed on the log scale in R if present (no underflow), else from P."""
    if r.get('neg_log10_p', '') != '':
        return num(r['neg_log10_p'])
    return -math.log10(num(r['p']))


def validate(primary, models, flow, alpha):
    """Fail-fast checks for every input except MOCK (FROZEN, DRAFT and SYNTHETIC)."""
    def check_stat(where, est, lo, hi, p, n):
        if not all(math.isfinite(v) for v in (est, lo, hi, p)):
            raise InputError(f'{where}: non-finite estimate, CI or P')
        if not 0 <= p <= 1:
            raise InputError(f'{where}: P = {p} outside [0, 1]')
        if not lo <= est <= hi:
            raise InputError(f'{where}: estimate outside its CI')
        if not n.isdigit():
            raise InputError(f'{where}: n missing')

    names = [r['protein'] for r in primary]
    if len(names) != len(set(names)):
        raise InputError('primary.csv: duplicate proteins')
    for r in primary:
        check_stat(f"primary {r['protein']}", num(r['estimate']), num(r['ci_low']), num(r['ci_high']), num(r['p']), r['n'])
    for r in primary:
        if not 0 <= num(r['q']) <= 1:
            raise InputError(f"primary {r['protein']}: q outside [0, 1]")
        if num(r['p']) == 0 and r.get('neg_log10_p', '') == '':
            raise InputError(f"primary {r['protein']}: P = 0 and no neg_log10_p column")
    selected = {r['protein'] for r in primary if is_selected(r, alpha)}
    keys = [(r['model'], r['protein']) for r in models]
    if len(keys) != len(set(keys)):
        raise InputError('models.csv: duplicate (model, protein) rows')
    expected = {(m, p) for m in MODELS for p in selected}
    if set(keys) != expected:
        raise InputError(f'models.csv must hold every model for every selected protein; '
                         f'missing {sorted(expected - set(keys))[:5]}, extra {sorted(set(keys) - expected)[:5]}')
    for r in models:
        check_stat(f"{r['model']} {r['protein']}", num(r['beta']), num(r['ci_low']), num(r['ci_high']), num(r['p']), r['n'])
    counts = {'analysed': len(primary), 'bonferroni': len(selected),
              'fdr': sum(is_fdr(r, flow.get('fdr_q', 0.05)) for r in primary)}
    for label, n, _ in flow['proteins']:
        key = next((k for k in counts if f'[{k}]' in label), None)
        if key and str(counts[key]) != n.replace(',', ''):
            raise InputError(f'flow.json "{label}" = {n}, but results give {counts[key]}')


WATERMARKS = {'MOCK': 'LAYOUT MOCK - thesis appendix values - not for publication',
              'DRAFT': 'DRAFT - analysis plan not frozen - not for publication',
              'SYNTHETIC': 'SYNTHETIC DATA - pipeline test - not real results'}


def watermark(fig, kind):
    if kind in WATERMARKS:
        fig.text(0.5, 0.5, WATERMARKS[kind], ha='center', va='center',
                 rotation=20, fontsize=11, color='#d0342c', alpha=0.25, zorder=100)


def save(fig, out, name):
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(out, f'{name}.{ext}'), dpi=300, bbox_inches='tight')
    plt.close(fig)


# ---------- Figure 1: flow ----------
def figure1(flow, out, kind):
    height = max(110, 11 * len(flow['participants']))    # keep each box tall enough for two lines of text
    fig, ax = plt.subplots(figsize=(174 * MM, height * MM))
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    def box(x, y, w, h, text, sub, strong):
        ax.add_patch(FancyBboxPatch((x, y - h / 2), w, h, boxstyle='round,pad=0.004,rounding_size=0.01',
                                    fc='white', ec=INK if strong else MUTED, lw=0.8 if strong else 0.6))
        ax.text(x + 0.012, y + (0.2 * h if sub else 0), text, va='center', fontsize=6.5, color=INK)
        if sub:                                          # offsets scale with the box, so text stays inside
            ax.text(x + 0.012, y - 0.22 * h, sub, va='center', fontsize=6.5, color=INK2)

    # participants: main boxes on a spine at the left; each exclusion box hangs off the spine with a connector
    ax.text(0.02, 0.97, 'Participants', fontsize=8, fontweight='bold', color=INK)
    rows = flow['participants']
    step = 0.88 / len(rows)
    ys = [0.90 - i * step for i in range(len(rows))]
    spine_x = 0.035
    for (label, value), y in zip(rows, ys):
        excl = label.startswith('Excluded')
        if excl:
            box(0.08, y, 0.38, step * 0.8, label, value, False)
            ax.plot([spine_x, 0.08], [y, y], color=INK2, lw=0.6)
        else:
            box(0.02, y, 0.40, step * 0.8, label, value, True)
    mains = [y for (label, _), y in zip(rows, ys) if not label.startswith('Excluded')]
    for a, b in zip(mains, mains[1:]):
        ax.annotate('', xy=(spine_x, b + step * 0.4), xytext=(spine_x, a - step * 0.4),
                    arrowprops=dict(arrowstyle='-|>', lw=0.6, color=INK2, mutation_scale=6))

    # proteins: one box per step, count in the right-hand column
    ax.text(0.50, 0.97, 'Proteins', fontsize=8, fontweight='bold', color=INK)
    ax.text(0.97, 0.97, 'n', fontsize=6.5, color=INK2, ha='right')
    prow = flow['proteins']
    ys = [0.90 - i * 0.13 for i in range(len(prow))]
    for (label, n, note), y in zip(prow, ys):
        excl = label.startswith('Excluded')
        box(0.50, y, 0.48, 0.09, label.split(' [')[0], note or None, not excl)
        ax.text(0.97, y, n, va='center', ha='right', fontsize=7, color=INK, fontweight='bold')
    watermark(fig, kind)
    save(fig, out, 'Figure1_flow')


# ---------- Figure 2: primary analysis, all proteins ----------
def figure2(primary, alpha, q_thr, out, kind, xlabel):
    fig, ax = plt.subplots(figsize=(120 * MM, 85 * MM))
    sel = [r for r in primary if is_selected(r, alpha)]
    fdr = [r for r in primary if not is_selected(r, alpha) and is_fdr(r, q_thr)]
    rest = [r for r in primary if r not in sel and r not in fdr]
    groups = [('Other proteins', rest, dict(s=6, c=MUTED, lw=0, marker='o')),
              (f'FDR < {100 * q_thr:g}% only', fdr, dict(s=14, facecolors='white', edgecolors=INK2, lw=0.7, marker='o')),
              (f'P < 0.05/{len(primary)}', sel, dict(s=18, c=ACCENT, edgecolors=INK, lw=0.5, marker='D'))]
    for label, rows, kw in groups:
        ax.scatter([num(r['estimate']) for r in rows], [nlp(r) for r in rows],
                   label=f'{label} (n = {len(rows)})', zorder=3, **kw)
    ax.axvline(0, color=GRID, lw=0.6, zorder=0)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r'$-\log_{10}\,P$')
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)

    # direct labels: inverse associations in an evenly spaced column left of the data, positive ones above their point
    pts = [(r['protein'], num(r['estimate']), nlp(r)) for r in sel]
    x0, x1 = ax.get_xlim()
    span = x1 - x0
    ythr = -math.log10(alpha)
    ax.set_xlim(x0 - 0.30 * span, x1)
    ax.hlines(ythr, x0, x1, color=INK2, lw=0.6, ls=(0, (3, 2)), zorder=1)
    ax.text(x1, ythr + 0.12, f'P = 0.05/{len(primary)}', ha='right', va='bottom', fontsize=6, color=INK2)
    neg = sorted([p for p in pts if p[1] < 0], key=lambda p: -p[2])
    if neg:
        top, bottom = max(p[2] for p in pts) + 1.0, ythr + 0.4
        ax.set_ylim(top=top + 0.4)
        room = int((top - bottom) / 0.45) + 1                # label only as many as fit legibly (smallest P first)
        if len(neg) > room:
            ax.text(x0 - 0.04 * span, bottom - 0.45, f'+{len(neg) - room + 1} more (Table 2)', ha='right',
                    va='center', fontsize=6, color=INK2)
            neg = neg[:room - 1]
        step = (top - bottom) / max(len(neg) - 1, 1)
        for i, (name, x, y) in enumerate(neg):
            ax.annotate(name, xy=(x, y), xytext=(x0 - 0.04 * span, top - i * step), ha='right', va='center',
                        fontsize=6, color=INK, arrowprops=dict(arrowstyle='-', lw=0.4, color=INK2, shrinkA=1, shrinkB=3))
    for name, x, y in (p for p in pts if p[1] >= 0):
        ax.annotate(name, xy=(x, y), xytext=(x, y + 0.9), ha='center', va='bottom', fontsize=6, color=INK,
                    arrowprops=dict(arrowstyle='-', lw=0.4, color=INK2, shrinkA=1, shrinkB=3))
    ticks = [t for t in ax.get_xticks() if x0 - 1e-9 <= t <= x1 + 1e-9]
    ax.set_xticks(ticks)
    ax.spines['bottom'].set_bounds(ticks[0], ticks[-1])
    ax.legend(loc='upper center', frameon=False, fontsize=6, handletextpad=0.3, borderaxespad=0.2,
              bbox_to_anchor=(0.5, 1.12), ncol=3, columnspacing=1.0)
    watermark(fig, kind)
    save(fig, out, 'Figure2_primary')


# ---------- Figure 3: forest plot, one panel per model ----------
def figure3(primary, models, alpha, out, kind, beta_unit):
    selected = {r['protein']: num(r['estimate']) for r in primary if is_selected(r, alpha)}
    if not selected:                                    # nothing passed the primary threshold: say so, do not crash
        fig, ax = plt.subplots(figsize=(120 * MM, 30 * MM))
        ax.set_axis_off()
        ax.text(0.5, 0.5, 'No assay met the primary (Bonferroni) threshold;\nno secondary models were fitted.',
                ha='center', va='center', fontsize=8, color=INK)
        watermark(fig, kind)
        save(fig, out, 'Figure3_forest')
        return
    order = sorted(selected, key=lambda p: selected[p], reverse=True)   # y index 0 = bottom: most negative on top
    by = {(r['model'], r['protein']): r for r in models}
    vals = [num(by[k][c]) for k in by if k[1] in selected for c in ('ci_low', 'ci_high', 'beta')]
    vals = [v for v in vals if math.isfinite(v)]
    lo, hi = min(vals), max(vals)
    pad = 0.1 * (hi - lo)
    fig, axes = plt.subplots(1, len(MODELS), figsize=(174 * MM, (40 + 4.2 * len(order)) * MM), sharey=True)
    for ax, mname in zip(axes, MODELS):
        ns = []
        for i, p in enumerate(order):
            r = by.get((mname, p))
            if r is None:                                   # MOCK only: legacy pipeline did not fit this cell
                ax.text(0, i, 'not fitted', va='center', ha='center', fontsize=5.5, color=MUTED)
                continue
            b, l, h = (num(r[k]) for k in ('beta', 'ci_low', 'ci_high'))
            if r['n'].isdigit():
                ns.append(int(r['n']))
            if math.isfinite(l) and math.isfinite(h):
                ax.plot([l, h], [i, i], color=INK2, lw=0.8, solid_capstyle='round', zorder=2)
            # Secondary models describe change in proteins already selected by the primary analysis, so no
            # significance encoding here (re-applying the threshold would read as replication).
            ax.scatter([b], [i], s=18, marker='o', zorder=3, lw=0.7, facecolors=INK, edgecolors=INK)
        ax.axvline(0, color=INK2, lw=0.6, ls=(0, (3, 2)), zorder=1)
        ax.set_xlim(lo - pad, hi + pad)
        ax.set_title(mname, fontsize=7, color=INK, pad=12)
        n_txt = ('[MISSING]' if not ns else f'{min(ns):,}' if min(ns) == max(ns) else f'{min(ns):,}–{max(ns):,}')
        ax.text(0.5, 1.01, f'n = {n_txt}', transform=ax.transAxes, ha='center', va='bottom', fontsize=6, color=INK2)
        ax.xaxis.grid(True, color=GRID, lw=0.5)
        ax.set_axisbelow(True)
        ax.tick_params(axis='y', length=0)
        for s in ('top', 'right', 'left'):
            ax.spines[s].set_visible(False)
    axes[0].set_yticks(range(len(order)), order)
    axes[0].set_ylim(-0.7, len(order) - 0.3)
    fig.supxlabel(f'β (NPX per {beta_unit} higher BD-PRS) with 95% CI', fontsize=7, y=-0.01)
    fig.subplots_adjust(wspace=0.12)
    watermark(fig, kind)
    save(fig, out, 'Figure3_forest')


def main(src, out):
    with open(os.path.join(src, 'PROVENANCE.txt'), encoding='utf-8') as f:
        provenance = f.readline().strip()
    kind = provenance.split(' ', 1)[0]
    if kind not in ('FROZEN', 'DRAFT', 'SYNTHETIC', 'MOCK'):
        raise InputError('PROVENANCE.txt must start with FROZEN, DRAFT, SYNTHETIC or MOCK')
    primary = read_csv(os.path.join(src, 'primary.csv'), ['protein', 'estimate', 'ci_low', 'ci_high', 'p', 'q', 'n'])
    models = read_csv(os.path.join(src, 'models.csv'), ['model', 'protein', 'beta', 'ci_low', 'ci_high', 'p', 'n'])
    # sample: 'common' (Basic, Extended 1-3) or 'available' (Sensitivity); described in the figure legend
    with open(os.path.join(src, 'flow.json'), encoding='utf-8') as f:
        flow = json.load(f)
    alpha = flow.get('bonferroni_alpha', 0.05) / len(primary)    # the threshold R used (flow.json)
    q_thr = flow.get('fdr_q', 0.05)
    if kind != 'MOCK':
        validate(primary, models, flow, alpha)
    os.makedirs(out, exist_ok=True)
    xlabel = flow.get('primary_xlabel', 'β (NPX per SD higher BD-PRS), basic model')
    figure1(flow, out, kind)
    figure2(primary, alpha, q_thr, out, kind, xlabel)
    figure3(primary, models, alpha, out, kind, flow.get('beta_unit', '1-SD'))
    print(f'{provenance} -> {out}')


if __name__ == '__main__':
    try:
        main(sys.argv[1], sys.argv[2])
    except InputError as e:
        sys.exit(f'input check failed: {e}')
