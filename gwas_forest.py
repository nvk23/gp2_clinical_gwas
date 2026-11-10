"""Forest plots of per-study effect sizes, with or without a meta-analysis row.

Per-study estimates come from the pipeline's `*_y_allresults.tsv` files; the
fixed- and random-effect summaries come from PLINK's `.meta` output.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import numpy as np
import pandas as pd
from scipy.stats import norm

# Columns read from each per-study results file.
NEEDED_COLS = ['ID', 'A1', 'A1_FREQ', 'OBS_CT', 'BETA', 'SE', 'P']

# chrom -> {pos: gene}. Pre-populate to title plots with gene names;
# otherwise the title falls back to chrom:pos.
GENE_LOOKUP: dict[str, dict[int, str]] = {}

STUDY_COLOR = '#1f5f9e'
SUB_COLOR = '#8fb3d9'
FIXED_COLOR = '#c0392b'
RANDOM_COLOR = '#27ae60'


def load_study_subset(filepath, snps) -> pd.DataFrame:
    """Rows for `snps` (a set of full chrom:pos:ref:alt IDs) from one study."""
    df = pd.read_csv(filepath, sep='\t', usecols=NEEDED_COLS,
                     dtype={'ID': str, 'A1': str})
    return df[df['ID'].isin(set(snps))].copy()


def load_locus(filepath, locus_prefix: str) -> pd.DataFrame:
    """Rows whose ID starts with `locus_prefix` (e.g. 'chr7:97841613:').

    Useful when the alleles at a locus differ between result files.
    """
    df = pd.read_csv(filepath, sep='\t', usecols=NEEDED_COLS,
                     dtype={'ID': str, 'A1': str})
    return df[df['ID'].str.startswith(locus_prefix)].copy()


def _format_p(p) -> str:
    if pd.isna(p):
        return 'NA'
    return f'{p:.2e}' if p < 1e-3 else f'{p:.3f}'


def _format_af(af) -> str:
    return f'{af:.3f}' if pd.notna(af) else 'NA'


def _column_layout(study_labels, plot_w=6.0):
    """Size the table columns from the labels they actually have to hold.

    Everything is laid out in inches first (study names, AF / N, beta CI and P
    columns, with a fixed gap between each), then converted to the
    axes-fraction anchors used with the blended transform. Long study names
    widen the figure instead of running into the next column.
    """
    char_w = 0.085          # approx width of one character at fontsize 10
    pad = 0.15              # outer margin
    gap = 0.40              # gap between columns

    labels = list(study_labels) + ['Study']
    w_study = char_w * max(len(s) for s in labels)
    w_afn = char_w * len('0.000 / 000000')
    w_beta = char_w * len('-0.000 (-0.000, 0.000)')
    w_p = char_w * len('0.00e-00')

    x_study_l = pad
    x_afn_r = x_study_l + w_study + gap + w_afn
    ax_l = x_afn_r + gap
    ax_r = ax_l + plot_w
    x_beta_r = ax_r + gap + w_beta
    x_p_r = x_beta_r + gap + w_p
    fig_w = x_p_r + pad

    def frac(x_in):
        return (x_in - ax_l) / plot_w

    return {
        'fig_w': fig_w,
        'left': ax_l / fig_w,
        'right': ax_r / fig_w,
        'X_STUDY': frac(x_study_l),
        'X_AFN': frac(x_afn_r),
        'X_BETACI': frac(x_beta_r),
        'X_P': frac(x_p_r),
        'X_SEPL': frac(x_study_l - 0.05),
        'X_SEPR': frac(x_p_r + 0.05),
    }


def _aligned_rows(data, snp, effect_allele, indent=False):
    """Per-study rows for one SNP, sign-flipped onto the chosen effect allele."""
    rows = []
    for study, df in (data or {}).items():
        sub = df[df['ID'] == snp]
        if sub.empty:
            continue
        r = sub.iloc[0]
        beta, se, af, a1 = r['BETA'], r['SE'], r['A1_FREQ'], r['A1']
        if pd.isna(beta) or pd.isna(se):
            continue
        if a1 != effect_allele:
            beta = -beta
            af = 1 - af if pd.notna(af) else af
        rows.append({'study': f'   {study}' if indent else study,
                     'beta': float(beta), 'se': float(se),
                     'af': float(af) if pd.notna(af) else np.nan,
                     'n': int(r['OBS_CT']), 'p': float(r['P']), 'sub': indent})
    return rows


def _se_from_p(beta, p):
    """SE implied by a reported beta and p-value.

    Deriving the summary CI from PLINK's p-value (rather than recomputing it
    from the per-study SEs) guarantees the drawn CI and the printed p-value
    cannot disagree.
    """
    if not np.isfinite(beta) or not np.isfinite(p) or p <= 0 or p >= 1:
        return float('nan')
    if abs(beta) < 1e-12:
        return float('nan')
    z = norm.ppf(1 - p / 2)
    return abs(beta) / z if z > 0 else float('nan')


def forest_plot(snp, study_data, meta_df=None, *, beta_label='Beta (log OR)',
                gene_name=None, save_path=None, se_thresh=100,
                effect_allele=None, sub_data=None, sub_parent=None):
    """Forest plot for `snp` across every study in `study_data`.

    `meta_df` is a PLINK `.meta` table indexed by SNP; when it is given the
    fixed- and random-effect summaries are drawn below the studies and the
    effect allele is taken from the meta row. Without it the plot shows the
    per-study estimates only, against `effect_allele` (default: the ALT allele
    in the SNP ID).

    Optional `sub_data` ({label: df}, same columns as `study_data`) is drawn as
    indented rows under `sub_parent`, e.g. a per-cohort breakdown of a pooled
    study. Those rows are display-only: excluded from k and from the x-axis
    scale, with whiskers that run off the axis ending in an arrow.

    Studies with |SE| > `se_thresh` are shown as grey diamonds at the edge and
    excluded from the scale and from the heterogeneity statistics.

    Returns a dict of summary statistics (meta fields are NaN when `meta_df`
    is not given), or None when the SNP could not be plotted.
    """
    chrom, pos, ref, alt = snp.split(':')

    m = None
    if meta_df is not None:
        if snp not in meta_df.index:
            print(f'  [skip] {snp}: no meta result')
            return None
        m = meta_df.loc[snp]
        if isinstance(m, pd.DataFrame):
            m = m.iloc[0]
        effect_allele = m['A1']
    elif effect_allele is None:
        effect_allele = alt

    if effect_allele == alt:
        arrow_str = f'{ref} → {alt}'
    elif effect_allele == ref:
        arrow_str = f'{alt} → {ref}'
    else:
        arrow_str = f'effect={effect_allele}'

    if gene_name is None:
        gene_name = GENE_LOOKUP.get(chrom, {}).get(int(pos))

    rows = _aligned_rows(study_data, snp, effect_allele)
    if not rows:
        print(f'  [skip] {snp}: no per-study results')
        return None
    sdf = pd.DataFrame(rows).sort_values('study').reset_index(drop=True)

    # Insert the breakdown rows directly under their parent study.
    sub_rows = _aligned_rows(sub_data, snp, effect_allele, indent=True) if sub_data else []
    if sub_rows and sub_parent in set(sdf['study']):
        subs = pd.DataFrame(sub_rows).sort_values('study')
        i = sdf.index[sdf['study'] == sub_parent][0] + 1
        sdf = pd.concat([sdf.iloc[:i], subs, sdf.iloc[i:]]).reset_index(drop=True)

    sdf['extreme'] = sdf['se'].abs() > se_thresh
    valid = sdf[~sdf['extreme'] & ~sdf['sub']].copy()
    k = len(valid)
    if k == 0:
        print(f'  [skip] {snp}: no studies under SE threshold')
        return None

    # Meta summaries. `Q` in the .meta file is already the p-value of Cochran's
    # Q test and `I` is the I-squared index, in percent.
    if m is not None:
        Q = float(m['Q']) if pd.notna(m['Q']) else float('nan')
        I2 = float(m['I']) if pd.notna(m['I']) else float('nan')
        beta_fixed, beta_random = float(m['BETA']), float(m['BETA(R)'])
        p_fixed, p_random = float(m['P']), float(m['P(R)'])
        se_fixed = _se_from_p(beta_fixed, p_fixed)
        se_random = _se_from_p(beta_random, p_random)
        ci_fixed = (beta_fixed - 1.96 * se_fixed, beta_fixed + 1.96 * se_fixed)
        ci_random = (beta_random - 1.96 * se_random, beta_random + 1.96 * se_random)
        if not np.isfinite(se_fixed):
            ci_fixed = (beta_fixed, beta_fixed)
        if not np.isfinite(se_random):
            ci_random = (beta_random, beta_random)
    else:
        Q = I2 = beta_fixed = beta_random = float('nan')
        p_fixed = p_random = se_fixed = se_random = float('nan')
        ci_fixed = ci_random = ()

    # x scale: non-extreme studies, plus the meta CIs when present
    lo_vals = [(valid['beta'] - 1.96 * valid['se']).min(), *ci_fixed[:1], *ci_random[:1]]
    hi_vals = [(valid['beta'] + 1.96 * valid['se']).max(), *ci_fixed[1:], *ci_random[1:]]
    x_lo, x_hi = min(lo_vals), max(hi_vals)
    span = x_hi - x_lo if x_hi > x_lo else 1.0
    x_lo -= 0.10 * span
    x_hi += 0.10 * span

    n_studies = len(sdf)
    fig_h = max(6.0, 0.42 * n_studies + 4.0) if m is not None \
        else max(5.0, 0.42 * n_studies + 2.5)

    extra_labels = ['Fixed effect', 'Random effect'] if m is not None else []
    lay = _column_layout(list(sdf['study']) + extra_labels)
    fig, ax = plt.subplots(figsize=(lay['fig_w'], fig_h))
    fig.subplots_adjust(left=lay['left'], right=lay['right'], top=0.88, bottom=0.12)

    blended = mtransforms.blended_transform_factory(ax.transAxes, ax.transData)
    X_STUDY, X_AFN = lay['X_STUDY'], lay['X_AFN']
    X_BETACI, X_P = lay['X_BETACI'], lay['X_P']
    X_SEPL, X_SEPR = lay['X_SEPL'], lay['X_SEPR']

    y_studies = np.arange(n_studies, 0, -1, dtype=float)
    y_meta_info, y_fixed, y_random = 0.2, -0.8, -1.8
    cap = 0.15

    # Header
    y_header = n_studies + 1.2
    for x, text, ha in ((X_STUDY, 'Study', 'left'), (X_AFN, 'AF / N', 'right'),
                        (X_BETACI, 'β (95% CI)', 'right'), (X_P, 'P', 'right')):
        ax.text(x, y_header, text, fontweight='bold', transform=blended,
                va='center', ha=ha)
    ax.plot([X_SEPL, X_SEPR], [y_header - 0.5, y_header - 0.5],
            transform=blended, color='gray', lw=0.6, clip_on=False)

    # Per-study rows, marker area proportional to inverse-variance weight
    max_w = float((1.0 / valid['se'] ** 2).max())
    for i, row in sdf.iterrows():
        y = y_studies[i]
        ax.text(X_STUDY, y, row['study'], transform=blended, va='center', ha='left',
                color='#555' if row['sub'] else 'black')
        ax.text(X_AFN, y, f"{_format_af(row['af'])} / {row['n']}",
                transform=blended, va='center', ha='right', color='#555')

        beta, se = row['beta'], row['se']
        ci_lo, ci_hi = beta - 1.96 * se, beta + 1.96 * se

        if row['extreme']:
            x_edge = x_hi - 0.04 * span if beta > 0 else x_lo + 0.04 * span
            ax.scatter([x_edge], [y], marker='D', color='lightgray', s=80,
                       zorder=5, clip_on=False)
            ax.text(X_BETACI, y, f'{beta:.2f} [extreme SE]',
                    transform=blended, va='center', ha='right', color='#888')
            ax.text(X_P, y, _format_p(row['p']),
                    transform=blended, va='center', ha='right', color='#888')
            continue

        if row['sub']:
            # Breakdown row: clipped to the axis, arrows where it runs off.
            lo, hi = max(ci_lo, x_lo), min(ci_hi, x_hi)
            ax.plot([lo, hi], [y, y], color=SUB_COLOR, lw=1.0, zorder=2)
            for edge, clipped, mk in ((lo, ci_lo < x_lo, '<'), (hi, ci_hi > x_hi, '>')):
                if clipped:
                    ax.scatter([edge], [y], marker=mk, s=25, color=SUB_COLOR, zorder=3)
                else:
                    ax.plot([edge, edge], [y - cap * 0.7, y + cap * 0.7],
                            color=SUB_COLOR, lw=1.0)
            if x_lo <= beta <= x_hi:
                size = min(200, 15 + 120 * (1.0 / se ** 2) / max_w)
                ax.scatter([beta], [y], s=size, marker='s', color=SUB_COLOR, zorder=3)
            color = '#555'
        else:
            ax.plot([ci_lo, ci_hi], [y, y], color=STUDY_COLOR, lw=1.2, zorder=2)
            ax.plot([ci_lo, ci_lo], [y - cap, y + cap], color=STUDY_COLOR, lw=1.2)
            ax.plot([ci_hi, ci_hi], [y - cap, y + cap], color=STUDY_COLOR, lw=1.2)
            ax.scatter([beta], [y], s=25 + 175 * (1.0 / se ** 2) / max_w,
                       marker='s', color=STUDY_COLOR, zorder=3)
            color = 'black'

        ax.text(X_BETACI, y, f'{beta:.3f} ({ci_lo:.3f}, {ci_hi:.3f})',
                transform=blended, va='center', ha='right', color=color)
        ax.text(X_P, y, _format_p(row['p']),
                transform=blended, va='center', ha='right', color=color)

    if m is not None:
        ax.plot([X_SEPL, X_SEPR], [0.6, 0.6], transform=blended,
                color='gray', lw=0.6, clip_on=False)
        ax.text(X_STUDY, y_meta_info,
                f'Meta-analysis (k={k}  I²={I2:.1f}%  Q={_format_p(Q)})',
                transform=blended, va='center', style='italic',
                color='#444', fontsize=9)

        for y, b, ci, p, name, color in (
                (y_fixed, beta_fixed, ci_fixed, p_fixed, 'Fixed effect', FIXED_COLOR),
                (y_random, beta_random, ci_random, p_random, 'Random effect', RANDOM_COLOR)):
            ax.plot([ci[0], ci[1]], [y, y], color=color, lw=1.4, zorder=4)
            ax.plot([ci[0], ci[0]], [y - cap, y + cap], color=color, lw=1.4)
            ax.plot([ci[1], ci[1]], [y - cap, y + cap], color=color, lw=1.4)
            ax.scatter([b], [y], s=140, marker='s', color=color, zorder=5)
            ax.text(X_STUDY, y, name, fontweight='bold', transform=blended,
                    va='center', ha='left', color=color)
            ax.text(X_BETACI, y, f'{b:.3f} ({ci[0]:.3f}, {ci[1]:.3f})',
                    transform=blended, va='center', ha='right', color=color)
            ax.text(X_P, y, _format_p(p), transform=blended,
                    va='center', ha='right', color=color)

    ax.axvline(x=0, color='gray', linestyle='--', lw=0.7, zorder=1)

    title = f'{chrom}:{pos}  {arrow_str}'
    ax.set_title(f'{gene_name}\n{title}' if gene_name else title, fontsize=13, pad=18)

    ax.set_xlim(x_lo, x_hi)
    n_sub = int(sdf['sub'].sum())
    n_extreme = int(sdf['extreme'].sum())
    if m is not None:
        y_bottom = y_random - (1.2 if n_extreme else 0.6) - (0.6 if n_sub else 0)
    else:
        y_bottom = 0.3
    ax.set_ylim(y_bottom, n_studies + 2.0)
    ax.set_xlabel(beta_label)
    ax.set_yticks([])
    for spine in ('left', 'top', 'right'):
        ax.spines[spine].set_visible(False)

    note_y = (y_random - 1.0) if m is not None else 0.5
    axes_data = mtransforms.blended_transform_factory(ax.transAxes, ax.transData)
    if n_extreme:
        excluded = ('x-axis scale and meta-analysis' if m is not None else 'x-axis scale')
        ax.text(0.5, note_y,
                f"* {n_extreme} {'study' if n_extreme == 1 else 'studies'} with "
                f'|SE| > {se_thresh} shown as ◇ but excluded from {excluded}',
                transform=axes_data, va='center', ha='center',
                style='italic', color='#888', fontsize=8)
    if n_sub:
        ax.text(0.5, note_y - (0.6 if n_extreme else 0.0),
                f'Indented rows: per-cohort breakdown of {sub_parent} '
                '(display only; not part of the meta-analysis)',
                transform=axes_data, va='center', ha='center',
                style='italic', color='#888', fontsize=8)

    if save_path:
        plt.savefig(save_path, dpi=120, bbox_inches='tight')
    plt.close(fig)

    return {'snp': snp, 'k': k, 'Q': Q, 'I2': I2,
            'beta_fixed': beta_fixed, 'se_fixed': se_fixed, 'p_fixed': p_fixed,
            'beta_random': beta_random, 'se_random': se_random, 'p_random': p_random}


def extreme_studies(study_data, se_thresh=100) -> pd.DataFrame:
    """Per-study rows whose |SE| exceeds the threshold, across all SNPs loaded."""
    rows = []
    for study, df in study_data.items():
        for _, r in df[df['SE'].abs() > se_thresh].iterrows():
            rows.append({'SNP': r['ID'], 'study': study, 'A1': r['A1'],
                         'BETA': r['BETA'], 'SE': r['SE'], 'A1_FREQ': r['A1_FREQ'],
                         'OBS_CT': r['OBS_CT'], 'P': r['P']})
    return pd.DataFrame(rows)
