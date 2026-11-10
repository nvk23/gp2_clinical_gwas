"""Beta-beta comparison of two sets of summary statistics.

`compare_effect()` is the entry point: it loads two result tables, merges them
on SNP, classifies each variant by where it reaches significance, and draws the
scatter of effect sizes with concordance statistics.

Both tables must use the column names `SNP`, `BETA`, `SE`, `P` (and optionally
`A1`, `A2`, `EAF`); `rename_summary_stats()` maps the METAL/locuszoom columns
onto them.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.path import Path as MplPath
from matplotlib.transforms import Bbox
from scipy import stats

# METAL output merged with locuszoom annotation -> the names used here.
METAL_COLUMNS = {
    'rsid': 'SNP', 'beta': 'BETA', 'stderr_beta': 'SE', 'Pvalue': 'P',
    'alt': 'A1', 'ref': 'A2', 'alt_allele_freq': 'EAF',
}


def rename_summary_stats(df, mapping=None) -> pd.DataFrame:
    """Rename summary-statistic columns onto SNP / BETA / SE / P / A1 / A2 / EAF."""
    mapping = METAL_COLUMNS if mapping is None else mapping
    return df.rename(columns={k: v for k, v in mapping.items() if k in df.columns})


def report_duplicates(df, name: str):
    """Describe duplicated SNPs before they are dropped. Returns the dupe rows."""
    if 'SNP' not in df.columns:
        print(f'{name}: no SNP column')
        return None

    dupes = df[df['SNP'].duplicated(keep=False)].copy()
    if dupes.empty:
        print(f'{name}: no duplicates')
        return None

    counts = dupes['SNP'].value_counts()
    print(f'{name}: {len(dupes):,} duplicate rows over {dupes["SNP"].nunique():,} variants '
          f'({(counts == 2).sum():,} x2, {(counts == 3).sum():,} x3, '
          f'{(counts >= 4).sum():,} x4+)')

    cols = [c for c in ['SNP', 'CHR', 'BP', 'A1', 'A2', 'BETA', 'SE', 'P', 'EAF']
            if c in dupes.columns]
    for var_id, count in counts.head(5).items():
        var_data = dupes[dupes['SNP'] == var_id]
        print(f'\n{var_id} ({count} copies):')
        print(var_data[cols].to_string(index=False))
        for col, tol in (('P', 0.01), ('BETA', 0.01)):
            if col in var_data.columns:
                spread = var_data[col].max() - var_data[col].min()
                if spread > tol:
                    print(f'  {col} values differ (range {spread:.3g})')
    return dupes


def load_and_prep(df_or_path, name: str, snplist=None, drop_duplicates=True) -> pd.DataFrame:
    """Load one summary-statistics table, filter to `snplist`, drop duplicates.

    Duplicates are resolved by keeping the lowest p-value, which is the row a
    reader of the plot would quote.
    """
    if isinstance(df_or_path, (str, bytes)) or hasattr(df_or_path, '__fspath__'):
        try:
            df = pd.read_csv(df_or_path, sep='\t', low_memory=False)
        except Exception:
            df = pd.read_csv(df_or_path, sep=',', low_memory=False)
    else:
        df = df_or_path.copy()
    print(f'{name}: loaded {len(df):,} variants')

    if snplist is not None:
        if 'SNP' not in df.columns:
            raise KeyError(f'{name}: no SNP column, cannot filter to a SNP list')
        before = len(df)
        df = df[df['SNP'].isin(snplist)].copy()
        print(f'{name}: {len(df):,} of {before:,} variants are in the SNP list')
        missing = set(snplist) - set(df['SNP'])
        if missing:
            shown = sorted(missing)[:10]
            print(f'{name}: {len(missing)} requested SNPs not found'
                  + (f' ({", ".join(shown)})' if len(missing) <= 10 else ''))

    if drop_duplicates and 'SNP' in df.columns:
        before = len(df)
        sort_by = 'P' if 'P' in df.columns else 'SNP'
        df = df.sort_values(sort_by).drop_duplicates(subset='SNP', keep='first')
        if before - len(df):
            print(f'{name}: removed {before - len(df):,} duplicates (kept lowest P)')

    df['COHORT'] = name
    return df


def merge_for_comparison(df1, df2, label1: str, label2: str) -> pd.DataFrame:
    """Inner-join two summary tables on SNP, suffixing columns with their label."""
    if 'SNP' not in df1.columns or 'SNP' not in df2.columns:
        raise ValueError('both tables need a SNP column')

    merged = df1.merge(df2, on='SNP', how='inner', suffixes=('_1', '_2'))
    print(f'{len(merged):,} overlapping variants')

    rename = {f'{c}_{i}': f'{c}_{label}'
              for c in ('BETA', 'SE', 'P', 'A1', 'A2', 'EAF')
              for i, label in ((1, label1), (2, label2))}
    return merged.rename(columns={k: v for k, v in rename.items() if k in merged.columns})


def classify_significance(merged, label1, label2, sig_level=5e-8) -> pd.DataFrame:
    """Tag each variant with where it reaches `sig_level`: either, both, neither."""
    p1, p2 = f'P_{label1}', f'P_{label2}'
    if p1 not in merged.columns or p2 not in merged.columns:
        merged['category'] = 'None'
        return merged

    sig_1, sig_2 = merged[p1] < sig_level, merged[p2] < sig_level
    merged['category'] = 'None'
    merged.loc[sig_1 & sig_2, 'category'] = 'Both'
    merged.loc[sig_1 & ~sig_2, 'category'] = label1
    merged.loc[~sig_1 & sig_2, 'category'] = label2

    print('\nSignificance classification:')
    print(merged['category'].value_counts().to_string())
    return merged


def concordance(merged, label1, label2) -> dict | None:
    """Correlation, regression slope and sign agreement between the two betas."""
    b1, b2 = f'BETA_{label1}', f'BETA_{label2}'
    if b1 not in merged.columns or b2 not in merged.columns:
        return None

    valid = merged[[b1, b2]].dropna()
    if valid.empty:
        return None

    r, p = stats.pearsonr(valid[b1], valid[b2])
    slope, intercept, _, _, _ = stats.linregress(valid[b1], valid[b2])
    return {'r': r, 'r2': r ** 2, 'p_value': p, 'slope': slope, 'intercept': intercept,
            'sign_concordance': (np.sign(valid[b1]) == np.sign(valid[b2])).mean(),
            'n': len(valid)}


def most_discordant(merged, label1, label2, n=5) -> list[str]:
    """SNPs whose effect sizes differ most between the two analyses."""
    b1, b2 = f'BETA_{label1}', f'BETA_{label2}'
    if b1 not in merged.columns or b2 not in merged.columns:
        return []
    merged = merged.assign(beta_diff=np.abs(merged[b1] - merged[b2]))
    return merged.nlargest(n, 'beta_diff')['SNP'].tolist()


def place_labels(ax, x, y, names, segments, s=60, fontsize=8):
    """Label every point without overlapping markers, lines, labels or the axes.

    Greedy search over candidate positions around each point, nearest first;
    crowded points are placed before isolated ones.
    """
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    to_px = ax.transData.transform
    px_per_pt = fig.dpi / 72

    pts = to_px(np.column_stack([x, y]))
    r = np.sqrt(s) / 2 * px_per_pt + 2  # marker radius plus a small pad, in pixels
    marker_boxes = [Bbox([[px - r, py - r], [px + r, py + r]]) for px, py in pts]
    seg_paths = [MplPath(to_px(np.asarray(seg, dtype=float))) for seg in segments]
    ax_box = ax.get_window_extent(renderer)

    sizes = []
    for name in names:
        t = ax.text(0, 0, name, fontsize=fontsize)
        bb = t.get_window_extent(renderer)
        sizes.append((bb.width + 6, bb.height + 4))   # padding for the label's box
        t.remove()

    angles = np.deg2rad(np.arange(0, 360, 15))
    dists = [8, 14, 22, 32, 45, 60, 80, 105]  # points

    dist_px = np.linalg.norm(pts[:, None, :] - pts[None, :, :], axis=-1)
    order = np.argsort(-(dist_px < 80 * px_per_pt).sum(axis=1), kind='stable')

    placed = []
    for i in order:
        px, py = pts[i]
        w, h = sizes[i]
        best, best_cost = None, np.inf
        for d in dists:
            for a in angles:
                cx = px + np.cos(a) * (d * px_per_pt + w / 2)
                cy = py + np.sin(a) * (d * px_per_pt + h / 2)
                box = Bbox.from_bounds(cx - w / 2, cy - h / 2, w, h)

                cost = d
                if not (ax_box.x0 < box.x0 and box.x1 < ax_box.x1
                        and ax_box.y0 < box.y0 and box.y1 < ax_box.y1):
                    cost += 1e5
                cost += 1e3 * sum(box.overlaps(b) for b in marker_boxes + placed)
                cost += 1e3 * sum(p.intersects_bbox(box) for p in seg_paths)
                # Discourage leader lines that cross other markers
                leader = MplPath([[px, py], [cx, cy]])
                cost += 50 * sum(leader.intersects_bbox(b)
                                 for j, b in enumerate(marker_boxes) if j != i)
                if cost < best_cost:
                    best, best_cost = (cx, cy, box), cost
            if best_cost < 1e3:  # nearest clean spot found
                break

        cx, cy, box = best
        placed.append(box)
        ax.annotate(
            names[i], xy=(x[i], y[i]),
            xytext=((cx - px) / px_per_pt, (cy - py) / px_per_pt),
            textcoords='offset points', ha='center', va='center', fontsize=fontsize,
            bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='none', alpha=0.85),
            arrowprops=dict(arrowstyle='-', color='0.35', lw=0.6,
                            shrinkA=0, shrinkB=np.sqrt(s) / 2),
            zorder=10)


def beta_beta_plot(merged, label1, label2, labels=None, sig_level=5e-8,
                   annotate_snps=None, annotate_discordant=0, annotate_all=False,
                   plot_error_bars=True, color_map=None, alpha=0.85, s=60,
                   figsize=(12, 9), save=None, dpi=300):
    """Scatter of the two effect sizes, coloured by where each variant is significant."""
    b1, b2 = f'BETA_{label1}', f'BETA_{label2}'
    se1, se2 = f'SE_{label1}', f'SE_{label2}'

    if labels is None:
        labels = [label1, label2, 'Both', 'None']
    label_map = {label1: labels[0], label2: labels[1],
                 'Both': labels[2], 'None': labels[3]}
    legend_names = {labels[0]: f'{labels[0]} only', labels[1]: f'{labels[1]} only',
                    labels[2]: labels[2],
                    labels[3]: 'Neither' if labels[3] == 'None' else labels[3]}

    if color_map is None:
        color_map = {labels[0]: '#87CEEB', labels[1]: '#E64B35',
                     labels[2]: '#00A087', labels[3]: '#808080'}
    marker_map = {labels[0]: 'o', labels[1]: '^', labels[2]: 's', labels[3]: 'o'}

    merged = merged.copy()
    merged['display_label'] = merged['category'].map(label_map)

    fig, ax = plt.subplots(figsize=figsize)
    show_se = plot_error_bars and se1 in merged.columns and se2 in merged.columns
    valid = merged.dropna(subset=[b1, b2])

    # Symmetric limits shared by both axes so y = x is the true diagonal
    extent = np.abs(valid[[b1, b2]].values)
    if show_se:
        extent = np.concatenate([np.abs(valid[b1]) + valid[se1].fillna(0),
                                 np.abs(valid[b2]) + valid[se2].fillna(0)])
    lim = 1.15 * np.nanmax(extent)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)

    # Reference lines, kept as segments so the labels can route around them
    segments = [[(-lim, 0), (lim, 0)], [(0, -lim), (0, lim)], [(-lim, -lim), (lim, lim)]]
    ax.axhline(0, color='0.75', linewidth=0.8, zorder=1)
    ax.axvline(0, color='0.75', linewidth=0.8, zorder=1)
    ax.plot([-lim, lim], [-lim, lim], color='0.5', linewidth=1, linestyle=':',
            label='y = x', zorder=2)

    for cat in (labels[3], labels[1], labels[0], labels[2]):
        data = valid[valid['display_label'] == cat]
        if data.empty:
            continue
        if show_se:
            ax.errorbar(data[b1], data[b2], xerr=data[se1], yerr=data[se2],
                        fmt='none', ecolor=color_map[cat], alpha=0.6,
                        linewidth=0.9, zorder=2)
        ax.scatter(data[b1], data[b2], c=color_map[cat], marker=marker_map[cat], s=s,
                   label=f'{legend_names[cat]} (n={len(data):,})', alpha=alpha,
                   edgecolors='black', linewidth=0.6, zorder=3)

    if show_se:
        for _, row in valid.iterrows():
            bx, by = row[b1], row[b2]
            if pd.notna(row[se1]):
                segments.append([(bx - row[se1], by), (bx + row[se1], by)])
            if pd.notna(row[se2]):
                segments.append([(bx, by - row[se2]), (bx, by + row[se2])])

    metrics = concordance(merged, label1, label2)
    if metrics is not None:
        x_range = np.array([valid[b1].min(), valid[b1].max()])
        y_pred = metrics['slope'] * x_range + metrics['intercept']
        ax.plot(x_range, y_pred, 'k--', alpha=0.6, linewidth=1.5,
                label='Regression line', zorder=2)
        segments.append(list(zip(x_range, y_pred)))

        ax.text(1.03, 0.0,
                f"r = {metrics['r']:.3f} (P = {metrics['p_value']:.2g})\n"
                f"R² = {metrics['r2']:.3f}\n"
                f"Slope = {metrics['slope']:.3f}\n"
                f"Sign concordance = {metrics['sign_concordance']:.1%}\n"
                f"n = {metrics['n']:,}",
                transform=ax.transAxes, verticalalignment='bottom',
                horizontalalignment='left', fontsize=10,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    snps_to_annotate = []
    if annotate_all:
        snps_to_annotate += valid['SNP'].tolist()
    if annotate_snps:
        snps_to_annotate += list(annotate_snps)
    if annotate_discordant > 0:
        discordant = most_discordant(merged, label1, label2, n=annotate_discordant)
        snps_to_annotate += discordant
        print('\nMost discordant variants (by |BETA difference|):')
        for snp in discordant:
            row = merged[merged['SNP'] == snp].iloc[0]
            print(f'  {snp}: b1={row[b1]:.3f}, b2={row[b2]:.3f}, '
                  f'|d|={abs(row[b1] - row[b2]):.3f}')

    ax.set_xlabel(f'Effect size (β) - {label1}', fontsize=12, fontweight='bold')
    ax.set_ylabel(f'Effect size (β) - {label2}', fontsize=12, fontweight='bold')
    ax.set_title(f'{label1} vs {label2} effect sizes', fontsize=14,
                 fontweight='bold', pad=15)
    ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)

    # Legend outside the axes, categories before the reference lines
    handles, handle_labels = ax.get_legend_handles_labels()
    order = sorted(range(len(handles)), key=lambda k: isinstance(handles[k], plt.Line2D))
    ax.legend([handles[k] for k in order], [handle_labels[k] for k in order],
              loc='upper left', bbox_to_anchor=(1.02, 1.0),
              title=f'Significant (P < {sig_level:.0e}) in', frameon=True,
              fontsize=9, title_fontsize=9, framealpha=0.9)
    ax.set_aspect('equal', adjustable='box')
    plt.tight_layout()

    # Annotate after the layout is final, so the positions are exact
    snps_to_annotate = [snp for snp in dict.fromkeys(snps_to_annotate)
                        if snp in valid['SNP'].values]
    if snps_to_annotate:
        rows = valid.drop_duplicates('SNP').set_index('SNP').loc[snps_to_annotate]
        place_labels(ax, rows[b1].values, rows[b2].values, snps_to_annotate,
                     segments, s=s)

    if save:
        plt.savefig(save, dpi=dpi, bbox_inches='tight')
        print(f'\nsaved {save}')

    return fig, ax, metrics


def compare_effect(data1, data2, snplist=None, label=None, drop_duplicates=True,
                   sig_level=5e-8, investigate_dupes=False, return_data=False, **plot_kw):
    """Load, merge, classify and plot two sets of summary statistics.

    `data1` / `data2` are DataFrames or paths; `label` is
    [name1, name2, 'Both', 'None']. Extra keyword arguments go to
    `beta_beta_plot()`.
    """
    label = label or ['Cohort 1', 'Cohort 2', 'Both', 'None']
    name1, name2 = label[0], label[1]

    if investigate_dupes:
        for data, name in ((data1, name1), (data2, name2)):
            df = data if isinstance(data, pd.DataFrame) else pd.read_csv(
                data, sep='\t', low_memory=False)
            report_duplicates(df, name)

    df1 = load_and_prep(data1, name1, snplist=snplist, drop_duplicates=drop_duplicates)
    df2 = load_and_prep(data2, name2, snplist=snplist, drop_duplicates=drop_duplicates)

    merged = merge_for_comparison(df1, df2, name1, name2)
    merged = classify_significance(merged, name1, name2, sig_level=sig_level)

    fig, ax, metrics = beta_beta_plot(merged, name1, name2, labels=label,
                                      sig_level=sig_level, **plot_kw)
    return (fig, ax, metrics, merged) if return_data else (fig, ax, metrics)
