"""Manhattan and QQ plots for meta-analysis and per-study summary statistics."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import beta, chi2

GENOME_WIDE = 5e-8
SUGGESTIVE = 1e-6

# Alternating shades of one hue, so the Manhattan bands read as position
# rather than as categories.
BAND_COLORS = ['#2a78d6', '#0d366b']
LINE_COLOR = '#9a9993'
ALERT_COLOR = '#d03b3b'
TEXT_COLOR = '#52514e'


def clean_pvals(p) -> np.ndarray:
    """Finite p-values in (0, 1]; everything else is dropped."""
    p = np.asarray(p, dtype=float)
    return p[np.isfinite(p) & (p > 0) & (p <= 1)]


def genomic_lambda(p) -> float:
    """Genomic inflation factor: observed median chi-square (1 df) over the null."""
    return np.median(chi2.isf(p, 1)) / chi2.ppf(0.5, 1)


def manhattan_plot(df, title=None, thin_below=1.3, keep_frac=0.1, gap=2e7,
                   p_col='P', chr_col='CHR', bp_col='BP', snp_col='SNP',
                   label_hits=True, ax=None, save_path=None):
    """Manhattan plot of a PLINK-style meta-analysis table.

    The dense, non-significant cloud is thinned so the figure stays quick to
    draw and small on disk; every point above `thin_below` is kept. Pass
    `thin_below=None` to plot every variant.
    """
    df = df[df[p_col].notna() & (df[p_col] > 0)].copy()
    df['-log10P'] = -np.log10(df[p_col])

    if thin_below is not None:
        n_before = len(df)
        low = df['-log10P'] < thin_below
        df = pd.concat([df[~low], df[low].sample(frac=keep_frac, random_state=0)])
        print(f'plotting {len(df):,} of {n_before:,} variants '
              f'(all with -log10P >= {thin_below}, plus {keep_frac:.0%} of the rest)')
    else:
        print(f'plotting all {len(df):,} variants')

    # Lay the chromosomes end to end on a single x axis.
    df = df.sort_values([chr_col, bp_col])
    chroms = sorted(df[chr_col].unique())
    offset, offsets = 0.0, {}
    for c in chroms:
        offsets[c] = offset
        offset += df.loc[df[chr_col] == c, bp_col].max() + gap
    df['x'] = df[bp_col] + df[chr_col].map(offsets)

    if ax is None:
        _, ax = plt.subplots(figsize=(14, 5))

    ticks, ticklabels = [], []
    for i, c in enumerate(chroms):
        sub = df[df[chr_col] == c]
        ax.scatter(sub['x'], sub['-log10P'], s=3, alpha=0.8, linewidths=0,
                   color=BAND_COLORS[i % 2], rasterized=True)
        ticks.append((sub['x'].min() + sub['x'].max()) / 2)
        ticklabels.append(str(c))

    ax.axhline(-np.log10(GENOME_WIDE), color=ALERT_COLOR, linestyle='--', linewidth=1,
               label=r'genome-wide ($5\times10^{-8}$)')
    ax.axhline(-np.log10(SUGGESTIVE), color=LINE_COLOR, linestyle='--', linewidth=1,
               label=r'suggestive ($1\times10^{-6}$)')

    x_lo, x_hi = df['x'].min(), df['x'].max()
    if label_hits:
        # Label the lead variant on each chromosome reaching significance.
        for _, grp in df[df[p_col] < GENOME_WIDE].groupby(chr_col):
            lead = grp.loc[grp[p_col].idxmin()]
            rel = (lead['x'] - x_lo) / (x_hi - x_lo)
            ha = 'left' if rel < 0.06 else 'right' if rel > 0.94 else 'center'
            dx = 4 if ha == 'left' else -4 if ha == 'right' else 0
            ax.annotate(lead[snp_col], xy=(lead['x'], lead['-log10P']),
                        xytext=(dx, 7), textcoords='offset points',
                        ha=ha, fontsize=7, color=TEXT_COLOR)

    ax.set_xticks(ticks)
    ax.set_xticklabels(ticklabels, fontsize=8)
    ax.set_xlabel('Chromosome')
    ax.set_ylabel(r'$-\log_{10}(P)$')
    if title:
        ax.set_title(title)
    ax.set_xlim(x_lo - gap, x_hi + gap)
    ax.set_ylim(0, max(8.0, df['-log10P'].max() * 1.15))
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(axis='x', length=0)
    # Legend above the axes so it never sits on top of the data.
    ax.legend(frameon=False, fontsize=8, ncol=2,
              loc='lower right', bbox_to_anchor=(1.0, 1.0))

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300)
    return ax


def qq_plot(p, ax=None, title=None, color=BAND_COLORS[0], n_points=20000, band=True):
    """QQ plot of -log10(P) with the 95% null band and lambda annotated.

    Lambda is computed on every variant; only the drawn points are thinned
    (log-spaced, so the tail keeps every point and the bulk is sampled).
    Returns lambda, or NaN when the input has no usable p-values.
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(5, 5))

    p = clean_pvals(p)
    n = len(p)
    if n == 0:
        ax.text(0.5, 0.5, 'no usable\nP values', ha='center', va='center',
                fontsize=8, color=ALERT_COLOR, transform=ax.transAxes)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.spines[['top', 'right']].set_visible(False)
        if title:
            ax.set_title(title, fontsize=10)
        return np.nan

    lam = genomic_lambda(p)

    obs = -np.log10(np.sort(p))
    rank = np.arange(1, n + 1)
    exp = -np.log10((rank - 0.5) / n)

    idx = np.unique(np.geomspace(1, n, num=min(n_points, n)).astype(int)) - 1
    e, o, r = exp[idx], obs[idx], rank[idx]

    if band:   # 95% pointwise null band from the order statistics
        lo = -np.log10(beta.ppf(0.975, r, n - r + 1))
        hi = -np.log10(beta.ppf(0.025, r, n - r + 1))
        ax.fill_between(e, lo, hi, color='#d7d6d1', alpha=0.6, linewidth=0, zorder=1)

    lim = max(e.max(), o.max()) * 1.05
    ax.plot([0, lim], [0, lim], color=LINE_COLOR, linewidth=1, zorder=2)
    ax.scatter(e, o, s=6, color=color, alpha=0.85, linewidths=0, zorder=3, rasterized=True)

    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_aspect('equal')
    ax.spines[['top', 'right']].set_visible(False)
    ax.annotate(f'$\\lambda$ = {lam:.3f}\nn = {n:,}', xy=(0.04, 0.96),
                xycoords='axes fraction', va='top', ha='left',
                fontsize=8, color=TEXT_COLOR)
    if title:
        ax.set_title(title, fontsize=10)
    return lam


def qq_by_study(labelled_paths, title=None, ncols=5, p_col='P', save_path=None):
    """One QQ panel per input file, plus a table of genomic inflation factors.

    Reads only the p-value column of each file, one at a time, so this stays
    usable on the full per-chromosome result sets.
    """
    items = list(labelled_paths.items())
    nrows = int(np.ceil(len(items) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.1 * ncols, 3.1 * nrows))
    axes = np.atleast_1d(axes).ravel()

    rows = []
    for ax, (name, path) in zip(axes, items):
        p = pd.read_csv(path, sep='\t', usecols=[p_col])[p_col]
        n_usable = len(clean_pvals(p))
        lam = qq_plot(p, ax, title=name)
        rows.append({'study': name, 'lambda': lam,
                     'n_usable': n_usable, 'n_rows': len(p)})
        print(f'{name:28} lambda = {lam:6.3f}   usable = {n_usable:>9,} / {len(p):,}')

    for ax in axes[len(items):]:
        ax.axis('off')

    fig.supxlabel(r'Expected $-\log_{10}(P)$')
    fig.supylabel(r'Observed $-\log_{10}(P)$')
    if title:
        fig.suptitle(title, fontsize=12)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=200)

    lambdas = pd.DataFrame(rows).sort_values('lambda', ascending=False)

    # Flag studies that look broken rather than merely inflated.
    bad = lambdas[lambdas['lambda'].isna()
                  | (lambdas['lambda'] < 0.9) | (lambdas['lambda'] > 1.2)]
    if len(bad):
        print('\ncheck these before trusting the meta-analysis:')
        for _, r in bad.iterrows():
            why = ('no usable P values' if pd.isna(r['lambda'])
                   else 'deflated' if r['lambda'] < 0.9 else 'inflated')
            print(f"  {r['study']:28} lambda = {r['lambda']:6.3f}  ({why})")

    return lambdas
