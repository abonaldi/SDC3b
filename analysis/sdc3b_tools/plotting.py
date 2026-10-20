"""
Corner plots and comparison metrics for gridded SDC3b posteriors.

Ported from the sdc3b_results_comparison notebook. The grid coordinates are
passed in (from io.load_posterior) rather than assumed, and 1-D credible
intervals are read from the CDF at cell edges.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import pandas as pd
from scipy.linalg import inv

import corner
from getdist import plots, MCSamples

from .io import LABELS


class GridCornerPlotter:
    def __init__(self, coords, engine='corner', true_values=None, labels=None):
        """
        Initialize the plotter.
        
        Parameters:
        - coords: 1-D cell-centre coordinates shared by all axes (from io.load_posterior)
        - engine: 'corner' or 'getdist'
        - true_values: List of values [x1, x2, x3]
        - labels: List of LaTeX strings (no $ signs)
        """
        self.engine = engine.lower()
        self.true_values = true_values
        self.labels = labels or LABELS
        self.posteriors = [] 
        self.coords_1d = np.asarray(coords)
        
        # 30 visually distinct colours: tab20 (20) + first 10 of tab20b
        _cm1 = plt.get_cmap('tab20')
        _cm2 = plt.get_cmap('tab20b')
        self.fallback_cycle = [_cm1(i) for i in range(20)] + [_cm2(i) for i in range(10)]
        self.marker_cycle = ['o', 's', '^', 'v', 'D', 'p', 'h', '*', 'P', 'X', '<', '>', 'd', '8']
        
    def add_posterior(self, grid, label, color=None, marker=None):
        """Add a 3D probability grid with an optional color and marker."""
        self.posteriors.append({
            'grid': grid,
            'label': label,
            'color': color,
            'marker': marker
        })
        
    def _prepare_weighted_samples(self, grid):
        """Flatten 3D grid into weighted points."""
        X1, X2, X3 = np.meshgrid(self.coords_1d, self.coords_1d, self.coords_1d, indexing='ij')
        points = np.vstack([X1.flatten(), X2.flatten(), X3.flatten()]).T
        weights = grid.flatten()
        return points, weights

    def plot_corner(self, **kwargs):
        if not self.posteriors:
            raise ValueError("No posteriors added. Call .add_posterior() first.")
        if self.engine == 'corner':
            return self._plot_corner_engine(**kwargs)
        elif self.engine == 'getdist':
            return self._plot_getdist_engine(**kwargs)

    def _plot_corner_engine(self, **kwargs):
        formatted_labels = [f"${l}$" for l in self.labels]
        fig = None
        legend_handles = []
        
        for i, post in enumerate(self.posteriors):
            points, weights = self._prepare_weighted_samples(post['grid'])
            
            # Determine color from user input or the C0-C9 cycle
            current_color = post['color'] if post['color'] else self.fallback_cycle[i % len(self.fallback_cycle)]
            
            opts = {
                'labels': formatted_labels,
                'bins': 40,
                'plot_datapoints': False,
                'fill_contours': True,
                'levels': (0.68, 0.95),
                'range': [(0, 1), (0, 1), (0, 1)],
                'color': current_color,
                'smooth': 1.0,    # Added smoothing as requested
                'smooth1d': 1.0,  # Added 1D smoothing as requested
            }
            
            if fig is not None:
                opts['fig'] = fig
            
            opts.update(kwargs)
            fig = corner.corner(points, weights=weights, **opts)
            
            # Create legend handle for the team
            legend_handles.append(mlines.Line2D([], [], color=current_color, label=post['label'], lw=2))

        # corner fixes each 1-D panel's y-range from the first posterior; refit to all of them
        ndim = len(self.labels)
        for p in range(ndim):
            ax = fig.axes[p * (ndim + 1)]
            ymax = max((np.max(line.get_ydata()) for line in ax.get_lines()), default=1.0)
            ax.set_ylim(0, 1.1 * ymax)

        if self.true_values is not None:
            # Mark true values in gray (lines and points) without adding to legend
            corner.overplot_lines(fig, self.true_values, color="gray", linestyle="--", alpha=0.6)
            corner.overplot_points(fig, np.expand_dims(self.true_values, 0), 
                                   color="gray", marker="s", markersize=5)

        # Place legend in the top-right empty space of the corner plot
        fig.legend(handles=legend_handles, loc='upper right', bbox_to_anchor=(0.95, 0.95), fontsize=12)
            
        return fig

    def _plot_getdist_engine(self, **kwargs):
        samples_list = []
        names = ["x1", "x2", "x3"]
        for i, post in enumerate(self.posteriors):
            points, weights = self._prepare_weighted_samples(post['grid'])
            s = MCSamples(samples=points, weights=weights, names=names, 
                          labels=self.labels, label=post['label'], ranges={n: (0, 1) for n in names})
            samples_list.append(s)
            
        g = plots.get_subplot_plotter()
        
        # Sync colors with the same C0-C9 cycle for GetDist
        user_colors = [p['color'] for p in self.posteriors]
        final_colors = [user_colors[i] if user_colors[i] else self.fallback_cycle[i % len(self.fallback_cycle)] 
                        for i in range(len(user_colors))]

        g.triangle_plot(samples_list, filled=True, contour_colors=final_colors,
                        markers=self.true_values, marker_args={'color': 'gray', 'ls': '--', 'alpha': 0.6},
                        **kwargs)
        return g


class ComparisonPlotter(GridCornerPlotter):
    def __init__(self, coords, engine='corner', true_values=None, labels=None, true_ranges=None):
        super().__init__(coords, engine=engine, true_values=true_values, labels=labels)
        self.true_ranges = true_ranges  # list of [lo, hi] per parameter (defines the truth cube)
        if self.true_values is None:
            print('Warning: true_values not provided. Metrics will not be calculated.')

    def add_posterior(self, grid, label, color=None, marker=None):
        """Overrides base to calculate metrics immediately when data is added."""
        super().add_posterior(grid, label, color, marker)
        self.posteriors[-1]['stats'] = self._calculate_metrics(grid)

    def _marginal_credible_intervals(self, grid):
        """
        Asymmetric credible intervals from 1-D marginal posteriors.

        Marginalises the 3-D grid along the two non-target axes and reads
        percentiles from the CDF - no Gaussianity assumed. The CDF is sampled
        at cell edges, so each cell's mass is spread across its width.

        Returns a list of dicts (one per parameter p) with keys:
            median, lo1, hi1  -- equal-tailed 68% CI (16th / 84th percentile)
            lo2, hi2          -- equal-tailed 95% CI (2.5th / 97.5th percentile)
            cdf_at            -- callable: CDF evaluated at an arbitrary value
        """
        c = self.coords_1d
        dx = c[1] - c[0]
        x = np.append(c - dx / 2, c[-1] + dx / 2)
        cis = []
        for axis in range(3):
            other_axes = tuple(i for i in range(3) if i != axis)
            marginal = np.sum(grid, axis=other_axes)
            marginal = marginal / np.sum(marginal)
            cdf = np.append(0.0, np.cumsum(marginal))

            def quantile(q, _cdf=cdf, _x=x):
                return float(np.interp(q, _cdf, _x))

            def cdf_at(val, _cdf=cdf, _x=x):
                return float(np.interp(val, _x, _cdf))

            cis.append({
                'median': quantile(0.500),
                'lo1':    quantile(0.160),
                'hi1':    quantile(0.840),
                'lo2':    quantile(0.025),
                'hi2':    quantile(0.975),
                'cdf_at': cdf_at,
            })
        return cis

    def _compute_sdc3b_score(self, grid):
        """SDC3b official score: posterior mass within the true parameter cube.

        The cube is defined by self.true_ranges — the xHI range spanned by each
        frequency band (from xhi_true_min_max).  Higher score = better.
        """
        if self.true_ranges is None:
            return np.nan
        grid_norm = grid / np.sum(grid)
        x = self.coords_1d
        slices = []
        for lo, hi in self.true_ranges:
            i_lo = np.searchsorted(x, min(lo, hi), side='left')
            i_hi = np.searchsorted(x, max(lo, hi), side='right')
            slices.append(slice(i_lo, i_hi))
        return float(np.sum(grid_norm[tuple(slices)]))

    def _calculate_metrics(self, grid):
        """Calculate weighted stats, credible intervals, Z-scores, Mahalanobis, PIT, coverage, and SDC3b score."""
        X1, X2, X3 = np.meshgrid(self.coords_1d, self.coords_1d, self.coords_1d, indexing='ij')
        points  = np.vstack([X1.flatten(), X2.flatten(), X3.flatten()]).T
        weights = grid.flatten() / np.sum(grid)

        means  = np.average(points, axis=0, weights=weights)
        cov    = np.cov(points.T, aweights=weights)
        sigmas = np.sqrt(np.diag(cov))
        cis    = self._marginal_credible_intervals(grid)

        metrics = {'means': means, 'sigmas': sigmas, 'cov': cov, 'cis': cis}

        if self.true_values is not None:
            tv = np.asarray(self.true_values)
            metrics['z_scores']  = np.abs(means - tv) / sigmas

            delta = means - tv
            try:
                metrics['mahalanobis'] = float(np.sqrt(delta @ inv(cov) @ delta))
            except np.linalg.LinAlgError:
                metrics['mahalanobis'] = np.nan

            metrics['pit'] = np.array([cis[p]['cdf_at'](tv[p]) for p in range(len(tv))])

            metrics['cover_68'] = np.array([cis[p]['lo1'] <= tv[p] <= cis[p]['hi1'] for p in range(len(tv))])
            metrics['cover_95'] = np.array([cis[p]['lo2'] <= tv[p] <= cis[p]['hi2'] for p in range(len(tv))])

        metrics['sdc3b_score'] = self._compute_sdc3b_score(grid)

        return metrics

    def plot_forest(self):
        """
        Horizontal whisker plot of 68% and 95% credible intervals for each team.
        Centre marker = marginal median.  Use plot_colour_legend() to identify teams.
        """
        num_teams  = len(self.posteriors)
        num_params = len(self.labels)
        fig, axes  = plt.subplots(1, num_params, figsize=(14, num_teams * 0.4 + 2), sharey=True)
        if num_params == 1:
            axes = [axes]

        for p in range(num_params):
            ax = axes[p]
            for i, post in enumerate(self.posteriors):
                ci     = post['stats']['cis'][p]
                med    = ci['median']
                color  = post['color']  or self.fallback_cycle[i % len(self.fallback_cycle)]
                marker = post['marker'] or self.marker_cycle[i % len(self.marker_cycle)]
                ax.errorbar(med, i, xerr=[[med - ci['lo2']], [ci['hi2'] - med]],
                            fmt='none', color=color, lw=1, alpha=0.4)
                ax.errorbar(med, i, xerr=[[med - ci['lo1']], [ci['hi1'] - med]],
                            fmt=marker, color=color, lw=3, markersize=8)
            if self.true_values is not None:
                ax.axvline(self.true_values[p], color='gray', ls='--', alpha=0.7)
            ax.set_xlabel(f'${self.labels[p]}$', fontsize=15)
            if p == 0:
                ax.set_yticks(range(num_teams))
                ax.set_yticklabels([post['label'] for post in self.posteriors])
            ax.invert_yaxis()
        fig.tight_layout()
        return fig

    def plot_boxplot(self):
        """
        Horizontal box plot of 68% (box) and 95% (whiskers) credible intervals.
        The team's marker is placed at the marginal median, matching plot_forest.
        Use plot_colour_legend() to identify teams.
        """
        num_teams  = len(self.posteriors)
        num_params = len(self.labels)
        fig, axes  = plt.subplots(1, num_params, figsize=(14, num_teams * 0.4 + 2), sharey=True)
        if num_params == 1:
            axes = [axes]

        box_h = 0.45

        for p in range(num_params):
            ax = axes[p]
            for i, post in enumerate(self.posteriors):
                ci     = post['stats']['cis'][p]
                med    = ci['median']
                color  = post['color']  or self.fallback_cycle[i % len(self.fallback_cycle)]
                marker = post['marker'] or self.marker_cycle[i % len(self.marker_cycle)]

                # 95% whisker line
                ax.plot([ci['lo2'], ci['hi2']], [i, i], color=color, lw=3, alpha=0.8)
                # Whisker end caps
                for x_cap in [ci['lo2'], ci['hi2']]:
                    ax.plot([x_cap, x_cap], [i - box_h / 4, i + box_h / 4],
                            color=color, lw=3, alpha=0.8)

                # 68% filled box
                rect = plt.Rectangle(
                    (ci['lo1'], i - box_h / 2), ci['hi1'] - ci['lo1'], box_h,
                    facecolor=color, alpha=0.5, edgecolor=color, lw=1.5
                )
                ax.add_patch(rect)

                # Median marker (team's assigned marker shape)
                ax.plot(med, i, marker=marker, color=color, markersize=8, zorder=5,
                        markeredgecolor='k', markeredgewidth=0.5)

            if self.true_values is not None:
                ax.axvline(self.true_values[p], color='gray', ls='--', lw=1.5, alpha=1.0)
            if self.true_ranges is not None:
                ax.axvspan(self.true_ranges[p][0], self.true_ranges[p][1], color='gray', alpha=0.25)

            ax.set_xlabel(f'${self.labels[p]}$', fontsize=15)
            if p == 0:
                ax.set_yticks(range(num_teams))
                ax.set_yticklabels([post['label'] for post in self.posteriors])
            ax.invert_yaxis()
        fig.tight_layout()
        return fig

    def plot_precision_accuracy(self, markersize=80):
        """
        Scatter: Bias (median - truth) vs. 68% CI width for each parameter.
        Use plot_colour_legend() to identify teams.
        """
        num_params = len(self.labels)
        fig, axes  = plt.subplots(1, num_params, figsize=(15, 5))
        for p in range(num_params):
            ax = axes[p]
            for i, post in enumerate(self.posteriors):
                ci    = post['stats']['cis'][p]
                bias  = ci['median'] - self.true_values[p]
                width = ci['hi1'] - ci['lo1']
                color  = post['color']  or self.fallback_cycle[i % len(self.fallback_cycle)]
                marker = post['marker'] or self.marker_cycle[i % len(self.marker_cycle)]
                ax.scatter(bias, width, color=color, marker=marker, edgecolors='k', s=markersize, zorder=3)
            ax.axvline(0, color='gray', lw=1.5, ls='--')
            ax.set_xlabel(f'Bias  ($\\tilde{{\\theta}}_{p+1} - \\theta_{{\\mathrm{{truth}},{p+1}}}$)', fontsize=15)
            ax.set_ylabel(r'68% CI width  $\Delta_{68}$', fontsize=15)
            ax.set_title(f'${self.labels[p]}$', fontsize=15)
        fig.tight_layout()
        return fig

    # ----------------------------------------------------------------
    # Metric labels for axes
    _METRIC_LABELS = {
        'Z':           lambda p : r'$Z_'+str(p)+r' = |\tilde{x}_'+str(p)+r' - x_{\mathrm{truth},'+str(p)+r'}|\,/\,\sigma_'+str(p)+r'$',
        'PIT':         lambda p : r'$F_'+str(p)+r'(x_{\mathrm{truth},'+str(p)+r'})$',
        'Bias':        lambda p : r'$\tilde{x}_'+str(p)+r' - x_{\mathrm{truth},'+str(p)+r'}$',
        'CI68':        lambda p : r'$\Delta_{68,'+str(p)+r'}$',
        'Mahalanobis': lambda p : r'$D_M$',
        'Cover_68':    lambda p : r'Cover 68%',
        'Cover_95':    lambda p : r'Cover 95%',
        'Score':       lambda p : r'$S_{\mathrm{SDC3b}}$',
    }
    _PER_PARAM = {'Z', 'PIT', 'Bias', 'CI68'}

    _IDEAL_VALUES = {
        'Z':           0.0,
        'PIT':         0.5,
        'Bias':        0.0,
        'CI68':        0.0,   # narrower = better, but no strict ideal value
        'Mahalanobis': 1.0,
        'Cover_68':    3,
        'Cover_95':    3,
        'Score':       1.0,   # all posterior mass inside the truth cube = perfect
    }

    def _get_metric_value(self, post, metric, p=None):
        s = post['stats']
        tv = np.asarray(self.true_values)
        if metric == 'Z':           return s['z_scores'][p]
        if metric == 'PIT':         return s['pit'][p]
        if metric == 'Bias':        return s['cis'][p]['median'] - tv[p]
        if metric == 'CI68':        return s['cis'][p]['hi1'] - s['cis'][p]['lo1']
        if metric == 'Mahalanobis': return s.get('mahalanobis', np.nan)
        if metric == 'Cover_68':    return int(np.sum(s['cover_68']))
        if metric == 'Cover_95':    return int(np.sum(s['cover_95']))
        if metric == 'Score':       return s.get('sdc3b_score', np.nan)
        raise ValueError(f'Unknown metric: {metric!r}. '
                         f'Choose from {sorted(self._METRIC_LABELS)}')

    def plot_metric_combinations(self, metric1, metric2, markersize=80):
        """
        Scatter plot comparing two metrics across all teams.

        Per-parameter metrics (three panels, one per parameter):
            'Z', 'PIT', 'Bias', 'CI68'
        Scalar metrics (single panel):
            'Mahalanobis', 'Cover_68', 'Cover_95', 'Score'

        Use plot_colour_legend() to identify teams.

        Example:
            comp.plot_metric_combinations('Z', 'PIT')
            comp.plot_metric_combinations('Score', 'Mahalanobis')
        """
        per_param = (metric1 in self._PER_PARAM) or (metric2 in self._PER_PARAM)
        n_panels  = len(self.labels) if per_param else 1

        fig, axes = plt.subplots(1, n_panels, figsize=(5 * n_panels, 5))
        if n_panels == 1:
            axes = [axes]

        ideal1 = self._IDEAL_VALUES.get(metric1)
        ideal2 = self._IDEAL_VALUES.get(metric2)

        for panel in range(n_panels):
            ax = axes[panel]
            p1 = panel if metric1 in self._PER_PARAM else None
            p2 = panel if metric2 in self._PER_PARAM else None

            for i, post in enumerate(self.posteriors):
                x = self._get_metric_value(post, metric1, p1)
                y = self._get_metric_value(post, metric2, p2)
                color  = post['color']  or self.fallback_cycle[i % len(self.fallback_cycle)]
                marker = post['marker'] or self.marker_cycle[i % len(self.marker_cycle)]
                ax.scatter(x, y, color=color, marker=marker, edgecolors='k', s=markersize, zorder=3)

            if ideal1 is not None:
                ax.axvline(ideal1, color='gray', ls='--', lw=1.5, alpha=1.0, zorder=1)
            if ideal2 is not None:
                ax.axhline(ideal2, color='gray', ls='--', lw=1.5, alpha=1.0, zorder=1)
            # if ideal1 is not None and ideal2 is not None:
            #     ax.scatter(ideal1, ideal2, marker='*', s=markersize * 2, color='gold',
            #                edgecolors='k', lw=0.8, zorder=5, label='ideal')

            ax.set_xlabel(self._METRIC_LABELS.get(metric1, metric1)(panel+1), fontsize=15)
            ax.set_ylabel(self._METRIC_LABELS.get(metric2, metric2)(panel+1), fontsize=15)
            if per_param:
                ax.set_title(f'${self.labels[panel]}$', fontsize=15)

        fig.tight_layout()
        return fig

    def plot_colour_legend(self, ncol=2):
        """Standalone colour legend panel for team identification."""
        handles = [
            mlines.Line2D([], [],
                          color=post['color']  or self.fallback_cycle[i % len(self.fallback_cycle)],
                          marker=post['marker'] or self.marker_cycle[i % len(self.marker_cycle)],
                          linestyle='none', markersize=8, label=post['label'])
            for i, post in enumerate(self.posteriors)
        ]
        n     = len(handles)
        nrows = int(np.ceil(n / ncol))
        fig, ax = plt.subplots(figsize=(ncol * 2.5, nrows * 0.4 + 0.5))
        ax.legend(handles=handles, ncol=ncol, loc='center', frameon=False, fontsize=9)
        ax.axis('off')
        fig.tight_layout()
        return fig

    def get_metrics_table(self, sort_metric=None):
        """
        DataFrame of Z-scores, Mahalanobis distance, and SDC3b score for each team.
        Rows follow the order the posteriors were added (paper order); pass a
        column name as sort_metric to sort by it.
        """
        rows = []
        for post in self.posteriors:
            s   = post['stats']
            row = {'Team': post['label']}
            for p, z in enumerate(s.get('z_scores', [])):
                row[f'Z_{p+1}'] = round(z, 3)
            row['Mahalanobis'] = round(s.get('mahalanobis', np.nan), 3)
            row['Score']       = round(s.get('sdc3b_score', np.nan), 4)
            rows.append(row)
        df = pd.DataFrame(rows)
        if sort_metric in df.columns:
            df = df.sort_values(sort_metric)
        return df
