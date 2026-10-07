"""Plot frozen scientific records; never import or run the scientific solver.

Run with Python, NumPy and Matplotlib. All arithmetic here is extraction,
unit scaling, error composition already specified by the scientific methods, or a
clearly labelled visual guide. Numerical allowances are not confidence intervals.
"""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import NullFormatter

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'figures'
OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'pdf.fonttype': 42, 'savefig.dpi': 240})
BLUE, ORANGE, GREEN, GRAY = '#23658a', '#bb5729', '#38795b', '#616971'
data = json.loads((ROOT/'scientific_data.json').read_text())
shifted = data['shifted_angular']
baseline = data['baseline_angular']
linear = data['linear_response']
endpoint = data['endpoint']
topology = data['branch_topology']
witness = data['perturbation']
settings = data['numerical_settings']
assert len(shifted['rows']) == len(baseline['rows']) == 93
assert len({tuple(r['configuration']) for r in shifted['rows']}) == 31
eep = witness['Eep']
eps = np.array([r['epsilon'] for r in shifted['complete_results']])
c = np.array([r['coefficient'] for r in shifted['complete_results']])
ei = np.array([r['EI'] for r in shifted['complete_results']])
ci, eli = linear['independent_coefficient'], linear['complete_physical_EI']
fine = settings['finest']
lookup = {(tuple(r['configuration']), r['epsilon']): r for r in shifted['rows']}
for e in eps:
    assert (tuple(fine), e) in lookup


def save(fig, name):
    fig.savefig(OUT / f'{name}.pdf', bbox_inches='tight')
    fig.savefig(OUT / f'{name}.png', bbox_inches='tight')
    plt.close(fig)


def decorate(ax):
    ax.grid(alpha=.16)
    ax.set_axisbelow(True)


# 1. Actual coexistence nodes and intervening stationary points; no fitted branch.
fig, axs = plt.subplots(1, 2, figsize=(9.4, 3.5))
tr = topology['rows']
t = np.array([r['fraction_Tcep'] for r in tr])
qlo = [min(r['minima'][0]['x'][2], r['minima'][1]['x'][2]) for r in tr]
qhi = [max(r['minima'][0]['x'][2], r['minima'][1]['x'][2]) for r in tr]
qmid = [r['intervening']['x'][2] for r in tr]
axs[0].scatter(t, qlo, color=BLUE, label='Lower-q uniform minimum', s=32)
axs[0].scatter(t, qhi, color=ORANGE, label='Upper-q uniform minimum', s=32)
axs[0].scatter(t, qmid, color=GRAY, marker='x', label='Intervening stationary point', s=32)
axs[0].scatter([1], [endpoint['fine']['x'][2]], color=GREEN, marker='*', s=110, label='Uniform endpoint')
axs[0].set(xlabel=r'$T/T_{\mathrm{CEP}}$', ylabel=r'$q$', title='(a) Local uniform branch coalescence')
axs[0].legend(frameon=False, fontsize=7.6)
axs[1].scatter(1-t, [r['barrier_per_N0'] for r in tr], color=BLUE, s=40)
axs[1].set(xscale='log', yscale='log', xlabel=r'$1-T/T_{\mathrm{CEP}}$',
           ylabel=r'Uniform barrier $\Delta F/(N_0 T_{c0}^2)$', title='(b) Measured coexistence barriers')
axs[1].set_xticks([.02,.05,.1,.2], ['0.02','0.05','0.10','0.20'])
axs[1].xaxis.set_minor_formatter(NullFormatter())
for ax in axs: decorate(ax)
fig.tight_layout(w_pad=2)
save(fig, 'fig01_endpoint_topology')

# 2. Independent linear qualification: actual frequency refinement at finest angle.
lr = sorted([r for r in linear['rows'] if r['Ntheta'] == 2048], key=lambda r:r['Nw'])
nw = np.array([r['Nw'] for r in lr])
fig, axs = plt.subplots(1, 2, figsize=(9.4, 3.5))
axs[0].plot(nw, [(r['coefficient']-ci)*1e9 for r in lr], 'o-', color=BLUE,
            label=r'$N_\theta=2048$; lines guide the eye')
axs[0].set(xscale='log', xlabel=r'Matsubara cutoff $N_\omega$',
           ylabel=r'$(c(N_\omega)-c_I)\times10^9$', title='(a) Resolved linear coefficient')
axs[0].axhline(0, color=GRAY, lw=.8)
axs[0].legend(frameon=False, fontsize=7.5)
for key, label, color in [('numerical_error','Total row numerical allowance', BLUE),
                         ('tail_remainder','Tail allowance', ORANGE),
                         ('primitive_error','Empirical arithmetic allowance', GREEN)]:
    values = [r.get(key, r.get('tail', {}).get('remainder')) for r in lr]
    assert all(v is not None for v in values), key
    axs[1].plot(nw, values, 'o-', color=color, label=label)
axs[1].axhline(eep, color=GRAY, ls='--', label=r'Source/endpoint allowance $E_{\mathrm{ep}}$')
axs[1].set(xscale='log', yscale='log', xlabel=r'$N_\omega$', ylabel='Coefficient allowance',
           title='(b) Tail refinement and retained source floor')
axs[1].legend(frameon=False, fontsize=7.5)
for ax in axs: decorate(ax)
fig.tight_layout(w_pad=2)
save(fig, 'fig02_linear_qualification')

# 3. The abscissa epsilon squared makes the leading finite-amplitude correction visible.
fig, ax = plt.subplots(figsize=(6.5, 4.0))
x = eps**2 * 1e5
ax.axhline(0, color=GRAY, lw=.8)
ax.axhspan((ci-eli-eep)*1e5, (ci+eli+eep)*1e5, color=GREEN, alpha=.12,
           label=r'Linear $c_I\pm(E_I+E_{\mathrm{ep}})$')
ax.axhline(ci*1e5, color=GREEN, ls='--', label='Independent quadratic limit')
ax.errorbar(x, c*1e5, yerr=(ei+eep)*1e5, fmt='o', capsize=4, color=BLUE,
            label=r'Shifted origin: $C(\epsilon)\pm(E_I(\epsilon)+E_{\mathrm{ep}})$')
ax.errorbar(x, c*1e5, yerr=ei*1e5, fmt='none', capsize=2, color=BLUE, lw=2)
cb = np.array([r['coefficient'] for r in baseline['complete_results']])
ax.scatter(x, cb*1e5, marker='x', color=ORANGE, zorder=4, label='Baseline angular origin')
ax.set(xlabel=r'$\epsilon^2\times10^5$', ylabel=r'$C(\epsilon)\times10^5$',
       title='Finite-amplitude quotient approaches a negative limit')
ax.legend(frameon=False, fontsize=8, loc='best'); decorate(ax)
fig.tight_layout(); save(fig, 'fig03_nonlinear_quotient')

# 4. Quotient times epsilon squared: purely algebraic reconstruction of observed descent.
fig, ax = plt.subplots(figsize=(6.5, 4.0))
ax.axhline(0, color=GRAY, lw=.8)
ax.errorbar(eps*1e3, eps**2*c*1e9, yerr=eps**2*(ei+eep)*1e9,
            fmt='o', capsize=4, color=BLUE, label='Mean finite-amplitude energy change')
guide = np.linspace(0, eps.max(), 200)
ax.plot(guide*1e3, guide**2*ci*1e9, '--', color=GREEN,
        label=r'Quadratic guide $c_I\epsilon^2$ (not measurements)')
ax.set(xlabel=r'$\epsilon\times10^3$', ylabel=r'Mean $\Delta F/(N_0 T_{c0}^2)\times10^9$',
       title='Finite perturbations lower the free energy')
ax.text(.03, .04, 'Bars: numerical + source/endpoint allowances; no statistical CI',
        transform=ax.transAxes, fontsize=7.7)
ax.legend(frameon=False, fontsize=8); decorate(ax)
fig.tight_layout(); save(fig, 'fig04_energy_descent')

# 5. Actual deterministic origin differences versus recorded comparison allowance.
origin = shifted['origin_comparisons']
fig, axs = plt.subplots(1, 2, figsize=(9.4, 3.5))
pos = np.arange(3)
axs[0].bar(pos-.17, [r['difference'] for r in origin], .34, color=BLUE, label='Absolute origin difference')
axs[0].bar(pos+.17, [r['allowed'] for r in origin], .34, color=ORANGE, label=r'Allowed $E_I^{\mathrm{shift}}+E_I^{\mathrm{base}}$')
axs[0].set(yscale='log', ylabel='Coefficient difference / allowance', title='(a) Shifted-origin agreement')
axs[0].legend(frameon=False, fontsize=7.5)
ratios = np.array([r['difference']/r['allowed'] for r in origin])
axs[1].bar(pos, ratios, color=BLUE)
axs[1].set(yscale='log', ylim=(1e-5, 1.8), ylabel='Difference / allowed difference',
           title='(b) Agreement far inside the comparison gate')
axs[1].axhline(1, color=ORANGE, ls='--', label='Comparison threshold')
axs[1].legend(frameon=False, fontsize=7.5)
for ax in axs:
    ax.set_xticks(pos, [f'{e:g}' for e in eps]); ax.set_xlabel(r'$\epsilon$'); decorate(ax)
fig.tight_layout(w_pad=2); save(fig, 'fig05_origin_agreement')

# 6. Five empirical axes and all ten mixed checks, from the saved 93-row schedule.
fig, axs = plt.subplots(2, 3, figsize=(10.2, 6.2))
labels = [r'Profile degree $J$', r'Matsubara cutoff $N_\omega$',
          r'Angular nodes $N_\theta$', r'Force half-degree $P$', r'Path nodes $N_t$']
colors = [ORANGE, BLUE, GREEN]
for a, ax in enumerate(axs.flat[:5]):
    for e, color in zip(eps, colors):
        levels = settings['levels'][a]; vals = []
        for level in levels:
            cfg = fine.copy(); cfg[a] = level
            vals.append((lookup[(tuple(cfg), e)]['coefficient']-lookup[(tuple(fine), e)]['coefficient'])*1e6)
        ax.plot(levels, vals, 'o-', color=color, label=rf'$\epsilon={e:g}$', ms=4)
    ax.set(xlabel=labels[a], ylabel=r'$(C-C_{\mathrm{fine}})\times10^6$',
           title=f'({chr(97+a)}) Empirical refinement')
    ax.set_xticks(settings['levels'][a]); ax.ticklabel_format(axis='y', style='sci', scilimits=(-3, 3)); decorate(ax)
    if a == 0: ax.legend(frameon=False, fontsize=7)
pairs = settings['mixed_pairs']; ratios = []
for h in shifted['hierarchy']:
    col = []
    for m in h['mixed']:
        d0,d1 = map(abs,m['contrasts']); n0,n1 = m['noise']
        bound1 = .75*d0+n0+n1
        bound2 = .25*sum(a['error'] for a in h['axes'])+n1
        col.append(max(d1/bound1, d1/bound2))
    ratios.append(col)
mixed_ratios = np.array(ratios).T
ax = axs.flat[5]
im=ax.imshow(np.maximum(mixed_ratios,1e-12), cmap='viridis',
             norm=LogNorm(vmin=1e-12, vmax=1), aspect='auto')
for i in range(10):
    for j in range(3):
        if mixed_ratios[i,j] == 0:
            ax.text(j,i,'0',color='white',ha='center',va='center',fontsize=6)
ax.set_xticks(range(3), [f'{e:g}' for e in eps]); ax.set_yticks(range(10), [f'{a+1},{b+1}' for a,b in pairs])
ax.set(xlabel=r'$\epsilon$', ylabel='Axis pair (panels a–e = 1–5)', title='(f) Mixed gate utilization')
fig.colorbar(im, ax=ax, label='Gate ratio (color floor $10^{-12}$; 0 marked)', fraction=.08)
fig.tight_layout(w_pad=1.4, h_pad=1.5); save(fig, 'fig06_axis_mixed_convergence')

# 7. Exact sign tests: Eep is retained once, not hidden in a numerical-only bar.
rich = shifted['epsilon_reduction']
labels = ['Primary', 'Independent\nlinear', r'$\epsilon=0.005$', r'$\epsilon=0.0025$', 'Richardson\n12']
centres = np.array([witness['cP'], ci, c[1], c[2], rich['Richardson12']])
numerical = np.array([witness['EP'], eli, ei[1], ei[2], rich['complete_Richardson_error']])
source = np.full(5, eep)
upper = centres + numerical + source
assert np.all(upper < 0)
fig, axs = plt.subplots(1, 2, figsize=(9.8, 3.8))
y = np.arange(5)
axs[0].barh(y, numerical*1e6, color=BLUE, label=r'$E_P$, $E_I$, or Richardson allowance')
axs[0].barh(y, source*1e6, left=numerical*1e6, color=ORANGE, label=r'Retained $E_{\mathrm{ep}}$')
axs[0].scatter(abs(centres)*1e6, y, color=GREEN, marker='|', s=130, label='Available negative magnitude')
axs[0].set_yticks(y, labels); axs[0].invert_yaxis()
axs[0].set(xlabel=r'Coefficient magnitude $\times10^6$', title='(a) Complete error budget versus negative signal')
axs[0].legend(frameon=False, fontsize=7.2, loc='lower right')
axs[1].errorbar(centres*1e5, y, xerr=(numerical+source)*1e5,
               fmt='o', color=BLUE, capsize=4)
axs[1].scatter(upper*1e5,y, color=ORANGE, marker='s', s=23, label='Conservative upper endpoint')
axs[1].axvline(0,color=GRAY,lw=.9); axs[1].set_yticks(y, labels); axs[1].invert_yaxis()
axs[1].set(xlabel=r'Coefficient with complete allowance $\times10^5$', title='(b) Every required upper endpoint remains negative')
axs[1].legend(frameon=False, fontsize=7.2, loc='lower right')
for ax in axs: decorate(ax)
fig.tight_layout(w_pad=1.7); save(fig, 'fig07_error_margin')

# Small extracted data, deliberately excluding checkpoint histories and raw primitives.

print('Generated seven scientific figures from scientific_data.json.')
