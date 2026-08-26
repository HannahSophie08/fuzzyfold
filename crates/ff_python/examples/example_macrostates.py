import fuzzyfold as ff
import argparse, re
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

seq = "AUAUUAGAUAUUAGUCAUAUGACUGACGGAAGUGGAGUUACCACAUGAAGUAUGACUAGGCAUAUUAUCUUAUAUGCCACAAAAA"

ssa = ff.Simulator(k0=1e5)

#results = ssa.simulate_timecourse(seq, None, t_ext=40, t_end=40, num_sims=50)

occupancy = ssa.simulate_macrostates(
    seq, 
    t_ext=0.02,
    t_end=30,
    t_lin=168,
    t_log=1000,
    num_sims=100,
    macrostates=["../../examples/pfl-riboswitch/pfl-IH1.ms", "../../examples/pfl-riboswitch/pfl-IH1_P2.ms", "../../examples/pfl-riboswitch/pfl-P1_P2.ms", "../../examples/pfl-riboswitch/pfl-P1_P2_linker.ms" ], 
)

def linlog_x(t, t_split, t_end, frac=0.7):
    t = np.asarray(t, dtype=float)
    lin = t / t_split * frac
    log = frac + np.log(np.maximum(t, t_split) / t_split) / np.log(t_end / t_split) * (1 - frac)
    return np.where(t <= t_split, lin, log)

def linlog_ticks(t_split, t_end, frac=0.7, n_lin=5, n_log=5):
    lin = np.linspace(0, t_split, n_lin + 1)
    log = np.geomspace(t_split, t_end, n_log + 1)[1:]  # geometric spacing, skip t_split
    ticks_t = np.concatenate([lin, log])
    ticks_x = linlog_x(ticks_t, t_split, t_end, frac)
    labels   = [f'{t:.3g}' for t in ticks_t]
    return ticks_x, labels


def plot_occupancy(occupancy, t_split, t_end, title, path):

    times = [t for t, _ in occupancy]

    fig, ax = plt.subplots(figsize=(14, 5.5))
    fig.suptitle(title, fontsize=15)
    t_sim = t_split + t_end 
    x = linlog_x(times, t_split, t_sim)

    macrostates = occupancy[0][1].keys()   # macrostate names, same keys at every timepoint

    for macrostate in macrostates:
        y = [fractions[macrostate] for _, fractions in occupancy]
        ax.plot(x, y, label=macrostate)

    tick_x, tick_labels = linlog_ticks(t_split, t_sim)
    ax.set_xticks(tick_x)
    ax.set_xticklabels(tick_labels, fontsize=10)
    ax.set(xlim=(0, 1), ylim=(0, 1), xlabel='', ylabel='occupancy')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_color('lightgrey')
    ax.spines['right'].set_linewidth(0.5)
    ax.grid(True, color='lightgrey', lw=0.5)

    for (x0, x1), lbl in [((0.0, 0.7), 'linear time [s]'), ((0.7, 1.0), 'logarithmic time [s]')]:
        mid = (x0 + x1) / 2
        kw = dict(transform=ax.transAxes, color='k', lw=1, clip_on=False)
        ax.plot([x0, x1], [-0.12, -0.12], **kw)          # horizontal bar
        ax.plot([x0, x0], [-0.12, -0.09], **kw)           # left serif
        ax.plot([x1, x1], [-0.12, -0.09], **kw)           # right serif
        ax.text(mid, -0.16, lbl, ha='center', va='top',
                transform=ax.transAxes, fontsize=10)

    # Subtle marker at the linear/log boundary
    ax.axvline(0.7, color='k', alpha=0.5, lw=1)

    ax.legend(fontsize=9, framealpha=0.85)
    fig.savefig(path, bbox_inches='tight', dpi=300);  plt.close(fig);  print(f'Wrote {path}')


t_split = (len(seq) - 1) * 0.02
plot_occupancy(occupancy, t_split, 30, "Macro-state occupancy", "../../examples/pfl-riboswitch/pfl_occupancy" )