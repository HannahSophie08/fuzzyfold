import fuzzyfold as ff
import numpy as np
import matplotlib.pyplot as plt 

seq = "AUAUUAGAUAUUAGUCAUAUGACUGACGGAAGUGGAGUUACCACAUGAAGUAUGACUAGGCAUAUUAUCUUAUAUGCCACAAAAA"

ssa = ff.Simulator(k0=1e5)
num_sims = 1000
structures = ssa.simulate_timecourse(seq, start=None, t_ext=0.02, t_end=0.02, t_lin=84, t_log=1, num_sims= num_sims)


#for t, counts in structures: 
    #for structure, n in counts.items():
        #print(f"t={t:.3f} {structure} {n}")

accessibility = []

for t, counts in structures: 
    struct_len = len(next(iter(counts)))
    counts_t = [0] * struct_len
    for structure, n in counts.items():
        for i, char in enumerate(structure): 
            if char == '.': 
                counts_t[i] += n
    
    accessibility.append([c / num_sims for c in counts_t])

times = [t for t, _ in structures]

max_len = max(len(row) for row in accessibility)
n_t = len(accessibility)

matrix = np.full((n_t, max_len), np.nan)
for i, row in enumerate(accessibility):
    matrix[i, :len(row)] = row
cmap = plt.cm.viridis.copy()
cmap.set_bad("white")
fig, ax = plt.subplots(figsize=(10, 8))
im = ax.imshow(matrix, aspect='auto', origin='upper', cmap=cmap, vmin=0, vmax=1)

ax.set_xlabel("sequence position")
ax.set_ylabel('transcript length')

fig.colorbar(im, ax=ax, label='accessibility')
path = '../../examples/pfl-riboswitch/accessibility_heatmap.png'
fig.savefig(path, bbox_inches='tight', dpi=300);  plt.close(fig);  print(f'Wrote {path}')



