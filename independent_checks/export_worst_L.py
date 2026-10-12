# Writes worst_states_L.json: the stage-L worst pendulum states (results/l1_polish_*.json) as plain Gaussian mixtures
# (weights, means, covariances). Only the parametrisation of the study code is used (h2_adversarial.unpack).
import json, sys, numpy as np
study = sys.argv[1] if len(sys.argv) > 1 else '../study'; sys.path.insert(0, study)
from h2_adversarial import unpack
out = []
for f in ['l1_polish_sigma0.01.json', 'l1_polish_sigma0.001.json']:
    for r in json.load(open(f'{study}/results/{f}')):
        ws, mus, Ss = unpack(np.array(r['x']), r['K'], r['sigma'])
        out.append(dict(sigma=r['sigma'], K=r['K'], gap=r['gap_search'], ws=[float(w) for w in ws],
                        mus=[[float(v) for v in m] for m in mus], Ss=[np.asarray(S).tolist() for S in Ss]))
json.dump(out, open('worst_states_L.json', 'w'))
print([(o['sigma'], o['K'], round(o['gap'], 6)) for o in out])
