import sys
import numpy as np
import validation_tools as vt

ref_dir, data_dir = vt.processArgv(sys.argv[1:])

results = []

physicals = ['rayleigh', 'prandtl', 'velocity', 'temperature']

out = vt.checkHdf5("visState0000.hdf5", ref_dir, data_dir, physicals, 0)
results.append(out[0])

tids = [0, 1, 2, 3]
tols = [100, 176, 564, 564]
datasets = ['temperature/temperature', 'velocity/velocity_r', 'velocity/velocity_theta', 'velocity/velocity_phi']
for tid, t, ds in  zip(tids, tols, datasets):
    results.append(vt.hdf5Test(*out[1], ds, tid, tol = t, perrow = True, percol = True, threshold = 2e-14))

# Output test summary
vt.printSummary(results, tids, reftol = tols)
