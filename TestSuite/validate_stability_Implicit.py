import sys
import numpy as np
import validation_tools as vt

ref_dir, data_dir = vt.processArgv(sys.argv[1:])

results = []

cond = True

# Check simulation setup
cond = cond and vt.check_setup('OUT_stdout', ref_dir, data_dir, 'truncation', 6)
cond = cond and vt.check_setup('OUT_stdout', ref_dir, data_dir, 'physical', 8)
cond = cond and vt.check_setup('OUT_stdout', ref_dir, data_dir, 'fixed_temp', 3)

# Check marginal.log
cond = cond and vt.check_setup('marginal.log', ref_dir, data_dir, 'converged to the bracket', 2)

if cond:
    print(vt.stability_success_str)
