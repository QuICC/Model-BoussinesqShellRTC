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

# Check ULP of Rac
def checkRa(dlines, rlines):
    msg = f'Checking error tolerance'
    # Set reference Rac
    t = rlines[1].split(" ")
    Rac_ref = (float(t[0])+float(t[-1]))/2.
    # Set data Rac
    t = dlines[1].split(" ")
    Rac_dat = (float(t[0])+float(t[-1].strip()))/2.
    # Compute ulp
    ulp = vt.compute_ulp(Rac_ref)
    err_ulp = int(np.abs(Rac_dat - Rac_ref)/ulp)

    cond = err_ulp < tol
    vt.printResult(cond, msg)
    if not cond:
        details = f'(Ra: {Rac_dat}, tol: {tol:.0f}, ulp: {err_ulp:.3e})'
        print(f'\t{details}')

    return cond

tol = 1e6
# Check marginal.log 
cond = cond and vt.check_setup('marginal.log', ref_dir, data_dir, 'converged to the bracket', 2, checkRa)

if cond:
    print(vt.stability_success_str)
