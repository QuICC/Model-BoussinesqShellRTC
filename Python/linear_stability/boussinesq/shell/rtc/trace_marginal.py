"""Script to run a marginal curve trace for the Boussinesq rotating thermal convection in a spherical shell (Toroidal/Poloidal formulation)"""

# Command line options used: -st_type sinvert -eps_target 0 -eps_target_real

import numpy as np

import quicc_solver.model.boussinesq.shell.rtc.implicit.physical_model as mod
import quicc.linear_stability.marginal_curve as MarginalCurve

# Create the model and activate linearization
model = mod.PhysicalModel()
model.linearize = True
model.use_galerkin = False

# Set resolution, parameters, boundary conditions
Rac = None
mc = None
rescale = False

# SF/SF, FT/FT, differential heating, small gap
#bc_vel = 1; bc_temp = 0; heating = 1; rratio = 0.99; Pr = 1; rescale = True
#res = [64, 92, 0]
#Ta = 1e12; Ek = Ta**-0.5; Rac = 7.4758215298779; mc = 245

# SF/SF, FT/FT, differential heating
#bc_vel = 1; bc_temp = 0; heating = 1; rratio = 0.35; Pr = 1; rescale = True
#res = [96, 96, 0]
#Ta = 1e10; Ek = Ta**-0.5; Rac = 40; mc = 10

# SF/SF, FT/FT, internal heating
#bc_vel = 1; bc_temp = 0; heating = 0; rratio = 0.35; Pr = 1; rescale = True
#Ta = 1e6; Ek = Ta**-0.5
#res = [32, 32, 0]
#Ta = 1e7; Ek = Ta**-0.5
#res = [32, 32, 0]
#Ta = 1e8; Ek = Ta**-0.5; Rac = 31.957658460641; mc = 7
#res = [48, 48, 0]
#Ta = 1e9; Ek = Ta**-0.5; Rac = 43.034758690274; mc = 10
#res = [48, 48, 0]
#Ta = 1e10; Ek = Ta**-0.5; Rac = 59.975071459666; mc = 13
#res = [64, 64, 0]
#Ta = 1e11; Ek = Ta**-0.5; Rac = 85.363356944817; mc = 20
#res = [64, 64, 0]
#Ta = 1e12; Ek = Ta**-0.5; Rac = 122.69214718393; mc = 30
#res = [128, 128, 0]
#Ta = 1e13; Ek = Ta**-0.5; Rac = 177.55422348123; mc = 44
#res = [192, 192, 0]
#Ta = 1e14; Ek = Ta**-0.5; Rac = 258.13410447601; mc = 65
#res = [512, 256, 0]
#Ta = 1e15; Ek = Ta**-0.5; Rac = 376.44742717745; mc = 95
#res = [512, 384, 0]
#Ta = 1e16; Ek = Ta**-0.5
#res = [768, 512, 0]
#Ta = 1e17; Ek = Ta**-0.5
#res = [768, 512, 0]
#Ta = 1e18; Ek = Ta**-0.5
#res = [1024, 384, 0]

# NS/NS, FT/FT, internal heating
bc_vel = 0; bc_temp = 0; heating = 0; rratio = 0.35; Pr = 1; rescale = True
#Ta = 1e6; Ek = Ta**-0.5
#res = [32, 32, 0]
#Ta = 1e7; Ek = Ta**-0.5
#res = [32, 32, 0]
Ta = 1e8; Ek = Ta**-0.5; Rac = 31.534088376364; mc = 6
res = [48, 48, 0]
#Ta = 1e9; Ek = Ta**-0.5; Rac = 42.219154540505; mc = 9
#res = [48, 48, 0]
#Ta = 1e10; Ek = Ta**-0.5; Rac = 59.124359856967; mc = 13
#res = [96, 96, 0]
#Ta = 1e11; Ek = Ta**-0.5; Rac = 84.487326687693; mc = 20
#res = [128, 128, 0]
#Ta = 1e12; Ek = Ta**-0.5; Rac = 121.87739395205; mc = 30
#res = [192, 192, 0]
#Ta = 1e13; Ek = Ta**-0.5; Rac = 176.79656879674; mc = 44
#res = [256, 256, 0]
#Ta = 1e14; Ek = Ta**-0.5; Rac = 257.45628575047; mc = 65
#res = [384, 256, 0]
#Ta = 1e15; Ek = Ta**-0.5; Rac = 375.86277729259; mc = 95
#res = [512, 384, 0]
#Ta = 1e16; Ek = Ta**-0.5
#res = [768, 512, 0]
#Ta = 1e17; Ek = Ta**-0.5
#res = [768, 768, 0]
#Ta = 1e18; Ek = Ta**-0.5
#res = [1024, 1024, 0]

# NS/NS, FT/FT, differential heating
#bc_vel = 0; bc_temp = 0; heating = 1; rratio = 0.35; Pr = 1; rescale = True
#Ta = 1e6; Ek = Ta**-0.5
#res = [32, 32, 0]
#Ta = 1e7; Ek = Ta**-0.5
#res = [32, 32, 0]
#Ta = 1e8; Ek = Ta**-0.5; Rac = 28.93487228774; mc = 5
#res = [48, 48, 0]
#Ta = 1e9; Ek = Ta**-0.5; Rac = 32.817417129811; mc = 7
#res = [48, 48, 0]
#Ta = 1e10; Ek = Ta**-0.5; Rac = 39.148927020559; mc = 9
#res = [96, 96, 0]
#Ta = 1e11; Ek = Ta**-0.5; Rac = 84.487326687693; mc = 20
#res = [128, 128, 0]
#Ta = 1e12; Ek = Ta**-0.5; Rac = 105; mc = 25
#res = [192, 192, 0]
#Ta = 1e13; Ek = Ta**-0.5; Rac = 176.79656879674; mc = 44
#res = [256, 256, 0]
#Ta = 1e14; Ek = Ta**-0.5; Rac = 257.45628575047; mc = 65
#res = [384, 256, 0]
#Ta = 1e15; Ek = Ta**-0.5; Rac = 375.86277729259; mc = 95
#res = [512, 384, 0]
#Ta = 1e16; Ek = Ta**-0.5
#res = [768, 512, 0]
#Ta = 1e17; Ek = Ta**-0.5
#res = [768, 768, 0]
#Ta = 1e18; Ek = Ta**-0.5
#res = [1024, 1024, 0]

# NS/NS, FT/FT, internal heating, ri/ro =  0.4
#bc_vel = 0; bc_temp = 0; heating = 0; rratio = 0.4; Pr = 1; rescale = False
#Ek = 3e-3; Rac = 27.254; mc = 4
#res = [32, 32, 0]
#Ek = 1e-3; Rac = 23.7692; mc = 5
#res = [32, 32, 0]
#Ek = 5e-4; Rac = 25.6031; mc = 6
#res = [64, 64, 0]
#Ek = 1e-4; Rac = 36.4955; mc = 10
#res = [48, 48, 0]
#Ek = 1e-5; Rac = 71.55899; mc = 19
#res = [64, 64, 0]

# NS/NS, FT/FT, internal heating, ri/ro =  0.2
#bc_vel = 0; bc_temp = 0; heating = 0; rratio = 0.2; Pr = 1; rescale = False
#Ek = 1e-4; Rac = 46.482; mc = 6

# NS/NS, FT/FT, internal heating, ri/ro =  0.6
#bc_vel = 0; bc_temp = 0; heating = 0; rratio = 0.6; Pr = 1; rescale = False
#Ek = 1e-4; Rac = 27.0031; mc = 18

# Convert Ekman to Taylor number
Ta = Ek**-2

# Create parameters (rescaling to proper nondimensionalisation)
ro = model.automatic_parameters({'ekman':Ek, 'r_ratio':rratio})['upper1d']
if mc is None:
    m = np.int_(0.3029*Ta**(1./6.)) # Asymptotic prediction for minimum
else:
    m = mc
if Rac is None:
    Ra = (4.1173*Ta**(2./3.) + 17.7815*Ta**(0.5))*(1.0+rratio)/(2.0*ro**2*Ta**0.5) # Asymptotic prediction for critical Rayleigh number
else:
    Ra = Rac

res = [res[0], res[1]+m, 0] # Extend harmonic degree by harmonic order (fixed number of modes)
if rescale:
    Ek = (Ta*(1.0-rratio)**4)**-0.5
else:
    Ek = Ta**-0.5
eq_params = {'ekman':Ek, 'prandtl':Pr, 'rayleigh':Ra, 'r_ratio':rratio, 'heating':heating}
eq_params.update(model.automatic_parameters(eq_params))
bcs = {'bcType':model.SOLVER_HAS_BC, 'velocity':bc_vel, 'temperature':bc_temp}

# Wave number function from single "index" (k perpendicular)
def wave(m):
    return [float(m)]

eigs = wave(1)

# Collect GEVP setup parameters into single dictionary
gevp_opts = {'model':model, 'res':res, 'eq_params':eq_params, 'eigs':eigs, 'bcs':bcs, 'wave':wave}

# Setup computation, visualization and IO
marginal_options = MarginalCurve.default_options()
marginal_options['evp_tol'] = 1e-10
marginal_options['geometry'] = 'shell'
marginal_options['curve'] = True
marginal_options['minimum'] = True
marginal_options['minimum_int'] = True
marginal_options['plot_curve'] = True
marginal_options['solve'] = True
marginal_options['point_k'] = m
marginal_options['plot_point'] = True
marginal_options['viz_mode'] = 0
marginal_options['show_spectra'] = True
marginal_options['show_physical'] = True
marginal_options['save_physical'] = True
marginal_options['save_spectra'] = True
marginal_options['impose_symmetry'] = False
marginal_options['use_spherical_evp'] = False
marginal_options['curve_points'] = np.arange(max(1,m-2), m+2, 1)

# Compute
MarginalCurve.compute(gevp_opts, marginal_options)
