"""Module provides the functions to generate the Boussinesq rotating thermal convection in a spherical shell (Toroidal/Poloidal formulation)"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp

import quicc.base.utils as utils
import quicc.geometry.spherical.shell as geo
import quicc.base.base_model as base_model
from quicc.geometry.spherical.shell_boundary import no_bc

class PhysicalModel(base_model.BaseModel):
    """Class to setup the Boussinesq rotating thermal convection in a spherical shell (Toroidal/Poloidal formulation)"""

    def periodicity(self):
        """Get the domain periodicity"""

        return [False, False, False]

    def nondimensional_parameters(self):
        """Get the list of nondimensional parameters"""

        return ["ekman", "prandtl", "rayleigh", "r_ratio", "heating","alpha","beta"]

        # alpha = 1 yields a linear gravity in r, whereas alpha !=1 leads to a alpha*r+(1-alpha)/r^2 type gravity profile
        # beta = 1 yields a linear dT/dr in r, whereas beta !=1 leads to a beta*r+(1-beta)/r^2 type dT/dr profile
        # heating=0: internal heating;   
        # heating=1: differential heating; 
        # heating=2: differential heating + internal heating (proportional to r)
        # heating=3: differential heating + internal heating (general form)

    def automatic_parameters(self, eq_params):
        """Extend parameters with automatically computable values"""

        # CFL parameters
        E = eq_params['ekman']
        d = {
                "cfl_inertial":0.1*E
                }

        rratio = eq_params['r_ratio']
        # Unit gap width
        if True:
            gap = {
                    "lower1d":rratio/(1.0 - rratio),
                    "upper1d":1.0/(1.0 - rratio)
                    }
        # Unit radius
        else:
            gap = {
                    "lower1d":rratio,
                    "upper1d":1.0
                    }

        d.update(gap)

        return d

    def config_fields(self):
        """Get the list of fields that need a configuration entry"""

        return ["velocity", "temperature"]

    def stability_fields(self):
        """Get the list of fields needed for linear stability calculations"""

        fields =  [("velocity","tor"), ("velocity","pol"), ("temperature","")]

        return fields

    def implicit_fields(self, field_row):
        """Get the list of coupled fields in solve"""

        fields =  [("velocity","tor"), ("velocity","pol"), ("temperature","")]

        return fields

    def explicit_fields(self, timing, field_row):
        """Get the list of fields with explicit linear dependence"""

        # Explicit linear terms
        if timing == self.EXPLICIT_LINEAR:
            if field_row == ("temperature",""):
                fields = [("velocity","pol")]
            else:
                fields = []

        # Explicit nonlinear terms
        elif timing == self.EXPLICIT_NONLINEAR:
            if field_row == ("temperature",""):
                fields = [("temperature","")]
            else:
                fields = []

        # Explicit update terms for next step
        elif timing == self.EXPLICIT_NEXTSTEP:
            fields = []

        return fields

    def block_size(self, res, eigs, bcs, field_row):
        """Create block size information"""

        tau_n = res[0]
        if self.use_galerkin:
            if field_row == ("velocity","tor") or field_row == ("temperature",""):
                shift_r = 2
            elif field_row == ("velocity","pol"):
                shift_r = 4
            else:
                shift_r = 0

            gal_n = (res[0] - shift_r)

        else:
            gal_n = tau_n
            shift_r = 0

        block_info = (tau_n, gal_n, (shift_r,0,0), 1)
        return block_info

    def stencil(self, res, eq_params, eigs, bcs, field_row, make_square):
        """Create the galerkin stencil"""

        assert(eigs[0].is_integer())

        m = int(eigs[0])

        ri, ro = (self.automatic_parameters(eq_params)['lower1d'], self.automatic_parameters(eq_params)['upper1d'])

        # Get boundary condition
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_row)
        return geo.stencil(res[0], ri, ro, res[1], m, bc, make_square)

    def stability_sizes(self, res, eigs, bcs):
        """Get the block sizes in the stability calculation matrix"""

        assert(eigs[0].is_integer())

        m = int(eigs[0])

        # Block sizes
        blocks = []
        for f in self.stability_fields():
            blocks.append(self.block_size(res, eigs, bcs, f)[1]*(res[1]-m))

        # Invariant size (local dimension in spectral space, no restriction)
        invariant = (res[0],)*len(self.stability_fields())

        # Index shift
        shift = m

        return (blocks, invariant, shift)

    def equation_info(self, res, field_row):
        """Provide description of the system of equation"""

        # Matrix operator is complex except for vorticity and mean temperature
        is_complex = True

        # Index mode: SLOWEST_SINGLE_RHS, SLOWEST_MULTI_RHS, MODE, SINGLE
        index_mode = self.SLOWEST_SINGLE_RHS

        return self.compile_equation_info(res, field_row, is_complex, index_mode)

    def convert_bc(self, eq_params, eigs, bcs, field_row, field_col):
        """Convert simulation input boundary conditions to ID"""

        # Solver: no tau boundary conditions
        if bcs["bcType"] == self.SOLVER_NO_TAU and not self.use_galerkin:
            bc = no_bc()

        # Solver: tau and Galerkin
        elif bcs["bcType"] == self.SOLVER_HAS_BC or bcs["bcType"] == self.SOLVER_NO_TAU:
            bc = no_bc()
            bcId = bcs.get(field_col[0], -1)
            if bcId == 0:
                if self.use_galerkin:
                    if field_col == ("velocity","tor"):
                        bc = {0:-20, 'rt':0}
                    elif field_col == ("velocity","pol"):
                        bc = {0:-40, 'rt':0}
                    elif field_col == ("temperature",""):
                        bc = {0:-20, 'rt':0}

                else:
                    if field_row == ("velocity","tor") and field_col == ("velocity","tor"):
                            bc = {0:20}
                    elif field_row == ("velocity","pol") and field_col == ("velocity","pol"):
                            bc = {0:40}
                    elif field_row == ("temperature","") and field_col == ("temperature",""):
                            bc = {0:20}

            elif bcId == 1:
                if self.use_galerkin:
                    if field_col == ("velocity","tor"):
                        bc = {0:-22, 'rt':0}
                    elif field_col == ("velocity","pol"):
                        bc = {0:-41, 'rt':0}

                else:
                    if field_row == ("velocity","tor") and field_col == ("velocity","tor"):
                        bc = {0:22}
                    elif field_row == ("velocity","pol") and field_col == ("velocity","pol"):
                        bc = {0:41}

            # Set LHS galerkin restriction
            if self.use_galerkin:
                if field_row == ("velocity","tor"):
                    bc['rt'] = 2
                elif field_row == ("velocity","pol"):
                    bc['rt'] = 4
                elif field_row == ("temperature",""):
                    bc['rt'] = 2

        # Stencil:
        elif bcs["bcType"] == self.STENCIL:
            if self.use_galerkin:
                bcId = bcs.get(field_col[0], -1)
                if bcId == 0:
                    if field_col == ("velocity","tor"):
                        bc = {0:-20, 'rt':2}
                    elif field_col == ("velocity","pol"):
                        bc = {0:-40, 'rt':4}
                    elif field_col == ("temperature",""):
                        bc = {0:-20, 'rt':2}

                elif bcId == 1:
                    if field_col == ("velocity","tor"):
                        bc = {0:-22, 'rt':2}
                    elif field_col == ("velocity","pol"):
                        bc = {0:-41, 'rt':4}

        # Field values to RHS:
        elif bcs["bcType"] == self.FIELD_TO_RHS:
            bc = no_bc()
            if self.use_galerkin:
                if field_row == ("velocity","tor"):
                    bc['rt'] = 2
                elif field_row == ("velocity","pol"):
                    bc['rt'] = 4
                elif field_row == ("temperature",""):
                    bc['rt'] = 2

        else:
            bc = no_bc()

        return bc

    def explicit_block(self, res, eq_params, eigs, bcs, field_row, field_col, restriction = None):
        """Create matrix block for explicit linear term"""

        assert(eigs[0].is_integer())

        m = int(eigs[0])

        Ra_eff, bg_eff = self.nondimensional_factors(eq_params)

        ri, ro = (self.automatic_parameters(eq_params)['lower1d'], self.automatic_parameters(eq_params)['upper1d'])

        mat = None
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_col)
        if field_row == ("temperature","") and field_col == ("velocity","pol"):
            if eq_params["heating"] == 0:
                mat = geo.i2r2(res[0], ri, ro, res[1], m, bc, -bg_eff, with_sh_coeff = 'laplh', restriction = restriction)
            elif eq_params["heating"] == 1:
                mat = geo.i2(res[0], ri, ro, res[1], m, bc, -bg_eff, with_sh_coeff = 'laplh', restriction = restriction)
            elif eq_params["heating"] == 2:
                beta = eq_params['beta']
                # assumes length-scale is the depth of the shell
                c1 = beta*bg_eff/ro # internal heating contribution
                c2 = ro**2*(1-beta*bg_eff) # differential heating contribution
                if beta==1: # same as for heating=0
                    mat = geo.i2r2(res[0], ri, ro, res[1], m, bc, c1, with_sh_coeff = 'laplh', restriction = restriction)
                else:
                    mat = geo.i2r3(res[0], ri, ro, res[1], m, bc, c1, with_sh_coeff = 'laplh', restriction = restriction) + geo.i2(res[0], ri, ro, res[1], m, bc, c2, with_sh_coeff = 'laplh', restriction = restriction)

        if mat is None:
            raise RuntimeError("Equations are not setup properly!")

        return mat

    def nonlinear_block(self, res, eq_params, eigs, bcs, field_row, field_col, restriction = None):
        """Create the explicit nonlinear operator"""

        assert(eigs[0].is_integer())

        m = int(eigs[0])

        ri, ro = (self.automatic_parameters(eq_params)['lower1d'], self.automatic_parameters(eq_params)['upper1d'])

        mat = None
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_col)
        if field_row == ("temperature","") and field_col == field_row:
            if eq_params["heating"] == 0:
                mat = geo.i2r2(res[0], ri, ro, res[1], m, bc, restriction = restriction)
            elif eq_params["heating"] == 1:
                mat = geo.i2r3(res[0], ri, ro, res[1], m, bc, restriction = restriction)
            elif eq_params["heating"] == 2:
                beta = eq_params['beta']
                if beta==1: # same as for heating=0
                    mat = geo.i2r2(res[0], ri, ro, res[1], m, bc, restriction = restriction)
                else:
                    mat = geo.i2r3(res[0], ri, ro, res[1], m, bc, restriction = restriction)

        if mat is None:
            raise RuntimeError("Equations are not setup properly!")

        return mat

    def implicit_block(self, res, eq_params, eigs, bcs, field_row, field_col, restriction = None):
        """Create matrix block linear operator"""

        assert(eigs[0].is_integer())

        Pr = eq_params['prandtl']
        T = 1.0/eq_params['ekman']
        Ra_eff, bg_eff = self.nondimensional_factors(eq_params)

        # Rescale time for eigensolver
        if self.linearize:
            c_dt = self.rescale_time(eq_params)
        else:
            c_dt = 1.0

        m = int(eigs[0])

        ri, ro = (self.automatic_parameters(eq_params)['lower1d'], self.automatic_parameters(eq_params)['upper1d'])

        mat = None
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_col)
        if field_row == ("velocity","tor"):
            if field_col == ("velocity","tor"):
                mat = geo.i2r2lapl(res[0], ri, ro, res[1], m, bc, 1.0/c_dt, l_zero_fix = 'zero', restriction = restriction)
                bc[0] = min(bc[0], 0)
                mat = mat + geo.i2r2(res[0], ri, ro, res[1], m, bc, 1j*m*T/c_dt, with_sh_coeff = 'invlaplh', l_zero_fix = 'zero', restriction = restriction)

            elif field_col == ("velocity","pol"):
                mat = geo.i2r2coriolis(res[0], ri, ro, res[1], m, bc, -T/c_dt, with_sh_coeff = 'invlaplh', l_zero_fix = 'zero', restriction = restriction)

            elif field_col == ("temperature",""):
                mat = geo.zblk(res[0], ri, ro, res[1], m, bc)

        elif field_row == ("velocity","pol"):
            if field_col == ("velocity","tor"):
                mat = geo.i4r4coriolis(res[0], ri, ro, res[1], m, bc, T/c_dt, with_sh_coeff = 'invlaplh', l_zero_fix = 'zero', restriction = restriction)

            elif field_col == ("velocity","pol"):
                mat = geo.i4r4lapl2(res[0], ri, ro, res[1], m, bc, 1.0/c_dt, l_zero_fix = 'zero', restriction = restriction)
                bc[0] = min(bc[0], 0)
                mat = mat + geo.i4r4lapl(res[0], ri, ro, res[1], m, bc, 1j*m*T/c_dt, with_sh_coeff = 'invlaplh', l_zero_fix = 'zero', restriction = restriction)

            elif field_col == ("temperature",""):
                mat = geo.i4r4(res[0], ri, ro, res[1], m, bc, -Ra_eff/c_dt**2, l_zero_fix = 'zero', restriction = restriction)

        elif field_row == ("temperature",""):
            if field_col == ("velocity","tor"):
                mat = geo.zblk(res[0], ri, ro, res[1], m, bc)

            elif field_col == ("velocity","pol"):
                if self.linearize:
                    if eq_params["heating"] == 0:
                        mat = geo.i2r2(res[0], ri, ro, res[1], m, bc, bg_eff/Pr, with_sh_coeff = 'laplh', restriction = restriction)
                    elif eq_params["heating"] == 1:
                        mat = geo.i2(res[0], ri, ro, res[1], m, bc, bg_eff/Pr, with_sh_coeff = 'laplh', restriction = restriction)
                    elif eq_params["heating"] == 2:
                        beta = eq_params['beta']
                        # assumes length-scale is the depth of the shell
                        c1 = beta*bg_eff/ro # internal heating contribution
                        c2 = ro**2*(1-beta*bg_eff) # differential heating contribution
                        if beta==1: # same as for heating=0
                            mat = geo.i2r2(res[0], ri, ro, res[1], m, bc, c1/Pr, with_sh_coeff = 'laplh', restriction = restriction)
                        else:
                            mat = geo.i2r3(res[0], ri, ro, res[1], m, bc, c1/Pr, with_sh_coeff = 'laplh', restriction = restriction) + geo.i2(res[0], ri, ro, res[1], m, bc, c2/Pr, with_sh_coeff = 'laplh', restriction = restriction)

                else:
                    mat = geo.zblk(res[0], ri, ro, res[1], m, bc)

            elif field_col == ("temperature",""):
                if eq_params["heating"] == 0:
                    mat = geo.i2r2lapl(res[0], ri, ro, res[1], m, bc, 1.0/(Pr*c_dt), restriction = restriction)
                elif eq_params["heating"] == 1:
                    mat = geo.i2r3lapl(res[0], ri, ro, res[1], m, bc, 1.0/(Pr*c_dt), restriction = restriction)
                elif eq_params["heating"] == 2:
                    beta = eq_params['beta']
                    if beta==1: # same as for heating=0
                        mat = geo.i2r2lapl(res[0], ri, ro, res[1], m, bc, 1.0/(Pr*c_dt), restriction = restriction)
                    else:
                        mat = geo.i2r3lapl(res[0], ri, ro, res[1], m, bc, 1.0/(Pr*c_dt), restriction = restriction)

        if mat is None:
            raise RuntimeError("Equations are not setup properly!")

        return mat

    def time_block(self, res, eq_params, eigs, bcs, field_row, restriction = None):
        """Create matrix block of time operator"""

        assert(eigs[0].is_integer())

        m = int(eigs[0])

        ri, ro = (self.automatic_parameters(eq_params)['lower1d'], self.automatic_parameters(eq_params)['upper1d'])

        mat = None
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_row)
        if field_row == ("velocity","tor"):
            mat = geo.i2r2(res[0], ri, ro, res[1], m, bc, l_zero_fix = 'set', restriction = restriction)

        elif field_row == ("velocity","pol"):
            mat = geo.i4r4lapl(res[0], ri, ro, res[1], m, bc, l_zero_fix = 'set', restriction = restriction)

        elif field_row == ("temperature",""):
            if eq_params["heating"] == 0:
                mat = geo.i2r2(res[0], ri, ro, res[1], m, bc, restriction = restriction)
            elif eq_params["heating"] == 1:
                mat = geo.i2r3(res[0], ri, ro, res[1], m, bc, restriction = restriction)
            elif eq_params["heating"] == 2:
                beta = eq_params['beta']
                if beta == 1: # same as for heating=0
                    mat = geo.i2r2(res[0], ri, ro, res[1], m, bc, restriction = restriction)
                else:
                    mat = geo.i2r3(res[0], ri, ro, res[1], m, bc, restriction = restriction)

        if mat is None:
            raise RuntimeError("Equations are not setup properly!")

        return mat

    def boundary_block(self, res, eq_params, eigs, bcs, field_row, field_col, restriction = None):
        """Create matrix block of boundary operator"""

        assert(eigs[0].is_integer())

        m = int(eigs[0])

        ri, ro = (self.automatic_parameters(eq_params)['lower1d'], self.automatic_parameters(eq_params)['upper1d'])

        mat = None
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_col)
        if field_row in [("velocity","tor"), ("velocity","pol")] and field_row == field_col:
            mat = geo.zblk(res[0], ri, ro, res[1], m, bc, l_zero_fix = 'zero', restriction = restriction)
        else:
            mat = geo.zblk(res[0], ri, ro, res[1], m, bc, restriction = restriction)

        if mat is None:
            raise RuntimeError("Equations are not setup properly!")

        return mat

    def rescale_time(self, eq_params):
        """Rescale time for linear stability calculation"""

        T = 1.0/(eq_params['ekman']*(1.0-eq_params['r_ratio'])**2)
        #c_dt = T**(2./3.)*(0.4715 - 0.6089*T**(-1/3.))
        c_dt = T**(2./3.)*(0.3144 - 0.6089*T**(-1./3.)) + T**(1./3.)*0.5186*int(0.3029*T**(1./3.))
        c_dt = c_dt*(1.0-eq_params['r_ratio'])**2
        c_dt = 1.0
        return c_dt

    def nondimensional_factors(self, eq_params):
        """Compute the effective Rayleigh number and background depending on nondimensionalisation"""

        Ra = eq_params['rayleigh']
        ri, ro = (self.automatic_parameters(eq_params)['lower1d'], self.automatic_parameters(eq_params)['upper1d'])
        rratio = eq_params['r_ratio']
        T = 1.0/eq_params['ekman']

        # Easy switch from nondimensionalistion by R_o (Dormy) and (R_o - R_i) (Christensen)
        # Parameters match as:  Dormy   Christensen
        #                       Ra      Ra/(R_o*R_i*Ta^0.5)
        #                       Ta      Ta*(1-R_i/R_o)^4
        if ro == 1.0:
            # R_o rescaling
            Ra_eff = Ra
            bg_eff = 1.0
        elif eq_params['heating'] == 0:
            # (R_o - R_i) rescaling
            Ra_eff = (Ra*T/ro)
            bg_eff = 2.0/(ro*(1.0 + rratio))
        elif eq_params['heating'] == 1:
            # (R_o - R_i) rescaling
            Ra_eff = (Ra*T/ro)
            bg_eff = ro**2*rratio       
        # Rescale Rayleigh by E^{-4/3}
        # Ra_eff = Ra_eff*T**(1./3.)
        elif eq_params['heating'] == 2:
            Ra_eff = (Ra*T/ro)
            bg_eff = 1


        return (Ra_eff, bg_eff)
