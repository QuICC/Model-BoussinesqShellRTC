# Boussinesq rotating thermal convection(RTC) in a spherical shell


The model equations are

$$
(\partial_t -\Delta){\bf u}= {\bf u}\times(\nabla\times{\bf u})- \frac{1}{E}\hat{\bf z}\times{\bf u}  -\nabla \Pi + \frac{Ra}{E} \left(\frac{d}{r_o}\right)\Theta{\bf r}\\
$$

$$
(\partial_t -\frac{1}{Pr}\Delta)\Theta = -{\bf u}\cdot\nabla\Theta + \mathcal{H}{\bf u}\cdot{\bf r}
$$

with the parameters defined as

$$
Pr = \frac{\nu}{\kappa};\quad E=\frac{\nu}{2\Omega d^2};\quad Ra=\frac{ \alpha g_o \mathcal{T} d}{2\nu\Omega}
$$

where $\nu$, $\kappa$, $\alpha$ are the kinematic viscosity, thermal conductivity, and thermal expansion coefficents; $\Omega$ is the rotation rate; $r_o$ is the outer radius of the spherical shell and $d$ is its thickness; $g_o$ is defined via the gravitational acceleration, ${\bf g} = -g_o {\bf r}/r_o$. $\mathcal{T}$ is a characteristic temperature and $\mathcal{H}$ is the effective temperature background. The spherical shell has an inner radius of $r_i$, outeer radius $r_o$ and thickness $d=r_o-r_i$.

The model equations are non-dimensionalised using $d$ as the characteristic length-scale and the viscous time-scale $\nu/d^2$. Note that in the model equation, the factor $d/r_o$ as written in the buoyancy term, involves dimensional quantities. If the outer radius $r_o$ is already non-dimensional, that factor is simply $1/r_o$.

The definition of both $\mathcal{T}$ and $\mathcal{H}$ depend on the mode of thermal driving (differential or internal heating, controlled by the parameter `heating`) and on the boundary conditions (fixed-temperature or fixed-flux):

- For fixed temperature boundary conditions, the background temperature is as given in Eq. 1.6 of Dormy et al., 2004 (DOI: 10.1017/S0022112003007316). Then $\mathcal{T} = \Delta T = T_i - T_o$ (where $T$ is the temperature of the resulting reference state) and:
    + $\mathcal{H}=$`bg_eff` for uniform internal heating (`heating=0`)
    + $\mathcal{H}=$ `bg_eff` $/r^3$ for differential heating (`heating=1`)  
where `bg_eff` is a pre-factor, function of $r_o$ and $d$.

- **The fixed-flux case is not currently implemented cleanly**. For fixed-flux at either one of the boundaries, the relevant thermal scale is defined by a thermal gradient $\beta$, so that the characteristic temperature is $\mathcal{T}=\beta d$. The exact definition of $\beta$ (and therefore of the Rayleigh number) depends on the mode of heating, but still:
    + $\mathcal{H}=$`bg_eff` for uniform internal heating (`heating=0`)
    + $\mathcal{H}=$ `bg_eff` $/r^3$ for differential heating (`heating=1`)  
**Currently the correct** `bg_eff` **is not implemented for fixed-flux conditions, but the code can still be used to calculate the onset of convection. The resulting Rayleigh number needs, however, to be corrected to account for this discrepancy.**


**At present, only uniform heating is implemented in the internally heated case.**

The system is solved using a Toroidal/Poloidal decomposition of the velocity

$$
{\bf u} = \nabla\times T {\bf r} + \nabla\times\nabla\times P {\bf r} .
$$


## General background profiles

**The following functionalities are currently implemented only for the Explicit models, and not for the SplitEqaution formalism.**

The parameters $\alpha$ (`alpha`) and $\beta$ (`beta`) can be used to set up general gravity and background temperature profiles. The following are currently implemented:

### Gravity:
The following general gravity profile is implemented:
$$
{\bf g} = -g_o\left[\left(\frac{r}{r_o}\right)\alpha + \left(\frac{r_o}{r}\right)^2(1-\alpha)\right]\hat{\bf r}
$$
where, with $\alpha=1$ a linear profile is obtained, adequate for an inner core of the same density as the outer core.

**This is the feature not implemented in the SplitEqaution formalism.**

### Temperature:

General background temperature profiles, $T$, containing a mix of both internal and differential heating are implemented.

There are two options, chosed via the parameter `heating`:

- `heating`=2 implements a purely linear internal heating component, resulting in:
$$
T = T'_o\left[\left(\frac{r}{r_o}\right)\beta + \left(\frac{r_o}{r}\right)^2(1-\beta)\right]
$$

- `heating`=3 implements a fully general form for the  internal heating component, resulting in:
$$
\frac{dT}{dr} = T'_o\left[\frac{1}{1-\chi^3}\left(\frac{r}{r_o}\right)\beta + \left(\frac{r_o}{r}\right)^2\left(1-\frac{\beta}{1-\chi^3}\right)\right]
$$

In both cases, tempreature is nondimensionalised via $T'_o$, the temperature gradient at the outer boundary, $r_o$. In the above definition of the Rayleigh number, therefore: $\mathcal{T} = D T'_o$.