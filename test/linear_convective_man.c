#include "solver.h"
#include "error_norms.h"

#define FREE_FLUID_PERMEABILITY ((Real)1e30)
static Real amplitude(Real time)
{
    return (Real)0.5 + (Real)0.25 * sin((Real)40 * time);
}

static Real amplitude_derivative(Real time)
{
    return (Real)10 * cos((Real)40 * time);
}

static Real linear_velocity(Real x, Real y, Real z, Real time, int component)
{
    const Real a = amplitude(time);
    switch (component) {
    case 0:
        return a * sin(x) * sin(x) * sin((Real)2 * y) * sin(z) * sin(z);
    case 1:
        return -a * sin((Real)2 * x) * sin(y) * sin(y) * sin(z) * sin(z);
    default:
        return (Real)0;
    }
}

static Real zero_pressure(Real x, Real y, Real z, Real time)
{
    (void)x;
    (void)y;
    (void)z;
    (void)time;
    return (Real)0;
}

static Real free_fluid_permeability(Real x, Real y, Real z, Real time,
                                    int component)
{
    (void)x;
    (void)y;
    (void)z;
    (void)time;
    (void)component;
    return FREE_FLUID_PERMEABILITY;
}

static Real linear_forcing(Real x, Real y, Real z, Real time, int component)
{
    const Real a = amplitude(time);
    const Real ax = amplitude_derivative(time);
    const Real brinkman = (Real)NU / FREE_FLUID_PERMEABILITY;
    const Real sx = sin(x), sy = sin(y), sz = sin(z);
    const Real s2x = sin((Real)2*x), s2y = sin((Real)2*y);
    const Real c2x = cos((Real)2*x), c2y = cos((Real)2*y), c2z = cos((Real)2*z);
    const Real u0 = sx*sx*s2y*sz*sz;
    const Real v0 = -s2x*sy*sy*sz*sz;
    switch (component) {
    case 0: {
        const Real convection = a*a*(u0*s2x*s2y*sz*sz +
            v0*(Real)2*sx*sx*c2y*sz*sz);
        const Real laplacian = a*((Real)2*c2x*s2y*sz*sz -
            (Real)4*u0 + (Real)2*sx*sx*s2y*c2z);
        return ax*u0 + convection - (Real)NU*laplacian + brinkman*a*u0;
    }
    case 1: {
        const Real convection = a*a*(u0*(Real)-2*c2x*sy*sy*sz*sz +
            v0*(-s2x*s2y*sz*sz));
        const Real laplacian = a*((Real)4*s2x*sy*sy*sz*sz -
            (Real)2*s2x*c2y*sz*sz - (Real)2*s2x*sy*sy*c2z);
        return ax*v0 + convection - (Real)NU*laplacian + brinkman*a*v0;
    }
    default:
        return (Real)0;
    }
}

int main(void)
{
    Data data = {
        .name = "Linear convective manufactured solution",
        .bc_velocity = linear_velocity,
        .forcing_fn = linear_forcing,
        .porosity_fn = free_fluid_permeability,
        .porosity_time_dependent = 0,
        .velocity_fn = linear_velocity,
        .pressure_fn = zero_pressure,
    };
    SolverMemState state;
    SolverStats stats = {0};
    const Real velocity_time = (Real)STEPS * (Real)DT;
    const Real pressure_time = velocity_time - (Real)DT / (Real)2;

    solver_init(&state, &data, NULL);
    solver_solve(&state, &data, &stats, 0);

    const SolverErrorNorms errors =
        compute_solver_error_norms(&state, &data,
                                   velocity_time, pressure_time);
    print_solver_error_norms(&errors, velocity_time, pressure_time);
    return 0;
}
