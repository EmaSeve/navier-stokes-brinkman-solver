#include "solver.h"
#include "error_norms.h"

#ifdef USE_FLOAT
#define REAL_SIN sinf
#define REAL_COS cosf
#else
#define REAL_SIN sin
#define REAL_COS cos
#endif

#define FREE_FLUID_PERMEABILITY ((Real)1e30)

/* Manufactured solution used for the convective convergence study. */
static Real manufactured_velocity(Real x, Real y, Real z, Real t,
                                  int component)
{
    switch (component) {
    case 0:
        return REAL_SIN(x) * REAL_COS(t + y) * REAL_SIN(z);
    case 1:
        return REAL_COS(x) * REAL_SIN(t + y) * REAL_SIN(z);
    case 2:
        return (Real)2 * REAL_COS(x) * REAL_COS(t + y) * REAL_COS(z);
    default:
        return (Real)0;
    }
}

static Real manufactured_pressure(Real x, Real y, Real z, Real t)
{
    return (Real)3 * (Real)NU * REAL_COS(x) * REAL_COS(y) *
           REAL_COS(z) * REAL_COS(t);
}

static Real free_fluid_permeability(Real x, Real y, Real z, Real t,
                                    int component)
{
    (void)x;
    (void)y;
    (void)z;
    (void)t;
    (void)component;
    return FREE_FLUID_PERMEABILITY;
}

/*
 * f = du/dt + (u . grad)u - NU lap(u) + grad(p) + NU/K u.
 */
static Real manufactured_forcing(Real x, Real y, Real z, Real t,
                                 int component)
{
    const Real sin_x = REAL_SIN(x);
    const Real cos_x = REAL_COS(x);
    const Real sin_ty = REAL_SIN(t + y);
    const Real cos_ty = REAL_COS(t + y);
    const Real sin_z = REAL_SIN(z);
    const Real cos_z = REAL_COS(z);
    const Real sin_y = REAL_SIN(y);
    const Real cos_y = REAL_COS(y);
    const Real cos_t = REAL_COS(t);
    const Real u = sin_x * cos_ty * sin_z;
    const Real v = cos_x * sin_ty * sin_z;
    const Real w = (Real)2 * cos_x * cos_ty * cos_z;
    const Real conv_x =
        u * cos_x * cos_ty * sin_z -
        v * sin_x * sin_ty * sin_z +
        w * sin_x * cos_ty * cos_z;
    const Real conv_y =
        -u * sin_x * sin_ty * sin_z +
        v * cos_x * cos_ty * sin_z +
        w * cos_x * sin_ty * cos_z;
    const Real conv_z =
        -u * (Real)2 * sin_x * cos_ty * cos_z -
        v * (Real)2 * cos_x * sin_ty * cos_z -
        w * (Real)2 * cos_x * cos_ty * sin_z;
    const Real brinkman = (Real)NU / FREE_FLUID_PERMEABILITY;

    switch (component) {
    case 0:
        return -sin_x * sin_ty * sin_z + conv_x +
               (Real)3 * (Real)NU * u +
               -(Real)3 * (Real)NU * sin_x * cos_y * cos_z * cos_t +
               brinkman * u;
    case 1:
        return cos_x * cos_ty * sin_z + conv_y +
               (Real)3 * (Real)NU * v +
               -(Real)3 * (Real)NU * cos_x * sin_y * cos_z * cos_t +
               brinkman * v;
    case 2:
        return -(Real)2 * cos_x * sin_ty * cos_z + conv_z +
               (Real)3 * (Real)NU * w +
               -(Real)3 * (Real)NU * cos_x * cos_y * sin_z * cos_t +
               brinkman * w;
    default:
        return (Real)0;
    }
}

int main(void)
{
    Data data = {
        .name = "Convective manufactured solution",
        .bc_velocity = manufactured_velocity,
        .forcing_fn = manufactured_forcing,
        .porosity_fn = free_fluid_permeability,
        .porosity_time_dependent = 0,
        .velocity_fn = manufactured_velocity,
        .pressure_fn = manufactured_pressure,
    };
    SolverMemState state;
    SolverStats stats = {0};
    const Real velocity_time = (Real)STEPS * (Real)DT;
    const Real pressure_time = velocity_time - (Real)DT / (Real)2;

    solver_init(&state, &data, NULL);
    solver_solve(&state, &data, &stats, 0);

    const SolverErrorNorms errors = compute_solver_error_norms(
        &state, &data, velocity_time, pressure_time);
    print_solver_error_norms(&errors, velocity_time, pressure_time);
    return 0;
}
