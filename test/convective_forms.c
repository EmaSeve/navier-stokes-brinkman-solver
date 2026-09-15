#include "field.h"
#include "momentum.h"
#include "physics.h"
#include "solver.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef USE_FLOAT
#define TEST_EPSILON ((Real)1e-5)
#else
#define TEST_EPSILON ((Real)1e-12)
#endif

static size_t index_of(int i, int j, int k)
{
    return (size_t)k * (size_t)WIDTH * (size_t)HEIGHT +
           (size_t)j * (size_t)WIDTH + (size_t)i;
}

static void assert_close(const char *name, Real actual, Real expected)
{
    if (fabs((double)(actual - expected)) > (double)TEST_EPSILON) {
        fprintf(stderr, "%s: expected %.17g, got %.17g\n",
                name, (double)expected, (double)actual);
        exit(EXIT_FAILURE);
    }
}

static void assert_true(const char *name, int condition)
{
    if (!condition) {
        fprintf(stderr, "%s: assertion failed\n", name);
        exit(EXIT_FAILURE);
    }
}

static void alloc_velocity_history(SolverMemState *state)
{
    vectorField_alloc(&state->u);
    vectorField_alloc(&state->u_prev);
    vectorField_alloc(&state->u_star);
}

static void free_velocity_history(SolverMemState *state)
{
    free(state->u.v_x);
    free(state->u.v_y);
    free(state->u.v_z);
    free(state->u_prev.v_x);
    free(state->u_prev.v_y);
    free(state->u_prev.v_z);
    free(state->u_star.v_x);
    free(state->u_star.v_y);
    free(state->u_star.v_z);
}

static void test_valid_convective_point(void)
{
    assert_true("x component interior", valid_convective_point(1, 1, 1, 0));
    assert_true("y component interior", valid_convective_point(1, 1, 1, 1));
    assert_true("z component interior", valid_convective_point(1, 1, 1, 2));
    assert_true("lower x face rejected", !valid_convective_point(0, 1, 1, 0));
    assert_true("lower y face rejected", !valid_convective_point(1, 0, 1, 1));
    assert_true("lower z face rejected", !valid_convective_point(1, 1, 0, 2));
    assert_true("upper x face rejected",
                !valid_convective_point(WIDTH - 1, 1, 1, 0));
    assert_true("upper y face rejected",
                !valid_convective_point(1, HEIGHT - 1, 1, 1));
    assert_true("upper z face rejected",
                !valid_convective_point(1, 1, DEPTH - 1, 2));
    assert_true("invalid component rejected", !valid_convective_point(1, 1, 1, 3));
}

static void test_compute_extrapolated_velocity(void)
{
    SolverMemState state = {0};
    const size_t interior = index_of(2, 2, 2);
    const size_t boundary = index_of(0, 2, 2);
    alloc_velocity_history(&state);

    for (size_t p = 0; p < GRID_CELLS; ++p) {
        state.u.v_x[p] = (Real)10 + (Real)p;
        state.u.v_y[p] = (Real)20 + (Real)p;
        state.u.v_z[p] = (Real)30 + (Real)p;
        state.u_prev.v_x[p] = (Real)1 + (Real)p;
        state.u_prev.v_y[p] = (Real)2 + (Real)p;
        state.u_prev.v_z[p] = (Real)3 + (Real)p;
    }

    compute_extrapolated_velocity(&state, 0);
    assert_close("step zero copies u x", state.u_star.v_x[interior], state.u.v_x[interior]);
    assert_close("step zero copies u y", state.u_star.v_y[interior], state.u.v_y[interior]);
    assert_close("step zero copies u z", state.u_star.v_z[interior], state.u.v_z[interior]);

    compute_extrapolated_velocity(&state, 1);
    assert_close("interior extrapolated x", state.u_star.v_x[interior],
                 (Real)1.5 * state.u.v_x[interior] -
                 (Real)0.5 * state.u_prev.v_x[interior]);
    assert_close("interior extrapolated y", state.u_star.v_y[interior],
                 (Real)1.5 * state.u.v_y[interior] -
                 (Real)0.5 * state.u_prev.v_y[interior]);
    assert_close("interior extrapolated z", state.u_star.v_z[interior],
                 (Real)1.5 * state.u.v_z[interior] -
                 (Real)0.5 * state.u_prev.v_z[interior]);
    assert_close("boundary extrapolated x", state.u_star.v_x[boundary],
                 (Real)1.5 * state.u.v_x[boundary] -
                 (Real)0.5 * state.u_prev.v_x[boundary]);
    assert_close("boundary extrapolated y", state.u_star.v_y[boundary],
                 (Real)1.5 * state.u.v_y[boundary] -
                 (Real)0.5 * state.u_prev.v_y[boundary]);
    assert_close("boundary extrapolated z", state.u_star.v_z[boundary],
                 (Real)1.5 * state.u.v_z[boundary] -
                 (Real)0.5 * state.u_prev.v_z[boundary]);

    free_velocity_history(&state);
}

static void test_swap_velocity_buffers(void)
{
    SolverMemState state = {0};
    Real *u_x;
    Real *u_y;
    Real *u_z;
    Real *u_prev_x;
    Real *u_prev_y;
    Real *u_prev_z;
    alloc_velocity_history(&state);

    u_x = state.u.v_x;
    u_y = state.u.v_y;
    u_z = state.u.v_z;
    u_prev_x = state.u_prev.v_x;
    u_prev_y = state.u_prev.v_y;
    u_prev_z = state.u_prev.v_z;
    state.u.v_x[0] = (Real)1;
    state.u_prev.v_x[0] = (Real)2;

    swap_velocity_buffers(&state);
    assert_true("x buffers swapped", state.u.v_x == u_prev_x && state.u_prev.v_x == u_x);
    assert_true("y buffers swapped", state.u.v_y == u_prev_y && state.u_prev.v_y == u_y);
    assert_true("z buffers swapped", state.u.v_z == u_prev_z && state.u_prev.v_z == u_z);
    assert_close("x values not copied", state.u.v_x[0], (Real)2);
    assert_close("previous x values not copied", state.u_prev.v_x[0], (Real)1);

    free_velocity_history(&state);
}

static void fill_constant(VectorField *vel)
{
    for (size_t p = 0; p < GRID_CELLS; ++p) {
        vel->v_x[p] = (Real)2;
        vel->v_y[p] = (Real)-3;
        vel->v_z[p] = (Real)4;
    }
}

static void fill_divergent_linear(VectorField *vel)
{
    for (int k = 0; k < DEPTH; ++k) {
        for (int j = 0; j < HEIGHT; ++j) {
            for (int i = 0; i < WIDTH; ++i) {
                size_t p = index_of(i, j, k);
                vel->v_x[p] = staggered_physical_coord(i, 0);
                vel->v_y[p] = (Real)0;
                vel->v_z[p] = (Real)0;
            }
        }
    }
}

static void fill_incompressible_linear(VectorField *vel)
{
    for (int k = 0; k < DEPTH; ++k) {
        for (int j = 0; j < HEIGHT; ++j) {
            for (int i = 0; i < WIDTH; ++i) {
                size_t p = index_of(i, j, k);
                vel->v_x[p] = centered_physical_coord(j, 1); /* u = y */
                vel->v_y[p] = centered_physical_coord(i, 0); /* v = x */
                vel->v_z[p] = (Real)0;
            }
        }
    }
}

static Real temporal_velocity(Real x, Real y, Real z, Real t, int component)
{
    switch (component) {
    case 0:
        return (Real)sin((double)x) * (Real)cos((double)(t + y)) *
               (Real)sin((double)z);
    case 1:
        return (Real)cos((double)x) * (Real)sin((double)(t + y)) *
               (Real)sin((double)z);
    case 2:
        return (Real)2 * (Real)cos((double)x) *
               (Real)cos((double)(t + y)) * (Real)cos((double)z);
    default:
        return (Real)0;
    }
}

static void fill_velocity_at_time(VectorField *vel, Real time)
{
    for (int k = 0; k < DEPTH; ++k) {
        for (int j = 0; j < HEIGHT; ++j) {
            for (int i = 0; i < WIDTH; ++i) {
                const size_t p = index_of(i, j, k);
                const Real x = centered_physical_coord(i, 0);
                const Real y = centered_physical_coord(j, 1);
                const Real z = centered_physical_coord(k, 2);

                vel->v_x[p] = temporal_velocity(
                    staggered_physical_coord(i, 0), y, z, time, 0);
                vel->v_y[p] = temporal_velocity(
                    x, staggered_physical_coord(j, 1), z, time, 1);
                vel->v_z[p] = temporal_velocity(
                    x, y, staggered_physical_coord(k, 2), time, 2);
            }
        }
    }
}

static Real convective_midpoint_error(const VectorField *extrapolated,
                                      const VectorField *exact_midpoint)
{
    Real sum = (Real)0;
    size_t count = 0;

    for (int k = 1; k < DEPTH - 1; ++k) {
        for (int j = 1; j < HEIGHT - 1; ++j) {
            for (int i = 1; i < WIDTH - 1; ++i) {
                for (int component = 0; component < 3; ++component) {
                    const Real difference =
                        convective_term(i, j, k, component, extrapolated,
                                        NULL, (Real)0) -
                        convective_term(i, j, k, component, exact_midpoint,
                                        NULL, (Real)0);
                    sum += difference * difference;
                    ++count;
                }
            }
        }
    }

    return (Real)sqrt((double)(sum / (Real)count));
}

static Real temporal_extrapolation_error(Real delta)
{
    SolverMemState state = {0};
    VectorField midpoint;
    const Real time = (Real)0.4;
    Real error;

    alloc_velocity_history(&state);
    vectorField_alloc(&midpoint);
    fill_velocity_at_time(&state.u, time);
    fill_velocity_at_time(&state.u_prev, time - delta);
    fill_velocity_at_time(&midpoint, time + (Real)0.5 * delta);
    compute_extrapolated_velocity(&state, 1);
    error = convective_midpoint_error(&state.u_star, &midpoint);

    free(midpoint.v_x);
    free(midpoint.v_y);
    free(midpoint.v_z);
    free_velocity_history(&state);
    return error;
}

static void test_temporal_convective_extrapolation(void)
{
    const Real coarse_error = temporal_extrapolation_error((Real)0.10);
    const Real fine_error = temporal_extrapolation_error((Real)0.05);

    assert_true("temporal extrapolation has nonzero error", coarse_error > (Real)0);
    assert_true("temporal extrapolation is second order",
                fine_error / coarse_error < (Real)0.27);
}

int main(void)
{
    const int i = 3, j = 3, k = 3;
    VectorField vel;
    vectorField_alloc(&vel);

    test_valid_convective_point();
    test_compute_extrapolated_velocity();
    test_swap_velocity_buffers();
    test_temporal_convective_extrapolation();

    fill_constant(&vel);
    for (int component = 0; component < 3; ++component) {
        assert_close("constant advective", advective_form(i, j, k, component, &vel,
                                                            NULL, (Real)0), 0);
        assert_close("constant divergence", divergence_form(i, j, k, component, &vel,
                                                              NULL, (Real)0), 0);
    }

    fill_divergent_linear(&vel);
    Real x_face = staggered_physical_coord(i, 0);
    assert_close("linear advective", advective_form(i, j, k, 0, &vel,
                                                      NULL, (Real)0),
                 (Real)0.5 * x_face);
    assert_close("linear divergence", divergence_form(i, j, k, 0, &vel,
                                                        NULL, (Real)0),
                 x_face);

    fill_incompressible_linear(&vel);
    Real y_face = staggered_physical_coord(j, 1);
    assert_close("incompressible x advective", advective_form(i, j, k, 0, &vel,
                                                                NULL, (Real)0),
                 (Real)0.5 * x_face);
    assert_close("incompressible x divergence", divergence_form(i, j, k, 0, &vel,
                                                                  NULL, (Real)0),
                 (Real)0.5 * x_face);
    assert_close("incompressible y advective", advective_form(i, j, k, 1, &vel,
                                                                NULL, (Real)0),
                 (Real)0.5 * y_face);
    assert_close("incompressible y divergence", divergence_form(i, j, k, 1, &vel,
                                                                  NULL, (Real)0),
                 (Real)0.5 * y_face);

    free(vel.v_x);
    free(vel.v_y);
    free(vel.v_z);
    puts("convective_forms: all tests passed");
    return EXIT_SUCCESS;
}
