#include "physics.h"
#include "solver.h"

#include <stdbool.h>

static Real spacing_from_component(int component) {
    switch (component) {
        case 0:
            return (Real)DX;
        case 1:
            return (Real)DY;
        case 2:
            return (Real)DZ;
        default:
            fprintf(stderr, "Invalid vector component: %d\n", component);
            exit(1);
    }
}

Real beta_from_k(Real k) {
    return 1.0 + (DT * NU) / (2.0 * k);
}

Real gamma_from_k(Real k) {
    Real beta = beta_from_k(k);
    return (DT * NU) / (2.0 * beta);
}

Real time_physical_coord(Real t_step) {
    return t_step * (Real)DT;
}

Real centered_physical_coord(int index, int component) {
    return (Real)index * spacing_from_component(component);
}

Real staggered_physical_coord(int index, int component) {
    return ((Real)index + 0.5) * spacing_from_component(component);
}

static inline Real boundary_increment(VectorFunction bc_velocity,
                                      Real x, Real y, Real z,
                                      Real t, int t_step, int component) {
    Real current = bc_velocity(x, y, z, t, component);
    if (t_step == 0) {
        return current;
    }
    return current - bc_velocity(x, y, z, t - (Real)DT, component);
}

/*
 * Return the boundary increment u(t_step) - u(t_step - 1).
 * On a lower face, the normal component is reconstructed from the
 * divergence-free constraint; tangential components are sampled directly.
 * At lower edges and corners the prescribed staggered value has priority.
 */
Real bc_left(VectorFunction bc_velocity,
             int i, int j, int k, int t_step, int component) {
    Real t = time_physical_coord(t_step);
    Real x = (Real)i * (Real)DX;
    Real y = (Real)j * (Real)DY;
    Real z = (Real)k * (Real)DZ;
    Real vx = x + (Real)DX / 2.0;
    Real vy = y + (Real)DY / 2.0;
    Real vz = z + (Real)DZ / 2.0;
    int lower_face_count = (i == 0) + (j == 0) + (k == 0);

    if ((unsigned int)component > 2U) {
        fprintf(stderr, "Invalid vector component: %d\n", component);
        exit(1);
    }
    if (lower_face_count == 0) {
        return 0.0;
    }

#define BC_INCREMENT(px, py, pz, comp) \
    boundary_increment(bc_velocity, (px), (py), (pz), \
                       t, t_step, (comp))

    if (lower_face_count > 1) {
        switch (component) {
            case 0:
                return BC_INCREMENT(vx, y, z, 0);
            case 1:
                return BC_INCREMENT(x, vy, z, 1);
            default:
                return BC_INCREMENT(x, y, vz, 2);
        }
    }

    if (i == 0) {
        if (component == 0) {
            Real divergence_y =
                (BC_INCREMENT(0.0, vy, z, 1) -
                 BC_INCREMENT(0.0, vy - (Real)DY, z, 1)) *
                (Real)DY_INVERSE;
            Real divergence_z =
                (BC_INCREMENT(0.0, y, vz, 2) -
                 BC_INCREMENT(0.0, y, vz - (Real)DZ, 2)) *
                (Real)DZ_INVERSE;
            return BC_INCREMENT(0.0, y, z, 0) -
                   ((Real)DX / 2.0) * (divergence_y + divergence_z);
        }
        if (component == 1) {
            return BC_INCREMENT(0.0, vy, z, 1);
        }
        return BC_INCREMENT(0.0, y, vz, 2);
    }

    if (j == 0) {
        if (component == 0) {
            return BC_INCREMENT(vx, 0.0, z, 0);
        }
        if (component == 1) {
            Real divergence_x =
                (BC_INCREMENT(vx, 0.0, z, 0) -
                 BC_INCREMENT(vx - (Real)DX, 0.0, z, 0)) *
                (Real)DX_INVERSE;
            Real divergence_z =
                (BC_INCREMENT(x, 0.0, vz, 2) -
                 BC_INCREMENT(x, 0.0, vz - (Real)DZ, 2)) *
                (Real)DZ_INVERSE;
            return BC_INCREMENT(x, 0.0, z, 1) -
                   ((Real)DY / 2.0) * (divergence_x + divergence_z);
        }
        return BC_INCREMENT(x, 0.0, vz, 2);
    }

    if (component == 0) {
        return BC_INCREMENT(vx, y, 0.0, 0);
    }
    if (component == 1) {
        return BC_INCREMENT(x, vy, 0.0, 1);
    }

    {
        Real divergence_x =
            (BC_INCREMENT(vx, y, 0.0, 0) -
             BC_INCREMENT(vx - (Real)DX, y, 0.0, 0)) *
            (Real)DX_INVERSE;
        Real divergence_y =
            (BC_INCREMENT(x, vy, 0.0, 1) -
             BC_INCREMENT(x, vy - (Real)DY, 0.0, 1)) *
            (Real)DY_INVERSE;
        return BC_INCREMENT(x, y, 0.0, 2) -
               ((Real)DZ / 2.0) * (divergence_x + divergence_y);
    }

#undef BC_INCREMENT
}

/*
 * Upper faces use the prescribed value at the physical wall.  When a point
 * belongs to more than one upper face, Z has priority over Y, then X, matching
 * the boundary overwrite order.  Lower faces retain priority at mixed edges.
 */
Real bc_right(VectorFunction bc_velocity,
              int i, int j, int k, int t_step, int component) {
    if (i == 0 || j == 0 || k == 0) {
        return bc_left(bc_velocity, i, j, k, t_step, component);
    }

    if ((unsigned int)component > 2U) {
        fprintf(stderr, "Invalid vector component: %d\n", component);
        exit(1);
    }

    Real t = time_physical_coord(t_step);
    Real x = (Real)i * (Real)DX;
    Real y = (Real)j * (Real)DY;
    Real z = (Real)k * (Real)DZ;
    Real vx = x + (Real)DX / 2.0;
    Real vy = y + (Real)DY / 2.0;
    Real vz = z + (Real)DZ / 2.0;

#define BC_INCREMENT(px, py, pz, comp) \
    boundary_increment(bc_velocity, (px), (py), (pz), \
                       t, t_step, (comp))

    if (k == DEPTH - 1) {
        switch (component) {
            case 0:
                return BC_INCREMENT(vx, y, vz, 0);
            case 1:
                return BC_INCREMENT(x, vy, vz, 1);
            default:
                return BC_INCREMENT(x, y, vz, 2);
        }
    }

    if (j == HEIGHT - 1) {
        switch (component) {
            case 0:
                return BC_INCREMENT(vx, vy, z, 0);
            case 1:
                return BC_INCREMENT(x, vy, z, 1);
            default:
                return BC_INCREMENT(x, vy, vz, 2);
        }
    }

    if (i == WIDTH - 1) {
        switch (component) {
            case 0:
                return BC_INCREMENT(vx, y, z, 0);
            case 1:
                return BC_INCREMENT(vx, vy, z, 1);
            default:
                return BC_INCREMENT(vx, y, vz, 2);
        }
    }

#undef BC_INCREMENT
    return 0.0;
}

static inline Real interior_second_derivative(const Real *restrict field,
                                              size_t index,
                                              size_t stride,
                                              Real inverse_spacing_square) {
    return (field[index - stride] -
            2.0 * field[index] +
            field[index + stride]) * inverse_spacing_square;
}

static inline Real upper_second_derivative(const Real *restrict field,
                                           size_t index,
                                           size_t stride,
                                           Real boundary_value,
                                           Real inverse_spacing_square) {
    return (field[index - stride] -
            3.0 * field[index] +
            2.0 * boundary_value) * inverse_spacing_square;
}

/* Advective form = (u . grad) U_component. */

bool valid_convective_point(int i, int j, int k, int component) {
    if (i < 1 || j < 1 || k < 1) {
        return false;
    }

    switch (component) {
        case 0:
            return i < WIDTH - 1 && j < HEIGHT && k < DEPTH;
        case 1:
            return i < WIDTH && j < HEIGHT - 1 && k < DEPTH;
        case 2:
            return i < WIDTH && j < HEIGHT && k < DEPTH - 1;
        default:
            return false;
    }
}

static Real convective_component_value(const VectorField *restrict vel,
                                       int component, int i, int j, int k,
                                       VectorFunction bc_velocity, Real time)
{
    const Real *field;
    int upper_index;
    int normal_component;
    int *coordinate;
    Real x;
    Real y;
    Real z;

    switch (component) {
        case 0:
            field = vel->v_x;
            break;
        case 1:
            field = vel->v_y;
            break;
        case 2:
            field = vel->v_z;
            break;
        default:
            return (Real)0;
    }

    if (i >= 0 && i < WIDTH && j >= 0 && j < HEIGHT &&
        k >= 0 && k < DEPTH) {
        return field[(size_t)k * (size_t)WIDTH * (size_t)HEIGHT +
                     (size_t)j * (size_t)WIDTH + (size_t)i];
    }

    if (bc_velocity == NULL || i < 0 || j < 0 || k < 0 ||
        i > WIDTH || j > HEIGHT || k > DEPTH ||
        ((i == WIDTH) + (j == HEIGHT) + (k == DEPTH) != 1)) {
        return (Real)0;
    }

    x = ((Real)i + (component == 0 ? (Real)0.5 : (Real)0)) * (Real)DX;
    y = ((Real)j + (component == 1 ? (Real)0.5 : (Real)0)) * (Real)DY;
    z = ((Real)k + (component == 2 ? (Real)0.5 : (Real)0)) * (Real)DZ;

    if (i == WIDTH) {
        upper_index = WIDTH - 1;
        normal_component = component == 0;
        coordinate = &i;
        x = ((Real)WIDTH - (Real)0.5) * (Real)DX;
    } else if (j == HEIGHT) {
        upper_index = HEIGHT - 1;
        normal_component = component == 1;
        coordinate = &j;
        y = ((Real)HEIGHT - (Real)0.5) * (Real)DY;
    } else {
        upper_index = DEPTH - 1;
        normal_component = component == 2;
        coordinate = &k;
        z = ((Real)DEPTH - (Real)0.5) * (Real)DZ;
    }

    *coordinate = normal_component ? upper_index - 1 : upper_index;
    return (Real)2 * bc_velocity(x, y, z, time, component) -
           field[(size_t)k * (size_t)WIDTH * (size_t)HEIGHT +
                 (size_t)j * (size_t)WIDTH + (size_t)i];
}

Real advective_form(int i, int j, int k, int component,
                    const VectorField *restrict vel,
                    VectorFunction bc_velocity, Real time) {
    const Real inv_4dx = (Real)0.25 * (Real)DX_INVERSE;
    const Real inv_4dy = (Real)0.25 * (Real)DY_INVERSE;
    const Real inv_4dz = (Real)0.25 * (Real)DZ_INVERSE;

    if (!valid_convective_point(i, j, k, component)) {
        return (Real)0;
    }

    /* The vast majority of points have a complete in-domain stencil.  Keep
     * their hot path free of ghost-coordinate checks and boundary callbacks. */
    if (i < WIDTH - 1 && j < HEIGHT - 1 && k < DEPTH - 1) {
        const size_t plane = (size_t)WIDTH * (size_t)HEIGHT;
        const size_t index = (size_t)k * plane + (size_t)j * WIDTH + (size_t)i;
        const Real *restrict u = vel->v_x;
        const Real *restrict v = vel->v_y;
        const Real *restrict w = vel->v_z;

        switch (component) {
            case 0: {
                const Real v_bar = (Real)0.25 *
                    (v[index] + v[index + 1] + v[index - WIDTH] +
                     v[index + 1 - WIDTH]);
                const Real w_bar = (Real)0.25 *
                    (w[index] + w[index + 1] + w[index - plane] +
                     w[index + 1 - plane]);
                return u[index] * (u[index + 1] - u[index - 1]) * inv_4dx +
                       v_bar * (u[index + WIDTH] - u[index - WIDTH]) * inv_4dy +
                       w_bar * (u[index + plane] - u[index - plane]) * inv_4dz;
            }
            case 1: {
                const Real u_bar = (Real)0.25 *
                    (u[index] + u[index + WIDTH] + u[index - 1] +
                     u[index - 1 + WIDTH]);
                const Real w_bar = (Real)0.25 *
                    (w[index] + w[index + WIDTH] + w[index - plane] +
                     w[index + WIDTH - plane]);
                return u_bar * (v[index + 1] - v[index - 1]) * inv_4dx +
                       v[index] * (v[index + WIDTH] - v[index - WIDTH]) * inv_4dy +
                       w_bar * (v[index + plane] - v[index - plane]) * inv_4dz;
            }
            case 2: {
                const Real u_bar = (Real)0.25 *
                    (u[index] + u[index + plane] + u[index - 1] +
                     u[index - 1 + plane]);
                const Real v_bar = (Real)0.25 *
                    (v[index] + v[index + plane] + v[index - WIDTH] +
                     v[index - WIDTH + plane]);
                return u_bar * (w[index + 1] - w[index - 1]) * inv_4dx +
                       v_bar * (w[index + WIDTH] - w[index - WIDTH]) * inv_4dy +
                       w[index] * (w[index + plane] - w[index - plane]) * inv_4dz;
            }
            default:
                return (Real)0;
        }
    }

    switch (component) {
        case 0: {
            Real u_center = convective_component_value(
                vel, 0, i, j, k, bc_velocity, time);
            Real dudx = convective_component_value(vel, 0, i + 1, j, k,
                                                    bc_velocity, time) -
                        convective_component_value(vel, 0, i - 1, j, k,
                                                    bc_velocity, time);
            Real dudy = convective_component_value(vel, 0, i, j + 1, k,
                                                    bc_velocity, time) -
                        convective_component_value(vel, 0, i, j - 1, k,
                                                    bc_velocity, time);
            Real dudz = convective_component_value(vel, 0, i, j, k + 1,
                                                    bc_velocity, time) -
                        convective_component_value(vel, 0, i, j, k - 1,
                                                    bc_velocity, time);

            Real v_bar = (Real)0.25 *
                (convective_component_value(vel, 1, i, j, k, bc_velocity, time) +
                 convective_component_value(vel, 1, i + 1, j, k, bc_velocity, time) +
                 convective_component_value(vel, 1, i, j - 1, k, bc_velocity, time) +
                 convective_component_value(vel, 1, i + 1, j - 1, k, bc_velocity, time));
            Real w_bar = (Real)0.25 *
                (convective_component_value(vel, 2, i, j, k, bc_velocity, time) +
                 convective_component_value(vel, 2, i + 1, j, k, bc_velocity, time) +
                 convective_component_value(vel, 2, i, j, k - 1, bc_velocity, time) +
                 convective_component_value(vel, 2, i + 1, j, k - 1, bc_velocity, time));

            return u_center * dudx * inv_4dx +
                   v_bar * dudy * inv_4dy +
                   w_bar * dudz * inv_4dz;
        }

        case 1: {
            Real v_center = convective_component_value(
                vel, 1, i, j, k, bc_velocity, time);
            Real dvdx = convective_component_value(vel, 1, i + 1, j, k,
                                                    bc_velocity, time) -
                        convective_component_value(vel, 1, i - 1, j, k,
                                                    bc_velocity, time);
            Real dvdy = convective_component_value(vel, 1, i, j + 1, k,
                                                    bc_velocity, time) -
                        convective_component_value(vel, 1, i, j - 1, k,
                                                    bc_velocity, time);
            Real dvdz = convective_component_value(vel, 1, i, j, k + 1,
                                                    bc_velocity, time) -
                        convective_component_value(vel, 1, i, j, k - 1,
                                                    bc_velocity, time);

            Real u_bar = (Real)0.25 *
                (convective_component_value(vel, 0, i, j, k, bc_velocity, time) +
                 convective_component_value(vel, 0, i, j + 1, k, bc_velocity, time) +
                 convective_component_value(vel, 0, i - 1, j, k, bc_velocity, time) +
                 convective_component_value(vel, 0, i - 1, j + 1, k, bc_velocity, time));
            Real w_bar = (Real)0.25 *
                (convective_component_value(vel, 2, i, j, k, bc_velocity, time) +
                 convective_component_value(vel, 2, i, j + 1, k, bc_velocity, time) +
                 convective_component_value(vel, 2, i, j, k - 1, bc_velocity, time) +
                 convective_component_value(vel, 2, i, j + 1, k - 1, bc_velocity, time));

            return u_bar * dvdx * inv_4dx +
                   v_center * dvdy * inv_4dy +
                   w_bar * dvdz * inv_4dz;
        }

        case 2: {
            Real w_center = convective_component_value(
                vel, 2, i, j, k, bc_velocity, time);
            Real dwdx = convective_component_value(vel, 2, i + 1, j, k,
                                                    bc_velocity, time) -
                        convective_component_value(vel, 2, i - 1, j, k,
                                                    bc_velocity, time);
            Real dwdy = convective_component_value(vel, 2, i, j + 1, k,
                                                    bc_velocity, time) -
                        convective_component_value(vel, 2, i, j - 1, k,
                                                    bc_velocity, time);
            Real dwdz = convective_component_value(vel, 2, i, j, k + 1,
                                                    bc_velocity, time) -
                        convective_component_value(vel, 2, i, j, k - 1,
                                                    bc_velocity, time);

            Real u_bar = (Real)0.25 *
                (convective_component_value(vel, 0, i, j, k, bc_velocity, time) +
                 convective_component_value(vel, 0, i, j, k + 1, bc_velocity, time) +
                 convective_component_value(vel, 0, i - 1, j, k, bc_velocity, time) +
                 convective_component_value(vel, 0, i - 1, j, k + 1, bc_velocity, time));
            Real v_bar = (Real)0.25 *
                (convective_component_value(vel, 1, i, j, k, bc_velocity, time) +
                 convective_component_value(vel, 1, i, j, k + 1, bc_velocity, time) +
                 convective_component_value(vel, 1, i, j - 1, k, bc_velocity, time) +
                 convective_component_value(vel, 1, i, j - 1, k + 1, bc_velocity, time));

            return u_bar * dwdx * inv_4dx +
                   v_bar * dwdy * inv_4dy +
                   w_center * dwdz * inv_4dz;
        }

        default:
            return (Real)0;
    }
}

Real divergence_form(int i, int j, int k, int component,
                     const VectorField *restrict vel,
                     VectorFunction bc_velocity, Real time) {
    const Real inv_2dx = (Real)0.5 * (Real)DX_INVERSE;
    const Real inv_2dy = (Real)0.5 * (Real)DY_INVERSE;
    const Real inv_2dz = (Real)0.5 * (Real)DZ_INVERSE;

    if (!valid_convective_point(i, j, k, component)) {
        return (Real)0;
    }

    if (i < WIDTH - 1 && j < HEIGHT - 1 && k < DEPTH - 1) {
        const size_t plane = (size_t)WIDTH * (size_t)HEIGHT;
        const size_t index = (size_t)k * plane + (size_t)j * WIDTH + (size_t)i;
        const Real *restrict u = vel->v_x;
        const Real *restrict v = vel->v_y;
        const Real *restrict w = vel->v_z;

        switch (component) {
            case 0: {
                const Real u_i = (Real)0.5 * (u[index - 1] + u[index]);
                const Real u_ip1 = (Real)0.5 * (u[index] + u[index + 1]);
                const Real u_jp = (Real)0.5 * (u[index] + u[index + WIDTH]);
                const Real u_jm = (Real)0.5 * (u[index] + u[index - WIDTH]);
                const Real v_jp = (Real)0.5 * (v[index] + v[index + 1]);
                const Real v_jm = (Real)0.5 *
                    (v[index - WIDTH] + v[index + 1 - WIDTH]);
                const Real u_kp = (Real)0.5 * (u[index] + u[index + plane]);
                const Real u_km = (Real)0.5 * (u[index] + u[index - plane]);
                const Real w_kp = (Real)0.5 * (w[index] + w[index + 1]);
                const Real w_km = (Real)0.5 *
                    (w[index - plane] + w[index + 1 - plane]);
                return (u_ip1 * u_ip1 - u_i * u_i) * inv_2dx +
                       (u_jp * v_jp - u_jm * v_jm) * inv_2dy +
                       (u_kp * w_kp - u_km * w_km) * inv_2dz;
            }
            case 1: {
                const Real u_ip = (Real)0.5 * (u[index] + u[index + WIDTH]);
                const Real u_im = (Real)0.5 *
                    (u[index - 1] + u[index - 1 + WIDTH]);
                const Real v_ip = (Real)0.5 * (v[index] + v[index + 1]);
                const Real v_im = (Real)0.5 * (v[index] + v[index - 1]);
                const Real v_j = (Real)0.5 * (v[index - WIDTH] + v[index]);
                const Real v_jp1 = (Real)0.5 * (v[index] + v[index + WIDTH]);
                const Real v_kp = (Real)0.5 * (v[index] + v[index + plane]);
                const Real v_km = (Real)0.5 * (v[index] + v[index - plane]);
                const Real w_kp = (Real)0.5 * (w[index] + w[index + WIDTH]);
                const Real w_km = (Real)0.5 *
                    (w[index - plane] + w[index + WIDTH - plane]);
                return (u_ip * v_ip - u_im * v_im) * inv_2dx +
                       (v_jp1 * v_jp1 - v_j * v_j) * inv_2dy +
                       (v_kp * w_kp - v_km * w_km) * inv_2dz;
            }
            case 2: {
                const Real u_ip = (Real)0.5 * (u[index] + u[index + plane]);
                const Real u_im = (Real)0.5 *
                    (u[index - 1] + u[index - 1 + plane]);
                const Real w_ip = (Real)0.5 * (w[index] + w[index + 1]);
                const Real w_im = (Real)0.5 * (w[index] + w[index - 1]);
                const Real v_jp = (Real)0.5 * (v[index] + v[index + plane]);
                const Real v_jm = (Real)0.5 *
                    (v[index - WIDTH] + v[index - WIDTH + plane]);
                const Real w_jp = (Real)0.5 * (w[index] + w[index + WIDTH]);
                const Real w_jm = (Real)0.5 * (w[index] + w[index - WIDTH]);
                const Real w_k = (Real)0.5 * (w[index - plane] + w[index]);
                const Real w_kp1 = (Real)0.5 * (w[index] + w[index + plane]);
                return (u_ip * w_ip - u_im * w_im) * inv_2dx +
                       (v_jp * w_jp - v_jm * w_jm) * inv_2dy +
                       (w_kp1 * w_kp1 - w_k * w_k) * inv_2dz;
            }
            default:
                return (Real)0;
        }
    }

    switch (component) {
        case 0: {
            Real u_i = (Real)0.5 *
                (convective_component_value(
                     vel, 0, i - 1, j, k, bc_velocity, time) +
                 convective_component_value(
                     vel, 0, i, j, k, bc_velocity, time));
            Real u_ip1 = (Real)0.5 *
                (convective_component_value(
                     vel, 0, i, j, k, bc_velocity, time) +
                 convective_component_value(
                     vel, 0, i + 1, j, k, bc_velocity, time));
            Real u_jp = (Real)0.5 *
                (convective_component_value(
                     vel, 0, i, j, k, bc_velocity, time) +
                 convective_component_value(
                     vel, 0, i, j + 1, k, bc_velocity, time));
            Real u_jm = (Real)0.5 *
                (convective_component_value(
                     vel, 0, i, j, k, bc_velocity, time) +
                 convective_component_value(
                     vel, 0, i, j - 1, k, bc_velocity, time));
            Real v_jp = (Real)0.5 *
                (convective_component_value(
                     vel, 1, i, j, k, bc_velocity, time) +
                 convective_component_value(
                     vel, 1, i + 1, j, k, bc_velocity, time));
            Real v_jm = (Real)0.5 *
                (convective_component_value(
                     vel, 1, i, j - 1, k, bc_velocity, time) +
                 convective_component_value(
                     vel, 1, i + 1, j - 1, k, bc_velocity, time));
            Real u_kp = (Real)0.5 *
                (convective_component_value(
                     vel, 0, i, j, k, bc_velocity, time) +
                 convective_component_value(
                     vel, 0, i, j, k + 1, bc_velocity, time));
            Real u_km = (Real)0.5 *
                (convective_component_value(
                     vel, 0, i, j, k, bc_velocity, time) +
                 convective_component_value(
                     vel, 0, i, j, k - 1, bc_velocity, time));
            Real w_kp = (Real)0.5 *
                (convective_component_value(
                     vel, 2, i, j, k, bc_velocity, time) +
                 convective_component_value(
                     vel, 2, i + 1, j, k, bc_velocity, time));
            Real w_km = (Real)0.5 *
                (convective_component_value(
                     vel, 2, i, j, k - 1, bc_velocity, time) +
                 convective_component_value(
                     vel, 2, i + 1, j, k - 1, bc_velocity, time));

            return (u_ip1 * u_ip1 - u_i * u_i) * inv_2dx +
                   (u_jp * v_jp - u_jm * v_jm) * inv_2dy +
                   (u_kp * w_kp - u_km * w_km) * inv_2dz;
        }

        case 1: {
            Real u_ip = (Real)0.5 *
                (convective_component_value(
                     vel, 0, i, j, k, bc_velocity, time) +
                 convective_component_value(
                     vel, 0, i, j + 1, k, bc_velocity, time));
            Real u_im = (Real)0.5 *
                (convective_component_value(
                     vel, 0, i - 1, j, k, bc_velocity, time) +
                 convective_component_value(
                     vel, 0, i - 1, j + 1, k, bc_velocity, time));
            Real v_ip = (Real)0.5 *
                (convective_component_value(
                     vel, 1, i, j, k, bc_velocity, time) +
                 convective_component_value(
                     vel, 1, i + 1, j, k, bc_velocity, time));
            Real v_im = (Real)0.5 *
                (convective_component_value(
                     vel, 1, i, j, k, bc_velocity, time) +
                 convective_component_value(
                     vel, 1, i - 1, j, k, bc_velocity, time));
            Real v_j = (Real)0.5 *
                (convective_component_value(
                     vel, 1, i, j - 1, k, bc_velocity, time) +
                 convective_component_value(
                     vel, 1, i, j, k, bc_velocity, time));
            Real v_jp1 = (Real)0.5 *
                (convective_component_value(
                     vel, 1, i, j, k, bc_velocity, time) +
                 convective_component_value(
                     vel, 1, i, j + 1, k, bc_velocity, time));
            Real v_kp = (Real)0.5 *
                (convective_component_value(
                     vel, 1, i, j, k, bc_velocity, time) +
                 convective_component_value(
                     vel, 1, i, j, k + 1, bc_velocity, time));
            Real v_km = (Real)0.5 *
                (convective_component_value(
                     vel, 1, i, j, k, bc_velocity, time) +
                 convective_component_value(
                     vel, 1, i, j, k - 1, bc_velocity, time));
            Real w_kp = (Real)0.5 *
                (convective_component_value(
                     vel, 2, i, j, k, bc_velocity, time) +
                 convective_component_value(
                     vel, 2, i, j + 1, k, bc_velocity, time));
            Real w_km = (Real)0.5 *
                (convective_component_value(
                     vel, 2, i, j, k - 1, bc_velocity, time) +
                 convective_component_value(
                     vel, 2, i, j + 1, k - 1, bc_velocity, time));

            return (u_ip * v_ip - u_im * v_im) * inv_2dx +
                   (v_jp1 * v_jp1 - v_j * v_j) * inv_2dy +
                   (v_kp * w_kp - v_km * w_km) * inv_2dz;
        }

        case 2: {
            Real u_ip = (Real)0.5 *
                (convective_component_value(
                     vel, 0, i, j, k, bc_velocity, time) +
                 convective_component_value(
                     vel, 0, i, j, k + 1, bc_velocity, time));
            Real u_im = (Real)0.5 *
                (convective_component_value(
                     vel, 0, i - 1, j, k, bc_velocity, time) +
                 convective_component_value(
                     vel, 0, i - 1, j, k + 1, bc_velocity, time));
            Real w_ip = (Real)0.5 *
                (convective_component_value(
                     vel, 2, i, j, k, bc_velocity, time) +
                 convective_component_value(
                     vel, 2, i + 1, j, k, bc_velocity, time));
            Real w_im = (Real)0.5 *
                (convective_component_value(
                     vel, 2, i, j, k, bc_velocity, time) +
                 convective_component_value(
                     vel, 2, i - 1, j, k, bc_velocity, time));
            Real v_jp = (Real)0.5 *
                (convective_component_value(
                     vel, 1, i, j, k, bc_velocity, time) +
                 convective_component_value(
                     vel, 1, i, j, k + 1, bc_velocity, time));
            Real v_jm = (Real)0.5 *
                (convective_component_value(
                     vel, 1, i, j - 1, k, bc_velocity, time) +
                 convective_component_value(
                     vel, 1, i, j - 1, k + 1, bc_velocity, time));
            Real w_jp = (Real)0.5 *
                (convective_component_value(
                     vel, 2, i, j, k, bc_velocity, time) +
                 convective_component_value(
                     vel, 2, i, j + 1, k, bc_velocity, time));
            Real w_jm = (Real)0.5 *
                (convective_component_value(
                     vel, 2, i, j, k, bc_velocity, time) +
                 convective_component_value(
                     vel, 2, i, j - 1, k, bc_velocity, time));
            Real w_k = (Real)0.5 *
                (convective_component_value(
                     vel, 2, i, j, k - 1, bc_velocity, time) +
                 convective_component_value(
                     vel, 2, i, j, k, bc_velocity, time));
            Real w_kp1 = (Real)0.5 *
                (convective_component_value(
                     vel, 2, i, j, k, bc_velocity, time) +
                 convective_component_value(
                     vel, 2, i, j, k + 1, bc_velocity, time));

            return (u_ip * w_ip - u_im * w_im) * inv_2dx +
                   (v_jp * w_jp - v_jm * w_jm) * inv_2dy +
                   (w_kp1 * w_kp1 - w_k * w_k) * inv_2dz;
        }

        default:
            return (Real)0;
    }
}

Real convective_term(int i, int j, int k, int component,
                     const VectorField *restrict vel,
                     VectorFunction bc_velocity, Real time) {
    return advective_form(i, j, k, component, vel, bc_velocity, time) +
           divergence_form(i, j, k, component, vel, bc_velocity, time);
}

/*
 * Momentum source at (t_step - 1/2) DT:
 *
 *   g = f - grad(p*) - (NU / K) u
 *       + NU (Dxx(eta) + Dyy(zeta) + Dzz(u)).
 *
 * The support excludes lower faces and the upper face normal to the velocity
 * component.  Tangential upper faces use a Dirichlet ghost value.
 */
Real g_value(int i, int j, int k, int t_step, Real k_i,
             const SolverMemState *solver_mem_state,
             const Data *data, int component) {
    size_t plane_size = (size_t)WIDTH * HEIGHT;
    size_t index = (size_t)k * plane_size +
                   (size_t)j * WIDTH +
                   (size_t)i;
    Real forcing_time = ((Real)t_step - 0.5) * (Real)DT;
    Real velocity_time = ((Real)t_step - 1.0) * (Real)DT;
    Real upper_x = ((Real)WIDTH - 0.5) * (Real)DX;
    Real upper_y = ((Real)HEIGHT - 0.5) * (Real)DY;
    Real upper_z = ((Real)DEPTH - 0.5) * (Real)DZ;
    Real x = (Real)i * (Real)DX;
    Real y = (Real)j * (Real)DY;
    Real z = (Real)k * (Real)DZ;
    const Real *restrict eta;
    const Real *restrict zeta;
    const Real *restrict velocity;
    const Real *restrict pressure = solver_mem_state->pressure_star.v;
    Real pressure_gradient;
    Real convection;
    Real laplacian_x;
    Real laplacian_y;
    Real laplacian_z;

    switch (component) {
        case 0:
            if (i < 1 || i >= WIDTH - 1 ||
                j < 1 || j >= HEIGHT ||
                k < 1 || k >= DEPTH) {
                return 0.0;
            }

            eta = solver_mem_state->eta.v_x;
            zeta = solver_mem_state->zeta.v_x;
            velocity = solver_mem_state->u.v_x;
            x += (Real)DX / 2.0;

            pressure_gradient =
                (pressure[index + 1] - pressure[index]) *
                (Real)DX_INVERSE;
            laplacian_x = interior_second_derivative(
                eta, index, 1, (Real)DX_INVERSE_SQUARE);
            laplacian_y = (j == HEIGHT - 1)
                ? upper_second_derivative(
                      zeta, index, WIDTH,
                      data->bc_velocity(x, upper_y, z,
                                        velocity_time, component),
                      (Real)DY_INVERSE_SQUARE)
                : interior_second_derivative(
                      zeta, index, WIDTH, (Real)DY_INVERSE_SQUARE);
            laplacian_z = (k == DEPTH - 1)
                ? upper_second_derivative(
                      velocity, index, plane_size,
                      data->bc_velocity(x, y, upper_z,
                                        velocity_time, component),
                      (Real)DZ_INVERSE_SQUARE)
                : interior_second_derivative(
                      velocity, index, plane_size,
                      (Real)DZ_INVERSE_SQUARE);
            break;

        case 1:
            if (i < 1 || i >= WIDTH ||
                j < 1 || j >= HEIGHT - 1 ||
                k < 1 || k >= DEPTH) {
                return 0.0;
            }

            eta = solver_mem_state->eta.v_y;
            zeta = solver_mem_state->zeta.v_y;
            velocity = solver_mem_state->u.v_y;
            y += (Real)DY / 2.0;

            pressure_gradient =
                (pressure[index + WIDTH] - pressure[index]) *
                (Real)DY_INVERSE;
            laplacian_x = (i == WIDTH - 1)
                ? upper_second_derivative(
                      eta, index, 1,
                      data->bc_velocity(upper_x, y, z,
                                        velocity_time, component),
                      (Real)DX_INVERSE_SQUARE)
                : interior_second_derivative(
                      eta, index, 1, (Real)DX_INVERSE_SQUARE);
            laplacian_y = interior_second_derivative(
                zeta, index, WIDTH, (Real)DY_INVERSE_SQUARE);
            laplacian_z = (k == DEPTH - 1)
                ? upper_second_derivative(
                      velocity, index, plane_size,
                      data->bc_velocity(x, y, upper_z,
                                        velocity_time, component),
                      (Real)DZ_INVERSE_SQUARE)
                : interior_second_derivative(
                      velocity, index, plane_size,
                      (Real)DZ_INVERSE_SQUARE);
            break;

        case 2:
            if (i < 1 || i >= WIDTH ||
                j < 1 || j >= HEIGHT ||
                k < 1 || k >= DEPTH - 1) {
                return 0.0;
            }

            eta = solver_mem_state->eta.v_z;
            zeta = solver_mem_state->zeta.v_z;
            velocity = solver_mem_state->u.v_z;
            z += (Real)DZ / 2.0;

            pressure_gradient =
                (pressure[index + plane_size] - pressure[index]) *
                (Real)DZ_INVERSE;
            laplacian_x = (i == WIDTH - 1)
                ? upper_second_derivative(
                      eta, index, 1,
                      data->bc_velocity(upper_x, y, z,
                                        velocity_time, component),
                      (Real)DX_INVERSE_SQUARE)
                : interior_second_derivative(
                      eta, index, 1, (Real)DX_INVERSE_SQUARE);
            laplacian_y = (j == HEIGHT - 1)
                ? upper_second_derivative(
                      zeta, index, WIDTH,
                      data->bc_velocity(x, upper_y, z,
                                        velocity_time, component),
                      (Real)DY_INVERSE_SQUARE)
                : interior_second_derivative(
                      zeta, index, WIDTH, (Real)DY_INVERSE_SQUARE);
            laplacian_z = interior_second_derivative(
                velocity, index, plane_size,
                (Real)DZ_INVERSE_SQUARE);
            break;

        default:
            return 0.0;
    }

    convection = convective_term(i, j, k, component,
                                 &solver_mem_state->u_star,
                                 data->bc_velocity, forcing_time);

    return data->forcing_fn(x, y, z, forcing_time, component) -
           pressure_gradient -
           convection -
           ((Real)NU / k_i) * velocity[index] +
           (Real)NU * (laplacian_x + laplacian_y + laplacian_z);
}
