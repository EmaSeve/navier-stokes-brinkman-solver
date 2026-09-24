# Navier–Stokes–Brinkman Solver

The `serial` and `convective` branches use the
same build procedure and produce a serial executable.

## Build

### Requirements

- a C11 compiler;
- GNU Make.

To enable the explicit SIMD momentum kernels, compile with:

```sh
make solver SIMD=1 ZETA_SIMD_VECTORS=16 U_SIMD_VECTORS=16
```

`ZETA_SIMD_VECTORS` and `U_SIMD_VECTORS` control the number of SIMD vectors
processed in each block and can be tuned independently.

Compile and run all programs in `test/` with:

```sh
make tests SIMD=1
./build/tests/moving_sphere
./build/tests/channel_obstacle
```
