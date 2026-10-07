# MathRix

MathRix is an open-source, header-only C++17 library for the mathematics that simulation and engineering code needs every day.
It has no dependencies: copy the `include/` folder or add it with CMake.
It is the math layer of [Rocket-Up](https://github.com/KaanBahaSever/Rocket-Up), the rocket flight simulator.

## Features

| Header | Contents |
|---|---|
| `Vec3.hpp` | 3-D vectors: arithmetic, `dot`, `cross`, `norm`, `normalized` |
| `Mat3.hpp` | 3×3 matrices: product, transpose, determinant, inverse |
| `Quat.hpp` | unit quaternions: `rotate`, `inverseRotate`, `fromAxisAngle`, `fromTwoVectors`, `derivative(ω)` for attitude integration, `toMatrix` |
| `Matrix.hpp` | dynamic M×N matrices: addition, subtraction, multiplication, transpose, trace, determinant (LU with pivoting), `inverse`, `solve` |
| `Polynomial.hpp` | evaluation (Horner), derivative, antiderivative, definite integral, `+ − ×`, root finding |
| `Interpolation.hpp` | `Interpolator1D` with step, linear, natural cubic spline, Akima and monotone PCHIP schemes. Each provides a derivative, an exact integral and three extrapolation modes. `Table2D` covers tabulated 2-D data |
| `Ode.hpp` | fixed-step Euler, midpoint, Heun and RK4, plus an adaptive Dormand-Prince 5(4) step with error control |
| `Constants.hpp` | π, degree/radian conversion, `clamp`, `smoothstep`, `lerp` |

All cubic interpolation schemes share one Hermite representation, so evaluation, derivatives and integrals use one code path.

## Usage

```cpp
#include <mathrix/mathrix.hpp>
using namespace mathrix;

Quat q = Quat::fromAxisAngle(Vec3::unitZ(), deg2rad(90.0));
Vec3 v = q.rotate(Vec3::unitX());                               // (0, 1, 0)

Interpolator1D thrust(t, F, Interp::Pchip, Extrapolation::Zero); // no overshoot
double impulse = thrust.integrate(0.0, t.back());

// Linear system A x = b with n unknowns (Gaussian elimination with partial pivoting)
Matrix A{{ 2,  1, -1},
         {-3, -1,  2},
         {-2,  1,  2}};
Matrix b = Matrix::column({8, -11, -3});
Matrix x = A.solve(b);                                          // [2, 3, -1]

Polynomial p{1, 0, -3, 1};                                      // x^3 - 3x^2 + 1
double r = p.rootIn(0.0, 1.0);

using S = VecN<2>;
auto f = [](double, const S& y) { S d; d[0] = y[1]; d[1] = -y[0]; return d; };
S y; y[0] = 1.0;
y = fixedStep(Integrator::RK4, f, 0.0, y, 0.01);
```

### CMake

```cmake
# as a subdirectory / git submodule
add_subdirectory(external/MathRix)
target_link_libraries(my_app PRIVATE MathRix::MathRix)

# or fetched
include(FetchContent)
FetchContent_Declare(MathRix GIT_REPOSITORY https://github.com/KaanBahaSever/MathRix.git GIT_TAG v0.1.0)
FetchContent_MakeAvailable(MathRix)

# or installed: find_package(MathRix REQUIRED)
```

See [`examples/linear_system.cpp`](examples/linear_system.cpp) for a complete n-unknown example, including the residual check
and singular systems.

## Building the tests

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build
```

## Roadmap

- Complex-number helpers on top of `std::complex`
- Eigenvalues of symmetric matrices, Cholesky and QR decompositions
- Fixed-size `MatN<R, C>` for small embedded problems

## License

MIT. See [LICENSE](LICENSE).
