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
| `Complex.hpp` | complex numbers: `+ − × ÷`, conjugate, modulus, argument, polar form, `exp`, `log`, `sqrt`, `pow`, `sin`, `cos`, n-th roots, `3.0_i` literal, `std::complex` interop |
| `Polynomial.hpp` | evaluation (real and complex), derivative, antiderivative, definite integral, `+ − ×`, all real and complex roots |
| `Interpolation.hpp` | `Interpolator1D` with step, linear, natural cubic spline, Akima and monotone PCHIP schemes. Each provides a derivative, an exact integral and three extrapolation modes. `Table2D` covers tabulated 2-D data |
| `Ode.hpp` | fixed-step Euler, midpoint, Heun and RK4, plus an adaptive Dormand-Prince 5(4) step with error control |
| `Constants.hpp` | π, degree/radian conversion, `clamp`, `smoothstep`, `lerp` |

```cpp
#include <mathrix/mathrix.hpp>
using namespace mathrix;
```

## Linear systems: $A\mathbf{x} = \mathbf{b}$

A system of $n$ linear equations in $n$ unknowns, for example

```math
\left\{
\begin{aligned}
 2x + y - z &= 8 \\
-3x - y + 2z &= -11 \\
-2x + y + 2z &= -3
\end{aligned}
\right.
```

is written in matrix form as

```math
\underbrace{\begin{bmatrix} 2 & 1 & -1 \\ -3 & -1 & 2 \\ -2 & 1 & 2 \end{bmatrix}}_{A}
\underbrace{\begin{bmatrix} x \\ y \\ z \end{bmatrix}}_{\mathbf{x}}
=
\underbrace{\begin{bmatrix} 8 \\ -11 \\ -3 \end{bmatrix}}_{\mathbf{b}}
```

`Matrix::solve` uses Gaussian elimination with partial pivoting and finds

```math
\mathbf{x} = A^{-1}\mathbf{b} = \begin{bmatrix} 2 \\ 3 \\ -1 \end{bmatrix}
\qquad\Longrightarrow\qquad x = 2,\; y = 3,\; z = -1
```

```cpp
Matrix A{{ 2,  1, -1},
         {-3, -1,  2},
         {-2,  1,  2}};
Matrix b = Matrix::column({8, -11, -3});
Matrix x = A.solve(b);            // [2, 3, -1]
double d = A.determinant();       // -1, non-zero so the solution is unique
Matrix Ainv = A.inverse();
```

The same call works for any $n$. If $\det A = 0$ the system has no unique solution and `solve` throws `std::domain_error`.
[`examples/linear_system.cpp`](examples/linear_system.cpp) solves a 6-unknown system and checks the residual
$\lVert A\mathbf{x} - \mathbf{b} \rVert \approx 10^{-15}$.

## Complex numbers

A complex number $z = a + bi$ with $i^2 = -1$ also has a polar form $z = r e^{i\theta}$, where $r = |z| = \sqrt{a^2 + b^2}$ and
$\theta = \arg z$. For $z_1 = 3 + 4i$ and $z_2 = 1 - 2i$:

```math
\begin{aligned}
z_1 + z_2 &= 4 + 2i, &\qquad z_1 z_2 &= (3 + 4i)(1 - 2i) = 11 - 2i, \\
\frac{z_1}{z_2} &= \frac{(3 + 4i)(1 + 2i)}{(1 - 2i)(1 + 2i)} = -1 + 2i, &\qquad |z_1| &= \sqrt{3^2 + 4^2} = 5
\end{aligned}
```

The elementary functions follow Euler's formula $e^{i\theta} = \cos\theta + i \sin\theta$:

```math
e^{i\pi} = -1, \qquad \sqrt{-4} = 2i, \qquad i^{\,i} = e^{-\pi/2} \approx 0.2079, \qquad
\sqrt[n]{z} = |z|^{1/n}\, e^{\,i(\arg z + 2\pi k)/n}, \quad k = 0, \dots, n-1
```

```cpp
using namespace mathrix::literals;

Complex z1 = 3.0 + 4.0_i;
Complex z2{1, -2};
Complex p = z1 * z2;                       // 11 - 2i
Complex q = z1 / z2;                       // -1 + 2i
double r = z1.abs(), theta = z1.arg();     // 5, 0.927 rad
Complex e = chop(exp(Complex(0, kPi)));    // -1   (chop removes 1e-16 round-off)
auto w = nthRoots(Complex(16), 4);         // 2, 2i, -2, -2i
std::complex<double> s = z1.toStd();       // interop with the standard library
```

## Polynomials

Polynomials are stored by ascending coefficients, $p(x) = c_0 + c_1 x + \dots + c_n x^n$. `roots()` returns all $n$ roots,
real and complex (Durand-Kerner iteration with Newton polishing):

```math
p(x) = x^3 - 2x^2 + 4x - 8 = (x - 2)(x^2 + 4)
\quad\Longrightarrow\quad x_1 = 2,\; x_{2,3} = \pm 2i
```

```cpp
Polynomial p{-8, 4, -2, 1};               // x^3 - 2x^2 + 4x - 8
auto roots = p.roots();                    // -2i, 2i, 2
Polynomial dp = p.derivative();            // 3x^2 - 4x + 4
double area = p.integrate(0.0, 2.0);       // definite integral, exact
Complex v = p(Complex(0, 2));              // 0
```

## Interpolation and ODEs

```cpp
// Monotone interpolation of a thrust curve and its exact integral (total impulse).
Interpolator1D thrust(t, F, Interp::Pchip, Extrapolation::Zero);
double impulse = thrust.integrate(0.0, t.back());

// Harmonic oscillator x'' = -x with RK4.
using S = VecN<2>;
auto f = [](double, const S& y) { S d; d[0] = y[1]; d[1] = -y[0]; return d; };
S y; y[0] = 1.0;
y = fixedStep(Integrator::RK4, f, 0.0, y, 0.01);

// Attitude: rotate the body x axis 90 degrees about z.
Quat q = Quat::fromAxisAngle(Vec3::unitZ(), deg2rad(90.0));
Vec3 v = q.rotate(Vec3::unitX());          // (0, 1, 0)
```

## Examples

| Program | Shows |
|---|---|
| [`examples/linear_system.cpp`](examples/linear_system.cpp) | $A\mathbf{x} = \mathbf{b}$ with 3 and 6 unknowns, residual check, singular systems |
| [`examples/complex.cpp`](examples/complex.cpp) | complex arithmetic, Euler's identity, n-th roots, complex polynomial roots |
| [`examples/example.cpp`](examples/example.cpp) | quaternions, interpolation, RK4, polynomials in a few lines |

## CMake

```cmake
# as a subdirectory / git submodule
add_subdirectory(external/MathRix)
target_link_libraries(my_app PRIVATE MathRix::MathRix)

# or fetched
include(FetchContent)
FetchContent_Declare(MathRix GIT_REPOSITORY https://github.com/KaanBahaSever/MathRix.git GIT_TAG v0.2.0)
FetchContent_MakeAvailable(MathRix)

# or installed: find_package(MathRix REQUIRED)
```

## Building the tests

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build
```

## Roadmap

- Complex matrices and linear systems
- Eigenvalues of symmetric matrices, Cholesky and QR decompositions
- Fixed-size `MatN<R, C>` for small embedded problems

## License

MIT. See [LICENSE](LICENSE).
