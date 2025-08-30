# ieee754_seq

A C++ class template defining a sequence of IEEE-754 floating point values.

## Description

The following code

```cpp
ieee754_seq<float, N, Emin, Emax> s;
```

defines a sequence `s` of radix-2 IEEE-754 floats

$$
x_i=m_i \times 2^{E_i}, \quad i=0, 1, \dots, M= 2^N\left(E_{max} - E_{min}\right)
$$

where $m_i$ is the $N$-bit "mantissa"

$$
m_i = 1 + (i \bmod 2^N) / 2^N
$$

and

$$
E_i = E_{min} + (i \div 2^N).
$$

The range of $x_i$ is 

$$
x_0=2^{E_{min}} \leq x_i \leq x_M = 2^{E_{max}}.
$$

For example, the sequence `ieee754_seq<float,2,0,2>` is

| $i$   | 0   | 1    | 2   | 3    | 4   | 5   | 6   | 7   | 8   |
| ----- | --- | ---- | --- | ---- | --- | --- | --- | --- | --- |
| $x_i$ | 1   | 1.25 | 1.5 | 1.75 | 2   | 2.5 | 3   | 3.5 | 4   |

The floating-pont numbers are not stored in the `ieee754_seq` object. They are generated on the fly by fast bit-manipulation operations on the integer index $i$.

Conversely, a floating-point number $x$, within the range, can be converter by bit manipulation to the index $i$, such that $x_i \leq x < x_{i+1}$.

These bit-wise conversions are possible because the IEEE-754 standard is used in most computer systems, thus, the internal representation of a floating-point number is essentially very similar to the one shown above.

## Purpose

The `ieee754_seq` number sequence is quasi log spaced. Thus it can be employed for optimized log-interpolation schemes with fast lookup and fewer calls to `log()`.

This was originally proposed by Yuan et al NIMB83(1993) p.413 for the interpolation of scattering cross sections tabulated over a wide range of projectile energy.

The implentation here is based on the program Corteo, written by Francois Schiettekatte and released under the GNU GPL.
The original corteo source code can be found here:
http://www.lps.umontreal.ca/~schiette/index.php?n=Recherche.Corteo

## Usage

Just include the single header file in your c++ code
```cpp
#include <"ieee754_seq.h">

ieee754_seq<float, 3, -10, 10> s; // 160 points from ~1e-3 to ~1e3
```

Iteration over the number sequence can be done in various ways:
```cpp
// for-loop with int
for (int i = 0; i < s.size(); ++i)
        cout << i << ' ' << s[i] << endl;

// for-loop with iterator
for (auto i = s.begin(); i < s.end(); ++i)
    cout << i << ' ' << *i << ' ' << i.log2v() << endl;

// range for
for (const float &x : s)
    cout << x << endl;
```
When using an `ieee754_seq::iterator` the dereference operator `*` returns the floating point value. Further, the function `i.log2v()` returns the base-2 logarithm without calling `std::log2()`.

### Interpolation

The `ieee754_seq` class defines 2 interpolator objects for log-lin and log-log interpolation.

```cpp
typedef ieee754_seq<float, 3, 0, 2> seq_t;

seq_t s;
seq_t::log_log_interp interp;

std::vector<float> y(s.size());

for (int i = 0; i < y.size(); ++i) y[i] = some_func(s[i]);

float x = 1.3333f;
float y_i = interp(x); // interpolate y at x

```


