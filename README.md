# Optimized N-dimensional Hilbert curve generation in C++

## Features

* *Speed*: Able to compute `Curve()` at `> 1GB/s` (tested on my
  faithful 7 year-old i5-2500K system with `N=2`, `K=15`,
  `Int=int16_t`).
* *Portability*: All the functions are standalone and should be easily
  portable to most other languages.
* *Small code size*: Currently, `Curve()` only has 25 non-empty
  lines. The functions have also been made non-recursive.
* *Compile-time evaluation*: All functions are `constexpr`. See
  example usage in `example.cpp`.
* *Full git history*: Everything has been left in, no matter how
  embarassing my mistakes may be.

## Installation

Simply add `hilbert.hpp` to your project.

## Usage

There are 3 functions for manipulating Hilbert curves:

1. `void Curve(UInt N, UInt K, Int curve[])`: Computes the `K`'th step
   of an `N` dimensional Hilbert curve.  The result will be stored in
   `curve`, which must have space for `2^(N\*K)` `N`-dimensional
   vectors, for a total of `N\*2^(N*K)` integers.
2. `void IToV(UInt N, UInt K, UInt i, Int v[])`: Computes the `i`'th
   vector in `Curve(N, K)`.  The result will be stored in `v`.
   Requires `i < 2^(N*K)`.
3. `UInt VToI(UInt N, UInt K, Int v[])`: Computes the index that `v`
   would have in the `K`'th step of an `N` dimensional Hilbert curve.
   Behavior is undefined if `v` is not a point on the curve.  `v` will
   be zeroed when this function returns.

`Curve()`, `IToV()`, and `VToI()` all operate on Hilbert curves
centered at the origin with points separated a distance of 2.  For
example, the 2nd iteration of a 1D Hilbert curve would have points at
`[{-3}, {-1}, {1}, {3}]`.  This data may be more useful based at 0
with a distance 1 between points.  In this format, the curve would
have points at `[{0}, {1}, {2}, {3}]`.  `OffsetV()` and `CenterV()`
converts between these formats.

1. `void OffsetV(UInt N, UInt K, const Int cv[], Int ov[])`: Scales
   `cv` down by a factor of 2 and shifts it to lie in the first
   orthant.  Stores the result in `ov`.  `cv` may point to the same
   vector as `ov`.
2. `void CenterV(UInt N, UInt K, const Int ov[], Int cv[])`: Shifts
   `ov` to be centered at the origin and scales it up by a factor
   of 2.  Stores the result in `cv`.  `ov` may point to the same
   vector as `cv`.

### Calling the functions

The above functions live in a static class called `Hilbert`. To call
`Curve()`, for example, you might use:

```
// Compute the K'th iteration of an N-dimensional Hilbert curve.
int curve[1 << N * K][N];
Hilbert<>::Curve(N, K, curve[0]);
```

### Template versions of each function

Each function takes the number of dimensions `N` and the iteration `K`
as parameters.  If `N` and/or `K` are known at compile time, template
functions are available which can help the compiler generate more
optimized code. For example, the 4 versions of `Curve()` are listed
below.

```
// Compute the K'th iteration of an N-dimensional Hilbert curve.
constexpr int N = 2;
constexpr int K = 3;
int curve[1 << N * K][N];

// Untemplated version.
// Use this version if N and K are not known at compile time.
Hilbert<>::Curve(N, K, curve[0]);

// Use this version if N is known at compile time.
Hilbert<>::CurveN<N>(K, curve[0]);

// Use this version if K is known at compile time.
Hilbert<>::CurveK<K>(N, curve[0]);

// Use this version if N and K are known at compile time.
Hilbert<>::Curve<N, K>(curve[0]);
```

If `N` is known at compile time, the compiler can generate much faster
code. There's a much smaller marginal speedup if `K` is known at
compile time, and only when `N` is large. If binary size is more
important than runtime, you might want to stick to using the
non-templated versions.

### Using different integer types

The `Hilbert` class itself has two template parameters: `Int` and
`UInt` which are defaulted to `int` and `std::size_t`. If you know
that `K` is small, you can use a narrower integer like `int16_t` if `K
< 16` or `int8_t` if `K < 8`. This will also increase performance of
`Curve()` since it's mostly limited by memory.

It's recommended to leave `UInt` set to `std::size_t`. Some day I may
add support for bit vectors in case you want to compute individual
points of the curve for very large `N` or `K`, but it is mostly a
placeholder until then.

## Examples

See `example.cpp`.

## Building

The tests and example can be built with `cmake`:

```
cmake .
make
./hilbert_test
./example
```
