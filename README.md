# Optimized N-dimensional Hilbert curve generation in C++

## Features

* *Speed*: Able to compute `IsToVs()` at `2.93 GiB/s` as indicated by
  running `benchmark` in `Release` mode with `clang++` on an i5-2500K
  system.
* *Small code size*: Currently the implementation is only about 200
  lines.  The functions are non-recursive and do not call other
  functions.  They should be trivially portable to most other
  languages.
* *Compile-time evaluation*: All functions are `constexpr`. See
  example usage in `example.cpp`.
* *No dependencies*. `hilbert.hpp` has 0 `include`s.
* Functions are also defined for `N = 0` and `K = 0`.
* *Full git history*: Everything has been left in, no matter how
  embarassing my mistakes may be.

## Installation

Simply add `hilbert.hpp` to your project.

## Usage

There are 4 functions for manipulating Hilbert curves:

1. `void IsToVs(ITy N, ITy K, ViTy curve[])`: Computes the `K`'th
   step of an `N` dimensional Hilbert curve.  The result will be
   stored in `curve`, which must have space for `2^(N*K)`
   `N`-dimensional vectors, for a total of `N*2^(N*K)` integers.
2. `void VsToIs(STy N, STy K, ITy is[])`: Computes the inverse of the
   `K`'th step of an `N` dimensional Hilbert curve.  The result will
   be stored in `is`, which must have space for `2^(N*K)` integers.
   `is` maps vectors to indices.  The vector is encoded in the index
   into `is`.  For example, if `N=3` and `K=2`, then `is` could be an
   array of type `ViTy[2][2][2]`, and to get the index for a vector
   `v` of type `ITy[3]`, you would use `is[v[2], v[1], v[0]]`.  If
   `is` is a flat array, you wold use `is[(v[0] << 0*2) + (v[1] <<
   1*2) + (v[2] << 2*2)]`.
2. `void IToV(ITy N, ITy K, ITy i, ViTy v[])`: Computes the `i`'th
   vector in `IsToVs(N, K)`.  The result will be stored in `v`.
   Requires `i < 2^(N*K)`.
3. `ITy VToI(ITy N, ITy K, ViTy v[])`: Computes the index that `v`
   would have in the `K`'th step of an `N` dimensional Hilbert curve.
   Behavior is undefined if `v` is not a point on the curve.  `v` will
   be zeroed when this function returns.

### The integer types `ITy`, `ViTy`, and `STy`

The `Hilbert` class itself has three tempalte parameters: `ITy`,
`ViTy`, and `STy` which are all defaulted to `unsigned int`.

`ViTy` is the type of the component of a vector.  If you know that `K`
is small, you can use a narrower integer for `ViTy` like `uint16_t` if
`K <= 16`. This will increase performance of `IsToVs()` since it's
mostly limited by memory.

`ITy` is the type of an index of a vector in the curve.  If you known
`N*K` is small, you can use a narrower integer for `ITy` like
`uint32_t` if `N*K <= 32`.  This will also increase performance of
`VsToIs()` since it's limited by memory.

`STy` is used as the type of the variables `N` and `K`.  There's not
many situations where you'd want to change it from the default.

### Calling the functions

The above functions live in a static class called `Hilbert`. To call
`IsToVs()`, for example, you might use:

```
// Compute the K'th iteration of an N-dimensional Hilbert curve.
unsigned int curve[1U << N * K][N];
Hilbert<>::IsToVs(N, K, curve[0]);
```

### Template versions of each function

Each function takes the number of dimensions `N` and the iteration `K`
as parameters.  If `N` and/or `K` are known at compile time, template
functions are available which can help the compiler generate more
optimized code. For example, the 4 versions of `IsToVs()` are listed
below.

```
// Compute the K'th iteration of an N-dimensional Hilbert curve.
constexpr int N = 2;
constexpr int K = 3;
unsigned int curve[1U << N * K][N];

// Untemplated version.
// Use this version if N and K are not known at compile time.
Hilbert<>::IsToVs(N, K, curve[0]);

// Use this version if N is known at compile time.
Hilbert<>::IsToVsN<N>(K, curve[0]);

// Use this version if K is known at compile time.
Hilbert<>::IsToVsK<K>(N, curve[0]);

// Use this version if N and K are known at compile time.
Hilbert<>::IsToVs<N, K>(curve[0]);
```

If `N` is known at compile time, the compiler can generate much faster
code. There's a much smaller marginal speedup if `K` is known at
compile time, and only when `N` is large. If binary size is more
important than runtime, you might want to stick to using the
non-templated versions.

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

A compiler with C++17 support is needed for `constexpr` usage.  If
`constexpr` is removed, the code should be portable to C++98.
