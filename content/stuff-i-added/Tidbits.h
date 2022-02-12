// for stuff that isn't long enough to include in it's own file but still important

/* PRAGMAS */
// #pragma GCC optimize("O3") - use for floating point cause Ofast uses fast-math which could mess up floating point precision
#pragma GCC optimize("Ofast")
#pragma GCC target("avx2")

/* RANDOM */
mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());
// uniform_int_distribution<int>(a, b)(rng) generates a random integer between a and b inclusive

/* PYTHON (if we need to avoid overflow) */
import sys
input, print = sys.stdin.readline, sys.stdout.write
// input will now contain a trailing newline character
// print will no longer automatically add a newline
