/**
 * Description:
 * Author:
 */

// for stuff that isn't long enough to include in it's own file but still important

/* FAST I/O */
ios_base::sync_with_stdio(false);
cin.tie(NULL);

/* MACROS */
#define int long long
#define pb push_back
#define mp make_pair
#define pi pair<int, int>
#define endl "\n"

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

/* STRESS TESTING (Windows) */
@echo off
set i=1

:loop
gen %i% > test.in
echo %i%
a < test.in > a.out
b < test.in > b.out
fc /b a.out b.out > nul
if errorlevel 1 (
goto :eof
)
set /a i += 1
goto loop
// for generator to take in command line arguments in C++, change main signature to int main(int argc, char* argv[])

/* STRESS TESTING (Linux) */
#!/bin/bash
for ((i = 1; ; i++))
do
    echo $i
    ./gen $i > test.in
    ./a < test.in > a.out
    ./b < test.in > b.out
    diff -w a.out b.out || break
done
