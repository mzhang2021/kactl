# Math Reference Sheet

For semi-frequent stuff that I don't remember off the top of my head.

## BEST Theorem

Counts number of Eulerian circuits in directed graphs.

$$
\text{ec}(G) = t_w(G) \prod_{v \in V} (\deg(v) - 1)!
$$

$t_w(G)$ denotes the number of arborescences (directed trees pointing towards the root) rooted at any arbitrary node $w$ and is calculated via Kirchhoff's matrix tree theorem.

## Kirchhoff's Matrix Tree Theorem

Counts number of spanning trees of a graph. The algorithm builds the Laplacian matrix defined as the following for a simple graph:

$$
L_{i, j} := \begin{cases}
\deg(v_i) & \text{if } i = j \\
-1 & \text{if } i \neq j \text{ and } v_i \text{ is adjacent to } v_j \\
0 & \text{otherwise}
\end{cases}
$$

The answer is any cofactor of this matrix (e.g. the determinant after deleting the last row and column).

- Cayley's formula: the number of spanning trees of a complete graph of size $n$ is $n^{n-2}$
- For multigraphs, $L_{i, j}$ equals $-m$, where $m$ is the number of edges between $i$ and $j$, self-loops are excluded
- For directed multigraphs, $L_{i, j}$ equals $-m$, where $m$ is the number of edges from $i$ to $j$, and $L_{i, i}$ equals the indegree of $i$ minus the number of loops at $i$
- Removing the $i$th row and column and taking the determinant gives the number of oriented spanning trees rooted at (pointing towards) vertex $i$.

## Stirling Numbers of the First Kind

Counts number of permutations of length $n$ with $k$ cycles.

$$
dp[n+1][k] = n \cdot dp[n][k] + dp[n][k-1] \\
dp[0][0] = 1 \\
dp[n][0] = dp[0][k] = 0
$$

Explanation: if you have $n$ elements split into $k$ cycles and are inserting a new element, you can either create a new cycle or insert it directly behind any of the previous $n$ elements in an existing cycle.

The generating function for signed Stirling numbers of the first kind can be computed for a fixed $n$ in $\mathcal O(n \log^2 n)$:

$$
\sum_{k=0}^n s(n, k) x^k = x(x - 1) \dots (x - (n - 1))
$$

The unsigned Stirling numbers of the first kind are similar:

$$
\newcommand{\stirlingi}{\genfrac{[}{]}{0pt}{}}
\sum_{k=0}^n \stirlingi{n}{k} = x(x + 1) \dots (x + n - 1)
$$

The unsigned Stirling numbers of the first kind can be computed for a fixed $k$ in $\mathcal O(n \log n)$:

$$
\newcommand{\stirlingi}{\genfrac{[}{]}{0pt}{}}
\stirlingi{n}{k} = \frac{n!}{k!} [x^n] \left(\sum_{n=1}^\infin \frac{x^n}{n}\right)^k
$$

## Stirling Numbers of the Second Kind

Counts number of ways to partition $n$ labeled objects into $k$ non-empty unlabeled subsets.

$$
dp[n+1][k] = k \cdot dp[n][k] + dp[n][k-1] \\
dp[0][0] = 1 \\
dp[n][0] = dp[0][k] = 0
$$

Explanation: the $n + 1$th object is either a singleton or not. If it is a singleton, distribute the remaining $n$ objects among $k - 1$ groups, otherwise insert $n + 1$ into one of the $k$ groups and distribute the remaining $n$ objects in $k$ groups as well.

The generating function for Stirling numbers of the second kind can be computed for a fixed $n$ in $\mathcal O(n \log n)$:

$$
\sum_{k=0}^n S(n, k) x^k = \left(\sum_{i=0}^n \frac{(-1)^i}{i!}\right) \left(\sum_{i=0}^n \frac{i^n}{i!}\right) \mod x^{n+1}
$$

They can also be computed for a fixed $k$ in $\mathcal O(n \log n)$:

$$
S(n, k) = n! [x^n] \frac{(e^x - 1)^k}{k!} \\
e^x = \sum_{n=0}^\infin \frac{x^n}{n!}
$$

## Partition Function

Counts number of ways to partition $n$ into non-negative integer parts. Partition of $n$ into $k$ parts follows the following recurrence:

$$
dp[n][k] = dp[n-1][k-1] + dp[n-k][k] \\
dp[0][0] = 1 \\
dp[n][0] = dp[0][k] = 0
$$

Explanation: there are two possibilities. If we include a $1$ in the partition, we simply partition the remaining of $n - 1$ into $k - 1$ parts. Otherwise, each part has size greater than $1$, so we subtract $1$ from each part and solve recursively.

$p(n)$ denotes the number of partitions of $n$ and its generating function can be computed in $\mathcal O(n \log n)$:
$$
\sum_{n=0}^\infin p(n) x^n = \prod_{k=1}^\infin \frac{1}{1 - x^k}
$$

The denominator can be computed in $\mathcal O(n)$ time with the pentagonal number theorem:

$$
\prod_{n=1}^\infin (1 - x^n) = 1 + \sum_{k=1}^\infin (-1)^k \left(x^{k(3k+1)/2} + x^{k(3k-1)/2}\right)
$$

## Derangements

Counts the number of permutations where $p_i \neq i$ for all $i$.

$$
dp[n] = (n - 1)(dp[n-1] + dp[n-2]) \\
dp[0] = 1 \\
dp[1] = 0
$$

Explanation: With $n$ elements, consider index $1$ which may receive $n - 1$ different values. There are two cases. Either index $1$ swaps with another index $i$, so we count the number of derangements among the remaining $n - 2$ indices. Alternatively, index $1$ receives value $i$ but index $i$ does not receive value $1$, so this is equivalent to counting the number of derangements of $n - 1$ indices as we can renumber value $1$ as value $i$.

The probability of getting a derangement from a random shuffle is $\frac{1}{e}$, which means the number of shuffles needed to get a derangement is effectively constant.

## Catalan Numbers

Shows up in numerous counting problems, such as number of correct bracket sequences of length $2n$, number of full binary trees with $n + 1$ leaves, etc.

$$
\{1, 1, 2, 5, 14, 42, 132, 429, \dots\} \\
C_n = \frac{1}{n + 1} \binom{2n}{n} \\
C_{n+1} = \sum_{i=0}^n C_i C_{n-i} \\
C_0 = 1
$$

Proof for combinatorial definition of Catalan numbers is based on reflection argument for grid paths.