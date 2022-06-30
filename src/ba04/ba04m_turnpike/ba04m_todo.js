/*
Rosalind: BA4M
Solve the Turnpike Problem

If A = (a1 = 0, a2, ... , an) is a set of n points on a line segment
in increasing order (a1 < a2 < · · · < an),
then ΔA denotes the collection of all pairwise differences between points in A.
For example, if A = (0, 2, 4, 7), then

ΔA = (7, 5, 4, 3, 2, 2, 0, 0, 0, 0, 2, 2, 3, 4, 5, 7).
The following problem asks us to reconstruct A from ΔA.

Turnpike Problem
Given all pairwise distances between points on a line segment,
reconstruct the positions of those points.

Given: A collection of integers L.

Return: A set A such that ∆A = L.

Sample Dataset
-10 -8 -7 -6 -5 -4 -3 -3 -2 -2 0 0 0 0 0 2 2 3 3 4 5 6 7 8 10

Sample Output
0 2 4 7 10

═════════════════════════════════════════════════

Plan 1.

References:
-
*/
// #!/usr/bin/env javascript

function main() {
    const fs = require('fs');
    const lines = fs.readFileSync('/home/wsl/rosalind/data/ba04m.txt').toString().split("\n");

    const startTime = performance.now();
    console.log(`${performance.now() - startTime} milliseconds`);
}

// execute main function
main()