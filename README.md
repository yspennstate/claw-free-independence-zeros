# Two questions of Jerrum and Patel on zeros of the independence polynomial

This repository holds the paper `main.tex` (compiled as `main.pdf`) and the script `check.py`.

The paper is archived on Zenodo: https://doi.org/10.5281/zenodo.24434503 (all versions).

The paper answers two questions asked by Jerrum and Patel about the parabolic regions
R(K) = {z : (Im z)^2 < Re z / K + 1/(4K^2)}. Every claw-free graph of maximum degree at most Delta has the
multivariate zero-free region R(2500(Delta+1)), and every claw-free graph of clique number at most k has
R(2500 r(3,k)), where r(3,k) is the Ramsey number. Up to the constant, both are best possible. A family of
graphs of maximum degree 6 with no induced S_{1,2,2} has zeros at arguments pi/m + O(1/m^2), so no region
containing a sector |arg z| < theta is zero-free for all such graphs.

## check.py

Requires Python 3 and mpmath. Run `python check.py` (a few seconds). It checks:

1. for the graphs G_{m,h} of Section 4: the number of vertices and the maximum degree, the transfer-matrix
   formula against a direct evaluation of the independence polynomial, the absence of an induced S_{1,2,2}
   in G_{2,2} (with a positive control on S_{1,2,2} itself), and the zeros near 2^m e^{i pi/m} for
   m = 8, 12, 16, 20, 30;
2. on small random claw-free graphs with weights in R(2500(Delta+1)), that the independence polynomial is
   nonzero and the ratios Z_H/Z_B of the proof have real part greater than 0.97;
3. the examples of Proposition 1.2: K_{Delta+1} and the complements of C_5 and of the Wagner graph, with
   their independence polynomials computed from the graphs.
