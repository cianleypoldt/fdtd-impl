Maxwell Solver
==============

Finite-difference time-domain (FDTD) solver for Maxwell's equations in C.

Compile and run:
    gcc fdtd3d.c -o maxwell_solver -lm
    ./maxwell_solver

References:
-----------
[1] K. S. Yee, "Numerical solution of initial boundary value problems
    involving Maxwell's equations in isotropic media," IEEE TAP, vol. 14,
    no. 3, pp. 302-307, May 1966.
    DOI: 10.1109/TAP.1966.1138693
    URL: https://home.cc.umanitoba.ca/~lovetrij/cECE7810/Papers/Yee%201966%20HiRes.pdf

[2] W. C. Chew, "Lecture 37: Finite Difference Method, Yee Algorithm,"
    ECE 604, Purdue University, 2020.
    URL: https://engineering.purdue.edu/wcchew/ece604s20/Lecture%20Notes/Lect37.pdf

Notes:
------
Developed for educational purposes at age 17. Closely mirrors the referenced
methods with little regard for performance, extensibility, or advanced
features (e.g., absorbing boundaries, material variation).
