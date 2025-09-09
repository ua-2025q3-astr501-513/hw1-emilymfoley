#!/usr/bin/env python3
#
# Please look for "TODO" in the comments, which indicate where you
# need to write your code.
#
# Part 3: Implement a Numerically Stable Quadratic Equation Solver (1 point)
#
# * Objective:
#   Implement a numerically stable quadratic equation solver that does
#   not catastrophic cancellation.
# * Details:
#   The description of the problem and the solution template can be
#   found in `hw1/p3.py`.
#
# From lecture `01w`, we learned about catastrophic cancellation---the
# significant loss of precision that occurs when subtracting two
# nearly equal numbers.
# This problem actually appeared in CK's research!
# While solving for the initial conditions of (unstable) spherical
# photon orbits around a black hole for the convergence test of
# [GRay2](https://ui.adsabs.harvard.edu/abs/2018ApJ...867...59C),
# catastrophic cancellation introduced errors so severe that photons
# would not remain on their spherical orbits for an radian.
# CK suspected a bug in the integrator and spent an entire month
# debugging the wrong part of the code.
# At the end, he realized the real problem.
# The standard quadratic formula we all learn in high school was
# simply not accurate enough for reliable numerical computation.
#
# Here, let's implement a numerically stable quadratic equation solver
# to overcome catastrophic cancellation.
# Please make sure that you take care of all the special cases.

import math

def quadratic(a, b, c):
    if a == 0:
        raise ValueError("Coefficient 'a' cannot be zero for a quadratic equation.")

    disc = b**2 - 4*a*c

    if disc < 0:
        return [None, None]     # no real roots
    elif disc == 0:
        x = -b / (2*a)
        return [x, x]           # repeated root
    else:
        sqrt_disc = math.sqrt(disc)
        if b >= 0:
            x1 = (-b - sqrt_disc) / (2*a)
        else:
            x1 = (-b + sqrt_disc) / (2*a)

        x2 = (c / a) / x1

        if x1 > x2:
            x1, x2 = x2, x1

        return [x1, x2]
