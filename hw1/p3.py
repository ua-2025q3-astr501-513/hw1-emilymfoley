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

import math

def quadratic(a, b, c):
    """Numerically stable quadratic equation solver.

    Solves a x^2 + b x + c = 0 avoiding catastrophic cancellation.

    Args:
        a, b, c: coefficients of the quadratic equation.

    Returns:
        (x1, x2): tuple of two roots.
                  If there are two real roots, x1 < x2.
                  If there is one real root, x2 is None.
                  If there are no real roots, x1 and x2 are None.
    """
    if a == 0:
        raise ValueError("Coefficient 'a' cannot be zero.")

    discriminant = b**2 - 4*a*c

    if discriminant < 0:
        # No real roots
        return None, None
    elif discriminant == 0:
        # One real root
        x = -b / (2*a)
        return x, None
    else:
        sqrt_disc = math.sqrt(discriminant)
        # Numerically stable computation
        if b >= 0:
            x1 = (-b - sqrt_disc) / (2*a)
        else:
            x1 = (-b + sqrt_disc) / (2*a)

        # Compute the other root using the conjugate trick
        x2 = (c / a) / x1

        # Ensure x1 < x2
        if x1 > x2:
            x1, x2 = x2, x1

        return x1, x2
