"""
Solve a general cubic equation using Cardano's method.

This implementation follows the classical analytical approach:
the cubic is first transformed into its depressed form, after which
Cardano's formula is applied to obtain the roots.

The solver supports complex coefficients and explicitly computes
cube roots in polar form to handle all branches correctly. A pairing
condition (αβ = A/3) is used to match corresponding cube roots,
ensuring the correct construction of all three solutions.

High-precision arithmetic (mpmath) is used to reduce numerical error,
and small residual components are cleaned for stable output.

This implementation focuses on mathematical transparency rather than
relying on built-in polynomial solvers.
"""

import mpmath as mp

mp.mp.dps = 50  # precision


def solve_cubic(a, b, c, d, tol=mp.mpf('1e-12')):
    if abs(a) == 0:
        raise ValueError("Not a cubic equation.")

    # Convert to mpmath complex
    a, b, c, d = mp.mpc(a), mp.mpc(b), mp.mpc(c), mp.mpc(d)
    scale = max(abs(a), abs(b), abs(c), abs(d))

    # Depressed cubic coefficients
    A = (1/3)*(b/a)**2 - c/a
    B = b*c/(3*a**2) - d/a - (2/27)*(b/a)**3

    # Discriminant
    D = B**2 - 4*A**3/27

    # Cube roots (polar form)
    alpha3 = (B + mp.sqrt(D)) / 2
    beta3 = (B - mp.sqrt(D)) / 2

    alpha_r, alpha_w = abs(alpha3), mp.arg(alpha3)
    beta_r, beta_w = abs(beta3), mp.arg(beta3)

    alpha = [
        mp.root(alpha_r, 3) * mp.e**(1j*(alpha_w + 2*k*mp.pi)/3)
        for k in range(3)
    ]

    beta = [
        mp.root(beta_r, 3) * mp.e**(1j*(beta_w - 2*k*mp.pi)/3)
        for k in range(3)
    ]

    # Match correct pairs
    idx = 0
    for i in range(3):
        if abs(alpha[0] * beta[i] - A/3) < tol*scale:
            idx = i
            break

    roots = [
        alpha[k] + beta[(k + idx) % 3] - b/(3*a)
        for k in range(3)
    ]

    f = lambda x: a*x**3 + b*x**2 + c*x + d
    df = lambda x: 3*a*x**2 + 2*b*x + c
    for _ in range(2):
        for i in range(len(roots)):
            if abs(df(roots[i])) < tol*scale:
                continue
            roots[i] = roots[i] - f(roots[i])/df(roots[i])

    # Clean numerical noise
    def clean(z):
        real = 0 if abs(mp.re(z)) < tol*scale else mp.re(z)
        imag = 0 if abs(mp.im(z)) < tol*scale else mp.im(z)
        return mp.mpc(real, imag)

    return [clean(r) for r in sorted(roots, key=lambda z: (abs(z), mp.re(z), mp.im(z)))]


# -----------------------
# Example usage
# -----------------------

if __name__ == "__main__":
    print("Enter coefficients a b c d (real or complex, e.g. 1+2j 0 -3 4j):")
    try:
        a, b, c, d = [mp.mpc(complex(x)) for x in input().split()]
    except:
        raise ValueError("Invalid input format. Use numbers like 1, -2.5, 3+4j")

    roots = solve_cubic(a, b, c, d)

    print("\nRoots:")
    for i, r in enumerate(roots, 1):
        print(f"r{i} = {mp.nstr(r, 10)}")