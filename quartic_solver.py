import mpmath as mp
from cubic_cardano_solver import solve_cubic

mp.mp.dps = 50 # precision

def solve_quartic(a, b, c, d, e, tol = mp.mpf(1e-25)):
    if abs(a) == 0:
        raise ValueError("Not a quartic")
    
    shift = b/(4*a)
    
    A = c/a - 6*(shift)**2
    B = d/a - 4*(shift)**3 - 2*(shift)*(c/a - 6*(shift)**2)
    C = e/a - (shift)**4 + (shift)**2*(c/a - 6*(shift)**2) - (shift)*(d/a - 4*(shift)**3)

    t_values = solve_cubic(8, -4*A, -8*C, 4*A*C-B**2)
    best_roots = None
    best_error = mp.inf

    for t in t_values:
        for R in [mp.sqrt(2*t - A), -mp.sqrt(2*t - A)]:
            if abs(R) < tol:
                roots = [
                    mp.sqrt(t+t**2-C)-shift,
                    -mp.sqrt(t+t**2-C)-shift,
                    mp.sqrt(t-t**2+C)-shift,
                    -mp.sqrt(t-t**2+C)-shift
                ]
            else:
                D1 = R**2 -4*(t-B/(2*R))
                D2 = R**2 -4*(t+B/(2*R))

                roots = [
                    (-R + mp.sqrt(D1))/2 - shift,
                    (-R - mp.sqrt(D1))/2 - shift,
                    ( R + mp.sqrt(D2))/2 - shift,
                    ( R - mp.sqrt(D2))/2 - shift,
                ]

            err = sum(abs(a*r**4 + b*r**3 + c*r**2 + d*r +e) for r in roots)

            if err < best_error:
                best_error = err
                best_roots = roots

    roots = best_roots

    f = lambda x: a*x**4 + b*x**3 + c*x**2 + d*x + e
    df = lambda x: 4*a*x**3 + 3*b*x**2 + 2*c*x + d
    for _ in range(2):
        for i in range(len(roots)):
            if abs(df(roots[i])) < tol:
                continue
            roots[i] = roots[i] - f(roots[i])/df(roots[i])
    
    def clean(z):
        real = 0 if abs(mp.re(z)) < tol else mp.re(z)
        imag = 0 if abs(mp.im(z)) < tol else mp.im(z)
        return mp.mpc(real, imag)

    return [clean(r) for r in sorted(roots, key=lambda z: (abs(z), mp.re(z), mp.im(z)))]

# -----------------------
# Example usage
# -----------------------

if __name__ == "__main__":
    print("Enter coefficients a b c d e (real or complex, e.g. 1+2j 0 -3 4j):")
    try:
        a, b, c, d, e = [mp.mpc(complex(x)) for x in input().split()]
    except:
        raise ValueError("Invalid input format. Use numbers like 1, -2.5, 3+4j")

    roots = solve_quartic(a, b, c, d, e)

    print("\nRoots:")
    for i, r in enumerate(roots, 1):
        print(f"r{i} = {mp.nstr(r, 10)}")

    print([mp.nstr(abs(a*r**4 + b*r**3 + c*r**2 + d*r +e), 10) for r in roots])