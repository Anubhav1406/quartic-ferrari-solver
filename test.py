from quartic_solver import solve_quartic
import mpmath as mp

mp.mp.dps = 50

if __name__ == "__main__":
    n = int(input())
    for _ in range(n):
        line = input()
        print(line)
        a, b, c, d, e = [mp.mpc(complex(x)) for x in line.split()]
        if abs(a) == 0:
            print("Not a quartic\n")
            continue
        roots = solve_quartic(a, b, c, d, e)
        for i, r in enumerate(roots, 1):
            print(f"r{i} = {mp.nstr(r, 10)}")
        print('\n')
