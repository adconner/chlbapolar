from sage.all import *
import bapolar

T = [
    matrix(ZZ, 4, 4, {(0, 0): 1, (1, 2): 1}),
    matrix(ZZ, 4, 4, {(2, 0): 1, (3, 2): 1}),
    matrix(ZZ, 4, 4, {(0, 1): 1, (1, 3): 1}),
    matrix(ZZ, 4, 4, {(2, 1): 1, (3, 3): 1}),
]

reps = [
    (
        [
            vector(ZZ, [1, -1, 0]),
            vector(ZZ, [1, 1, 0]),
            vector(ZZ, [-1, -1, 0]),
            vector(ZZ, [-1, 1, 0]),
        ],
        [
            matrix(QQ, 4, 4, {(0, 2): 1, (1, 3): 1}),
            matrix(QQ, 4, 4, {(1, 0): -1, (3, 2): -1}),
            matrix(QQ, 4, 4, {}),
        ],
    ),
    (
        [
            vector(ZZ, [0, 1, -1]),
            vector(ZZ, [0, 1, 1]),
            vector(ZZ, [0, -1, -1]),
            vector(ZZ, [0, -1, 1]),
        ],
        [
            matrix(QQ, 4, 4, {}),
            matrix(QQ, 4, 4, {(0, 2): 1, (1, 3): 1}),
            matrix(QQ, 4, 4, {(1, 0): -1, (3, 2): -1}),
        ],
    ),
    (
        [
            vector(ZZ, [-1, 0, 1]),
            vector(ZZ, [1, 0, 1]),
            vector(ZZ, [-1, 0, -1]),
            vector(ZZ, [1, 0, -1]),
        ],
        [
            matrix(QQ, 4, 4, {(1, 0): -1, (3, 2): -1}),
            matrix(QQ, 4, 4, {}),
            matrix(QQ, 4, 4, {(0, 2): 1, (1, 3): 1}),
        ],
    ),
]

res = bapolar.border_apolarity_cycl_inv(T, reps, 6)

print()
print(len(res[1]), "110 candidates")
print(len(res[0]), "111 candidates")

# vim: ft=python
