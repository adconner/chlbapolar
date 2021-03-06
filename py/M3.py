
# Add current directory to path to import bapolar
import os
import sys
sys.path.append(os.getcwd())

from sage.all import *
import bapolar

T = [
    matrix(QQ, 9, 9, {(0, 0): 1, (1, 3): 1, (2, 6): 1}),
    matrix(QQ, 9, 9, {(3, 0): 1, (4, 3): 1, (5, 6): 1}),
    matrix(QQ, 9, 9, {(6, 0): 1, (7, 3): 1, (8, 6): 1}),
    matrix(QQ, 9, 9, {(0, 1): 1, (1, 4): 1, (2, 7): 1}),
    matrix(QQ, 9, 9, {(3, 1): 1, (4, 4): 1, (5, 7): 1}),
    matrix(QQ, 9, 9, {(6, 1): 1, (7, 4): 1, (8, 7): 1}),
    matrix(QQ, 9, 9, {(0, 2): 1, (1, 5): 1, (2, 8): 1}),
    matrix(QQ, 9, 9, {(3, 2): 1, (4, 5): 1, (5, 8): 1}),
    matrix(QQ, 9, 9, {(6, 2): 1, (7, 5): 1, (8, 8): 1}),
]

reps = [
    (
        [
            vector(ZZ, [1, 0, -1, 0, 0, 0]),
            vector(ZZ, [1, 0, 1, -1, 0, 0]),
            vector(ZZ, [1, 0, 0, 1, 0, 0]),
            vector(ZZ, [-1, 1, -1, 0, 0, 0]),
            vector(ZZ, [-1, 1, 1, -1, 0, 0]),
            vector(ZZ, [-1, 1, 0, 1, 0, 0]),
            vector(ZZ, [0, -1, -1, 0, 0, 0]),
            vector(ZZ, [0, -1, 1, -1, 0, 0]),
            vector(ZZ, [0, -1, 0, 1, 0, 0]),
        ],
        [
            matrix(QQ, 9, 9, {(0, 3): 1, (1, 4): 1, (2, 5): 1}),
            matrix(QQ, 9, 9, {(3, 6): 1, (4, 7): 1, (5, 8): 1}),
            matrix(QQ, 9, 9, {(1, 0): -1, (4, 3): -1, (7, 6): -1}),
            matrix(QQ, 9, 9, {(2, 1): -1, (5, 4): -1, (8, 7): -1}),
            matrix(QQ, 9, 9, {}),
            matrix(QQ, 9, 9, {}),
        ],
    ),
    (
        [
            vector(ZZ, [0, 0, 1, 0, -1, 0]),
            vector(ZZ, [0, 0, 1, 0, 1, -1]),
            vector(ZZ, [0, 0, 1, 0, 0, 1]),
            vector(ZZ, [0, 0, -1, 1, -1, 0]),
            vector(ZZ, [0, 0, -1, 1, 1, -1]),
            vector(ZZ, [0, 0, -1, 1, 0, 1]),
            vector(ZZ, [0, 0, 0, -1, -1, 0]),
            vector(ZZ, [0, 0, 0, -1, 1, -1]),
            vector(ZZ, [0, 0, 0, -1, 0, 1]),
        ],
        [
            matrix(QQ, 9, 9, {}),
            matrix(QQ, 9, 9, {}),
            matrix(QQ, 9, 9, {(0, 3): 1, (1, 4): 1, (2, 5): 1}),
            matrix(QQ, 9, 9, {(3, 6): 1, (4, 7): 1, (5, 8): 1}),
            matrix(QQ, 9, 9, {(1, 0): -1, (4, 3): -1, (7, 6): -1}),
            matrix(QQ, 9, 9, {(2, 1): -1, (5, 4): -1, (8, 7): -1}),
        ],
    ),
    (
        [
            vector(ZZ, [-1, 0, 0, 0, 1, 0]),
            vector(ZZ, [1, -1, 0, 0, 1, 0]),
            vector(ZZ, [0, 1, 0, 0, 1, 0]),
            vector(ZZ, [-1, 0, 0, 0, -1, 1]),
            vector(ZZ, [1, -1, 0, 0, -1, 1]),
            vector(ZZ, [0, 1, 0, 0, -1, 1]),
            vector(ZZ, [-1, 0, 0, 0, 0, -1]),
            vector(ZZ, [1, -1, 0, 0, 0, -1]),
            vector(ZZ, [0, 1, 0, 0, 0, -1]),
        ],
        [
            matrix(QQ, 9, 9, {(1, 0): -1, (4, 3): -1, (7, 6): -1}),
            matrix(QQ, 9, 9, {(2, 1): -1, (5, 4): -1, (8, 7): -1}),
            matrix(QQ, 9, 9, {}),
            matrix(QQ, 9, 9, {}),
            matrix(QQ, 9, 9, {(0, 3): 1, (1, 4): 1, (2, 5): 1}),
            matrix(QQ, 9, 9, {(3, 6): 1, (4, 7): 1, (5, 8): 1}),
        ],
    ),
]

# for interactive use, turn verbose = True to show progress
res = bapolar.border_apolarity_cycl_inv(T, reps, 16, verbose = False)

print()
print(len(res[1]), "110 candidates")
print(len(res[0]), "111 candidates")
