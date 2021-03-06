
# Add current directory to path to import bapolar
import os
import sys
sys.path.append(os.getcwd())

from sage.all import *
import bapolar

T = [
    matrix(QQ, 9, 9, {(4, 8): 1, (5, 7): -1, (7, 5): -1, (8, 4): 1}),
    matrix(QQ, 9, 9, {(3, 8): -1, (5, 6): 1, (6, 5): 1, (8, 3): -1}),
    matrix(QQ, 9, 9, {(3, 7): 1, (4, 6): -1, (6, 4): -1, (7, 3): 1}),
    matrix(QQ, 9, 9, {(1, 8): -1, (2, 7): 1, (7, 2): 1, (8, 1): -1}),
    matrix(QQ, 9, 9, {(0, 8): 1, (2, 6): -1, (6, 2): -1, (8, 0): 1}),
    matrix(QQ, 9, 9, {(0, 7): -1, (1, 6): 1, (6, 1): 1, (7, 0): -1}),
    matrix(QQ, 9, 9, {(1, 5): 1, (2, 4): -1, (4, 2): -1, (5, 1): 1}),
    matrix(QQ, 9, 9, {(0, 5): -1, (2, 3): 1, (3, 2): 1, (5, 0): -1}),
    matrix(QQ, 9, 9, {(0, 4): 1, (1, 3): -1, (3, 1): -1, (4, 0): 1}),
]

reps = [
    (
        [
            vector(ZZ, [1, 0, 1, 0]),
            vector(ZZ, [1, 0, -1, 1]),
            vector(ZZ, [1, 0, 0, -1]),
            vector(ZZ, [-1, 1, 1, 0]),
            vector(ZZ, [-1, 1, -1, 1]),
            vector(ZZ, [-1, 1, 0, -1]),
            vector(ZZ, [0, -1, 1, 0]),
            vector(ZZ, [0, -1, -1, 1]),
            vector(ZZ, [0, -1, 0, -1]),
        ],
        [
            matrix(QQ, 9, 9, {(0, 3): 1, (1, 4): 1, (2, 5): 1}),
            matrix(QQ, 9, 9, {(3, 6): 1, (4, 7): 1, (5, 8): 1}),
            matrix(QQ, 9, 9, {(0, 1): 1, (3, 4): 1, (6, 7): 1}),
            matrix(QQ, 9, 9, {(1, 2): 1, (4, 5): 1, (7, 8): 1}),
        ],
    ),
    (
        [
            vector(ZZ, [1, 0, 1, 0]),
            vector(ZZ, [1, 0, -1, 1]),
            vector(ZZ, [1, 0, 0, -1]),
            vector(ZZ, [-1, 1, 1, 0]),
            vector(ZZ, [-1, 1, -1, 1]),
            vector(ZZ, [-1, 1, 0, -1]),
            vector(ZZ, [0, -1, 1, 0]),
            vector(ZZ, [0, -1, -1, 1]),
            vector(ZZ, [0, -1, 0, -1]),
        ],
        [
            matrix(QQ, 9, 9, {(0, 3): 1, (1, 4): 1, (2, 5): 1}),
            matrix(QQ, 9, 9, {(3, 6): 1, (4, 7): 1, (5, 8): 1}),
            matrix(QQ, 9, 9, {(0, 1): 1, (3, 4): 1, (6, 7): 1}),
            matrix(QQ, 9, 9, {(1, 2): 1, (4, 5): 1, (7, 8): 1}),
        ],
    ),
    (
        [
            vector(ZZ, [1, 0, 1, 0]),
            vector(ZZ, [1, 0, -1, 1]),
            vector(ZZ, [1, 0, 0, -1]),
            vector(ZZ, [-1, 1, 1, 0]),
            vector(ZZ, [-1, 1, -1, 1]),
            vector(ZZ, [-1, 1, 0, -1]),
            vector(ZZ, [0, -1, 1, 0]),
            vector(ZZ, [0, -1, -1, 1]),
            vector(ZZ, [0, -1, 0, -1]),
        ],
        [
            matrix(QQ, 9, 9, {(0, 3): 1, (1, 4): 1, (2, 5): 1}),
            matrix(QQ, 9, 9, {(3, 6): 1, (4, 7): 1, (5, 8): 1}),
            matrix(QQ, 9, 9, {(0, 1): 1, (3, 4): 1, (6, 7): 1}),
            matrix(QQ, 9, 9, {(1, 2): 1, (4, 5): 1, (7, 8): 1}),
        ],
    ),
]

# for interactive use, turn verbose = True to show progress
res = bapolar.border_apolarity_cycl_inv(T, reps, 16, verbose = False)

print()
print(len(res[1]), "110 candidates")
print(len(res[0]), "111 candidates")
