T = [matrix(ZZ, 4, 4, {(0,0):1, (1,2):1}), matrix(ZZ, 4, 4, {(2,0):1, (3,2):1}), matrix(ZZ, 4, 4, {(0,1):1, (1,3):1}), matrix(ZZ, 4, 4, {(2,1):1, (3,3):1})]

reps = [ [ [
    matrix(QQ, 4, 4, {(0,2):1, (1,3):1}), matrix(QQ, 4, 4, {(1,0):-1, (3,2):-1}), matrix(QQ, 4, 4, {})], [
    matrix(QQ, 4, 4, {(2,0):1, (3,1):1}), matrix(QQ, 4, 4, {(0,1):-1, (2,3):-1}), matrix(QQ, 4, 4, {})], [
    matrix(QQ, 4, 4, {(0,0):1, (1,1):1, (2,2):-1, (3,3):-1}), matrix(QQ, 4, 4, {(0,0):-1, (1,1):1, (2,2):-1, (3,3):1}), matrix(QQ, 4, 4, {})]], [
[
    matrix(QQ, 4, 4, {}), matrix(QQ, 4, 4, {(0,2):1, (1,3):1}), matrix(QQ, 4, 4, {(1,0):-1, (3,2):-1})], [
    matrix(QQ, 4, 4, {}), matrix(QQ, 4, 4, {(2,0):1, (3,1):1}), matrix(QQ, 4, 4, {(0,1):-1, (2,3):-1})], [
    matrix(QQ, 4, 4, {}), matrix(QQ, 4, 4, {(0,0):1, (1,1):1, (2,2):-1, (3,3):-1}), matrix(QQ, 4, 4, {(0,0):-1, (1,1):1, (2,2):-1, (3,3):1})]], [
[
    matrix(QQ, 4, 4, {(1,0):-1, (3,2):-1}), matrix(QQ, 4, 4, {}), matrix(QQ, 4, 4, {(0,2):1, (1,3):1})], [
    matrix(QQ, 4, 4, {(0,1):-1, (2,3):-1}), matrix(QQ, 4, 4, {}), matrix(QQ, 4, 4, {(2,0):1, (3,1):1})], [
    matrix(QQ, 4, 4, {(0,0):-1, (1,1):1, (2,2):-1, (3,3):1}), matrix(QQ, 4, 4, {}), matrix(QQ, 4, 4, {(0,0):1, (1,1):1, (2,2):-1, (3,3):-1})]]]

C = matrix(QQ,[[2, 0, 0], [0, 2, 0], [0, 0, 2]])

# x^1_2 y^1_1 (36,34) -> x^2_1 y^2_2
# x^1_1y^1_2-x^1_2y^2_2 (35,35) 
# x^2_2 y^1_2 (34,36) -> x^1_1 y^2_1

def new_from_old(reps):
    wtss = []
    repsn = []
    for xs,ys,hs in reps:
        repsn.append(xs)
        assert all(h.is_diagonal() for h in hs)
        wtss.append([vector(ZZ,h.diagonal()) for h in hs])
    return wtss, repsn

# vim: ft=python
 