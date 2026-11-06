# coding: utf-8
from ase import Atoms
from ase.visualize import view
import numpy as np
import math

# List of tools :
"""
=================Vector tools==================
vector_def(A,B)     :   returns simply the vector translating A to B
vector_norm(V)      :   return euclidian norm of a vector
vector_trslt(P,V)   :   return image of translation of P from vector V
throw_away (P,O,d)  :   returns P' so as O,P and P' are aligned, and OP'=d. So it translate P from O with distance d to direction OP
angle_vector(u,v)   :   returns the value of the convex angle defined by vectors u and v, in radians.
ColinearTest(u,v)   :   returns 1 if u and v colinear, and 0 if not.

===========Distance and Proximity tools========
distance(a,b)                   :   calculate distance between 2 points
search_nearest(point,set,d)     :   search the point in set wich is the nearest from point. We must know that the min distance is d
point_set_proximity(point, set) :   returns the min distance between the point and the set of points
set_set_proximity(S1,S2)        :   returns minimal distance between each points of each sets.
longest_dist(set)               :   returns the longest distance between the points in the set
shortest_dist(set)              :   returns the shortest distance between the points in the set
dist_point_plan(Pt,Plan)        :   From Pt=[x,y,z] and Plan=[a,b,c,d], return distance beetween Pt and Plan

===============Construction tools===============
Isobarycenter(set)          :    Calculate isobarycenter of a set of points. Returns his coordinates    
Invert_Coord(set,O,r)       :    Apply circular inversion to every point in the set, excepted for origin, remaining origin.
Midpoint(P1,P2)             :    Computes the middle point of P1 and P2
rot3D(P,A,u)                :    Returns P' image of P by rotation from angle A around unity vector u
===================Data tools===================
commonlist(L1,L2)           :    Returns a list of elements common to L1 and L2, ie output is L1 intercection with L2
cleanlist(list)             :    Returns a list without repeated elements (each element in cleaned list appears only once)
"""


# ========================Vector Tools=============================
def vector_def(A, B):
    # returns simply the vector translating A to B
    l = len(A)
    if l != len(B):
        print("Major Error : The 2 given points aren't same dimensioned :\nA = {}\nB={}".format(A, B))
    V = []
    for i in range(0, l):
        V.append(B[i] - A[i])
    return V


# ==========================================
def vector_norm(V):
    # return euclidian norm of a vector
    l = (len(V))
    Origin = np.ndarray.tolist(np.zeros(l))
    return distance(V, Origin)


# ==========================================
def vector_trslt(P, V):
    x, y, z = P
    a, b, c = V
    return [x + a, y + b, z + c]


# ==========================================
def throw_away(P, O, d):
    # returns P' so as O,P and P' are aligned, and OP'=d. So it translate P from O with distance d to direction OP
    OP = tool.vector_def(O, P)
    OP = np.ndarray.tolist(np.multiply(OP, 1. / tool.vector_norm(OP)))
    output = np.ndarray.tolist(np.multiply(OP, d))
    return output


# ==========================================
def angle_vector(u, v):
    # returns the value in radian of the convex angle defined by vectors u and v.
    U = vector_norm(u)
    V = vector_norm(v)
    UV = np.dot(u, v)
    # print("U : {}\nV : {}\n u.v : {}".format(U,U,UV,))
    if UV == 0:
        # print("u.v = 0, on a donc un angle de {}, soit {}°".format(math.pi/2,math.pi/2*180/math.pi))
        return math.pi / 2
    else:
        # print("L'angle est donc donc arccos({}) = {} = {}°".format(UV / (U * V),np.arccos(UV/(U*V)),np.arccos(UV/(U*V))*180/math.pi))
        return np.arccos(UV / (U * V))


# ==========================================
def ColinearTest(u, v):
    # Tests if u and v colinear
    L = len(u)
    k = 999999999999  # initial value, to be sure first ktest will become k
    if L != len(v):
        print(
        "Error : u and v ant same dimension : \nu={}\nv={}\nSo we return u and v not aligned... but verify consequences".format(
            u, v))
    else:
        for i in range(0, L):
            if u[i] == 0 or v[i] == 0:
                if u[i] != 0 or v[i] != 0:
                    return 1
            else:
                ktest = u[i] / v[i]
                if k == 999999999999:
                    k = ktest
                elif np.abs(k - ktest) < 0.0000001:  # We accept almost colinearité at 1/10^6
                    return 1
    return 0


# ==========================================

# ==========================================
# ==========================================
# ==============Distance and Proximity Tools=======================
def distance(a, b):
    # Calculate distance between 2 points
    d = 0.0
    dim = len(a)
    for i in range(0, dim):
        d += (b[i] - a[i]) ** 2
    d = np.sqrt(d)
    return d


# ==========================================
def search_nearest(point, set, d):
    # search the point in set wich is the nearest from point. We must know that the min distance is d
    for p in set:
        if distance(point, p) == d:
            return p


# ==========================================
# ==========================================
def point_set_proximity(point, set):
    # returns the min distance between the point and the set of points
    d = distance(point, set[0])
    for p in set[1:]:
        d = min(d, distance(point, p))
    return d


# ==========================================
def set_set_proximity(S1, S2):
    # returns minimal distance between each points of each sets.
    d = point_set_proximity(S1[0], S2)
    for s in S1[1:]:
        d = min(d, point_set_proximity(s, S2))
    return d


# ==========================================
def longest_dist(set):
    # returns the longest distance between the points in the set
    L = len(set)
    if L < 2:
        print("Major ERROR : the set has become a 1-UPLET")
        quit()
    else:
        d = distance(set[L - 2], set[L - 1])  # The last ones
        for i in range(0, L - 2):
            for j in range(i + 1, L - 2):
                d = max(d, distance(set[i], set[j]))
    return d


# ==========================================
def shortest_dist(set):
    # returns the shortest distance between the points in the set
    L = len(set)
    if L < 2:
        print("Major ERROR : the set has become a 1-UPLET")
        quit()
    else:
        d = distance(set[L - 2], set[L - 1])  # The last ones
        for i in range(0, L - 2):
            for j in range(i + 1, L - 2):
                d = min(d, distance(set[i], set[j]))
    return d


# ==========================================
def dist_point_plan(Pt, Plan):
    # From Pt=[x,y,z] and Plan=[a,b,c,d] corresponding to plan's equation ax + by + cz + d = 0, give the distance beetween Pt and the Plan
    x, y, z = Pt
    a, b, c, d = Plan
    dist = np.abs(a * x + b * y + c * z + d) / np.sqrt(a ** 2 + b ** 2 + c ** 2)
    return dist


# ======================Construction tools=========================
def Isobarycenter(set):
    # Calculate isobarycenter of a set of points. Returns his coordinates
    x = y = z = 0.0
    l = len(set)
    for point in set:
        x += point[0] / l
        y += point[1] / l
        z += point[2] / l
    # print("New Centroid of cluter:",[x,y,z])
    return [x, y, z]


# ==========================================
def Invert_Coord(set, O, r):
    # Apply circular inversion to every point in the set, excepted for origin, remaining origin.
    output = []
    for point in set:
        if point == O:
            output.append(point)
        else:
            D = distance(O, point)
            OP = vector_def(O, point)
            OP = np.multiply(OP, 1. / D)  # set OP to unity vector
            add = np.ndarray.tolist(np.array(O) + np.multiply(OP, r ** 2 / D))
            output.append(add)
    return output


# ==========================================
def Midpoint(P1, P2):
    # From points coordinates P1 and P2 : construct the middle point of [P1,P2]
    output = np.ndarray.tolist(np.multiply((np.array(P1) + np.array(P2)), 0.5))
    return output


# ==========================================
def rot3D(P, A, u):
    # Verify first if u is unity vector :
    if vector_norm(u) != 1:
        u = np.ndarray.tolist(np.multiply(u, 1. / vector_norm(u)))
    # From P in R3, A in R and u in R3, returns the image from P's rotation around u from angle A
    x, y, z = P
    ux, uy, uz = u
    c = np.cos(A)
    s = np.sin(A)
    # We compute directly the result : see the 3D rotation matrix
    X = x * (ux ** 2 * (1 - c) + c) + y * (ux * uy * (1 - c) - uz * s) + z * (ux * uz * (1 - c) + uy * s)
    Y = x * (ux * uy * (1 - c) + uz * s) + y * (uy ** 2 * (1 - c) + c) + z * (uy * uz * (1 - c) - ux * s)
    Z = x * (ux * uz * (1 - c) - uy * s) + y * (uy * uz * (1 - c) + ux * s) + z * (uz ** 2 * (1 - c) + c)

    return [X, Y, Z]


# =========================Data Tools==============================
def commonlist(L1, L2):
    # Returns a list of elements common to L1 and L2, ie output is L1 intercection with L2
    output = []
    for i in L1:
        if i in L2:
            output.append(i)
    return output


# ==========================================
def cleanlist(list):
    # Returns a list without repeated elements (each element in cleaned list appears only once)
    ls = len(list)
    Index = range(ls)
    Output = []

    for iS in range(ls):
        if iS in Index:
            for iOS in range(iS + 1, ls):
                S = list[iS]
                OS = list[iOS]
                if S == OS:
                    Index.remove(iOS)  # S and OS are same coord or almost : we remove the last : OS
                    # else : iS correspond to coord needed to be deleted, so we don't add it

    for I in Index:
        Output.append(list[I])

    return Output

# ==========================================

