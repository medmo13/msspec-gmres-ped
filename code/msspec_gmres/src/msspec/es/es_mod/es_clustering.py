# coding: utf-8
import unittest
from subprocess import call
import numpy as np
from numpy import sqrt, dot, cross
from numpy.linalg import norm
from ase import Atoms
from ase.io import read, write
from ase.visualize import view
import empty_spheres as esph
import es_tools as tool

# ===================================================================
# List of routines :
"""
=================
Tangent_Fourth_Sphere(Spheres_data, r, inout=1)         :   From 3 tangent spheres, returns the fourth, tangent to others, with radius r.
                                                            inout determines the side of the coordinate.
Spheres_Data_Structure_Extractor (Structure,list)       :   From Structure, returns data as [[S1][S2]...] where Si = [[xi,yi,zi],ri]
                                                            used in Tangent_Fourth routine. List determines the radius we will use.                                                      

Spheres_Data_XYZ_Extractor (name,list)                  :   From name, given as "_____.xyz", returns data as [[S1][S2]...] where Si = [[xi,yi,zi],ri]
                                                            used in Tangent_Fourth routine. List determines the radius we will use.    

Tetrahedron_ES (Spheres_data)                           :   From 4 spheres forming a tetrahedron,returns the sphere tangent to others, with radius r.

Triangle_ES (Spheres_data,R)                            :   Returns the 2 solutions of tangent of 3 spheres with radius R.                                                                                                         
=================
"""


# ===================================================================

def Tangent_Fourth_Sphere(Spheres_data, r, inout=1):
    # From Spheres_data build like : [[S1],[S2],[S3]] containing spheres informations Si=[[xi,yi,zi],ri], the wanted radius r,
    # and with the inout variable (equal to 1 or -1, depends on the facet, to build empty sphere in or out the hull)
    # computes the center of the fourth sphere tangent to S1, S2 and S3, and with radius r
    # WARNING : Require S1, S2 and S3 are tangent each to another !
    # ================
    # Read information
    S1, S2, S3 = Spheres_data

    d1 = S1[1] + r
    d2 = S2[1] + r
    d3 = S3[1] + r
    # Solve the problem in the right space : using u,v,t as base to solution
    u = tool.vector_def(S2[0], S1[0])
    v = tool.vector_def(S3[0], S1[0])
    unorm = tool.vector_norm(u)
    vnorm = tool.vector_norm(v)
    u = np.multiply(u, 1. / unorm)
    v = np.multiply(v, 1. / vnorm)
    w = np.multiply(S1[0], -2)
    # =====
    a = (d2 ** 2 - d1 ** 2 + S1[0][0] ** 2 - S2[0][0] ** 2 + S1[0][1] ** 2 - S2[0][1] ** 2 + S1[0][2] ** 2 -
         S2[0][2] ** 2) / (2 * unorm)
    b = (d3 ** 2 - d1 ** 2 + S1[0][0] ** 2 - S3[0][0] ** 2 + S1[0][1] ** 2 - S3[0][1] ** 2 + S1[0][2] ** 2 -
         S3[0][2] ** 2) / (2 * vnorm)
    c = d1 ** 2 - S1[0][0] ** 2 - S1[0][1] ** 2 - S1[0][2] ** 2
    # =====
    t = np.ndarray.tolist(np.cross(u, v))
    tnorm = tool.vector_norm(t)
    t = np.multiply(t, 1 / tnorm)
    #
    alpha = (a - b * np.dot(u, v)) / (1 - np.dot(u, v) ** 2)
    beta = (b - a * np.dot(u, v)) / (1 - np.dot(u, v) ** 2)
    d = alpha ** 2 + beta ** 2 + 2 * alpha * beta * np.dot(u, v) + alpha * np.dot(u, w) + beta * np.dot(v, w) - c
    theta = (-1 * np.dot(w, t) + inout * np.sqrt(np.dot(w, t) ** 2 - 4 * d)) / 2
    #
    solution = np.ndarray.tolist(np.multiply(u, alpha) + np.multiply(v, beta) + np.multiply(t, theta))
    return solution


# ===================================================================
def Spheres_Data_Structure_Extractor(Structure, list):
    # From Structure, returns data as [[S1][S2]...] where Si = [[xi,yi,zi],ri]
    # used in Tangent_Fourth routine. List determines the radius we will use.
    # list determines wich radius we take : 1 for covalent

    set_pos = np.ndarray.tolist(Structure.positions)
    set_nb = np.ndarray.tolist(Structure.numbers)
    data = []
    for i in range(0, len(set_nb)):
        data.append([set_pos[i], esph.Atom_Radius(Structure, i, list)])
    return data


# ===================================================================
def Spheres_Data_XYZ_Extractor(name, list):
    # From name, given as "_____.xyz", returns data as [[S1][S2]...] where Si = [[xi,yi,zi],ri]
    # used in Tangent_Fourth routine. List determines the radius we will use.
    # list determines wich radius we take : 1 for covalent
    set = read(name)
    set_pos = np.ndarray.tolist(set.positions)
    set_nb = np.ndarray.tolist(set.numbers)
    data = []
    for i in range(0, len(set_nb)):
        data.append([set_pos[i], esph.Atom_Radius(set, i, list)])
    return data


# ===================================================================
def Tetrahedron_ES(Spheres_data):
    # print "Spheres datas :",Spheres_data
    S1, S2, S3, S4 = Spheres_data
    x1, y1, z1 = S1[0]
    r1 = S1[1]
    x2, y2, z2 = S2[0]
    r2 = S2[1]
    x3, y3, z3 = S3[0]
    r3 = S3[1]
    x4, y4, z4 = S4[0]
    r4 = S4[1]

    # Solve this problem equals to found V=[x,y,z] and r as : MV = S + r P with :
    M = np.array([[x1 - x2, y1 - y2, z1 - z2], [x2 - x3, y2 - y3, z2 - z3], [x3 - x4, y3 - y4, z3 - z4]])
    P = np.array([r2 - r1, r3 - r2, r4 - r3])
    S = np.array([(r1 ** 2 - r2 ** 2) + (x2 ** 2 - x1 ** 2) + (y2 ** 2 - y1 ** 2) + (z2 ** 2 - z1 ** 2),
                  (r2 ** 2 - r3 ** 2) + (x3 ** 2 - x2 ** 2) + (y3 ** 2 - y2 ** 2) + (z3 ** 2 - z2 ** 2),
                  (r3 ** 2 - r4 ** 2) + (x4 ** 2 - x3 ** 2) + (y4 ** 2 - y3 ** 2) + (z4 ** 2 - z3 ** 2)])
    S = np.multiply(S, 0.5)

    """
    print("M= {}".format(M))
    print("P= {}".format(P))
    print("S= {}".format(S))
    #"""

    # MV = S + r P <=> V = Minv S + r MinvP, rewrite as V = Mbar + r Pbar
    if np.linalg.det(M) == 0:
        print("Singular matrix on the Tetrahedron Fourth SPhere Problem : So we cancel routine")
        return 999

    Minv = np.linalg.inv(M)
    Sbar = np.matmul(Minv, S)
    Pbar = np.matmul(Minv, P)
    """
    print("Minv= {}\n So verify inversion : Minv * M = \n{}".format(Minv,np.matmul(M,Minv)))
    print("Pbar= {}".format(Pbar))
    print("Sbar= {}".format(Sbar))
    #"""
    # V = [x,y,z] depends on r : We need to solve r first : as a root of a 2 degree polygon
    Vi = np.array(S1[0])
    ri = np.array(S1[1])
    D = Sbar - Vi
    a = np.linalg.norm(Pbar) ** 2 - 1
    b = 2 * (np.dot(D, P) - ri)
    c = np.linalg.norm(Sbar) ** 2 - 2 * np.dot(Vi, Sbar) - ri ** 2

    delta = b ** 2 - 4 * a * c  # Theorically delta >=0
    if delta < 0:
        return 999
    rsol = (-b - np.sqrt(delta)) / (2 * a)  # radius of sphere internally tangent
    # print("rsol found : {}".format(r))
    # print("Verify solution : ar² + br + c = {}".format(a * rsol**2 + b * rsol + c))

    # Now we can compute the solution V = x,y,z
    V = -1 * (Sbar + np.multiply(Pbar, rsol))
    r = tool.distance(V, Vi) - ri
    # print("r that should really be : {}".format(r))
    # print("Verify solution : ar² + br + c = {}".format(a * r ** 2 + b * r + c))
    return [np.ndarray.tolist(V), r]


# ===================================================================

def Triangle_ES(Spheres_data, R):
    # Routine posted by Andrew Wagner, Thanks to him.
    # Implementaton based on Wikipedia Trilateration article, about intercection of 3 spheres
    # This problem solves sphere tangent to 3 spheres : and returns the 2 solutions
    P1 = np.array(Spheres_data[0][0])
    r1 = Spheres_data[0][1] + R
    P2 = np.array(Spheres_data[1][0])
    r2 = Spheres_data[1][1] + R
    P3 = np.array(Spheres_data[2][0])
    r3 = Spheres_data[2][1] + R
    temp1 = P2 - P1
    if norm(temp1) == 0:
        return 666, 666
    e_x = temp1 / norm(temp1)
    temp2 = P3 - P1
    i = dot(e_x, temp2)
    temp3 = temp2 - i * e_x
    if norm(temp3) == 0:
        return 666, 666
    e_y = temp3 / norm(temp3)
    e_z = cross(e_x, e_y)
    d = norm(P2 - P1)
    j = dot(e_y, temp2)
    x = (r1 * r1 - r2 * r2 + d * d) / (2 * d)
    y = (r1 * r1 - r3 * r3 - 2 * i * x + i * i + j * j) / (2 * j)
    temp4 = r1 * r1 - x * x - y * y
    if temp4 < 0:
        # print("The three spheres do not intersect!")
        return 666, 666
    else:
        z = sqrt(temp4)
        p_12_a = P1 + x * e_x + y * e_y + z * e_z
        p_12_b = P1 + x * e_x + y * e_y - z * e_z

    p_12_a = np.ndarray.tolist(p_12_a)
    p_12_b = np.ndarray.tolist(p_12_b)

    return p_12_a, p_12_b

# ===================================================================



