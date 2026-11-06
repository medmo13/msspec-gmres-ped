# coding: utf-8

import unittest
from subprocess import call
import numpy as np
from ase import Atoms
from ase.io import read,write
import empty_spheres as esph
import es_tools as tool

#===================================================================

def Apollonius_CCC(circles_data) :
    # From the circle_data list defined like : [[list of centers coordinates 3D],[list of radius]]
    # defines the Apollonius point, equidistant to the 3 circles
    #___________________________________________________________
    # Regroup data, and define the case
    O1=circles_data[0][0]
    O2=circles_data[0][1]
    O3=circles_data[0][2]
    r1 = circles_data[1][0]
    r2 = circles_data[1][1]
    r3 = circles_data[1][2]
    rmin=min(circles_data[1])
    nbmin=circles_data[1].count(rmin)#So it is the number of circles with the smallest radius
    if nbmin==1 : #then we have 1 little circle, we go to Apollonius CCP
        pass
    elif nbmin==2: #then we have 2 little circles, we go to Apollonius CPP
        if r1!=rmin :
            data=[[O1,r1],O2,O3]
        elif r2!=rmin :
            data = [[O2, r2], O1, O3]
        elif r3!=rmin :
            data = [[O3, r3], O1, O2]
        Apollo=Apollonius_CPP(data)

    elif nbmin==3 :#then the 3 circles have the same radius, so we search simply the centroid, form Apollonius PPP
        data=[O1,O2,O3]
        Apollo=tool.Isobarycenter(data)

    return Apollo

# ===================================================================
def Apollonius_CPP(set) :
    #From a set define like this : [[Coordinates of circle center,Circle radius],Coord P1, Coord P2]
    # define the point equidistant to P1, P2 and the circle. We will use circular inversion method
    Apollo=0
    return Apollo
# ===================================================================
def Circular_Inversion (P,O,R) :
    #Computes P' the image of circular inversion of P by the circle (O,R). It means OP*OP'=R²
    OP=tool.vector_def(O,P)
    nOP=tool.vector_norm(OP)
    OPprim = (OP * R ** 2) / (nOP ** 2)
    Apollo = np.ndarray.tolist(np.array(OPprim)-np.array(O))
    return Apollo
# ===================================================================




