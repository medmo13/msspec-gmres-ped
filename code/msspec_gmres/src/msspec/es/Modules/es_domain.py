# coding: utf-8

import unittest
import numpy
import delaunay.core as delc
from ase import Atoms

def Define_Domain(set,ncover,ratio) :
    #defines the sphere englobing the molecule, such as we will have at least ncover empty spheres around the molecule,
    #and this empty sphere have a size around ratio x {global size of muffin tin atom size}


def Distance(A,B)
    #calculate distance between A and B
    try :
        assert len(A) == len(B)
        d=0
        for l in range(0,len(A))
            d+=(B[l]-A[l])**2
        d=sqrt(d)
    except AssertionError :
        print("Error : distance calculated from two points in different spaces")

