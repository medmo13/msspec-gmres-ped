# coding: utf-8

import unittest
import delaunay.core as delc
from subprocess import call
import numpy as np
from scipy.spatial import ConvexHull
from ase import Atoms
from ase.io import read,write
from ase.visualize import view
import empty_spheres as esph
import es_tools as tool
import math
#===================================================================
# List of routines :
"""
=================
sym_analyse(Cluster)                        :    Convert set into xyz folder, then finds all symetries using. Uses compare_sym and read_sym_file
major_plane(set,[multiple])                 :    Search the major(s) plane with max nb of points in set and returns all of his points, and his equation
                                                 if multiple = True, returns all Plane equations of planes containing enough points
Cluster_flatness_informations(set)          :    Returns set's hull's total volume and area
Cluster_search_hollow(set,tol)              :    Returns the hollow datas : hollow=[[list,center,volume]...] where list is hollow' vertice's 
Cluster_emptyness_informations(Structure)   :    Returns the % of volume of hull occupied by his spheres
                                                 index,[center,volume] his hollow's center and volume. tol defines hollow's diagonale
Vertice_Sphere_Proportion(O,hull)           :    Returns the proportion of the sphere centered in Oth pt in the hull (proportion in ]O,1])
hull_search_neighbors(O,Simplices)          :    Returns list including O and neighbors (list contains only index, no coordinates)
convex_base(Neil,Simplices)                 :    Returns all Neil[0] neighbors so as ConvBase is the base of pyramid from top Neil[0]
Neighbors_List(numAtom,Structure)           :    Returns index list in Structure of neighbors of Atom indexed numAtom. numatom not included
Hull_Tetracut(O,Structure,hull)             :    Returns index list in Structure of centers of spheres defining the terahedron cuting the vertice O
facet_cap(O,Structure,hull)                 :    Return the proportion of sphere centered in O out from the hull 
                                                 (ie return the cap proportion of the sphere defined by cuting sphere with hull facets)
=================
"""
#===================================================================
def sym_analyse(Cluster) :
    #This function convert the given set into a xyz folder, then uses modified clus_geom executable to find
    # the set symmetries. Then, it computes and returns the list of symmetries the list into "sym"
    #Sym begins with the integer number of symetries, then with all symmetry-names on strings
    sym=[]
    write('Cluster.xyz',Cluster) #Sym_Analys folder contains actually the cluster_geom executables
    call(["./es_mod/Sym_Analys/proc_geom"])
    sym=read_sym_file('sym_analysis.lis')
    #For the moment, the file isn't removed...
    return sym
# ===================================================================
def compare_sym(Nsym,Osym):
    # Compare the new list of symmetries Nsym to the old one (Osym).
    # If the lists are not the same, it tells user to fusion more empty-spheres, to keep sym unchanged
    # Situation will be defined with output : 1 if symmetry conserved, 0 if not.
    # If there was no symmetry except identity : (-1). And if symmetries have been gained : 2.
    N=Nsym[0]
    O=Osym[0]
    print("Old number of sym :",O)
    print("New number of sym :", N)
    if N>O :
        output=2
        print("Symmetries gained. Symmetry-list updated")
    elif O==1 :
        output=-1
        print("No symmetry existing to help at decision")
    elif N == O:
        output=1
        print("Symmetry has been conserved.")
    elif N<O :
        output=0
        print("Symmetries Lost: Fusion-Operation on last empty-spheres advised")
    return output
# ===================================================================
def read_sym_file(symfile) :
    #The function reads the .lis file. The first argument is a 2space-int, after this,
    # we have 5 space-string. Each time, we have to delete the \n readed
    sym = []
    file=open(symfile,"r")
    NbSym=int(file.read(2))
    sym.append(NbSym)
    delete=file.read(1)#to delete the \n
    for i in range(0,NbSym) :
        Namesym=str(file.read(5))#read the 5 sym-name char, + the 2 empty spaces as char
        delete=file.read(1)#to delete the \n
        sym.append(Namesym)
    file.close()
    return sym
# ===================================================================
def major_plane(set,multiple=False) :
    #This function search the major plane where we can find the maximum of set-points.
    #In case of a draw, it takes the plane where the points are the closest (to be conform to molecules planes)
    l=len(set)
    PlanesLen=[]
    #We will use this lists :
    AllPlanes=[]
    Plane_list=[]#List of all planes, defined by every points in this plane
    Plane_eq=[] #List of 4-uplet corresponding to plane-equation for each Plane. Linked to Plane_list
    #
    #Regroup every unordered combinaisons of 3-uplet from set :
    Pts_Combi=combinaison_3_noorder_norepeat(set)
    #print("Pts_Combi={}".format(Pts_Combi))

    #For every Plane in Pts_Combi, define plane equation from the 3 pts defining P
    for P in Pts_Combi :
        #print("3 points Combi studied : {}".format(P))
        #Check caracteristics of the plane P, then update Plane_list and Plane_eq
        v1=tool.vector_def(P[0],P[1])
        n1=tool.vector_norm(v1)
        if n1 == 0 :
            print("n1 = 0 on this Pt Combi :{}\nVerify no double coordinates given".format(P))
        u1=np.ndarray.tolist(np.multiply(v1,1./n1))
        v2 = tool.vector_def(P[0],P[2])
        n2 = tool.vector_norm(v2)
        if n2 == 0 :
            print("n2 = 0 on this Pt Combi :{}\nVerify no double coordinates given".format(P))
        u2 = np.ndarray.tolist(np.multiply(v2, 1. / n2))
        #print("The 2vectors for the potential plane : \nv1={}, of norm{}, normed to {}\nv2={},of norm{}, normed to {}".format(v1,n1,u1,v2,n2,u2))

        if tool.ColinearTest(u1,u2) != 1 :#IE the 3 points are not aligned
            norm=np.cross(v1,v2) #normal vector to the plane
            #print("u1 :{}\nu2 :{}\n       So norm is {} ".format(u1,u2,norm))
            a, b, c = norm
            d=np.dot(norm,P[2])
            Equation=[a,b,c,d]

            """print("Equation of plane found : {}\nFrom Combi {}\n".format(Equation,P))
            if abs(a)+abs(b)+abs(c)==0 :
                print("More Details :\n\nv1 = {}\nv2={}\nSo norm ={}\n\n".format(v1,v2,norm))
            #"""

            # Add points on Plane list

            if Equation not in Plane_eq :
                Plane_pts=[]
                Plane_eq.append(Equation)
                for pt in set :
                    x, y, z = pt
                    #print("Pt : {}, Eq : {}\nSo ax +by + cz +d = {}".format(pt, Equation, a * x + b * y + c * z + d))
                    if abs(a * x + b * y + c * z + d) < 0.02:
                        Plane_pts.append(pt)
                Plane_list.append(Plane_pts)
                PlanesLen.append(len(Plane_pts))

            #Else : Case already found

    #Find now the plane with the most points in :

    #print("Plane_list = {}".format(Plane_list))


    if multiple == False :

        Plane = Plane_list[0]

        for P in Plane_list[1:]:
            if len(P) > len(Plane):
                Plane = P
            elif len(P) == len(Plane):
                if tool.longest_dist(P) < tool.longest_dist(Plane):
                    Plane = P
            # Else, Plane remains the one with the most and the nearest points
        #
        index = Plane_list.index(Plane)
        Plane_eq = Plane_eq[index]
        return [Plane, Plane_eq]

    elif multiple == True :
        #print("We have Planelist at {} element and Planeeq at {} elements".format(len(Plane_list),len(Plane_eq)))
        CleanLen=tool.cleanlist(PlanesLen)
        print("Planes include number of spheres in {}".format(sorted(CleanLen)))
        planesize = input("Please select minimal size of plane we should take :     ")
        AllPlanes = []

        for P in Plane_list :
            if len(P) >= planesize :
                index=Plane_list.index(P)
                PEQ=Plane_eq[index]

                #AllPlanes.append(PEQ)
                AllPlanes.append([P,PEQ])
        """
        print("In total we have {} planes".format(len(AllPlanes)))
        for P in AllPlanes:
            print("P n°{} : Equation : {}\nPoints : {}\n".format(AllPlanes.index(P) + 1, P[1],P[0]))

        #"""

        #Clean Double Planes : ax + by + cz + d = 0 <=> -ax -by -cz -d =0

        AllPlanesIndex = range(len(AllPlanes))
        Output = []

        for iP in range(len(AllPlanes)):
            if iP in AllPlanesIndex:
                for iOP in range(iP + 1, len(AllPlanes)):
                    Plist = AllPlanes[iP][0]
                    OPlist = AllPlanes[iOP][0]
                    if Plist == OPlist:
                        AllPlanesIndex.remove(iOP)  # P and OP are same Planes as they include same points
        for I in AllPlanesIndex:
            Output.append(AllPlanes[I])

        return Output


# ===================================================================

def combinaison_3_noorder_norepeat(set) :
    #returns every 3-uplet combinaison of the set, with no repetition, and no order variation (ie :{1,2,3}={1,3,2})
    l=len(set)
    Combi = []
    for i in range(0, l - 2):
        for j in range(i + 1, l - 1):
            for k in range(j + 1, l):
                Combi.append([set[i], set[j], set[k]])
    return Combi
# ===================================================================
def Cluster_flatness_informations(set) :
    #returns the total volume and area of the set's hull, then
    hull= ConvexHull(set)
    set_area=hull.area
    set_volume=hull.volume
    return[set_volume,set_area]
# ===================================================================
def Cluster_search_hollow(set,tol) :
    #Returns the hollow datas : hollow=[[[list],center,volume]...] where list is hollow vertices, and center and volume respectively
    #his center and volume.'Tol' can be 2 (big) or sqrt(3) to consider cubes as hollow, or even sqrt(2) to consider terahedrons as hollow
    #Notice that no hollows are detected fo a pentagone, due to his particular alignment of vertices...
    #Compute Centroid of Cluster :
    dmin = tool.shortest_dist(set)
    #Search for hollows
    hollow = []
    L=len(set)
    for i in range(0,L-1):
        for j in range(i+1,L) :
            #Data print (debug) :
            P1=set[i]
            P2=set[j]
            print("i,j : {},{}\nP1 :{}\nP2 :{}\ndistance(P1,P2):{}\n dmin={} so tol becomes {}\n".format(i+1,j+1,P1,P2,tool.distance(P1,P2),dmin,tol*0.99*dmin))
            if tool.distance(P1,P2) > 0.99*tol * dmin :#Then we may have find a hole center between P1 and P2
                hcenter=tool.Midpoint(P1,P2)
                d = tool.distance(P1,P2)/2
                print("hcenter is distant form set to {}...\nMust be > {}".format(tool.point_set_proximity(hcenter, set),0.98*d))
                if tool.point_set_proximity(hcenter, set) > 0.98*d :#We have a hollow there :
                    hollow_list=[]
                    hollowinvset=tool.Invert_Coord(set,hcenter,10)#10 is just decent : the volumewill be multiplied by 1000
                    hollowhull=ConvexHull(hollowinvset)
                    hollow_index =np.ndarray.tolist(hollowhull.vertices)
                    hollow_volume=hollowhull.volume/1000
                    hollow.append([hollow_index,hcenter,hollow_volume])
                    print("Hollow founded :{}\n".format(hcenter))
                    #
                #else : no hollow finally
            # else : P1 and P2 pretty neighboors
        #
    #
    return hollow
# ===================================================================
def Cluster_emptyness_informations(Structure) :
    #returns the % of volume of hull occupied by his spheres
    Allproportions=[]
    set=Structure.positions
    hull= ConvexHull(set)
    #esph.lookhull(Structure,hull)
    set_volume=hull.volume
    spheres_volume=0
    for pt in range(0,len(set)) : #Study all set point from index.
        print("Pt n°{}".format(pt))
        if pt in hull.vertices :
            #pt is a boundary point. We need to know the proportion of the sphere in the hull.
            proportion=Vertice_Sphere_Proportion(pt,hull,Structure)
        else : #We must verify sphere is completely is the cluster or not : it can have a cap out :
            NeibPt=Neighbors_List(pt,Structure)
            NeibInHull=tool.commonlist(NeibPt,hull.vertices)
            if len(NeibInHull) != 0 :#We must search the cap cutting the sphere.
                cap=facet_cap(pt,Structure,hull)
                proportion = 1 - cap
            else :
                proportion = 1

        previous = spheres_volume
        spheres_volume += proportion * (4 * math.pi / 3) * (esph.Atom_Radius(Structure,pt,1)**3)
        print("So sphere n°{} is at proportion {} in hull : We add {} to total sphere volume, being at {} now".format(pt,proportion,spheres_volume-previous,spheres_volume))
        raw_input("\nPress Enter to continue ...\n")

        Allproportions.append(proportion)
    print("AllProportions : ",Allproportions)

    print("Hull volume = {}".format(set_volume))
    return spheres_volume / set_volume * 100
# ===================================================================

# ===================================================================
def Vertice_Sphere_Proportion(O,hull,Structure) :
    #Returns the proportion of the sphere centered in O in the hull (proportion in ]O,1]).
    #O is the index of the center in hull vertices list.
    R = esph.Atom_Radius(Structure, O, 1) # Radius of O
    Vertices_list = np.ndarray.tolist(hull.vertices)
    Simplices_list = np.ndarray.tolist(hull.simplices)
    Point_list = np.ndarray.tolist(hull.points)
    Norm_list = np.ndarray.tolist(hull.equations)
    #print("Vertices : {}\nSimplices :{}\nPoints : {}\nO :{}".format(Vertices_list,Simplices_list,Point_list,O))
    Proportion = 0    #initialisation
    Neighbors=Hull_Neighbors_List(O,Structure,hull)
    print("Neighbors :{}".format(Neighbors))
    #Neighbors=tool.commonlist(Neighbors,Vertices_list)#We delete all neighbors not included in the hull
    Neighbors.insert(0,O)
    #print("Neighbors in hull :{}".format(Neighbors))

    if len(Neighbors) == 3 : #O is on a ridge
        Pt=[]
        for pt in Neighbors :
            Pt.append(hull.points[pt])
        Ox = tool.vector_def(Pt[0],Pt[1])
        Oy = tool.vector_def(Pt[0], Pt[2])
        Alpha=tool.angle_vector(Ox,Oy)
        Proportion=Alpha / (2 * math.pi)

    else :

        if len(Neighbors) < 3 : #There is an error here
            print("Hull is strangely defined : you can check on view")
        #Test if O is in a big facet :
        norms=[]
        for Sim in Simplices_list :
            if O in Sim :
                indx=Simplices_list.index(Sim)
                norm=Norm_list[indx][:3]
                if norm not in norms :
                    norms.append(norm)
        if len(norms)==1 : #Only one norm for all facets containing O : O in center of a big facet
            Proportion = 0.5
        else : #Here the work begins : O top of a pyramid

            Triag = delc.Triangulation(Structure.positions)
            for tetra in Triag.indices :
                if O in tetra: #A tetrahedron containing O : O will be counted as top
                    print("Tetraedron found with {} in : {}".format(O,tetra))
                    tetra.remove(O)
                    P1 = Point_list[tetra[0]]
                    P2 = Point_list[tetra[1]]
                    P3 = Point_list[tetra[2]]
                    #View tetrahedron :
                    Tetraview=Atoms("XH3",positions=[Point_list[O],P1,P2,P3])
                    view(Structure+Tetraview)
                    V1 = tool.vector_def(Point_list[O], P1)
                    V2 = tool.vector_def(Point_list[O], P2)
                    V3 = tool.vector_def(Point_list[O], P3)
                    print("Vectors from top of pyramid to base points" , [V1, V2, V3])
                    a = tool.angle_vector(V1, V2)
                    b = tool.angle_vector(V2, V3)
                    c = tool.angle_vector(V3, V1)

                    # """____________________________________________________________
                    # L'Huilier Theorem :
                    p = (a + b + c) / 2
                    #print("p = {} , or {}° , Angles : {}, {}, and  {}".format(p,p*180/math.pi,a*180/math.pi,b*180/math.pi,c*180/math.pi))
                    S = 4 * np.arctan(
                        np.sqrt(np.tan(p / 2) * np.tan((p - a) / 2) * np.tan((p - b) / 2) * np.tan((p - c) / 2)))
                    #print("L'huilier pocess :\nS   =4arctan( Sqrt( tan(p/2)*tan(p-a/2)* tan(p-b/2) * tan(p-c/2) ) )\n    =4arctan(sqrt(tan({})*tan({})*tan({})*tan({})))\n    =4arctan(sqrt({}*{}*{}*{}))\n    =4arctan({})\n    = {} \n\n Sphere Volume = {}, so proportion = {}".format(p/2,(p-a)/2,(p-b)/2,(p-c)/2,np.tan(p/2),np.tan((p-a)/2),np.tan((p-b)/2),np.tan((p-c)/2),np.tan(p/2)*np.tan((p-a)/2)*np.tan((p-b)/2)*np.tan((p-c)/2),S,(4 * math.pi),S / (4 * math.pi)))
                    print("S with Huilier = {}".format(S))
                    # """____________________________________________________________

                    # """____________________________________________________________
                    # Spherical area computing : Sinus formula + Girard's formula
                    cosT = (np.cos(c) - np.cos(a)*np.cos(b)) / (np.sin(a)*np.sin(b))
                    Theta = np.arccos(cosT)
                    sinT = np.sin(Theta)
                    sinA = np.sin(a) * sinT / np.sin(c)
                    sinB = np.sin(b) * sinT / np.sin(c)
                    print("sinB detail : b = {}, sin(b) = {}, sinT = {},c = {} sin(c) = {}\So : sinB = {}".format(b,np.sin(b),sinT,c,np.sin(c),sinB))
                    if sinA >=1 :
                        Alpha=math.pi/2
                    else :
                        Alpha = np.arcsin(sinA)
                    if sinB >=1 :
                        Beta=math.pi/2
                    else :
                        Beta = np.arcsin(sinB)
                    S = (Alpha + Beta + Theta - math.pi)
                    print("S with sinus and Girard's formulaes = {}".format(S))
                    #"""____________________________________________________________



                    Proportion += S / (4 * math.pi)
                    print("We had {} to actual {} Proportion, being now at {}".format(S / (4*math.pi), Proportion - S / (4 * math.pi), Proportion))
                #else the tetraedron in disconnected to the center of sphere : no need to study him

    print("proportion : {}".format(Proportion, 1. / Proportion))

    """Little thing to see exactly what happens
    print("proportion : {} , or 1/{}".format(Proportion, 1./Proportion))
    Lset=len(Neighbors)-1
    set = []
    Lhul=len(np.ndarray.tolist(hull.points))
    namehull= "C" +str(Lhul)
    HullView=Atoms(namehull,positions=np.ndarray.tolist(hull.points))
    for pt in Neighbors:
        set.append(Point_list[pt])
    set.remove(Point_list[O])
    nameset= "H" + str(Lset)
    ViewNeig=Atoms(nameset,positions=set)
    Ostr=Atoms(positions=[Point_list[O]])
    view(HullView+ViewNeig+Ostr)
    #"""

    return Proportion

# ===================================================================

def hull_search_neighbors(O,Simplices) :
    #Returns the list including O in fist position and all of his neighbors (list contains only index, no coordinates)
    Neil=[O]
    print("O in the search neighbors routine :",O)
    for facet in Simplices :
        if O in facet :
            for i in facet :
                if i not in Neil :
                    Neil.append(i)
                #else : i still in neighbors list
    return Neil
# ===================================================================

def convex_base(Neil,Simplices) :
    #Order the neighbors list by listing them so as  each consecutives elements in the list are themselves neighbors.
    #So it returns all Neil[0] neighbors so as ConvBase is the base of pyramid from top Neil[0]
    ConvBase=[]
    ConvBase.append(Neil[1])
    Rest=Neil[2:]
    print("Initial Neil : ",Neil)
    while len(Rest)>0 :
        for i in Rest :
            print("Convbase :{}\nRest :{}\nactual i:{}".format(ConvBase, Rest,i))
            Last=ConvBase[-1]
            print(Last)
            Neighborsi=hull_search_neighbors(i,Simplices)
            print("The last elements :{}\nHis neighbors : {}".format(ConvBase[-1],Neighborsi))
            if ConvBase[-1] in Neighborsi :
                ConvBase.append(i)
                Rest.remove(i)

    return ConvBase

# ===================================================================

def Neighbors_List(numAtom,Structure) :
    # Returns the list of all atom's indexes in Structure wich are neighbors of Atom indexed numAtom (numAtom must be int, 0 included)
    set=np.ndarray.tolist(Structure.positions)
    radAtom=esph.Atom_Radius(Structure,numAtom,1)
    posAtom=set[numAtom]
    Neilist=[]
    for i in range(0,len(set)) :
        int(i)
        radi=esph.Atom_Radius(Structure,i,1)
        D = tool.distance(set[i],posAtom) - radi - radAtom
        if D < 0.1 :
            Neilist.append(i)
    Neilist.remove(numAtom)
    return Neilist

# ===================================================================
def Hull_Neighbors_List(O,Structure,hull)  :
    # Returns the list of all atom's indexes in hull wich are neighbors of Atom indexed numAtom (numAtom must be int, 0 included)
    Allfacets = np.ndarray.tolist(hull.simplices)
    Hull_Neilist = []
    print(Allfacets)
    for facet in Allfacets :
        if O in facet :
            for index in facet :
                if index not in Hull_Neilist :
                    Hull_Neilist.append(index)
    Hull_Neilist.remove(O)
    return Hull_Neilist
# ===================================================================

def facet_cap(O,Structure,hull) :
    #From the hull, the Structure and the index of the sphere center O (int), return the proportion of sphere cuted by hull facets
    #(ie the caps out of the hull). cap is returned proportionnal to total sphere volume (ie cap in [0,1[)
    cap = 0
    Center = Structure.positions[O]
    Radius = esph.Atom_Radius(Structure,O,1)
    AllFacets = np.ndarray.tolist(hull.equations)
    hullplan=tool.cleanlist(AllFacets)
    cutplan = [] # List of all cuting plan equations
    cuthigh = [] # List of cuting high, corelated to cutplan
    cutpoint = [] #List containing for each plan one point of the hull included in this plan (used for overlap routine)
    for facet_plan in AllFacets:
        d=tool.dist_point_plan(Center,facet_plan)
        if d<Radius : #Then we can have a cap out : but first verify that the facet isn't to far away
            plan_index = AllFacets.index(facet_plan)
            FacetPtsIndex = np.ndarray.tolist(hull.simplices[plan_index])
            FacetPtsCoord=[]
            for facetpt in FacetPtsIndex :
                FacetPtsCoord.append(np.ndarray.tolist(hull.points[facetpt]))
            Facetnorm=facet_plan[:3]
            GoOnPlan=np.multiply(Facetnorm,d/tool.vector_norm(Facetnorm))
            OrthoProjO = np.ndarray.tolist(hull.points[O] + GoOnPlan) #OrthoProjO is the orthogonal projection of O on facet_plan
            # Sphere cut by facet only if OrtoProjO is in the facet : sum of angles with facet points must be 2*pi
            V0 = tool.vector_def(OrthoProjO, FacetPtsCoord[0])
            V1 = tool.vector_def(OrthoProjO, FacetPtsCoord[1])
            V2 = tool.vector_def(OrthoProjO, FacetPtsCoord[2])
            SumAngle= np.abs(tool.angle_vector(V0,V1))+np.abs(tool.angle_vector(V0,V2))+np.abs(tool.angle_vector(V1,V2))
            print("Sum of angle : {} = {}°".format(SumAngle,SumAngle*180/math.pi))

            if SumAngle > 1.99 * math.pi :#considering a little error here from 1% : because soon OrthoProjO is in a facet ridge
                if facet_plan not in cutplan :#If OrthoProjO in a ridge, the plan can be counted twice or more...
                    h=Radius-d
                    cutplan.append(facet_plan)
                    cuthigh.append(h)
                    cutpoint.append(FacetPtsCoord[0]) # Just one point is enough :

                # Draw to see better :===================
                Draw = FacetPtsCoord
                Draw.append(np.ndarray.tolist(hull.points[O]))
                Draw.append(OrthoProjO)
                DrawView = Atoms("X3CH", positions=Draw)
                view(DrawView)
                # =======================================



            """________________________________________________________________
            # Little thing to see exactly what happens
            hullview = []
            for i in hull.vertices:
                hullview.append(Structure.positions[i])
            L = len(hullview)
            Lookhull = Atoms(positions=hullview)
            set = []
            for ffp in AllFacets :
                if ffp == cutplan :
                    FacetPts = np.ndarray.tolist(hull.simplices[plan_index])
                    for pt in FacetPts:
                        if pt not in set :
                            set.append(np.ndarray.tolist(hull.points[pt]))
            nameset = "H" + str(len(set))
            ViewFacet = Atoms(nameset, positions=set)
            Ostr = Atoms("C", positions=[np.ndarray.tolist(hull.points[O])])
            view(Lookhull + ViewFacet + Ostr)
            # __________________________________________________________________"""

    if len(cutplan) <1 :
        print("Finally no plan cutting our sphere...\n\n" )

    if len(cutplan) == 1 :#Only one plan cuting sphere : One cap to be calculated :
        h = cuthigh[0]
        print("Sphere (radius {}) is cuted by one plan at cap from high {}\n Plan equation :{} ".format(Radius, h,
                                                                                                              facet_plan))

        Vtot = 4 * math.pi * (Radius ** 3) / 3
        Vcap = math.pi * (h ** 2) / 3 * (3 * Radius - h)
        cap += Vcap / Vtot

    if len(cutplan) > 1:
        print("We have {} plans cutting the sphere.... We have to calculate overlap... ".format(len(cutplan)))

    return cap
