# coding: utf-8

import unittest
import math
import numpy as np
import delaunay.core as delc
import es_sym_analys as essyma
import es_tools as tool
import es_clustering as esclus
from scipy.spatial import ConvexHull, Voronoi

# ===========
from ase import Atoms
from ase.visualize import view
from ase.data import covalent_radii

# ===================================================================
# List of routines :
"""
=================
 ES_AllRadius_Updater(NewES,Structure,[list])  :   Update ES_AllRadius global variable with new radius of empty spheres 
                                                   given as NewES

 Voronoi_Vertex(Structure)          :   Computes Voronoi Vertices of Structure
 Delaunay_Tetrahedral_ES(Structure,[minsize],[maxsize],[tol]) :   Creates a tetrehedral mesh from the structure,
                                                            then returns for each center the perfect sphere going in.
 Convex_Hull_Cover(Structure,[es_radius],[tol],[missing])     :   Finds the exterior Hull from the set, create triangular 
                                                                  mesh then returns cover coordinates. tol=0 => no fusion
 Select_Int_Ext(Centroid,E1,E2,IE)  :   Clean the Cover, taking only internal or external
 Internal_Hull_Cover(Structure,[es_radius],[tol],[missing])   :   Finds the interior Hull from the set, create triangular 
                                                                  mesh then returns cover coordinates
 Internal_Hull_Centers(set)         :   Finds the interior Hull from the set, create triangular mesh then returns centers coordinates
 ES_Fusion(set, structure, size)    :   Change the set by clustering spheres near from size to each other. No size given => take shortest
                                        Maintain distances with structure, compared with the ancient set.
 Fusion_Overlap(Spheres_Data,tol)   :   Find Spheres touching each other, and fusions them. Don't return  radius : only final coordinates  
 Flat_Covering(Structure,[R],[tol],[Curved])    :   For flat (or almost) set : Searchs major plane, triangulates, 
                                                    and cover the 2 sides.
 Plane_Triangulation(Plane3D,Plane_eq): Return triangulation of a 3D Plane (convert into 2D, uses Delny)
 Atom_Radius(Structure,n,list)      :   Returns radius of n°th atom in Structure (Angstrom). Regroup different radius lists.
 Convex_Hull_InterCover(set)        :   Return list of internal cover using ConvexHull : Different from Delaunay_Intersphere : 
                                        made for empty clusters

=================
"""

ES_AllRadius = []  # Global Variable taking all empty spheres radius


# ===================================================================
def ES_AllRadius_Updater(NewES, Structure, list=1):
    # Update ES_AllRadius global variable with new radius of empty spheres given as NewES
    global ES_AllRadius
    StrPos = np.ndarray.tolist(Structure.positions)
    for ES in np.ndarray.tolist(NewES.positions):
        # print("ES = {}".format(ES))
        Tempo = []
        for A in StrPos:
            # print("A = {}".format(A))
            d = tool.distance(ES, A)
            Tempo.append(d - Atom_Radius(Structure, StrPos.index(A), list))
        ES_AllRadius.append(min(Tempo))
    return ES_AllRadius


# ===================================================================
def Voronoi_Vertex(Structure):
    # Uses Delaunay triangulation to create a tetrahedral mesh from the set of
    # atoms-centers coordinates, then computes each tetrahedron's inerty center
    # Returns list of tetrahedrons-centers coordinates
    # Be careful on tetrahedralisation done : the cube for example will not be correctly tesselated (6tetrahedrons in)
    """
        Triag=delc.Triangulation(set)
        tetracenters=[]
        for tetra in Triag.indices:
            x = y = z = 0.0
            for vertex in tetra:
                x += set[vertex][0] / 4.0
                y += set[vertex][1] / 4.0
                z += set[vertex][2] / 4.0
            tetracenters.append((x, y, z))
        """
    struct = np.ndarray.tolist(Structure.positions)
    Vor = Voronoi(struct)
    return np.ndarray.tolist(Vor.vertices)


# ===================================================================
def Delaunay_Tetrahedral_ES(Structure, minsize=0, maxsize=999, tol=0.6):
    # Uses Delaunay triangulation to create a tetrahedral mesh from the Structure of
    # atoms-centers coordinates, then adds empty sphere in each tetrahedron, equidistant to all Atoms
    es_data = []
    allradius = []
    set = Structure.positions
    All_Spheres_Data = esclus.Spheres_Data_Structure_Extractor(Structure, 1)
    # print All_Spheres_Data
    Triag = delc.Triangulation(set)
    for tetra in Triag.indices:
        Data = []
        for vertex in tetra:
            Data.append(All_Spheres_Data[vertex])
        tetra_es = esclus.Tetrahedron_ES(Data)
        if tetra_es == 999:
            EV = []
            for vertex in tetra:
                EV.append(vertex)
            print(
            "Error by computing tangent solution to tetrahedron {} : delta < 0 for radius determination ".format(EV))
        elif tetra_es != 999:  # 999 is returned when the problem has no solution due to singular matrix
            if tetra_es[1] <= maxsize:
                if tetra_es[1] >= minsize:
                    es_data.append(tetra_es)
                    # output.append(tetra_es[0])

            """Verify result : ________________
            print("E_S created at position {}, in center of tetrahedron {}".format(output[0],tetra))
            viewer = [output[-1]]
            for vertex in tetra :
                d=tool.distance(set[vertex],tetra_es[0])
                print("Distance with vertex {} : {}\nAtom radius = {}, ES radius = {}, so the sum is : {}".format(vertex,d,All_Spheres_Data[vertex][1],tetra_es[1],All_Spheres_Data[vertex][1]+tetra_es[1]))
                #viewer.append(set[vertex])
            #View=Atoms("XC4",positions = viewer)
            #view(View)
            #raw_input("\nPress Enter to continue ...\n")
            #"""  # _______________________________

    # print "allradius : \n",allradius
    output = Fusion_Overlap(es_data, tol)
    return output


# ===================================================================
def Convex_Hull_Cover(Structure, radius_es=0, tol=0.6, missing=False, Int_Ext=0):
    # Uses ConvexHull to find the triangular Hull from the set, then generate a covering by adding
    # an empty sphere on each hull triangle.
    # Default tol value is used for Fusion Overlap
    # Returns list of cover coordinates
    if Int_Ext == 0:
        print("Select wich covering you wish :\n0: Internal Cover\n1: External Cover\n2: Both Cover")
        CovChoice = input()

    Cover_Data = []
    set = np.ndarray.tolist(Structure.positions)
    if radius_es == 0:
        radius_es = input("Please select the radius of empty spheres you desire :    ")
    hull = ConvexHull(set)
    xc = yc = zc = 0
    lh = len(hull.vertices)
    for hpt in hull.vertices:
        xc += set[hpt][0] / lh
        yc += set[hpt][1] / lh
        zc += set[hpt][2] / lh
    Centroid = [xc, yc, zc]

    if missing == False:  # It means we don't care if some hull points are not implemented because they are in facet
        # Computes the centroïd of the hull

        AllData = esclus.Spheres_Data_Structure_Extractor(Structure, 1)
        for facet in hull.simplices:  # hull.simplices contains index of simplices vertex, grouped by 3.
            Data = []
            for vertex in facet:
                Data.append(AllData[vertex])
            ES1, ES2 = esclus.Triangle_ES(Data, radius_es)
            if ES1 == 666:
                print("For radius {}, no solution of tangent to triangle {} ".format(radius_es, facet))
            if ES1 != 666:
                ES = Select_Int_Ext(Centroid, ES1, ES2, 1)  # Last parameter != 0 => external cover
                if tool.distance(Centroid, ES) < 1000000000:
                    Cover_Data.append([ES, radius_es])



    elif missing == True:
        Data = esclus.Spheres_Data_Structure_Extractor(Structure, 1)
        HullPlane = tool.cleanlist(np.ndarray.tolist(hull.equations))
        for HP in HullPlane:
            a, b, c, d = HP
            HPList = []
            HPIndex = []
            for pt in set:
                x, y, z = pt
                if abs(a * x + b * y + c * z + d) < 0.01:  # Then pt is in Plane HP
                    HPList.append(pt)
                    HPIndex.append(set.index(pt))

            Tria = Plane_Triangulation(HPList, HP)
            for tria in Tria.indices:
                Spheres_data = [Data[HPIndex[tria[0]]], Data[HPIndex[tria[1]]], Data[HPIndex[tria[2]]]]
                ES1, ES2 = esclus.Triangle_ES(Spheres_data, radius_es)
                if ES1 == 666:
                    print("For radius {}, no solution of tangent to triangle {} ".format(radius_es, [HPIndex[tria[0]],
                                                                                                     HPIndex[tria[1]],
                                                                                                     HPIndex[tria[2]]]))
                if ES1 != 666:
                    # ES = Select_Int_Ext(Centroid, ES1, ES2, 1)  # Last parameter != 1 => external cover

                    # Cover_Data.append([ES, radius_es])
                    if CovChoice > 1:
                        Cover_Data.append([ES1, radius_es])
                        Cover_Data.append([ES2, radius_es])
                    else:
                        ES = Select_Int_Ext(Centroid, ES1, ES2, CovChoice)  # 1 : external; 0 : internal
                        Cover_Data.append([ES, radius_es])

        """Verify result : ________________
        print("E_S created at position {}, defined on triangle {}".format(ES, facet))
        viewer = [ES,ES2]
        for vertex in facet:
            d = tool.distance(set[vertex], ES)
            print("Distance with vertex {} : {}\nAtom radius = {}, ES radius = {}, so the sum is : {}".format(vertex, d,AllData[vertex][1],radius_es,AllData[vertex][1]+radius_es))
            viewer.append(set[vertex])
        View=Atoms("XH4",positions = viewer)
        view(View)
        # raw_input("\nPress Enter to continue ...\n")
        # """  # _______________________________

        # Fusion overlapping spheres

    Output = Fusion_Overlap(Cover_Data, tol)  # tol at 1% : means fusion if d < (r1 + r2) * 0.01

    return Output


# ===================================================================
def Select_Int_Ext(Centroid, E1, E2, IE):
    # Returns only internal or external part of cover, using the fact hull is convex : so using his centroid
    # IE = 0 : take internal part, IE != 0 : take external part
    d1 = tool.distance(Centroid, E1)
    d2 = tool.distance(Centroid, E2)
    if IE == 0:  # Internal part : the closest to centroid
        if d1 < d2:
            return E1
    else:  # External part : the farest from centroïd
        if d2 < d1:
            return E1
    # Excepted this 2 double conditions : we have to take E2
    return E2


# ===================================================================


def Internal_Hull_Cover(Structure, radius_es=0, tol=0.6, missing=False):
    # Uses ConvexHull to find the triangular Hull from the set, then generate a covering by adding
    # an empty sphere on each hull triangle.
    # Default tol value is used for Fusion Overlap
    # Returns list of cover coordinates
    Cover_Data = []
    set = np.ndarray.tolist(Structure.positions)

    if radius_es == 0:
        radius_es = input("Please select the radius of empty spheres you desire :  ")

    hull = ConvexHull(set)

    xc = yc = zc = 0
    lh = len(hull.vertices)
    for hpt in hull.vertices:
        xc += set[hpt][0] / lh
        yc += set[hpt][1] / lh
        zc += set[hpt][2] / lh
    Centroid = [xc, yc, zc]

    invset = tool.Invert_Coord(set, Centroid, 10)
    hull = ConvexHull(invset)

    if missing == False:  # It means we don't care if some hull points are not implemented because they are in facet
        # Computes the centroïd of the hull

        AllData = esclus.Spheres_Data_Structure_Extractor(Structure, 1)
        for facet in hull.simplices:  # hull.simplices contains index of simplices vertex, grouped by 3.
            Data = []
            for vertex in facet:
                Data.append(AllData[vertex])
            ES1, ES2 = esclus.Triangle_ES(Data, radius_es)
            if ES1 == 666:
                print("For radius {}, no solution of tangent to triangle {} ".format(radius_es, facet))
            if ES1 != 666:
                ES = Select_Int_Ext(Centroid, ES1, ES2, 0)  # Last parameter = 0 => internal cover

                Cover_Data.append([ES, radius_es])

    elif missing == True:
        Data = esclus.Spheres_Data_Structure_Extractor(Structure, 1)
        HullPlane = tool.cleanlist(np.ndarray.tolist(hull.equations))
        for HP in HullPlane:
            a, b, c, d = HP
            HPList = []
            HPIndex = []
            for pt in set:
                x, y, z = pt
                if abs(a * x + b * y + c * z + d) < 0.01:  # Then pt is in Plane HP
                    HPList.append(pt)
                    HPIndex.append(set.index(pt))

            Tria = Plane_Triangulation(HPList, HP)
            for tria in Tria.indices:
                Spheres_data = [Data[HPIndex[tria[0]]], Data[HPIndex[tria[1]]], Data[HPIndex[tria[2]]]]
                ES1, ES2 = esclus.Triangle_ES(Spheres_data, radius_es)
                if ES1 == 666:
                    print("For radius {}, no solution of tangent to triangle {} ".format(radius_es, [HPIndex[tria[0]],
                                                                                                     HPIndex[tria[1]],
                                                                                                     HPIndex[tria[2]]]))
                if ES1 != 666:
                    ES = Select_Int_Ext(Centroid, ES1, ES2, 0)  # Last parameter != 0 => internal cover

                    Cover_Data.append([ES, radius_es])

        """Verify result : ________________
        print("E_S created at position {}, defined on triangle {}".format(ES, facet))
        viewer = [ES,ES2]
        for vertex in facet:
            d = tool.distance(set[vertex], ES)
            print("Distance with vertex {} : {}\nAtom radius = {}, ES radius = {}, so the sum is : {}".format(vertex, d,AllData[vertex][1],radius_es,AllData[vertex][1]+radius_es))
            viewer.append(set[vertex])
        View=Atoms("XH4",positions = viewer)
        view(View)
        # raw_input("\nPress Enter to continue ...\n")
        # """  # _______________________________

        # Fusion overlapping spheres

    Output = Fusion_Overlap(Cover_Data, tol)  # tol at 1% : means fusion if d < (r1 + r2) * 0.01

    return Output


# ===================================================================
def Internal_Hull_Centers(set):
    # Uses ConvexHull to find the intern triangular Hull from the set, then computes all centers of simplices
    # Returns list of centers coordinates
    invset = tool.Invert_Coord(set, [0, 0, 0], 10)
    hull = ConvexHull(invset)
    output = []
    for facet in hull.simplices:
        x = y = z = 0.0
        for vertex in facet:
            x += set[vertex][0] / 3.0
            y += set[vertex][1] / 3.0
            z += set[vertex][2] / 3.0
        facet_center = [x, y, z]  # Center of the triangular facet
        output.append(facet_center)
    output = np.array(output).tolist()
    return output


# ===================================================================

# ===================================================================


# ===================================================================
def ES_Fusion(set, structure, size=0):
    # study the given set, and fusion some empty spheres to create a better set. Size is used for the partitionnal
    # clustering, and structure assures the clustering will not reduce the distances. It it set basically to the min distance in the set.
    if size == 0:
        size = tool.shortest_dist(set)
        # print("initial size :", size)
        size = size * 11. / 10  # we add 10%, to include the very similar cases
        # print("with error correction  size :", size)
    fusion_set = []  # output : defined like set
    dmin = tool.set_set_proximity(set, structure)
    hull = ConvexHull(structure)
    # simplice_centers=Convex_Hull_Centers(structure+set)
    while len(set) > 0:
        # Define a new cluster, add as much empty spheres as possible, and regroup around the centroid
        cluster = [set[0]]
        centroid = cluster[0]  # Initialisation of the next cluster (it may rest one element, to progress until set=void
        set.pop(0)
        reroll = 1
        while reroll == 1:
            d0_error = 0
            reroll = 0  # We must scan the set everytime we change centroid, or we could miss some ES in set
            for ES in set:
                # Other possible condition :   tool.point_set_proximity(ES, cluster)<=size
                if tool.distance(ES, centroid) < size:  # We can fusion to the cluster
                    reroll = 1  # Centroid will be updated, so reroll the scan of set
                    print("Fusionned a sphere to the cluster")
                    cluster.append(ES)
                    set.remove(ES)  # It is in the cluster, so remove from set : we studied it
                    centroid = tool.Isobarycenter(cluster)
                    if tool.point_set_proximity(centroid,
                                                structure) < dmin:  # We have to put centroid farer to balance fusion
                        Nearest = tool.search_nearest(centroid, structure,
                                                      tool.point_set_proximity(centroid, structure))
                        V = np.ndarray.tolist(
                            np.array(centroid) - np.array(Nearest))  # So we need nearest structure point
                        d = tool.distance(centroid, Nearest)
                        if d == 0:  # it means the centroid came right on existing structure
                            print("Cluster centered exactly in the structure. Size must be revised : cluster cancelled")
                            d0_error = 1
                            fusion_set += cluster
                            reroll = 0
                            break
                        else:
                            V = np.multiply(V, 1. / d)  # set V as norm 1
                            V = np.multiply(V, dmin)
                            #
                            print("\n\n We put away form:\n{}\ndmin={}\n".format(tool.vector_norm(V), dmin))
                            #
                            centroid = np.ndarray.tolist(
                                np.array(Nearest) + V)  # get centroid away from Nearest to dmin
                            #
                            #
                            #
        if d0_error == 0:
            fusion_set.append(centroid)
            #
    #
    return fusion_set


# ===================================================================
def Fusion_Overlap(Spheres_Data, tol):
    # Find Spheres touching each other, and fusions them. Don't return  radius : only final coordinates
    Output = []
    ls = len(Spheres_Data)
    Index = range(ls)
    for iS in range(ls):
        if iS in Index:  # The we have to treat this case
            FusionIndex = [iS]
            for iOS in range(iS + 1, ls):
                if iOS in Index:
                    S = Spheres_Data[iS][0]
                    OS = Spheres_Data[iOS][0]
                    rS = Spheres_Data[iS][1]
                    rOS = Spheres_Data[iOS][1]
                    # print("S : {}\nOS : {}".format(S,OS))
                    if tool.distance(S, OS) < (rS + rOS) * tol:
                        # print("Overlap detected : d= {}, r1 ={}, r2 = {}".format(tool.distance(S, OS),rS,rOS))
                        Index.remove(iOS)  # S and OS are same coord or almost : we remove the last : OS
                        FusionIndex.append(iOS)
            lf = len(FusionIndex)
            x = y = z = 0
            for i in FusionIndex:
                x += Spheres_Data[i][0][0] / lf
                y += Spheres_Data[i][0][1] / lf
                z += Spheres_Data[i][0][2] / lf
            Output.append([x, y, z])
            # else : iS correspond to coord already fusionned
    return Output


# ===================================================================

# ===================================================================

# ===================================================================

# ===================================================================

# ===================================================================
# ===================================================================
def Flat_Covering(Structure, R=0, tol=0.6, Curved=False):
    # Designed for quite flat set : Searchs major plane, triangulates it, and cover both sides with empty spheres wich radius=size.

    if R != 0:
        NoAsk = 1
    else:
        NoAsk = 0

    if Curved == False:

        flatness = input("Please describe cluster :\nOnly one major plane : Enter 0\nMore major planes : Enter 1\n")

        set = np.ndarray.tolist(Structure.positions)
        struct = set
        FlatCover = []
        # Search major plane(s)
        if flatness != 0:
            AllPlanes = essyma.major_plane(struct, multiple=True)

        else:
            [Plane3D, Plane_eq] = essyma.major_plane(struct)
            AllPlanes = [[Plane3D, Plane_eq]]

        """
        for P in AllPlanes :
            PlaneView = Atoms(positions=P[0])
            #view(Structure+PlaneView)
            view(PlaneView)
            print("Plane n°{} : \nContains : {}\n PlaneEq = {}".format(AllPlanes.index(P)+1,P[0],P[1]))

        #"""

        # Build empty spheres for all major planes :

        for AP in AllPlanes:
            Plane3D, Plane_eq = AP

            if NoAsk == 0:
                Pview = Atoms(positions=Plane3D)
                view(Pview)
                print("Please select the radius of empty spheres you desire :    ")
                R = input("(See the view of current plane to get help)\n")

            Index = []
            for Ppt in Plane3D:
                Index.append(struct.index(Ppt))

            """ Show Details on Plane
            #print("Plane : Equation is {}x+{}y+{}z+{}=0\n Plane norm is {}".format(Plane_eq[0],Plane_eq[1],Plane_eq[2],Plane_eq[3],Norm))
            Lset = len(set)
            name = "C" + str(Lset)
            Structure = Atoms(name, positions=struct)
            Plane3DView = Atoms(positions=Plane3D)
            view(Structure + Plane3DView)
            #"""  # =====

            Triang = Plane_Triangulation(Plane3D, Plane_eq)
            # print("Triangulation 2D : {}".format(Triang.indices))
            # Extract DataforTangent_Fourth_Sphere fromStructure
            Data = esclus.Spheres_Data_Structure_Extractor(Structure, 1)
            simplicenters = []
            for tria in Triang.indices:
                Spheres_data = [Data[Index[tria[0]]], Data[Index[tria[1]]], Data[Index[tria[2]]]]
                """
                h =  (Spheres_data[0][1] + Spheres_data[1][1] + Spheres_data[2][1]) #h is set as 3 times the average covalent radii of atoms
                xc = 1. / 3 * (Spheres_data[0][0][0] + Spheres_data[1][0][0] + Spheres_data[2][0][0])
                yc = 1. / 3 * (Spheres_data[0][0][1] + Spheres_data[1][0][1] + Spheres_data[2][0][1])
                zc = 1. / 3 * (Spheres_data[0][0][2] + Spheres_data[1][0][2] + Spheres_data[2][0][2])
                Addpoint = np.ndarray.tolist(np.array([xc,yc,zc])+np.multiply(Norm,h))
                Spheres_data.append([Addpoint,1]) # We add one sphere with radius=0, that must be tangent to solution
                """

                P1, P2 = esclus.Triangle_ES(Spheres_data, R)
                if P1 == 666:
                    print("For radius {}, no solution of tangent to triangle {} ".format(R, [Index[tria[0]],
                                                                                             Index[tria[1]],
                                                                                             Index[tria[2]]]))
                if P1 != 666:
                    FlatCover.append([P1, R])
                    FlatCover.append([P2, R])
                """
                ViewPos=[Spheres_data[0][0],Spheres_data[1][0],Spheres_data[2][0],Addpoint,Cov1[0]]
                print "preview :",ViewPos
                Viewer=Atoms("H2OClX",positions=[Spheres_data[0][0],Spheres_data[1][0],Spheres_data[2][0],Addpoint,Cov1[0]])
                view(Viewer)
                Spheres_data.remove([Addpoint, 1])
                Addpoint = np.ndarray.tolist(np.array([xc, yc, zc]) - np.multiply(Norm, h))
                Spheres_data.append([Addpoint, 1])  # Our
                print Spheres_data
                Cov2 = esclus.Tetrahedron_ES (Spheres_data)
                ViewPos = [Spheres_data[0][0], Spheres_data[1][0], Spheres_data[2][0], Addpoint, Cov2]
                print "preview :", ViewPos
                Viewer = Atoms("H2OClX", positions=[Spheres_data[0][0], Spheres_data[1][0], Spheres_data[2][0], Addpoint, Cov2[0]])
                view(Viewer)
                FlatCover.append(Cov1[0])
                FlatCover.append(Cov2[0])

            """
            """
            for tria in Triang.indices: #compute classic simplice center :
                x=y=z=0
                for vertex in tria:
                    x += set[vertex][0] / 3.0
                    y += set[vertex][1] / 3.0
                    z += set[vertex][2] / 3.0
                simcen=[x,y,z]
                for simcen in simplicenters :
                    C1=np.ndarray.tolist(np.array(simcen) + np.array(Norm))
                    C2=np.ndarray.tolist(np.array(simcen) - np.array(Norm))
                    FlatCover.append(C1)
                    FlatCover.append(C2)
            #"""  #



    elif Curved == True:  # With curve : Complicated to define planes... So we need to do as if it was pure 2D.
        # In Case of Curved routine, ALL spheres will be included in Cover iteration.
        FlatCover = []
        Plane3D = Structure.positions
        Plane2D = []
        Projection_selected = 0
        while Projection_selected == 0:
            PProj = input(
                "To create mesh, we need to project on a plane : please select him.\n1 for x=0\n2 for y=0\n3 for z=0\n")
            if PProj == 1:
                for pt in Plane3D:
                    Plane2D.append(pt[1:])
                Projection_selected = 1
            elif PProj == 3:
                for pt in Plane3D:
                    Plane2D.append(pt[:2])
                Projection_selected = 1
            elif PProj == 2:
                for pt in Plane3D:
                    Plane2D.append([pt[0], pt[2]])
                Projection_selected = 1
        Triang = delc.Triangulation(Plane2D)
        if R == 0:
            R = input("Please select the radius of empty spheres you desire :    ")

        Data = esclus.Spheres_Data_Structure_Extractor(Structure, 1)
        simplicenters = []
        for tria in Triang.indices:
            Spheres_data = [Data[tria[0]], Data[tria[1]], Data[tria[2]]]
            P1, P2 = esclus.Triangle_ES(Spheres_data, R)
            if P1 == 666:
                print("For radius {}, no solution of tangent to triangle {} ".format(R, [tria[0], tria[1],
                                                                                         tria[2]]))
            if P1 != 666:
                FlatCover.append([P1, R])
                FlatCover.append([P2, R])

    FlatCover = Fusion_Overlap(FlatCover, tol)

    return FlatCover


# ====================================================================
def Plane_Triangulation(Plane3D, Plane_eq):
    # Transform plane into same z
    Norm = Plane_eq[:3]
    a, b, c, d = Plane_eq
    Plane2D = []
    NormZ0 = [0, 0, 1]
    Alpha = tool.angle_vector(Norm, NormZ0)  # Angle of rotation

    if Alpha % math.pi == 0:  # Plane already at z=K : no rotation needed
        for pt in Plane3D:
            Plane2D.append(pt[:2])
    else:  # We must see the axis of rotation, then rotate all points :
        if a == 0:
            u = [1, 0, 0]
            M = [-d / b, 0, 0]
        else:
            u = [-b / a, 1, 0]
            M = [-d / a, 0, 0]
        # print("Norm of Plane = {}, so : \nAlpha = {}° \nu={} and M = {}".format(Norm,Alpha * 180 / math.pi, u,M))
        # Plane def by ax + by + cz + d = 0. Intercect z=0 at array ax + by + d = 0, so vector directing intercection is (-b,a,0)
        nu = tool.vector_norm(u)
        u = np.ndarray.tolist(np.multiply(u, 1. / nu))
        # print("Intercection is line defined by vector u={} and point p={}\nAngle of rotation will be {}°".format(u,M,Alpha*180/math.pi))
        # Now rotate all points with rotation angled Alpha round u.
        Rview = []
        for pt in Plane3D:  # Translate point until axis of rotation include origin : for easier rotation
            x, y, z = pt
            Mm = [-M[0], -M[1], -M[2]]
            pt = tool.vector_trslt(pt, Mm)
            rpt = tool.rot3D(pt, Alpha, u)
            rpt = tool.vector_trslt(rpt, M)

            # print("pt = {}, rotate to rpt = {} (Verify same z !)".format(pt,rpt))
            Rview.append(rpt)
            Plane2D.append(rpt[:2])

    # Triangulate with Delny : indices will be the same as original Plane : no need to invert transformation
    Triang = delc.Triangulation(Plane2D)
    return Triang


# ====================================================================

# ====================================================================
def Atom_Radius(Structure, n, list):
    # Returns the radius of the n°th atom in Structure. 0 are placed to keep information. Unit = Angstrom
    # List variable determines wich information we need
    global ES_AllRadius

    if 0 in Structure.numbers:
        FirstES = np.ndarray.tolist(Structure.numbers).index(0)  # Number of the first Empty_Spheres
    N = Structure.numbers[n]
    if N == 0:  # Atom is X : empty sphere
        Atom_Radius = ES_AllRadius[n - FirstES]

    if list == 1:  # Covalent Radii
        Atom_Radius_List = covalent_radii
    else:
        print("No list found, verify list variable. Routine returns 0")
        return 0
    if Atom_Radius_List[N] == 0:
        print ("No information on this atom, or false n° given. Routine returns 0")
    return Atom_Radius_List[N]


# ===================================================================

def Convex_Hull_InterCover(set):
    # Uses ConvexHull to find the triangular Hull from the set, then generate a covering by adding
    # an empty sphere on each hull triangle.
    # Returns list of cover coordinates
    reverset = tool.Invert_Coord(set, [0, 0, 0], 2)
    hull = ConvexHull(reverset)
    cover_coord = []
    counter = 0
    for facet in hull.simplices:  # hull.simplices contains index of simplices vertex, grouped by 3.
        x = y = z = 0.0
        for vertex in facet:
            x += reverset[vertex][0] / 3.0
            y += reverset[vertex][1] / 3.0
            z += reverset[vertex][2] / 3.0
        facet_center = [x, y, z]  # Center of the triangular facet
        normal_facet = hull.equations[counter][:3]  # The exterior normal of the facet
        addpoint = facet_center + normal_facet
        cover_coord.append(addpoint)
        counter += 1
    cover_coord = np.array(cover_coord).tolist()
    cover_coord = tool.Invert_Coord(cover_coord, [0, 0, 0], 2)
    return cover_coord


# ===================================================================
def lookhull(Structure, hull):
    # view the structure and his hull
    hullview = []
    for i in hull.vertices:
        hullview.append(Structure.positions[i])
    L = len(hullview)
    Lookhull = Atoms(positions=hullview)
    view(Lookhull)
    view(Structure + Lookhull)
    return 0