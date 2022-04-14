import matplotlib.pyplot as plt
import numpy as np
import gemmi
from mpl_toolkits.mplot3d import axes3d
from sklearn.decomposition import PCA
from scipy.spatial import ConvexHull, convex_hull_plot_2d
from scipy.spatial import Delaunay

#This function recieves a PDB file as input and returns the X,Y and Z coordinates of the structure as arrays
def getProtein(file):
    stc = gemmi.read_pdb(file)
    X = []
    Y = []
    Z = []
    model = stc[0]
    for chain in model:
        for res in chain:
            for atom in res:
                coord = atom.pos
                X.append(coord[0])
                Y.append(coord[1])
                Z.append(coord[2])
    return X,Y,Z

#This function centers a set of coordinates (X,Y or Z) on the origin (0,0,0)
def centerProtein(CoordSet):
    N = len(CoordSet)
    meanCoord = sum(CoordSet)/N
    newCoordSet = []
    for coord in CoordSet:
        new_coord = coord - meanCoord
        newCoordSet.append(new_coord)
    return newCoordSet

#This funcion returns the rotation matrix that aligns two vectors, where vec1 is "moved" to vec2
def rotation_matrix_from_vectors(vec1, vec2):
    """
    Reference
    ---------
    See https://stackoverflow.com/questions/45142959/calculate-rotation-matrix-to-align-two-vectors-in-3d-space
    """
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    return rotation_matrix

#This function formats coordinates from [X1,...,Xi,...,Xn],[Y1,...,Yi,...,Yn],[Z1,...,Zi,...,Zn]
#to ([X1,Y1,Z1],...[Xi,Yi,Zi],..., [Xn,Yn,Zn])
def format_coord(X,Y,Z):
    data = list()
    for i in range(len(X)):
        data.append([X[i],Y[i],Z[i]])
    return data

#This function checks if a point (x,y,z) belong in an Ellipsoid with axis a,b,c
def eqEllipsoid(x,y,z,a,b,c):
    eq = ((x**2)/a**2) + ((y**2)/b**2) + ((z**2)/c**2)
    return eq <= 1

#This function generates half of an ellipse of axis a,b,c, returning
#arrays for coordinates X,Y and Z of pseudoatoms with ray r
def generateHalf(r,x,y,z,a,b,c):
    X = []
    Y = []
    Z = []
    while eqEllipsoid(x,y,z,a,b,c):
        X.append(x)
        Y.append(y)
        Z.append(z)
        x += 2*r
        X.append(-x)
        Y.append(y)
        Z.append(z)
    return X,Y,Z

#This function generates an ellipse of axis a,b,c for a given z coordinate,
#returning arrays for coordinates X,Y and Z of pseudoatoms with ray r
def generateSlice(r,Pz,a,b,c):
    Px,Py = 0,0
    X = []
    Y = []
    Z = []
    while eqEllipsoid(Px,Py,Pz,a,b,c):
        P = generateHalf(r,Px,Py,Pz,a,b,c)
        for i in P[0]: X.append(i)
        for i in P[1]: Y.append(i)
        for i in P[2]: Z.append(i)
        Py += 2*r
        for i in P[0]: X.append(i)
        for i in P[1]: Y.append(-i)
        for i in P[2]: Z.append(i)
    return X,Y,Z

#This function generates an ellipsoid of axis a,b,c, returning arrays
#for coordinates X,Y,Z of pseudoatoms with ray r
def generateEllipsoid(r,a,b,c):
    Px,Py,Pz = 0,0,0
    X = []
    Y = []
    Z = []
    while eqEllipsoid(Px,Py,Pz,a,b,c):
        P = generateSlice(r,Pz,a,b,c)
        for i in P[0]: X.append(i)
        for i in P[1]: Y.append(i)
        for i in P[2]: Z.append(i)
        Pz += 2*r
        for i in P[0]: X.append(i)
        for i in P[1]: Y.append(i)
        for i in P[2]: Z.append(-i)
    return X,Y,Z

#This function checks if a point is inside a convex hull
def in_hull(poly, point):
    """
    Reference
    ---------
    See https://stackoverflow.com/questions/29311682/finding-if-point-is-in-3d-poly-in-python 
    """
    hull = ConvexHull(poly)					#creates a hull with the initial points
    new_hull = ConvexHull(np.concatenate((poly, [point])))	#adds the point to the initial array and generates a new hull
    return np.array_equal(new_hull.vertices, hull.vertices)	#compares the vertices of the two hulls

#This function gets a list of all atom types in the protein
def getData(file):
    stc = gemmi.read_pdb(file)
    data = []
    model = stc[0]
    for chain in model:
        for res in chain:
            for atom in res:
                data.append(atom.name)
    return data
