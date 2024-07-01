"""
Implementation of Antonietti's algorithm 2
https://doi.org/10.1007/s10915-018-0802-y

Hats off to Thijs van Putten for coding this
"""
import numpy as np
from math import sqrt
from numba import jit

def integrate_monomials_polyhedron(node_coords, elem_faces, orders, restrict_complete):
    """
    Integrate all monomials up to specified orders (3 orders in u,v,w) over a polyhedron
    in parametrized u,v,w space (range -1 to 1).
    Polyhedron is defined by node coordinates, and faces specified as a 2D array of indices into node_coords
    If restrict complete is set; restrict monomials to complete basis only
    """
    node_monomials = evaluate_monomials(node_coords, orders, restrict_complete)
    num_subgeom, center, euclid_dist, index_offset, node_indices = hierarchical_geometry_polyhedron(node_coords, elem_faces)
    return integrate_monomials_recursive(3, 0, orders, restrict_complete, num_subgeom, center, euclid_dist, index_offset, node_indices, node_monomials)

def integrate_monomials_triangle(node_coords, orders, restrict_complete):
    """
    Integrate all monomials up to specified orders (3 orders in u,v,w) over a triangle
    in parametrized u,v,w space (range -1 to 1).
    Triangle is defined by 3 sets of node coordinates
    If restrict complete is set; restrict monomials to complete basis only
    """
    node_monomials = evaluate_monomials(node_coords, orders, restrict_complete)
    num_subgeom, center, euclid_dist, index_offset, node_indices = hierarchical_geometry_triangle(node_coords)
    result = integrate_monomials_recursive(2, 0, orders, restrict_complete, num_subgeom, center, euclid_dist, index_offset, node_indices, node_monomials)
    return result

def integrate_monomials_mesh(vertices, faces, orders, restrict_complete):
    """
    Integrate all monomials up to specified orders (3 orders in u,v,w) over a triangle mesh
    in parametrized u,v,w space (range -1 to 1).
    Each triangle is defined by a face which specifies 3 vertices
    If restrict complete is set; restrict monomials to complete basis only
    """
    vertex_monomials = evaluate_monomials(vertices, orders, restrict_complete) 
    result = np.zeros((len(faces), orders[0]+1, orders[1]+1, orders[2]+1))
    for iface in range(len(faces)):
        node_monomials = vertex_monomials[faces[iface],:,:,:]
        num_subgeom, center, euclid_dist, index_offset, node_indices = hierarchical_geometry_triangle(vertices[faces[iface]])
        result[iface] = integrate_monomials_recursive(2, 0, orders, restrict_complete, num_subgeom, center, euclid_dist, index_offset, node_indices, node_monomials)
    return result

@jit(nopython=True, cache=True, nogil=True)
def evaluate_monomials(points, orders, restrict_complete):
    """
    evaluate all monomials up to supplied orders in u,v,w at supplied points
    If restrict complete is set; restrict monomials to complete basis only
    """
    out = np.zeros((points.shape[0], orders[0]+1, orders[1]+1, orders[2]+1))
    for i in range(orders[0]+1):
        j_max = orders[1]-restrict_complete*i
        for j in range(j_max+1):
            k_max = orders[2]-restrict_complete*(i+j)
            for k in range(k_max+1):
                out[:,i,j,k] = points[:,0]**i * points[:,1]**j * points[:,2]**k
    return out

@jit(nopython=True, cache=True, nogil=True)
def hierarchical_geometry_polyhedron(node_coords, elem_faces):
    """
    Construct hierarchical geometry representation of an element defined by nodes node_coords
    (in paramterized u,v,w coordinates range -1 to 1) and faces elem_faces defined as lists of
    indices into the nodes
    """
    # set number of subgeometries per dimension
    nfaces = len(elem_faces)
    num_subgeom = np.array([0,2,3,nfaces], np.int32)
    # create arrays for center, euclidean distance and node indices
    center = np.zeros((nfaces + 3*nfaces, 3)) # unneeded for 3D, nfaces coords for 2D, 3*nfaces coords for 1D, unneeded for 0D
    euclid_dist = np.zeros(nfaces + 3*nfaces + 2*3*nfaces) # unneeded for 3D, nfaces vals for 2D, 3*nfaces vals for 1D, 2*3*nfaces vals for 0D
    node_indices = np.zeros(2*3*nfaces,np.int32) # 2 nodes per edge, 3 edges per face
    index_offset = np.array([nfaces + 3*nfaces, nfaces, 0, 0], np.int32) # offset into center and euclid_dist per dimension
    # loop through the faces
    for iface, facenodes in enumerate(elem_faces):
        fill_hierarchical_geometry_triangle(node_coords, iface, facenodes, center, euclid_dist, node_indices, index_offset)
    return num_subgeom, center, euclid_dist, index_offset, node_indices

@jit(nopython=True, cache=True, nogil=True)
def hierarchical_geometry_triangle(node_coords):
    """
    Construct hierarchical geometry representation of a triangle defined by nodes node_coords
    (in paramterized u,v,w coordinates range -1 to 1)
    """
    # set number of subgeometries per dimension
    num_subgeom = np.array([0,2,3], np.int32)
    # create arrays for center, euclidean distance and node indices
    center = np.zeros((1 + 3, 3)) # 1 coord for 2D, 3 coords for 1D
    euclid_dist = np.zeros(1 + 3 + 2*3) # 1 val for 2D, 3 vals for 1D, 2*3 vals for 0D
    node_indices = np.array([0, 1, 1, 2, 2, 0], np.int32) # node inidices of the two nodes per edge
    index_offset = np.array([4, 1, 0], np.int32) # offset into center and euclid_dist per dimension
    # for single triangle, iface is 0 and facenodes are just 0,1,2
    iface = 0
    facenodes = [0, 1, 2]
    fill_hierarchical_geometry_triangle(node_coords, iface, facenodes, center, euclid_dist, node_indices, index_offset)
    return num_subgeom, center, euclid_dist, index_offset, node_indices

@jit(nopython=True)
def fill_hierarchical_geometry_triangle(node_coords, iface, tri_nodes, center, euclid_dist, node_indices, index_offset):
    """
    Fill hierarchical geometry arrays for a triangle with index iface, and nodes tri_nodes referreing to coordinates node_coords,
    Arrays to be filled are center, euclid_dist and node_indices
    """
    # properties of the triangle
    face_center = (node_coords[tri_nodes[0]] + node_coords[tri_nodes[1]] + node_coords[tri_nodes[2]]) / 3
    face_normal = normalize(np.cross(node_coords[tri_nodes[1]] - node_coords[tri_nodes[0]], node_coords[tri_nodes[2]] - node_coords[tri_nodes[1]]))
    distance = np.dot(face_center, face_normal)
    center[index_offset[2] + iface] = face_center
    euclid_dist[index_offset[2] + iface] = distance

    # loop through the edges
    for iedge in range(3):
        edgenode0 = tri_nodes[iedge]
        edgenode1 = tri_nodes[iedge+1] if iedge < 2 else tri_nodes[0]
        edge_center = (node_coords[edgenode0] + node_coords[edgenode1]) / 2
        edge_normal = normalize(np.cross(node_coords[edgenode1] - node_coords[edgenode0], face_normal))
        distance = np.dot(edge_center - face_center, edge_normal)
        edge_index = iface * 3 + iedge
        center[index_offset[1] + edge_index] = edge_center
        euclid_dist[index_offset[1] + edge_index] = distance
        # loop through the nodes
        for ipoint, node in enumerate([edgenode0, edgenode1]):
            node_normal = normalize(node_coords[node] - edge_center)
            distance = np.dot(node_coords[node] - edge_center, node_normal)
            point_index = edge_index * 2 + ipoint
            euclid_dist[index_offset[0] + point_index] = distance
            node_indices[point_index] = node

@jit(nopython=True)
def normalize(vec):
    norm = sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)
    return vec/norm

@jit(nopython=True, cache=True, nogil=True)
def integrate_monomials_recursive(dim, geom_index, orders, restrict_complete, num_subgeom, center, euclid_dist, index_offset, node_indices, node_monomials):
    """
    Integrate monomials of given orders (in u,v,w) over a hierarchical geometry structure
    Evaluation of the monomials at the nodes should be precomputed and given as node_monomials
    This is algorithm 2 of Antonietti et al
    Arguments:
    - dim: dimension of the current geometry
    - geom_index: index of the current geometry in lists of center & euclidean distance
    - orders (u,v,w): orders of the monomials being integrated
    - restrict_complete: restrict monomials to complete basis only
    - num_subgeom: number of subgeometries per dimension
    - center: center of each geometry
    - euclid_dist: euclidean distance of each geometry to the cetner of the parent
    - index_offset: offset into center and euclid_dist per dimensiont
    - node_index: index in node_monomials of all 0D geometries (ie: nodes)
    - node_monomials: values of the monomials precomputed at all nodes
    """
    # if this is a node, return the precomputed integral
    if dim == 0:
        return node_monomials[node_indices[geom_index]]

    # collect the integrals of the geometry one level lower (ie: if this is an element, get the face integrals)
    sub_integral = np.zeros((num_subgeom[dim], orders[0]+1, orders[1]+1, orders[2]+1))
    for isub in range(num_subgeom[dim]):
        subgeom_index = geom_index * num_subgeom[dim] + isub
        sub_integral[isub] = integrate_monomials_recursive(dim-1, subgeom_index, orders, restrict_complete, num_subgeom, center, euclid_dist, index_offset, node_indices, node_monomials)

    # combine to the integral of the current geometry
    integral = np.zeros((orders[0]+1, orders[1]+1, orders[2]+1))
    for i in range(orders[0]+1):
        j_max = orders[1]-restrict_complete*i
        for j in range(j_max+1):
            k_max = orders[2]-restrict_complete*(i+j)
            for k in range(k_max+1):
                # sum up integrals of subgeometry times euclidean distance
                value = 0
                for isub in range(num_subgeom[dim]):
                    subgeom_index = geom_index * num_subgeom[dim] + isub
                    value += euclid_dist[index_offset[dim-1] + subgeom_index] * sub_integral[isub,i,j,k]
                # for faces & edges, recursively apply stokes law
                if dim < 3:
                    # recursive derivative contributions from u/v/w-monomials
                    if i > 0:
                        value += center[index_offset[dim] + geom_index][0] * i * integral[i-1,j,k]
                    if j > 0:
                        value += center[index_offset[dim] + geom_index][1] * j * integral[i,j-1,k]
                    if k > 0:
                        value += center[index_offset[dim] + geom_index][2] * k * integral[i,j,k-1]
                integral[i,j,k] = 1 / (dim + i + j + k) * value                
    return integral
