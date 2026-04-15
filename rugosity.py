# requirements: shapely
# requirements: trimesh
# requirements: numpy


import Rhino.Geometry as rg
import trimesh
import numpy as np
from shapely.geometry import Polygon
from shapely.ops import unary_union
from collections import defaultdict
from trimesh.transformations import translation_matrix

def rugosity(mesh, voxel_size):
    tmesh = trimesh.Trimesh( 
        [[p.X, p.Y, p.Z] for p in mesh.Vertices],
        [[p.A, p.B, p.C] for p in mesh.Faces] )

    def mesh_faces_to_grid(mesh, cell_size):
        """
        Assign mesh faces to a 3D voxel grid.

        Parameters
        ----------
        mesh : trimesh.Trimesh
        cell_size : float

        Returns
        -------
        grid : dict
            {(i,j,k): [face_indices]}
        bbox_min : (3,) array
            grid origin
        grid_shape : (3,) array
            number of cells in each dimension
        """

        # mesh bounds
        bbox_min = mesh.vertices.min(axis=0)
        bbox_max = mesh.vertices.max(axis=0)

        grid_shape = np.ceil((bbox_max - bbox_min) / cell_size).astype(int)

        triangles = mesh.triangles  # (F,3,3)

        # triangle bounding boxes
        tri_min = triangles.min(axis=1)
        tri_max = triangles.max(axis=1)

        # convert to grid coordinates
        grid_min = np.floor((tri_min - bbox_min) / cell_size).astype(int)
        grid_max = np.floor((tri_max - bbox_min) / cell_size).astype(int)

        grid = defaultdict(list)

        for f in range(len(triangles)):
            gmin = grid_min[f]
            gmax = grid_max[f]

            for i in range(gmin[0], gmax[0] + 1):
                for j in range(gmin[1], gmax[1] + 1):
                    for k in range(gmin[2], gmax[2] + 1):
                        grid[(i, j, k)].append(f)

        return grid, (bbox_min, bbox_max), grid_shape


    def voxel_index_to_box(voxel_index, bbox_min, bbox_max, grid_shape):
        """
        Create a trimesh box for a voxel given its index in the grid.
        
        Parameters
        ----------
        voxel_index : tuple of 3 ints
            (i, j, k) voxel index in the grid
        bbox_min : array-like, shape (3,)
            Minimum coordinates of the bounding box
        bbox_max : array-like, shape (3,)
            Maximum coordinates of the bounding box
        grid_shape : tuple of 3 ints
            Number of cells along each axis (nx, ny, nz)
            
        Returns
        -------
        box : trimesh.Trimesh
            A box mesh representing the voxel
        """
        
        voxel_index = np.array(voxel_index)
        bbox_min = np.array(bbox_min)
        bbox_max = np.array(bbox_max)
        grid_shape = np.array(grid_shape)
        
        # compute cell size in each dimension
        cell_size = (bbox_max - bbox_min) / grid_shape
        
        # compute voxel center
        center = bbox_min + (voxel_index + 0.5) * cell_size
        
        # create box
        box = trimesh.creation.box(
            extents=cell_size,
            transform=translation_matrix(center)
        )
        
        return box

    geometry = []
    grid, bounds, shape = mesh_faces_to_grid( tmesh, voxel_size )

    O = bounds[0]
    D = bounds[1] - bounds[0]

    def triangle_area(a, b, c):
        return 0.5 * np.linalg.norm(np.cross(b - a, c - a))

    def project_to_plane(points, origin, normal):
        normal = normal / np.linalg.norm(normal)
        dist = np.dot(points - origin, normal)
        return points - np.outer(dist, normal)

    def plane_basis(normal):
        normal = normal / np.linalg.norm(normal)

        if abs(normal[0]) < 0.9:
            t = np.array([1,0,0])
        else:
            t = np.array([0,1,0])

        u = np.cross(normal, t)
        u /= np.linalg.norm(u)

        v = np.cross(normal, u)

        return u, v

    def to_2d(points, origin, u, v):
        d = points - origin
        x = np.dot(d, u)
        y = np.dot(d, v)
        return np.column_stack((x,y))

    rugosity_dict = {}
    global_area = 0
    global_projected = 0

    for cell, faces in grid.items():

        faces_idx = tmesh.faces[faces]
        vertices = tmesh.vertices

        # total triangle area
        total_area = 0.0
        for f in faces_idx:
            a,b,c = vertices[f]
            total_area += triangle_area(a,b,c)

        # find best fit plane
        points = vertices[np.unique(faces_idx.flatten())]
        origin = points.mean(axis=0)
        u_pts = points - origin
        eigvals, eigvecs = np.linalg.eigh(u_pts.T @ u_pts)
        normal = eigvecs[:,0]
        u,v = plane_basis(normal)

        projected_triangles = []

        # project triangles
        for f in faces_idx:
            tri = vertices[f]
            proj = project_to_plane(tri, origin, normal)
            tri2d = to_2d(proj, origin, u, v)
            projected_triangles.append(tri2d)

        # compute union area
        polys = [Polygon(tri) for tri in projected_triangles]
        union_poly = unary_union(polys)
        projected_area = union_poly.area

        # calculate rugosity
        if projected_area > 0:
            rugosity = total_area / projected_area
        else:
            rugosity = 0

        print(cell, rugosity)    
        rugosity_dict[cell] = rugosity
        global_area += total_area
        global_projected += projected_area

    global_rugosity = global_area / global_projected

    # calculate rugosity for each mesh vertex for visualisation 
    face_values = np.zeros(len(tmesh.faces))
    for cell, faces in grid.items():
        face_values[faces] = rugosity_dict[cell]

    vertex_values = np.zeros(len(tmesh.vertices)) 
    vertex_counts = np.zeros(len(tmesh.vertices))

    for face_idx, face in enumerate(tmesh.faces):
        for v in face: 
            vertex_values[v] += face_values[face_idx]
            vertex_counts[v] += 1

    vertex_values /= np.maximum(vertex_counts, 1)

    return global_rugosity, vertex_values
