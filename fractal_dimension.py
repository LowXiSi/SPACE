# requirements: numpy
# requirements: trimesh

import numpy as np
import trimesh
import time
import random
import Rhino.Geometry as rg

def fractal_dimension(mesh, scaling_factor, max_scale_coef, min_scale, resolution):

    geometry = []

    random.seed(42)

    # -----------------------------
    # Convert Rhino mesh to Trimesh
    # -----------------------------
    def RhinoToLists(mesh, ftype=np.float32, itype=np.int32):
        points = np.zeros((mesh.Vertices.Count, 3), dtype=ftype)
        for idx, v in enumerate(mesh.Vertices):
            points[idx, :] = [v.X, v.Y, v.Z]

        faces = np.zeros((mesh.Faces.Count, 3), dtype=itype)
        for idx, f in enumerate(mesh.Faces):
            faces[idx, :] = [f.A, f.B, f.C]

        return points, faces

    def RhinoToTrimesh(mesh):
        points, faces = RhinoToLists(mesh)
        return trimesh.Trimesh(points, faces)

    m = RhinoToTrimesh(mesh)

    # -----------------------------
    # Box-counting function with random offsets
    # -----------------------------
    def count_boxes_random(voxel_matrix, k, trials=10):
        sx, sy, sz = voxel_matrix.shape
        counts = []
        
        for _ in range(trials):
            # get random grid offsets
            ox = random.randint(0, k - 1)
            oy = random.randint(0, k - 1)
            oz = random.randint(0, k - 1)
            
            # compute padding
            pad_x_left = k - ox
            pad_y_left = k - oy
            pad_z_left = k - oz
            
            pad_x_right = k - ((sx + pad_x_left) % k)
            pad_y_right = k - ((sy + pad_y_left) % k)
            pad_z_right = k - ((sz + pad_z_left) % k)
            
            # apply padding
            padded = np.pad(voxel_matrix,
                         ((pad_x_left, pad_x_right),
                          (pad_y_left, pad_y_right),
                          (pad_z_left, pad_z_right)),
                         mode='constant',
                         constant_values=0)
            
            # reshape into blocks
            sx2, sy2, sz2 = padded.shape
            reshaped = padded.reshape(sx2//k, k,
                              sy2//k, k,
                              sz2//k, k)
            
            pooled = reshaped.max(axis=(1,3,5))
            counts.append(np.sum(pooled))
        
        return np.mean(counts)

    # -----------------------------
    # Fractal dimension calculation
    # -----------------------------
    t0 = time.time()

    # Largest scale: cube that fits the entire structure
    bbox = m.bounds  # (min corner, max corner)
    bbox_sizes = bbox[1] - bbox[0]
    largest_scale = np.max(bbox_sizes)
    if max_scale_coef:
        largest_scale = largest_scale * max_scale_coef

    # Generate scales: divide by scaling factor each time until min_scale

    max_k = int(largest_scale / resolution)

    k_values = []
    k = max_k
    while k >= 1:
        k_values.append(k)
        k = int(k / scaling_factor)

    k_values = np.unique(k_values)
    scales = k_values * resolution

    mask = scales >= min_scale
    scales = scales[mask]
    k_values = k_values[mask]

    # Voxelize at given resolution
    voxels = m.voxelized(resolution)
    voxel_matrix = voxels.matrix

    # Count boxes at each scale using random offsets
    counts = np.array([count_boxes_random(voxel_matrix, k, trials=10)
                       for k in k_values])

    # Only keep scales with >1 box
    mask = counts > 1
    scales = scales[mask]
    counts = counts[mask]

    # Log-log regression for fractal dimension
    log_e = np.log(1 / scales)          # log(1/box_size)
    log_c = np.log(counts)              # log(N)

    fractal_dim, intercept = np.polyfit(log_e, log_c, 1)

    # Store geometry for plotting (optional)
    points = [rg.Point3d(e, c, 0) for e, c in zip(log_e, log_c)]
    geometry.extend(points)

    def f(x):
        y = fractal_dim * x + intercept 
        return y

    line_pts = [rg.Point3d(x, f(x), 0) for x in log_e]
    geometry.append(rg.Polyline(line_pts))

    t_elapsed = time.time() - t0
    print(f'Fractal dimension: {fractal_dim:.4f}')
    print(f'Elapsed time: {t_elapsed:.3f} sec')
    print(f'Scales used: {scales}')

    from sklearn.metrics import r2_score
    x_values = log_e
    y_values = log_c
    predict = np.poly1d([fractal_dim, intercept])
    predicted_y_values = predict(x_values)
    R2 = r2_score(y_values, predicted_y_values)
    print(f"R-squared value: {R2}")

    print(messages)
    
    return fractal_dim, scales, counts

