"""
The test code to the three classes
A square area with a wire perpendicular to the palen at the center is generated.
"""
import time
import numpy as np
from jinja2.filters import do_title

import FEM_SciPy_gmsh as FEM_SciPy
from Visualization_motor import Visualization
from readMesh_mod import ReadMesh
import os

def meshMap(ky):
    eleInfo = {"Steel": {"mu_r": 4000.0, "J": 0.0}, "Wire": {"mu_r": 4000.0, "J": 10000.0}}
    return eleInfo[ky]


def color(ky):
    colInfo = {"Steel": "#55ff55", "Wire": "#ff5555", "Boundary": "purple", "Inner": "cyan"}
    return colInfo[ky]


if __name__ == "__main__":
    """
    n = 55

    def generate_element(poles, i, j, k):
        J = 0.0
        mu = 100.0
        J_ = 10.0
        mu_ = 1.0
        r = 10
        if (r >= poles[i][0] >= -r and r >= poles[i][1] >= -r and
            r >= poles[j][0] >= -r and r >= poles[j][1] >= -r and
            r >= poles[k][0] >= -r and r >= poles[k][1] >= -r):
            J = J_

        return FEM_SciPy.Element(poles, [i, j, k], J=J, mu_r=mu)


    poles = []
    for i in range(n - 1, -n, -1):
        for j in range(1 - n, n):
            if i % 2 == j % 2:
                poles.append([j, i])

    e = []
    b = []
    for i in range(len(poles)):
        # print(f"{i}: {poles[i]}")
        J = 0.0
        if poles[i][0] % 2 == 0 and poles[i][0] < n - 1 and poles[i][1] > 1 - n:
            e.append(generate_element(poles, i, i + n, i + 1))
            e.append(generate_element(poles, i, i + 2 * n - 1, i + n))
        if poles[i][0] % 2 == 1:
            e.append(generate_element(poles, i, i + n, i - n + 1))
            e.append(generate_element(poles, i, i + n - 1, i + n))

        if (poles[i][0] == n - 1 or poles[i][0] == 1 - n) and poles[i][1] > 1 - n:
            b.append([i, i + 2 * n - 1])
            # print(f"{i} {poles[i]}, {i + 2 * n - 1} {poles[i + 2 * n - 1]}")
        if (poles[i][1] == n - 1 or poles[i][1] == 1 - n) and poles[i][0] < n - 1:
            b.append([i, i + 1])
            # print(f"{i} {poles[i]}, {i + 1} {poles[i + 1]}")

    poles = np.array(poles)
    e = np.array(e)
    b = np.array(b)
    """

    dir = os.getcwd()
    mesh_filename = "square.msh"

    mesh_reader = ReadMesh(dir, mesh_filename)
    mesh_reader.get_mesh()
    bd = mesh_reader.get_boundary_info()
    rg = mesh_reader.get_region_info()

    p = list(zip(mesh_reader.x_coord, mesh_reader.y_coord))
    e = []
    b = bd["Boundary"]
    print(b)

    print("Mesh Done")
    start_time = time.time()
    M1 = FEM_SciPy.Mesh(mesh_reader, meshMap, [], 0, True, name = "M1")
    # M1.print_S()
    # M1.print_T()
    end_time = time.time()  # Record the end time
    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time:.4f} seconds")

    Vis = Visualization(M1)

    print(M1)
    print(M1.get_B(0, 0))
    print(M1.get_B(95, 95))

    Vis.mesh_nodes_triangulation(color, linewidth=0.5, grid=False)
    Vis.mesh_colormap("Area", "Area", unit="$\\mathrm{m^2}$")
    Vis.mesh_colormap("J", "current density", unit="$\\mathrm{J·m^{-2}}$")
    Vis.mesh_colormap("mu", "magnetic permeability",unit="$\\mathrm{H·m^{-1}}$")
    # Vis.mesh_elements_triangulation(linewidth=0.5)
    Vis.mesh_colormap(shading="gouraud")
    Vis.mesh_contourf()
    Vis.mesh_contours()
    Vis.mesh_colormap("B", "magnetic flux density", unit="$\\mathrm{T}$")
    Vis.mesh_colormap2("B", "magnetic flux density", unit="$\\mathrm{T}$", shading="gouraud" )
    Vis.mesh_field_quiver(grid_num=(50, 50), scale = 40)
    Vis.mesh_3d2("A", "magnetic potential", unit="$\\mathrm{T·m}$")
    Vis.mesh_scatter3d("B", "magnetic flux density", unit="$\\mathrm{T}$", s=0.3)
