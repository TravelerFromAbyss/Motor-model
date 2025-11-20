"""
Show three phases of the motor at θ = 0, 60, and 120 degree.
"""
import time
import numpy as np
import FEM_SciPy_gmsh as FEM_SciPy
import os
import matplotlib.pyplot as plt
from Visualization_motor import Visualization
from Motor_model import MotorModel
from readMesh_mod import ReadMesh

def meshMap(ky, phi, SI_amp, RI):
    """
    The meshMap function as the parameter of Mesh class.
    :param ky: The region name.
    :type ky: str
    :param phi: The phase angle of stator current. Unit: rad.
    :type phi: float
    :return: A dict containing the permittivity and current density.
    :rtype: dict
    """
    eleInfo = {"Stator": (4000.0, 0.0), "Rotor": (4000.0, 0.0), "Air": (1.0, 0.0), "IAir": (1.0, 0.0)}
    if ky[0:4] != "Wire":
        return {"mu_r": eleInfo[ky][0], "J": eleInfo[ky][1]}
    else:
        if ky[4] == "S":
            Phase = int(ky[5])
            Pole = 1 if ky[6] == "+" else -1
            return {"mu_r": 1.0, "J": Pole * SI_amp * np.sin(np.deg2rad(phi) - Phase * 2 * np.pi / 3)}
        elif ky[4] == "R":
            Pole = 1 if ky[5] == "+" else -1
            return {"mu_r": 1.0, "J": Pole * RI}
        return {"mu_r": 1.0, "J": 0.0}


def color(ky):
    """
    The color of each region and boundary as the parameter Visualization.mesh_nodes_triangulation().

    :param ky: The region name or "Boundary".
    :return: The color.
    :rtype: str
    """
    colInfo = {"Stator": "#55ff55", "Rotor": "#55ff55", "Air": "#555555", "IAir": "#555555", "Boundary": "cyan"}
    if ky[0:4] != "Wire":
        return colInfo[ky]
    else:
        if ky[4] == "S":
            return "red" if ky[6] == "+" else "blue"
        elif ky[4] == "R":
            return "red" if ky[5] == "+" else "blue"
        return "None"


def airField(M, phi):
    """
    Draw the quiver at the air gap.
    :param M: The motor model.
    :type M: MotorModel
    :param phi: the phase angle of stator current. Unit: rad.
    :type phi: float
    """
    grid_num = (50, 50)
    X = []; Y = []; U = []; V = []; value = []
    for the in np.linspace(0, 2 * np.pi, 100):
        X.append(M.R_rotor * np.cos(the)); Y.append(M.R_rotor * np.sin(the))
        X.append(M.R_rotor * np.cos(the)); Y.append(M.R_rotor * np.sin(the))
        B = M.get_B(X[-1], Y[-1])
        # (By·cosθ - Bx·sinθ)(cos(pi/2 + θ), sin(pi/2 + θ))
        U.append((B[1] * np.cos(the) - B[0] * np.sin(the)) * np.cos(np.pi / 2 + the))
        V.append((B[1] * np.cos(the) - B[0] * np.sin(the)) * np.sin(np.pi / 2 + the))
        # (Bx·cosθ + By·sinθ)(cos(θ), sin(θ))
        U.append((B[0] * np.cos(the) + B[1] * np.sin(the)) * np.cos(the))
        V.append((B[0] * np.cos(the) + B[1] * np.sin(the)) * np.sin(the))
        value.append(B[2])
        value.append(B[2])

    fig, ax = plt.subplots()
    Q = ax.quiver(X, Y, U, V, value, cmap="jet", scale=13)
    text = ax.text(0.95, 0.95, f"$\\varphi$={phi}\n$\\tau$={M.torque(0.3):.4e}", transform=ax.transAxes, ha='right',
                   va='top', fontsize=10, color='black')
    fig.colorbar(Q, ax=ax, label="Magnetic flux density/$T$")  # Add a colorbar
    ax.set_title("Magnetic Field")
    ax.set_xlabel("x/m")
    ax.set_ylabel("y/m")

    ax.set_xlim(-0.5, 0.5)
    ax.set_ylim(-0.5, 0.5)
    ax.set_aspect("equal")
    # plt.get_current_fig_manager().window.showMaximized()
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    statorWirer = 0.04
    rotorWirer = 0.05

    SI_amp = 800 / (np.pi * statorWirer**2)
    RI = 0 # 200 / (np.pi * rotorWirer**2)

    # when θ = 0 degree
    theta = 0
    model1 = MotorModel(ThetaR=theta + 90, RotorR=0.247, AirW=0.003, StatorR=0.30, InnerR=0.15,
                        statorWireR=0.275, statorWirer=0.02, rotorWireR=0.2, rotorWirer=0.025,
                        modelName="Three_phases.msh")

    direction = os.getcwd()
    mesh_filename = str(model1)
    print("Mesh Done")

    start_time = time.time()
    M1 = FEM_SciPy.Mesh(ReadMesh(direction, mesh_filename), meshMap, [theta, SI_amp, RI], name = "M1")
    end_time = time.time()  # Record the end time
    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time:.4f} seconds")

    # when θ = 60 degree
    theta = 60
    model2 = MotorModel(ThetaR=theta + 90, RotorR=0.247, AirW=0.003, StatorR=0.30, InnerR=0.15,
                        statorWireR=0.275, statorWirer=0.02, rotorWireR=0.2, rotorWirer=0.025,
                        modelName="Three_phases.msh")
    direction = os.getcwd()
    mesh_filename = str(model2)
    print("Mesh Done")

    start_time = time.time()
    M2 = FEM_SciPy.Mesh(ReadMesh(direction, mesh_filename), meshMap, [theta, SI_amp, RI], name = "M2")
    end_time = time.time()  # Record the end time
    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time:.4f} seconds")

    # when θ = 120 degree
    theta = 120
    model3 = MotorModel(ThetaR=theta + 90, RotorR=0.247, AirW=0.003, StatorR=0.30, InnerR=0.15,
                        statorWireR=0.275, statorWirer=0.02, rotorWireR=0.2, rotorWirer=0.025,
                        modelName="Three_phases.msh")
    direction = os.getcwd()
    mesh_filename = str(model3)
    print("Mesh Done")

    start_time = time.time()
    M3 = FEM_SciPy.Mesh(ReadMesh(direction, mesh_filename), meshMap, [theta, SI_amp, RI], name = "M3")
    end_time = time.time()  # Record the end time
    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time:.4f} seconds")

    Vis1 = Visualization(M1)
    Vis2 = Visualization(M2)
    Vis3 = Visualization(M3)
    print(M1, M2, M3)

    print(M1.torque(0.3))
    print(M2.torque(0.3))
    print(M3.torque(0.3))

    """
    airField(M1, 0)
    airField(M2, 60)
    airField(M3, 120)
    Vis1.mesh_colormap2(shading="gouraud")
    Vis2.mesh_colormap2(shading="gouraud")
    Vis3.mesh_colormap2(shading="gouraud")
    """

    Vis1.mesh_nodes_triangulation(color, linewidth=0.5, grid=False)
    Vis1.mesh_field_quiver(False, grid_num=(60, 60), scale = 40, vrange=(0, 5))
    Vis2.mesh_field_quiver(False, grid_num=(60, 60), scale = 40, vrange=(0, 5))
    Vis3.mesh_field_quiver(False, grid_num=(60, 60), scale = 40, vrange=(0, 5))
