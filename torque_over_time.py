"""
Calculate the change in torque over time
"""
import time
import numpy as np
import matplotlib.pyplot as plt
import FEM_SciPy_gmsh as FEM_SciPy
import os
from Visualization_motor import Visualization
from Motor_model import MotorModel
from readMesh_mod import ReadMesh
from matplotlib.animation import FuncAnimation


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


statorWirer = 0.04
rotorWirer = 0.05

SI_amp = 800 / (np.pi * statorWirer**2)
RI = 200 / (np.pi * rotorWirer**2)

f = 50
omega = 60 * f / 1

n = 61
Theta = np.linspace(0, 360, n) # (-180, 180, n)
T = Theta / omega

torque = []
for i in range(n):
    model1 = MotorModel(ThetaR=Theta[i] + 90, RotorR=0.247, AirW=0.003, StatorR=0.30, InnerR=0.15,
                        statorWireR=0.275, statorWirer=0.02, rotorWireR=0.2, rotorWirer=0.025,
                        modelName="Over_time.msh")
    direction = os.getcwd()
    mesh_filename = str(model1)
    print("Mesh Done")

    start_time = time.time()
    M = FEM_SciPy.Mesh(ReadMesh(direction, mesh_filename), meshMap, [Theta[i], SI_amp, RI], name="M1") # -90
    end_time = time.time()  # Record the end time
    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time:.4f} seconds")
    torque.append(M.torque(0.3))
    print(f"{i}: φ={Theta[i]}, τ={torque[-1]:.4f}")

    # Vis = Visualization(M)
    # Vis.mesh_field_quiver(False, grid_num=(50, 50), scale = 40)

plt.rcParams['font.family'] = "Arial"
plt.rcParams['font.size'] = 20

plt.plot(Theta, torque, marker="o")
plt.xlabel("$\\varphi$ in $\\mathrm{deg}$") #("Angle displacement between the rotor and stator field / $\\mathrm{deg}$") #
plt.ylabel("Torque / $\\mathrm{N\\cdot m}$")
#plt.title("The relationship between torque and angle displacement") # ("The relationship between torque and angle position") #
plt.show()

plt.plot(T, torque, marker="o")
plt.xlabel("time / $s$")
plt.ylabel("Torque / $\\mathrm{N\\cdot m}$")
#plt.title("The torque over time")
plt.show()
