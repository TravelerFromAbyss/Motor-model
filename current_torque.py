"""
Calculate the torque according to the change current.
The calculation is performed at several rotor angle (each 10 between 0 and 180 degree).
"""
import time
import numpy as np
import FEM_SciPy_gmsh as FEM_SciPy
import os
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as m_colors
from Visualization_motor import Visualization
from Motor_model import MotorModel
from readMesh_mod import ReadMesh
import matplotlib.font_manager as fm


plt.rcParams['font.family'] = "Arial"
plt.rcParams['font.size'] = 20

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

Theta = range(0, 180, 10) # The simple angles are taken
cmap = plt.get_cmap("jet")
norm = m_colors.Normalize(vmin=min(Theta), vmax=max(Theta))

fitCoe = [np.array([]), np.array([]), np.array([])]

for theta in Theta:
    model1 = MotorModel(ThetaR=theta + 90, RotorR=0.247, AirW=0.003, StatorR=0.30, InnerR=0.15,
                        statorWireR=0.275, statorWirer=0.02, rotorWireR=0.2, rotorWirer=0.025,
                        modelName="Torque_calculation.msh")
    direction = os.getcwd()
    mesh_filename = str(model1)
    print(f"{theta}: Mesh Done")

    tau_SI = []
    I = np.linspace(0, 60, 16) # np.array([2])

    for i in I:
        SI_amp = i * 1e2 / (np.pi * statorWirer**2) # 8e2
        RI = 2e2 / (np.pi * rotorWirer**2) # 2e2

        start_time = time.time()
        M = FEM_SciPy.Mesh(ReadMesh(direction, mesh_filename), meshMap, [theta, SI_amp, RI], name="M")
        end_time = time.time()  # Record the end time
        tau_SI.append(M.torque(0.3))
        elapsed_time = end_time - start_time
        print(f"{len(tau_SI)}: I = {RI:.3e}, torque = {tau_SI[-1]:.3e};\tElapsed time: {elapsed_time:.4f} seconds")

        """
        This code was used to 
        if i == 2:
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
            cf = ax.quiver(X, Y, U, V, value, cmap="jet", scale=10)

            fig.colorbar(cf, ax=ax, label="Magnetic flux density/$T$")  # Add a colorbar
            ax.set_title("Magnetic Field")
            ax.set_xlabel("x/m")
            ax.set_ylabel("y/m")
            ax.set_aspect("equal")
            plt.grid(False)

            plt.show()
        """

    c = np.polyfit(I * 1e2, tau_SI, 1)
    fitCoe[0] = np.append(fitCoe[0], c[0])
    fitCoe[1] = np.append(fitCoe[1], c[1])
    print(c)
    plt.plot(I * 1e2, tau_SI, marker='o', label=f"$\\varphi$={theta}", color=cmap(norm(theta)))


plt.xlabel("Amplitude of the stator current / $\\mathrm{A}$") # ("The rotor current / $\\mathrm{A}$") #
plt.ylabel("Torque / $\\mathrm{N\\cdot m}$")
# plt.title("The relationship between the torque and the rotor current") # ("The relationship between the torque and \nthe amplitude of the stator current") #
plt.legend(loc="upper left", prop=fm.FontProperties(size=10))
plt.show()

plt.plot(Theta, fitCoe[0] * 5, label="$a\\times 5$") # r5 s50
plt.plot(Theta, fitCoe[1] * 1e-1, label="$b\\times 10^{-1}$") # r-10 s-10
plt.xlabel("$\\varphi$ / $\\mathrm{deg}$")
plt.ylabel("Fitting coefficients")
# plt.title("Fitting coefficients for $\\tau$-$I_r$ relationship")
plt.legend(prop=fm.FontProperties(size=12))
plt.show()

