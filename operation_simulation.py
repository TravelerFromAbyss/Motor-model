"""
Simulate the rotation of an EESM motor and make an animation of the change in magnetic field in the air gap.
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


def meshMap(ky, phi):
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

model = MotorModel(ThetaR=90, RotorR=0.247, AirW=0.003, StatorR=0.30, InnerR=0.15,
                   statorWireR=0.275, statorWirer=0.02, rotorWireR=0.2, rotorWirer=0.025,
                   modelName="Simulation.msh")
direction = os.getcwd()
mesh_filename = str(model)
print("Mesh Done")

start_time = time.time()
M = FEM_SciPy.Mesh(ReadMesh(direction, mesh_filename), meshMap, [0], name ="M1")
end_time = time.time()  # Record the end time
elapsed_time = end_time - start_time
print(f"Elapsed time: {elapsed_time:.4f} seconds")

# Create the initial frame
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

# The quiver of the firast frame
fig, ax = plt.subplots()
Q = ax.quiver(X, Y, U, V, value, cmap="jet", scale=13)
text = ax.text(0.95, 0.95, f"$\\varphi$={0}\n$\\tau$={M.torque(0.3):.4e}", transform=ax.transAxes, ha='right', va='top', fontsize=10, color='black')
fig.colorbar(Q, ax=ax, label="Magnetic flux density/$T$")  # Add a colorbar
ax.set_title("Magnetic Field")
ax.set_xlabel("x/m")
ax.set_ylabel("y/m")

ax.set_xlim(-0.5, 0.5)
ax.set_ylim(-0.5, 0.5)
ax.set_aspect("equal")
plt.tight_layout()
# plt.get_current_fig_manager().window.showMaximized()
# plt.show()
q = []
tau = []
fs = 72
# Simulate each frame during the rotation and record the vector field.
for i in range(fs):
    print(i)
    model = MotorModel(ThetaR=i * 5 + 90, RotorR=0.247, AirW=0.003, StatorR=0.30, InnerR=0.15,
                       statorWireR=0.275, statorWirer=0.02, rotorWireR=0.2, rotorWirer=0.025,
                       modelName="Simulation.msh")
    print("Mesh Done")

    direction = os.getcwd()
    mesh_filename = str(model)
    # Calculate new U and V components based on the frame number
    M = FEM_SciPy.Mesh(ReadMesh(direction, mesh_filename), meshMap, [i * 5], name="M1")
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

    q.append((U, V, value)) # Update the quiver vectors
    tau.append(M.torque(0.3))

def update(frame):
    print(frame - 1)
    Q.set_UVC(*q[frame - 1]) # Update the quiver vectors
    text.set_text(f"$\\varphi$={frame * 5}\n$\\tau$={tau[frame - 1]:.4e}")
    return Q, # Return the updated Quiver object


plt.tight_layout()
#plt.get_current_fig_manager().window.showMaximized()
ani = FuncAnimation(fig, update, frames=fs, interval=10, blit=False)
ani.save('field in air_rotor.gif', writer='pillow', fps=8)
plt.show()
