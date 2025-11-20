"""
The final version of finite element method algorithm`
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse, integrate
from matplotlib.tri import Triangulation
from readMesh_mod import ReadMesh

class Element:
    def __init__(self, nodes, points, region, J=0.0, mu_r=1.0, mu0=4e-7 * np.pi):
        """
        This class will store and process each single element in the mesh.

        :param nodes: The x and y coordinates of the all nodes.
        :param points: The id of the three angles of this element.
        :param region: The region name which this node is in.
        :param J: The current density in this node. Unit: Ampere·meter^(-2). Default 0.0.
        :param mu_r: The relative permittivity in this node: Unit: none. Default 1.0.
        :param mu0: The reference permittivity of vacuum. Unit: Faraday/meter. Default 4e-7 π.
        """
        self.NumberOfNodes = len(nodes)
        self.nodes = nodes
        self.points = np.array(points)
        self.region = region

        self.p1, self.p2, self.p3 = points
        self.com = ((nodes[self.p1][0] + nodes[self.p2][0] + nodes[self.p3][0]) / 3,
                    (nodes[self.p1][1] + nodes[self.p2][1] + nodes[self.p3][1]) / 3)
        # The coordinates of the center of the element.

        self.J = J
        self.mu0 = mu0
        self.mu_r = mu_r
        self.mu = self.mu0 * self.mu_r
        self.Area = 0.5 * abs(nodes[self.p1][0] * (nodes[self.p2][1] - nodes[self.p3][1]) +
                              nodes[self.p2][0] * (nodes[self.p3][1] - nodes[self.p1][1]) +
                              nodes[self.p3][0] * (nodes[self.p1][1] - nodes[self.p2][1]))
        self.a = [nodes[self.p2][0] * nodes[self.p3][1] - nodes[self.p3][0] * nodes[self.p2][1],
                  nodes[self.p3][0] * nodes[self.p1][1] - nodes[self.p1][0] * nodes[self.p3][1],
                  nodes[self.p1][0] * nodes[self.p2][1] - nodes[self.p2][0] * nodes[self.p1][1]]
        self.b = [nodes[self.p2][1] - nodes[self.p3][1],
                  nodes[self.p3][1] - nodes[self.p1][1],
                  nodes[self.p1][1] - nodes[self.p2][1]]
        self.c = [nodes[self.p3][0] - nodes[self.p2][0],
                  nodes[self.p1][0] - nodes[self.p3][0],
                  nodes[self.p2][0] - nodes[self.p1][0]]

        self.S = [] # sparse.dok_matrix((self.NumberOfNodes, self.NumberOfNodes))
        self.set_S()
        # the stiffness matrix shown in article
        self.T = [] # sparse.dok_matrix((self.NumberOfNodes, 1))
        self.set_T()
        # the load vector shown in article

        self.Bx = None
        self.By = None
        self.B = None

    def set_S(self):
        """
        Calculate the stiffness matrix and set it to self.S

        :return: the calculated S
        :rtype: list
        """
        self.S.append([self.p1, self.p1,
                       (self.b[0] * self.b[0] + self.c[0] * self.c[0]) / (4 * self.Area * self.mu)])
        self.S.append([self.p1, self.p2,
                       (self.b[0] * self.b[1] + self.c[1] * self.c[0]) / (4 * self.Area * self.mu)])
        self.S.append([self.p1, self.p3,
                       (self.b[0] * self.b[2] + self.c[2] * self.c[0]) / (4 * self.Area * self.mu)])
        self.S.append([self.p2, self.p2,
                       (self.b[1] * self.b[1] + self.c[1] * self.c[1]) / (4 * self.Area * self.mu)])
        self.S.append([self.p2, self.p3,
                       (self.b[1] * self.b[2] + self.c[2] * self.c[1]) / (4 * self.Area * self.mu)])
        self.S.append([self.p3, self.p3,
                       (self.b[2] * self.b[2] + self.c[2] * self.c[2]) / (4 * self.Area * self.mu)])
        return self.S

    def set_T(self):
        """
        Calculate the load vector and set it to self.T.
        :return: the calculated T.
        """
        self.T.append([self.p1, self.J / 3 * self.Area])
        self.T.append([self.p2, self.J / 3 * self.Area])
        self.T.append([self.p3, self.J / 3 * self.Area])
        return self.T

    def set_B(self, A):
        """
        Calculate the magnetic field in the element by taking the gradient of the magnetic potential A.
        :param A: a tuple of three number. The magnetic potential at the three nodes.
                  Unit: Volt·s·meter^(-1).
        :return: a tuple of the magnetic field magnitude in this element, its x-component, and its y-component.
                 Unit: Tesla.
        """
        # B = (ci*Bi + cj*Bj + ck*Bk, -bi*Ai - bj*Aj - bk*Ak)/(2*Area)
        self.Bx = 1 / (2 * self.Area) * (self.c[0] * A[0] + self.c[1] * A[1] + self.c[2] * A[2])
        self.By = -1 / (2 * self.Area) * (self.b[0] * A[0] + self.b[1] * A[1] + self.b[2] * A[2])
        self.B = np.sqrt(self.Bx**2 + self.By**2)
        return self.B, self.Bx, self.By

class Mesh:
    def __init__(self, mesh: ReadMesh, meshMap, mapPara, A0=0, printing=False, name="None", mu0=4e-7 * np.pi, N=50):
        """
        This class stores all the information of the mesh and derive the magnetic potential A and field B.

        :param mesh: A ReadMesh object.
        :type mesh: ReadMesh
        :param meshMap: a function that return the values of a certain type of region.
            The function should have parameters in form of:
                key: str. The name of the region type.
                phi: float. the phase angle of the stator current. Unit: rad.
                SI_amp: float. the stator current amplitude. Unit: Ampere.
                RI: float. the rotor current. Unit: Ampere.
            and return type of:
                a dict containing two entries {"mu_r": float, "J": float}:
                    "mu_r": The relative permittivity. Unit: none.
                    "J": The current density. Unit: Ampere·meter^(-2).
        :type meshMap: Callable[[str, float, float, float]] -> dict
        :param mapPara: The list of phi, SI_amp, and RI used in meshMap.
        :type mapPara: iterable of three elements
        :param A0: The reference magnetic potential at the edge of the mesh. Unit: Volt·s·meter^(-1).
        :type A0: float
        :param printing: Whether the running stage of processing should be printed.
        :type printing: bool
        :param name: The name of the Mesh.
        :type name: str
        :param mu0: The reference permittivity of vacuum. Unit: Faraday/meter. Default 4e-7 π.
        :type mu0: float
        :param N: The number of coils in the solenoid. Unit: none. Default 50.
        :type N: int
        """
        self.name = name
        self.mu_0 = mu0
        self.N = N

        self.mesh = mesh
        mesh.get_mesh()

        bd = mesh.get_boundary_info() # Get the boundary nodes of the mesh
        rg = mesh.get_region_info() # Get which elements are contained in each type of region.

        self.nodes = np.array(list(zip(mesh.x_coord, mesh.y_coord)))

        self.elements = []
        self.rg_elements = {}
        self.R_rotor = 0
        # Create the Element object for each element in the mesh.
        for ky in rg.keys():
            self.rg_elements[ky] = []
            for ele in rg[ky]:
                # if ky == "IAir":
                #     continue
                e = Element(self.nodes, ele, ky, mu_r=meshMap(ky, *mapPara)["mu_r"], J=meshMap(ky, *mapPara)["J"])
                self.rg_elements[ky].append(e)
                self.elements.append(e)
                if ky == "Air":
                    self.R_rotor += np.sqrt(e.com[0]**2 + e.com[1]**2)

        self.R_rotor /= len(rg["Air"]) if "Air" in rg.keys() else 1
        # The radius of rotor is taken as the average distance from all air gap element to the center
        # It will be used in torque calculation.
        self.elements = np.array(self.elements)
        if printing:
            print("Elements done")

        self.boundary = np.array(bd["Boundary"])

        self.NumberOfNodes = len(self.nodes)
        self.NumberOfElements = len(self.elements)

        boundaryNodes = self.boundary # Extracting the nodes on the boundary

        # Assembling S. The rows and columns of nodes in boundaryNodes are set to 0,
        # with the diagonal elements_info to 1
        self.S = sparse.dok_matrix((self.NumberOfNodes, self.NumberOfNodes))
        for e in self.elements:
            for s in e.S:
                if s[0] in boundaryNodes or s[1] in boundaryNodes:
                    self.S[s[0], s[1]] = 0 if s[0] != s[1] else 1
                    continue
                self.S[s[0], s[1]] += s[2]; self.S[s[1], s[0]] += s[2] if s[0] != s[1] else 0
        if printing:
            print("S summed")

        # Assembling T. The elements_info of the nodes in boundaryNodes are set to A0
        self.T = sparse.dok_matrix((self.NumberOfNodes, 1))
        for e in self.elements:
            for t in e.T:
                if t[0] in boundaryNodes:
                    self.T[t[0]] = A0
                    continue
                self.T[t[0], 0] += t[1]
        if printing:
            print("T summed")

        # Solve the equation SA=T
        self.A = self.N * sparse.linalg.spsolve(sparse.csr_array(self.S), sparse.csr_array(self.T))
        if printing:
            print("Solving done")

        self.B = np.zeros((self.NumberOfElements, 3))
        self.set_B()
        if printing:
            print("B set")

        self.ele, self.com = self.get_Elements()
        self.tri_nodes = Triangulation(self.nodes[:, 0], self.nodes[:, 1], self.ele)
        self.tri_elements = Triangulation(self.com[:, 0], self.com[:, 1])
        self.triFinder = self.tri_nodes.get_trifinder()

    def set_B(self):
        """
        Set the magnetic field of the elements.
        :return: The calculated array of B
        :rtype: numpy array
        """
        for i in range(self.NumberOfElements):
            e: Element = self.elements[i]
            e.set_B((self.A[e.p1], self.A[e.p2], self.A[e.p3]))
            self.B[i] = np.array([e.Bx, e.By, e.B])
        return self.B

    def print_S(self):
        """
        Print the calculated stiffness matrix of the elements.
        :return: None
        """
        for i in range(self.NumberOfNodes):
            for j in range(self.NumberOfNodes):
                print("%.4f\t" % float(self.S[i, j]), end="")
            print()

    def print_T(self):
        """
        Print the calculated load vector of the elements.
        :return: None
        """
        for i in range(self.NumberOfNodes):
            print("%.4f\t" % float(self.T[i, 0]))

    def get_Elements(self):
        """
        Return the id and the COM of all the elements.
        :return: The id and the COM of all the elements.
        :rtype: tuple of two iterables: (id, COM coordinates)
        """
        ele = []
        com = []
        for e in self.elements:
            ele.append(e.points)
            com.append(e.com)
        return np.array(ele), np.array(com)

    def get_Elements_Area(self):
        """
        Return the area of all the elements.
        :return: The area of all the elements.
        :rtype: list
        """
        area = []
        for e in self.elements:
            area.append(e.Area)
        return area

    def get_Elements_J(self):
        """
        Return the current density of all the elements.
        :return: The current density of all the elements.
        :rtype:list
        """
        j = []
        for e in self.elements:
            j.append(e.J)
        return j

    def get_Elements_mu(self):
        """
        Return the relative permittivity of all the elements.
        :return: The relative permittivity of all the elements.
        :rtype: list
        """
        mu = np.array([])
        for e in self.elements:
            mu = np.append(mu, e.mu_r)
        return mu

    def get_B(self, x, y):
        """
        Return the magnetic field at a position of the mesh.
        :param x: The x-coordinate.
        :type x: float
        :param y: The y-coordinate.
        :type y: float
        :return: The magnetic field magnitude, x-component, and y-component.
        :rtype: tuple of three
        """
        eIndex = self.triFinder(x, y)
        if eIndex != -1:
            e: Element = self.elements[eIndex]
            return np.array((e.Bx, e.By, e.B))
        else:
            return np.array((0, 0, 0))

    def torque(self, l_ef, r=-1.0, **kwargs):
        """
        Calculate the torque.
        :param l_ef: the axis-length of the motor. Unit: meter
        :type l_ef: float
        :param r:
        :param kwargs:
        :return: The torque. Unit: Newton·meter
        :rtype: float
        """
        """
        r = self.R_rotor if r == -1.0 else r

        def getB_theta(theta):
            X = r * np.cos(theta); Y = r * np.sin(theta)
            Bx, By, Bm = self.get_B(X, Y)
            Br = Bx * np.cos(theta) + By * np.sin(theta)
            Bt = By * np.cos(theta) - Bx * np.sin(theta)
            # print(f"{np.rad2deg(theta): .2f}: r {Br: .4f}, t {Bt: .4f}, {Br * Bt: .4f}")
            return Br * Bt

        return l_ef * r ** 2 / self.mu_0 * (integrate.quad(getB_theta, 0, 2 * np.pi)[0])
        N = kwargs["N"] if "N" in kwargs.keys() else len(self.mesh.get_region_info()["Air"])

        Ang = np.linspace(0, 2 * np.pi, N)
        pre_tau = 0
        for theta in Ang[1:]:
            X = r * np.cos(theta); Y = r * np.sin(theta)
            Bx, By, Bm = self.get_B(X, Y)
            Br = Bx * np.cos(theta) + By * np.sin(theta)
            Bt = By * np.cos(theta) - Bx * np.sin(theta)
            pre_tau += Br * Bt
        return l_ef * r ** 2 * pre_tau * 2 * np.pi / N / self.mu_0
        """

        pre_tau = 0
        for e in self.rg_elements["Air"]:
            com = e.com
            theta = np.atan(com[1] / com[0]) + (0 if com[1] >= 0 else np.pi)
            r2 = com[0] ** 2 + com[1] ** 2
            # print(com[0], com[1], np.rad2deg(theta))

            Br = e.Bx * np.cos(theta) + e.By * np.sin(theta)
            Bt = e.By * np.cos(theta) - e.Bx * np.sin(theta)
            pre_tau += Br * Bt * r2

        return l_ef * pre_tau * 2 * np.pi / len(self.rg_elements["Air"]) / self.mu_0




    def __str__(self):
        return f"Mesh name: {self.name}, with {self.NumberOfNodes} nodes, {self.NumberOfElements} elements_info"

