import gmsh
import numpy as np

class MotorModel:
    def __init__(self, meshSize=0.03, StatorR=0.3, RotorR = 0.247,
                 AirW=0.003, InnerR=0.15, ThetaR=90.0,  poles=2, phases=3,
                 statorWireR=0.275, statorWirer=0.02, rotorWireR=0.2, rotorWirer=0.025,
                 modelName="motor.msh", printing=False):
        """
        Create the model of an electrically excited synchronous motor.

        :param meshSize: The relative size of each element. Default: 0.03
        :type meshSize: float
        :param StatorR: The radius of the stator. Unit: meter. Default: 0.3
        :type StatorR: float
        :param RotorR: The radius of the rotor. Unit: meter. Default: 0.247
        :type RotorR: float
        :param AirW: The width of the air gap. Unit: meter. Default: 0.003
        :type AirW: float
        :param InnerR: The radius of the air at the center. Unit: meter. Default: 0.15
        :rtype InnerR: float
        :param ThetaR: The angle of the rotor clockwise with respect to horizontal. Unit: degree. Default: 90.0
        :type ThetaR: float
        :param poles: The number of poles. Unit: none. Default: 2
        :type poles: int
        :param phases: The number of phases. Unit: none. Default: 3
        :type phases: int
        :param statorWireR: The distance from the motor center to the stator wire center. Unit: meter. Default: 0.275
        :type statorWireR: float
        :param statorWirer: The radius of the stator wire. Unit: meter. Default: 0.02
        :type statorWirer: float
        :param rotorWireR: The distance from the motor center to the rotor wire center. Unit: meter. Default: 0.2
        :type rotorWireR: float
        :param rotorWirer: The radius of the rotor wire. Unit: meter. Default: 0.025
        :type rotorWirer: float
        :param modelName: The name of the model file to be saved. The filename extension should be ".msh".
                          Default: "motor.msh"
        :type modelName: str
        """
        Model = gmsh.model
        Geo = Model.geo
        gmsh.initialize()
        gmsh.model.add("EESM_motor")
        if not printing:
            gmsh.option.setNumber("General.Terminal", 0)

        self.mesh_size = meshSize

        self.StatorRO = StatorR
        self.RotorR = RotorR
        self.AirW = AirW
        self.StatorRI = self.RotorR + self.AirW
        self.InnerR = InnerR

        # The angle of magnetic field
        self.rTheta = ThetaR

        OAng = np.deg2rad(self.rTheta - 90.0)
        IAng = OAng
        printing = False

        Ori = Geo.add_point(0, 0 ,0, self.mesh_size)
        # Create the outer circle of the stator
        Rso_left = Geo.add_point(-self.StatorRO * np.cos(OAng), -self.StatorRO * np.sin(OAng), 0, self.mesh_size)
        Rso_right = Geo.add_point(self.StatorRO * np.cos(OAng), self.StatorRO * np.sin(OAng), 0, self.mesh_size)
        Cso_up = Geo.add_circle_arc(Rso_right, Ori, Rso_left)
        Cso_dn = Geo.add_circle_arc(Rso_left, Ori, Rso_right)
        Curve_so = Geo.add_curve_loop([Cso_up, Cso_dn])
        if printing:
            print("Stator outer done")

        # Create the inner circle of the stator
        Rsi_left = Geo.add_point(-self.StatorRI * np.cos(OAng), -self.StatorRI * np.sin(OAng), 0, self.mesh_size / 2.5)
        Rsi_right = Geo.add_point(self.StatorRI * np.cos(OAng), self.StatorRI * np.sin(OAng), 0, self.mesh_size / 2.5)
        Csi_up = Geo.add_circle_arc(Rsi_right, Ori, Rsi_left)
        Csi_dn = Geo.add_circle_arc(Rsi_left, Ori, Rsi_right)
        Curve_si = Geo.add_curve_loop([Csi_up, Csi_dn])
        if printing:
            print("Stator Inner done")

        # Create the circle of the rotor
        Rr_left = Geo.add_point(-self.RotorR * np.cos(IAng), -self.RotorR * np.sin(IAng), 0, self.mesh_size / 2.5)
        Rr_right = Geo.add_point(self.RotorR * np.cos(IAng), self.RotorR * np.sin(IAng), 0, self.mesh_size / 2.5)
        Cr_up = Geo.add_circle_arc(Rr_right, Ori, Rr_left)
        Cr_dn = Geo.add_circle_arc(Rr_left, Ori, Rr_right)
        Curve_r = Geo.add_curve_loop([Cr_up, Cr_dn])
        if printing:
            print("Rotor outer done")

        # Create the inner hole
        Ir_left = Geo.add_point(-self.InnerR * np.cos(OAng), -self.InnerR * np.sin(OAng), 0, self.mesh_size)
        Ir_right = Geo.add_point(self.InnerR * np.cos(OAng), self.InnerR * np.sin(OAng), 0, self.mesh_size)
        Ci_up = Geo.add_circle_arc(Ir_right, Ori, Ir_left)
        Ci_dn = Geo.add_circle_arc(Ir_left, Ori, Ir_right)
        Curve_i = Geo.add_curve_loop([Ci_up, Ci_dn])
        if printing:
            print("Rotor inner done")

        self.Poles = poles
        self.Phases = phases

        self.RC_wireS = statorWireR
        self.r_wireS = statorWirer

        # Create the stator wires
        Center_wireS = []
        Sides_wireS = []
        Circle_wireS = []
        Curve_wireS = []
        Plane_wireS = []
        Phy_wires = []
        for i in range(self.Poles):
            Phi = np.deg2rad(i * (360 / self.Poles))
            for j in range(self.Phases):
                # The angle of the  slot
                Theta = np.deg2rad(j * (360 / self.Phases))
                c = self.RC_wireS * np.array([np.cos(Theta + Phi), np.sin(Theta + Phi)])
                r = self.r_wireS * np.array([-np.sin(Theta + Phi), np.cos(Theta + Phi)])

                # Add the centers and the points on circles.
                Center_wireS.append(Geo.add_point(c[0], c[1], 0, self.mesh_size / 2))
                Sides_wireS.append([])
                Sides_wireS[-1].append(Geo.add_point((c - r)[0], (c - r)[1], 0, self.mesh_size / 2))
                Sides_wireS[-1].append(Geo.add_point((c + r)[0], (c + r)[1], 0, self.mesh_size / 2))

                # Add circles
                Circle_wireS.append([])
                Circle_wireS[-1].append(Geo.add_circle_arc(Sides_wireS[-1][0], Center_wireS[-1], Sides_wireS[-1][1]))
                Circle_wireS[-1].append(Geo.add_circle_arc(Sides_wireS[-1][1], Center_wireS[-1], Sides_wireS[-1][0]))
                # print(Circle_wireS[i])

                # Add curves
                Curve_wireS.append(Geo.add_curve_loop([Circle_wireS[-1][0], Circle_wireS[-1][1]]))
                # print(Curve_wireS)

                # Add plane
                Plane_wireS.append(Geo.add_plane_surface([Curve_wireS[-1]]))
                # print(Plane_wireS)

                Phase = str(j)
                Pole = "+" if i % 2 == 0 else "-"
                Geo.synchronize()
                Phy_wires.append(Model.add_physical_group(2, [Plane_wireS[-1]], -1, "WireS" + Phase + Pole + str(i)))
        if printing:
            print("Stator wires done")

        # The rotor wire
        self.rTheta = np.deg2rad(self.rTheta)
        self.R_rWire = rotorWireR
        self.r_wire = rotorWirer

        M = self.R_rWire * np.array([-np.sin(self.rTheta), np.cos(self.rTheta)])
        r = self.r_wire * np.array([np.cos(self.rTheta), np.sin(self.rTheta)])

        Points_WireRU = [Geo.add_point(M[0], M[1], 0, self.mesh_size / 2),
                         Geo.add_point((M - r)[0], (M - r)[1], 0, self.mesh_size / 2),
                         Geo.add_point((M + r)[0], (M + r)[1], 0, self.mesh_size / 2)]

        Circle_wireRU = [Geo.add_circle_arc(Points_WireRU[1], Points_WireRU[0], Points_WireRU[2]),
                         Geo.add_circle_arc(Points_WireRU[2], Points_WireRU[0], Points_WireRU[1])]

        Points_WireRD = [Geo.add_point(-M[0], -M[1], 0, self.mesh_size / 2),
                         Geo.add_point((-M + r)[0], (-M + r)[1], 0, self.mesh_size / 2),
                         Geo.add_point((-M - r)[0], (-M - r)[1], 0, self.mesh_size / 2)]

        Circle_wireRD = [Geo.add_circle_arc(Points_WireRD[2], Points_WireRD[0], Points_WireRD[1]),
                         Geo.add_circle_arc(Points_WireRD[1], Points_WireRD[0], Points_WireRD[2])]

        Curve_wireRU = Geo.add_curve_loop(Circle_wireRU)
        Curve_wireRD = Geo.add_curve_loop(Circle_wireRD)
        Plane_RU = Geo.add_plane_surface([Curve_wireRU])
        Plane_RD = Geo.add_plane_surface([Curve_wireRD])
        if printing:
            print("Rotor wires done")

        #Other plane surfaces
        Plane_stator = Geo.add_plane_surface([Curve_so, Curve_si] + Curve_wireS)
        Plane_air = Geo.add_plane_surface([Curve_si, Curve_r])
        Plane_rotor = Geo.add_plane_surface([Curve_r, Curve_wireRU, Curve_wireRD, Curve_i])#
        Plane_inner = Geo.add_plane_surface([Curve_i])
        if printing:
            print("Plane surface done")

        Geo.synchronize()
        if printing:
            print("Synchronization done")

        Model.mesh.generate(2)
        if printing:
            print("Mesh generation done")

        Phy_wirer = [Model.add_physical_group(2, [Plane_RU], -1, "WireR+"),
                     Model.add_physical_group(2, [Plane_RD], -1, "WireR-")]
        Phy_Boundary = Model.add_physical_group(1, [Cso_up, Cso_dn], -1, "Boundary")
        Stator = Model.add_physical_group(2, [Plane_stator], -1, "Stator")
        Air = Model.add_physical_group(2, [Plane_air], -1, "Air")#
        Rotor = Model.add_physical_group(2, [Plane_rotor], -1, "Rotor")
        IAir = Model.add_physical_group(2, [Plane_inner], -1, "IAir")
        if printing:
            print("Physical groups done")

        self.modelName = modelName
        gmsh.write(self.modelName)
        gmsh.finalize()

    def __str__(self):
        return self.modelName

if __name__ == "__main__":
    M = MotorModel()