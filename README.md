This project use finite element method (FEM) to derive the electromagnetic field in an motor.

The functions of main files:
Motor_model: Contains a MotorModel class. It creates the geometry and meshes it using gmsh according to the information given by the user. A .msh is saved.
readMesh_mod: Contains a ReadMeshread class. It need a file directary and file name to read the mesh, and it tells you the information about each region and boundary information.
FEM_SciPy_gmsh: Contains an Element class and a Mesh class. The Mesh class perform FEM on given model. It also calculates the torque produced by the motor.
Visualization_motor: Contains a Visualization class. It receive the Mesh object and preoduce diagram according to the user's needs.
