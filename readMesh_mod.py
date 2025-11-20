"""
A modified version of readMesh.py to be compatible with our model.
"""
import meshio
import numpy as np
import matplotlib.pyplot as plt
import time
import os


class ReadMesh:
    def __init__(self, dir, mesh_file, printing=False):
        """
        Read the .msh file and get the region and boundary information
        :param dir: The directory of the .msh file.
        :type dir: str
        :param mesh_file: The name of the .msh file. It should end with extension name ".msh".
        :type mesh_file: str
        :param printing: Whether to print the physical regions and boundaries. Default: False.
        :type printing: bool
        """
        filename = dir+'\\'+mesh_file
        self.filename = filename       # Complete filename
        self.mesh = None               # The mesh object returned by meshio.read()
        self.x_coord = None            # The x-coordinates of all nodes
        self.y_coord = None            # The y-coordinates of all nodes
        self.elementcon = None         # The nodes for all triangular elements
        self.linecon = None            # The nodes on the boundaries
        self.region = None             # The names of physical planes
        self.boundary = None           # The names of the boundaries
        self.boundary_info = None      # Boundary information dict
        self.region_info = None        # Physical plane information dict
        self.printing = printing

    def get_mesh(self):
        """
        Read the mesh, print existing physical planes and boundaries.
        """
        tic = time.time()    
        # Read the mesh
        if self.printing:
            print('---------------------------------------------------------')
            print('--------------------- Read mesh -------------------------')
            print('---------------------------------------------------------')
        mesh = meshio.read(self.filename)
        if self.printing:
            print("- Mesh file has been read. ")
        
        # Get the node coordinates
        x_coord = mesh.points[:,0] # x-coordinates in meter
        y_coord = mesh.points[:,1] # y-coordinates in meter
        if self.printing:
            print("- Node coordinate has been assembled.")
        
        # Get the element nodes
        elementcon = mesh.cells_dict['triangle']
        if self.printing:
            print("- Element connection has been assembled.")
        
        # Get the boundary node
        linecon = mesh.cells_dict['line']
        if self.printing:
            print("- Line connection has been assembled.\n")
        toc = time.time()-tic
        if self.printing:
            print('Elapsed time is:',toc,'seconds\n')
        
        # Save the result
        self.mesh = mesh
        self.x_coord = x_coord
        self.y_coord = y_coord
        self.elementcon = elementcon
        self.linecon = linecon
        
        # Print the physical plane names and the boundary names
        if self.printing:
            print('---------------------------------------------------------')
            print('--------- List physical regions and boundaries ----------')
            print('---------------------------------------------------------')
        tic = time.time()  
        key = self.mesh.field_data.keys()
        key_list = list(key)
        key_regions = []
        key_boundaries = []
        for key in key_list:
            if self.printing:
                print(key)
            key_name = list(self.mesh.cell_sets_dict[key].keys())[0]     # physical plane or boundary key name, "triangle" and "line"
            if key_name == 'triangle':
                key_regions.append(key)
            elif key_name == 'line':
                key_boundaries.append(key)

        if self.printing:
            print("-- Physical regions:")
            for key in key_regions:
                print('- '+ key)

            print("-- Boundaries:")
            for key in key_boundaries:
                print('- '+ key)
            print("\n")
        toc = time.time()-tic

        if self.printing:
            print('Elapsed time is:',toc,'seconds')
            print("\n")
        
        # Save the result
        self.region = key_regions
        self.boundary = key_boundaries

    def get_boundary_info(self):
        """
        Get the boundary information from the mesh.
        :return: The boundary information
        :rtype: dict
        """
        tic = time.time()
        if self.printing:
            print(self.printing)
            print('---------------------------------------------------------')
            print('---------- Assemble boundary information ----------------')
            print('---------------------------------------------------------')
        boundary_info = {} # Empty dict
        for key in self.boundary:                           # Traverse the boundary list
            subline = self.mesh.cell_sets_dict[key]
            line_id = subline.values()                      # Get the values
            line_id_list = list(line_id)                    # value to list
            linecon_part = self.linecon[line_id_list[0]]    # Read the nodes
            linecon_part_array = linecon_part.reshape(-1,1) # Make the nodes list a vertical list
            node_id= np.unique(linecon_part_array)          # Eliminate the repetitive nodes
            values = node_id
            boundary_info[key] = values                     # Save the nodes list to the dict.
        if self.printing:
            #print(boundary_info)
            print("- Boundary information has been assemblied.")
            print('\n')
        toc = time.time()-tic
        if self.printing:
            print('Elapsed time is:',toc,'seconds')
            print("\n")

        # returns
        self.boundary_info = boundary_info
        return boundary_info

    def get_region_info(self):
        """
        Get the region information from the mesh.
        :return: the region information
        :rtype: dict
        """
        tic = time.time()
        if self.printing:
            print(self.printing)
            print('---------------------------------------------------------')
            print('----------- Assemble region information -----------------')
            print('---------------------------------------------------------')
        region_info = {} # Empty dict
        for key in self.region:                                   # Traverse the physical plane list
            submesh = self.mesh.cell_sets_dict[key]
            element_id = submesh.values()                         # Get the values
            element_id_list = list(element_id)                    # value to list
            elementcon_part = self.elementcon[element_id_list[0]] # Read the nodes
            values = elementcon_part                              # Make the nodes list a vertical list
            region_info[key] = values                             # Save the nodes list to the dict.
        if self.printing:
            #print(region_info)
            print("- Region information has been assemblied.")
            print('\n')
        toc = time.time()-tic
        if self.printing:
            print('Elapsed time is:',toc,'seconds')
            print("\n")
        
        # returns
        self.region_info = region_info
        return region_info

    def save_mesh(self):
        """
        Produce and save the mesh plot.
        It draws the diagrams of the planes and boundaries one by one
        """
        tic = time.time()
        print('---------------------------------------------------------')
        print('------------------ Save mesh figure ---------------------')
        print('---------------------------------------------------------')
        #plt.figure().set_size_inches(30,30)
        plt.figure()
        
        # The mesh of physical planes
        for key in self.region_info:
            region_elementcon = self.region_info[key]                                       # Get the nodes
            plt.triplot(self.x_coord, self.y_coord, region_elementcon, linewidth=0.5) # Plot
            # plt.axis('equal')
            # plt.axis('off')
       
        # The mesh of boundaries
        for key in self.boundary_info:
            boundary_node_id = self.boundary_info[key]                                       # Get the nodes
            plt.scatter(self.x_coord[boundary_node_id], self.y_coord[boundary_node_id], s=1) # Plot
            # plt.axis('equal')
            # plt.axis('off')

        print("- Mesh figure has been saved.")
        plt.savefig('Mesh.jpeg', dpi=300,  bbox_inches = 'tight', pad_inches = 0.1) # Save the Image
        print('\n')
        toc = time.time()-tic
        print('Elapsed time is:',toc,'seconds')
        print("\n")
        
        #plt.show()
        
        
if __name__ == "__main__":
    dir = os.getcwd()
    mesh_filename = "mesh1.msh"

    mesh_reader = ReadMesh(dir,mesh_filename)
    mesh_reader.get_mesh()    
    bd = mesh_reader.get_boundary_info()
    rg = mesh_reader.get_region_info()
    mesh_reader.save_mesh()
    print(mesh_reader.x_coord)
    print(mesh_reader.y_coord)
    print(mesh_reader.region_info.keys())
    # print(mesh_reader.region_info)
    # print(mesh_reader.boundary_info)
    plt.title("Triangulation of the area with Gmsh")
    plt.xlabel("x/m")
    plt.ylabel("y/m")
    plt.gca().set_aspect("equal")
    # plt.gca().set_xlim(-60, 60)  # Sets the x-axis limits on the Axes object
    # plt.gca().set_ylim(-60, 60) # Sets the y-axis limits on the Axes object
    plt.show()

    
    
    
    