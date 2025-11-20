"""
A modified version of readMesh.py to be compatible with our model.
"""
import meshio
import numpy as np
import matplotlib.pyplot as plt
import time
import os

# 输入的参数：
# dir——文件夹
# mesh_filename——文件名

class ReadMesh:
    def __init__(self, dir, mesh_file):
        """初始化属性filename"""
        filename = dir+'\\'+mesh_file  # 文件名
        self.filename = filename       # 文件名
        self.mesh = None               # 返回mesh对象
        self.x_coord = None            # x坐标，单位m，按节点号1、2、3...排列
        self.y_coord = None            # y坐标，单位m，按节点号1、2、3...排列
        self.elementcon = None         # 单元由节点表示的连接规则，按单元号1、2、3...排列
        self.linecon = None            # 边界由节点表示的连接规则，按边界标签1、2、3...排列
        self.region = None             # 物理面名称
        self.boundary = None           # 边界名称
        self.boundary_info = None      # 边界信息，字典
        self.region_info = None        # 物理面信息，字典
        self.printing = False
        
        
    # 读取剖分，打印物理面和线域，显示剖分
    def get_mesh(self):
        tic = time.time()    
        # 读取剖分文件
        if self.printing:
            print('---------------------------------------------------------')
            print('--------------------- Read mesh -------------------------')
            print('---------------------------------------------------------')
        mesh = meshio.read(self.filename)
        if self.printing:
            print("- Mesh file has been read. ")
        
        # 组装单元坐标
        x_coord = mesh.points[:,0]                         # x坐标，单位m
        y_coord = mesh.points[:,1]                         # x坐标，单位m
        if self.printing:
            print("- Node coordinate has been assembled.")
        
        # 组装单元连接
        elementcon = mesh.cells_dict['triangle']
        if self.printing:
            print("- Element connection has been assembled.")
        
        # 组装线连接
        linecon = mesh.cells_dict['line']
        if self.printing:
            print("- Line connection has been assembled.\n")
        toc = time.time()-tic
        if self.printing:
            print('Elapsed time is:',toc,'seconds\n')
        
        # 返回参数
        self.mesh = mesh
        self.x_coord = x_coord
        self.y_coord = y_coord
        self.elementcon = elementcon
        self.linecon = linecon
        
        # 打印全物理面和线域
        if self.printing:
            print('---------------------------------------------------------')
            print('--------- List physical regions and boundaries ----------')
            print('---------------------------------------------------------')
        tic = time.time()  
        key = self.mesh.field_data.keys()                 # 获取所有的key
        key_list = list(key)                              # 把key转换为列表
        key_regions = []                                  # 空列表储存物理面
        key_boundaries = []                               # 空列表储存边界
        for key in key_list:
            if self.printing:
                print(key)
            key_name = list(self.mesh.cell_sets_dict[key].keys())[0]     # 获物理面/线键的名称，有两种：“triangle”和"line"
            if key_name == 'triangle':
                key_regions.append(key)                                  # 储存所有的物理面名称
            elif key_name == 'line':
                key_boundaries.append(key)                               # 储存所有的边界条件名称

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
        
        # 返回参数
        self.region = key_regions
        self.boundary = key_boundaries
            
    # 组装边界条件信息
    def get_boundary_info(self):
        tic = time.time()
        if self.printing:
            print(self.printing)
            print('---------------------------------------------------------')
            print('---------- Assemble boundary information ----------------')
            print('---------------------------------------------------------')
        boundary_info = {}                                                               # 先构造一个空字典
        for key in self.boundary:                                                        # 循环每个边界
            subline = self.mesh.cell_sets_dict[key]
            line_id = subline.values()                                                   # 获取所有的value
            line_id_list = list(line_id)                                                 # 把value转换为list
            linecon_part = self.linecon[line_id_list[0]]                                 # 读取某一线的节点号，格式为头到尾
            linecon_part_array = linecon_part.reshape(-1,1)                              # 将每个线段头尾展成一列向量 
            node_id= np.unique(linecon_part_array)                                       # 剔除线段连接时候的重复点
            values = node_id                                                             # 将节点id的array变成列表的一个元素
            boundary_info[key] = values                                                  # 更新字典
        if self.printing:
            #print(boundary_info)
            print("- Boundary information has been assemblied.")
            print('\n')
        toc = time.time()-tic
        if self.printing:
            print('Elapsed time is:',toc,'seconds')
            print("\n")
        
        # 返回参数
        self.boundary_info = boundary_info
        return boundary_info                                                             # 按节点序号排好
    
    # 组装物理面信息
    def get_region_info(self):
        tic = time.time()
        if self.printing:
            print(self.printing)
            print('---------------------------------------------------------')
            print('----------- Assemble region information -----------------')
            print('---------------------------------------------------------')
        region_info = {}                                                                 # 先构造一个空字典        
        for key in self.region:                                                          # 循环每个物理面
            submesh = self.mesh.cell_sets_dict[key]
            element_id = submesh.values()                                                # 获取所有的value
            element_id_list = list(element_id)                                           # 把value转换为list
            elementcon_part = self.elementcon[element_id_list[0]]                        # 读取某一物理面的单元号
            values = elementcon_part                                                     # 将节点id的array变成列表的一个元素
            region_info[key] = values                                                    # 更新字典
        if self.printing:
            #print(region_info)
            print("- Region information has been assemblied.")
            print('\n')
        toc = time.time()-tic
        if self.printing:
            print('Elapsed time is:',toc,'seconds')
            print("\n")
        
        # 返回参数
        self.region_info = region_info
        return region_info                                                              # 按物理域类型排列
    
    # 画出全剖分图，逐个物理面和线来画
    def save_mesh(self):
        tic = time.time()  
        # 画出全剖分图，逐个物理面和线来画
        print('---------------------------------------------------------')
        print('------------------ Save mesh figure ---------------------')
        print('---------------------------------------------------------')
        #plt.figure().set_size_inches(30,30)                                             # 第一位为宽，第二位为高
        plt.figure()
        
        # 画物理面
        for key in self.region_info:                                                     # 画每个物理面的剖分
            region_elementcon = self.region_info[key]                                    # 该物理域的节点连接关系
            plt.triplot(self.x_coord, self.y_coord, region_elementcon, linewidth=0.5)                   # 画这一部分的剖分
            # plt.axis('equal')
            # plt.axis('off')
       
        # 画边界上的点  
        for key in self.boundary_info:                                                   # 画每个边界的剖分
            boundary_node_id = self.boundary_info[key]                                   # 该边界上节点号                                   
            plt.scatter(self.x_coord[boundary_node_id], self.y_coord[boundary_node_id], s=1)  # 画这一部分的边界
            # plt.axis('equal')
            # plt.axis('off')

        print("- Mesh figure has been saved.")
        plt.savefig('Mesh.jpeg', dpi=300,  bbox_inches = 'tight', pad_inches = 0.1)      # 保存图片
        print('\n')
        toc = time.time()-tic
        print('Elapsed time is:',toc,'seconds')
        print("\n")
        
        #plt.show()                                                                      # 用plt.show()会线程阻塞
        
        
if __name__ == "__main__":
    # 定义剖分文件夹和名称
    dir = os.getcwd()                                                              # 文件夹
    mesh_filename = 'mesh1.msh'                                                    # 文件名

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

    
    
    
    