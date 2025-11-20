"""
The final version of our visualization tool
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib import cm


class Visualization:
    def __init__(self, mesh, fontname="Arial", fontsize=20, title=False):
        """
        Load the mesh into the visualization tool.

        :param mesh: The Mesh object to visualize.
        :type mesh: Mesh
        :param fontname: The name of font to use. Default: "Arial".
        :type fontname: str
        :param fontsize: The size of font to use. Default: 20
        :type fontsize: int
        :param title: Whether to show the title. Default: False
        :type title: bool
        """
        plt.rcParams['font.family'] = fontname
        plt.rcParams['font.size'] = fontsize

        self.do_title = title

        self.mesh = mesh
        self.meshR = mesh.mesh
        self.nodes = self.mesh.nodes

        self.ele = self.mesh.elements
        self.com = self.mesh.com
        self.tri_nodes = self.mesh.tri_nodes
        self.tri_elements = self.mesh.tri_elements

        self.Value={"A": self.mesh.A,                                # The magnetic potential
                    "B": self.mesh.B[:, 2],                          # The magnitude of magnetic field
                    "Bx": self.mesh.B[:, 0],                         # The x-component of the magnetic field
                    "By": self.mesh.B[:, 1],                         # The y-component of the magnetic field
                    "mu": np.array(self.mesh.get_Elements_mu()),     # The relative permittivity
                    "J": np.array(self.mesh.get_Elements_J()),       # The current density
                    "Area": np.array(self.mesh.get_Elements_Area())} # The area of the elements

    @staticmethod
    def get_value_name_dict():
        """
        The symbol used in value_name parameter and its corresponding value.
        :return: a dict of symbol and its corresponding value.
        :rtype: dict
        """
        return {"A": "The magnetic potential",
                "B": "The magnitude of magnetic field",
                "Bx": "The x-component of the magnetic field",
                "By": "The y-component of the magnetic field",
                "mu": "The relative permittivity",
                "J": "The current density",
                "Area": "The area of the elements"}

    def mesh_nodes_triangulation(self, color_ele=lambda p: "black",
                                 show_fig=True, grid=False, style="k-", linewidth=1.0, s=1.0):
        """
        Demonstrate the mesh of the area.

        :param color_ele: A function that take the region name and return the color. Default: black for all.
        :type color_ele: Callable[[str], str]
        :param show_fig: Whether to show the figure directly. Default: True.
        :type show_fig: bool
        :param grid: Whether to show the grid lines.
        :type grid: bool
        :param style: The stype of the mesh line.
                      The style parameter for the plt.triplot() method. Default: "k-".
        :type style: str
        :param linewidth: The width of mesh line.
                          The linewidth parameter for the plt.triplot() method. Default: 1.0.
        :type linewidth: float
        :param s: The point size of the edge elements.
                  The s parameter for the plt.scatter() method. Default: 1.0.
        :type s: float
        :return: Four parameters. tp: the list of triplot of all region; bd: the list of scatter of te boundary plots;
                 fig, ax: the return value of plt.subplot.
        :rtype: tuple
        """

        fig, ax = plt.subplots()
        tp = []; bd = []
        for key in self.meshR.region_info:
            region_element = self.meshR.region_info[key]
            triTemp = Triangulation(self.meshR.x_coord, self.meshR.y_coord, region_element)
            tp.append(ax.triplot(triTemp, style, linewidth=linewidth, color=color_ele(key)))

        for key in self.meshR.boundary:
            boundary_node_id = self.meshR.boundary_info[key]
            bd.append(ax.scatter(self.meshR.x_coord[boundary_node_id],
                                 self.meshR.y_coord[boundary_node_id],
                                 s=s, color=color_ele(key)))

        """
        This code was used to label the six stator wires. 
        d = 0.33
        tx = 0.02
        ty = 0.01
        ax.text(d - tx, 0 - ty, "A+")
        ax.text(-d - tx, 0 - ty, "A-")
        ax.text(d * np.cos(2 * np.pi / 3) - tx, d * np.sin(2 * np.pi / 3) - ty, "B+")
        ax.text(-d * np.cos(2 * np.pi / 3) - tx, -d * np.sin(2 * np.pi / 3) - ty, "B-")
        ax.text(d * np.cos(4 * np.pi / 3) - tx, d * np.sin(4 * np.pi / 3) - ty, "C+")
        ax.text(-d * np.cos(4 * np.pi / 3) - tx, -d * np.sin(4 * np.pi / 3) - ty, "C-")
        ax.set_xlim(-0.38, 0.38)
        ax.set_ylim(-0.38, 0.38)
        """

        if self.do_title:
            ax.set_title("Triangulation of the region")
        ax.set_xlabel("$x$/m")
        ax.set_ylabel("$y$/m")
        plt.grid(grid)
        ax.set_aspect('equal')

        if show_fig:
            plt.show()
        return tp, bd, fig, ax

    def mesh_elements_triangulation(self, show_fig=True, grid=False, style="k-", linewidth=1.0):
        """
        Demonstrate the mesh of the area. However, in this method, the nodes are the center of mass of all
        triangle element, and the edges care connecting adjacent centers.

        :param show_fig: Whether to show the figure directly. Default: True.
        :type show_fig: bool
        :param grid: Whether to show the grid lines.
        :type grid: bool
        :param style: The stype of the mesh line.
                      The style parameter for the plt.triplot() method. Default: "k-".
        :type style: str
        :param linewidth: The width of mesh line.
                          The linewidth parameter for the plt.triplot() method. Default: 1.0.
        :type linewidth: float
        :return: Three parameters. me: the triplot return. fig, ax: the return value of plt.subplot.
        """

        fig, ax = plt.subplots()
        me = ax.triplot(self.tri_elements, style, linewidth=linewidth)
        ax.set_xlabel("$x$/m")
        ax.set_ylabel("$y$/m")
        plt.grid(grid)
        ax.set_aspect('equal')

        if show_fig:
            plt.show()
        return me, fig, ax

    def mesh_colormap(self, value_name="A", value_text="magnetic potential",
                      unit="$\\mathrm{T\\cdot m}$", show_fig=True, cmap="jet",
                      vrange=None, shading="flat"):
        """
        Show a color map of a specific value at the nodes.

        :param value_name: The symbol of the value to show. Default: "A".
                           Learn about the supporting symbols and their meaning by call Visualization.get_value_name_dict().
        :type value_name: str
        :param value_text: The complete name or other texts shown in the title. Default: "magnetic potential".
        :type value_text: str
        :param unit: The unit text shown in the colorbar. Default: "$\\mathrm{T\\cdot m}$".
        :type unit: str
        :param show_fig: Whether to show the figure directly. Default: True.
        :type show_fig: bool
        :param cmap: The cmap type used. The cmap parameter of plt.tripcolor(). Default: "jet".
        :type cmap: str
        :param vrange: The range of values that corresponds to the colormap.
                       The vmin and vmax parameters of plt.tripcolor(). Default: None.
        :type vrange: iterable of two
        :param shading: The shading type used.
                        The shading parameter of plt.tripcolor(). Default: "flat".
        :type shading: str
        :return: Three parameters. me: the tripcolor return. fig, ax: the return value of plt.subplot.
        :rtype: tuple
        """

        value = self.Value[value_name]
        if vrange is None:
            vrange = [None, None]

        fig, ax = plt.subplots()
        tpc = ax.tripcolor(self.tri_nodes, value, cmap=cmap, shading=shading, vmin=vrange[0], vmax=vrange[1])
        fig.colorbar(tpc, ax=ax, label=value_text + " of each element / " + unit)
        if self.do_title:
            ax.set_title("The " + value_text + " of the elements in the region")
        ax.set_xlabel('$x$/m')
        ax.set_ylabel('$y$/m')
        ax.set_aspect('equal')

        if show_fig:
            plt.show()
        return tpc, fig, ax

    def mesh_3d(self, value_name="A", value_text="magnetic potential",
                unit="$\\mathrm{T\\cdot m}$", show_fig=True, cmap="jet", vrange=None):
        """
        Show a 3D diagram of a specific value at the nodes. The z-axis represents the value.

        :param value_name: The symbol of the value to show. Default: "A".
                           Learn about the supporting symbols and their meaning by call Visualization.get_value_name_dict().
        :type value_name: str
        :param value_text: The complete name or other texts shown in the title. Default: "magnetic potential".
        :type value_text: str
        :param unit: The unit text shown in the colorbar. Default: "$\\mathrm{T\\cdot m}$".
        :type unit: str
        :param show_fig: Whether to show the figure directly. Default: True.
        :type show_fig: bool
        :param cmap: The cmap type used. The cmap parameter of plt.plot_trisurf(). Default: "jet".
        :type cmap: str
        :param vrange: The range of values that corresponds to the colormap.
                       The vmin and vmax parameters of plt.plot_trisurf(). Default: None.
        :type vrange: iterable of two
        :return: Three parameters. me: the plot_trisurf return. fig, ax: the return value of plt.subplot.
        :rtype: tuple
        """
        if vrange is None:
            vrange = [None, None]
        value = self.Value[value_name]

        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        me = ax.plot_trisurf(self.tri_nodes if value_name == "A" else self.tri_elements,
                             value, cmap=cmap, vmin=vrange[0], vmax=vrange[1])
        fig.colorbar(me, ax=ax, label=value_text + " of each element / " + unit)
        ax.set_xlabel('$x$/m')
        ax.set_ylabel('$y$/m')
        if self.do_title:
            ax.set_title("The " + value_text +" of the elements in the region")
        if show_fig:
            plt.show()
        return me, fig, ax

    def mesh_contourf(self, value_name="A", value_text="magnetic potential",
                      unit="$\\mathrm{T\\cdot m}$", show_fig=True, cmap="jet",
                      vrange=None, levels=25):
        """
        Show the contour map of a certain value in the mesh using tricontourf.

        :param value_name: The symbol of the value to show. Default: "A".
                           Learn about the supporting symbols and their meaning by call Visualization.get_value_name_dict().
        :type value_name: str
        :param value_text: The complete name or other texts shown in the title. Default: "magnetic potential".
        :type value_text: str
        :param unit: The unit text shown in the colorbar. Default: "$\\mathrm{T\\cdot m}$".
        :type unit: str
        :param show_fig: Whether to show the figure directly. Default: True.
        :type show_fig: bool
        :param cmap: The cmap type used. The cmap parameter of plt.tricontourf(). Default: "jet".
        :type cmap: str
        :param vrange: The range of values that corresponds to the colormap.
                       The vmin and vmax parameters of plt.tricontourf(). Default: None.
        :type vrange: iterable of two
        :param levels: Determines the number and positions of the contour lines / regions.
                       The levels parameter of plt.tricontourf(). Default: 25.
        :type levels: int
        :return: Three parameters. cf: the tricontourf return. fig, ax: the return value of plt.subplot.
        :rtype: tuple
        """
        value = self.Value[value_name]

        if vrange is None:
            vrange = [None, None]
        fig, ax = plt.subplots()
        cf = ax.tricontourf(self.tri_nodes if value_name == "A" else self.tri_elements,
                            value, levels=levels, cmap=cmap, vmin=vrange[0], vmax=vrange[1])
        fig.colorbar(cf, ax=ax, label=value_text + " of each element / " + unit)
        if self.do_title:
            ax.set_title("The " + value_text +" of the elements in the region")
        ax.set_xlabel('$x$/m')
        ax.set_ylabel('$y$/m')
        ax.set_aspect('equal')
        if show_fig:
            plt.show()
        return cf, fig, ax

    def mesh_contours(self, value_name="A", value_text="magnetic potential",
                                unit="$\\mathrm{T\\cdot m}$", show_fig=True, cmap="jet",
                                vrange=None, levels=25):
        """
        Show the contour map of a certain value in the mesh using tricontour.

        :param value_name: The symbol of the value to show. Default: "A".
                                   Learn about the supporting symbols and their meaning by call Visualization.get_value_name_dict().
        type value_name: str
        :param value_text: The complete name or other texts shown in the title. Default: "magnetic potential".
        :type value_text: str
        :param unit: The unit text shown in the colorbar. Default: "$\\mathrm{T\\cdot m}$".
        :type unit: str
        :param show_fig: Whether to show the figure directly. Default: True.
        :type show_fig: bool
        :param cmap: The cmap type used. The cmap parameter of plt.tricontour(). Default: "jet".
        :type cmap: str
        :param vrange: The range of values that corresponds to the colormap.
                       The vmin and vmax parameters of plt.tricontour(). Default: None.
        :type vrange: iterable of two
        :param levels: Determines the number and positions of the contour lines / regions.
                       The levels parameter of plt.tricontourf(). Default: 25.
        :type levels: int
        :return: Three parameters. mcs: the tricontour return. fig, ax: the return value of plt.subplot.
        :rtype: tuple
        """
        value = self.Value[value_name]

        if vrange is None:
            vrange = [None, None]

        fig, ax = plt.subplots()
        cs = ax.tricontour(self.tri_nodes if value_name == "A" else self.tri_elements,
                           value, levels=levels, cmap=cmap, vmin=vrange[0], vmax=vrange[1])
        cs.clabel(inline=True, colors='k')
        ax.set_aspect("equal")
        ax.set_xlabel('$x$/m')
        ax.set_ylabel('$y$/m')
        fig.colorbar(cs, ax=ax, label=value_text + " of each element / " + unit)
        if self.do_title:
            ax.set_title("The " + value_text +" of the elements in the region")
        if show_fig:
            plt.show()
        return cs, fig, ax

    def mesh_scatter3d(self, value_name="A", value_text="magnetic potential",
                       unit="$\\mathrm{T\\cdot m}$", show_fig=True, cmap="jet",
                       vrange=None, s=plt.rcParams["lines.markersize"] ** 2, marker="o"):
        """
        Show a 3D colored scatter of a specific value of the mesh.

        :param value_name: The symbol of the value to show. Default: "A".
                           Learn about the supporting symbols and their meaning by call Visualization.get_value_name_dict().
        :type value_name: str
        :param value_text: The complete name or other texts shown in the title. Default: "magnetic potential".
        :type value_text: str
        :param unit: The unit text shown in the colorbar. Default: "$\\mathrm{T\\cdot m}$".
        :type unit: str
        :param show_fig: Whether to show the figure directly. Default: True.
        :type show_fig: bool
        :param cmap: The cmap type used. Default: "jet".
        :type cmap: str
        :param vrange: The range of values that corresponds to the colormap. Default: None.
        :type vrange: iterable of two
        :param s: The size of scatter points.
                  The s parameter of plt.scatter(). Default: plt.rcParams["lines.markersize"] ** 2.
        :type s: float
        :return: Three parameters. sc: the tripcolor return. fig, ax: the return value of plt.subplot.
        :rtype: tuple
        """
        value = self.Value[value_name]

        if vrange is None:
            norm = plt.Normalize(value.min(), value.max())
        else:
            norm = plt.Normalize(vrange[0], vrange[1])
        colors = cm.get_cmap(cmap)(norm(value))

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        X = self.nodes[:, 0] if value_name == "A" else self.com[:, 0]
        Y = self.nodes[:, 1] if value_name == "A" else self.com[:, 1]
        sc = ax.scatter(X, Y, value, c=colors, s=s, marker=marker)

        mappable = cm.ScalarMappable(norm=norm, cmap=cmap)
        mappable.set_array(value)
        fig.colorbar(mappable, ax=ax, shrink=0.6, label=value_text + " of each element / " + unit)
        ax.set_xlabel('$x$/m')
        ax.set_ylabel('$y$/m')
        if self.do_title:
            ax.set_title("The " + value_text +" of the elements in the region")
        if show_fig:
            plt.show()
        return sc, fig, ax

    def mesh_3d2(self, value_name="B", value_text="magnetic flux density",
                 unit="$\\mathrm{T}$", show_fig=True, cmap="jet", vrange=None):
        """
        Show a 3D diagram of a specific value of the elements.
        This method should not be used to show the information at nodes, such as the magnetic potential.
        If you try so, mesh_3d() will be called.

        :param value_name: The symbol of the value to show. Default: "B".
                           Learn about the supporting symbols and their meaning by call Visualization.get_value_name_dict().
        :type value_name: str
        :param value_text: The complete name or other texts shown in the title. Default: "magnetic flux density".
        :type value_text: str
        :param unit: The unit text shown in the colorbar. Default: "$\\mathrm{T}$".
        :type unit: str
        :param show_fig: Whether to show the figure directly. Default: True.
        :type show_fig: bool
        :param cmap: The cmap type used. Default: "jet".
        :type cmap: str
        :param vrange: The range of values that corresponds to the colormap. Default: None.
        :type vrange: iterable of two
        :return: Three parameters. mesh: a Poly3DCollection. fig, ax: the return value of plt.subplot.
        :rtype: tuple
        """
        if value_name == "A":
            self.mesh_3d()
            return
        value = self.Value[value_name]

        if vrange is None:
            norm = plt.Normalize(value.min(), value.max())
        else:
            norm = plt.Normalize(vrange[0], vrange[1])
        colors = cm.get_cmap(cmap)(norm(value))

        vertices = []
        for e, val in zip(self.mesh.get_Elements()[0], value):
            vertices.append([(self.mesh.nodes[i, 0], self.mesh.nodes[i, 1], val) for i in e])

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        mesh = Poly3DCollection(vertices, facecolors=colors)
        ax.add_collection3d(mesh)

        mappable = cm.ScalarMappable(norm=norm, cmap=cmap)
        mappable.set_array(value)
        fig.colorbar(mappable, ax=ax, shrink=0.6, label=value_text + " of each element / " + unit)
        ax.set_xlabel('$x$/m')
        ax.set_ylabel('$y$/m')
        if self.do_title:
            ax.set_title("The " + value_text + " of the elements in the region")

        if show_fig:
            plt.show()
        return mesh, fig, ax

    def mesh_colormap2(self, value_name="B", value_text="magnetic flux density",
                       unit="$\\mathrm{T}$", show_fig=True, cmap="jet",
                       vrange=None, shading="flat"):
        """
        Show a color map of a specific value at the nodes.
        This method should not be used to show the information at nodes, such as the magnetic potential.
        If you try so, mesh_colormap() will be called.

        :param value_name: The symbol of the value to show. Default: "A".
                           Learn about the supporting symbols and their meaning by call Visualization.get_value_name_dict().
        :type value_name: str
        :param value_text: The complete name or other texts shown in the title. Default: "magnetic potential".
        :type value_text: str
        :param unit: The unit text shown in the colorbar. Default: "$\\mathrm{T\\cdot m}$".
        :type unit: str
        :param show_fig: Whether to show the figure directly. Default: True.
        :type show_fig: bool
        :param cmap: The cmap type used. The cmap parameter of plt.tripcolor(). Default: "jet".
        :type cmap: str
        :param vrange: The range of values that corresponds to the colormap.
                       The vmin and vmax parameters of plt.tripcolor(). Default: None.
        :type vrange: iterable of two
        :param shading: The shading type used.
                                The shading parameter of plt.tripcolor(). Default: "flat".
        :type shading: str
        :return: Three parameters. tpc: the tripcolor return. fig, ax: the return value of plt.subplot.
        :rtype: tuple
        """
        if value_name == "A":
            return self.mesh_colormap()
        value = self.Value[value_name]
        if vrange is None:
            vrange = [None, None]

        fig, ax = plt.subplots()
        tpc = ax.tripcolor(self.tri_elements, value, cmap=cmap, shading=shading, vmin=vrange[0], vmax=vrange[1])
        fig.colorbar(tpc, ax=ax, label=value_text + " of each element / " + unit)
        if self.do_title:
            ax.set_title("The " + value_text + " of the elements in the region")
        ax.set_xlabel('$x$/m')
        ax.set_ylabel('$y$/m')
        ax.set_aspect('equal')

        if show_fig:
            plt.show()
        return tpc, fig, ax

    def mesh_field_quiver(self, grid=False, cons=lambda x,y: True,
                          grid_size=(5, 5), grid_num=None,
                          show_length=False, show_fig=True,
                          cmap="jet", vrange=None, scale=None):
        """
        Show the magnetic vector field quiver.

        :param grid: Whether to show the grid lines. Default: False.
        :type grid: bool
        :param cons: A function that tells whether to show the vector at (x, y). Default: True for all.
        :type cons: Callable[[float, float], bool]
        :param grid_size: the distance between two sampling points on x and y direction. Default: (5, 5).
        :type grid_size: iterable of two
        :param grid_num: the number of sampling points on x and y direction. Default: taken according to the grid size.
        :type grid_num: iterable of two
        :param show_length: Whether to use the length of the vector to represent its value. Default: False.
        :type show_length: bool
        :param show_fig: Whether to show the figure directly. Default: True.
        :type show_fig: bool
        :param cmap: The cmap type used. The cmap parameter of plt.quiver(). Default: "jet".
        :type cmap: str
        :param vrange: The range of values that corresponds to the colormap. Default: None.
        :type vrange: iterable of two
        :param scale: The size of the quiver arrows. Default: None.
        :type scale: float
        :return: Three parameters. cf: the quiver return. fig, ax: the return value of plt.subplot.
        """
        Xmin = min(self.nodes[:, 0])
        Xmax = max(self.nodes[:, 0])
        Ymin = min(self.nodes[:, 1])
        Ymax = max(self.nodes[:, 1])

        grid_num = grid_num if grid_num else ((Xmax-Xmin) // grid_size[0], (Ymax-Ymin) // grid_size[1])
        if vrange is None:
            vrange = [None, None]

        X = []; Y = []; U = []; V = []; value = []
        for i in np.linspace(Xmin, Xmax, grid_num[0]):
            for j in np.linspace(Ymin, Ymax, grid_num[1]):
                if cons(i, j):
                    X.append(i); Y.append(j)
                    B = self.mesh.get_B(i, j)
                    U.append(B[0])
                    V.append(B[1])
                    value.append(B[2])
                    if (not show_length) and value[-1]:
                        U[-1] /= value[-1]
                        V[-1] /= value[-1]

        if vrange is not [None, None]:
            value = np.clip(value, vrange[0], vrange[1])

        fig, ax = plt.subplots()
        cf = ax.quiver(X, Y, U, V, value, cmap=cmap, scale=scale)

        fig.colorbar(cf, ax=ax, label="Magnetic flux density / $\\mathrm{T}$")  # Add a colorbar
        if self.do_title:
            ax.set_title("Magnetic Field")
        ax.set_xlabel('$x$/m')
        ax.set_ylabel('$y$/m')
        ax.set_aspect("equal")
        plt.grid(grid)
        if show_fig:
            plt.show()

        return cf, fig, ax
