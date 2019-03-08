import numpy as np
import matplotlib.pyplot as plt
import math
from operator import itemgetter


def prepare_xy_data_from_file(filename, columns=2):
    """Prepares two lists based on data in txt file"""
    with open(filename) as f:
        lines = f.readlines()
        # print(lines)
        x = [line.split()[0] for line in lines]
        try:
            y = [line.split()[1] for line in lines]
        except:
            y = [line.replace("\n", "") for line in lines]
            y = [line.split(" ")[1] for line in lines]
        if columns == 3:
            y2 = [line.split()[2] for line in lines]
    for i in range(len(x)):
        x[i] = float(x[i])

    for i in range(len(y)):
        try:
            y[i] = float(y[i])
        except:
            print(y[i])
    if columns == 3:
        for i in range(len(y2)):
            y2[i] = float(y2[i])
        return (x, y, y2)
    return (x, y)


def sort_depth_value(x, y):
    """Sorts values according to the depth, x ... value, y ... depth"""
    assert len(x) == len(y)
    list1 = [[x[i], y[i]] for i in range(len(x))]

    list1 = sorted(list1, key=itemgetter(0))
    x = [item[0] for item in list1]
    y = [item[1] for item in list1]

    return (x, y)


class MyPlot:
    """Class for automated plotting of xy-scatter data"""
    def __init__(self):
        self.x_list = []
        self.y_list = []
        self.c_list = []
        self.label_list = []
        self.linewidth = []
        pass

    def load_data(self, filepath, sort=False):
        x, y = prepare_xy_data_from_file(filepath)
        if sort:
            x, y = sort_depth_value(x, y)
        return (x, y)

    def append_data(self, x, y, c, label, linewidth=2.0):
        self.x_list.append(x)
        self.y_list.append(y)
        self.c_list.append(c)
        self.label_list.append(label)
        self.linewidth.append(linewidth)

    def construct_plot(
        self,
        title,
        xlabel,
        ylabel,
        save=False,
        xymin=False,
        xymax=False,
        figsize=(7, 5),
    ):
        fig = plt.figure(figsize=figsize)
        ax1 = fig.add_subplot(111)
        ax1.set_title(title)
        ax1.set_xlabel(xlabel)
        ax1.set_ylabel(ylabel)
        for i in range(len(self.x_list)):
            ax1.plot(
                self.x_list[i],
                self.y_list[i],
                c=self.c_list[i],
                label=self.label_list[i],
                linewidth=self.linewidth[i],
            )

        # ax1.plot(y,x, c='r', label='Mean T depth profile, Viscoplasticity - case 1', linewidth=1.0)
        leg = ax1.legend()
        if xymin:
            ax1.set_xlim(xmin=xymin[0])
            ax1.set_ylim(ymin=xymin[1])
        if xymax:
            ax1.set_xlim(xmax=xymax[0])
            ax1.set_ylim(ymax=xymax[1])
        if save:
            fig.savefig("img\\" + save, format="png",
                        dpi=200, bbox_inches="tight")
        return (ax1, fig)


def logarithmic(x):
    """recalculates whole list of values to decadic logarithm - used in plots"""
    for i in range(len(x)):
        if x[i] <= 0:
            x[i] = 1e-10
        x[i] = math.log10(x[i])
    return x

dir1 = "ALE_results_c1\\"
dir2 = "Blankenbach_results_c1\\"

dir3 = ("40x40\\Viscoplastic_results_c2_8\\")
#dir4 = "40x40\\Viscoplastic_results_c2_8\\"
#dir5 = "40x40\\Viscoplastic_results_c3_8\\"
#dir6 = "40x40\\Viscoplastic_results_c4_8\\"
#dir7 = "40x40\\Viscoplastic_results_c5_8\\"
case_no = "2"


plot2 = MyPlot()
x, y = prepare_xy_data_from_file(dir1 + "h1.dat")
x2, y2 = prepare_xy_data_from_file(dir2 + "h1.dat")
x3, y3 = prepare_xy_data_from_file("klara_results\\ale_top_bndry.dat")
x4, y4 = prepare_xy_data_from_file("klara_results\\dtopo_top_bndry.dat")
plot2.append_data(
    x, y, 'r', 'Arbitrary Lagrangian Eulerian method', linewidth=2.0)
plot2.append_data(x2, y2, 'k', 'Dynamic topography', linewidth=2.0)
plot2.append_data(x3,y3, c='g', label='Klara - ALE')
plot2.append_data(x4,y4, c='b', label='Klara - Dynamic topography')
plot2.construct_plot("Displacement", "Time", "$\Delta h$",
                     save="klara_dynamic_topography.png")



plot3 = MyPlot()
x, y = prepare_xy_data_from_file(dir2 + "Nusselt.dat")
x2, y2 = prepare_xy_data_from_file(dir2 + "Nusselt_benchmark.dat")
#x3, y3 = prepare_xy_data_from_file("klara_results\\ale_top_bndry.dat")
#x4, y4 = prepare_xy_data_from_file("klara_results\\dtopo_top_bndry.dat")
plot3.append_data(
    x, y, 'r', 'Nusselt number', linewidth=2.0)
plot3.append_data(x2, y2, 'k', 'Benchmark nusselt number', linewidth=2.0)
#plot2.append_data(x3,y3, c='g', label='Klara - ALE')
#plot2.append_data(x4,y4, c='b', label='Klara - Dynamic topography')
plot3.construct_plot("Nusselt number", "Time", "$Nu$",
                     save="blankenbach_nusselt.png",xymin=[0,0],xymax=[0.30,7.5])


plot4 = MyPlot()
x, y = prepare_xy_data_from_file(dir2 + "Rmsvel.dat")
x2, y2 = prepare_xy_data_from_file(dir2 + "Rmsvel_benchmark.dat")
#x3, y3 = prepare_xy_data_from_file("klara_results\\ale_top_bndry.dat")
#x4, y4 = prepare_xy_data_from_file("klara_results\\dtopo_top_bndry.dat")
plot4.append_data(
    x, y, 'r', 'Rms velocity', linewidth=2.0)
plot4.append_data(x2, y2, 'k', 'Benchmark rms velocity', linewidth=2.0)
#plot2.append_data(x3,y3, c='g', label='Klara - ALE')
#plot2.append_data(x4,y4, c='b', label='Klara - Dynamic topography')
plot4.construct_plot("Rms velocity", "Time", "$v_{rms}$",
                     save="blankenbach_rmsvel.png",xymin=[0,0],xymax=[0.30,100])







"""
plot1 = MyPlot()
x, y = prepare_xy_data_from_file(dir1 + "h1.dat")
x2, y2 = prepare_xy_data_from_file(dir1 + "h2.dat")
x3, y3 = prepare_xy_data_from_file(dir1 + "h3.dat")
x4, y4 = prepare_xy_data_from_file(dir1 + "h4.dat")
x5, y5 = prepare_xy_data_from_file(dir1 + "h5.dat")
plot1.append_data(x, y, 'r', '(0.0,1.0)', linewidth=2.0)
plot1.append_data(x2, y2, 'g', '(0.25,1.0)', linewidth=2.0)
plot1.append_data(x3, y3, 'c', '(0.5,1.0)', linewidth=2.0)
plot1.append_data(x4, y4, 'b', '(0.75,1.0)', linewidth=2.0)
plot1.append_data(x5, y5, 'm', '(1.0,1.0)', linewidth=2.0)
plot1.construct_plot("Displacement", "Time", "$\Delta h$", save="h12345.png")


plot2 = MyPlot()
x, y = prepare_xy_data_from_file(dir1 + "h1.dat")
x2, y2 = prepare_xy_data_from_file(dir2 + "h1.dat")
x3, y3 = prepare_xy_data_from_file("klara_results\\ale_top_bndry.dat")
x4, y4 = prepare_xy_data_from_file("klara_results\\dtopo_top_bndry.dat")
plot2.append_data(
    x, y, 'r', 'Arbitrary Lagrangian Eulerian method', linewidth=2.0)
plot2.append_data(x2, y2, 'k', 'Dynamic topography', linewidth=2.0)
#plot2.append_data(x3,y3, c='g', label='Klara - ALE')
#plot2.append_data(x4,y4, c='b', label='Klara - Dynamic topography')
plot2.construct_plot("Displacement", "Time", "$\Delta h$",
                     save="dynamic_topography.png")


plot3 = MyPlot()
x, y = prepare_xy_data_from_file(dir1 + "Nusselt.dat")
x2, y2 = prepare_xy_data_from_file(dir2 + "Nusselt.dat")
x3, y3, y6 = prepare_xy_data_from_file(
    "klara_results\\ale_nusselt_vrms.dat", 3)
x4, y4, y5 = prepare_xy_data_from_file(
    "klara_results\\dtopo_nusselt_vrms.dat", 3)
plot3.append_data(
    x, y, 'r', 'Arbitrary Lagrangian Eulerian method', linewidth=2.0)
plot3.append_data(x2, y2, 'k', 'Fixed domain', linewidth=2.0)
#plot3.append_data(x3,y3, 'g', 'Klara - ALE')
#plot3.append_data(x4,y4, 'b', 'Klara - Dynamic topography')
plot3.construct_plot("Nusselt number", "Time", "$Nu$", save="nusselt.png")


plot4 = MyPlot()
x, y = prepare_xy_data_from_file(dir1 + "Rmsvel.dat")
x2, y2 = prepare_xy_data_from_file(dir2 + "Rmsvel.dat")
x3, y3, y6 = prepare_xy_data_from_file(
    "klara_results\\ale_nusselt_vrms.dat", 3)
x4, y4, y5 = prepare_xy_data_from_file(
    "klara_results\\dtopo_nusselt_vrms.dat", 3)
plot4.append_data(
    x, y, 'r', 'Arbitrary Lagrangian Eulerian method', linewidth=2.0)
plot4.append_data(x2, y2, 'k', 'Fixed domain', linewidth=2.0)
#plot4.append_data(x3,y6, 'g', 'Klara - ALE')
#plot4.append_data(x4,y5, 'b', 'Klara - Dynamic topography')
plot4.construct_plot("Rms velocity", "Time", "$v_{rms}$", save="rmsvel.png")
"""

# VISCOPLASTICITY ###################x


# for dir3 in [dir3,dir4,dir5,dir6,dir7]
plot5 = MyPlot()
x, y = plot5.load_data(dir3 + "Rmsvel.dat")
# x2,y2=plot5.load_data(dir4+"Rmsvel.dat")
plot5.append_data(
    x,
    y,
    "r",
    "Viscoplasticity - case " + case_no + ", Viscoplasticity - case " + case_no,
)
# plot5.append_data(x2,y2,'b', 'Viscoplasticity - case 3, Viscoplasticity - case 3')
plot5.construct_plot(
    "Rms velocity",
    "Time",
    "$v_{rms}$",
    save="vp_rmsvel_c" + case_no + ".png",
    xymin=[0, 0],
    xymax=[0.5, 200],
)


plot6 = MyPlot()
x, y = plot6.load_data(dir3 + "nus_bot.dat")
x2, y2 = plot6.load_data(dir3 + "nus_top.dat")
plot6.append_data(
    x, y, "r", "Nusselt bottom, Viscoplasticity - case " + case_no)
plot6.append_data(
    x2, y2, "b", "Nusselt top, Viscoplasticity - case " + case_no)
# x3,y3=plot6.load_data(dir4+"nus_bot.dat")
# x4,y4=plot6.load_data(dir4+"nus_top.dat")
# plot6.append_data(x3,y3,'g', 'Nusselt bottom, Viscoplasticity - case 3')
# plot6.append_data(x4,y4,'k', 'Nusselt top, Viscoplasticity - case 3')
plot6.construct_plot(
    "Nusselt number",
    "Time",
    "$Nu$",
    save="vp_nusselt_c" + case_no + ".png",
    xymin=[0, 1],
)

# print(x[10],x2[10],x3[10],x4[10])

plot7 = MyPlot()
x, y = plot7.load_data(dir3 + "mean_t.dat")
plot7.append_data(x, y, "r", "Mean T, Viscoplasticity - case " + case_no)

# x2,y2=plot7.load_data(dir4+"mean_t.dat")
# plot7.append_data(x2,y2,'k', 'Mean T, Viscoplasticity - case 3')
plot7.construct_plot(
    "Mean T", "Time", "$<T>$", save="vp_mean_T_c" + case_no + ".png", xymin=[0, 0.5]
)


plot8 = MyPlot()
x, y = plot8.load_data(dir3 + "T_depth.dat", sort=True)
x = x[::2]
y = y[::2]
# x2,y2=plot8.load_data(dir4+"T_depth.dat",sort=True)
plot8.append_data(
    y, x, "r", "Mean T depth profile, Viscoplasticity - case " + case_no)
# plot8.append_data(y2,x2,'g', 'Mean T depth profile, Viscoplasticity - case 3')
plot8.construct_plot(
    "Mean T depth profile",
    "Temperature",
    "Depth",
    save="vp_t_depth_profile_c" + case_no + ".png",
    figsize=(5, 5),
    xymin=[0,0],
    xymax=[1,1]
    
)


plot9 = MyPlot()
x, y = plot9.load_data(dir3 + "eta_depth.dat", sort=True)
y = logarithmic(y)
x = x[::2]
y = y[::2]
plot9.append_data(
    y, x, "r", "$\eta$ depth profile, Viscoplasticity - case " + case_no)

# x2,y2=plot9.load_data(dir4+"eta_depth.dat",sort=True)
# y2=logarithmic(y2)
# plot9.append_data(y2,x2,'g', '$\eta$ depth profile, Viscoplasticity - case 3')

plot9.construct_plot(
    "$\eta$ depth profile",
    "$log_{10}\eta$",
    "Depth",
    save="vp_eta_depth_profile_c" + case_no + ".png",
    xymin=[-5, 0],
    xymax=[0, 1],
    figsize=(5, 5),
)


plot10 = MyPlot()
x, y = plot10.load_data(dir3 + "rms_depth.dat", sort=True)
x = x[::2]
y = y[::2]
plot10.append_data(
    y, x, "r", "$v_{rms}$ depth profile, Viscoplasticity - case " + case_no
)
# x2,y2=plot10.load_data(dir4+"rms_depth.dat",sort=True)
# plot10.append_data(y2,x2,'g', '$v_rms$ depth profile, Viscoplasticity - case 3')

plot10.construct_plot(
    "$v_{rms}$ depth profile",
    "$v_{rms}$",
    "Depth",
    save="vp_rms_depth_profile_c" + case_no + ".png",
    xymin=[90, 0],
    xymax=[220, 1],
    figsize=(5, 5),
)


plot11 = MyPlot()
x, y, z = prepare_xy_data_from_file(dir3 + "eta_minmax.dat", columns=3)
y = logarithmic(y)
z = logarithmic(z)
plot11.append_data(
    x, y, "r", "$\eta_{min}$, Viscoplasticity - case " + case_no)
plot11.append_data(
    x, z, "g", "$\eta_{max}$, Viscoplasticity - case " + case_no)
plot11.construct_plot(
    "$\eta_{min/max}$",
    "$Time$",
    "$log_{10}\eta$",
    save="vp_eta_minmax_c" + case_no + ".png",
    xymin=[0, -8],
    figsize=(5, 5),
)

"""
plot12 = MyPlot()
x, y, z = prepare_xy_data_from_file(dir3 + "eta_minmax.dat", columns=3)
y = logarithmic(y)
z = logarithmic(z)
plot12.append_data(x, y, 'r', '$\eta_{max}$, Viscoplasticity - case 3')
plot12.append_data(x, z, 'g', '$\eta_{min}$, Viscoplasticity - case 3')
plot12.construct_plot("$\eta_{min/max}$", "$Time$", "$log_{10}\eta$",
                      save="vp_eta_minmax_c3.png", xymin=[0, -8], figsize=(5, 5))
"""

plot13 = MyPlot()
x, y, z = prepare_xy_data_from_file(dir3 + "top_v_x.dat", columns=3)
plot13.append_data(x, y, "r", "$v^{surf}$, Viscoplasticity - case " + case_no)
plot13.append_data(
    x, z, "g", "$v_{max}^{surf}$, Viscoplasticity - case " + case_no)
plot13.construct_plot(
    "$v^{surf}$",
    "$Time$",
    "$v^{surf}$",
    save="top_v_x_c" + case_no + ".png",
    xymin=[0, 0],
    figsize=(7, 5),
)

"""
plot14 = MyPlot()
x, y, z = prepare_xy_data_from_file(dir4 + "top_v_x.dat", columns=3)
plot14.append_data(x, y, 'r', '$v^{surf}$, Viscoplasticity - case 3')
plot14.append_data(x, z, 'g', '$v_{max}^{surf}$, Viscoplasticity - case 3')
plot14.construct_plot("$v^{surf}$", "$Time$", "$v^{surf}$", save="top_v_x_c3.png", xymin=[
                      0, 0], xymax=[0.2, 1500], figsize=(7, 5))
"""

plot15 = MyPlot()
x, y = plot15.load_data(dir3 + "Rmsvel.dat")
x2, y2 = plot15.load_data(dir3 + "nus_top.dat")
plot15.append_data(y, y2, "r", "$v^{surf}$, Viscoplasticity - case " + case_no)
plot15.append_data(
    x2, y2, "g", "$v_{max}^{surf}$, Viscoplasticity - case " + case_no)
plot15.construct_plot(
    "$v^{surf}$",
    "$Time$",
    "$v^{surf}$",
    save="top_v_x_c" + case_no + ".png",
    xymin=[30, 2.5],
    xymax=[100, 7.5],
    figsize=(7, 5),
)
