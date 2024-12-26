import matplotlib as mpl
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib import cm
from scipy.ndimage import uniform_filter1d
import re
import os
import json
import numpy as np
import pandas as pd
import json
import sys
from vtk import vtkUnstructuredGridReader
from vtk.util import numpy_support as VN
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib import cm
from multiprocessing import Pool

def elastic_mat(e,nu):
    fac = (0.5 * (1 - (2 * nu)))
    return (e / ((1 + nu) * (1 - (2 * nu)))) * np.array(
        [[(1-nu),nu,nu,0,0,0],
            [nu,(1-nu),nu,0,0,0],
            [nu,nu,(1-nu),0,0,0],
            [0,0,0,fac,0,0],
            [0,0,0,0,fac,0],
            [0,0,0,0,0,fac]])

def trace(sig):
    return sig[0]+sig[1]+sig[2]

def get_pq(eps_xx,eps_yy,eps_xy,dt,dc,ds):
    de = elastic_mat(1e9,0.24)
    sig = np.matmul(de,np.array([eps_xx,eps_yy,0,0,0,eps_xy]))
    p = trace(sig)/3
    s = sig-np.array([p,p,p,0,0,0])
    if p>0:
        p = p * (1-dt)
    else:
        p = p * (1-dc)
    q = np.sqrt(np.matmul(s,s.transpose())) * (1-ds)
    return -p,q

def damage_map(p,q,damage_t,damage_c,damage_s):
    pr = 0
    pr = np.copy(p)
    pr[p>0] = p[p>0] * (1-damage_c)
    pr[p<=0] = p[p<=0] * (1-damage_t)
    qr = q * (1-damage_s)
    return pr,qr

def plot_yield_surface_mc(name,d_t,d_s,d_c):
    scale = 100
    xmin, xmax, ymin, ymax = -20e3*scale, 30e3*scale, 0e3, 30e3*scale
    ticks_frequency = 1
    theta = 50 * np.pi/180
    c = 26e3 * 100
    A = c
    B = np.tan(theta)
    sig_m = np.linspace(0,10000e3,100) - A
    tau_m = np.linspace(0,10000e3,100) * B
    # plt.plot(p, q,label="Undamaged yield surface")
    #_damaged,q_damaged = damage_map(p,q,d_t,d_c,d_s)
    plt.plot(sig_m,tau_m,label="d = {}".format(name))
    plt.legend()
    ax = plt.gca()
    ax.set(xlim=(xmin-1, xmax+1), ylim=(ymin, ymax+1), aspect='equal')
    #ax.set(aspect='equal')
    ax.spines['bottom'].set_position('zero')
    ax.spines['left'].set_position('zero')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel('$x$', size=14, labelpad=-24, x=1.02)
    ax.set_ylabel('$y$', size=14, labelpad=-21, y=1.02, rotation=0)
    ax.grid(which='both', color='grey', linewidth=1, linestyle='-', alpha=0.2)

def plot_yield_varsurface_mc(name,angle,scale):
    xmin, xmax, ymin, ymax = -26e3*scale, 40e3*scale, 0e3, 50e3*scale
    ticks_frequency = 1
    theta = angle * np.pi/180
    c = 26e3 * scale
    A = c
    B = np.tan(theta)
    sig_m = np.linspace(0,-10000e3,100) + A/B
    tau_m = sig_m * B + A
    # plt.plot(p, q,label="Undamaged yield surface")
    #_damaged,q_damaged = damage_map(p,q,d_t,d_c,d_s)
    plt.plot(sig_m,tau_m,label="d = {}".format(name))
    plt.legend()
    ax = plt.gca()
    ax.set(xlim=(xmin-1, xmax+1), ylim=(ymin, ymax+1), aspect='equal')
    #ax.set(aspect='equal')
    ax.spines['bottom'].set_position('zero')
    ax.spines['left'].set_position('zero')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel('$x$', size=14, labelpad=-24, x=1.02)
    ax.set_ylabel('$y$', size=14, labelpad=-21, y=1.02, rotation=0)
    ax.grid(which='both', color='grey', linewidth=1, linestyle='-', alpha=0.2)

def plot_yield_surface_mc_damaged(name):
    scale = 100
    xmin, xmax, ymin, ymax = -20e3, 300e3, 0e3, 300e3
    ticks_frequency = 1
    theta = 30 * np.pi/180
    c = 26e3 * 0
    A = c * np.cos(theta)
    B = np.sin(theta)
    sig_m = np.linspace(0,10000e3,100) - A
    tau_m = np.linspace(0,10000e3,100) * B
    # plt.plot(p, q,label="Undamaged yield surface")
    #_damaged,q_damaged = damage_map(p,q,d_t,d_c,d_s)
    plt.plot(sig_m,tau_m,label="d = {}".format(name))
    plt.legend()
    ax = plt.gca()
    ax.set(xlim=(xmin-1, xmax+1), ylim=(ymin, ymax+1), aspect='equal')
    #ax.set(aspect='equal')
    ax.spines['bottom'].set_position('zero')
    ax.spines['left'].set_position('zero')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel('$x$', size=14, labelpad=-24, x=1.02)
    ax.set_ylabel('$y$', size=14, labelpad=-21, y=1.02, rotation=0)
    ax.grid(which='both', color='grey', linewidth=1, linestyle='-', alpha=0.2)


def plot_yield_surface(name,d_t,d_s,d_c):
    scale = 100
    xmin, xmax, ymin, ymax = -20e3, 500e3, 0e3, 200e3
    ticks_frequency = 1
    theta = 50 * np.pi/180
    c = 26e3 * 1
    A = (3 * c * np.cos(theta))/(np.sqrt(9 + (3 * np.power(np.sin(theta),2))))
    B = np.sin(theta)/(np.sqrt(9 + (3 * np.power(np.sin(theta),2))))
    p = np.linspace(0,10000e3,100) - A
    q = np.linspace(0,10000e3,100) * B
    # plt.plot(p, q,label="Undamaged yield surface")
    p_damaged,q_damaged = damage_map(p,q,d_t,d_c,d_s)
    plt.plot(p_damaged, q_damaged,label="d = {}".format(name))
    plt.legend()
    ax = plt.gca()
    ax.set(xlim=(xmin-1, xmax+1), ylim=(ymin, ymax+1), aspect='equal')
    #ax.set(aspect='equal')
    ax.spines['bottom'].set_position('zero')
    ax.spines['left'].set_position('zero')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel('$x$', size=14, labelpad=-24, x=1.02)
    ax.set_ylabel('$y$', size=14, labelpad=-21, y=1.02, rotation=0)
    ax.grid(which='both', color='grey', linewidth=1, linestyle='-', alpha=0.2)

def plot_yield_surface_undamaged(name,d_t,d_s,d_c):
    scale = 100
    xmin, xmax, ymin, ymax = -200e3, 40000e3, 0e3, 20000e3
    ticks_frequency = 1
    theta = 50 * np.pi/180
    c = 26e3 * 100
    A = (3 * c * np.cos(theta))/(np.sqrt(9 + (3 * np.power(np.sin(theta),2))))
    B = np.sin(theta)/(np.sqrt(9 + (3 * np.power(np.sin(theta),2))))
    p = np.linspace(0,100000e3,100) - A
    q = np.linspace(0,100000e3,100) * B
    # plt.plot(p, q,label="Undamaged yield surface")
    p_damaged,q_damaged = damage_map(p,q,d_t,d_c,d_s)
    plt.plot(p_damaged, q_damaged,label="d = {}".format(name))
    plt.legend()
    ax = plt.gca()
    ax.set(xlim=(xmin-1, xmax+1), ylim=(ymin, ymax+1), aspect='equal')
    #ax.set(aspect='equal')
    ax.spines['bottom'].set_position('zero')
    ax.spines['left'].set_position('zero')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel('$x$', size=14, labelpad=-24, x=1.02)
    ax.set_ylabel('$y$', size=14, labelpad=-21, y=1.02, rotation=0)
    ax.grid(which='both', color='grey', linewidth=1, linestyle='-', alpha=0.2)


def get_data(filename):
    global data_name
    reader = vtkUnstructuredGridReader()
    reader.SetFileName(filename)
    reader.ReadAllVectorsOn()
    reader.ReadAllScalarsOn()
    reader.Update()

    data = reader.GetOutput()

    vtk_points = data.GetPoints()
    xyz3d = vtk_to_numpy( vtk_points.GetData() )
    xy = xyz3d[:,0:2]
    scalar_names = [reader.GetScalarsNameInFile(i) for i in range(0, reader.GetNumberOfScalarsInFile())]
    scalar_data = data.GetPointData()
    #scalar_names = scalar_data.GetArrayNames()
    def GetScalar(scalar_name):
        return vtk_to_numpy(scalar_data.GetArray(scalar_names.index(scalar_name)))
    lx = GetScalar("size_x")
    ly = GetScalar("size_y")
    damage = GetScalar(data_name)
    #damage = GetScalar("plastic_strain")
    #damage = GetScalar("damage")
    return pd.DataFrame({"coord_x":xy[:,0], "coord_y":xy[:,1],"lx":lx,"ly":ly,
                         "damage":damage,
                         "damage-shear":GetScalar("damage-shear"),
                         "damage-tension":GetScalar("damage-tension"),
                         "damage-compression":GetScalar("damage-compression"),
                         "s_1":GetScalar("s_1"),
                         "s_3":GetScalar("s_3"),
                         "su_1":GetScalar("su_1"),
                         "su_3":GetScalar("su_3"),
                         "unique-id":GetScalar("unique-id"),
                         "p-undamaged":GetScalar("p-undamaged"),
                         "p-damaged":GetScalar("p-damaged"),
                         "q-undamaged":GetScalar("q-undamaged"),
                         "q-damaged":GetScalar("q-damaged")
                         })

def get_data_all(folder,frame_number):
    regex = re.compile(r'sim(_\d+)?_{}'.format(frame_number))
    files = list(filter(regex.search,os.listdir(folder)))
    subframes = [get_data(folder + "/" + f) for f in files]
    df = pd.concat(subframes)
    return df



import subprocess

plt.rc('font', family='serif', serif='Times')
# plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=8)
plt.rc('ytick', labelsize=8)
plt.rc('axes', labelsize=8)
width = 3.487
height = width / 1.618


#output_regex = re.compile("output-\d+")
#output_list = list(filter(output_regex.match,os.listdir()))
#print(output_list)
# output_dir = "./{}/".format(output_list[int(input())])
chalk_dir ="./" 
output_regex = re.compile("output*")
output_list = list(filter(output_regex.match,os.listdir(chalk_dir)))
for i,out in enumerate(output_list):
    print("{}: {}".format(i,out))
output_dir = chalk_dir + "./{}/".format(output_list[int(input())]) 

#output_dir = "..//output/"
#output_dir = "./output/"
#xlim = [0.06,0.12+0.08]
#ylim = [0,0.10]
with open(output_dir+"settings.json") as f:
    json_settings = json.load(f)
    #print("Water level:{}".format(json_settings["OCEAN-HEIGHT"]))
    #print("Domain size:{}".format(json_settings["DOMAIN-SIZE"]))
    h = json_settings["RESOLUTION"]
    xlim = [0,json_settings["DOMAIN-SIZE"][0]]
    ylim = [0,json_settings["DOMAIN-SIZE"][1]]
    #ylim[0] = 20

water_height = 0

ice_height = 200

plt.close("all")

files = os.listdir(output_dir)
#Grab all files that are unique to first rank
finalcsv = re.compile("sim(_0+)?_\d*\.vtk")
files_csvs = list(filter(finalcsv.match,files))

#finalcsv = re.compile("sim_0+_.*\.vtk")

framenumber_regex = re.compile("\d+")
files_csvs = list(map(lambda x: framenumber_regex.findall(x)[-1],files_csvs))
files_csvs.sort(key=int)
#files_csvs = list(map(lambda x: "sim_{}.vtk".format(x), files_csvs))
print("files: {}".format(files_csvs))
dt = 1e4/60
time = []
timesteps = pd.read_csv(output_dir+"timesteps.csv")
max_stress = []
damage = []
full_data = []


def get_nearest_mp(df,x,y):
    p = np.array([df["coord_x"], df["coord_y"]]).transpose()
    x = np.array([x,y])
    dist = np.power(p-x,2).sum(axis=-1)
    mp_index = np.argmin(dist)
    return mp_index,df["unique-id"].values[mp_index]

mp_index = 0
def get_plot(i):
    global mp_index
    fname = files_csvs[i]
    plt.clf()
    df = get_data_all(output_dir,fname)
    print("Plot frame {}".format(i),flush=True)
    ax = fig.add_subplot(111,aspect="equal")
    #df = full_data[i]
    patch_list=[]
    # patch = Rectangle(xy=(0,0) ,width=xlim[1], height=water_height,color="blue")
    # patch_sea = [patch]
    # ps = PatchCollection(patch_sea)
    # ax.add_collection(ps)

    for a_x, a_y,lx,ly,damage in zip(df["coord_x"],
                                     df["coord_y"],
                                     df["lx"],
                                     df["ly"],
                                     df["damage"]):
        patch = Rectangle(
            xy=(a_x-lx/2, a_y-ly/2) ,width=lx, height=ly)
        patch_list.append(patch)
    p = PatchCollection(patch_list, cmap=cm.plasma, alpha=1)
    colour_data = np.zeros(len(df["damage"])) * 0
    mp_index,uid = get_nearest_mp(df,15.5,4)
    print(mp_index)
    print(uid)
    print("Loc: {} {}".format(df["coord_x"].values[mp_index],df["coord_y"].values[mp_index]))
    colour_data[mp_index]=1.0
    p.set_array(colour_data)
    p.set_clim([0,1.0])
    ax.add_collection(p)
    #fig.colorbar(p,location="bottom",label=data_name)

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    plt.title("Frame {} - T {:.2f}s - {} - E: {} - OOBF: {}".format(i,
        timesteps["time"].iloc[i],
        timesteps["step-type"].iloc[i],
        timesteps["energy"].iloc[i],
        timesteps["oobf"].iloc[i]
        ))
    d_t = df["damage-tension"].values[mp_index]
    d_s = df["damage-shear"].values[mp_index]
    d_c = df["damage-compression"].values[mp_index]
    print("Damages: T {} S {} C {}".format(d_t,
                                           d_s,
                                           d_c))
    plt.figure()
    plot_yield_surface("undamaged",0,0,0)
    plot_yield_surface("damaged",d_t,d_s,d_c)
    p,q = -df["p-damaged"].values[mp_index],df["q-damaged"].values[mp_index]
    plt.scatter(p,q,label="Current stress")
    print("Current stress, P {} Q {}".format(p,q))
    print("Residual yield: {}".format(
        (180/np.pi)*np.arctan(((1-d_s)/(1-d_c))*np.tan(30*np.pi/180))))
    # plt.savefig("outframes/frame_{:05}.png".format(i))
    # plt.clf()
    mp_index = uid


def get_data_over_time(mp_uindex,max_time=-1):
    df_0 = get_data_all(output_dir,files_csvs[0])
    length = len(files_csvs)
    if not max_time == -1:
        length = min(max_time,length)
    else:
        max_tme = length

    p_undamaged = np.zeros(length)
    q_undamaged = np.zeros(length)
    p_damaged = np.zeros(length)
    q_damaged = np.zeros(length)
    coord_x = np.zeros(length)
    coord_y = np.zeros(length)
    damage = np.zeros(length)
    s_1 = np.zeros(length)
    s_3 = np.zeros(length)
    su_1 = np.zeros(length)
    su_3 = np.zeros(length)
    for i,fname in enumerate(files_csvs):
        if i > length-1:
            break
        df = get_data_all(output_dir,fname)
        print(i)
        mp_index = np.where(df["unique-id"]==mp_uindex)[0]
        p_undamaged[i] = df["p-undamaged"].values[mp_index]
        q_undamaged[i] = df["q-undamaged"].values[mp_index]
        p_damaged[i] = df["p-damaged"].values[mp_index]
        q_damaged[i] = df["q-damaged"].values[mp_index]
        coord_x[i] = df["coord_x"].values[mp_index]
        coord_y[i] = df["coord_y"].values[mp_index]
        s_1[i] = df["s_1"].values[mp_index]
        s_3[i] = df["s_3"].values[mp_index]
        su_1[i] = df["su_1"].values[mp_index]
        su_3[i] = df["su_3"].values[mp_index]
        damage[i] = df["damage"].values[mp_index]
    filter_width = 5
    return pd.DataFrame({"damage":     uniform_filter1d(damage     ,filter_width),
                         "p_undamaged":uniform_filter1d(p_undamaged,filter_width),
                         "q_undamaged":uniform_filter1d(q_undamaged,filter_width),
                         "p_damaged":  uniform_filter1d(p_damaged  ,filter_width),
                         "q_damaged":  uniform_filter1d(q_damaged  ,filter_width),
                         "coord_x":  uniform_filter1d(coord_x  ,filter_width),
                         "coord_y":  uniform_filter1d(coord_y  ,filter_width),
                         "s_1":  uniform_filter1d(s_1,filter_width),
                         "s_3":  uniform_filter1d(s_3,filter_width),
                         "su_1":  uniform_filter1d(su_1,filter_width),
                         "su_3":  uniform_filter1d(su_3,filter_width),
                         })





current_frame = 0
max_frame = len(files_csvs)-1
def replot():
    get_plot(current_frame,files_csvs[current_frame])
    fig.canvas.draw()



data_name = "damage"
data_map = {
        "p":"plastic_strain",
        "d":"damage",
        "y":"damage-ybar-scaled",
        "u":"sig_xx",
        "o":"mpi-index",
        "q":"split-depth",
        "x":"disp_x"
        }

fig = plt.figure()#plt.figure(figsize=(16,9),dpi=200)
get_plot(90)
plt.show()
#mp_index = 13137.0 #15468
df_time = get_data_over_time(mp_index,100)


plt.figure()
plot_yield_varsurface_mc("undamaged",50,1)
#plot_yield_surface_mc_damaged("damaged")
# filter_size = 10
tau_m = 0.5*(df_time["s_1"].values - df_time["s_3"].values) 
sig_m = 0.5*(df_time["s_1"].values + df_time["s_3"].values)
plt.scatter(-sig_m,tau_m,c=df_time["damage"].values)#

#plt.figure()
#plot_yield_surface("undamaged",0,0,0)
#filter_size = 10
#plt.plot(-df_time["p_damaged"].values,df_time["q_damaged"].values)#,c=df_time["damage"].values
#
#plt.figure()
#plt.title("Undamage ss")
#plot_yield_surface_undamaged("undamaged",0,0,0)
#plt.plot(-df_time["p_undamaged"].values,df_time["q_undamaged"].values)#,c=df_time["damage"].values
#
#
fig = plt.figure()
ax = plt.gca()
def animate(i):
    ax.cla() # clear the previous image
    print(i)
    plot_yield_surface("undamaged",0,0,0)
    #plt.plot(-uniform_filter1d(df_time["p_damaged"].values[:i],filter_size),uniform_filter1d(df_time["q_damaged"].values[:i],filter_size))
    plt.scatter(-df_time["p_damaged"].values[:i],df_time["q_damaged"].values[:i],c=df_time["damage"].values[:i])
   #plt.plot(-df_time["p_undamaged"].values[:i],df_time["q_undamaged"].values[:i])#,c=df_time["damage"].values
#
#def animate_undamaged(i):
#    ax.cla() # clear the previous image
#    print(i)
#    plot_yield_surface_undamaged("undamaged",0,0,0)
#    plt.scatter(-df_time["p_undamaged"].values[:i],df_time["q_undamaged"].values[:i],c=df_time["damage"].values[:i])
#
def animate_coord(i):
    ax.cla() # clear the previous image
    plt.xlim([0,30])
    plt.ylim([0,16])
    plt.scatter(df_time["coord_x"].values[:i],df_time["coord_y"].values[:i],c=df_time["damage"].values[:i])
#
def animate_mc_undamaged(i):
    ax.cla() # clear the previous image
    print(i)
    plot_yield_surface_mc("undamaged",0,0,0)
    # plot_yield_surface_mc_damaged("damaged")
    tau_m = 0.5*(df_time["su_1"].values - df_time["su_3"].values) 
    sig_m = 0.5*(df_time["su_1"].values + df_time["su_3"].values)
    plt.scatter(-sig_m[:i],tau_m[:i],c=df_time["damage"].values[:i])#
#anim = animation.FuncAnimation(fig, animate_mc_undamaged, frames = len(df_time) + 1, interval = 1, blit = False)

def animate_mc(i):
    ax.cla() # clear the previous image
    print(i)
    plot_yield_varsurface_mc("undamaged",50,1)
    # plot_yield_surface_mc_damaged("damaged")
    tau_m = 0.5*(df_time["s_1"].values - df_time["s_3"].values) 
    sig_m = 0.5*(df_time["s_1"].values + df_time["s_3"].values)
    plt.scatter(-sig_m[:i],tau_m[:i],c=df_time["damage"].values[:i])#
#anim = animation.FuncAnimation(fig, animate_mc, frames = len(df_time) + 1, interval = 1, blit = False)

#anim = animation.FuncAnimation(fig, animate, frames = len(df_time) + 1, interval = 1, blit = False)
#anim = animation.FuncAnimation(fig, animate_coord, frames = len(df_time) + 1, interval = 1, blit = False)
plt.show()
