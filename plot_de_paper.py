import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import pandas as pd
import json
import os,re

plt.style.use("seaborn-paper")
plt.rc('font', family='serif', serif='Times')
# plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=8)
plt.rc('ytick', labelsize=8)
plt.rc('axes', labelsize=8)
scale = 2
width = 3.487*scale
# height = width / 1.400
height = width / 1.618

# chalk_dir ="./-paper-boundary/"
# output_regex = re.compile("output*")
# output_list = list(filter(output_regex.match,os.listdir(chalk_dir)))
# for i,out in enumerate(output_list):
#     print("{}: {}".format(i,out))

# output_dir = chalk_dir + "./{}/".format(output_list[int(input())]) 
#chalk_dir ="./"
chalk_dir ="/nobackup/rmvn14/ham-chalk/" 
output_dir =chalk_dir+"./output-paper-strong/"
df = pd.read_csv(output_dir+"timesteps.csv")

with open(output_dir+"./settings.json") as f:
    json_settings = json.load(f)
    h = json_settings["RESOLUTION"]
    xlim = [0,json_settings["DOMAIN-SIZE"][0]]
    ylim = [0,json_settings["DOMAIN-SIZE"][1]]
    thresh_energy = json_settings["CRITERIA-ENERGY"]
    thresh_oobf = json_settings["CRITERIA-OOBF"]
    tau = 1e1
    #ylim[0] = 20

thresh_scale = 1
thresh = thresh_scale
#thresh_energy = 1e-1
#thresh_oobf = 1e-1

fig = plt.figure(figsize=(width,height),dpi=200)
ax = fig.gca()
print(df)
offset = 0.1
#plt.xlim([00,100])

time = df["time"].values/tau
#time = df["steps"].values
energy = df["energy"].values/thresh_energy
oobf = df["oobf"].values/thresh_oobf
prop_cycle = plt.rcParams['axes.prop_cycle']
colours = prop_cycle.by_key()['color']

ax.plot(time,energy,label="Energy",c=colours[0],ls="-")
ax.plot(time,oobf,label="OOBF",c=colours[1],ls="-")
ax.set_ylabel("Normalised dynamic criteria")
ax.set_xlabel("Normalised time")
#ax.set_xlabel("Step")
ax.set_ylim(0,thresh*4)
ax.axhline(thresh_scale,c="black",ls="--")

#plt.plot(df["time"].values,df["damage"].values/df["damage"].max() ,label="Damage")
#plt.plot(df["time"].values,df["energy"].values ,label="Energy")
#plt.plot(df["time"].values,df["oobf"].values ,label="OOBF")

quasi_point = None
for i in range(len(df)-1):
    if df["step-type"].iloc[i] == "COLLAPSE":
        if quasi_point == None:
            quasi_point = time[i]
    else:
        if quasi_point:
            plt.axvspan(quasi_point,time[i],alpha=0.25,color="black")
            quasi_point = None
if quasi_point:
    plt.axvspan(quasi_point,time[i],alpha=0.25,color="black")

#total_mass = df["mass"].values[0]#800 * 800 * 918
total_mass = 1
ax_damage = ax.twinx()
ax_damage.plot(time,df["damage"]/total_mass,label="Damage",c="black")
ax_damage.set_ylim(bottom=0)

# if True:
#     for i in range(len(df)-1):
#         if df["step-type"].iloc[i] != df["step-type"].iloc[i+1]:
#             ##Transition found
#             x = df["time"].values[i+1]
#             plt.axvline(x)
#             plt.text(x+offset,0,'Transition to {}'.format(df["step-type"].values[i+1]),rotation=90)
lines, labels = ax.get_legend_handles_labels()
lines2, labels2 = ax_damage.get_legend_handles_labels()
ax_damage.legend(lines + lines2, labels + labels2, loc=0)
#plt.legend(["Energy","Force","Damage"])
ax.xaxis.set_major_locator(MultipleLocator(1))
plt.xlim(left=0)
plt.xlim(right=10)
plt.xlabel("Time (s)")
plt.ylabel("Mass-Damage (kg)")
plt.tight_layout()
plt.savefig("de.pdf")
plt.show()
