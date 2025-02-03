import matplotlib.pyplot as plt
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
chalk_dir ="./"
output_dir ="/mnt/c/Temp/output-paper/"
df = pd.read_csv(output_dir+"timesteps.csv")

with open(output_dir+"./settings.json") as f:
    json_settings = json.load(f)
    h = json_settings["RESOLUTION"]
    xlim = [0,json_settings["DOMAIN-SIZE"][0]]
    ylim = [0,json_settings["DOMAIN-SIZE"][1]]
    thresh_energy = json_settings["CRITERIA-ENERGY"]
    thresh_oobf = json_settings["CRITERIA-OOBF"]
    #ylim[0] = 20

thresh_scale = 0.1
#thresh_energy = 1e-1
#thresh_oobf = 1e-1

fig = plt.figure(figsize=(width,height),dpi=200)
print(df)
offset = 0.1
plt.xlim([00,100])
plt.plot(df["time"].values,df["damage"].values/df["damage"].max() ,label="Damage")
plt.plot(df["time"].values,thresh_scale*df["energy"].values/thresh_energy ,label="Energy")
plt.plot(df["time"].values,thresh_scale*df["oobf"].values/thresh_oobf ,label="OOBF")
plt.axhline(thresh_scale,c="green",ls="--")
# if True:
#     for i in range(len(df)-1):
#         if df["step-type"].iloc[i] != df["step-type"].iloc[i+1]:
#             ##Transition found
#             x = df["time"].values[i+1]
#             plt.axvline(x)
#             plt.text(x+offset,0,'Transition to {}'.format(df["step-type"].values[i+1]),rotation=90)
plt.xlabel("Time (s)")
plt.legend()
plt.tight_layout()
plt.savefig("de.pdf")
plt.show()
