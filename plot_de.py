import matplotlib.pyplot as plt
import pandas as pd
import json
import os,re

chalk_dir ="/nobackup/rmvn14/ham-chalk/" 
output_regex = re.compile("output*")                                 
output_list = list(filter(output_regex.match,os.listdir(chalk_dir)))
output_list.sort()
for i,out in enumerate(output_list):                                    
    print("{}: {}".format(i,out))                                       

output_dir = chalk_dir + "./{}/".format(output_list[int(input())]) 
df = pd.read_csv(output_dir+"timesteps.csv")

with open(output_dir+"./settings.json") as f:
    json_settings = json.load(f)
    h = json_settings["RESOLUTION"]
    xlim = [0,json_settings["DOMAIN-SIZE"][0]]
    ylim = [0,json_settings["DOMAIN-SIZE"][1]]
    thresh_energy = [0,json_settings["DOMAIN-SIZE"][1]]
    thresh_energy = json_settings["CRITERIA-ENERGY"]
    thresh_oobf = json_settings["CRITERIA-OOBF"]
    #ylim[0] = 20

thresh_scale = 0.1
#thresh_energy = 1e-1
#thresh_oobf = 1e-1

print(df)
plt.plot(df["time"],df["damage"]/df["damage"].max() ,label="Damage")
plt.plot(df["time"],thresh_scale*df["energy"]/thresh_energy ,label="Energy")
plt.plot(df["time"],thresh_scale*df["oobf"]/thresh_oobf ,label="OOBF")
#plt.plot(df["time"],df["mass"]/thresh_oobf ,label="Mass")
plt.axhline(thresh_scale,c="green",ls="--")
for i in range(len(df)-1):
    if df["step-type"].iloc[i] != df["step-type"].iloc[i+1]:
        ##Transition found
        x = df["time"].values[i+1]
        plt.axvline(x)
        plt.text(x+1,0,'Transition to {}'.format(df["step-type"].values[i+1]),rotation=90)
plt.legend()
plt.savefig("de.pdf")
plt.show()
