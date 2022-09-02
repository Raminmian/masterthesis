
import pandas as pd

import os
data = pd.read_csv('padloc_summary/system_info.txt', sep='\t')
df=data[['system',"yaml.name"]]
systems=list(set(df.iloc[:,0]))

Dict = {}
for i in range(len(df)):
    Dict[df.iloc[i,1]]=df.iloc[i,0]
    
def read_output(file,spename):
    with open(file) as f:
        pro= f.readlines()
    pro1=[pro[i].split("\n")[0] for i in range(len(pro))]
    pro1.append("GCF_lastline")
    index=[i for i in range(len(pro1)) if "GCF" in pro1[i]]
    alldf= pd.DataFrame(columns = systems)
    
    zero=[0]*alldf.shape[1]
    for k in range(len(index)-1):
        name=pro1[index[k]]
        output=pro1[index[k]+1:index[k+1]]
        uout=list(set(output))
        sysout=[uout[i].split(",")[1] for i in range(len(uout))]
        sys=[Dict[sysout[i]] for i in range(len(sysout))]
        alldf.loc[name]=zero
        for j in alldf.columns:
            alldf.loc[name][j]=sys.count(j)
    alldf['spename']=[spename]*alldf.shape[0]
    return alldf


#allaa=read_output('padloc_summary/Acinetobacterbaumannii_summary/systems.txt')

spe1=snakemake.input[0].split("_summary")[1].split("/")[1]
pse=read_output("{}/systems.txt".format(snakemake.input[0]),spe1)
pse = pse.reindex(sorted(pse.columns), axis=1)
pse.to_csv(snakemake.output[0])