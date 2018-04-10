import csv
import numpy as np
from scipy import stats


dim = 1000
# functiList = ['sphere','cigar','discus','Ellipse','Rosenbrock','Bohachevsky','Griewank','Step','Tablet']
savecat = "F:/Fireworks/BBFWA/mapping_reslut/Tablet_result/Tablet_result_d"

del_allresult ,ran_allresult,mir_allresult ,ori_allresult,bound_allresult= [],[],[],[],[]
namelist = ["del_allresult", "ori_allresult", "ran_allresult", "mir_allresult", "bound_allresult"]
savelist = [del_allresult, ori_allresult, ran_allresult, mir_allresult, bound_allresult]

filename = savecat + str(dim) + '/' + namelist[0] + '.csv'
dcount,ocount,rcount,mcount,bcount=0,0,0,0,0

def finddupl(lst):
    """找出重复的项
    """
    exists, dupl = [], []
    for item in lst:
        if item[0] in exists:
            dupl.append(item[0])
        else:
            exists.append(item[0])

    return dupl

for i in range(len(savelist)):
    filename = savecat + str(dim) + '/' + namelist[i] + '.csv'
    with open(filename, 'r') as f:
        print(filename)
        reader = csv.reader(f)
        for row in reader:
            temrow = []
            for j in range(len(row)):
                if row[j].find('e') != -1:
                    row[j] = row[j].replace('e', 'E')
                if row[j].find(']') != -1 and row[j].find('[') != -1:
                    row[j] = row[j].lstrip('[')
                    row[j] = row[j].strip(']')
                temrow.append(float(row[j]))
            savelist[i].append(temrow)

del_meanresult = np.mean(np.array(del_allresult),axis=0)
ori_meanresult = np.mean(np.array(ori_allresult),axis=0)
ran_meanresult = np.mean(np.array(ran_allresult),axis=0)
mir_meanresult = np.mean(np.array(mir_allresult),axis=0)
bound_meanresult = np.mean(np.array(bound_allresult),axis=0)


pvalue = stats.kruskal(del_meanresult,ori_meanresult,ran_meanresult,mir_meanresult,bound_meanresult)
print(pvalue)


dic_del_meanresult = [(del_meanresult[i], 0) for i in range(len(del_meanresult))]
dic_ori_meanresult = [(ori_meanresult[i],1) for i in range(len(ori_meanresult))]
dic_ran_meanresult = [(ran_meanresult[i], 2) for i in range(len(ran_meanresult))]
dic_mir_meanresult = [(mir_meanresult[i],3) for i in range(len(mir_meanresult))]
dic_bound_meanresult = [(bound_meanresult[i],4) for i in range(len(bound_meanresult))]

all_meanresule = dic_del_meanresult + dic_mir_meanresult + dic_ori_meanresult + dic_ran_meanresult + dic_bound_meanresult
all_meanresule = sorted(all_meanresule, key=lambda tup: tup[0])
rep = finddupl(all_meanresule)
print(rep,len(rep))

for i in range(len(all_meanresule)):
    if all_meanresule[i][1] ==0:
        dcount += i+1
    if all_meanresule[i][1] == 1:
        ocount += i+1
    if all_meanresule[i][1] == 2:
        rcount += i+1
    if all_meanresule[i][1] == 3:
        mcount += i+1
    if all_meanresule[i][1] == 4:
        bcount += i+1

print(dcount,ocount,rcount,mcount,bcount)

H = 12/(5001*5000) * (dcount*dcount/1000 + ocount*ocount/1000  + rcount*rcount/1000+ mcount*mcount/1000+ bcount*bcount/1000)- 3 *(5001)
Ti = 20*20*20 - 20
Hc = H/(1-Ti/(5000*5000-5000))
print(H,Hc)





