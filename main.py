import numpy as np
from numpy import *
import matplotlib.pyplot as plt
from bbfwa_class import BBFWA
import math
import csv

dim = 30
maxeval = 30000
n = 30; Cr = 0.9; Ca = 1.2 #parameters
plot_episode = 1001
run_episode = 50
plot_begin = 1

savecat = "F:/Fireworks/BBFWA/mapping_reslut/cigar_result/cigar_result_d"

dellastfx,del_allresult ,ranlastfx,ran_allresult= [],[],[],[] #Save the value of each episode,the last value of 1000 run time.
mirlastfx,mir_allresult,orilastfx = [],[],[]
ori_allresult,boundlastfx,bound_allresult,mixlastfx,mix_allresult= [],[],[],[],[]

del_meanresult,ran_meanresult,mir_meanresult,ori_meanresult,bound_meanresult,no_meanresult= [],[],[],[],[],[]
mixfx, mixresult = [],[]

episode = 0

delbbfwa = BBFWA(dim,maxeval,n,Cr,Ca,'delete',plot_episode,run_episode)
oribbfwa = BBFWA(dim, maxeval, n, Cr, Ca, 'origin', plot_episode, run_episode)
ranbbfwa = BBFWA(dim,maxeval,n,Cr,Ca,'random',plot_episode,run_episode)
mirbbwfa = BBFWA(dim,maxeval,n,Cr,Ca,'mirror',plot_episode,run_episode)
boundbbfwa = BBFWA(dim,maxeval,n,Cr,Ca,'boundary',plot_episode,run_episode)

while episode < run_episode:

    lb = -100 * ones((dim, 1))
    ub = 100 * ones((dim, 1))
    r = random.rand(dim,1) * (ub - lb) + lb # offset sphere function

    delfx, delresult = delbbfwa.optimal(r)
    orifx, oriresult = oribbfwa.optimal(r)
    ranfx, ranresult = ranbbfwa.optimal(r)
    mirfx, mirresult = mirbbwfa.optimal(r)
    boundfx, boundresult = boundbbfwa.optimal(r)

    dellastfx.append(delfx)  # last value of 1000
    del_allresult.append(delresult)
    orilastfx.append(orifx)
    ori_allresult.append(oriresult)
    ranlastfx.append(ranfx)
    ran_allresult.append(ranresult)
    mirlastfx.append(mirfx)
    mir_allresult.append(mirresult)
    boundlastfx.append(boundfx)
    bound_allresult.append(boundresult)

    episode += 1


namelist = ["del_allresult", "ori_allresult", "ran_allresult", "mir_allresult", "bound_allresult"]
savelist = [del_allresult, ori_allresult, ran_allresult, mir_allresult, bound_allresult]
for i in range(len(savelist)):
    filename = savecat + str(dim) + '/' + namelist[i] + '.csv'
    with open(filename, 'wb') as f:
        writer = csv.writer(f)
        writer.writerows(savelist[i])


del_stdresult = np.std(np.array(del_allresult), axis=0)
del_meanresult = np.mean(np.array(del_allresult),axis=0)
ori_stdresult = np.std(np.array(ori_allresult), axis=0)
ori_meanresult = np.mean(np.array(ori_allresult),axis=0)
ran_stdresult = np.std(np.array(ran_allresult), axis=0)
ran_meanresult = np.mean(np.array(ran_allresult),axis=0)
mir_stdresult = np.std(np.array(mir_allresult), axis=0)
mir_meanresult = np.mean(np.array(mir_allresult),axis=0)
bound_stdresult = np.std(np.array(bound_allresult), axis=0)
bound_meanresult = np.mean(np.array(bound_allresult),axis=0)


fig = plt.figure()
plt.xlabel('Iterations')
plt.ylabel('Best Fitness')

plt.plot(arange(plot_begin,plot_episode),del_meanresult,label='delete',color = '#0048c3')  #blue
plt.plot(arange(plot_begin,plot_episode),ori_meanresult,label='modular',color = "#ffa500")  # orange
plt.plot(arange(plot_begin,plot_episode),ran_meanresult,label='random',color = '#ff0000') # red
plt.plot(arange(plot_begin,plot_episode),mir_meanresult,label='mirror',color = '#008000') # green
plt.plot(arange(plot_begin,plot_episode),bound_meanresult,label='boundary',color = '#800080') # purple

plt.legend(loc='upper right')
filename = 'one'+'_d' + str(dim) + '_n' + str(n)

for ext in ["png", "pdf", "svg", "eps"]:
    print("saving picture.%s" % (ext,))
    plt.savefig(savecat + str(dim) + '/' + filename + ".%s" % (ext,)) #bbox_inches="tight"



fig = plt.figure()
plt.xlabel('Iterations')
plt.ylabel('Best Fitness')

plt.plot(arange(plot_begin,plot_episode),del_meanresult,label='delete',color = '#0048c3')  #blue
plt.plot(arange(plot_begin,plot_episode),ori_meanresult,label='modular',color = "#ffa500")  # orange
plt.plot(arange(plot_begin,plot_episode),ran_meanresult,label='random',color = '#ff0000') # red
plt.plot(arange(plot_begin,plot_episode),mir_meanresult,label='mirror',color = '#008000') # green
plt.plot(arange(plot_begin,plot_episode),bound_meanresult,label='boundary',color = '#800080') # purple


plt.fill_between(arange(plot_begin,plot_episode),np.array(del_meanresult + del_stdresult),
                 np.array(del_meanresult - del_stdresult),color = "#7fa3e1")#f5dfdd
plt.fill_between(arange(plot_begin,plot_episode),np.array(ori_meanresult + ori_stdresult),
                 np.array(ori_meanresult - ori_stdresult),color = "#ffe4b2")
plt.fill_between(arange(plot_begin,plot_episode),np.array(ran_meanresult + ran_stdresult),
                 np.array(ran_meanresult - ran_stdresult),color = "#ffb2b2")
plt.fill_between(arange(plot_begin,plot_episode),np.array(mir_meanresult + mir_stdresult),
                 np.array(mir_meanresult - mir_stdresult),color = "#99cc99")
plt.fill_between(arange(plot_begin,plot_episode),np.array(bound_meanresult + bound_stdresult),
                 np.array(bound_meanresult - bound_stdresult),color = "#d8b2d8")


plt.legend(loc='upper right')
filename = 'd' + str(dim) + '_n' + str(n)

for ext in ["png", "pdf", "svg", "eps"]:
    print("saving picture.%s" % (ext,))
    plt.savefig(savecat + str(dim) + '/' + filename + ".%s" % (ext,)) #bbox_inches="tight"


del_logallresult = map(log,del_allresult)
del_stdlogresult = np.std(np.array(del_logallresult), axis=0)
del_meanlogresult = np.mean(np.array(del_logallresult),axis=0)
delup_stdlogresult = list(map(lambda  x:x[0]+x[1],zip(del_meanlogresult,del_stdlogresult)))
deldown_stdlogresult = list(map(lambda  x:x[0]-x[1],zip(del_meanlogresult,del_stdlogresult)))

ori_logallresult = map(log,ori_allresult)
ori_stdlogresult = np.std(np.array(ori_logallresult), axis=0)
ori_meanlogresult = np.mean(np.array(ori_logallresult),axis=0)
oriup_stdlogresult = list(map(lambda  x:x[0]+x[1],zip(ori_meanlogresult,ori_stdlogresult)))
oridown_stdlogresult = list(map(lambda  x:x[0]-x[1],zip(ori_meanlogresult,ori_stdlogresult)))


ran_logallresult = map(log,ran_allresult)
ran_stdlogresult = np.std(np.array(ran_logallresult), axis=0)
ran_meanlogresult = np.mean(np.array(ran_logallresult),axis=0)
ranup_stdlogresult = list(map(lambda  x:x[0]+x[1],zip(ran_meanlogresult,ran_stdlogresult)))
randown_stdlogresult = list(map(lambda  x:x[0]-x[1],zip(ran_meanlogresult,ran_stdlogresult)))


mir_logallresult = map(log,mir_allresult)
mir_stdlogresult = np.std(np.array(mir_logallresult), axis=0)
mir_meanlogresult = np.mean(np.array(mir_logallresult),axis=0)
mirup_stdlogresult = list(map(lambda  x:x[0]+x[1],zip(mir_meanlogresult,mir_stdlogresult)))
mirdown_stdlogresult = list(map(lambda  x:x[0]-x[1],zip(mir_meanlogresult,mir_stdlogresult)))


bound_logallresult = map(log,bound_allresult)
bound_stdlogresult = np.std(np.array(bound_logallresult), axis=0)
bound_meanlogresult = np.mean(np.array(bound_logallresult),axis=0)
boundup_stdlogresult = list(map(lambda  x:x[0]+x[1],zip(bound_meanlogresult,bound_stdlogresult)))
bounddown_stdlogresult = list(map(lambda  x:x[0]-x[1],zip(bound_meanlogresult,bound_stdlogresult)))


fig = plt.figure()
plt.xlabel('Iterations')
plt.ylabel('Log Best Fitness')

plt.plot(arange(plot_begin,plot_episode),del_meanlogresult,label='delete',color = '#0048c3')
plt.plot(arange(plot_begin,plot_episode),ori_meanlogresult,label='modular',color = '#ffa500')
plt.plot(arange(plot_begin,plot_episode),ran_meanlogresult,label='random',color = '#ff0000')
plt.plot(arange(plot_begin,plot_episode),mir_meanlogresult,label='mirror',color = '#008000')
plt.plot(arange(plot_begin,plot_episode),bound_meanlogresult,label='boundary',color = '#800080')

plt.legend(loc='upper right')

filename = 'one_log' + '_d' + str(dim) + '_n' + str(n)
for ext in ["png", "pdf", "svg", "eps"]:
    print("saving picture.%s" % (ext,))
    plt.savefig(savecat + str(dim) + '/' + filename + ".%s" % (ext,))


fig = plt.figure()
plt.xlabel('Iterations')
plt.ylabel('Log Best Fitness')

plt.plot(arange(plot_begin,plot_episode),del_meanlogresult,label='delete',color = '#0048c3')
plt.plot(arange(plot_begin,plot_episode),ori_meanlogresult,label='modular',color = '#ffa500')
plt.plot(arange(plot_begin,plot_episode),ran_meanlogresult,label='random',color = '#ff0000')
plt.plot(arange(plot_begin,plot_episode),mir_meanlogresult,label='mirror',color = '#008000')
plt.plot(arange(plot_begin,plot_episode),bound_meanlogresult,label='boundary',color = '#800080')


plt.fill_between(arange(plot_begin,plot_episode),delup_stdlogresult,
                 deldown_stdlogresult,color = "#7fa3e1")
plt.fill_between(arange(plot_begin,plot_episode),oriup_stdlogresult,
                 oridown_stdlogresult,color = "#ffe4b2")
plt.fill_between(arange(plot_begin,plot_episode),ranup_stdlogresult,
                 randown_stdlogresult,color = "#ffb2b2")
plt.fill_between(arange(plot_begin,plot_episode),mirup_stdlogresult,
                 mirdown_stdlogresult,color = "#99cc99")
plt.fill_between(arange(plot_begin,plot_episode),boundup_stdlogresult,
                 bounddown_stdlogresult,color = "#d8b2d8")

plt.legend(loc='upper right')

filename = 'log' + '_d' + str(dim) + '_n' + str(n)
for ext in ["png", "pdf", "svg", "eps"]:
    print("saving picture.%s" % (ext,))
    plt.savefig(savecat + str(dim) + '/' + filename + ".%s" % (ext,))

plt.show()

