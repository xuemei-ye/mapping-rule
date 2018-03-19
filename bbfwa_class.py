
# coding: utf-8
# References: Bare bones fireworks algorithm by Junzhi Li 
# https://www.sciencedirect.com/science/article/pii/S1568494617306609
# ljz@pku.edu.cn

import numpy as np
from numpy import *
import matplotlib.pyplot as plt

class BBFWA:

    def __init__(self,dim,maxeval,n,Cr,Ca,mapping_rule,plot_episode,run_episode):
        self.dim = dim
        self.maxeval = maxeval
        self.n = n
        self.Cr = Cr
        self.Ca = Ca
        self.mapping_rule = mapping_rule
        self.plot_episode = plot_episode
        self.run_episode = run_episode

    # sphere function
    def sphere(self,x,r):
        return sum((x - r) * (x - r), 0)

    def cigar(self,individual,r):
        individual = individual - r
        return individual[0] ** 2 + 1e6 * sum(gene * gene for gene in individual)

    def discus(self,individual,r):
        individual = individual - r
        return 1e6 * (individual[0] ** 2) + sum(gene * gene for gene in individual)

    def schaffer(self,individual, r):
        individual = individual - r
        return sum(((x ** 2 + x1 ** 2) ** 0.25 * ((sin(50 * (x ** 2 + x1 ** 2) ** 0.1)) ** 2 + 1.0)
                    for x, x1 in zip(individual[:-1], individual[1:])), 0)

    def rosenbrock(self,individual, r):
        individual = individual - r
        return sum(100 * (x * x - y) ** 2 + (1. - x) ** 2 \
                   for x, y in zip(individual[:-1], individual[1:]))

    def schwefel(self,individual,r):   # ub = 500,lb = -500  two places
        individual = individual - r
        N = len(individual)
        return 418.9828872724339 * N - sum(x * sin(sqrt(abs(x))) for x in individual)

    def step(self,x,r):
        x = x - r
        return sum((floor(x + 0.5)) ** 2)

    def optimal(self,r):
        lb = -100 * ones((self.dim,1))  #将其变为30维
        ub = 100 * ones((self.dim,1))
        A = ub - lb
        x = random.rand(self.dim, 1) * (ub - lb) + lb
        fx = self.cigar(x,r)
        result = []
        result.append(fx[0])
        eval = 1
        while eval < self.maxeval:
            s = (random.rand(self.dim,self.n) * 2 - 1) * tile(A,(1,self.n)) + tile(x,(1,self.n))
            # boundary handling
            if self.mapping_rule == 'delete':
                temn = self.n
                j = 0
                for i in range(self.dim):
                    while j < temn:
                        if s[i,j] > ub[0] or s[i,j] < lb[0] :
                            s = delete(s,j,1)  # 删除第j列
                            temn -= 1
                            j -= 1
                        j += 1
            else:
                for i in range(self.dim):
                    index = logical_or(s[i, :] > ub[i],s[i, :] < lb[i])  # print(len(index),len(ub[i])) 300 * 1  ub[i] 100
                    ubdex = s[i, :] > ub[i]
                    lbdex = s[i, :] < lb[i]
                    if self.mapping_rule == 'random':
                        s[i, index] = random.rand(1, sum(index)) * (ub[i] - lb[i]) + lb[i]
                    if self.mapping_rule == 'origin':
                        s[i, index] = abs(s[i, index]) % (ub[i] - lb[i]) + lb[i]
                    if self.mapping_rule == 'boundary':
                        s[i, ubdex] = ub[i]
                        s[i, lbdex] = lb[i]
                    if self.mapping_rule == 'mirror':
                        s[i, ubdex] = ub[i] - (s[i, ubdex] - ub[i])
                        s[i, lbdex] = lb[i] + (lb[i] - s[i, lbdex])
                    if self.mapping_rule == 'nothing':
                        pass
                    if self.mapping_rule == 'mixture':
                        if eval < 3000:
                            s[i, index] = random.rand(1, sum(index)) * (ub[i] - lb[i]) + lb[i]
                        else:
                            s[i, ubdex] = ub[i] - (s[i, ubdex] - ub[i])
                            s[i, lbdex] = lb[i] + (lb[i] - s[i, lbdex])


            fs = ones((len(s[0])))  #300 * 1
            for i in range(len(s[0])):  # 对每个点计算 x - r，一列为一个点
                tems = s[:, i]
                tems = tems.reshape(self.dim, 1) # 转为列向量
                fs[i] = self.cigar(tems, r)

            eval = eval + self.n
            if min(fs) < fx:
                x = s[:, argmin(fs)].reshape(self.dim, 1) # x:30*1
                fx = min(fs.tolist())
                A = A * self.Ca
            else:
                A = A * self.Cr
            if self.plot_episode == 1001:
                    if eval <= 30000:  # plot_episode change to 502 存前500次的结果
                        result.append(fx)
            else:
                result.append(fx) # 存每轮迭代结果
        print(fx)
        return fx,result




