
# coding: utf-8

# Bare bones fireworks algorithm
# https://www.sciencedirect.com/science/article/pii/S1568494617306609
# by Junzhi Li
# ljz@pku.edu.cn

import numpy as np
from numpy import *
import matplotlib.pyplot as plt
import matplotlib.patches as patches


class BBFWA:

    def __init__(self,dim,maxeval,n,Cr,Ca,mapping_rule,plot_episode,run_episode,s):
        self.dim = dim
        self.maxeval = maxeval
        self.n = n
        self.Cr = Cr
        self.Ca = Ca
        self.mapping_rule = mapping_rule
        self.plot_episode = plot_episode
        self.run_episode = run_episode
        self.s = s

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

    def rastrigin(self,individual,r):
        individual = individual - r
        return 10 * len(individual) + sum(gene * gene - 10 * \
                                          cos(2 * pi * gene) for gene in individual)
    def drawpoint(self,s,x,title):

        fig, ax = plt.subplots()
        ax.scatter(s[0, :], s[1, :], color='#79b220')
        ax.scatter(x[0], x[1], color="#e06064")
        # ax.scatter(88,80,color = '#ffdb19')
        # ax.set_title(title)
        plt.plot((100, 100), (100, -100), color='#a54e9c')
        plt.plot((-100, -100), (100, -100), color='#a54e9c')
        plt.plot((-100, 100), (100, 100), color='#a54e9c')
        plt.plot((-100, 100), (-100, -100), color='#a54e9c')
        # plt.xlim(-300, 300)
        # plt.ylim(-300, 300)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        plt.axis('off')
        for ext in ["png", "pdf","eps"]:
            print("saving picture.%s")
            plt.savefig('F:/Fireworks/BBFWA/mapping_reslut/point_image/' + title + ".%s" % (ext,))
        plt.show()

    def optimal(self,r):
        lb = -100 * ones((self.dim,1))
        ub = 100 * ones((self.dim,1))
        A = 50 * ones((self.dim,1)) #ub - lb
        x = np.array([[80],[80]])
        isdraw = True
        fx = self.sphere(x,r)
        result = []
        result.append(fx[0])
        eval = 1
        while eval < self.maxeval:
            # s = (random.rand(self.dim,self.n) * 2 - 1) * tile(A,(1,self.n)) + tile(x,(1,self.n))
            temps = (random.rand(self.dim, self.n) * 2 - 1) * tile(A, (1, self.n)) + tile(x, (1, self.n))
            # boundary handling
            if self.mapping_rule == 'delete':
                temn = self.n
                j = 0
                if isdraw:
                    self.drawpoint(self.s, x, 'Before boundary handling')
                    isdraw = False
                for i in range(self.dim):
                    while j < temn:

                        if self.s[i,j] > ub[0] or self.s[i,j] < lb[0] :
                            self.s = delete(self.s,j,1)  # 删除第j列
                            temn -= 1
                            j -= 1
                        j += 1
                    j = 0
                self.drawpoint(self.s, x, 'After Delete Mapping')
                isdraw = True

            # 随机映射到爆炸范围内，开始时A很大，需要多几次才能填补第一次爆炸范围外的点，如此做，跟随机的效果区别不大， 所以没有用到文章中。
            if self.mapping_rule == 'explodeIn':
                k, m, l, l2, k2, k3 = 0, 0, 0, 0, 0,0
                incount, incount2, incount3 = [], [],[]
                temps = (random.rand(self.dim, self.n) * 2 - 1) * tile(A, (1, self.n)) + tile(x, (1, self.n))
                temps2 = (random.rand(self.dim, self.n) * 2 - 1) * tile(A, (1, self.n)) + tile(x, (1, self.n))
                temps3 = (random.rand(self.dim, self.n) * 2 - 1) * tile(A, (1, self.n)) + tile(x, (1, self.n))
                while k < self.n:
                    temp_colindex = logical_or(temps[:, k] > ub[0], temps[:, k] < lb[0])
                    if not temp_colindex.any():
                        incount.append(k)
                    k += 1

                while k2 < self.n:
                    temp_colindex2 = logical_or(temps2[:, k2] > ub[0], temps2[:, k2] < lb[0])
                    if not temp_colindex2.any():
                        incount2.append(k2)
                    k2 += 1

                while k3 < self.n:
                    temp_colindex3 = logical_or(temps3[:, k3] > ub[0], temps3[:, k3] < lb[0])
                    if not temp_colindex3.any():
                        incount3.append(k3)
                    k3 += 1

                for j in range(self.n):
                    colindex = logical_or(self.s[:, j] > ub[0], self.s[:, j] < lb[0])
                    if colindex.any():
                        if isdraw:
                            self.drawpoint(self.s, x, 'before_explodeIn')
                            isdraw = False
                        if m < len(incount):
                            self.s[:,j] = temps[:,incount[m]]
                            m += 1
                        elif l < len(incount2):
                            self.s[:, j] = temps2[:, incount2[l]]
                            l += 1
                        elif l2 <len(incount3):
                            self.s[:,j] = temps3[:,incount3[l2]]
                            l2 += 1

                self.drawpoint(self.s, x, 'after_explodeIn')
                isdraw = True

            for i in range(self.dim):
                index = logical_or(self.s[i, :] > ub[i],self.s[i, :] < lb[i])  # print(len(index),len(ub[i])) 300 * 1  ub[i] 100
                ubdex = self.s[i, :] > ub[i]
                lbdex = self.s[i, :] < lb[i]
                if self.mapping_rule == 'random':
                    if isdraw:
                        if ubdex.any() or lbdex.any():
                            self.drawpoint(self.s,x,'before_random')
                            isdraw = False
                    self.s[i, index] = random.rand(1, sum(index)) * (ub[i] - lb[i]) + lb[i]
                    if i == self.dim - 1:
                        self.drawpoint(self.s,x,'After Random Mapping')
                        isdraw = True

                if self.mapping_rule == 'origin':
                    if isdraw :
                        if ubdex.any() or lbdex.any():
                            self.drawpoint(self.s,x,'before_origin')
                            isdraw = False
                    self.s[i, index] = abs(self.s[i, index]) % (ub[i] - lb[i]) + lb[i]
                    if i == self.dim - 1:
                        self.drawpoint(self.s,x,'after_origin')
                        isdraw = True

                if self.mapping_rule == 'boundary':
                    if isdraw :
                        if ubdex.any() or lbdex.any():
                            self.drawpoint(self.s,x,'before_boundary')
                            isdraw = False
                    self.s[i, ubdex] = ub[i]
                    self.s[i, lbdex] = lb[i]
                    if i == self.dim - 1:
                        self.drawpoint(self.s,x,'After Boundary Mapping')
                        isdraw = True

                if self.mapping_rule == 'mirror':
                    if isdraw :
                        if ubdex.any() or lbdex.any():
                            self.drawpoint(self.s,x,'before_mirror')
                            isdraw = False
                    self.s[i, ubdex] = ub[i] - (self.s[i, ubdex] - ub[i])
                    self.s[i, lbdex] = lb[i] + (lb[i] - self.s[i, lbdex])
                    if i == self.dim - 1:
                        self.drawpoint(self.s,x,'after_mirror')
                        isdraw = True

                # if self.mapping_rule == 'mixture': 混合随机镜像，效果不佳，未采用。

            fs = ones((len(self.s[0])))  #300 * 1
            for i in range(len(self.s[0])):  # 对每个点计算 x - r，一列为一个点
                tems = self.s[:, i]
                tems = tems.reshape(self.dim, 1)
                fs[i] = self.sphere(tems, r)

            eval = eval + self.n
            if min(fs) < fx:
                x = self.s[:, argmin(fs)].reshape(self.dim, 1) # x:30*1
                fx = min(fs.tolist())
                A = A * self.Ca
            else:
                A = A * self.Cr
            if self.plot_episode == 1001:
                    if eval <= 30000:
                        result.append(fx)
            else:
                result.append(fx) # 存每轮迭代结果
        print(fx)
        return fx,result




