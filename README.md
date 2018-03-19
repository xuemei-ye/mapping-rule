# Which Mapping Rule in the Fireworks Algorithm is Better for Large Scale Optimization

IEEE CEC 2018, by Xuemei Ye, Junzhi Li, Bo Xu and Ying Tan

This is the code for implementing the BBFWA algorithm to investigate the preformence of five mapping rule presented in the paper: Which Mapping Rule in the Fireworks Algorithm is Better for Large Scale Optimization. bbfwa_drawPoint.py is used to draw abridged general view of five mapping rules.

### Experiment

The experiments designed 9 groups experiments on the selected 9 benchmark functions, the dimensions of each experiment is 100, 400, 700, 1000 respectively.The parameters of each run time: n = 30, Cr = 0.9, Ca = 1.2, Evaluation times: 30000.  Each line in the figure is the mean of 50 runs, the range of each line presents the confidence interval. e.g. Cigar function'result :

![Cigar Function](https://github.com/xuemei-ye/mapping-rule/blob/master/Cigar.PNG)

### Five Mapping Rules

![Mapping Rule](https://github.com/xuemei-ye/mapping-rule/blob/master/Mapping%20Rule.PNG)

**Read more :**

[Fireworks Algorithm](http://www.cil.pku.edu.cn/research/fa/)

[Note](https://zhuanlan.zhihu.com/p/34699945)

[Relevant paper and code](http://www.cil.pku.edu.cn/research/fwa/resources/index.html)
