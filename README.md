**Dynamic-Movement-Primitives-and-Imitation-Learning**

Dynamic movement primitives (DMPs) are a method of trajectory control/planning from Prof.Stefan Schaal’s lab. 

Complex movements have long been thought to be composed of sets of primitive action ‘building blocks’ executed in sequence and \ or in parallel, and DMPs are a proposed mathematical formalization of these primitives. The difference between DMPs and previously proposed building blocks is that each DMP is a nonlinear dynamical system. The basic idea is that you take a dynamical system with well specified, stable behavior and add another term that makes it follow some interesting trajectory as it goes about its business. 

The DMP differential equations (Transformation System, Canonical System, Non-linear Function) realize a general way of generating point-to-point movements.

Imitation learning using linear regression is performed to compute the weight factor W from a demonstrated trajectory dataset, given by a teacher. The quality of the imitation is evaluated by comparing the training data with the data generated by the DMP.



![Test Image 4](https://github.com/suryakiranmg/Dynamic-Movement-Primitives-and-Imitation-Learning-Robotics/blob/master/imitation_learn.png)


**References:**

http://www-clmc.usc.edu/Research/MotorPrimitives 

http://www-clmc.usc.edu/Research/ImitationLearning

https://studywolf.wordpress.com/2013/11/16/dynamic-movement-primitives-part-1-the-basics/


