# -------------------------------------------------------------------------
#
# PYTHON for FEM DUMMIES 18-19
# Projet "tsunami"
#
# Script de test
#  Vincent Legat
#
# -------------------------------------------------------------------------
# 
#import mutex
#from threading import Thread, Lock
#import time
import numpy as np
import tsunamiBoucle as tsunami

filename = "../../Mesh/Tiny.txt"
R = 6371220
theTsunami = tsunami.Tsunami(filename, np.zeros((199,3)), np.zeros((199,3)), np.zeros((199,3)))




Ainverse = np.array([18.0,-6.0,-6.0,-6.0,18.0,-6.0,-6.0,-6.0,18.0])
#Ainverse = Ainverse.reshape((1,9))
Zeros = np.array([[0,0,0]])
#print(Ainverse)
x = np.bincount([0,0,0],  weights=np.ones(3), minlength = 9)
Ainverse[:] += x
print(x)
#Ainverse += x

#for i in range(3):
#    Ainverse[Zeros[i]] += 1
#print(Ainverse)

