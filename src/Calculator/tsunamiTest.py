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

import numpy as np
#import tsunamiExperimental as tsunami
import tsunamiBoucle as tsunami
#import tsunami as tsunami
#import tsunamiv2 as tsunami
#import SauvegardeV0 as tsunami
#
#
# -1- Lecture des données
#
R = 6371220
#theMeshFile = "/Users/romaingraux/Library/Mobile Documents/com~apple~CloudDocs/Professionel/EPL/Q4/MAP/Elements finis/Projet/Mesh/Tiny.txt"
theMeshFile = "/Users/romaingraux/Library/Mobile Documents/com~apple~CloudDocs/Professionel/EPL/Q4/MAP/Elements finis/Projet/Mesh/PacificTriangleFine.txt"
[nNode,X,Y,H,nElem,elem] = tsunami.readMesh(theMeshFile)
print(" == Number of elements : %d " % nElem)
print(" == Number of nodes    : %d " % nNode)

#
# -2- On impose la condition initiale et on sauvegarde les élévations 
#     linéaires discontinues dans un fichier de résultats :-)
#
#     Observer qu'on a créé des coordonnées discontinues pour pouvoir
#     évaluer la condition initiale.
#
#theMesh = tsunami.Mesh(theMeshFile,R)
x = np.zeros([nElem,3])
y = np.zeros([nElem,3])
for iElem in range(nElem):
  nodes  = elem[iElem]
  x[iElem][:] = X[nodes]
  y[iElem][:] = Y[nodes] 
E = tsunami.initialConditionOkada(x,y)

theResultFiles = "test-%06d.txt"
tsunami.writeResult(theResultFiles,0,E)

#
# -3- Calcul du tsunami en relisant les conditions initiales 
#     dans le fichier qu'on vient juste d'écrire :-)
#

U = np.zeros([nElem,3])
V = np.zeros([nElem,3])
E = tsunami.readResult(theResultFiles,0,nElem)
dt = 5; nIter = 2; nSave = 1000
[U,V,E] = tsunami.compute(theMeshFile,theResultFiles,U,V,E,dt,nIter,nSave)

#for iElem in [27,28] :
#  print(" == Elevations for element %d : %14.7e %14.7e %14.7e " % (iElem,*E[iElem][:]) )
