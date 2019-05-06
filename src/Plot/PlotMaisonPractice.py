#
# PYTHON for FEM DUMMIES 18-19
# Projet "tsunami"
#
# Visualisation des résultats : OpenGL old-fashion
# Une implémentation avec des shaders sera présentée en S11
#
#  Vincent Legat
#
# -------------------------------------------------------------------------
# 

from pygame import *
import pygame
from OpenGL.GL   import *
from OpenGL.GLU  import *
from OpenGL.GLUT import *

import numpy as np
import tsunami
import sys
from tsunami import *
#theFile = "PacificTriangleTiny.txt"
theFile = "stupid.txt"
R = 6371220
#theMesh = Mesh(theFile, R)
#theEdges = Edges(theMesh)
#vertices = theMesh.xyz
#edges = theEdges.edges[0:2]
#vertices = ((-1,-1,0),(1,-1,0),(1,1,0),(-1,1,0))
#edges = ((0,1),(1,2),(2,3),(3,0))
vertices= (
    (1, -1, -1),
    (1, 1, -1),
    (-1, 1, -1),
    (-1, -1, -1),
    (1, -1, 1),
    (1, 1, 1),
    (-1, -1, 1),
    (-1, 1, 1)
    )
edges = (
    (0,1),
    (0,3),
    (0,4),
    (2,1),
    (2,3),
    (2,7),
    (6,3),
    (6,4),
    (6,7),
    (5,1),
    (5,4),
    (5,7)
    )



def Cube():
    glBegin(GL_LINES)
    for edge in edges:
        for vertex in edge:
            glVertex3fv(vertices[vertex])
    glEnd()
    
def reshape(width, height):   
  global theRatio
  
  glViewport(0,0,width,height)
  theRatio = width / float(height)
  glutPostRedisplay()
# -------------------------------------------------------------------------
def special(symbol,x,y):
  global translationHorizontal, translationVertical
  if symbol == GLUT_KEY_LEFT : 
    translationHorizontal += 3
  elif symbol == GLUT_KEY_RIGHT : 
    translationHorizontal -= 3
  elif symbol == GLUT_KEY_UP :        
    translationVertical -= 3
  elif symbol == GLUT_KEY_DOWN :
    translationVertical += 3
  else:
    return
  glutPostRedisplay()
  
def draw():  
    global translationHorizontal, translationVertical
    glClearColor(0.9, 0.9, 0.5, 1.0) # du jaune, 1.0 est la transparence
    glClear(GL_COLOR_BUFFER_BIT)
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    glTranslatef(0.2*translationHorizontal,14.0,0.0);
    glTranslatef(0.0,-1.0,0.2*translationVertical);
    glBegin(GL_TRIANGLES);
    glutSwapBuffers() 

# -------------------------------------------------------------------------


translationHorizontal = 0; translationVertical= 0; theRatio = 0.8  
glutInit(sys.argv)
glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH)
# Full screen
glutInitWindowSize(glutGet(GLUT_SCREEN_WIDTH), glutGet(GLUT_SCREEN_HEIGHT))
window = glutCreateWindow("My Glut")
glutDisplayFunc(draw)
#glutIdleFunc(draw)
glutMainLoop()





















