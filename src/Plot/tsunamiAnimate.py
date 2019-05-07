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


from OpenGL.GL   import *
from OpenGL.GLU  import *
from OpenGL.GLUT import *

import sys
sys.path.insert(0, '/Users/romaingraux/Library/Mobile Documents/com~apple~CloudDocs/Professionel/EPL/Q4/MAP/Elements finis/Projet/src/Calculator/')

import tsunamiBoucle as tsunami

import numpy as np
import time

# -------------------------------------------------------------------------

def draw():  
  global E,theFlagBathymetry,translationHorizontal,translationVertical,zoom,theRatio, degreesHorizontal, degreesVertical, view

  glClearColor( 0.9, 0.9, 0.8, 0.0 );
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(65.0,theRatio,1.0,100.0);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(0.0,1.0,0.0,0.0,20.0,0.0,0.0,0.0,1.0);  
  glTranslatef(0.2*translationHorizontal,3.0,0.2*translationVertical);
  glTranslatef(0.0,0.25*zoom,0.0);
#  glRotatef(degreesVertical,np.cos(degreesHorizontal * np.pi / 180),np.sin(degreesHorizontal * np.pi / 180),0.0);
  if(view == 0):
      glTranslatef(0.0,14.0,0.0);
      glRotatef(degreesHorizontal,0.0,0.0,-1.0);      
#      glRotatef(degreesVertical,np.cos(degreesHorizontal * np.pi / 180),np.sin(degreesHorizontal * np.pi / 180),0.0);
  if(view == 1):
      glTranslatef(0.0,7.0,-4.0);
      glRotatef(116,0.0,0.0,1.0);
  if(view == 2):
      glTranslatef(0.0,9.0,-3.0);
      glRotatef(116,0.0,0.0,1.0);
#  glRotatef(degreesVertical, 0, 1, 0)
#  print(degreesVertical)
#  glRotatef(degreesVertical, np.cos(degreesHorizontal * np.pi / 180),0,np.sin(degreesHorizontal * np.pi / 180) )
#  glRotatef(0.7*translationHorizontal,0.0,0.0,1.0);
#  glRotatef(0.7*translationVertical,1.0,1.0,0.0);
  
  quadratic = gluNewQuadric();         
  gluQuadricNormals(quadratic, GLU_SMOOTH); 
  glColor3f(1.0,1.0,1.0);
  gluSphere(quadratic,5.95,400,200);
 
  n = 9*nElem
  index  = np.ravel(elem)
  colors = np.zeros(n) 
  coord  = np.zeros(n)

  if (theFlagBathymetry == False) :
    value = 10* np.ravel(E)
  else :
    value = H[index]/BathMax
     
#  print(np.clip(value,0,1))
#  value = np.asarray(value)
  value = np.clip(value,a_min = 0.0,a_max = 1.0)        
  colors[0:n:3] = 3.5*(value)*(value)
  colors[1:n:3] = (1-value)*(value)*3.5
  colors[2:n:3] = (1-value)*(1-value)
  
  x = X[index]
  y = Y[index]
  factor = (4*R*R + x*x + y*y)*(R/6)
  coord[0:n:3] = 4*R*R * x / factor
  coord[1:n:3] = 4*R*R * y / factor
  coord[2:n:3] = (4*R*R - x*x - y*y)*R / factor   
 
  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_COLOR_ARRAY);
  glEnableClientState(GL_NORMAL_ARRAY);
  glVertexPointer(3,GL_FLOAT,0,coord);
  glNormalPointer(GL_FLOAT,0,coord);
  glColorPointer(3,GL_FLOAT,0,colors);
  glDrawArrays(GL_TRIANGLES,0,nElem*3);
  glDisableClientState(GL_NORMAL_ARRAY);    
  glDisableClientState(GL_COLOR_ARRAY);
  glDisableClientState(GL_VERTEX_ARRAY);      
    
#  string = "COUCOU"
#  glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, string[0]);
  coord *= 1.001;

  glColor3f(0.0, 0.0, 0.0);
  glEnableClientState(GL_VERTEX_ARRAY);    
  glVertexPointer(3,GL_FLOAT,0,coord);
  for i in range(nElem):
    glDrawArrays(GL_LINE_LOOP,3*i,3)
  glDisableClientState(GL_VERTEX_ARRAY) 
  glutSwapBuffers() 
 
# -------------------------------------------------------------------------
    
def reshape(width, height):   
  global theRatio
  
  glViewport(0,0,width,height)
  theRatio = width / float(height)
  glutPostRedisplay()

# -------------------------------------------------------------------------
 
def keyboard(key,x,y):
  global theFlagBathymetry, zoom, degreesHorizontal, degreesVertical, view, tab, run, paused, sleeper
  
  key = key.decode()
  if ord(key) == 27: # Escape
    sys.exit(0)
  elif key == 'v':
    view = tab[view - 1]
  elif key == 'd':
    degreesHorizontal += 4.0
  elif key == 'b':
    theFlagBathymetry = not(theFlagBathymetry)
  elif key == 'q':
    degreesHorizontal -= 4.0
  elif key == 'z':
    degreesVertical += 4.0
  elif key == 's':
    degreesVertical -= 4.0
  elif key == 'e':
    zoom -= 3 
  elif key == 'a':
    zoom += 3     
  elif key == 'r':
    run(1) 
  elif key == 'p':
    sleeper += 10**(-4)
  elif key == 'm':
    sleeper -= 10**(-4)
    if(sleeper < 0):
        sleeper = 0
  elif ord(key) == 32:
    paused = not(paused) 
  else:
    return
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
  
# -------------------------------------------------------------------------

def printf(x, y, Text):
    font = GLUT_BITMAP_9_BY_15
    blending = False 
    if glIsEnabled(GL_BLEND) :
        blending = True

    #glEnable(GL_BLEND)
    glColor3f(1,1,1)
    glRasterPos2f(x,y)
    for ch in Text :
        glutBitmapCharacter( font , ctypes.c_int( ord(ch) ) )


    if not blending :
        glDisable(GL_BLEND) 

# -------------------------------------------------------------------------
def idle():
  global iter,delta,E,theResultFiles, paused, sleeper
  
  try :
    if (paused == False):
        time.sleep(sleeper)
        iter += delta   
        E = tsunami.readResult(theResultFiles,iter,nElem)
        glutPostRedisplay()
  except FileNotFoundError: 
    pass

# -------------------------------------------------------------------------
  
iter = 0; delta = 500;
R = 6371220;
BathMax = 9368;
theMeshFile = "/Users/romaingraux/Library/Mobile Documents/com~apple~CloudDocs/Professionel/EPL/Q4/MAP/Elements finis/Projet/Mesh/PacificTriangleFine.txt"
theResultFiles = "/Users/romaingraux/Library/Mobile Documents/com~apple~CloudDocs/Professionel/EPL/Q4/MAP/Elements finis/Projet/perso-results/compute-%06d.txt"
theFlagBathymetry = False
paused = False
translationHorizontal = 0
translationVertical = 0
zoom = 0.0
degreesHorizontal = -116.0
degreesVertical = 0.0
theRatio = 0.8
view = 0
sleeper = 0
tab = [0, 1, 2]
[nNode,X,Y,H,nElem,elem] = tsunami.readMesh(theMeshFile)
E = np.zeros([nElem,3])
glutInit(sys.argv)
glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH)
# Full screen
glutInitWindowSize(glutGet(GLUT_SCREEN_WIDTH), glutGet(GLUT_SCREEN_HEIGHT))
glutCreateWindow("1")
matSpecular   = [1.0,1.0,1.0,0.0]
matShininess  = [50.0]
lightPosition = [8.0,8.0,8.0,0.0]
lightRadiance = [1.0,1.0,1.0,1.0]
glMaterialfv(GL_FRONT,GL_SPECULAR, matSpecular)
glMaterialfv(GL_FRONT,GL_SHININESS,matShininess)
glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,GL_TRUE)
glLightfv(GL_LIGHT0,GL_POSITION,lightPosition)
glLightfv(GL_LIGHT0,GL_DIFFUSE, lightRadiance)
glLightfv(GL_LIGHT0,GL_SPECULAR,lightRadiance)
glEnable(GL_LIGHTING)
glEnable(GL_LIGHT0)
glDepthFunc(GL_LEQUAL)
glEnable(GL_DEPTH_TEST)
glEnable(GL_COLOR_MATERIAL)
glEnable(GL_NORMALIZE)	

glutDisplayFunc(draw)
glutKeyboardFunc(keyboard)
glutSpecialFunc(special)
glutReshapeFunc(reshape)
glutIdleFunc(idle)
glClearColor( 0.9, 0.9, 0.8, 0.0 );
glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
glMatrixMode(GL_PROJECTION);
glLoadIdentity();
gluPerspective(65.0,1.0,1.0,100.0);
if "-info" in sys.argv:
  print("GL_RENDERER   = ",glGetString(GL_RENDERER).decode())
  print("GL_VERSION    = ",glGetString(GL_VERSION).decode())
  print("GL_VENDOR     = ",glGetString(GL_VENDOR).decode())

print('======================================')  
print(' b       : show bathymetry/elevation ')
print(' a | e   : zoom ')
print(' r       : restart ')
print(' p | m   : speed of the elevation ')
print(' q | d   : rotate the Earth ')
print(' UP/DOWN/LEFT/RIGHT : translate the Earth ')
print(' SPACE   : pause ')
print(' ESC     : exit ')
print('======================================')
# -------------------------------------------------------------------------

def run(boole):
    global iter, E, delta, R, BathMax, theMeshFile, theResultFiles, theFlagBathymetry, translationHorizontal, translationVertical, zoom , degreesHorizontal, degreesVertical, theRatio, nNode,X,Y,H,nElem,elem
    iter = 0; delta = 2;
    
     
    try :
      E = tsunami.readResult(theResultFiles,0,nElem)
    except FileNotFoundError:
      E = np.zeros([nElem,3])
    
    glutMainLoop()
    
run(0)