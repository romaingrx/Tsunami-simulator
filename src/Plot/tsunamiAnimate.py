#
# Projet "tsunami maker"
#
# Visualisation des résultats : OpenGL old-fashion
#
# Romain Graux & Lucas Delbecque
#
# -------------------------------------------------------------------------


from OpenGL.GL   import *
from OpenGL.GLU  import *
from OpenGL.GLUT import *
from PIL import Image

import sys
sys.path.insert(0, '/Users/romaingraux/Library/Mobile Documents/com~apple~CloudDocs/Professionel/EPL/Q4/MAP/Elements finis/Projet/src/Calculator/')

import tsunami as tsunami

import numpy as np
import time

# -------------------------------------------------------------------------

# Fonction qui dessine tout 
def draw():  
  global E,theFlagBathymetry,translationHorizontal,translationVertical,zoom,theRatio, degreesHorizontal, degreesVertical, view, iter, oneClick, textID
  
  lightness = 1.0
  glClearColor( *(3 * [lightness]), 0.0 );
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(65.0,theRatio,1.0,100.0);
  
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(0.0,1.0,0.0,
            0.0,20.0,0.0,
            0.0,0.0,1.0);  
            
  # === Affichage des vues si v est pressé sinon régit selon les paramètres === #            
  if(oneClick == False):
      glTranslatef(0.2*translationHorizontal,3.0,0.2*translationVertical);
      glTranslatef(0.0,0.25*zoom,0.0);
      glTranslatef(0.0,14.0,0.0);
      glRotatef(degreesHorizontal,0.0,0.0,-1.0);
  if(view == 0 and oneClick==True):
      degreesHorizontal = -130.0
      translationHorizontal = -20
      translationVertical = 0
      zoom = -15.0
      glTranslatef(0.2*translationHorizontal,3.0,0.2*translationVertical);
      glTranslatef(0.0,0.25*zoom,0.0);
      glTranslatef(0.0,14.0,0.0);
      glRotatef(degreesHorizontal,0.0,0.0,-1.0);  
      oneClick = False
     
  if(view == 1 and oneClick==True):
      degreesHorizontal = -140.0
      translationHorizontal = -10
      translationVertical = -15
      zoom = -34.0
      glTranslatef(0.2*translationHorizontal,3.0,0.2*translationVertical);
      glTranslatef(0.0,0.25*zoom,0.0);
      glTranslatef(0.0,14.0,0.0);
      glRotatef(degreesHorizontal,0.0,0.0,-1.0);
      oneClick = False
      
  if(view == 2 and oneClick==True):
      degreesHorizontal = -140.0
      translationHorizontal = -20
      translationVertical = -15
      zoom = -20.0
      glTranslatef(0.2*translationHorizontal,3.0,0.2*translationVertical);
      glTranslatef(0.0,0.25*zoom,0.0);
      glTranslatef(0.0,14.0,0.0);
      glRotatef(degreesHorizontal,0.0,0.0,-1.0);  
      oneClick = False
      
    
  # === Easter === #
  mult = 2
  verticies = (
        (-np.cos(degreesHorizontal* np.pi / 180)*mult*2.16,-np.sin(degreesHorizontal* np.pi / 180)*mult*2.16,-mult*1.44),
        (-np.cos(degreesHorizontal* np.pi / 180)*mult*2.16,-np.sin(degreesHorizontal* np.pi / 180)*mult*2.16, mult*1.44),
        ( np.cos(degreesHorizontal* np.pi / 180)*mult*2.16,np.sin(degreesHorizontal* np.pi / 180)*mult*2.16, mult*1.44),
        ( np.cos(degreesHorizontal* np.pi / 180)*mult*2.16,np.sin(degreesHorizontal* np.pi / 180)*mult*2.16,-mult*1.44))
  quadratic = gluNewQuadric();         
  gluQuadricNormals(quadratic, GLU_SMOOTH); 
  gluQuadricTexture(quadratic, GL_TRUE)
  glEnable(GL_TEXTURE_2D)
  glBindTexture(GL_TEXTURE_2D, textID)
  glBegin(GL_QUADS)
  glTexCoord2f(0,1)
  glVertex3f(*verticies[0])
  glTexCoord2f(0,0)
  glVertex3f(*verticies[1])
  glTexCoord2f(1,0)
  glVertex3f(*verticies[2])
  glTexCoord2f(1,1)
  glVertex3f(*verticies[3])
  glEnd()
  glDisable(GL_TEXTURE_2D)  
  
  # === Affichage de la sous-sphere un rien plus petite que les cooronnées pour éviter des défauts d'affichage === #
  glColor3f(0.8,1.0,0.7);
  gluSphere(quadratic,5.98,400,200);
  
  # === Affichage des coordonnées des triangles et prise en compte de la couleur avec la hauteur liée à E === #
  n = 9*nElem
  index  = np.ravel(elem)
  colors = np.zeros(n) 
  coord  = np.zeros(n)

  if (theFlagBathymetry == False) :
    value = 10* np.ravel(E)
  else :
    value = H[index]/BathMax
     
  value = np.asarray(value)
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
    
  coord *= 1.001;

  glColor3f(0.0, 0.0, 1.0);
  glEnableClientState(GL_VERTEX_ARRAY);    
  glVertexPointer(3,GL_FLOAT,0,coord);
  for i in range(nElem):
    glDrawArrays(GL_LINE_LOOP,3*i,3)
  glDisableClientState(GL_VERTEX_ARRAY) 
  
  # === Affichage des commandes dans la fenêtre === #
  
  # Prise en compte de la taille de la fenêtre pour le choix des différents paramètres
  Xmax = glutGet(GLUT_WINDOW_WIDTH); Ymax = glutGet(GLUT_WINDOW_HEIGHT);
  if(Xmax < 2000):
      xBegin = 72.5; yBegin = 92.5;
      font = 'medium'
      dx = 4.5; dy = 2.5
  else :
      xBegin = 80; yBegin = 92.5;
      font = 'big'
      dx = 3.0; dy = 2.5
  
  printf(xBegin,yBegin,"Listes des commandes",'black','big')
  
  # Listes les commandes en rouge sur la fenêtre
  printf(xBegin,yBegin-1*dy ,"esc"    ,'red',font)
  printf(xBegin,yBegin-2*dy ,'pause'  ,'red',font)
  printf(xBegin,yBegin-3*dy ,'r'      ,'red',font)
  printf(xBegin,yBegin-4*dy ,'b'      ,'red',font)
  printf(xBegin,yBegin-5*dy ,'p / m'  ,'red',font)
  printf(xBegin,yBegin-6*dy ,'a / e'  ,'red',font)
  printf(xBegin,yBegin-7*dy ,'q / d'  ,'red',font)
  printf(xBegin,yBegin-8*dy ,'flèches','red',font)
  printf(xBegin,yBegin-9*dy ,'v'      ,'red',font)
  
  # Listes les explications en noir sur la fenêtre
  printf(xBegin+dx,yBegin-1*dy ,': quitter le programme'        ,'black',font)
  printf(xBegin+dx,yBegin-2*dy ,': mettre pause à la simulation','black',font)
  printf(xBegin+dx,yBegin-3*dy ,': restart la simulation'       ,'black',font)
  printf(xBegin+dx,yBegin-4*dy ,': afficher la bathymétrie/simulation','black',font)
  printf(xBegin+dx,yBegin-5*dy ,': accélérer/ralentir la simulation','black',font)
  printf(xBegin+dx,yBegin-6*dy ,': zoom + / zoom -','black',font)
  printf(xBegin+dx,yBegin-7*dy ,': rotation horizontale autour de la terre','black',font)
  printf(xBegin+dx,yBegin-8*dy ,': translation','black',font)
  printf(xBegin+dx,yBegin-9*dy ,': 3 vues prédéfinie','black',font)
  
  
  # Affichage du temps sur la fenêtre
  printf(5,yBegin,'Temps : ','black','big')
  printf(5+1.5*dx,yBegin,str(convert(iter)),'yellow','big')
  
  glutSwapBuffers() 
 
# -------------------------------------------------------------------------

# Fonction qui gère le reshape de la fenêtre  
def reshape(width, height):   
  global theRatio
  
  glViewport(0,0,width,height)
  theRatio = width / float(height)
  glutPostRedisplay()

# -------------------------------------------------------------------------

# Fonction qui gère les entrées du clavier 
def keyboard(key,x,y):
  global theFlagBathymetry, zoom, degreesHorizontal, degreesVertical, view, tab, run, paused, speed, delta, oneClick
  
  key = key.decode()
  if ord(key) == 27: # Escape
    sys.exit(0)
  elif key == 'v':
    view = tab[view - 1]
    oneClick = True
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
    run() 
  elif key == 'p':
    speed += 1
    delta = mindelta*2*speed
  elif key == 'm':
    speed -= 1
    if(speed < 1):
        speed = 1
    delta = mindelta*2*speed
  elif ord(key) == 32:
    paused = not(paused) 
  else:
    return
  glutPostRedisplay()

# -------------------------------------------------------------------------

# Fonction qui gère les entrées des flèches
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

# Fonction qui print de manière fixe sur la fenêtre
def printf(xPourcent, yPourcent, Text, color, font):
    # color = black or white or red or yellow
    # font  = big or medium or small
    
    Xmax = glutGet(GLUT_WINDOW_WIDTH); Ymax = glutGet(GLUT_WINDOW_HEIGHT);
    x = int(xPourcent * Xmax / 100); y = int(yPourcent * Ymax /100)
    fonts = {
            'big'    : GLUT_BITMAP_TIMES_ROMAN_24,
            'medium' : GLUT_BITMAP_9_BY_15,
            'small'  : GLUT_BITMAP_HELVETICA_10,}
    colors = {
            'black'  : (0.0, 0.0, 0.0),
            'white'  : (1.0, 1.0, 1.0),
            'red'    : (1.0, 0.0, 0.0),
            'yellow' : (0.0, 0.0, 0.0),}
              
    blending = False 
    if glIsEnabled(GL_BLEND) :
        blending = True

    glColor3f(*colors[color])
    glWindowPos2i(x, y);
    for ch in Text :
        glutBitmapCharacter( fonts[font] , ctypes.c_int( ord(ch) ) )
        
    if not blending :
        glDisable(GL_BLEND) 

# -------------------------------------------------------------------------
        
# Fonction qui transforme un fichier photo en texture       
def read_texture(filename):
    img = Image.open(filename)
    img_data = np.array(list(img.getdata()), np.int8)
    textID = glGenTextures(1)
    glBindTexture(GL_TEXTURE_2D, textID)
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1)
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP)
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP)
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT)
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT)
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST)
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST)
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL)
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, img.size[0], img.size[1], 0, GL_RGB, GL_UNSIGNED_BYTE, img_data)
    return textID

# -------------------------------------------------------------------------

# Fonction qui transforme des secondes en heures, minutes, secondes    
def convert(seconds): 
    seconds = seconds % (24 * 3600) 
    hour = seconds // 3600
    seconds %= 3600
    minutes = seconds // 60
    seconds %= 60
      
    return "%d h %02d\'%02d\"" % (hour, minutes, seconds)

# -------------------------------------------------------------------------
    
# Fonction qui calcul les valeurs de E en tout temps
def idle():
  global iter,delta,E,theResultFiles, paused, speed  
  try :
    if (paused == False):
        iter += delta   
        E = tsunami.readResult(theResultFiles,iter,nElem)
        glutPostRedisplay()
  except FileNotFoundError: 
    pass

# -------------------------------------------------------------------------

#=== Initialisation de toutes les variables ===#
speed = 0
iter = 0; mindelta = 50
delta = mindelta;
R = 6371220;
BathMax = 9368;
theMeshFile = "/Users/romaingraux/Library/Mobile Documents/com~apple~CloudDocs/Professionel/EPL/Q4/MAP/Elements finis/Projet/Mesh/PacificTriangleFine.txt"
theResultFiles = "/Users/romaingraux/Library/Mobile Documents/com~apple~CloudDocs/Professionel/EPL/Q4/MAP/Elements finis/Projet/results/FineResults/compute-%06d.txt"
theFlagBathymetry = False # Affiche ou non la bathymétrie
paused = False # Met pause à la simulation
theRatio = 0.8 # Ratio de l'écran
oneClick = False # Se met à True si 'v' est pressé, évite de rester bloquer sur une vue dans la fonction 'draw'
view = 2 # L'initialise à la vue 2
tab = [0, 1, 2] # Toutes les valeurs de 'view'
degreesHorizontal = -130.0 # Initialise la rotation à -130°
translationHorizontal = -20 # Initialise la translation horizontale à -20
translationVertical = 0 # Initialise la translation verticale à 0
zoom = -15.0 # Initialise le zoom à -15
[nNode,X,Y,H,nElem,elem] = tsunami.readMesh(theMeshFile)
E = np.zeros([nElem,3])
Texture = "/Users/romaingraux/Library/Mobile Documents/com~apple~CloudDocs/Professionel/EPL/Q4/MAP/Elements finis/Projet/Textures/fem.jpg" # Image liée à la texture
textID = 0

#=== Initialise une fenêtre openGL en Full screen ===#
glutInit(sys.argv)
glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH)
glutInitWindowSize(glutGet(GLUT_SCREEN_WIDTH), glutGet(GLUT_SCREEN_HEIGHT))
glutCreateWindow("Tsunami Maker")

#=== Imposition de la lumière 
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

#=== Rends la terre 'smooth' entre la sphere et les vertex liés aux coordonnées des points ===#
glDepthFunc(GL_LEQUAL)
glEnable(GL_DEPTH_TEST)
glEnable(GL_COLOR_MATERIAL)
glEnable(GL_NORMALIZE)	

#=== Affectation des fonctions du programme dans les fonctions de glut ===#
glutDisplayFunc(draw)
glutKeyboardFunc(keyboard)
glutSpecialFunc(special)
glutReshapeFunc(reshape)
glutIdleFunc(idle)



if "-info" in sys.argv:
  print("GL_RENDERER   = ",glGetString(GL_RENDERER).decode())
  print("GL_VERSION    = ",glGetString(GL_VERSION).decode())
  print("GL_VENDOR     = ",glGetString(GL_VENDOR).decode())

# -------------------------------------------------------------------------

def run():
    global iter, E, R, BathMax, theMeshFile, theResultFiles, theFlagBathymetry, translationHorizontal, translationVertical, zoom , degreesHorizontal, degreesVertical, theRatio, nNode,X,Y,H,nElem,elem, paused, textID
    iter = 0; 
    paused = False
    textID = read_texture(Texture) #Charge la texture une seule fois 
     
    try :
      E = tsunami.readResult(theResultFiles,0,nElem)
    except FileNotFoundError:
      E = np.zeros([nElem,3])
    
    glutMainLoop() # Démarrage de l'animation
    

run() # Wouhou