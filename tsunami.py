#
# PYTHON for FEM DUMMIES 18-19
# Projet "tsunami"
#
# Canevas de départ
#  Vincent Legat
#
# -------------------------------------------------------------------------
# 

import numpy as np

_gaussTri3Xsi    = np.array([0.166666666666667,0.666666666666667,0.166666666666667])
_gaussTri3Eta    = np.array([0.166666666666667,0.166666666666667,0.666666666666667])
_gaussTri3Weight = np.array([0.166666666666667,0.166666666666667,0.166666666666667])

_gaussEdg2Xsi    = np.array([-0.5773502691896257, 0.5773502691896257])
_gaussEdg2Weight = np.array([1.0,1.0])

        
def xStar(x,y,R):
    return 4*R*R*x / (4*R*R + x*x + y*y)
def yStar(x,y,R):
    return  4*R*R*y / (4*R*R + x*x + y*y)
def zStar(x,y,R):
    return R*(4*R*R - x*x - y*y) / (4*R*R + x*x + y*y)
def longitude(xStar,yStar):
    return np.arctan2(yStar,xStar)*180/np.pi
def latitude(zStar,R):
    return np.arcsin(zStar/R)*180/np.pi

# -------------------------------------------------------------------------
    
def readMesh(fileName) :
  with open(fileName,"r") as f :
    nNode = int(f.readline().split()[3])
    xyz   = np.array(list(list(float(w) for w in f.readline().split()[2:]) for i in range(nNode)))
    nElem = int(f.readline().split()[3])
    elem  = np.array(list(list(int(w)   for w in f.readline().split()[2:]) for i in range(nElem)))
  X = xyz[:,0]
  Y = xyz[:,1]
  H = xyz[:,2] 
  return [nNode,X,Y,H,nElem,elem]

# -------------------------------------------------------------------------
  
def readResult(fileBaseName,iter,nElem) :
  fileName = fileBaseName % iter
  with open(fileName,"r") as f :
    nSize = int(f.readline().split()[3])
    if (nElem != nSize) :
      print(" ==== Error : incoherent sizes : %d != %d" % (nElem,nSize))     
    E = np.array(list(list(float(w) for w in f.readline().split()[2:5]) for i in range(nElem)))
    print(" === iteration %6d : reading %s ===" % (iter,fileName))
  return E

# -------------------------------------------------------------------------

def writeResult(fileBaseName,iter,E) :
  fileName = fileBaseName % iter
  nElem = E.shape[0]  
  with open(fileName,"w") as f :
    f.write("Number of elements %d\n" % nElem)
    for i in range(nElem):
      f.write("%6d : %14.7e %14.7e %14.7e\n" % (i,*E[i,:]))
    print(" === iteration %6d : writing %s ===" % (iter,fileName))

# -------------------------------------------------------------------------    

def initialConditionOkada(x,y) :
        R = 6371220;
        x3d = 4*R*R*x / (4*R*R + x*x + y*y);
        y3d = 4*R*R*y / (4*R*R + x*x + y*y);
        z3d = R*(4*R*R - x*x - y*y) / (4*R*R + x*x + y*y);
        lat = np.arcsin(z3d/R)*180/np.pi;
        lon = np.arctan2(y3d,x3d)*180/np.pi;
        lonMin = 142;
        lonMax = 143.75;
        latMin = 35.9;
        latMax = 39.5;
        olon = (lonMin+lonMax)/2;
        olat = (latMin+latMax)/2;
        angle = -12.95*np.pi/180; 
        lon2 = olon + (lon-olon)*np.cos(angle) + (lat-olat)*np.sin(angle);
        lat2 = olat - (lon-olon)*np.sin(angle) + (lat-olat)*np.cos(angle);
        return np.all([lon2 <= lonMax,lon2 >= lonMin,lat2 >= latMin,lat2 <= latMax],axis=0).astype(int)
    
# -------------------------------------------------------------------------

def compute(theFile,theResultFiles,U,V,E,dt,nIter,nSave):
    theTsunami = Tsunami(theFile,U,V,E);
    theTsunami.initialize();
    for t in range(nIter):
        iterCompute(theTsunami);
        if ((t % dt) == 0):
            theTsunami.writeFile(theResultFiles % t)
            
    return [theTsunami.U,theTsunami.V,theTsunami.E]

# -------------------------------------------------------------------------
    
def iterCompute(Tsunami):
    return 1

# -------------------------------------------------------------------------
    
def computeElem(Tsunami):
    theMesh = Tsunami.mesh;
    theRule = Tsunami.rule2D;
    
    U       = Tsunami.U;
    V       = Tsunami.V;
    E       = Tsunami.E;
    iterU   = Tsunami.iterU;
    iterV   = Tsunami.iterV;
    iterE   = Tsunami.iterE;
    X       = theMesh.x;
    Y       = theMesh.y;
    Z       = theMesh.zStar;
    R       = Tsunami.R; 
    omega   = Tsunami.omega;
    gamma   = Tsunami.gamma
    g       = Tsunami.g
    xsi    = theRule.xsi;
    eta    = theRule.eta;
    weight = theRule.weight;
    phi = np.array([1-xsi-eta,xsi,eta])
    dphidxsi = np.array([ -1.0, 1.0,0.0])
    dphideta = np.array([ -1.0, 0.0,1.0])
    
    for iElem in range(theMesh.nElem):
        mapCoord = theMesh.elem[iElem];
        mapElem  = [3*iElem + j for j in range(3)]
        
        xloc     = theMesh.x[mapCoord]
        yloc     = theMesh.y[mapCoord]
        
        dxdxsi   = xloc @ dphidxsi
        dxdeta   = xloc @ dphideta
        dydxsi   = yloc @ dphidxsi
        dydeta   = yloc @ dphideta
        
        jac = abs(dxdxsi*dydeta - dydxsi*dxdeta)
        dphidx = (dphidxsi * dydeta - dphideta * dydxsi) / jac;
        dphidy = (dphideta * dxdxsi - dphidxsi * dxdeta) / jac;
        
        
        
        xh = phi @ X[mapCoord]
        yh = phi @ Y[mapCoord]
        zh = phi @ Z[mapCoord]
        eh = phi @ E[iElem,:]
        uh = phi @ U[iElem,:]
        vh = phi @ V[iElem,:]
        
        
        lat = (4*R*R-x*x-y*y)/(4*R*R+x*x+y*y)
        term = (4*R*R+x*x+y*y)/(4*R*R)
        f = 2*omega*np.sin(np.arcsin(((4*R*R - x*x - y*y)*180) / ((4*R*R + x*x + y*y)*np.pi) ))
        iterE[iElem,:] += sum((np.outer(zh*uh*term,dphidx)+np.outer(zh*1*term,dphidy))) + (phi @ ((zh*(xh*uh+yh*vh))/(4*R*R)))
        iterU[iElem,:] += phi @ 
        
#        iterE[iElem,:] += (jac * weight) * (3 * (z*u*dphidx + z*v*dphidy) + (z*(x*u+y*v)/(R*R)) * phi @ np.ones(3)) 
#        iterU[iElem,:] += (jac * weight) * (3 * (e*g*term*dphidx) + (f*1-gamma*u + g*x*e/(2*R*R)) * phi @ np.ones(3))
    print(iterU[:,:])
        
        
        
    return 1
        

# -------------------------------------------------------------------------
    
def computeEdge(Tsunami):
    return 1

def inverseMatrix(Tsunami):
    return 1
# -------------------------------------------------------------------------
# |//////////////////////////---- CLASS ----//////////////////////////////|   
# -------------------------------------------------------------------------
class IntegrationRule(object):

  def __init__(self,elementType,n):
    if (elementType == "Triangle" and n == 3) :
      self.name = "Gauss with 3 points"
      self.n = n
      self.xsi    = _gaussTri3Xsi
      self.eta    = _gaussTri3Eta
      self.weight = _gaussTri3Weight
    elif (elementType == "Edge" and n == 2) :
      self.name = "Gauss with 2 points"
      self.n = n
      self.xsi    = _gaussEdg2Xsi
      self.weight = _gaussEdg2Weight

  def printf(self):
    print(" Integration rule : %s " % self.name)
    print(" Number of nodes = %d " % self.n)
    print(" xsi     = ",self.xsi)
    print(" eta     = ",self.eta)
    print(" weights = ",self.weight)
    
  def func(self,x):
    return x

# -------------------------------------------------------------------------
    
class Mesh(object):
    
    def __init__(self,fileName,R):
        with open(fileName,"r") as f :
            self.nNode = int(f.readline().split()[3])
            self.xyz   = np.array(list(list(float(w) for w in f.readline().split()[2:]) for i in range(self.nNode)))
            self.xyz[:,2] = self.xyz[:,2].clip(100)
            self.nElem = int(f.readline().split()[3])
            self.elem  = np.array(list(list(int(w)   for w in f.readline().split()[2:]) for i in range(self.nElem)))
            self.xStar = self.xyz[:,0]
            self.yStar = self.xyz[:,1]
            self.zStar = self.xyz[:,2] 
            self.xy    = np.array(list(list(np.array([2*R*self.xStar[i] / (R + self.zStar[i]),2*R*self.yStar[i] / (R + self.zStar[i])]) for i in range(self.nNode))))
            self.x     = self.xy[:,0]
            self.y     = self.xy[:,1] 
            self.lola  = np.array(list(list(np.array([np.arctan2(self.yStar[i],self.xStar[i])*180/np.pi,np.arcsin(self.zStar[i]/R)*180/np.pi]) for i in range(self.nNode))))
            self.longitude = self.lola[:,0]
            self.latitude = self.lola[:,1]            
            
    def printf(self):
        print("Number of nodes %d" % self.nNode)
        print("")
        print("Mesh (3D)     xStar          yStar          zStar")
        for i in range(self.nNode):
            print("%6d : %14.7e %14.7e %14.7e" % (i,*self.xyz[i,:]))
        print("")
        print("Mesh (2D)       x              y")
        for i in range(self.nNode):
            print("%6d : %14.7e %14.7e" % (i,*self.xy[i,:]))
        print("")
        print("            Longitude      Latitude")
        for i in range(self.nNode):
            print("%6d : %14.7e %14.7e" % (i,*self.lola[i,:]))
        print("")
        print("Number of triangles %d" % self.nElem)
        for i in range(self.nElem):
            print("%6d : %6d %6d %6d" % (i,*self.elem[i,:]))
            
# -------------------------------------------------------------------------
            
class Edges(object):

  def __init__(self,mesh):
    self.mesh = mesh
    self.nEdges = mesh.nElem * 3
    self.nBoundary = 0
    self.edges = [[0 for i in range(4)] for i in range(self.nEdges)]
    for i in range (mesh.nElem) :
      for j in range(3) :
        id = i*3 + j
        self.edges[id][0] = mesh.elem[i][j]
        self.edges[id][1] = mesh.elem[i][(j+1)%3]
        self.edges[id][2] = i
        self.edges[id][3] = -1
    self.edges.sort(key = lambda item : -(min(item[0:2])*self.nEdges)-max(item[0:2]))
    index = 0
    for i in range(self.nEdges) :
      if (self.edges[i][0:2] != self.edges[i-1][1::-1]) :
         self.edges[index] = self.edges[i]
         index += 1
      else :
         self.edges[index-1][3] = self.edges[i][2]
    del self.edges[index:]
    self.edges.sort(key = lambda item : item[3])
    self.nBoundary = 2*index - self.nEdges
    self.nEdges = index

  def printf(self):
    print("Number of edges %d" % self.nEdges)
    print("Number of boundary edges %d" % self.nBoundary)
    for i in range(self.nEdges):
      print("%6d : %4d %4d : %4d %4d" % (i,*self.edges[i]))

# -------------------------------------------------------------------------
      
class System(object):
    
    def __init__(self,n):
        self.A      = np.zeros((n,n));
        self.B      = np.zeros(n);
        self.n      = n;
    
    def printff(self):
        for i in range(self.n):
            for j in range(self.n):
                if (self.A[i,j] == 0) :
                    print("     ", end='')
                else :
                    print(" %+.1f" % self.A[i,j],end='')
            print(" :  %+.1f" % self.B[i]);

# -------------------------------------------------------------------------
        
class Tsunami(object):
    
    def __init__(self,fileName,U,V,E):
        self.R = 6371220
        self.gamma = 10e-7;
        self.g = 9.81;
        self.omega = 2*np.pi / 86400;
        self.mesh   = Mesh(fileName,self.R);
        self.size   = self.mesh.nElem;
        self.edges  = Edges(self.mesh);
        self.U      = U
        self.V      = V
        self.E      = E
        self.iterU  = np.zeros([self.size,3])
        self.iterV  = np.zeros([self.size,3])
        self.iterE  = np.zeros([self.size,3])
        self.rule1D = IntegrationRule("Edge",2);
        self.rule2D = IntegrationRule("Triangle",3);

    
    
    def initialize(self):
        for iElem in range(self.mesh.nElem):
            mapCoord = self.mesh.elem[iElem]
            for j in range(3):
                self.E[3*iElem + j] = self.initialConditionOkada(self.mesh.x[mapCoord[j]],self.mesh.y[mapCoord[j]])
         
    
    def writeFile(self,fileName):
        nElem = self.E.shape[0]; 
        with open(fileName,"w") as f :
            f.write("Number of elements %d\n" % nElem)
            for i in range(nElem):
                f.write("%6d : %14.7e %14.7e %14.7e\n" % (i,*(self.E[i,:])))
        print(" === iteration %6d : writing %s ===" % (iter,fileName))
        

# -------------------------------------------------------------------------
