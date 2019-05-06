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
#from threading import Thread, Lock

_gaussTri3Xsi    = np.array([0.166666666666667,0.666666666666667,0.166666666666667])
_gaussTri3Eta    = np.array([0.166666666666667,0.166666666666667,0.666666666666667])
_gaussTri3Weight = np.array([0.166666666666667,0.166666666666667,0.166666666666667])

_gaussEdg2Xsi    = np.array([-0.5773502691896257, 0.5773502691896257])
_gaussEdg2Weight = np.array([1.0,1.0])

#mutex = Lock()

        
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
    for n in range(nIter):
        print(n)
        iterCompute(theTsunami,dt);
        if ((n % nSave) == 0):
            theTsunami.writeFile(theResultFiles, n)
            
    return [theTsunami.U,theTsunami.V,theTsunami.E]

# -------------------------------------------------------------------------
    
def iterCompute(Tsunami,dt):
    theSize = Tsunami.mesh.nElem
    Tsunami.iterE = np.zeros([theSize,3])
    Tsunami.iterU = np.zeros([theSize,3])
    Tsunami.iterV = np.zeros([theSize,3])
    
#    threadElem = Thread(target = computeElem,kwargs=dict(Tsunami = Tsunami));
#    threadEdges = Thread(target = computeEdge,kwargs=dict(Tsunami = Tsunami));
#    threadElem.start()
#    threadEdges.start()
#    threadElem.join()
#    threadEdges.join()
#    await mutex.acquire()
#    Tsunami.iterE = Tsunami.iterE1 + Tsunami.iterE2
#    Tsunami.iterU = Tsunami.iterU1 + Tsunami.iterU2
#    Tsunami.iterV = Tsunami.iterV1 + Tsunami.iterV2
    
    computeElem(Tsunami);
    computeEdge(Tsunami);
    inverseMatrix(Tsunami);
    Tsunami.E += dt * Tsunami.iterE 
    print(Tsunami.E[27])
    Tsunami.U += dt * Tsunami.iterU
    Tsunami.V += dt * Tsunami.iterV
    return 

# -------------------------------------------------------------------------
        
def computeElem(Tsunami):
    theMesh = Tsunami.mesh;
    theRule = Tsunami.rule2D;
    
    U      = Tsunami.U;
    V      = Tsunami.V;
    E      = Tsunami.E;
    R      = Tsunami.R; 
    gamma  = Tsunami.gamma
    g      = Tsunami.g
    xsi    = theRule.xsi;
    eta    = theRule.eta;
    weight = theRule.weight;
    phi    = np.array([1-xsi-eta,xsi,eta])
    
    Eh = E @ phi
    Uh = U @ phi
    Vh = V @ phi
    
    Xh = theMesh.xh
    Yh = theMesh.yh
    Zh = theMesh.zh
    
    Jac  = theMesh.jac
    Term = theMesh.term
    F    = theMesh.f
    
    Dphidx        = theMesh.dphidx
    Dphidy        = theMesh.dphidy
    
    EtermDphidx   = sum((Zh * Uh * np.outer(Jac,np.ones(3)) * np.outer(np.ones(theMesh.nElem),weight) * Term).T)
    EtermDphidy   = sum((Zh * Vh * np.outer(Jac,np.ones(3)) * np.outer(np.ones(theMesh.nElem),weight) * Term).T)
    UVtermDphi    = sum((g * Eh * Term * np.outer(Jac,np.ones(3)) * np.outer(np.ones(theMesh.nElem),weight)).T)
    
    EtermPhi      = Zh * (Xh * Uh + Yh *Vh) * np.outer(Jac,np.ones(3)) * np.outer(np.ones(theMesh.nElem),weight) / (R*R)
    UtermPhi      = (F * Vh - gamma * Uh + (g * Eh * Xh / (2*R*R))) * np.outer(Jac,np.ones(3)) * np.outer(np.ones(theMesh.nElem),weight)
    VtermPhi      = (-F * Uh - gamma * Vh + (g * Eh * Yh / (2*R*R))) * np.outer(Jac,np.ones(3)) * np.outer(np.ones(theMesh.nElem),weight)
    
    Tsunami.iterU = UtermPhi @ phi + Dphidx * np.outer(UVtermDphi,np.ones(3))
    Tsunami.iterV = VtermPhi @ phi + Dphidy * np.outer(UVtermDphi,np.ones(3))
    Tsunami.iterE = EtermPhi @ phi + Dphidx * np.outer(EtermDphidx,np.ones(3)) + Dphidy * np.outer(EtermDphidy,np.ones(3))
    
    return

        

# -------------------------------------------------------------------------
    
def computeEdge(Tsunami):
    theMesh = Tsunami.mesh
    theRule = Tsunami.rule1D;
    theEdges= Tsunami.edges;
    nBoundary = theEdges.nBoundary
    
    U       = Tsunami.U;
    V       = Tsunami.V;
    E       = Tsunami.E;
    iterU   = Tsunami.iterU;
    iterV   = Tsunami.iterV;
    iterE   = Tsunami.iterE;
    g       = Tsunami.g
    xsi    = theRule.xsi;
    weight = theRule.weight;
    phi = np.asarray([1.0-xsi,1.0+xsi])/ 2
    
    Jac = theEdges.jac
    Nx  = theEdges.nx
    Ny  = theEdges.ny
    
    Zh  = theEdges.zh
    
    Edges = theEdges.edges
    MapEdgeLeft  = theEdges.mapEdgeLeft
    MapEdgeRight = theEdges.mapEdgeRight
    ElemLeft     = Edges[:,2]
    ElemRight    = Edges[:,3]
#    print(MapEdgeRight)
    
    UhLeft  = U.flatten()[(3*np.outer(ElemLeft,np.ones(2))+MapEdgeLeft).flatten().astype(int)].reshape((theEdges.nEdges,2)) @ phi
    VhLeft  = V.flatten()[(3*np.outer(ElemLeft,np.ones(2))+MapEdgeLeft).flatten().astype(int)].reshape((theEdges.nEdges,2)) @ phi
    EhLeft  = E.flatten()[(3*np.outer(ElemLeft,np.ones(2))+MapEdgeLeft).flatten().astype(int)].reshape((theEdges.nEdges,2)) @ phi
    UnLeft  = UhLeft * np.outer(Nx,np.ones(2)) + VhLeft * np.outer(Ny,np.ones(2))
    
    UhRight = U.flatten()[(3*np.outer(ElemRight,np.ones(2))+MapEdgeRight).flatten().astype(int)].reshape((theEdges.nEdges,2)) @ phi
    VhRight = V.flatten()[(3*np.outer(ElemRight,np.ones(2))+MapEdgeRight).flatten().astype(int)].reshape((theEdges.nEdges,2)) @ phi
    EhRight = E.flatten()[(3*np.outer(ElemRight,np.ones(2))+MapEdgeRight).flatten().astype(int)].reshape((theEdges.nEdges,2)) @ phi
    UnRight = UhRight * np.outer(Nx,np.ones(2)) + VhRight * np.outer(Ny,np.ones(2))
    
    UnRight[0:nBoundary,:] = - UnLeft[0:nBoundary,:]
    EhRight[0:nBoundary,:] =   EhLeft[0:nBoundary,:]
    
    EStar   = 0.5 * ((EhLeft + EhRight) + np.sqrt(Zh/g) * (UnLeft - UnRight)) 
    UnStar  = 0.5 * ((UnLeft + UnRight) + np.sqrt(g/Zh) * (EhLeft - EhRight)) 
    Term    = theEdges.term
    Eterm   = ((np.outer(np.ones(theEdges.nEdges),weight) * Zh * UnStar * Term) @ phi) * np.outer(Jac,np.ones(2)) / 2
    Uterm   = ((np.outer(np.ones(theEdges.nEdges),weight) * EStar * Term) @ phi) * np.outer(Nx,np.ones(2)) * g * np.outer(Jac,np.ones(2)) / 2
    Vterm   = ((np.outer(np.ones(theEdges.nEdges),weight) * EStar * Term) @ phi) * np.outer(Ny,np.ones(2)) * g * np.outer(Jac,np.ones(2)) / 2
    
    ElLeft   = ((3 * np.outer(ElemLeft, np.ones(2)).flatten()) + MapEdgeLeft.flatten()).astype(int)
    ElRight  = ((3 * np.outer(ElemRight,np.ones(2)).flatten()) + MapEdgeRight.flatten()).astype(int)
    
    goodOnes = np.ones(2 * theEdges.nEdges)
    goodOnes[0:2 * nBoundary] -= 1
#    print(len(ElLeft.flatten()))
#    testE           =  iterE.flatten()
#    testE[ElLeft]  -= Eterm.flatten()
#    testE[ElRight] += Eterm.flatten() * goodOnes
#    testE = testE.reshape((theMesh.nElem,3))
#    testE.reshape((theMesh.nElem,3))
#    print(theMesh.nElem)
    iterE = iterE.ravel()
    iterU = iterU.ravel()
    iterV = iterV.ravel()
    
    Eterm = Eterm.flatten()
    Uterm = Uterm.flatten()
    Vterm = Vterm.flatten()
    
    VO = goodOnes * Vterm
    UO = goodOnes * Uterm
    EO = goodOnes * Eterm
    
#    iterE[ElLeft]   -= np.array(Eterm.flatten())
#    iterE[ElRight]  -= goodOnes * np.array(Eterm.flatten())
#    iterU[ElLeft]   -= np.array(Uterm.flatten())
#    iterU[ElRight]  += goodOnes * np.array(Uterm.flatten())
#    iterV[ElLeft]   += np.array(Vterm.flatten())
#    iterV[ElRight]  += goodOnes * np.array(Vterm.flatten())
#    iterE[ElLeft] -= Eterm
#    iterE -= np.bincount(ElLeft,  weights = Eterm, minlength = theMesh.nElem)
#    iterU -= np.bincount(ElLeft,  weights = Uterm, minlength = theMesh.nElem)
#    iterV -= np.bincount(ElLeft,  weights = Vterm, minlength = theMesh.nElem)
#    iterE += np.bincount(ElRight, weights = EO,    minlength = theMesh.nElem)
#    iterU += np.bincount(ElRight, weights = UO,    minlength = theMesh.nElem)
#    iterV += np.bincount(ElRight, weights = VO,    minlength = theMesh.nElem)
    
    for iEdge in range(theEdges.nEdges * 2):              
        iterE[ElLeft[iEdge]]   -= Eterm[iEdge]
        iterU[ElLeft[iEdge]]   -= Uterm[iEdge]
        iterV[ElLeft[iEdge]]   -= Vterm[iEdge]
        iterV[ElRight[iEdge]]  += VO[iEdge]
        iterU[ElRight[iEdge]]  += UO[iEdge]
        iterE[ElRight[iEdge]]  += EO[iEdge]
        
        
    iterE = iterE.reshape((theMesh.nElem,3))  
    iterU = iterU.reshape((theMesh.nElem,3))  
    iterV = iterV.reshape((theMesh.nElem,3))  
#    print(iterE[101])
#    print(testE[101])      
    return   

def inverseMatrix(Tsunami):
    theMesh  = Tsunami.mesh
    Ainverse = Tsunami.Ainverse
    JacAdit  = np.outer(theMesh.jacAdit,np.ones(3))
    
    Tsunami.iterE = (Tsunami.iterE @ Ainverse) / JacAdit
    Tsunami.iterU = (Tsunami.iterU @ Ainverse) / JacAdit
    Tsunami.iterV = (Tsunami.iterV @ Ainverse) / JacAdit
    return 

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
    

# -------------------------------------------------------------------------
    
class Mesh(object):
    
    def __init__(self,fileName,R):
        self.R = R
        dphidxsi = np.array([ -1.0, 1.0,0.0])
        dphideta = np.array([ -1.0, 0.0,1.0])
        xsi = _gaussTri3Xsi
        eta = _gaussTri3Eta
        phi = np.array([1-xsi-eta,xsi,eta])
        self.omega = 2*np.pi / 86400
        with open(fileName,"r") as f :
            self.nNode = int(f.readline().split()[3])
            self.xyz   = np.array(list(list(float(w) for w in f.readline().split()[2:]) for i in range(self.nNode)))
            self.nElem = int(f.readline().split()[3])
            self.elem  = np.array(list(list(int(w)   for w in f.readline().split()[2:]) for i in range(self.nElem)))
            self.x = self.xyz[:,0]
            self.y = self.xyz[:,1]
            self.z = self.xyz[:,2]
            self.xElem    = self.x[self.elem]
            self.yElem    = self.y[self.elem]
            self.zElem    = self.z[self.elem]
            self.dxdxsi   = self.x[self.elem] @ dphidxsi
            self.dxdeta   = self.x[self.elem] @ dphideta
            self.dydxsi   = self.y[self.elem] @ dphidxsi
            self.dydeta   = self.y[self.elem] @ dphideta 
            self.xStar    = np.array((4*R*R*self.x) / (4*R*R + self.x*self.x + self.y*self.y))
            self.yStar    = np.array((4*R*R*self.y) / (4*R*R + self.x*self.x + self.y*self.y))
            self.xh       = self.x[self.elem] @ phi
            self.yh       = self.y[self.elem] @ phi
            self.zh       = self.z[self.elem] @ phi
            self.jac      = abs(self.dxdxsi*self.dydeta - self.dydxsi*self.dxdeta)
            self.dphidx   = (np.outer(np.ones(self.nElem),dphidxsi) * np.outer(self.dydeta,np.ones(3)) - np.outer(np.ones(self.nElem),dphideta) * np.outer(self.dydxsi,np.ones(3))) / np.outer(self.jac,np.ones(3));
            self.dphidy   = (np.outer(np.ones(self.nElem),dphideta) * np.outer(self.dxdxsi,np.ones(3)) - np.outer(np.ones(self.nElem),dphidxsi) * np.outer(self.dxdeta,np.ones(3))) / np.outer(self.jac,np.ones(3));
            self.sinLat   = (4*R*R - self.xh*self.xh - self.yh*self.yh) / (4*R*R + self.xh*self.xh + self.yh*self.yh)
            self.term     = (4*R*R+self.xh*self.xh+self.yh*self.yh)/(4*R*R)        
            self.f        = 2 * self.omega * self.sinLat  
            self.jacAdit  = abs((self.xElem[:,0]-self.xElem[:,1]) * (self.yElem[:,0]-self.yElem[:,2]) - (self.xElem[:,0]-self.xElem[:,2]) * (self.yElem[:,0]-self.yElem[:,1])) 
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
    xsi = _gaussEdg2Xsi 
    phi = np.asarray([1.0-xsi,1.0+xsi])/ 2
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
        self.edges[id][3] = 0
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
    self.edges = np.array(self.edges)
    self.nodesLeft    = self.mesh.elem[self.edges[:,2]]
    self.nodesRight   = self.mesh.elem[self.edges[:,3]]
    self.mapEdgeLeft = np.zeros((self.nEdges,2),dtype = int); self.mapEdgeRight = np.zeros((self.nEdges,2), dtype = int) ;
    for i in range(self.nEdges):
        self.mapEdgeLeft[i,:] = [np.nonzero(self.nodesLeft[i] == self.edges[i,j])[0][0] for j in range(2)]
        if (i >= self.nBoundary):
            self.mapEdgeRight[i,:] = np.array([np.nonzero(self.nodesRight[i] == self.edges[i,j])[0][0] for j in range(2)])
    self.dx = self.mesh.x[self.edges[:,1]] - self.mesh.x[self.edges[:,0]]
    self.dy = self.mesh.y[self.edges[:,1]] - self.mesh.y[self.edges[:,0]]
    self.jac = np.sqrt(self.dx*self.dx + self.dy*self.dy)
    self.nx  =   self.dy / self.jac
    self.ny  = - self.dx / self.jac
    self.xh  = self.mesh.x[self.edges[:,0:2]] @ phi
    self.yh  = self.mesh.y[self.edges[:,0:2]] @ phi
    self.zh  = self.mesh.z[self.edges[:,0:2]] @ phi
    self.term= (4*self.mesh.R*self.mesh.R + self.xh*self.xh + self.yh*self.yh)/(4*self.mesh.R*self.mesh.R)

    def printf(self):
        print("Number of edges %d" % self.nEdges)
        print("Number of boundary edges %d" % self.nBoundary)
        for i in range(self.nEdges):
            print("%6d : %4d %4d : %4d %4d" % (i,*self.edges[i]))

# -------------------------------------------------------------------------
        
class Tsunami(object):
    
    def __init__(self,fileName,U,V,E):
        self.Ainverse = np.array([[18.0,-6.0,-6.0],[-6.0,18.0,-6.0],[-6.0,-6.0,18.0]])
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
        self.mapEdgeLeft = np.zeros((self.edges.nEdges,3),dtype=np.int)
        self.mapEdgeRight = np.zeros((self.edges.nEdges,3),dtype=np.int)
    
    def writeFile(self,fileName,iter):
        fileName = fileName % iter
        nElem = self.mesh.nElem;
        with open(fileName,"w") as f :
            f.write("Number of elements %d\n" % nElem)
            for i in range(nElem):
                f.write("%6d : %14.7e %14.7e %14.7e\n" % (i,*(self.E[i,:])))
        print(" === iteration %6d : writing %s ===" % (iter,fileName))
        

# -------------------------------------------------------------------------
