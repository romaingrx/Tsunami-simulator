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

#import sys
#sys.path.insert(0, '/Users/romaingraux/Library/Mobile Documents/com~apple~CloudDocs/Professionel/EPL/Q4/MAP/Elements finis/Projet/src/Calculator/opt_einsum/optimize')
#from opt_einsum import contract
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
    theTsunami = Tsunami(theFile,U,V,E,dt);
    for n in range(nIter):
        print(n)
        iterCompute(theTsunami,dt);
        if ((n % nSave) == 0):
            theTsunami.writeFile(theResultFiles, n)
            
    return [theTsunami.U,theTsunami.V,theTsunami.E]

# -------------------------------------------------------------------------
    
def iterCompute(Tsunami,dt):
    theMesh       = Tsunami.mesh
    theEdges      = Tsunami.edges
    theSize       = Tsunami.mesh.nElem
    Tsunami.iterE = np.zeros([theSize,3])
    Tsunami.iterU = np.zeros([theSize,3])
    Tsunami.iterV = np.zeros([theSize,3])
    
    E        = Tsunami.E
    xElem    = theMesh.xElem
    yElem    = theMesh.yElem
    
    delta    = Tsunami.dt * 300
    
    Elements = np.unique(np.nonzero(E)[0])
#    print(Elements)
    
    XELEL    = xElem[Elements]
    YELEL    = yElem[Elements]
    xMax     = np.amax(XELEL) + delta
    xMin     = np.amin(XELEL) - delta
    yMax     = np.amax(YELEL) + delta
    yMin     = np.amin(YELEL) - delta
    
#    print(xStar([xMin,xMax],[yMin,yMax],Tsunami.R))
    
    xElemMax = np.einsum('ab->a',(xElem < xMax))
    xElemMin = np.einsum('ab->a',(xElem > xMin))
    yElemMax = np.einsum('ab->a',(yElem < yMax))
    yElemMin = np.einsum('ab->a',(yElem > yMin))
    
    xE = np.einsum('a,a->a',xElemMax,xElemMin)
    yE = np.einsum('a,a->a',yElemMax,yElemMin)
#    print(E)
    Tsunami.Elements = np.nonzero((xE + yE) > 0)[0]
    Tsunami.Edges    = np.nonzero(np.einsum('ab->a',np.isin(theEdges.edges[:,2:4],Elements)))[0]
#    print(len(Tsunami.Elements))
    
    
#    print(len(Tsunami.Elements))
    
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
    Tsunami.E += np.einsum(',ab->ab', dt, Tsunami.iterE)
#    print(sum(Tsunami.E.flatten()))
#    print(sum(Tsunami.E[Tsunami.Elements].flatten()))
    print(Tsunami.E[27])
    Tsunami.U += np.einsum(',ab->ab', dt, Tsunami.iterU)
    Tsunami.V += np.einsum(',ab->ab', dt, Tsunami.iterV)
#    print("Vraie valeur de E : %6f" %sum(Tsunami.E.flatten()))
#    print(sum(Tsunami.E[Tsunami.Elements].flatten()))
#    print("Valeur optimisée de E : %6f" % Tsunami.E[Tsunami.Elements])
    
    return 

# -------------------------------------------------------------------------
        
def computeElem(Tsunami):
    theMesh = Tsunami.mesh;
    theRule = Tsunami.rule2D;
    
    R      = Tsunami.R; 
    gamma  = Tsunami.gamma
    g      = Tsunami.g
    xsi    = theRule.xsi;
    eta    = theRule.eta;
    weight = theRule.weight;
    phi    = np.array([1-xsi-eta,xsi,eta])
    
    Elem = Tsunami.Elements
    
    Eh = np.einsum('ab,bc->ac',Tsunami.E[Elem], phi)
    Uh = np.einsum('ab,bc->ac',Tsunami.U[Elem], phi)
    Vh = np.einsum('ab,bc->ac',Tsunami.V[Elem], phi)
    
    Xh = theMesh.xh[Elem]
    Yh = theMesh.yh[Elem]
    Zh = theMesh.zh[Elem]
    
    Jac  = theMesh.jac[Elem]
    Term = theMesh.term[Elem]
    F    = theMesh.f[Elem]
    
    Dphidx        = theMesh.dphidx[Elem]
    Dphidy        = theMesh.dphidy[Elem] 
    
    outer         = np.einsum(      "a,b",       Jac, weight)
    outerTerm     = np.einsum("ab,ab->ab",     outer, Term)
    ZhouterTerm   = np.einsum("ab,ab->ab", outerTerm, Zh)
    
    EtermDphidx   = np.einsum( "ab,ab->a",        Uh, ZhouterTerm)
    EtermDphidx   = np.einsum( "ab,ab->a",        Uh, ZhouterTerm)
    EtermDphidy   = np.einsum( "ab,ab->a",        Vh, ZhouterTerm)
    UVtermDphi    = np.einsum(",ab,ab->a",   g,   Eh, outerTerm)
    
    EhgRR         = np.einsum(  ',ab->ab', g / (2*R*R),Eh)
    EtermPhi      = np.einsum("  ,ab->ab",   1 / (R*R), np.einsum("ab,ab->ab", np.einsum("ab,ab->ab", Zh, (np.einsum("ab,ab->ab", Xh, Uh) + np.einsum("ab,ab->ab", Yh, Vh))), outer))
    UtermPhi      = np.einsum('ab,ab->ab', (np.einsum('ab,ab->ab',  F, Vh) - np.einsum(',ab->ab', gamma, Uh) + (np.einsum('ab,ab->ab', EhgRR, Xh))), outer)
    VtermPhi      = np.einsum('ab,ab->ab', (np.einsum('ab,ab->ab', -F, Uh) - np.einsum(',ab->ab', gamma, Vh) + (np.einsum('ab,ab->ab', EhgRR, Yh))), outer)
#    print(len(VtermPhi))
    
    Tsunami.iterU[Elem] = np.einsum('ab,bc->ac', UtermPhi, phi) + np.einsum('ab,a->ab', Dphidx, UVtermDphi)
    Tsunami.iterV[Elem] = np.einsum('ab,bc->ac', VtermPhi, phi) + np.einsum('ab,a->ab', Dphidy, UVtermDphi)
    Tsunami.iterE[Elem] = np.einsum('ab,bc->ac', EtermPhi, phi) + np.einsum('ab,a->ab', Dphidx, EtermDphidx) + np.einsum('ab,a->ab', Dphidy, EtermDphidy)
    return
        

# -------------------------------------------------------------------------
    
def computeEdge(Tsunami):
    theMesh = Tsunami.mesh
    theRule = Tsunami.rule1D;
    theEdges= Tsunami.edges;
    
    g       = Tsunami.g
    xsi    = theRule.xsi;
    weight = theRule.weight;
    phi = np.asarray([1.0-xsi,1.0+xsi])/ 2

    EdgesIn = Tsunami.Edges 
#    print(Tsunami.Edges)
#    EdgesIn = np.arange(theEdges.nEdges)
    nEdges = len(EdgesIn)
    nBoundary = len(np.nonzero(np.clip(theEdges.edges[EdgesIn,-1],-theMesh.nElem - 10,0))[0])
    
    Jac = theEdges.jac[EdgesIn]
    Nx  = theEdges.nx[EdgesIn]
    Ny  = theEdges.ny[EdgesIn]
    
    Zh  = theEdges.zh[EdgesIn]
    
    Edges         = theEdges.edges[EdgesIn]
    MapEdgeLeft   = theEdges.mapEdgeLeft[EdgesIn]
    MapEdgeRight  = theEdges.mapEdgeRight[EdgesIn]
    ElemLeft      = Edges[:,2]
    ElemRight     = np.clip(Edges[:,3],0,theMesh.nElem+10)
    
    
    
    outerElLeft   = np.einsum('a,b',ElemLeft,np.ones(2))
    outerElRight  = np.einsum('a,b',ElemRight,np.ones(2))
    ElLeftMap     = (3*outerElLeft + MapEdgeLeft).flatten().astype(int)
    ElRightMap    = (3*outerElRight + MapEdgeRight).flatten().astype(int)    
    
    UEdgLeft = np.array([Tsunami.U[ElemLeft,MapEdgeLeft[:,0]],Tsunami.U[ElemLeft,MapEdgeLeft[:,1]]]).T
    VEdgLeft = np.array([Tsunami.V[ElemLeft,MapEdgeLeft[:,0]],Tsunami.V[ElemLeft,MapEdgeLeft[:,1]]]).T
    EEdgLeft = np.array([Tsunami.E[ElemLeft,MapEdgeLeft[:,0]],Tsunami.E[ElemLeft,MapEdgeLeft[:,1]]]).T
    
    UEdgRight = np.array([Tsunami.U[ElemRight,MapEdgeRight[:,0]],Tsunami.U[ElemRight,MapEdgeRight[:,1]]]).T
    VEdgRight = np.array([Tsunami.V[ElemRight,MapEdgeRight[:,0]],Tsunami.V[ElemRight,MapEdgeRight[:,1]]]).T
    EEdgRight = np.array([Tsunami.E[ElemRight,MapEdgeRight[:,0]],Tsunami.E[ElemRight,MapEdgeRight[:,1]]]).T
    
    UhLeft  = np.einsum('ab,bc->ac',UEdgLeft, phi)
    VhLeft  = np.einsum('ab,bc->ac',VEdgLeft, phi)
    EhLeft  = np.einsum('ab,bc->ac',EEdgLeft, phi)
    UnLeft  = np.einsum('ab,a->ab',UhLeft, Nx) + np.einsum('ab,a->ab',VhLeft, Ny)
    
    UhRight = np.einsum('ab,bc->ac',UEdgRight, phi)
    VhRight = np.einsum('ab,bc->ac',VEdgRight, phi)
    EhRight = np.einsum('ab,bc->ac',EEdgRight, phi)
    UnRight = np.einsum('ab,a->ab',UhRight, Nx) + np.einsum('ab,a->ab',VhRight, Ny)
    
    UnRight[0:nBoundary,:] = - UnLeft[0:nBoundary,:]
    EhRight[0:nBoundary,:] =   EhLeft[0:nBoundary,:]
    
    EStar        = ((EhLeft + EhRight) + np.einsum('ab,ab->ab', np.sqrt(Zh/g), (UnLeft - UnRight))) 
    UnStar       = ((UnLeft + UnRight) + np.einsum('ab,ab->ab', np.sqrt(g/Zh), (EhLeft - EhRight))) 
    Term         = theEdges.term[EdgesIn]
    
    weightTerm   = np.einsum('ab,b->ab',Term,weight)
    wTEStar      = np.einsum('ab,ab->ab', weightTerm, EStar)
    wTEStarphi   = np.einsum('ab,bc->ac', wTEStar, phi)
    wTESpJac     = np.einsum('ab,a->ab',wTEStarphi,Jac)
    OBO          = np.einsum(',ab->ab', 0.25 * g, wTESpJac)

    Eterm   = np.einsum('ab,a->ab',np.einsum('ab,bc->ac',np.einsum('ab,ab->ab',np.einsum('ab,ab->ab',weightTerm, UnStar), Zh), phi), Jac) / 4
    Uterm   = np.einsum('ab,a->ab',OBO, Nx)
    Vterm   = np.einsum('ab,a->ab',OBO, Ny)
    
    goodOnes2 = np.ones((nEdges,2))
    goodOnes2[:nBoundary] -= 1    
    
    goodOnes = np.ones(2 * nEdges)
    goodOnes[0:2 * nBoundary] -= 1
    
    Eterm2 = np.copy(Eterm)
    Uterm2 = np.copy(Uterm)
    Vterm2 = np.copy(Vterm)
    
    Eterm = Eterm.flatten()
    Uterm = Uterm.flatten()
    Vterm = Vterm.flatten()
    
    VO = np.einsum('a,a->a', goodOnes, Vterm)
    UO = np.einsum('a,a->a', goodOnes, Uterm)
    EO = np.einsum('a,a->a', goodOnes, Eterm)
    
#    VO2 = np.einsum('ab,ab->ab', goodOnes2, Vterm2)
#    UO2 = np.einsum('ab,ab->ab', goodOnes2, Uterm2)
#    EO2 = np.einsum('ab,ab->ab', goodOnes2, Eterm2)

    iterE = Tsunami.iterE
    iterU = Tsunami.iterU
    iterV = Tsunami.iterV
    
#    ETERM = Eterm2[np.nonzero(Eterm2[:,0]),0]
#    print(len(iterE[:,0]))
#    EbinLeft0 = np.bincount(ElemLeft,  weights = Eterm2[:,0])
#    ELeft0    = np.nonzero(EbinLeft0)[0]
#    print(ELeft0)
#    EbinLeft1 = np.bincount(ElemLeft,  weights = Eterm2[:,1])
#    ELeft1    = np.nonzero(EbinLeft1)[0]  
    
#    iterE[ELeft0,theEdges.mapEdgeLeft[ELeft0,0]] -= EbinLeft0[ELeft0]
#    iterE[ELeft1,theEdges.mapEdgeLeft[ELeft1,1]] -= EbinLeft1[ELeft1]
#    iterU[ElemLeft,MapEdgeLeft[:,0]] -= Uterm2[:,0]
#    iterU[ElemLeft,MapEdgeLeft[:,1]] -= Uterm2[:,1]
#    iterV[ElemLeft,MapEdgeLeft[:,0]] -= Vterm2[:,0]
#    iterV[ElemLeft,MapEdgeLeft[:,1]] -= Vterm2[:,1]
#    
#    iterE[ElemLeft,MapEdgeLeft[:,0]] += EO2[:,0]
#    iterE[ElemLeft,MapEdgeLeft[:,1]] += EO2[:,1]
#    iterU[ElemLeft,MapEdgeLeft[:,0]] += UO2[:,0]
#    iterU[ElemLeft,MapEdgeLeft[:,1]] += UO2[:,1]
#    iterV[ElemLeft,MapEdgeLeft[:,0]] += VO2[:,0]
#    iterV[ElemLeft,MapEdgeLeft[:,1]] += VO2[:,1]
    
    iterE = iterE.ravel()
    iterU = iterU.ravel()
    iterV = iterV.ravel()
    
    length = len(iterE)
#    print("EbinLeft0 : ",len(EbinLeft0))
#    print("Eterm[:,0] : ",len(Eterm2[:,0]))
    
#    print(len(Eterm2))
    iterE -= np.bincount(ElLeftMap,  weights = Eterm, minlength = length )
    iterU -= np.bincount(ElLeftMap,  weights = Uterm, minlength = length )
    iterV -= np.bincount(ElLeftMap,  weights = Vterm, minlength = length )
    iterE += np.bincount(ElRightMap, weights = EO,    minlength = length )
    iterU += np.bincount(ElRightMap, weights = UO,    minlength = length )
    iterV += np.bincount(ElRightMap, weights = VO,    minlength = length )
        
    iterE = iterE.reshape((theMesh.nElem,3))  
    iterU = iterU.reshape((theMesh.nElem,3))  
    iterV = iterV.reshape((theMesh.nElem,3))  
    
    return   

def inverseMatrix(Tsunami):
    theMesh  = Tsunami.mesh
    Ainverse = Tsunami.Ainverse
    JacInv = theMesh.jacInv
    
    Tsunami.iterE = np.einsum('a,ab->ab', JacInv, np.einsum('ab,bc->ac',Tsunami.iterE, Ainverse))
    Tsunami.iterU = np.einsum('a,ab->ab', JacInv, np.einsum('ab,bc->ac',Tsunami.iterU, Ainverse))
    Tsunami.iterV = np.einsum('a,ab->ab', JacInv, np.einsum('ab,bc->ac',Tsunami.iterV, Ainverse))
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
            self.dxdxsi   = np.einsum('ab,b->a',self.xElem, dphidxsi)
            self.dxdeta   = np.einsum('ab,b->a',self.xElem, dphideta)
            self.dydxsi   = np.einsum('ab,b->a',self.yElem, dphidxsi)
            self.dydeta   = np.einsum('ab,b->a',self.yElem, dphideta)
            self.xh       = np.einsum('ab,bc->ac',self.x[self.elem], phi)
            self.yh       = np.einsum('ab,bc->ac',self.y[self.elem], phi)
            self.zh       = np.einsum('ab,bc->ac',self.z[self.elem], phi)
            self.jac      = abs(self.dxdxsi*self.dydeta - self.dydxsi*self.dxdeta)
            self.jacInv   = 1 / self.jac
            self.dphidx   = np.einsum('a,ab->ab', self.jacInv , (np.einsum('a,b', self.dydeta, dphidxsi) - np.einsum('a,b', self.dydxsi, dphideta)));
            self.dphidy   = np.einsum('a,ab->ab', self.jacInv , (np.einsum('a,b', self.dxdxsi, dphideta) - np.einsum('a,b', self.dxdeta, dphidxsi)));
            self.sinLat   = (4*R*R - self.xh*self.xh - self.yh*self.yh) / (4*R*R + self.xh*self.xh + self.yh*self.yh)
            self.term     = (4*R*R+self.xh*self.xh+self.yh*self.yh)/(4*R*R)        
            self.f        = 2 * self.omega * self.sinLat  
            
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
    self.edges = np.zeros((self.nEdges,4)).tolist()
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
    
    def __init__(self,fileName,U,V,E,dt):
        self.Ainverse = np.array([[18.0,-6.0,-6.0],[-6.0,18.0,-6.0],[-6.0,-6.0,18.0]])
        self.R        = 6371220
        self.gamma    = 10e-7;
        self.g        = 9.81;
        self.omega    = 2 * np.pi / 86400;
        self.dt       = dt
        self.mesh     = Mesh(fileName,self.R);
        self.edges    = Edges(self.mesh);
        self.U        = U
        self.V        = V
        self.E        = E
        self.iterU    = 0
        self.iterV    = 0
        self.iterE    = 0
        self.rule1D   = IntegrationRule("Edge",2);
        self.rule2D   = IntegrationRule("Triangle",3);
        self.c        = np.sqrt(max(self.mesh.z) * self.g)
        Tsunami.Edges = 0
        Tsunami.Elements = 0
#        self.iterU1  = np.zeros([self.size,3])
#        self.iterV1  = np.zeros([self.size,3])
#        self.iterE1  = np.zeros([self.size,3])
#        self.iterU2  = np.zeros([self.size,3])
#        self.iterV2  = np.zeros([self.size,3])
#        self.iterE2  = np.zeros([self.size,3])
    
    def writeFile(self,fileName,iter):
        fileName = fileName % iter
        nElem = self.mesh.nElem;
        with open(fileName,"w") as f :
            f.write("Number of elements %d\n" % nElem)
            for i in range(nElem):
                f.write("%6d : %14.7e %14.7e %14.7e\n" % (i,*(self.E[i,:])))
        print(" === iteration %6d : writing %s ===" % (iter,fileName))
        

# -------------------------------------------------------------------------
