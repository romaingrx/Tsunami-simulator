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

_gaussEdg2Xsi    = np.array([0.5773502691896257, -0.5773502691896257])
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
    for n in range(nIter):
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
    computeElem(Tsunami);
#    computeEdge(Tsunami);
#    inverseMatrix(Tsunami);
#    print(Tsunami.U[0:10,:])
    Tsunami.E += dt * Tsunami.iterE 
    Tsunami.U += dt * Tsunami.iterU
    Tsunami.V += dt * Tsunami.iterV
#    print(Tsunami.U[0:5,:])
    return 

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
    X       = theMesh.xStar;
    Y       = theMesh.yStar;
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
        
        xloc     = X[mapCoord]
        yloc     = Y[mapCoord]
        zloc     = Z[mapCoord]
        
        eloc     = E[iElem,:]
        uloc     = U[iElem,:]
        vloc     = V[iElem,:]
        
        dxdxsiRomain   = xloc @ dphidxsi
        dxdetaRomain   = xloc @ dphideta
        dydxsiRomain   = yloc @ dphidxsi
        dydetaRomain   = yloc @ dphideta
        
        jacRomain = abs(dxdxsiRomain*dydetaRomain - dydxsiRomain*dxdetaRomain)
        dphidx = (dphidxsi * dydetaRomain - dphideta * dydxsiRomain) / jacRomain;
        dphidy = (dphideta * dxdxsiRomain - dphidxsi * dxdetaRomain) / jacRomain;
        
        eh = phi @ eloc
        uh = phi @ uloc
        vh = phi @ vloc
        
        xh = phi @ xloc
        yh = phi @ yloc
        zh = phi @ zloc
        
        sinLat = (4*R*R - xh*xh - yh*yh) / (4*R*R + xh*xh + yh*yh)
        term = (4*R*R+xh*xh+yh*yh)/(4*R*R)        
        f = 2 * omega * sinLat
        
        for k in range(theRule.n):
#            xsik = xsi[k] 
#            etak = eta[k]
            weightk = weight[k]
            phik = phi[k,:]
            
            dxdxsi   = 0
            dxdeta   = 0
            dydxsi   = 0
            dydeta   = 0
            for i in range(3):
                dxdxsi += xloc[i]*dphidxsi[i];
                dxdeta += xloc[i]*dphideta[i];
                dydxsi += yloc[i]*dphidxsi[i];
                dydeta += yloc[i]*dphideta[i];
            jac = abs(dxdxsi*dydeta - dydxsi*dxdeta)
            x=0; y=0; h=0;e=0;u=0;v=0
            for i in range(3):
                x += xloc[i]*phik[i]
                y += yloc[i]*phik[i]
                h += zloc[i]*phik[i]
                e += eloc[i]*phik[i]
                u += uloc[i]*phik[i]
                v += vloc[i]*phik[i]
#            if (iElem == 0):
#                print(dphidx,dphidxRomain)
            lat = ((4*R*R-x*x-y*y))/(4*R*R+x*x+y*y)
            coriolis = 2*omega*lat
            sphere = ((4*R*R+x*x+y*y)/(4*R*R))
#            if (iElem == 0):
#                print(dphi)
            
#            for i in range(3):
#                iterE[iElem,i] += ((dphidx[i]*zh[k]*uh[k]+dphidy[i]*zh[k]*vh[k])*term[k] + (phik[i]*(zh[k]*(xh[k]*uh[k]+yh[k]*vh[k])/(R*R))))            * jac*weightk;
#                iterU[iElem,i] += ((phik[i]*(f[k]*vh[k]-gamma*uh[k])+dphidx[i]*g*eh[k]*term[k]) + (phik[i]*g*x*e/(2*R*R)))    * jac*weightk;
#                iterV[iElem,i] += ((phik[i]*(-coriolis*u-gamma*v) + dphidy[i]*g*e*sphere) + (phik[i]*g*y*e/(2*R*R))) * jac*weightk;
        EtermDphidx = sum(zh*uh*jac*weight*term)
        EtermDphidy = sum(zh*vh*jac*weight*term)
        EtermPhi = zh * (xh * uh + yh *vh) * jac * weight / (R*R)
        UtermPhi = (f * vh - gamma * uh + (g * eh * xh / (2*R*R))) * jac * weight
        UVtermDphi = sum(g * eh * term * jac * weight)
        VtermPhi = (-f * uh - gamma * vh + (g * eh * yh / (2*R*R))) * jac *weight
#        if(iElem == 0):
#            print(phi[:,0]*(g * xh / (2*R*R)))
#            print(((f[0] * vh[0] - gamma * uh[0] + g * xh[0] / (2*R*R)) * jac * weight[0] * phi[0,:]) + ((f[1] * vh[1] - gamma * uh[1] + g * xh[1] / (2*R*R)) * jac * weight[1] * phi[1,:]) + ((f[2] * vh[2] - gamma * uh[2] + g * xh[2] / (2*R*R)) * jac * weight[2] * phi[2,:]))
        iterU[iElem,:] += (UtermPhi @ phi) + (dphidx * UVtermDphi)
        iterV[iElem,:] += (VtermPhi @ phi) + (dphidy * UVtermDphi)
        iterE[iElem,:] += (EtermPhi @ phi) + (dphidx * EtermDphidx + dphidy * EtermDphidy)
        
        
#        if (iElem < 3):
#            print((phi.T @ termPhi) + sum(np.outer(termDphi,dphidx)))
#        if(iElem == 2):
#            print(iterU[iElem,:])
    print(iterU[0:3,:])
    return

        

# -------------------------------------------------------------------------
    
def computeEdge(Tsunami):
    theMesh = Tsunami.mesh;
    theRule = Tsunami.rule1D;
    theEdges= Tsunami.edges;
    
    U       = Tsunami.U;
    V       = Tsunami.V;
    E       = Tsunami.E;
    iterU   = Tsunami.iterU;
    iterV   = Tsunami.iterV;
    iterE   = Tsunami.iterE;
    X       = theMesh.xStar;
    Y       = theMesh.yStar;
    Z       = theMesh.zStar;
    R       = Tsunami.R; 
    g       = Tsunami.g
    xsi    = theRule.xsi;
    weight = theRule.weight;
    phi = np.asarray([1.0-xsi,1.0+xsi])/ 2
    for iEdge in range(theEdges.nBoundary,theEdges.nEdges):
        nodes = theEdges.edges[iEdge][0:2]
        
        x = X[nodes]
        y = Y[nodes]
        z = Z[nodes]
        
        mapEdgeLeft  = Tsunami.mapEdgeLeft[iEdge][1:3]
        mapEdgeRight = Tsunami.mapEdgeRight[iEdge][1:3]
        
        iElemLeft = Tsunami.mapEdgeLeft[iEdge][0]
        iElemRight = Tsunami.mapEdgeRight[iEdge][1]
        
        dx      = x[1] - x[0]
        dy      = y[1] - y[0]
        jac     = np.sqrt(dx*dx+dy*dy)
        nx      = dy / jac
        ny      = -dx / jac
        jac     = 0.5 * jac
        
        x = phi @ x
        y = phi @ y
        z = phi @ z
        
        uhLeft = phi @ U[iElemLeft][mapEdgeLeft]
        vhLeft = phi @ V[iElemLeft][mapEdgeLeft]
        uhRight = phi @ U[iElemRight][mapEdgeRight]
        vhRight = phi @ V[iElemRight][mapEdgeRight]
        eLeft   = phi @ E[iElemLeft][mapEdgeLeft]
        eRight = phi @ E[iElemRight][mapEdgeRight]  
        
        unLeft  = uhLeft * nx + vhLeft * ny
        unRight = uhRight * nx + vhRight * ny    
             
        eStar   = 0.5 * ((eLeft + eRight)   + np.sqrt(z/g) * (unLeft - unRight))     
        unStar  = 0.5 * ((unLeft + unRight) + np.sqrt(g/z) * (eLeft - eRight))
        
        term = (4*R*R+x*x+y*y)/(4*R*R)
      
        
        
                
        iterE[iElemLeft][mapEdgeLeft]   -= ((weight * z * unStar * term) @ phi) * jac
        iterE[iElemRight][mapEdgeRight] += ((weight * z * unStar * term) @ phi) * jac
        iterU[iElemLeft][mapEdgeLeft]   -= ((weight * eStar * term) @ phi) * nx * g * jac
        iterU[iElemRight][mapEdgeRight] += ((weight * eStar * term) @ phi) * nx * g * jac
        iterV[iElemLeft][mapEdgeLeft]   -= ((weight * eStar * term) @ phi) * ny * g * jac
        iterV[iElemRight][mapEdgeRight] += ((weight * eStar * term) @ phi) * ny * g * jac
        
        
    for iEdge in range(theEdges.nBoundary):
        nodes = theEdges.edges[iEdge][0:2]
        x = X[nodes]
        y = Y[nodes]
        z = Z[nodes]
        
        mapEdgeLeft  = Tsunami.mapEdgeLeft[iEdge][1:3]
            
        iElemLeft = Tsunami.mapEdgeLeft[iEdge][0]
            
        dx      = x[1] - x[0]
        dy      = y[1] - y[0]
        jac     = np.sqrt(dx*dx+dy*dy)
        nx      = dy / jac
        ny      = -dx / jac
        jac     = 0.5 * jac
        
        x = phi @ x
        y = phi @ y
        z = phi @ z
        
        uhLeft = phi @ U[iElemLeft][mapEdgeLeft]
        vhLeft = phi @ V[iElemLeft][mapEdgeLeft]
        
        unLeft  = uhLeft * nx + vhLeft * ny
        unRight = -unLeft
        eLeft   = phi @ E[iElemLeft][mapEdgeLeft]
        eRight  = eLeft
        
        eStar   = 0.5 * ((eLeft + eRight)   + np.sqrt(z/g) * (unLeft - unRight))     
        unStar  = 0.5 * ((unLeft + unRight) + np.sqrt(g/z) * (eLeft - eRight))
        
        term = (4*R*R+x*x+y*y)/(4*R*R)
        
        
        iterE[iElemLeft][mapEdgeLeft]   -= (phi @ (weight * z * unStar * term)) * jac 
        iterU[iElemLeft][mapEdgeLeft]   -= (phi @ (weight * eStar * term * nx)) * g * jac 
        iterV[iElemLeft][mapEdgeLeft]   -= (phi @ (weight * eStar * term * ny)) * g * jac 
            
    return 

def inverseMatrix(Tsunami):
    theMesh = Tsunami.mesh
    iterU   = Tsunami.iterU;
    iterV   = Tsunami.iterV;
    iterE   = Tsunami.iterE;
    Ainverse = np.array([[18.0,-6.0,-6.0],[-6.0,18.0,-6.0],[-6.0,-6.0,18.0]])
    for iElem in range(theMesh.nElem) :
      nodes = theMesh.elem[iElem]
      x = theMesh.x[nodes]
      y = theMesh.y[nodes]
      jac = abs((x[0]-x[1]) * (y[0]-y[2]) - (x[0]-x[2]) * (y[0]-y[1]))    
      iterE[iElem] = Ainverse @ iterE[iElem] / jac
      iterU[iElem] = Ainverse @ iterU[iElem] / jac
      iterV[iElem] = Ainverse @ iterV[iElem] / jac
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
        self.mapEdgeLeft = np.zeros((self.edges.nEdges,3),dtype=np.int)
        self.mapEdgeRight = np.zeros((self.edges.nEdges,3),dtype=np.int)
        for iEdge in range(self.edges.nBoundary):
          myEdge = self.edges.edges[iEdge]
          elementLeft  = myEdge[2]
          nodesLeft    = self.mesh.elem[elementLeft]
          self.mapEdgeLeft[iEdge,0]  = elementLeft
          self.mapEdgeLeft[iEdge,1:3]  = [ np.nonzero(nodesLeft  == myEdge[j])[0][0] for j in range(2)]
        for iEdge in range(self.edges.nBoundary,self.edges.nEdges):
          myEdge = self.edges.edges[iEdge]
          elementLeft  = myEdge[2]
          elementRight = myEdge[3]
          nodesLeft    = self.mesh.elem[elementLeft]
          nodesRight   = self.mesh.elem[elementRight]
          self.mapEdgeLeft[iEdge,0]  = elementLeft
          self.mapEdgeRight[iEdge,0] = elementRight
          self.mapEdgeLeft[iEdge,1:3]  = [ np.nonzero(nodesLeft  == myEdge[j])[0][0] for j in range(2)]
          self.mapEdgeRight[iEdge,1:3] = [ np.nonzero(nodesRight == myEdge[j])[0][0] for j in range(2)]         
    
    def writeFile(self,fileName,iter):
        fileName = fileName % iter
        nElem = self.mesh.nElem;
        with open(fileName,"w") as f :
            f.write("Number of elements %d\n" % nElem)
            for i in range(nElem):
                f.write("%6d : %14.7e %14.7e %14.7e\n" % (i,*(self.E[i,:])))
        print(" === iteration %6d : writing %s ===" % (iter,fileName))
        

# -------------------------------------------------------------------------
