#!/usr/bin/env python
import sys,os
import numpy as np
import pylab as py
import lhapdf
import pandas as pd

class PDFCALC:

  def __init__(self,name,central_only=False,ismc=False):
    print 'loading ',name
    self.name=name
    self.ismc=ismc
    self.central_only=central_only
    if central_only==True:
      self.central=lhapdf.mkPDF(name,0)
    else:
      self.SETS=lhapdf.mkPDFs(name)

  def _get_xpdf(self,Set,flav,x,Q2):
    if   flav=='g': xpdf=Set.xfxQ2(21,x,Q2)
    elif flav=='u': xpdf=Set.xfxQ2(2,x,Q2)
    elif flav=='d': xpdf=Set.xfxQ2(1,x,Q2)
    elif flav=='s': xpdf= Set.xfxQ2(3,x,Q2)
    elif flav=='ub': xpdf= Set.xfxQ2(-2,x,Q2)
    elif flav=='db': xpdf= Set.xfxQ2(-1,x,Q2)
    elif flav=='sb': xpdf= Set.xfxQ2(-3,x,Q2)
    elif flav=='c': xpdf= Set.xfxQ2(4,x,Q2)
    elif flav=='b': xpdf= Set.xfxQ2(5,x,Q2)
    elif flav=='db+ub': xpdf= Set.xfxQ2(-2,x,Q2)+Set.xfxQ2(-1,x,Q2)
    elif flav=='db-ub': xpdf= Set.xfxQ2(-1,x,Q2)-Set.xfxQ2(-2,x,Q2)
    elif flav=='S': xpdf= 2*(Set.xfxQ2(-1,x,Q2)+Set.xfxQ2(-2,x,Q2)+Set.xfxQ2(-3,x,Q2)+Set.xfxQ2(-4,x,Q2)+Set.xfxQ2(-5,x,Q2))
    elif flav=='uv': xpdf= Set.xfxQ2(2,x,Q2)-Set.xfxQ2(-2,x,Q2)
    elif flav=='dv': xpdf= Set.xfxQ2(1,x,Q2)-Set.xfxQ2(-1,x,Q2)
    if np.isnan(xpdf): xpdf=0
    return xpdf

  def _get_xpdf_central(self,flav,x,Q2):
    if   flav=='g': xpdf= self.central.xfxQ2(21,x,Q2)
    elif flav=='u': xpdf= self.central.xfxQ2(2,x,Q2)
    elif flav=='d': xpdf= self.central.xfxQ2(1,x,Q2)
    elif flav=='s': xpdf= self.central.xfxQ2(3,x,Q2)
    elif flav=='c': xpdf= self.central.xfxQ2(4,x,Q2)
    elif flav=='b': xpdf= self.central.xfxQ2(5,x,Q2)
    elif flav=='db+ub': xpdf= self.central.xfxQ2(-2,x,Q2)+self.central.xfxQ2(-1,x,Q2)
    elif flav=='db-ub': xpdf= self.central.xfxQ2(-1,x,Q2)-self.central.xfxQ2(-2,x,Q2)
    elif flav=='ub': xpdf= self.central.xfxQ2(-2,x,Q2)
    elif flav=='db': xpdf= self.central.xfxQ2(-1,x,Q2)
    if np.isnan(xpdf): xpdf=0
    return xpdf

  def _get_symmetric_errors(self,OBS):
    n=len(OBS)-1
    feven=np.array([OBS[2*i] for i in range(1,n/2)])
    fodd=np.array([OBS[2*i-1] for i in range(1,n/2)])
    df=np.zeros(feven[0].size)
    for i in range(n/2-1):
      df+=(fodd[i]-feven[i])**2
    return df**0.5/2

  def _get_asymmetric_errors(self,OBS):
    n=len(OBS)-1
    f0=np.array(OBS[0])
    feven=np.array([OBS[2*i] for i in range(1,n/2)])
    fodd=np.array([OBS[2*i-1] for i in range(1,n/2)])
    dfeven=feven-f0
    dfodd=fodd-f0
    zeros=np.zeros(f0.size)
    dfP=np.zeros(f0.size)
    dfM=np.zeros(f0.size)
    for i in range(n/2-1):
      dfP+=np.amax([dfodd[i],dfeven[i],zeros],0)**2
      dfM+=np.amax([-dfodd[i],-dfeven[i],zeros],0)**2
    return dfP**0.5,dfM**0.5

  def get_xpdf(self,flav,X,Q2):
    D={}
    D['X']=X
    D['Q2']=Q2
    if self.central_only:
      D['xf0']=np.array([self._get_xpdf_central(flav,x,Q2) for x in X])
      D['dxf']=np.zeros(X.size)
      D['dxf+']=np.zeros(X.size)
      D['dxf-']=np.zeros(X.size)
    else:
      PDFS=[[self._get_xpdf(Set,flav,x,Q2) for x in X] for Set in self.SETS]
      if self.ismc==False:
        D['xf0']=np.array(PDFS[0])
        D['dxf']=self._get_symmetric_errors(PDFS)
        D['dxf+'],D['dxf-']=self._get_asymmetric_errors(PDFS)
      else:
        D['xf0']=np.mean(PDFS,axis=0)
        D['dxf+']=np.var(PDFS,axis=0)**0.5
        D['dxf-']=np.var(PDFS,axis=0)**0.5
    return D

  def get_d_over_u(self,X,Q2):
    D={}
    d=lambda Set,x,Q2: self._get_xpdf(Set,'d',x,Q2) 
    u=lambda Set,x,Q2: self._get_xpdf(Set,'u',x,Q2) 
    OBS=[[d(Set,x,Q2)/u(Set,x,Q2) for x in X] for Set in self.SETS]
    if self.ismc==False:
      D['central']=np.array(OBS[0])
      D['sym err']=self._get_symmetric_errors(OBS)
      D['asym err +'],D['asym err -']=self._get_asymmetric_errors(OBS)
    else:
      D['central']=np.mean(OBS,axis=0)
      D['asym err +']=np.var(OBS,axis=0)**0.5
      D['asym err -']=np.var(OBS,axis=0)**0.5
    return D

  def get_u_minus_d(self,X,Q2):
    D={'X':X}
    d=lambda Set,x,Q2: self._get_xpdf(Set,'d',x,Q2) 
    u=lambda Set,x,Q2: self._get_xpdf(Set,'u',x,Q2) 
    OBS=[[u(Set,x,Q2)-d(Set,x,Q2) for x in X] for Set in self.SETS]
    if self.ismc==False:
      D['central']=np.array(OBS[0])
      D['sym err']=self._get_symmetric_errors(OBS)
      D['asym err +'],D['asym err -']=self._get_asymmetric_errors(OBS)
    else:
      D['central']=np.mean(OBS,axis=0)
      D['asym err +']=np.var(OBS,axis=0)**0.5
      D['asym err -']=np.var(OBS,axis=0)**0.5
    return D

  def get_db_minus_ub(self,X,Q2):
    D={'X':X}
    db=lambda Set,x,Q2: self._get_xpdf(Set,'db',x,Q2) 
    ub=lambda Set,x,Q2: self._get_xpdf(Set,'ub',x,Q2) 
    OBS=[[db(Set,x,Q2)-ub(Set,x,Q2) for x in X] for Set in self.SETS]
    if self.ismc==False:
      D['central']=np.array(OBS[0])
      D['sym err']=self._get_symmetric_errors(OBS)
      D['asym err +'],D['asym err -']=self._get_asymmetric_errors(OBS)
    else:
      D['central']=np.mean(OBS,axis=0)
      D['asym err +']=np.var(OBS,axis=0)**0.5
      D['asym err -']=np.var(OBS,axis=0)**0.5
    return D

  def get_pvdis(self,X,Q2):

    sw2=0.2315
    eU=2*2/3.*(0.5-2*2/3.*sw2)
    eD=2*(-1/3.)*(-0.5+2*1/3.*sw2)
    eU2=4/9.
    eD2=1/9.

    D={}
    D['X']=X
    d=lambda Set,x,Q2: self._get_xpdf(Set,'d',x,Q2)+self._get_xpdf(Set,'db',x,Q2) 
    u=lambda Set,x,Q2: self._get_xpdf(Set,'u',x,Q2)+self._get_xpdf(Set,'ub',x,Q2) 
    s=lambda Set,x,Q2: self._get_xpdf(Set,'s',x,Q2)+self._get_xpdf(Set,'sb',x,Q2)
    F2gz=lambda Set,x,Q2:eU*u(Set,x,Q2)+eD*d(Set,x,Q2)+eD*s(Set,x,Q2)
    F2g=lambda Set,x,Q2:eU2*u(Set,x,Q2)+eD2*d(Set,x,Q2)+eD2*s(Set,x,Q2)

    OBS=[[F2gz(Set,x,Q2)/F2g(Set,x,Q2) for x in X] for Set in self.SETS]
    #OBS=[[F2g(Set,x,Q2) for x in X] for Set in self.SETS]

    if self.ismc==False:
      D['xf0']=np.array(OBS[0])
      D['dxf']=self._get_symmetric_errors(OBS)
      D['dxf+'],D['dxf-']=self._get_asymmetric_errors(OBS)
    else:
      D['xf0']=np.mean(OBS,axis=0)
      D['dxf+']=np.var(OBS,axis=0)**0.5
      D['dxf-']=np.var(OBS,axis=0)**0.5
    return D

  def get_Apv(self,X,Q2):
    D={}
    D['X']=X
    d=lambda Set,x,Q2: self._get_xpdf(Set,'d',x,Q2) 
    u=lambda Set,x,Q2: self._get_xpdf(Set,'u',x,Q2) 
    ub=lambda Set,x,Q2: self._get_xpdf(Set,'ub',x,Q2) 
    db=lambda Set,x,Q2: self._get_xpdf(Set,'db',x,Q2) 
    s=lambda Set,x,Q2: self._get_xpdf(Set,'s',x,Q2)
    sin2tW = 0.2315  # CJ code
    equ = 2/3.
    eqd = -1/3.
    gVu = +1/2. - 4/3.*sin2tW  # Cloet
    gVd = -1/2. + 2/3.*sin2tW  # Cloet
    c1u = 2*equ*gVu
    c1d = 2*eqd*gVd
    equ2 = equ*equ
    eqd2 = eqd*eqd
    OBS=[]
    for Set in self.SETS:
      OBSx=[]
      for x in X:
        xu=u(Set,x,Q2)
        xub=ub(Set,x,Q2)
        xd=d(Set,x,Q2)
        xdb=db(Set,x,Q2)
        xs=s(Set,x,Q2)
        xsb=s(Set,x,Q2)
        F2g = equ2*(xu + xub) + eqd2*(xd + xdb + xs + xsb)
        F2gZ = c1u*(xu + xub) + c1d*(xd + xdb + xs + xsb)
        apv = F2g#F2gZ#/F2g
        OBSx.append(apv)
      OBS.append(OBSx)

    if self.ismc==False:
      D['xf0']=np.array(OBS[0])
      D['dxf']=self._get_symmetric_errors(OBS)
      D['dxf+'],D['dxf-']=self._get_asymmetric_errors(OBS)
    else:
      D['xf0']=np.mean(OBS,axis=0)
      D['dxf+']=np.var(OBS,axis=0)**0.5
      D['dxf-']=np.var(OBS,axis=0)**0.5
    return D

  def get_pvdis_flav(self,X,Q2,flav):

    sw2=0.2315
    eU=2*2/3.*(0.5-2*2/3.*sw2)
    eD=2*(-1/3.)*(-0.5+2*1/3.*sw2)
    eU2=4/9.
    eD2=1/9.

    D={}
    D['X']=X
    d=lambda Set,x,Q2: self._get_xpdf(Set,'d',x,Q2)+self._get_xpdf(Set,'db',x,Q2) 
    u=lambda Set,x,Q2: self._get_xpdf(Set,'u',x,Q2)+self._get_xpdf(Set,'ub',x,Q2) 
    s=lambda Set,x,Q2: self._get_xpdf(Set,'s',x,Q2)+self._get_xpdf(Set,'sb',x,Q2)
    if flav=='u': F2gz_f=lambda Set,x,Q2:eU*u(Set,x,Q2)
    if flav=='d': F2gz_f=lambda Set,x,Q2:eD*d(Set,x,Q2)
    if flav=='s': F2gz_f=lambda Set,x,Q2:eD*s(Set,x,Q2)
    F2gz=lambda Set,x,Q2:eU*u(Set,x,Q2)+eD*d(Set,x,Q2)+eD*s(Set,x,Q2)

    OBS=[[F2gz_f(Set,x,Q2)/F2gz(Set,x,Q2) for x in X] for Set in self.SETS]

    if self.ismc==False:
      D['xf0']=np.array(OBS[0])
      D['dxf']=self._get_symmetric_errors(OBS)
      D['dxf+'],D['dxf-']=self._get_asymmetric_errors(OBS)
    else:
      D['xf0']=np.mean(OBS,axis=0)
      D['dxf+']=np.var(OBS,axis=0)**0.5
      D['dxf-']=np.var(OBS,axis=0)**0.5
    return D

  def get_F2g_flav(self,X,Q2,flav):

    sw2=0.2315
    eU=2*2/3.*(0.5-2*2/3.*sw2)
    eD=2*(-1/3.)*(-0.5+2*1/3.*sw2)
    eU2=4/9.
    eD2=1/9.

    D={}
    D['X']=X
    d=lambda Set,x,Q2: self._get_xpdf(Set,'d',x,Q2)+self._get_xpdf(Set,'db',x,Q2) 
    u=lambda Set,x,Q2: self._get_xpdf(Set,'u',x,Q2)+self._get_xpdf(Set,'ub',x,Q2) 
    s=lambda Set,x,Q2: self._get_xpdf(Set,'s',x,Q2)+self._get_xpdf(Set,'sb',x,Q2)
    if flav=='u': F2g_f=lambda Set,x,Q2:eU2*u(Set,x,Q2)
    if flav=='d': F2g_f=lambda Set,x,Q2:eD2*d(Set,x,Q2)
    if flav=='s': F2g_f=lambda Set,x,Q2:eD2*s(Set,x,Q2)
    F2g=lambda Set,x,Q2:eU2*u(Set,x,Q2)+eD2*d(Set,x,Q2)+eD2*s(Set,x,Q2)

    OBS=[[F2g_f(Set,x,Q2)/F2g(Set,x,Q2) for x in X] for Set in self.SETS]

    if self.ismc==False:
      D['xf0']=np.array(OBS[0])
      D['dxf']=self._get_symmetric_errors(OBS)
      D['dxf+'],D['dxf-']=self._get_asymmetric_errors(OBS)
    else:
      D['xf0']=np.mean(OBS,axis=0)
      D['dxf+']=np.var(OBS,axis=0)**0.5
      D['dxf-']=np.var(OBS,axis=0)**0.5
    return D

  def get_uv_dv(self,X,Q2):

    D={}
    D['X']=X
    dv=lambda Set,x,Q2: self._get_xpdf(Set,'d',x,Q2)-self._get_xpdf(Set,'db',x,Q2) 
    uv=lambda Set,x,Q2: self._get_xpdf(Set,'u',x,Q2)-self._get_xpdf(Set,'ub',x,Q2) 
    OBS=[[uv(Set,x,Q2)-dv(Set,x,Q2) for x in X] for Set in self.SETS]

    if self.ismc==False:
      D['xf0']=np.array(OBS[0])
      D['dxf']=self._get_symmetric_errors(OBS)
      D['dxf+'],D['dxf-']=self._get_asymmetric_errors(OBS)
    else:
      D['xf0']=np.mean(OBS,axis=0)
      D['dxf+']=np.var(OBS,axis=0)**0.5
      D['dxf-']=np.var(OBS,axis=0)**0.5
    return D


if __name__=="__main__":

  #pdf=PDFCALC('CJ15nlo')

  #X1=10**np.linspace(-4,-1,100)
  #X2=np.linspace(0.11,0.99,100)
  #X=np.append(X1,X2)
  #Q2=4.0
  #u_minus_d=pd.DataFrame(pdf.get_u_minus_d(X,Q2))
  #db_minus_ub=pd.DataFrame(pdf.get_db_minus_ub(X,Q2))

  #to_excel(u_minus_d,'u_minus_d')
  #to_excel(db_minus_ub,'db_minus_ub')
  #print db_minus_ub[['X','central',]] 
  #import pylab as py
  #py.plot(db_minus_ub['X'],db_minus_ub['central'])
  #py.plot(u_minus_d['X'],u_minus_d['central'])
  #py.semilogx()
  #py.show()


  pdf=lhapdf.mkPDF('CJ15nlo',0)
  x=0.5
  Q2=10.0
  print pdf.xfxQ2(21,x,Q2)



