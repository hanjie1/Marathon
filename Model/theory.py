#!/usr/bin/env python
import sys,os
import numpy as np
from scipy.integrate import quad,fixed_quad
import lhapdf

class DIS:
  
  def __init__(self):

    self.pdf=lhapdf.mkPDF("CJ15nlo", 0)
    self.mc=self.pdf.quarkThreshold(4)
    self.mb=self.pdf.quarkThreshold(5)
    self.TR=0.5
    self.CF=4./3.
    self.alfa=1/137.036
    self.M=0.93891897
    apU=4.0/9.0
    apD=1.0/9.0
    self.couplings={}
    self.couplings['p']={1:apD,2:apU,3:apD,4:apU,5:apD}
    self.couplings['n']={1:apU,2:apD,3:apD,4:apU,5:apD}
    self.fmap={}

    self.F2={'p':{},'n':{}}
    self.FL={'p':{},'n':{}}
 
  def integrator(self,f,xmin,xmax,method='gauss',n=100):
    f=np.vectorize(f)
    if method=='quad':
      return quad(f,xmin,xmax)[0]
    elif method=='gauss':
      return fixed_quad(f,xmin,xmax,n=n)[0]
    
  def log_plus(self,z,f,x):
    return np.log(1-z)/(1-z)*(f(x/z)/z-f(x)) + 0.5*np.log(1-x)**2*f(x)/(1-x)

  def one_plus(self,z,f,x):
    return 1/(1-z)*(f(x/z)/z-f(x))+ np.log(1-x)*f(x)/(1-x)

  def C2q(self,z,f,x):
    return self.CF*(2*self.log_plus(z,f,x)-1.5*self.one_plus(z,f,x)\
      +(-(1+z)*np.log(1-z)-(1+z*z)/(1-z)*np.log(z)+3+2*z)*f(x/z)/z\
      -(np.pi**2/3+4.5)*f(x)/(1-x))
    
  def C2g(self,z,f,x):
    return (((1-z)**2+z*z)*np.log((1-z)/z)-8*z*z+8*z-1)*f(x/z)/z
 
  def CLq(self,z,f,x):
    return 2*self.CF*z*f(x/z)/z #<--- note prefactor 2, instead of 4 used by MVV
    
  def CLg(self,z,f,x):
    return 4*z*(1-z)*f(x/z)/z

  def qplus(self,x,Q2):
    output=0
    for i in range(1,self.Nf+1):
      output+=self.couplings[self.tar][i]*(self.pdf.xfxQ2(i,x,Q2)/x+self.pdf.xfxQ2(-i,x,Q2)/x)
    return output

  def glue(self,x,Q2):
    output=0
    for i in range(1,self.Nf+1):
      output+=2*self.couplings[self.tar][i]
    return output*self.pdf.xfxQ2(21,x,Q2)/x
      
  def integrand_F2(self,x,z,Q2):
    return self.C2q(z,lambda y:self.qplus(y,Q2),x) + self.C2g(z,lambda y:self.glue(y,Q2),x)
    
  def integrand_FL(self,x,z,Q2):
    return self.CLq(z,lambda y:self.qplus(y,Q2),x) + self.CLg(z,lambda y:self.glue(y,Q2),x)
    
  def get_F2(self,x,Q2,tar):
    if (x,Q2) not in self.F2[tar]:
      self.tar=tar
      alphaS = self.pdf.alphasQ2(Q2)
      self.Nf=3
      if Q2>self.mc**2: self.Nf+=1
      if Q2>self.mb**2: self.Nf+=1
      LO=self.qplus(x,Q2)
      integrand=lambda z:self.integrand_F2(x,z,Q2)
      NLO=self.integrator(integrand,x,1)
      self.F2[tar][(x,Q2)]=x*(LO+alphaS/np.pi/2.0*NLO)
    return self.F2[tar][(x,Q2)]

  def get_FL(self,x,Q2,tar):
    if (x,Q2) not in self.FL[tar]:
      self.tar=tar
      alphaS = self.pdf.alphasQ2(Q2)
      self.Nf=3
      if Q2>self.mc**2: self.Nf+=1
      if Q2>self.mb**2: self.Nf+=1
      integrand=lambda z:self.integrand_FL(x,z,Q2)
      NLO=self.integrator(integrand,x,1)
      self.F2[tar][(x,Q2)]= x*alphaS/np.pi/2.0*NLO
    return self.FL[tar][(x,Q2)]

  def get_F1(self,x,Q2,tar):
    F2=get_F2(x,Q2,tar)
    FL=get_FL(x,Q2,tar)
    return ((1+4*self.M**2/Q2*x**2)*F2-FL)/(2*x)
 
  def get_dsigdxdQ2(self,x,y,Q2,target,precalc=False):
    if precalc==False: 
      return 4*np.pi*self.alfa**2/Q2**2/x*((1-y+y**2/2)*self.get_F2(x,Q2,target)-0*y**2/2*self.get_FL(x,Q2,target))
    else:
      return self.storage.retrieve([x,y,Q2,target])    


if __name__== "__main__":

  dis=DIS()
  for i in range(100):
    print i
    print dis.get_F2(0.5,100,'p')

  #M=0.9
  #E=12.0
  #S=M**2+2*M*E
  #x=0.5
  #Q2=10.0
  #y=Q2/x/S


  #x=0.13e-3
  #Q2=2.0
  #print dis.get_F2(x,Q2,'p')
  ##print dis.get_FL(x,Q2,'d')
  ##print dis.get_dsigdxdQ2(x,y,Q2,'d')

  #print dir(pdf)
























