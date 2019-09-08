        program BinCenter_Dp
        implicit none

        character*4 ratio
        character*80 infile,outfile
        real*8 xmin,xmax,xcenter,Q2,xavg,delx
        real*8 tmpX,tmp_xs
        real*8 XS_C1,XS_C2,xs_bin1,xs_bin2
        real*8 E0,Ep,cossq,sinsq,tansq
        real*8 Z1,Z2,A1,A2
        real*8 BC1,BC2,BC_fac
        integer kin,totN
        integer ii,jj

        real*8 xs_born
        external xs_born

        ratio='H3He'
        infile='../inputfile/'//trim(ratio)//'_x.dat'
        write(6,*) infile
        OPEN(UNIT=7,FILE=infile)

        outfile='OUT/model121/'//trim(ratio)//'_BCfac.dat'
        open(unit=66,file=outfile)

        if(ratio .eq. 'Dp') then
           Z1=1.0
           A1=2.0
           Z2=1.0
           A2=1.0
        endif

        if(ratio .eq. 'HeD') then
           Z1=2.0
           A1=3.0
           Z2=1.0
           A2=2.0
        endif

        if(ratio .eq. 'H3D') then
           Z1=1.0
           A1=3.0
           Z2=1.0
           A2=2.0
        endif

        if(ratio .eq. 'H3He') then
           Z1=1.0
           A1=3.0
           Z2=2.0
           A2=3.0
        endif

        do 99 ii=1,100
           READ(7,'(F7.4,2x,F7.4,2x,F10.5,2x,F7.3,2x,F10.6,2x,I3.2)',END=100)
     >     xmin,xmax,Q2,xcenter,xavg,kin
         
           XS_C1 = xs_born(xcenter,Q2,Z1,A1)
           XS_C2 = xs_born(xcenter,Q2,Z2,A2)

           delx=0.00001
           totN=int((xmax-xmin)/delx)-1;
           xs_bin1=0.0
           xs_bin2=0.0

           do 88 jj=0,totN
              tmpX=xmin+delx*jj
              tmp_xs = xs_born(tmpX,Q2,Z1,A1)
              xs_bin1=xs_bin1+tmp_xs
              tmp_xs = xs_born(tmpX,Q2,Z2,A2)
              xs_bin2=xs_bin2+tmp_xs     
88         continue            

           BC1=XS_C1/xs_bin1*(totN+1.0) 
           BC2=XS_C2/xs_bin2*(totN+1.0)
           BC_fac=BC1/BC2 
           write(66,'(5f10.5)')xavg,Q2,BC_fac,BC1,BC2

99      continue

100     end

c---------------------------------------------------------------------------------------------

        real*8 function xs_born(x,Q2,Z,A)
        real*8 x,Q2,Z,A
        real*8 Mp,E0,Ep,cossq,sinsq,tansq
        real*8 Wsq,nu,W1D,W2D,F1D,F2D
        real*8 F1,F2,W1,W2,R,DR
        real*8 sigmott,emccor,emciso,F2_NP,F_IS
        real*8 CJf2p,CJf2n
        integer model,D2_model,EMC_model,NP_model
        logical goodfit
        real*8 EMC_KP
        external EMC_KP
        real*8 emc_func_slac
        external emc_func_slac
        real*8 CJsfn
        external CJsfn
        real*8 F2NP_NMC
        external F2NP_NMC

        Mp=0.93827231
        E0=10.589

        nu=Q2/(2.0*Mp*x)
        Ep=E0-nu;
        Wsq=Mp**2+Q2*(1.0/x-1.0)
        sinsq=Q2/(4.0*E0*Ep)
        cossq=1.0-sinsq
        tansq=sinsq/cossq

        model=121
        D2_MODEL=MODEL/100
        EMC_MODEL=(MODEL-D2_MODEL*100)/10
        NP_MODEL=MODEL-D2_MODEL*100-EMC_MODEL*10

        if(A .eq. 1.0) then
           if(D2_model .eq. 1)then
              call ineft(Q2,sqrt(Wsq),W1,W2,dble(1.0))
              F2 = nu*W2
              F1 = Mp*W1
           endif
           if(D2_model .eq. 2)then
              call f2glob(X,Q2,'P',12,F2)
              call R1998(X,Q2,R,DR,goodfit)
              F1 = F2*(1.0+Q2/nu**2)/(2.0*x*(1.0+R))
           endif
         endif

        if(A .gt. 1.5) then
           if(D2_model .eq. 1)then
              call ineft(Q2,sqrt(Wsq),W1D,W2D,dble(2.0))
              F2D = nu*W2D
              F1D = Mp*W1D
           endif
           if(D2_model .eq. 2)then
              call f2glob(X,Q2,'D',12,F2D)
              call R1998(X,Q2,R,DR,goodfit)
              F1D = F2D*(1.0+Q2/nu**2)/(2.0*x*(1.0+R))
           endif

           emccor=1.0
           if(EMC_MODEL .EQ. 1) THEN 
              emccor=EMC_KP(x,Z,A)
           ENDIF
         
           if((EMC_MODEL .eq. 2) .and. (A.gt.2.5)) then
              emciso=emc_func_slac(x,A)

              if(NP_MODEL .eq. 1) then
                    F2_NP=1-0.8*x
              endif
              if(NP_MODEL .eq. 2) then
c                   F2_n/F2_p from CJ15
                    call setCJ(600)
                    CJf2p=CJsfn(1,x,sqrt(Q2))
                    CJf2n=CJsfn(2,x,sqrt(Q2))
                    F2_NP=CJf2n/CJf2p
              endif
              if(NP_MODEL .eq. 3) then
c                   F2_n/F2_p NMC paper 1992
                    F2_NP=F2NP_NMC(x,Q2)
              endif
              F_IS=(1+F2_NP)/(Z+(A-Z)*F2_NP)
              emccor=emciso/F_IS
           endif

           F2=F2D*emccor
           F1=F1D*emccor
         endif

        sigmott=(19732.0/(2.0*137.0388*E0*sinsq))**2*cossq/1.d6
        xs_born = 1d3*sigmott*(F2/nu+2.0*F1/Mp*tansq)
        return
        end

c---------------------------------------------------------------------------------------------

        real*8 function emc_KP(x,Z,A)
        real*8 x,Z,A

        emc_KP = 1.0
        if((Z .EQ. 1.0) .and. (A .eq. 3.0)) then
          emc_KP = 1.07251-1.61648*x+12.3626*x**2-65.5932*x**3+213.311*x**4
     >           -423.943*x**5+499.994*x**6-321.304*x**7+87.0596*x**8
          emc_KP = emc_KP*A/2.0
        endif

        if((Z .EQ. 2.0) .and. (A .eq. 3.0)) then
          emc_KP = 1.02967-0.135929*x+3.92009*x**2-21.2861*x**3
     >           +64.7762*x**4-129.928*x**5+169.609*x**6
     >           -127.386*x**7+41.0723*x**8
          emc_KP = emc_KP*A/2.0
        endif

        return
        end
c----------------------------------------------------------------------------------------------
c-------------------------------------------------------------------------------------------
        real*8 function emc_func_slac(x,A)
        real*8 x,A,atemp
        real*8 alpha,C
!       Javier EMC fit for isoscalar nuclei, Phys. Rev. D 49 (4348) 1994

        atemp = A
!       if(A.eq.4.or.A.eq.3) then  ! emc effect is more like C for these
!       2...
!          atemp = 12
!       endif

        alpha = -0.070+2.189*x - 24.667*x**2 + 145.291*x**3
     >         -497.237*x**4 + 1013.129*x**5 - 1208.393*x**6
     >         +775.767*x**7 - 205.872*x**8

        C = exp( 0.017 + 0.018*log(x) + 0.005*log(x)**2)
        
        emc_func_slac = C*atemp**alpha
        return
        end

c------------------------------------------------------------------------------------------------
        real*8 function F2NP_NMC(x,Q2)
        real*8 x,Q2,AX,BX

        AX=0.979-1.692*x+2.797*x**2-4.313*x**3+3.075*x**4
        BX=-0.171*x+0.244*x**2
        F2NP_NMC=AX*((Q2/20.0)**BX)*(1+x**2/Q2)

        return
        end

