      !File  snowdgtv22.f
      !Utah Energy Balance Snow Accumulation and Melt Model.
      !C
      !Authors
      !David G. Tarboton (dtarb@cc.usu.edu), Utah Water Research Laboratory, Utah State University
      !Charlie Luce cluce (cluce@fs.fed.us), USDA Forest Service, Rocky Mountain Research Station, Boise, ID
      !Graduate Students who have worked on this code
      !  Jinsheng You, PhD, 2004 (youjins@hotmail.com)
      !  Tanveer Chowdhury, MS, 1993
      !  Tom Jackson, PhD, 1994
      !c
      !Significant changes
      !5/4/04 Accommodate glacier surface runoff (version 2.1)
      !7/26/01 Jinsheng You dissertation work (version 2.0)
      !  - New surface temperature parameterizations comprising
      !     1) the Simple Gradient, almost the same as Original UEB,
      !     2) Simple Gradient approach with shallow snow correction.
      !     3) The classical force-restore approach.
      !     4) Modified force-restore approach.
      !  - New refreezing scheme
      !c
      !C
      !This program was written and prepared under agreement with and funding
      !by the U.S. Government , and therefore is in the public domain and not
      !subject to copyright.
      !C


! version modified by Chaopeng Shen (Michigan State Univ) to couple with the distributed hydrologic model
! in order to link, all REAL variables need to changed to REAL*8
! assembly input, sitev, statev, TsPrevday, taveprevday
! pay attention to the unit of dt, prcp, rad
! also uppack state variables so changes will be made onsite
MODULE SNOWUEB_MOD

CONTAINS
      SUBROUTINE SNOWUEB2(dt,nt,input,sitev,statev,                             &
     &   tsprevday, taveprevday, nstepday, param,iflag,                           &
     &   cump,cume,cummr,outv,mtime)
      USE COMM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      !    parameter(nsv=5,nxv=5,niv=7,nov=14)
      parameter(niv=7)

      real*8 MR
	  !real*8 lans, lang
      integer pflag,ounit,iflag(*)
      !data pi/3.141592654/ !define PI forangle conversion

      !Definitions
      !dt  Time step in hours
      !nt number of time steps
      !input  -- input forcing
      !  input(1,*) air temperature (C)
      !  input(2,*) precipitation (m/hr)
      !  input(3,*) wind speed  (m/s)
      !  input(4,*) relative humidity (on a 0-1 scale)
      !  input(5,*) incoming short wave  (kJ/m^2/h)
      !  input(6,*) net radiation  (kJ/m^2/h)
      !  input(7,*) Cosine of Zenith angle
      !SITEV -- site variables
      !      site variables (1-5)
      !      sitev(1)  forest cover fraction
      !      sitev(2)  drift factor  (No detailed information give 1)
      !      sitev(3)  air pressure (Pa)
      !      sitev(4) ground heat flux  Kj/m^2/hr (3.6 = 1 W/m^2)
      !      sitev(5) albedo extinction parameter (m)
      !STATEV
      !      statev(1)  Snow Energy Content  (KJ/m^2)
      !      statev(2)  Snow Water Equivalent (m) relative to T = 0 C solid phase
      !       statev(3)  Dimensionless age of snow surface (or albedo - depending on flag 4)
      !       statev(4)  Refreezing depth (m) used as refdepth
      !       statev(5)  Refreezing depth (m) used as totalrefdepth
      ! totalrefdepth is a misnomer.  These are the same quantity - the refreezing depth.  They are repeated because
      ! when refdepth exceeds the diurnal penetration depth it is set to 0 to tell the code to use the regular
      ! surface temperature functions while totalrefdepth is only reset to 0 when there is actual melt generated
      ! or energy content becomes negative indicating freezing of all liquid phase.  This ensures that the regular
      ! surface temperature functions persist when there is liquid water present but the refreezing front has penetrated
      ! to depth greater than diurnal temperature fluctuations.
      !       TsPrevday(1:nstepday)   Surface temperature over the last 24 hours
      !	   TavePrevday(1:nstepday)   Depth average temperature over the last 24 hours
      !c
      !PARAM  --  snowmelt model parameters (see below)
      !iflag  -- flags
      !   iflag(1) 0=radiation is shortwave in col 5 and longwave in col 6, else = net radiation in column 7
      !       iflag(2)        no 0 (/yes 1) printing
      !       iflag(3)  unit to which to print
      !       iflag(4)  how albedo calculations are done (a value 1 means albedo is calculated, otherwise statev(3) is albedo
      !	   iflag(5)  model option for surface temperature approximation
      !             1) the Simple Gradient, almost the same as Original UEB,
      !             2) Simple Gradient approach with shallow snow correction.
      !             3) The classical force-restore approach.
      !             4) Modified force-restore approach.
      !cump,cume,cummr  -- cumulative precipitation (with df), cumulative evaporation, cumulative melt over time step in m
      !outv  -- output variables
      !      outv(1)=prain   rain  m/hr
      !      outv(2)=ps     snow  m/hr
      !      outv(3)=a     albedo
      !      outv(4)=qh    sensible heat (kJ/m2/hr)
      !      outv(5)=qe    latent heat (kJ/m2/hr)
      !      outv(6)=e     sublimation m/hr
      !      outv(7)=mr    melt outflow m/hr
      !      outv(8)=qm    heat advected with melt
      !      outv(9)=q     Energy rate of change (kJ/m2/hr)
      !      outv(10)=fm   Mass rate of change (m/hr)
      !      outv(11)=tave  Average temperature (C)
      !      outv(12)=tsurf  Surface temperature C
      !      outv(13)=qnet  Net Radiation (kJ/m2/hr)
      ! outv(14)=smelt   Surface melt  m/hr
      !c
      !mtime   4 variable array with time information
      !   mtime(1)  year
      !   mtime(2)  month
      !   mtime(3)  day
      !   mtime(4)  hour (0-24 as a real number)
      !   mtime(5)  model element number
      !c
      ! Note: in this subroutine, the outv is an array which passes value to this subroutine and back to the snow
      ! drive program. The Outv(9) and outv(10) pass the daily average snowpack temperature and daily
      ! snow surface temperature to this subroutine but pass the Qin total and combined mass fluxes
      ! back.


      real*8 input(niv,nt)
      real*8 sitev(*)
      real*8 outv(*)
      real*8 statev(*)
      real*8 param(*)   ! inherit dimension from calling program
	  real*8 mtime(*)
	  real*8 tsprevday(*)
	  real*8 taveprevday(*)
	  integer iTsMethod       !yjs Add model time initialization 09/19/2000

	!common /ts_save/ ts_old, tave_old, Ts_Ave, Tave_ave

      !  Constant data set
      ! CP: declared in module comm

      !  End of constant declaration

      ! Parameters
       !c
      !  debugging
      !     if(mtime(1) .eq. 1948 .and. mtime(2) .eq. 11 .and.
      !    + mtime(3) .eq. 8 .and. mtime(4) .gt. 12 .and.
      !   + mtime(5) .eq. 169)then
      !   mtime(1)=mtime(1)
      !endif
      !    write(6,*)(input(j,1),j=1,7)
      ! State variables - These serve as initial conditions
      ub=statev(1)    ! Snow Energy Content  (KJ/m^2)
      w=statev(2)     ! Snow Water Equivalent (m) relative to T = 0 C solid phase
	  if(ub.le.0.0) then
	  refDepth = 0.0
	  totalRefDepth = 0.0
	  else
	  refDepth=statev(4)
	  totalRefDepth =statev(5)
      endif

      !	Save old Value 07/23/01
	ub_old = ub
	refDepth_old=refDepth
	qnetob = 0D0

      ! Site variables
      fc=sitev(1)     !  Forest cover fraction (0-1)
      df=sitev(2)     !  Drift factor
      pr=sitev(3)     !  Atmospheric Pressure (Pa)
      qg=sitev(4)     !  Ground heat flux (KJ/m^2/hr)  This is more logically an
      !             input variable,but is put here because it is never known at each
      !               time step. Usually it will be assigned a value 0.
      aep=sitev(5)   !  Albedo extinction parameter to smooth
      !    transition of albedo when snow is shallow. Depends on Veg. height (m)


      ! control flags
      iradfl=iflag(1)
      pflag=iflag(2)
      ounit=iflag(3)
	!iflag(4) albedo caculation
	  iTsMethod = iflag(5) ! the method to approximate the surface temperature
						! 1 normal old snow melt model
						! 2 revised direct gradient method (New method for ke) and qe
						! 3 force restore approach
						! 4 modified force restore approach

      ! model time step information
	  yy=mtime(1)
	  mm=mtime(2)
	  dd=mtime(3)
	  hr=mtime(4)

      ! Calculate constants that need only be calculated once.
       cd=k*k*hff/(log(z/zo)**2)*(1-0.8*FC)   ! factor in turbulent fluxes
      !  The FC correction is applied here to hit sensible and latent heat fluxes
      !  and evaporation rate in one place and maintain consistency in the surface
      !   temperature and stability calculations. It is consistent with wind speed
      !   reduction by a factor of 1-0.8 FC which is the most reasonable physical
      !   justification for this approach.
      !   I recognize that this is not a very good way to parameterize vegetation.
      !   This will be an area of future improvements (I hope).
      !  FC is also used below to adjust radiation inputs.
      RRHOI=RHOI/RHOW
      RRHO=RHO/RHOW
      RID=1.0/RRHO-1.0/RRHOI
	  rhom=lc*rho
      !c
      ! Loop for each time step
      !c
      do 2 i = 1,nt   ! DGT Time looping disabled as temperature averaging not handled here - must be handled outside
      ! Input variables
        ta=input(1,i)    ! Air temperature input (Degrees C)
        p=input(2,i)     ! Precipitation rate input (m/hr)
        ws=input(3,i)    ! Wind Speed (m/s)
        rh=input(4,i)    ! Relative humidity (fraction 0-1)
        if(iradfl.eq.0)then  ! input is incoming short and longwave
          qsi=input(5,i)   ! Incoming shortwave radiation (KJ/m^2/hr)
          qli=input(6,i)   ! Incoming longwave radiation (KJ/m^2/hr)
        else
          qnetob=input(5,i) ! Net allwave radiation (KJ/m^2/hr)
        endif
        coszen=input(7,i)   ! Cos(angle between direct sunlight and surface
      !                           normal).
      !       Representative value over time step used in albedo calculation.
      !       We neglect the difference between direct and diffuse albedo.
      ! Daily average temperatures handled internally so that multiple time steps will work
	  Ts_ave = daily_ave(tsprevday, nstepday, -100.D0)+tk ! (C)  Surface temperature average over last 24 hours
	  Tave_ave = daily_ave(taveprevday, nstepday, -100.D0)+tk ! (C)  Depth averaged temperature average over last 24 hours
        Ts_old = tsprevday(nstepday)+tk ! (C) Surface temperature from previous time step
        Tave_old = taveprevday(nstepday)+tk ! (C) Average temperature from previous time step

      !Separate rain and snow
        PS=PARTSNOW(P,TA,TR,TS)
        PRAIN = P - PS
      !Increase precipitation as snow by drift multiplication factor
        PS = PS*dF

      !      if(iflag(4).eq.1)then
      !Calculate albedo
      !        a=albedo(statev(3),coszen,w/rrho,aep,abg,avo,anir0)
      !Use of this albedo throughout time step neglects the
      !changes due to new snow within a time step.
      !      else
      !        a=statev(3)
      !      endif

      !  if there is any new snow without rain the snow age is change to 0.

      !     if(ps .gt. 0.0 .and. prain .le. 0.0) then
      !       statev(3) =0.0
      !  endif


        if(iflag(4).eq.1)then
      !Calculate albedo
          a=albedo(statev(3),coszen,w/rrho,aep,abg,avo,anir0)
      !Use of this albedo throughout time step neglects the
      !changes due to new snow within a time step.o
        else
          a=statev(3)
        endif
      cdh=0D0; cde=0D0

      ! Calculate neutral mass transfer coefficient
       rkn=cd*ws
      ! Adjust radiation inputs for the effect of forest canopy.
      qsi=qsi*(1.-fc)
      qli=qli*(1.-fc)
      qnetob=qnetob*(1.-fc)
      !   FC corrections are also in the following subroutines.
      !    QFM where outgoing longwave radiation is calculated.
      !    SNOTMP where outgoing longwave radiation is used in the surface
      !    temperature equilibrium approx.
      !Debugging
      !    if(mtime(1)== 1950 .and. mtime(2)==4 .and.
      !   + mtime(3)==1 .and. mtime(4) >12 .and. mtime(5)==169) then
      !   mtime(1)=mtime(1)
      !endif
      QE = 0D0; E = 0D0; Q = 0D0; TSurf = 0D0; Tave = 0D0
      ! Call predictor corrector subroutine to do all the work
        CALL PREDICORR(DT,UB,W,a,ta,prain,ps,ws,rh,qsi,                                     &
     &               qli,iradfl,rkn,                                                   &
     &               qnetob,cdh,cde,rid,param,sitev,iTsMethod,mtime,                          &
      !   Following variables are output
     &               QH,QE,E,                                                                   &
     &               MR,QM,Q,FM,TSURF,tave,qnet,refDepth, totalRefDepth,                          &
     &               smelt,gsurf)
      !c
      !DGT 4/2/05   Despite all checks in predicor It can (and does) occur that
      ! we still have ub so large that it results in tave greater than 0, which implies that all the
      ! snow is liquid.  In these cases - just force the snow to disappear and add the energy involved to Qm.
       tave=tavg(ub,w,rhow,cs,to,rhog,de,cg,hf)
      !if(tave < -1000)then
      !   tave=tave
      !endif
	 if(tave .gt. 0.)then   !  all is liquid so snow must disappear
		mr=mr+w/dt
		qm=qm+w/dt*rhow*hf
		q=q-w/dt*rhow*hf
		ub=ub-w/dt*rhow*hf
		w=0.
       endif
      ! 7/25/05   To guard against unreasonable UB when there is no snow do not allow bulk temperature to go above 10 C
       if(tave .gt. 10.)then
         ub=rhog*de*cg*10.
	 endif
      !cdgt 5/4/04 surface melt change
      !cyjsc   Update snow surface age based on snowfall in time step
      !       if(iflag(4).eq.1) call agesn(statev(3),dt,ps,tsurf,tk)
      ! Update snow surface age based on snowfall in time step
       if(iflag(4).eq.1) call agesn(statev(3),dt,ps,tsurf,tk,dNewS)
      !  accumulate for mass balance
       cump=cump+(ps+prain)*dt
       cume=cume+e*dt
       cummr=cummr+mr*dt
       tave=tavg(ub,w,rhow,cs,to,rhog,de,cg,hf)   !  this call
      ! necessary to keep track of average internal temperature used in some surface energy algorithms.
      !if(tave < -1000)then
      !   tave=tave
      !endif

      ! update the total depth of the refreezing depth in the snow pack according the
      ! the refreezing depth at time step and the positive energy input. 07/22/01
      !DGT's revised logic  1/13/05
      if(lc.gt.0.0) then
       if(refDepth.gt. 0.0) then
	    totalRefDepth=refDepth  ! if refreezing depth increases totalrefdepth keeps track
	 else  ! here refdepth has gone to 0 somehow
	  if(mr.gt.0.0 .or. (ub.gt.ub_old .and. ub .gt.0.0))                          &
      ! If there is melt or an increase in energy refdepth is reset
     &	    totalRefDepth = 0.0
	 endif
	elseif(mr.gt.0.0 .or. (ub.gt.ub_old .and. ub .gt.0.0)) then
      ! Here lc=0.  If there is melt or an increase in energy refdepth is reset
      ! This is likely redundant because if lc=0 then there is no meltwater to refreeze
	 totalRefDepth =0.0
	endif
      !Jinsheng's original logi      !
      !    if(lc.gt.0.0) then
      !     if((refDepth-refDepth_old) .gt. 0.0) then
      !    totalRefDepth=totalRefDepth +refDepth-refDepth_old
      ! else
      !  if(mr.gt.0.0)  totalRefDepth = 0.0 !totalRefDepth-mr*rhow/rhom
      !  if((ub-ub_old) .gt. 0.0 .and. ub. gt.0.0) then
      !    TotalRefDepth = totalRefDepth-(ub-ub_old)/(rhom*hf)
      !  endif
      ! endif
      !elseif(mr.gt.0.0 .or. (ub.gt.ub_old .and. ub .gt.0.0)) then
      ! totalRefDepth =0.0
      !endif
	if(totalRefDepth .lt. 0.0)  totalRefDepth=0.0
      !c
      ! update tsbackup and tavebackup
	 do 50 ii = 1 , nstepday-1
		tsprevday(ii)= tsprevday(ii+1)
	    taveprevday(ii)= taveprevday(ii+1)
 50	 continue
	   tsprevday(nstepday)= tsurf
	   taveprevday(nstepday)= tave


       IF(PFLAG.eq.1) THEN
            write(ounit,*) UB,W,statev(3),                                               &
     &                 prain,ps,a,qh,qe,e,mr,qm,q,fm,tave,tsurf                          &
     &                ,cump,cume,cummr,qnet,smelt,refdepth,totalrefdepth
      ! 5/4/04 surface melt smelt
       ENDIF
  2    continue

       statev(1)=ub
       statev(2)=w
	   statev(4)=refDepth
	   statev(5)=totalRefDepth
       outv(1)=prain
       outv(2)=ps
       outv(3)=a
       outv(4)=qh
       outv(5)=qe
       outv(6)=e
       outv(7)=mr
       outv(8)=qm
       outv(9)=q
       outv(10)=fm
       outv(11)=tave
       outv(12)=tsurf
       outv(13)=qnet
	   outv(14)=smelt   ! dgt 5/4/04
       RETURN
       END subroutine

      !*********************** PREDICORR () **************************
      !   Predictor-corrector scheme to update the state variables,
      !   U and W for each time step

      !       SUBROUTINE PREDICORR(DT,UB,W,a,ta,prain,ps,ws,rh,qsi,
      !     &                 qli,iradfl,rkn,
      !     &                 qnetob,cdh,cde,rid,param,sitev
      !   Following variables are output
      !     &                 QH,QE,E,
      !     &                 MR,QM,Q,FM,TSURF,qnet)

       SUBROUTINE PREDICORR(DT,UB,W,a,ta,prain,ps,ws,rh,qsi,                                &
     &               qli,iradfl,rkn,                                                         &
     &               qnetob,cdh,cde,rid,param,sitev,iTsMethod,mtime,                          &
      !   Following variables are output
     &               QH,QE,E,                                                                    &
     &               MR,QM,Q,FM,TSURF,tave,qnet,refDepth, totalRefDepth                          &
     &,smelt,gsurf)
      ! 5/4/04 surface melt smelt
      !  Constant data set
      USE comm
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      ! CP: declared in module comm

      !  End of constant declaration


      real*8 MR,mr1 !,Lans,LanG
      real*8 sitev(*)
      real*8 param(*)
	  real*8 mtime(*)  !yjs add model time

      data wtol,utol / 0.025D0,2000.D0/
      data ncall / 0/
      ncall=ncall+1

      CALL   QFM(ub,w,a,ta,prain,ps,ws,rh,qsi,qli,rkn,iradfl,qnetob,           &
     &           cdh,cde,rid,param,sitev,iTsMethod,mtime,                      &
     &           fm,q,qm,mr,qe,e,tsurf,tave,qh,qnet,dt,refDepth,               &
     &           totalRefDepth,smelt)     !yjs add average temperature
       ! 5/4/04 surface melt smelt
       !    PREDICTOR
       W1 = W + DT*FM
       IF(W1.LT.0.0) THEN
         w1=0.0
         CALL PREHELP(W1,W,DT,FM,0.D0,1.D0,PS,PRAIN,E,RHOW,HF,Q,QM,                 &
     &    MR,qe,hneu)
       ENDIF
       UB1 = UB + DT*Q

       Q1 = Q
       FM1 = FM
       ! save values so that they can be averaged for output
       qh1=qh
       qe1=qe
       e1=e
       mr1=mr
	   smelt1=smelt  !cdgt 5/4/04 surface melt smelt

       qm1=qm
       tsurf1=tsurf
       qnet1=qnet
       CALL QFM(ub1,w1,a,ta,prain,ps,ws,rh,qsi,qli,rkn,iradfl,qnetob,                     &
     &             cdh,cde,rid,param,sitev,iTsMethod,mtime,                               &
     &             fm,q,qm,mr,qe,e,tsurf,tave, qh,qnet,dt,refDepth,                       &
     &             totalRefDepth,smelt)  !yjs !yjs add average temperature
      ! 5/4/04 surface melt smelt

      !    CORRECTOR
       W2 = W + DT/2.0*(FM1 + FM)
       IF(W2.LT.0.0) THEN
         w2=0.0
         CALL PREHELP(W2,W,DT,FM,FM1,2.D0,PS,PRAIN,E,RHOW,HF,Q,QM,MR,qe,                 &
     &   hneu)
       ENDIF
       UB2 = UB + DT/2.0*(Q1 + Q)
      ! iterate to convergence to enhance stability
       niter=-3
       imax=1
  1    if((abs(w2-w1).gt.wtol .or. abs(ub2-ub1).gt.utol) .and.                          &
     &  niter .lt. imax)then
      !    write(6,10)ncall,"_ncall"
      !10   format(1x,i5,a6)
          w1=w2
          ub1=ub2
          CALL QFM(ub1,w1,a,ta,prain,ps,ws,rh,qsi,qli,rkn,iradfl,                        &
     &           qnetob,cdh,cde,rid,param,sitev,iTsMethod,mtime,                          &
     &           fm,q,qm,mr,qe,e,tsurf,tave,qh,qnet,dt,refDepth,                          &
     &           totalRefDepth,smelt)     !yjs add average temperature
      ! 5/4/04 surface melt smelt
      !  corrector again
         W2 = W + DT/2.0*(FM1 + FM)
         IF(W2.LT.0.0) THEN
           w2=0.0
           CALL PREHELP(W2,W,DT,FM,FM1,2.D0,PS,PRAIN,E,RHOW,HF,Q,QM,MR,                  &
     &    qe,hneu)
         ENDIF
         UB2 = UB + DT/2.0D0*(Q1 + Q)
         niter=niter+1
         if(niter .ge. 1)then  ! had * steps to converge now hit it.
      !What follows is a fix to numerical instability that results from
      !nonlinearity when the snowpack is shallow and melting rapidly.  If
      !convergence does not occur when the snowpack is not melting (a very
      !rare thing) I just accept the predictor corrector solution.
      !c
      !Modified by DGT 7/27/05 to add complete meltout condition
      !The quantities that this changes are w2, ub2, mr and qm
      !c
      !The fix first checks whether the added energy is sufficient to melt the snow
      !completely.  If this is the case then the snow disappears.
      !In other cases we assume that the liquid fraction of the snow remains constant.
      !This to some extent bypasses the melt outflow estimates.
      !ae is added energy during the time step.
          ae=(q1+q+qm1+qm)*.5*dt
      ! This fix is only physically sensible under melting conditions
      ! and when ae is positive and there is snow
          if(ub .gt. 0. .and. ae .gt. 0. .and. w .gt. 0.)then
	        e2=(e+e1)*.5   !  This is the average sublimation
      ! Check liquid fraction with added energy.  If this is greater than 1 then all snow melts
      ! Otherwise implement a solution assuming that the liquid fraction remains constant
	        rlf=(ub+ae)/(rhow*w*hf)
			if(rlf .ge. 1.)then
				mr=w/dt+(ps+prain-e2)   ! Here snow disappears
				if(mr .lt. 0.)then
					mr=0.   !  Force this to not go negative
      !    This can only occur if e2 is large compared to other terms.  Setting w2=0 implicitly reduces e2.
      !    There is a resulting energy discrepancy due to a limitation on sublimation and latent heat flux
      !    This is ignored because once the snow is gone energy computations are not pertinent.
                  endif
				qm=mr*rhow*hf
				w2=0.
				ub2=ub+ae-qm*dt
			else
      ! Determine the w/ub ratio at the beginning of the time step.
      ! Physically liquid fraction = ub/(rhow*w*hf) and since rhow and hf are constants
      ! keeping r=w/ub constant is equivalent to keeping liquid fraction constant.
      ! If liquid fraction is constant this is constant.
				r=w/ub
      ! Solve for ub2 and w2 that satisfy the three equations
      !          r=w2/ub2
      !          ub2=ub+ae-rhow*hf*mr*dt     Energy balance the last term being the energy advected by melt
      !          w2=w+(ps+prain-e2-mr)*dt    Mass balance
      ! The unknowns in the above are ub2, w2 and m and the equations are linear
      ! once the first eqn is multiplied by ub2
      ! The solution is
				ub2=(rhow*hf*(w+(ps+prain-e2)*dt)-ae-ub)/(rhow*hf*r-1)
				w2=r*ub2
				if(w2 .lt. 0.)then  ! Avoid negative W
					w2=0.
				endif
				mr=(w-w2)/dt-e2+ps+prain
				if(mr .lt. 0.)then   ! Avoid negative mr
					mr=0.
					w2=w+(ps+prain-e2)/dt
					if(w2 .lt. 0)then
						w2=0.
      !    This can only occur if e2 is large compared to other terms.  Setting w2=0 implicitly reduces e2.
      !    There is a resulting energy discrepancy due to a limitation on sublimation and latent heat flux
      !    This is ignored because once the snow is gone energy computations are not pertinent.
					endif
				endif
				qm=mr*rhow*hf
				ub2=ub+ae-qm*dt   ! redundant most of the time but recalc due to exceptions
			endif
      !  Check that nothing went wrong
              if(mr .lt. 0.D0)then
	            write(6,*)'Error - negative melt rate in snow'
			endif
	        if(w2 .lt. 0.D0)then
	            write(6,*)'Error - negative w2 in snow'
			endif
			q=ae/dt-qm
      ! Now set first pass values equal to final values to fake the averages below
			qm1=qm
			mr1=mr
			q1=q
          endif
         endif
         go to 1
       endif
       W = W2
       ub=ub2
      !average values from two time steps for output.  This is done for mr
      !and e to ensure mass balance and the others for better physical
      !comparisons
       qh=(qh+qh1)*.5D0
       qe=(qe+qe1)*.5D0
       e=(e+e1)*.5D0
       mr=(mr+mr1)*.5D0
       qm=(qm+qm1)*.5D0
       tsurf=(tsurf+tsurf1)*.5D0
       qnet=(qnet+qnet1)*.5D0
       q=(q+q1)*.5D0
	 smelt=(smelt+smelt1)*.5D0/(hf*rhow)  !cdgt 5/4/04 surface melt smelt
      !convert from energy KJ/m^2/hr to water depth m/hr of melt.
       RETURN
       END subroutine
      !************************ QFM () ********************************
      !   Calculates Energy and Mass Flux at any instant

      !       SUBROUTINE QFM(ub,w,a,ta,prain,ps,ws,rh,qsi,qli,rkn,iradfl,
      !     &                qnetob,cdh,cde,rid,param,sitev,
      !     &                fm,q,qm,mr,qe,e,tsurf,qh,qnet)
       SUBROUTINE QFM(ub,w,a,ta,prain,ps,ws,rh,qsi,qli,rkn,iradfl,                      &
     &           qnetob,cdh,cde,rid,param,sitev,iTsMethod,mtime,                         &
     &           fm,q,qm,mr,qe,e,tsurf,tave,qh,qnet,dt,refDepth,                          &
     &           totalRefDepth,smelt)     !yjs add average temperature
      ! 5/4/04 surface melt smelt
      USE COMM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      real*8 sitev(*)
      real*8 param(*)
	  real*8 mtime(*)   !yjs 09/19/2000 add the model time
      real*8 mr !,lc,ks,k,Lans,LanG

      !  Constant data set
      ! CP: declared in module comm


	  !common /ts_save/ ts_old, tave_old, Ts_ave, Tave_ave

      !  End of constant declaration

      ! gsurf added for modeling surface runoff from a glacier

      ! Site variables
      fc=sitev(1)     !  Forest cover fraction (0-1)
      !    df=sitev(2)    !  Drift factor
      pr=sitev(3)     !  Atmospheric Pressure (Pa)
      qg=sitev(4)     !  Ground heat flux (KJ/m^2/hr)  This is more logically an
      !    input variable, but is put here because it is never known at each
      !    time step.  Usually it will be assigned a value 0.

      QSN = QSI*(1.0D0-A)
      !   To ensure all temperatures in kelvin
      TAK = TA+TK
      TAVEK  = TAVE+TK
      DENSA  = PR/(RA*TAK)     ! Density of Air in kg/m3 if PR in Pa
      DENS = RHO
	  rhom=lc*rho

      !Cysj calculate the ks=k/ze
	  fKappaS = lans/(rho*cs)
	  ds=sqrt(2*fKappaS/w1day)
	  rs=Lans/(ds*rd1*rho*cs)


      ! save the old values
	tave_old=tave

       EA=SVPW(TA)*RH   ! The saturation vapour pressure over water is used
      !     because most instruments report relative humidity relative to water.
       QP=QPF(PRAIN,TA,TO,PS,RHOW,HF,CW,CS)
       tave=tavg(ub,w,rhow,cs,to,rhog,de,cg,hf)

      ! as described as below, the original UEB does not consider the refreezing. in this
      !change, the model will consider the refreezing effect starting from the top to the bottom.
      !Here is the predictor.

      ! the method is: refreezing can be modeled if the ub is greater than 0.0 without any precipiation,
      !meanwhile the total refreezing depth in the snowpack is less than the Refreezing depth times the daily damping of wet snow.

	  if(Ub .gt.0 .and. ps.le.0.0 .and. prain.le.0.0 .and.                               &
     &  totalRefDepth .le. rd1*ds .and. w.gt.0.0) then
	  qc1=QcEst(TK,RKN,WS,TAK,QP,DENSA,                                                   &
     &      PR,EA,DENS,RS,TAVEK,qsn,qli,fc,qnetob,iradfl,                                 &
     &      mtime,param,iTsMethod,w,dt)
	  qc2=QcEst(TK-0.01D0,RKN,WS,TAK,QP,DENSA,                                           &
     &      PR,EA,DENS,RS,TAVEK,qsn,qli,fc,qnetob,iradfl,                                &
     &      mtime,param,iTsMethod,w,dt)


	  call Grad(qc1,qc2,0.0D0,-0.01D0, var_a, var_b)
          x1=refDep(lans,var_a, var_b, hf, rhom, dt, refDepth) !refreezing depth
	  refDepth=x1


	else
        refDepth=0.0
      endif


      !     write(6,*)"Tave in qfm",tave
       TSURF = SRFTMP(QSI,A,QLI,QP,EA,TA,TAVE,TK,PR,RA,CP,RHO,                           &
     &            RKN,HNEU,ES,SBC,CS,RS,W,QNETOB,IRADFL,WS,Z,G,FC,                       &
     &            fstab,mtime,param,iTsMethod,dt,ub,refDepth,                            &
     &            smelt)   !yjs Change to add model control
      ! 5/4/04 surface melt smelt
      !     write(6,*)tsurf
       qle=(1-FC)*es*sbc*(tsurf+tk)**4
       qlnet=qli-qle

       CALL TURBFLUX(PR,RA,TA,TK,tsurf,z,g,cp,rkn,WS,EA,RHOW,                           &
     &                 HNEU,QH,QE,E,fstab) !yjs Add a fraction to reduce the evaporation after snow gone
       MR=FMELT(UB,RHOW,W,HF,LC,RID,KS,PRAIN)
      !      MR in m/hr
       QM=MR*RHOW*(HF+(tave-to)*cw)        !  Include advection of
      !       meltwater/rain that is at tave so that the model does not
      !       crash when there is no snow and it rains.
      !      QM in kj/m2/hr
      !dgt 5/4/04  Add surface melt that runs off to melt runoff and melt energy so that energy and mass balances
      ! are still correct
       mr=mr+smelt/(hf*rhow)*gsurf
	   qm=qm+smelt*gsurf
       IF(IRADFL.EQ.0) THEN
          QNET = QSI*(1.0D0-A)+QLNET
       ELSE
          QNET = QNETOB
       ENDIF
       Q = QNET+QP+QG+QH+QE-QM
       FM=PRAIN+PS-MR-E
       RETURN
       END subroutine

      !**************************** TAVG () ***************************
      !   Calculates the average temperature of snow and interacting
      !   soil layer

      real*8 FUNCTION  TAVG(UB,W,RHOW,CS,TO,RHOG,DE,CG,HF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SNHC = RHOW*W*CS
      SHC  = RHOG*DE*CG
      CHC  = SNHC+SHC

      !    SNHC = Snow heat capacity
      !    SHC  = Soil heat capacity
      !    CHC  = Combined heat capacity

      IF(UB.LE.0.) THEN
         TS=UB/CHC
      ELSE
         AL=UB/(RHOW*HF)
         IF(W.GT.AL) THEN
            TS=TO
         ELSE
            TS=(UB-W*RHOW*HF)/CHC
         ENDIF
      ENDIF
      TAVG=TS

      RETURN
      END function



      !**************************** SRFTMP () ***************************
      !   Computes the surface temperature of snow
      real*8 FUNCTION SRFTMP(QSI,A,QLI,QPIN,EA,TA,TAVE,TK,PR,RA,CP,RHO                          &
     &   ,RKN,HNEU,ES,SBC,CS,RS,W,qnetob,IRADFL,WS,Z,G,FC,fstab,mtime                          &
     &   ,param,iTsMethod,dt,ub,refDepth,smelt)	!yjs add three value to reflect model control changes
      ! 5/4/04 surface melt smelt
      !
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	real*8 param(*)
	real*8 mtime(*)
	integer iTsMethod
      ! This version written on 4/23/97 by Charlie Luce solves directly the
      ! energy balance equation using Newtons method - avoiding the linearizations
      ! used previously.  The derivative is evaluated numerically over the range
      ! ts to fff*ts  where fff = 0.999
      data tol,nitermax/0.05D0,10/
	  fff=(273D0-tol)/273.D0  ! Scale the difference used for derivatives by the tolerance
      QSN = QSI*(1.0D0-A)
      !   To ensure all temperatures in kelvin
      TAK = TA+TK
      TAVEK  = TAVE+TK
      DENSA  = PR/(RA*TAK)     ! Density of Air in kg/m3 if PR in Pa
      DENS = RHO
      !c
      qp=qpin            ! store input variable locally without changing global value
      !  debugging
      !    if(mtime(1) .eq. 1948 .and. mtime(2) .eq. 11
      !   + .and. mtime(3) .eq. 18 .and. mtime(4) .gt. 12 .and.
      !   + mtime(4) .lt. 13 .and. mtime(5) .eq. 39)then
      !	mtime(1)=mtime(1)
      !endif


      if ((w.le.0.) .and. (qp.gt.0.)) qp=0.
      !c
      !    ignore the effect of precip advected
      !    energy on the calculation of surface temperature when there is no snow.
      !    Without this ridiculously high temperatures can result as the model
      !    tries to balance outgoing radiation with precip advected energy.
      !c
       Ts = TaK                             ! first approximation
       ER=tol*2D0    ! so that it does not end on first loop
       niter = 0
 1     if(ER.gt.tol.and.niter.lt. nitermax)then
           Tslast = Ts
           F1 = surfeb(ts,RKN,WS,TAK,QP,DENSA,                                    &
     &      PR,EA,DENS,RS,TAVEK,qsn,qli,fc,qnetob,iradfl,                            &
     &	  mtime,param,iTsMethod,w,dt,ub,refDepth)
           f2 = surfeb(fff*ts,RKN,WS,TAK,QP,DENSA,                                &
     &      PR,EA,DENS,RS,TAVEK,qsn,qli,fc,qnetob,iradfl,                            &
     &	  mtime,param,iTsMethod,w,dt,ub,refDepth)
           Ts = Ts - ((1.-fff) * Ts * F1) / (F1 - F2)
	     if(ts .lt. tak - 50D0)go to 11
	     ER = abs(Ts - Tslast)
           niter=niter+1
	     go to 1   ! loop back and iterate
	 endif

	if(er.le.tol) goto 10  ! The solution has converged

      ! If still not converged use bisection method
 11    Tlb = TaK - 20.        ! First guess at a reasonable range
	 Tub = Tak + 10.
	 Flb = surfeb(tlb,RKN,WS,TAK,QP,DENSA,                                         &
     &      PR,EA,DENS,RS,TAVEK,qsn,qli,fc,qnetob,iradfl,                            &
     &	  mtime,param,iTsMethod,w,dt,ub,refDepth)
	 Fub = surfeb(tub,RKN,WS,TAK,QP,DENSA,                                         &
     &      PR,EA,DENS,RS,TAVEK,qsn,qli,fc,qnetob,iradfl,                            &
     &	  mtime,param,iTsMethod,w,dt,ub,refDepth)
	 ibtowrite=0
       if(Flb*fub .gt. 0.D0)then   ! these are of the same sign so the range needs to be enlarged
		Tlb= TaK - 150.    ! an almost ridiculously large range - solution should be in this if it exists
		Tub= Tak + 100.
		Flb = surfeb(tlb,RKN,WS,TAK,QP,DENSA,                                      &
     &      PR,EA,DENS,RS,TAVEK,qsn,qli,fc,qnetob,iradfl,                            &
     &	  mtime,param,iTsMethod,w,dt,ub,refDepth)
		Fub = surfeb(tub,RKN,WS,TAK,QP,DENSA,                                      &
     &      PR,EA,DENS,RS,TAVEK,qsn,qli,fc,qnetob,iradfl,                            &
     &	  mtime,param,iTsMethod,w,dt,ub,refDepth)
	    ibtowrite=1
          if(Flb*fub .gt. 0.)then   ! these are of the same sign so no bisection solution
			write(6,*)'Bisection solution failed with large range'
			write(6,*)'Date: ',mtime(1),mtime(2),mtime(3)
			write(6,*)'Time: ',mtime(4)
			write(6,*)'Model element: ',mtime(5)
			write(6,*)'A surface temperature of 273 K assumed'
	        ts=tk
			go to 10
	    else
			write(6,*) 'Bisection surface solution with large range'
			write(6,*)'Date: ',mtime(1),mtime(2),mtime(3)
			write(6,*)'Time: ',mtime(4)
			write(6,*)'Model element: ',mtime(5)
			write(6,*)                                                                          &
     &	'This is not a critical problem unless it happens frequently'
	         write(6,*)                                                                       &
     &     'and solution below appears incorrect'
	     endif
        else
      !	write(6,*)
      !   +    'Bisection surface temperature solution'
      !	write(6,*)'Date: ',mtime(1),mtime(2),mtime(3)
      !	write(6,*)'Time: ',mtime(4)
      !	write(6,*)'Model element: ',mtime(5)
      !    write(6,*)'This is not a critical problem'
        endif
      !   Here do the bisection
       niter = log((Tub-Tlb)/tol)/log(2.)   ! Number of iterations needed for temperature tolerance
	  do iter=1,niter
           ts = 0.5*(tub+tlb)
           F1 = surfeb(ts,RKN,WS,TAK,QP,DENSA,                                      &
     &      PR,EA,DENS,RS,TAVEK,qsn,qli,fc,qnetob,iradfl,                           & !25
     &	  mtime,param,iTsMethod,w,dt,ub,refDepth)	!yjs add three value to reflect model control changes))
		 if(f1.gt.0.0) then  ! This relies on the assumption (fact) that this is a monotonically decreasing function
			tlb=ts
	     else
			tub=ts
	     endif
       enddo
	  if(ibtowrite .eq. 1)then
		write(6,*)'Surface temperature: ',ts,' K'
		write(6,*)'Energy closure: ',f1
		write(6,*)'Iterations: ',niter
	  endif

 10    Ts = Ts - Tk
       IF(W.GT.0.D0.AND.TS.GT.0.D0) THEN
      ! 5/4/04 surface melt smelt
          SRFTMP = 0.0
	    smelt=surfeb(srftmp+Tk,RKN,WS,TAK,QP,DENSA,                                 &
     &      PR,EA,DENS,RS,TAVEK,qsn,qli,fc,qnetob,iradfl,                           & !25
     &	  mtime,param,iTsMethod,w,dt,ub,refDepth)
      ! 5/4/04 surface melt smelt is the energy not accommodated by conduction into the snow
      !so it results in surface melt which then infiltrates and adds energy to the snowpack
      !through refreezing
      !This change added to capture this quantity when the model is used for a glacier where
      !surface melt may run off rather than infiltrate.
      !For modeling of glaciers in Antarctica by Mike Gooseff
       ELSE
          SRFTMP = TS
	    smelt=0.0
      !No surface melt this case
       ENDIF
       RETURN
       END function
      !******************************************************************
       real*8 function surfeb(tss,RKN,WS,TAK,QP,DENSA,                          &
     &  PR,EA,DENS,RS,TAVEK,qsn,qli,fc,qnetob,iradfl,                               &
     &      mtime,param,iTsMethod,w,dt,ub,refDepth)	!yjs add three value to reflect model control changes)
      !    function to evaluate the surface energy balance for use in solving for
      !    surface temperature
      !    DGT and C Luce 4/23/97
      USE COMM
	  IMPLICIT DOUBLE PRECISION (a-z)
	  real*8 param(*)
	  real*8 mtime(*)
	  integer iTsMethod, iradfl
	  real*8 LanE_Ze, LanE_de, LanE_ze2,LanE_de2

	  !common /ts_save/ ts_old, tave_old, Ts_ave, Tave_ave

      !  Constant data set
      ! CP: declared in module comm


      !  End of constant declaration

      !Parameters

	  zs=w*rhow/rho; ze = 0D0

        RKIN=RKINST(RKN,WS,TAK,tss,Z,G,fstab)

      !  07/25/02   at Boise.  To make the UEB2 work under Unix/linux system the fancy stuff like "Select case" shall not be used
	    !select case(iTsMethod)			!Four choice of the surface temperature modeling
	    !case (1)						!1. Old snow, use calibrated snow surface heat conductance
          if(iTsMethod .eq. 1) then
			qcs = DENS*CS*RS*(tss-TAVEK) !2. Revised scheme LanE/Ze of the snow surface
										!3. Force restore approach
		!case (2)						!4. Modified force restore approach.
          elseif (iTsMethod .eq. 2) then
			fKappaS = LanS/(rho*cs)

			fKappaG = LanG/(rhog*cg)

			d1 = sqrt(2*fKappaS/w1day)
			if(zs .ge. rd1*d1) then

				LanE_Ze=LanS/(rd1*d1)
	        else

				LanE_Ze=LanE(LanS, LanG, Zs, rho, rhog, cs, cg, rd1,                           &
     &				ze, w1day) !Cyjs   call the subroutine to update the heat conductivity. LanE()
				LanE_Ze=LanE_Ze/ze
			end if

			qcs= LanE_Ze*(tss-TAVEK)
		elseif (iTsMethod .eq. 3) then !case (3)
			fKappaS = LanS/(rho*cs)
			fKappaG = LanG/(rhog*cg)

			d1 = sqrt(2*fKappaS/w1day)

			if(zs .ge. rd1*d1) then
				LanE_Ze=LanS/(rd1*d1)
				Ze=rd1*d1
	        else

				LanE_Ze=LanE(LanS, LanG, Zs, rho, rhog, cs, cg, rd1,                            &
     &				ze, w1day) !Cyjs   call the subroutine to update the heat conductivity. LanE()
				LanE_Ze=LanE_Ze/ze
			end if

			dee=ze/rd1
			LanE_de=LanE_ze/dee*ze
			qcs= LanE_de*(tss-Ts_old)/(w1day*dt)+LanE_Ze*(tss-TAVEK)


		else !case (4)   !change to all default cases. If not for model comparison
			fKappaS = LanS/(rho*cs)
			fKappaG = LanG/(rhog*cg)

			d1 = sqrt(2*fKappaS/w1day)
			dlf = sqrt(2*fKappaG/wlf)

			if(zs .ge. rd1*d1) then
				LanE_Ze=LanS/(rd1*d1)
				Ze=rd1*d1
			else
				LanE_Ze=LanE(LanS, LanG, Zs, rho, rhog, cs, cg, rd1,                            &
     &				ze, w1day) !Cyjs   call the subroutine to update the heat conductivity. LanE()
				LanE_Ze=LanE_Ze/ze
			end if

			if(zs .ge. rd1*dlf) then
				LanE_Ze2=LanS/(rd1*dlf)
				Ze2=rd1*dlf
			else
				LanE_Ze2=LanE(LanS, LanG, Zs, rho, rhog, cs, cg, rd1,                            &
     &				ze2, wlf) !Cyjs   call the subroutine to update the heat conductivity. LanE()
				LanE_Ze2=LanE_Ze2/ze2
			end if

			dee=ze/rd1
			LanE_de=LanE_ze/dee*ze
			de2=ze2/rd1
			LanE_de2=LanE_ze2/de2*ze2

			if(ub.le.0.0D0 .or. refDepth.le.0.0)then
			  qcs= LanE_de*(tss-Ts_old)/(w1day*dt)+LanE_Ze*(tss-Ts_Ave)+                          &
     &				LanE_de2*(Ts_ave-Tave_Ave)
	        elseif(refDepth .gt. rd1*d1) then
			  qcs= LanE_de*(tss-Ts_old)/(w1day*dt)+LanE_Ze*(tss-Ts_Ave)+                          &
     &				LanE_de2*(Ts_ave-Tave_Ave)
              else
                 qcs=lanE_ze*ze*(tss-tk)/refDepth
	        endif


		!End select
            endif


       	surfeb = QP + RKIN*DENSA*CP*(TAK - tss)                                          &
     &          +(HNEU*DENSA*0.622D0*RKIN)/PR*(EA-SVP(tss-Tk))                          &
     &          -qcs                !yjs This is a change. The original has only one selection of case 1.
	    If(Iradfl.eq.0)then
			surfeb = surfeb + QSN+QLI-(1.D0-FC)*ES*SBC*tss**4D0
      	ELSE
     			surfeb = surfeb + qnetob
      	ENDIF
        return
	  end function
      !******************************************************************
       real*8 function QcEst(tss,RKN,WS,TAK,QP,DENSA,                                &
     &      PR,EA,DENS,RS,TAVEK,qsn,qli,fc,qnetob,iradfl,                              &
     &      mtime,param,iTsMethod,w,dt)	!yjs add three value to reflect model control changes)
      !    function to estimate surface heat conduction for use in solving for
      !    surface temperature
      !    DGT and C Luce 4/23/97
       USE COMM
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)

	real*8 param(*)
	real*8 mtime(*)
	integer iTsMethod
	real*8 LanE_Ze, LanE_de, LanE_ze2,LanE_de2

      !  Constant data set
      ! CP: declared in module comm
      !  End of constant declaration

      !Parameters
	  zs=w*rhow/rho

        RKIN=RKINST(RKN,WS,TAK,tss,Z,G,fstab)

       	qcEst = QP + RKIN*DENSA*CP*(TAK - tss)                                          &
     &          +(HNEU*DENSA*0.622D0*RKIN)/PR*(EA-SVP(tss-Tk))

	    If(Iradfl.eq.0)then
			qcEst = qcEst + QSN+QLI-(1.D0-FC)*ES*SBC*tss**4D0
      	ELSE
     			qcEst = qcEst + qnetob
      	ENDIF
        return
	  end function
      !**************************** SRFTMPO () ***************************
      !   Computes the surface temperature of snow
      REAL*8 FUNCTION SRFTMPO(QSI,A,QLI,QPIN,EA,TA,TAVE,TK,PR,RA,CP,               &
     &  RHO,RKN,HNEU,ES,SBC,CS,RS,W,qnetob,IRADFL,WS,Z,G,FC,fstab)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 tub, tlb
      dimension tint(0:10,2)
      DATA ncall,tol/0D0,0.05D0/
      ncall = ncall+1
      NITER = 10
      TSTAR = TA
      qp=qpin    ! store input variable locally without changing global value
      if((w.le.0.) .and. (qp.gt.0.)) then
      qp=0.
      endif
         ! ignore the effect of precip advected
      !    energy on the calculation of surface temperature when there is no snow.
      !    Without this ridiculously high temperatures can result as the model
      !    tries to balance outgoing radiation with precip advected energy.
      TSURF = SNOTMP(TSTAR,QSI,A,QLI,QP,EA,TA,TAVE,PR,                             &
     &        RKN,rs,qnetob,IRADFL,WS,FC,mtime) !YJS pass modeling time
      !    The first calculation gets the solution on the right side of ta
      !    since it uses neutral stability and the linearization (equivalent
      !    to Newton-Rhapson) will move in the direction of a solution in the case
      !    of a well behaved function.  See Notebook 8/3/94 for further elaboration.
      !    Use this to place bounds on surface temperature for iteration.
      if(tsurf.gt.ta)then
        tlb=ta
        tub=ta+30.   !  Max upper bound 30 C warmer than surface
        if(tsurf.gt.tub)tsurf=(tlb+tub)*.5  ! sometimes tsurf is outside these bounds
      else
        tlb=ta-30.
        tub=ta
        if(tsurf.lt.tlb)tsurf=(tlb+tub)*.5
      endif
      tint(0,1)=tstar
      tint(0,2)=tsurf
      ! Now iterate
      tstar=tsurf
      DO 10 I=1,NITER
      TSURF = SNOTMP(TSTAR,QSI,A,QLI,QP,EA,TA,TAVE,PR,                                  &
     &        RKN,rs,qnetob,IRADFL,WS,FC,mtime) !YJS pass modeling time

      !    write(6,*)ta,tstar,tsurf
      tint(i,1)=tstar
      tint(i,2)=tsurf
      if(tlb.le.tsurf .and. tsurf .le. tub)then
         if(tsurf.gt.tstar)then
           tlb=tstar   ! increasing so can increase lower bound.
         else
           tub=tstar
         endif
       else if(tsurf .gt. tub)then   ! upper bound overshot
         tlb=tstar                   ! increase lower bound and
         tsurf=(tstar+tub)*.5D0        ! guess at solution halfway
       else    ! tsurf .lt. tlb  here, i.e. lower bound overshot
         tub=tstar
         tsurf=(tstar+tlb)*.5D0
       endif
      !  Check for convergence
       IF(ABS(TSURF-TSTAR).LT. tol) THEN
          GO TO 20
       ELSE
          TSTAR = TSURF
       ENDIF
 10    CONTINUE
      !     write(6,*) 'slipped through the loop in SRFTMP()',ncall
      !     do 15 i = 0,niter
      !15      write(6,*)tint(i,1),tint(i,2)
      !     write(6,*)"tsurf",tsurf
      !  Newton rhapson not converging so use bisection
       f1=sntmpb(tlb,QSI,A,QLI,QP,EA,TA,TAVE,TK,PR,RA,CP,RHO,                          &
     &         RKN,HNEU,ES,SBC,CS,RS,qnetob,IRADFL,WS,Z,G,FC,fstab)
       f2=sntmpb(tub,QSI,A,QLI,QP,EA,TA,TAVE,TK,PR,RA,CP,RHO,                          &
     &         RKN,HNEU,ES,SBC,CS,RS,qnetob,IRADFL,WS,Z,G,FC,fstab)
       tsurf=(tlb+tub)*.5D0
       if(f1*f2 .gt.0.)then
         write(6,*)'SRFTMP has failed to find a solution',ncall                          &
     &     ,tlb,tub,tsurf
       else
         nib=(log(tub-tlb)-log(tol))/log(2.D0)
         do 16 i=1,nib
      !       write(6,*)nib,tlb,tub,tsurf
         f=sntmpb(tsurf,QSI,A,QLI,QP,EA,TA,TAVE,TK,PR,RA,CP,RHO,                          &
     &         RKN,HNEU,ES,SBC,CS,RS,qnetob,IRADFL,WS,Z,G,FC,fstab)
         if(f*f1 .gt.0.)then
           tlb=tsurf
           f1=f
         else
           tub=tsurf
      !        f2=f
         endif
         tsurf=(tlb+tub)*.5D0
 16      continue
      !       write(6,*)"Bisection",nib,tsurf
       endif
 20     IF(W.GT.0..AND.TSURF.GT.0.D0) THEN
          SRFTMPO = 0.0
       ELSE
          SRFTMPO = TSURF
       ENDIF
      !     write(6,*)"Final",srftmp
       RETURN
       END function

      !**************************** SNOTMP () ******************************
      !   Function to compute surface temperature using Penman/surface
      !   resistance analogy and equilibrium approach

      !      FUNCTION SNOTMP(TSTAR,QSI,A,QLI,QP,EA,TA,TAVE,TK,PR,RA,CP,RHO,
      !     &        RKN,HNEU,ES,SBC,CS,RS,qnetob,IRADFL,WS,Z,G,FC,fstab)
      !REAL*8 FUNCTION SNOTMP(TSTAR,QSI,A,QLI,QP,EA,TA,TAVE,TK,PR,                          &
      !&    RA,CP,RHO,RKN,HNEU,ES,SBC,CS,RS,qnetob,IRADFL,WS,Z,G,FC,fstab,                      &
      !&    mtime,lanS,lanG) !add a modeltime
      REAL*8 FUNCTION SNOTMP(TSTAR,QSI,A,QLI,QP,EA,TA,TAVE,PR,                           &
     &    rkn,rs,qnetob,IRADFL,ws,fc,mtime) !add a modeltime
      use comm
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	  real*8 LanmdaE
      QSN = QSI*(1.0-A)
      !   To ensure all temperatures in kelvin
      TAK = TA+TK
      TSTARK = TSTAR+TK
      TAVEK  = TAVE+TK
      DENSA  = PR/(RA*TAK)
      !   DENSA in kg/m3 if PR in Pa
      !    CP1 = CP/1000.0D0
      !   CP1 in kj/kg/oK


      DENS = RHO
      RKIN=RKINST(RKN,WS,TAK,TSTARK,Z,G,fstab)

	LanmdaE=LanE(LanS, LanG, Zs, DENS, rhog, cs, cg, r, ze, w1day)

      !      UPPER = QP+(DENSA*CP*TAK)*RKIN
      !     &        -(HNEU*DENSA*0.622D0*RKIN)/PR*(SVPI(TSTAR)
      !     &        -EA-DELTA(TSTAR)*TSTARK)
      !     &        +DENS*CS*RS*TAVEK
      !      DEN = DENS*CS*RS+(DENSA*CP)*RKIN+(DELTA(TSTAR)*
      !     &       HNEU*DENSA*.622D0*RKIN)/PR

      UPPER = QP+(DENSA*CP*TAK)*RKIN                                                &
     &        -(HNEU*DENSA*0.622D0*RKIN)/PR*(SVPI(TSTAR)                            &
     &        -EA-DELTA(TSTAR)*TSTARK)                                               &
     &        +LanmdaE*TAVEK                     !change to combined layer (snow and soil)
      DEN = LanmdaE+(DENSA*CP)*RKIN+(DELTA(TSTAR)*                                  &
     &       HNEU*DENSA*.622*RKIN)/PR
      IF(IRADFL.EQ.0) THEN
         UPPER = UPPER + QSN+QLI+3.0D0*(1.-FC)*ES*SBC*TSTARK**4D0
         DEN = DEN+4.0D0*(1.-FC)*ES*SBC*TSTARK**3D0
      ELSE
         UPPER = UPPER + qnetob
      ENDIF

      SNOTMP = UPPER/DEN-TK

      RETURN
      END function
      !************************ sntmpb () ******************************
      !   Function to compute surface temperature using bisection

      REAL*8 FUNCTION sntmpb(TSTAR,QSI,A,QLI,QP,EA,TA,TAVE,TK,PR,RA,                          &
     &   CP,RHO,RKN,HNEU,ES,SBC,CS,RS,qnetob,IRADFL,WS,Z,G,FC,fstab)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      QSN = QSI*(1.0-A)
      !   To ensure all temperatures in kelvin
      TAK = TA+TK
      TSTARK = TSTAR+TK
      TAVEK  = TAVE+TK
      DENSA  = PR/(RA*TAK)
      !   DENSA in kg/m3 if PR in Pa
      !    CP1 = CP/1000.0
      !   CP1 in kj/kg/oK
      DENS = RHO
      RKIN=RKINST(RKN,WS,TAK,TSTARK,Z,G,fstab)

      UPPER = QP+(DENSA*CP*TAK)*RKIN                                              &
     &        -(HNEU*DENSA*0.622*RKIN)/PR*(SVPI(TSTAR)                            &
     &        -EA) +DENS*CS*RS*TAVEK
      DEN = DENS*CS*RS+(DENSA*CP)*RKIN
      IF(IRADFL.EQ.0) THEN
         UPPER = UPPER + QSN+QLI-(1.D0-FC)*ES*SBC*TSTARK**4D0
      ELSE
         UPPER = UPPER + qnetob
      ENDIF
      sntmpb=upper-den*tstark
      RETURN
      END function

      !**************************** RKINST() **********************
      REAL*8 function rkinst(rkn,ws,ta,ts,z,g,fstab)
      !   function to calculate no neutral turbulent transfer coefficient using the
      !   richardson number correction.
      IMPLICIT DOUBLE PRECISION (a-z)
      if(ws.le.0.)then
        rkinst=0.    !  No wind so no sensible or latent heat fluxes.
      else
        rich=g*(ta-ts)*z/(ws*ws*ta)    ! ta must be in K
        if(rich.ge.0.)then
        rkinst=rkn/(1.D0+10.D0*rich)
		! Cyjs PhiMH=1+10*rich    Change back
        else
      !          rkinst=rkn*(1D0-10.D0*rich)
		!PhiMH=(1-16*rich)**(-.75)
		rkinst=rkn*min(3.D0,(1D0-16D0*rich)**(.75D0))
        endif
      endif
      !Linear damping of stability correction through parameter fstab

      rkinst=rkn+fstab*(rkinst-rkn)
		!if(PhiMH .ne. 0.0) rkinst=rkn/PhiMH !!!!Change back to the original design purpose
      return
      end function


      !******************************* DELTA () ****************************
      !   Function to compute gradient of saturated vapour pressure,
      !   temperature function over ice
      !   Uses Lowe (1977) polynomial
      REAL*8 FUNCTION DELTA(T)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension a(7)
      data a/0.5030305237D0,0.0377325502D0,0.001267995369D0,                          &
     &   2.477563108D-5,3.005693132D-7,2.158542548D-9,                          &
     &   7.131097725D-12/
      delta=a(1)+t*(a(2)+t*(a(3)+t*(a(4)+t*(a(5)+t*(a(6)+t*a(7))))))
      delta=delta*100.    ! convert from mb to Pa
      RETURN
      END function

      !************************** PARTSNOW () ************************
      !   Partitioning of precipitation into rain and snow

      REAL*8 FUNCTION PARTSNOW(P,TA,TR,TS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IF(TA.LT.TS) THEN
        PARTSNOW=P
      ELSEIF(TA.GT.TR) THEN
        PARTSNOW=0.0D0
      ELSE
        PARTSNOW=P*(TR-TA)/(TR-TS)
      ENDIF

      RETURN
      END function

      !*************************** QPF () ***************************
      !   Calculates the heat advected to the snowpack due to rain

      REAL*8 FUNCTION QPF(PR,TA,TO,PS,RHOW,HF,CW,CS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IF(TA.GT.TO) THEN
         TRAIN=TA
         TSNOW=TO
      ELSE
         TRAIN=TO
         TSNOW=TA
      endif
      QPF=PR*RHOW*(HF+CW*(TRAIN-TO))+PS*RHOW*CS*(TSNOW-TO)
      RETURN
      END function

      !************************ TURBFLUX () ***************************
      !   Calculates the turbulent heat fluxes (sensible and latent
      !   heat fluxes) and condensation/sublimation.

      SUBROUTINE TURBFLUX(PR,RA,TA,TK,TS,z,g,cp,rkn,WS,EA,                          &
     &                    RHOW,HNEU,QH,QE,E,fstab)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

       real*8 RHOA,QH,ES,EA,QE,E

      tak=ta+tk
      tsk=ts+tk
      rkin=rkinst(rkn,ws,tak,tsk,z,g,fstab)
      RHOA=PR/(RA*(TAK))
      !   RHOA in kg/m3
      QH=RHOA*(TA-TS)*cp*rkin
      ES=SVPI(TS)
	QE=0.622D0*hneu/(RA*(TAK))*rkin*(EA-ES)
      E=-QE/(RHOW*HNEU)
      !   E in  m/hr

      RETURN
      END subroutine

      !*********************** FMELT () *************************
      !   Calculates the melt rate and melt outflow

      REAL*8 FUNCTION FMELT(UB,RHOW,W,HF,LC,RID,KS,PRAIN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      real*8 LC,KS

      !     write(*,*) 'I am in FMELT!!!!'


      UU=0.0
      IF(UB.LT.0.D0) THEN
         FMELT=0.0
      ELSEIF(W.LE.0.0D0) THEN
         FMELT=PRAIN
         IF(PRAIN.LE.0.D0) THEN
            FMELT=0.0
         ENDIF
      ELSE
         UU=UB/(RHOW*W*HF)
      !                            liquid fraction
         IF(UU.GT.0.99D0) THEN
            UU=0.99D0
      !                            TO ENSURE HIGH MELT RATE WHEN ALL LIQUID
      ENDIF
      IF((UU/(1-UU)).LE.LC) THEN
         SS=0.0D0
      ELSE
         SS=(UU/((1- UU)*RID)-LC/RID)/(1-LC/RID)
      ENDIF
      FMELT=KS*SS**3D0
      ENDIF
      IF(FMELT.LT.0.0) THEN
          WRITE(*,*)'FMELT is NEGATIVE!'
          STOP
       ENDIF

      !     write(*,*) ' I am leaving FMELT!!!'

       RETURN
       END function

      !*************************** SVP () *****************************
      !   Calculates the vapour pressure at a specified temperature over water or ice
      !   depending upon temperature.  Temperature is celsius here.
      REAL*8 FUNCTION SVP(T)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      if(t .ge. 0.D0)then
        svp=svpw(t)
      else
        svp=svpi(t)
      endif
      return
      end function
      !*************************** SVPI () *****************************
      !   Calculates the vapour pressure at a specified temperature over ice.
      !   using polynomial from Lowe (1977).
      REAL*8 FUNCTION SVPI(T)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SVPI=6.109177956D0 + t * (0.503469897D0 + t * (0.01886013408D0 +                           &
     & t * (0.0004176223716D0 + t * (5.82472028D-06 + t *                           &
     &   (4.838803174D-08 + t * 1.838826904D-10)))))
      SVPI=SVPI*100   ! convert from mb to Pa
      RETURN
      END function

      !*************************** SVPW () *****************************
      !   Calculates the vapour pressure at a specified temperature over water
      !   using polynomial from Lowe (1977).
      REAL*8 FUNCTION SVPW(T)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SVPW=6.107799961D0 + t * (0.4436518521D0 + t * (0.01428945805D0 +                           &
     & t * (0.0002650648471D0 + t * (3.031240936D-06 + t *                           &
     & (2.034080948D-08 + t * 6.136820929D-11)))))
      SVPW=SVPW*100D0   ! convert from mb to Pa
      RETURN
      END function

      !************************* PREHELP () ***************************
      !    Routine to correct energy and mass fluxes when
      !    numerical overshoots dictate that W was changed in
      !    the calling routine - either because W went negative
      !    or due to the liquid fraction being held constant.

       SUBROUTINE PREHELP(W1,W,DT,FM,FM1,fac,PS,PRAIN,E,RHOW,HF,Q,QM,                          &
     & MR,qe,hneu)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       real*8 MR,DT

      ! The next statement calculates the value of FM given
      !  the W and w1 values
       FM = (w1-W)/DT*fac-FM1
      ! The next statement calculates the changed MR and E due to the new FM.
      !  FM was = PRAIN+PS-MR-E
      !  Here the difference is absorbed in MR first and then E if mr < 0.
      !
       MR = MAX( 0.0D0 , (PS + PRAIN - FM - E))
       E = Ps + Prain - FM - MR
      !  Because melt rate changes the advected energy also changes.  Here
      !   advected and melt energy are separated,
       QOTHER = Q + QM - QE
      !   then the part due to melt recalculated
       QM = MR*RHOW*HF
      !   then the part due to evaporation recalculated
       QE = -E*RHOW*HNEU
      !   Energy recombined
       Q = QOTHER - QM + QE
       RETURN
       END subroutine

      !***************************** ALBEDO () *****************************
      !   Function to calculate Albedo
      !   BATS Albedo Model (Dickinson et. al P.21)
      REAL*8 FUNCTION ALBEDO(tausn,coszen,d,aep,abg,avo,airo)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      B = 2.0D0
      CS = 0.2D0
      !    AVO = 0.95
      CN = 0.5D0
      !    AIRO = 0.65

      FAGE = tausn/(1.0D0+tausn)

      IF(coszen.LT.0.5D0) THEN
         FZEN = 1.0D0/B*((B+1.0D0)/(1.0D0+2.0D0*B*coszen)-1.0D0)
      ELSE
         FZEN = 0.0D0
      ENDIF
      AVD = (1.0D0-CS*FAGE)*AVO
      AVIS = AVD+0.4D0*FZEN*(1.0D0-AVD)
      AIRD = (1.0D0-CN*FAGE)*AIRO
      ANIR = AIRD+0.4D0*FZEN*(1.0D0-AIRD)
      ALBEDO = (AVIS+ANIR)/2.0D0
      if(d.lt.aep)then   ! need to transition albedo to a bare ground value
        rr=(1.D0-d/aep)*exp(-d*.5D0/aep)
        albedo=rr*abg+(1.D0-rr)*albedo
      endif
      RETURN
      END function

      !******************************** AGESN () *****************************
      !   Function to calculate Dimensionless age of snow for use in
      !   BATS Albedo Model (Dickinson et. al P.21)

      Subroutine agesn(tausn,dt,ps,tsurf,tk,dNewS)
      !      Subroutine agesn(tausn,dt,ps,tsurf,tk)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      tsk=tsurf+tk  ! Express surface temperature in kelvin
      R1 = EXP(5000.0D0*(1.0D0/TK - 1.0D0/TSK))
      R2 = R1**10D0
      IF(R2.GT.1.0) THEN
         R2 = 1.0D0
      ENDIF
      R3 = 0.3D0
      ! Dickinson p 23 gives decay as DT*(R1+R2+R3)/tau0  with
      !  tau0 = 10**6 sec.  Here 0.0036 = 3600 s/hr * 10**-6 s**-1
      !  since dt is in hours.
      !      tausn = max((tausn+0.0036D0*(R1+R2+R3)*DT)*
      !     &            (1.0D0 - 100.0D0*PS*DT),0.D0)
      tausn = tausn+0.0036*(R1+R2+R3)*DT
	  if(ps.gt.0.0D0) then
	  if(dNewS .gt. 0.0) then
	    tausn = max(tausn*(1-ps*dt/dNewS),0.D0)
	  else
	    tausn = 0.0D0
	  endif
	  endif

      RETURN
      END subroutine


      !***************************** LanE () *************************************
      !Function to get the LanE which is the thermal conductivity by ze

	real*8 function LanE(LanS, LanG, Zs, rho, rhog, cs, cg, r, ze,                           &
     &     w1day)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	real*8 LanS, LanG, Zs, Dlf


	 fKappaS = LanS/(rho*cs)
	 fKappaG = LanG/(rhog*cg)

	 d1 = sqrt(2D0*fKappaS/w1day)
	 d2 = sqrt(2D0*fKappaG/w1day)

	 LanE=min(r,zs/d1)*d1/LanS + max(0.D0,(r-zs/d1))*d2/LanG

	 ze=min(r,zs/d1)*d1+max(0.D0,(r-zs/d1))*d2

	 if (LanE .ne. 0.D0) LanE = 1/LanE*ze

	 return

	end function
      !*************************** refDep() **************************************
      !   function to calculate the refreezing depth
	real*8 function refDep(flans, a, b, hf, rhom, dt, x1 )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	temp=flans*flans+2*b*(-a*flans*dt/(hf*rhom)+flans*x1+0.5*b*x1*x1)

	if(temp .lt. 0.0 .or. a.gt.0.0 .or. b.eq.0.0) then
	  refDep=0
	else
	  refDep=(-flans+sqrt(temp))/b

	endif

	return
	end function

      !************************** Grad () ********************************
      !Linear simple function to calculate the gradient of Q and T

      subroutine Grad(qc1,qc2,t1,t2, a, b)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	if((t2-t1) .ne. 0.0) then
	  b=(qc2-qc1)/(t1-t2)
	  a=qc1+b*t1
	endif

	return
	end subroutine

      !****************** Caculate the daily average  *********************
      !  Calculate the daily average value
      !  n number of records, a minimum value -100 or some
	real*8 function daily_ave(backup, n, a)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	real*8 backup(*)
	sum = 0.D0
	count=0.D0
	do 1 i = 1, n
	if(backup(i).gt.a) then
		sum=sum+backup(i)
		count=count+1.
	end if
 1	continue
	if(count.ne.0) then
		daily_ave=sum/count
	else
		daily_ave=0.D0
	end if
	return
	end function
      !File  snowx.f
      !Subroutines and function subprograms for the Utah Energy Balance
      !Snow Accumulation and Melt Model.
      !David G. Tarboton, Utah Water Research Laboratory, Utah State University
      !Version 2.1
      !Last Change 5/4/04 to accommodate glacier surface runoff.
      !C
      !This program was written and prepared under agreement with and funding
      !by the U.S. Government , and therefore is in the public domain and not
      !subject to copyright.

      !***********READ Bristow Campbell Parameter values **************

      SUBROUTINE bcparm(dtbar,bca,bcc,bcfile)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*200 bcfile
      real*8 dtbar(12)
      OPEN(8,FILE=bcfile,STATUS='OLD')
      READ(8,*)bca,bcc
 1    read(8,*,end=20)month,val
      dtbar(month)=val
      go to 1
 20   CLOSE(8)
      RETURN
      END subroutine

      !***********READ PARAMETER or state variable VALUES**************

      SUBROUTINE readvals(vals,vfile,n)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*200 vfile
      real*8 vals(1)
      OPEN(8,FILE=vfile,STATUS='OLD')
      READ(8,*)(vals(i),i=1,n)
      CLOSE(8)
      RETURN
      END subroutine

      !****************** READ THE INITIAL CONDITIONS *********************

      !    SUBROUTINE FINITIAL(statev,icfile,nxv,YEAR,MONTH,DAY,HOUR,DT)
      !    CHARACTER*20 icfile
      !    INTEGER YEAR,DAY
      !    real*8 statev(1)
      !    OPEN(8,FILE=icfile,STATUS='OLD')
      !    READ(8,*)(statev(i),i=1,nxv)
      !    READ(8,*)YEAR,MONTH,DAY,HOUR,DT
      !    CLOSE(8)
      !    RETURN
      !    END

      !****************** READ THE site variables  *********************

      SUBROUTINE readsv(sitev,svfile,nsv,slope,azi,lat)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*200 svfile
      real*8 sitev(1),lat
      OPEN(8,FILE=svfile,STATUS='OLD')
      READ(8,*)(sitev(i),i=1,nsv)
      READ(8,*)slope,azi,lat
      CLOSE(8)
      RETURN
      END subroutine

      !*********************** UPDATETIME () ***************************
      !               Update Time for each time step

      SUBROUTINE UPDATETIME(YEAR,MONTH,DAY,HOUR,DT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER YEAR,DAY,DMON(12),DM, MONTH,I       ! 30/03/2004 ITB
	  real*8  hour, dt  ! DGT Dec 10, 2004.  Fixing ITB errors

      DATA (DMON(I),I=1,12)/31,28,31,30,31,30,31,31,30,31,30,31/
      HOUR=HOUR+DT
	DM=DMON(MONTH)
      ! check for leap years
      if(month .eq. 2)dm=lyear(year)
 10   continue
      IF(HOUR.GE.24.0) THEN
        HOUR=HOUR-24.0
        DAY=DAY+1
	go to 10
      ENDIF
 20   continue
      IF(DAY.GT.DM) THEN
        DAY=day - dm
        MONTH=MONTH+1
        IF(MONTH.GT.12) THEN
          MONTH=1
          YEAR=YEAR+1
	    DM=DMON(MONTH)
          if(month .eq. 2)dm=lyear(year)
	  endif
	  go to 20
      ENDIF
      RETURN
      END subroutine

      !*********************** lyear () ***************************
      !   function to return number of days in February checking for
      !   leap years
      function lyear(year)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer year,lyear
      IF(MOD(YEAR,4).GT.0 .or.                                                     &
     &   (mod(year,100) .eq.0 .and. mod(year,400) .ne. 0)) THEN
      !Leap years are every 4 years
      !- except for years that are multiples of centuries (e.g. 1800, 1900)
      !- except again that when the century is divisible by 4 (e.g. 1600, 2000)
      !  then it is a leap year
         lyear=28
      ELSE
         lyear=29
      ENDIF
	return
	end function

      !************************* atf () ****************************
      !   to get the atmospheric transmissivity using the Bristow and Campbell
      !   (1984) approach

      Subroutine atf(atff,trange,month,dtbar,a,c)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION dtbar(12)
      b=0.036D0*exp(-0.154D0*dtbar(month))
      atff=a*(1D0-exp(-b*trange**c))
      !    write(6,*)trange,month,a,c,dtbar(month),atf
      RETURN
      END subroutine


      !*********************** hourlyRI () **********************
      !              To get hourly radiation index

      SUBROUTINE hyri(YEAR,MONTH,DAY,HOUR,DT,SLOPE,AZI,LAT,                     &
     &                    HRI,COSZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER YEAR,DAY
      real*8 LP,LAT1,LAT
      !lp= latitude of equivalent plane in radians
      !lat1 = latitude in radians
      !lat = latitude in degrees
      !a number that speaks for itself - every kissable digit
      PI=3.141592653589793238462643383279502884197169399375105820974944592308
      CRAD=PI/180.0
      ! crad = degree to radian conversion factor
      !   CONVERT TIMES TO RADIANS FROM NOON
      T=(HOUR-12.0)*PI/12.0
      DELT1=DT*PI/12.0
      !   CONVERT angles TO RADIANS
      SLOPE1=SLOPE*CRAD
      AZI1=AZI*CRAD
       LAT1=LAT*CRAD
	FJULIAN=FLOAT(JULIAN(year,MONTH,DAY))
      D=CRAD*23.5*SIN((FJULIAN-82.0)*0.017214206321)
      !0.017214206321 is 2 pi / 365
      !D is solar declination
      LP=ASIN(SIN(SLOPE1)*COS(AZI1)*COS(LAT1)                                    &
     &   +COS(SLOPE1)*SIN(LAT1))
      !LP is latitude of equivalent plane
      !    TD=ACOS(-TAN(LAT1)*TAN(D))  This formula abandoned 1/8/04
      !    to make the code work for polar conditions
      !TD is half day length, i.e. the time from noon to sunset.  Sunrise is at -TD
      tanprod=TAN(LAT1)*TAN(D)
	  if(tanprod .gt. 1.)then
	  td=pi  ! This is the condition for perpetual light
	  else if(tanprod .lt. -1.)then
	  td=0   ! The condition for perpetual night
	  else
	  td=acos(-tanprod)  ! The condition where there is a sunrise and set
	  endif
      ! Equivalent longitude offset.  Modified on 1/8/04
      ! so that it correctly accounts for shift in longitude if equivalent
      ! plane slope goes over a pole.  Achieved using atan2.
      !    DDT=ATAN(SIN(AZI1)*SIN(SLOPE1)/(COS(SLOPE1)*COS(LAT1)
      !   *    -COS(AZI1)*SIN(SLOPE1)*SIN(LAT1)))
      ddt=atan2(SIN(AZI1)*SIN(SLOPE1),                                                   &
     &(COS(SLOPE1)*COS(LAT1)-COS(AZI1)*SIN(SLOPE1)*SIN(LAT1)))

      ! Now similar logic as before needs to be repeated for equivalent plane
      ! but with times reflecting
      TPeqp=TAN(LP)*TAN(D)
      !Keep track of beginning and end of exposure of equiv plane to sunlight
      IF(tpeqp .gt. 1.0) THEN
          TPbeg=-pi   ! perpetual light
	    tpend=pi
      ELSEif(tpeqp .lt. -1.)then
          TPbeg=0.0  ! perpetual dark
          tpend=0.0
	  else
	    tpbeg = -acos(-tpeqp)-ddt
	    tpend = acos(-tpeqp)-ddt
      ENDIF

      ! Start and end times for integration of radiation exposure
      ! need to account for both horizon, slope and time step
      T1=max(T,tpbeg,-TD)
      T2=min(T+DELT1,TD,tpend)
      !    write(6,*)t1,t2
      IF(T2.LE.T1) THEN
        HRI=0.0
      ELSE
        HRI=(SIN(D)*SIN(LP)*(T2-T1)+COS(D)*COS(LP)*(SIN(T2+DDT)                          &
     &       -SIN(T1+DDT)))/(COS(SLOPE1)*DELT1)
      ! In the above the divide by cos slope normalizes illumination to per unit horizontal area
      ENDIF
      ! There is a special case if tpbeg is less than -pi that occurs in polar regions
      ! where a poleward facing slope may be illuminated at night more than the day.
      ! Add this in
      if(tpbeg .lt. -pi)then
	  T1=max(T,-tpbeg+2*pi,-TD)
	  T2=min(T+DELT1,TD)
	  if(t2 .gt. t1)then
	    hri=hri+(SIN(D)*SIN(LP)*(T2-T1)+COS(D)*COS(LP)*(SIN(T2+DDT)                          &
     &       -SIN(T1+DDT)))/(COS(SLOPE1)*DELT1)
	  endif
	  endif
      !for the purposes of calculating albedo we need a cosine of the illumination angle.  This
      !does not have slope correction so back this out again.  This is an average over the
      !time step
      COSZ = HRI*COS(SLOPE1)
      !    write(6,*)hri,coszen

      !RETURN
      END subroutine


      !************************** JULIAN () ****************************
      !            To convert the real date to julian date
      !  The Julian are change to a new version to take the Lean Yean into consideration
      !   in the old version, there are 365 days each year.
      !    FUNCTION JULIAN(MONTH,DAY)
      function julian(yy,mm,dd)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer yy,dd
      dimension mmstrt(1:12)
      data (mmstrt(i),i=1,12)/0,31,59,90,120,151,181,212,243,273,                          &
     &                        304,334/
      jday = mmstrt(mm) + dd
      ileap = yy - int(yy/4)*4
      if(ileap.eq.0.and.mm.ge.3) jday = jday + 1
      julian = jday
      return
      end function

      !***************************** QLIF () *******************************
      !   Computes the incoming longwave radiation using satterlund Formula
      !    Modified 10/13/94 to account for cloudiness.
      !    Emissivity of cloud cover fraction is assumed to be 1.
      !C
      subroutine qlif(qliff,TA,RH,TK,SBC,cf)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EA = SVPW(TA)*RH
      TAK = TA + TK
      EA1 = 1.08*(1.0-EXP(-(EA/100.0)**((TAK)/2016.0)))
      QLIFf =(cf+(1.-cf)*EA1)*SBC*TAK**4
      RETURN
      END subroutine



END MODULE