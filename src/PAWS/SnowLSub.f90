
	subroutine SnowLSub(DT,nStep,INPT,SITEV,STATEV    &
     &  ,PARAM,iflag,nstepday,                        &
     &  cump,cume,cummr,outv,tsbackup,tavebackup,     &
     &  ndepletionpoints,dfc)    ! 4/2/05  DGT added modelelement for debugging
      use comm
      USE SNOWUEB_MOD, Only : SNOWUEB2
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(niv=7)
      INTEGER YEAR,DAY,iflag(*),stepInDay  ! Pass a flag to control surface temperature modeling
      real*8 IO,LAT
      real*8 INPT(niv,1)
      real*8 sitev(*)
      real*8 outv(*)
      real*8 statev(*)
      real*8 param(*)
      real*8 dtbar(12)
	real*8 mtime(5)      ! YJS pass to reflect the change of Snow
	real*8 dfc(ndepletionpoints,2)
	real*8 hour
	real*8 tsbackup(*)
	real*8 tavebackup(*)

      !  Constant data set
      ! CP set in the module
      !  End of constant declaration
      IF (ComState .eq. 1) THEN
       ! CP: This should have already been intialized outside
      ELSE
      tr = Param(1)     !  Temperature above which all is rain (3 C)
      ts = Param(2)     !  Temperature below which all is snow (-1 C)
      es = Param(3)     !  emmissivity of snow (nominally 0.99)
      cg  = Param(4)    !  Ground heat capacity (nominally 2.09 KJ/kg/C)
      z = Param(5)      !  Nominal meas. height for air temp. and humidity (2m)
      zo = Param(6)     !  Surface aerodynamic roughness (m)
      rho = Param(7)    !  Snow Density (Nominally 450 kg/m^3)
      rhog = Param(8)   !  Soil Density (nominally 1700 kg/m^3)
      lc = Param(9)     !  Liquid holding capacity of snow (0.05)
      ks = Param(10)    !  Snow Saturated hydraulic conductivity (20 !160 m/hr)
      de = Param(11)    !  Thermally active depth of soil (0.1 m)
      abg = Param(12)   !  Bare ground albedo  (0.25)
      avo = Param(13)   !  Visual new snow albedo (0.95)
      anir0 = Param(14) !  NIR new snow albedo (0.65)
	lans =  Param(15) !  the thermal conductivity of fresh (dry) snow (0.0576 kJ/m/k/hr)
	lang =  Param(16) !  the thermal conductivity of soil (9.68 kJ/m/k/hr)
	wlf =  Param(17)  !  Low frequency fluctuation in deep snow/soil layer (1/4 w1 = 0.0654 radian/hr)
	rd1 =  Param(18)  !  Amplitude correction coefficient of heat conduction (1)
      fstab = Param(19) !  Stability correction control parameter 0 = no corrections, 1 = full corrections
	Tref = Param(20)  !  Reference temperature of soil layer in ground heat calculation input
	dNewS = Param(21) !  The threshold depth of for new snow (0.001 m)
	gsurf = Param(22) !  The fraction of surface melt that runs off (e.g. from a glacier)
	ComState = 1
	ENDIF

      !  Unravel state variables
      afrac=statev(7)
	baw=statev(2)
      Wmax = statev(6)
	meltflag=statev(8)
	Wtgt = statev(9)     !target snow water equivalent for snowfall hysteresis loops during melt
	Wtgtmax=statev(10)
	Aftgt=statev(11)
	irad=iflag(1)
	slope=sitev(6)
	!azi=sitev(7) ! these variables are no longer used as radiation are estimated outside
	!lat=sitev(8)
      ! debugging code to stop at a specific date
      !   if((year .eq. 1960) .and. (month .eq. 10) .and.
      !  + (day .eq. 27) .and. ! (hour .gt. 8) .and.
      ! + (modelelement .eq. 110))then
      !c	   year=year
      !c	endif

      do 2 iit=1,nStep    !  time loop
	IF(AFRAC.gt.0 .and. BAW .gt. 0) then
      !See discussion in "Depletion curve logic.doc"  ! DGT 7/25/05.  Revised to avoid spikes this was causing
      !  The idea here is to use a weighted average of the snow water equivalent depths over the old and
      !  new snow covered area as the representative snow water equivalent to give the point model while
      !  on an excursion.
		If(meltflag .eq.0 .and. Aftgt > 0                                       &
   &  	.and. Afrac > Aftgt .and. BAW > Wtgt)then  !  Here we are on an excursion
      !c			SCAW=Wtgt/Aftgt+BAW-Wtgt   ! DGT 7/25/05.  Initial approach superceded by the following
			AfNew=Afrac-Aftgt
			WoldSCA=Wtgt/Aftgt
			Wnew=BAW-Wtgt
			WnewSCA=Wnew/AfNew
			SCAW=(WnewSCA*Wnew+WoldSCA*Wtgt)/BAW
		else
		    SCAW = BAW/Afrac
		endif
      else  ! initialize for no snow
        SCAW = 0D0
	  Afrac=0D0
	  BAW=0D0
	  Wmax=0.000001D0   !  DGT 7/25/05 to avoid NaN when first used
	  Wtgtmax=0D0
	  meltflag=1D0
      endif
      statev(2)=SCAW
      ! DGT 7/26/05.  Statev(1) in the lumped model stores the average temperature except in the case
      ! when liquid water is present in which case it stores the liquid fraction.
      ! The assumption is that these quantities remain constant during depletion curve adjustments to W
      ! and are used to reinstate a physically realistic U.
      !cg =param(4)    !  Ground heat capacity (nominally 2.09 KJ/kg/C)
      !rhog=param(8)   !  Soil Density (nominally 1700 kg/m^3)
      !de=param(11)    !  Thermally active depth of soil (0.1 m)
	if(statev(1)< 0.)then
	    !  Here energy content state variable was temperature and needs to be converted to energy
	    if (statev(1)<-100) then
	    ELSE
		    statev(1)=statev(1)*(rhog*de*cg+rhow*SCAW*cs)
		  ENDIF
	else
		if(SCAW > 0.)then
	       !  Here energy content state variable was liquid fraction
			statev(1)=statev(1)*rhow*SCAW*hf
		else
		    if (statev(1)>1) then
			! CP added this IF since initial input was 11000
			else
			!  Here energy content variable was temperature.  SCAW is 0 so energy is only of soil
			   statev(1)=statev(1)*rhog*de*cg
			endif
		endif
	endif


      ! *****************This is start of point model
      !c
      ! From the  input list
      RH=inpt(4,1)
      QSIOBS=inpt(5,1)
	QNETOB=inpt(6,1)
	trange=inpt(7,1)
	TA=inpt(1,1)

	! the big chunk deleted by Chaopeng was computing short and long wave radiation
	! we will just pass in


      ! cumulative variables used for single time step over snow covered area
	cump=0D0
	cume=0D0
	cummr=0D0
	mtime = 0D0

       CALL SNOWUEB2(DT,1,INPT,SITEV,STATEV,                                                &
     &  tsbackup,tavebackup,nstepday,PARAM,iflag,                                &
     &  cump,cume,cummr,outv,mtime)   !Add a pass varible to snow mtime

      !c ************************ End of point model

      !   DeltaW = Change in water equivalent
      DeltaW=statev(2)-SCAW
      if(outv(2) .gt. 0. .and. DeltaW .gt. 0.)then   ! Here an increase in snow water equivalent due to snowfall
		If(meltflag .eq. 1)then
		  Meltflag=0   !  Start new excursion
		  Wtgt=BAW
		  Aftgt=Afrac
		  Wtgtmax=Wtgt  ! Establish new upper limit of new excursion
		endif
		BAW=BAW+DeltaW
		Wtgtmax = max(BAW,Wtgtmax)
		If(BAW .gt. Wmax)then
		 Wmax=BAW
		 Meltflag=1   ! end up on main curve if Wmax increases
		elseif(BAW > 3.D0* Wtgt)then   !  DGT 7/25/05   If new snow dominates old snow we
		!  reset the process to be on the depletion curve.
		!  Here domination is arbitrarily defined as greater than 3 times.
		    Meltflag=1
			Wmax=BAW
		endif
		Afracnew=1
		BaSwin=cummr    ! basin average surface water input
		Basub=cume   ! basin average sublimation
	Else if( DeltaW .gt. 0D0) then  ! here a rare case when W increases but there is no snowfall.
	!   This may be due to rainfall absorbing and freezing into the snow or condensation.
	!    This is assumed to only be effective over snow covered area fraction
		BAW=BAW+DeltaW*Afrac
		Wtgtmax = max(BAW,Wtgtmax)
		If(BAW .gt. Wmax)then
		 Wmax=BAW
		 Meltflag=1   ! end up on main curve if Wmax increases
		Endif
      ! DGT 7/25/05  Logically the dominates check above should be done here too, but it is not
      ! because if Wmax is reset Afrac should go to 1 to put the process at the top of a depletion curve
      ! and we are not increasing Afrac because there is no snowfall.  The implication is that BAW may
      ! go greater than 3 * Wtgt and the depletion curve only get reset next time snow falls.
      ! The 3 is arbitrary so this does not seem like a serious shortcoming.
		Afracnew=Afrac   ! (emphasize does not change)
      !       One could increase Wmax here to keep on main depletion curve.  However I decided to handle this issue below with the Min statements so that Af never increases when there is melt.  Excursions below the depletion curve will therefore return back to the same Wa before further reductions in Af occur.
          BaSwin=cummr*afrac+cump*(1-afrac)  !  All precipitation is surface water input over snow free area
		Basub=cume*Afrac  ! Sublimation is only counted over snow covered area.  Calling model has to handle ET from soil
	Else   ! Here (DeltaW < 0) so there is melt
	    If(meltflag .eq. 0D0 .and. (BAW-DeltaW*Afrac) < Wtgt)then
      ! Here reduction passes target so do in two steps to avoid overshooting and large steps towards end of season
			DW1=-(BAW-Wtgt)/Afrac
			DW2=-(DeltaW-DW1)
      !           now logically in two steps
      !              BAW=BAW+DW1*Afrac+Aftgt*DW2
      !           which simplifies to
			BAW=Wtgt+Aftgt*DW2
		  else
	        BAW=BAW+DeltaW*Afrac
	    endif
		If(BAW .lt. 0D0)then
      !    There is a mass balance issue here because BAW is insufficient to supply M and E that
      !    result in the reduction DeltaW.
             !r=1-(BAW/Afrac*DeltaW)   ! This assumes the shortfall will be evenly spread between M and E
             ! the above original UEB statements create mass balance errors
           r = 1 + BAW/(Afrac*(cume+cummr))
		   BAW=0.D0
			cummr=cummr*r
			cume=cume*r
      !     DGT 7/26/05  Since snow is disappearing here we should make sure energy behaves properly
      !     Set energy content to 0.  There was snow so any energy was in liquid fraction not temperature greater
      !     than 0 and that is now gone.
             if(statev(1) >= 0.)then
	           statev(1)=0.
      !c			else
      !c				statev(1)=statev(1)    ! redundant but debugging to see if this ever occurs
			endif
	    endif
		BaSwin=cummr*Afrac+cump*(1-Afrac)
		BaSub=cume*Afrac
		If(BAW .lt. Wtgt)meltflag=1  ! Back on main melt curve
		If(meltflag .eq. 1)then
			Afracnew=Min(Ablat(BAW/Wmax,ndepletionpoints,dfc),Afrac)   ! Regular depletion curve.  Make sure that Af does not increase in this step for special case where W may have increased to be below depletion curve by absorbing of precip or snow.  Assume Af remains the same until W returns to original Af
		Else
			Wref= Wmax-(Wtgtmax-BAW)/(Wtgtmax-Wtgt)*(Wmax-Wtgt)  ! See depletion curve logic.doc for derivation and picture
			Afracnew=Min(Ablat(Wref/Wmax,ndepletionpoints,dfc),Afrac)
		Endif
	endif

      !DGT 7/25/05.  Save in statev(1) the energy state variable that persists through depletion curve
      ! BAW and SCAW adjustments
      !
	if(statev(1) .le. 0.)then
	!  Here convert energy to temperature - using statev(2) output from the point model prior to
	!  depletion curve adjustments
	   statev(1)=statev(1)/(rhog*de*cg+rhow*statev(2)*cs)
	else
		if(BAW > 0.)then
	!  Here snow water equivalent is positive and energy content positive so
	!  energy content variable is the liquid fraction
	! Note - this is a strictly greater than.  A tolerance could be considered
	! Note also that the test is on BAW because depletion curve adjustments may result in this
	! being 0 when statev(2) is not so then the temperature inter
			if(statev(2) <= 0.)then
				if(BAW > 1e-6)then
	!  Due to arithmetic precision in the division and multiplication by afrac
	!   BAW > 0 can occur when statev(2) = 0.  Only consider this an error if > 1e-6
					write(6,*)"Error, BAW > 0, statev(2) <= 0."
					BAW=0.
					statev(1)=0.
				else
					BAW=0.
	!  Here energy content is positive and snow water equivalent is 0 so
	!  energy content state variable is the average soil temperature
					statev(1)=statev(1)/(rhog*de*cg)
				endif
			else
				statev(1)=statev(1)/(rhow*statev(2)*hf)
			endif
		else
			if(statev(2) > 0.)then
	!  Here BAW was 0 but statev(2) is positive.  This can occur due to depletion curve adjustments
	!  Treat this as a case where snow water equivalent is 0 but disregard the previous energy content
	!  because it is greater than 0 and reflects liquid water fraction which has been removed.
				statev(1)=0.
				write(6,*)'Liquid fraction energy content set to 0'
				write(6,*)'because depletion curve set W to 0'
			else
	!  Here energy content is positive and snow water equivalent is 0 so
	!  energy content state variable is the average soil temperature
				statev(1)=statev(1)/(rhog*de*cg)
			endif
		endif
	endif


      !Output variables
	outv(15)=BAW
	outv(16)=Afrac
	if(Afracnew .eq. 1)outv(16)=Afracnew  ! output old Afrac unless in increased because old Afrac was used in mass balances
	Afrac=Afracnew
	outv(17)=meltflag
	outv(18)=Wtgt
	outv(19)=Wtgtmax
	outv(20)=Wmax
	outv(21)=Aftgt
	outv(22)=BaSwin
	outv(23)=BaSub
	!if(iflag(2) .eq. 1)WRITE(iflag(3),*)(outv(i),i=15,23)
 2	continue   ! end of time step loop

      iflag(1)=irad  ! reinstate iflag(1) for next call
      statev(2) = BAW
	statev(6)=Wmax
	statev(7)=Afrac
	statev(8)=meltflag
	statev(9)=Wtgt
	statev(10)=Wtgtmax
	statev(11)=Aftgt
	return
      END


       function Ablat(DimBaw, ndfc, dfc)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       real*8 Ablat, DimBaw
	 real*8 dfc(ndfc,2)
       if(DimBaw .ge. 1.0) then
	   Afrac =1
	 else
	   if(DimBaw .le. 0.0) then
	     Afrac=0.0
	   else
           do 500 i=1 , ndfc
	       if(i .eq.1) then
		     Afrac=dfc(1,2)   ! DGT 7/25/05 to avoid Afrac ever being interpolated as 0 and avoid divide by 0 errors
	       else
               if(dfc(i,1).ge.DimBaw) then
	           Afrac = dfc(i-1,2)+(dfc(i,2)-dfc(i-1,2))*(DimBaw-                          &
     &		   dfc(i-1,1))/(dfc(i,1)-dfc(i-1,1))   ! assumes points not repeated in the first column of dfc
	           goto 1000
	         endif
	       endif
500        continue
         endif
	 endif

1000   if(Afrac .ge. 1.0) Afrac = 1.0

       Ablat = Afrac
       return
       end
