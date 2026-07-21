      MODULE comm
      ! in here stores temporary storage for computational subroutines that need to
      ! access a whole bunch of numbers are the same time
      ! like the snow model
      REAL*8,DIMENSION(:),POINTER:: SParm
      REAL*8 :: tr,ts,es,cg,z,zo,rho,rhog,lc,ks,de,abg,avo               &
     &   ,anir0,lans,lang,wlf,rd1,fstab,Tref,dNewS,gsurf
      REAL*8 :: ts_old, tave_old, Ts_ave, Tave_ave
      REAL*8 to,tk,sbc,hf,hneu,cw,cs,cp,ra,k,hff,rhoi,rhow               &
     &   ,grav,w1day,pi
      integer :: ComState
      data to /0.0D0/        !  Temperature of freezing (0 C)
      data tk /273.15D0/     !  Temperature to convert C to K (273.15)
      data sbc /2.0747D-7/ !  Stefan boltzman constant (2.0747e-7 KJ/m^2/hr/K)
      data hf /333.5D0/      !  Heat of fusion (333.5 KJ/kg)
      data hneu /2834.0D0/   !  Heat of Vaporization (Ice to Vapor, 2834 KJ/kg)
      data cw /4.18D0/       !  Water Heat Capacity (4.18 KJ/kg/C)
      data cs /2.09D0/       !  Ice heat capacity (2.09 KJ/kg/C)
      data cp /1.005D0/      !  Air Heat Capacity (1.005 KJ/kg/K)
      data ra /287.0D0/      !  Ideal Gas constant for dry air (287 J/kg/K)
      data k  /0.4D0/        !  Von Karmans constant (0.4)
      data hff /3600.0D0/    !  Factor to convert /s into /hr (3600)
      data rhoi /917.0D0/    !  Density of Ice (917 kg/m^3)
      data rhow /1000.0D0/   !  Density of Water (1000 kg/m^3)
      data grav    /9.81D0/     !  Gravitational acceleration (9.81 m/s^2)
	  data w1day /0.261799D0/!  Daily frequency (2pi/24 hr 0.261799 radians/hr)
	  data pi /3.141592654D0/!  Pi
      END MODULE comm
