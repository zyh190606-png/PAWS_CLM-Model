
!

      subroutine PenmanMonteith(nt,Rad,Hmd,e0,Wnd,T,Elev,tau,ref,rET,rETs,Rn,Hb,Hli,rss)
      !% (E) This function computes the reference ET using Penman Monteith method
      !% Steps of calculating ET:
      !% (1), reference ET for water, soil and vegetations (reference type)
      !% (2), spread solar radiation using Beer-Lambert law and LAI
      !% (3), compute total ET demand for each layer by considering crop coeff
      !% (4), compute actual ET by considering water stress and depth distribution
      !% (5), add up ET at each layer to become ET sink term for VDZ model
      !% inside this function, Rad, Hb, Hli, Rn are MJ/m2/day
      !% rET, rETs are mm/day
      !% RET is in mm/day
      !% Input:
      !% tau: air transmittance, ref: code for different reference basis
      !% (others self explanatory)
      !% Output:
      !% rET: reference ET; rETs: reference ET for bare soil, Rn: net radiation-
      !% Hb: outgoing long wave radiation, Hli: incoming long wave radiation
      implicit none
      integer nt
      real*8,dimension(nt)::Rad,Hmd,e0,Wnd,T,rET,rETs,Rn,Hb,Hli,rcc
      real*8::Elev,tau,rss(nt,2)
      integer::ref
      real*8::Cp,sigma,P,alb,Hmx,ct,cf,fcld,rc
      real*8,dimension(nt)::TK,rho,lambda,gamma,Delta,e,DE,Rs,emm,ccc,EA,EA1,c,a1,a2,ra
      real*8 b,m,pp,e_oc,den
      parameter(pp = 2, m = 8, b = 0.440D0, e_oc = 0.976D0)
      real*8,dimension(nt)::e_cs,eff

      TK = T + 273.15D0
      Cp = 1.013D-3
      sigma = 4.903D-9
      rho = 1710.0D0 - 6.85D0*T         !this is not air density, this is a combined term
      den = 1.2250D0
      P = 101.3D0-0.01152D0*Elev+0.544D-6*Elev**2       !% Atmospheric Pressure in kPa
      lambda = 2.501D0-2.361D-3*T       !% Specific Heat of evaporation MJ/kg
      gamma = (Cp/0.622D0)*(P/lambda)   !% psychrometric constant, kPa/C
      !%e0 = exp((16.78.*T-116.9)./(T+237.3)); % saturation vapor pressure in kPa
      !%Delta = 4098*(e0./TK.^2); % Slope of saturation Curve THIS IS A BUG!!!!
      Delta = 4098.0D0*(e0/(T+237.3D0)**2)
      e = Hmd*e0
      DE = e0 - e
      !alb = zeros(size(Elev)); alb(:) = 0.23;
      alb = 0.23D0
      !%alb(sno> 0.5)= 0.8;
      Rs = (1.0D0-alb)*Rad      !%SNOW COVER CHANGE ####################################
      !% long-wave radiation.
      emm = -(0.34D0 - 0.139D0*dsqrt(e))
      Hmx = maxval(Rad)
      !%fcld = Rs*(0.9/Hmx)+0.1;
      ct = tau/0.7D0
      cf = 1.0D0 - ct       !% cloud fraction. 0.7 is the max in Spokas setting
      fcld = 0.9D0*ct+0.1D0     !%fcld = tau/0.7 is the cloud cover factor
      !%fcld = 0.9*(tau)+0.1; % temporarily use this!!!!!!!!!!!!!
      ccc = sigma*TK**4

      Hb = fcld*emm*ccc     !% UNIT--MJ/m2/d
      Rn = Rs + Hb          !%Rn(Rn<0)=0; % THIS COULD BE A BUG!!!!!!!!!!!!!!!!!!!!!!
      !if any(Rn<0)
      !  4;
      !end
      EA = e0*Hmd
      !%EA1 = 1.08*(1.0-exp(-(EA/100.0).^(TK/2016.0))); % This is a bug, unit of
      !%EA is kPa, but it only impacts snowmelt
      !CP: the equation below is found to generate too high emissivity values, and may have caused plant temperature to be very high
      !metabolism breaks down, so less latent heat, then even higher temperature, therefore creating a positive feedback loop which results in too little ET
      !and too high soil temperature. 
      !EA1 = 1.08D0*(1.0D0-dexp(-(EA*10.0D0)**(TK/2016.0D0)))
      !Hli= (cf+ct*EA1)*ccc
      ! now using: Hli = (e_cs*ct^p + e_oc*(1-ct^p))*sigma*T**4 where e_oc is emmissivity under overcast, e_cs is emmissivity under clearsky
      ! Konzelmann et al. (1994): e_cs = 0.23+b*(e/TK)**(1/m). (e: Pa, in our calculations, it was in kPa)
      ! b = 0.484, m=8 (Konzelmann, via Sedlar, Cryosphere, 2008, On the use of incoming longwave radiation parameterizations in a glacier environment)
      ! e_oc = 0.952, p = 4
      
      e_cs = 0.23D0 + b*(e*1000D0/TK)**(1D0/m)
      eff  = (e_cs*ct**pp + e_oc*(1D0 - ct**pp))
      Hli = eff*ccc

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !% ref = 2;


      if (ref == 0) then    !% FAO reference grass, height 0.12, surface resistance 90 s/m, albedo=0.23
        c = (Delta+gamma*(1.0D0+0.34D0*Wnd))
        a1 = 0.408D0*Delta/c
        a2 = gamma*DE*Wnd/TK*900.0D0/c
        rET= a1*Rn + a2
        ra = 208.0D0/Wnd        !% ra of the grass reference
        !rETs = (Rn*Delta + rho*Cp*DE/ra)/(lambda*(Delta+gamma))     !% Soil Maximum ET. We can discuss this later
        ! above is bug code
        rETs = (Rn*Delta + rho*gamma*DE/ra)/(lambda*(Delta+gamma))
        !% 0.408 may be substituted by 1/lambda
        !% assume a 0 of Ground heat flux
        !% http://www.fao.org/docrep/X0490E/x0490e08.htm#chapter%204%20%20%20determination%20of%20eto
      elseif (ref == 1) then    !% alfalfa
        ra = 114.0D0/Wnd
        rc = 49.0D0
        !%CT = 1710 - 6.85*T;
        rET = (Delta*Rn/lambda+gamma*CT*DE/(ra*lambda))/(Delta+gamma*(1.0D0+rc/ra))
        !%E0 = (Delta.*Rn+ (1710-6.85*T).*DE./ra)./(Delta+gamma.*(1+rc/ra))
        !%WRONG
      elseif (ref == 2) then    !% add Priestley-Taylor as RETs
        rETs = (1.28D0*Delta/(Delta+gamma)*Rn)/lambda
        ra = 114.0D0/Wnd
        rc = 49.0D0
        !%CT = 1710 - 6.85*T;
        rET = (Delta*Rn+gamma*rho*DE/(ra))/(lambda*(Delta+gamma*(1.0D0+rc/ra)))
      elseif (ref == 3) then
      ! test ra and rs input
        ra = rss(:,1)
        rcc = rss(:,2)
        rET = (Rn*Delta + rho*gamma*DE/ra)/(lambda*(Delta+gamma*(1.0D0+rcc/ra)))
      endif

      !rET(rET<0)=0;
      !rETs(rETs<0)=0;
      where(rET<0)
        rET = 0.0D0
      endwhere
      where(rETs<0)
        rETs = 0.0D0
      endwhere


      end subroutine PenmanMonteith
