#include "fintrf.h"
!
      subroutine adhoc_Set()
      !% (E) These codes are ad-hoc settings that will be implemented in formal codes
      !% This include initialization of some variables and setting some parameter
      !% values

      !global Env
      use Vdata
      implicit none
      real*8, dimension(size(w%tData%dc,1),size(w%tData%dc,2)), target :: dc
      real*8, dimension(w%m(1),w%m(2)) :: botD,sum_temp,mask2D_dbl
      integer, dimension(size(w%CR)) :: CR

      integer, dimension(:), pointer :: rL=>NULL()
      real*8,dimension(w%m(1),w%m(2),w%m(3))::temp,mask3D_dbl

      integer,dimension(w%m(1),w%m(2))::SMAP
      integer,dimension(w%m(1)*w%m(2))::SMAP_1D

      logical*1,dimension(w%m(1),w%m(2),w%m(3))::mask3D

      integer :: nr, nc, nz, i, j, k
      real*8 :: fc, THES, THED, THE2

      g%VDZ%th%dfc = w%tData%dc
      !% snow sitev
      g%VDZ%th%Ssitev(:,:,2) = 1.0D0
      g%VDZ%th%Ssitev(:,:,3) = (101.3D0-0.01152D0*g%Topo%E+0.544D-6*g%Topo%E**2)*1000.0D0
      g%VDZ%th%Ssitev(:,:,4) = 3.6D0
      g%VDZ%th%Ssitev(:,:,5) = 0.1D0
      g%VDZ%th%Ssitev(:,:,6) = 0.0D0    !% slope, not used, but included for potential use
      nr = size(g%Veg%RPT,1)
      nc = size(g%Veg%RPT,2)
      nz = size(g%Veg%RPT,3)

      ! % comparison purposes::
      g%VDZ%th%Ssitev(:,:,3) = 98312.53D0
      g%VDZ%th%Sstatev(:,:,1) = 11000.0D0
      do i = 1, nr
        do j = 1, nc
          fc = 0.0D0
          do k = 1, nz
            if ((g%Veg%RPT(i,j,k)==3) .or. (g%Veg%RPT(i,j,k)==4)) then
              fc = fc + g%Veg%Frac(i,j,k)
            end if
            g%VDZ%th%Ssitev(i,j,6) = fc
          end do
        end do
      end do
      
      g%VDZ%th%tsa = -9999.0D0
      g%VDZ%th%tsa(:,:,24) = 0.0D0
      g%VDZ%th%tss = -9999.0D0
      g%VDZ%th%tss(:,:,24) = 0.0D0

      !%%
      g%VDZ%VParm = [0.08D0, 0.3D0, 1.0D0, 0.0D0, 1.46D0, 3.0D0, 0.0D0, 5D-2, 0.05D0, 5D-3, 1.0D0]
      !  % 1gamma, 2ST multiplier, 3frozen soil K coefficient--(fraction of KS at
      !  % T=-5C),4t0 (for freezing start), 5ho multiplier, 6h0 power,
      !  % 7KSdecay,8decayBottom,9anisotropic ratio for the botK between GW(1) and
      !  % GW(2),10hf for F_contribute, 11, multiplier for depressionS K
      !  % SHOULD USE Ko!!!!!!!!!
      g%VDZ%KS = g%VDZ%Ko
      call display( 'Using Ko instead of KS')

      g%OVN%FParm = [0.0D0, 15.0D0, 1.0D0]    !% [river bed rise, use 0.3, depS K divisor, riv K multiplier]
      nr = size(g%VDZ%h,1)
      nc = size(g%VDZ%h,2)
      nz = size(g%VDZ%h,3)
      !%g.GW(2).K(:)=0; % lower bed closed
      g%OVN%hf = 0.005D0
      g%OVN%OVF = 0.05D0   !% overland flow fraction
      g%Topo%S0p = 0.005D0
      g%OVN%dist = 506.0D0
      g%OVN%ho = 0.0D0
      g%OVN%Rf = 0.0D0
      g%GW(1)%MAPP(:,:) = w%DM%mask
      g%GW(2)%MAPP = w%DM%mask
      g%OVN%SATH = 0.0D0
      !g%VDZ%N = g%VDZ%N * 1.05
      !%g.VDZ.ALPHA = g.VDZ.ALPHA*(3/4);
      g%VDZ%KSSAVE = g%VDZ%KS + 0.0D0
      g%VDZ%WTL = g%VDZ%NZC + 0.0D0
      g%GW(1)%topBound = g%Topo%E + 15.0D0

      if (len_trim(w%msg)==0) then
        !%g.GW(1).K = g.GW(1).K;
        !%g.VDZ.KS = g.VDZ.KS*3;
        w%msg = 'adjusted'      !% flag for already modified data
        CR = w%CR
        !%g.VDZ.LAMBDA = g.VDZ.LAMBDA-2;
        !%g.VDZ.KS = g.VDZ.KS *1.5;
        !g%VDZ%ALPHA = g.VDZ.ALPHA ???
        g%Veg%hI = 0.0D0
      end if
      g%GW(1)%ST = 0.01D0
      g%GW(2)%ST = 5.9160D-7   !% find from previous optimization
      !%g.GW(1).Cond.K(:) = 0; % Block seepage to below!!!!!
      !%g.GW(1).Cond(1).K = g.GW(2).K(g.GW(1).Cond(1).indices);
      !% aquitard
      do j = 1, nc
        do i = 1, nr
          botD(i,j) = 0.0D0 + 2.0D0     ! 2 meters of aquitard
          if(botD(i,j) > g%GW(2)%DZ(i,j)/2.0D0) then
            botD(i,j) = g%GW(2)%DZ(i,j)/2.0D0
          end if
        end do
      end do
      do i = 1, 2
        g%GW(i)%E = 0.5D0*(g%GW(i)%EB+g%GW(i)%ES)
      end do
      g%GW(1)%botD = botD + 0.0D0
      g%GW(1)%botK = 0.0D0 + 8D-5     !% conductivity
      g%VDZ%GWL = g%GW(1)%h + 0.0D0

      !  %% fill missing temperature
      !  %{
      !  t = w.Stations(1).datenums;
      !  v = w.Stations(1).s0;
      !  v2= w.Stations(1).s1;
      !  v(t<datenum2(19990603))=v2(t<datenum2(19990603));
      !  loc = (isnan(v) | v<-90);
      !  v(loc) = interp1(t(~loc),v(~loc),t(loc));
      !  w.Stations(1).s0 = v;
      !  %}
      !%% precompute water deficit of unsat zone (to saturation)
      !% save a copy of initial groundwater head at year start
      w%Rec%hG = g%GW(1)%h + 0.0D0
      rL => w%CR
      nr = size(rL)
      w%Rec%acc%Qoc = 0.0D0
      w%Rec%acc%Qout = 0.0D0
      w%Rec%acc%trib = 0.0D0
      w%Rec%Mem = 0.0D0
      w%Rec%MemR = 0.0D0
      w%Rec%MemQgc = 0.0D0

      temp = g%VDZ%THES * g%VDZ%DZ
      sum_temp = sum(temp,3)
      THES = sum(sum_temp)/size(sum_temp)
      do i = 1, nz
        mask3D(:,:,i) = (.not. w%DM%mask)
      end do
      mask3D_dbl = 0.0D0
      where(mask3D)
        mask3D_dbl = 1.0D0
      end where
      !mask = repmat(~w.DM.mask,[1 1 nz]);
      temp = (g%VDZ%THES - g%VDZ%THE) * g%VDZ%DZ * mask3D_dbl
      THED = sum(temp(:,:,2))/size(temp(:,:,2))
      sum_temp = sum(temp,3)
      THE2 = sum(sum_temp)/size(sum_temp)
      w%Stats%THES     = THES
      w%Stats%THE_FULL = THE2
      w%Stats%THED     = THED

      !%%
      !where(g%VDZ%LAMBDA > 0)
      !  g%VDZ%LAMBDA = -g%VDZ%LAMBDA
      !end where
      ! using where occurs stupid stack overflow error
      do k = 1, size(g%VDZ%LAMBDA,3)
        do j = 1, size(g%VDZ%LAMBDA,2)
          do i = 1, size(g%VDZ%LAMBDA,1)
            if(g%VDZ%LAMBDA(i,j,k) > 0) then
              g%VDZ%LAMBDA(i,j,k) = -g%VDZ%LAMBDA(i,j,k)
            endif
          enddo
        enddo
      enddo
      !g.VDZ.LAMBDA(g.VDZ.LAMBDA>0)=-g.VDZ.LAMBDA(g.VDZ.LAMBDA>0);
      where(g%GW(2)%h < g%GW(1)%EB)
        g%GW(2)%h = g%GW(1)%h
      end where
      !loc = find(g.GW(2).h<g.GW(1).EB); g.GW(2).h(loc)=g.GW(1).h(loc);

      !loc = find(g.GW(2).h<g.GW(1).EB); g.GW(2).h(loc)=g.GW(1).h(loc);
      !  %% adjust some parameters
      !  %{
      !  g.GW(1).K=g.GW(1).K/5;
      !  for i=1:length(w.CR)
      !      c = w.CR(i);
      !      r(c).Riv.K = r(c).Riv.K/3;
      !      r(c).Riv.TA= sum(r(c).Riv.w(2:end-1).*r(c).Riv.DX);
      !  end

      !  CR = w(wid).CR;
      !  for i = 1:length(CR)vdzInitValues(gid)
      !      rid = CR(i);
      !      r(rid).Riv.Mann(:)=0.02;
      !      r(rid).Riv.K(:) = r(rid).Riv.K(:)/50; % Originally was too big!!
      !  end
      !  g.GW(1).ST = g.GW(1).ST*8;

      !  %}

      !g.GW(1).msg = struct('TF',0);
      g%GW(1)%msg%TF = 0.0D0
      !%g.OVN.Ds(:) = 2e-3; % 2mm depression storage
      !%g.OVN.Ds(:) = 0;

      !g.Wea.month =struct('prcp',[],'temp',[],'tmx',[],'tmn',[],'rad',[],'wnd',[],'hmd',[]);
      g%Wea%month%prcp = 0.0D0
      g%Wea%month%temp = 0.0D0
      g%Wea%month%tmx = 0.0D0
      g%Wea%month%tmn = 0.0D0
      g%Wea%month%rad = 0.0D0
      g%Wea%month%wnd = 0.0D0
      g%Wea%month%hmd = 0.0D0


      !SMAP = int32(zeros(size(g(gid).OVN.h)));
      SMAP = 0_4
      SMAP_1D = oneD_Ptr(size(SMAP),dble(SMAP))
      do i = 1, size(w%Stations)
        if(associated(w%Stations(i)%gIndices)) then
          SMAP_1D(w%Stations(i)%gIndices) = i
        endif
      end do
      SMAP = reshape(SMAP_1D,[size(SMAP,1),size(SMAP,2)])
      g%Veg%SMAP = SMAP
      g%Veg%SNOW = 0.0D0
      g%Veg%NewS = 0.0D0

      !g(1).Riv.K = g(1).Riv.Len;
      !% for i=1:length(g(1).Riv.K),g(1).Riv.K{i}(:)= 0;end
      do i = 1, size(g%Riv%K)
        g%Riv%K(i)%ptr(:)= 30.0D0
      end do
      !%Ex = E + E(i+1)
      !%% printout headers

      g%Veg%Gprcp = 0.0D0
      g%Veg%EvapG = 0.0D0

      !  %% display menu
      !  %{
      !  global Prob
      !  menu=struct('title','r(26).Riv.Eta','input',0,'mainobj','r','gid',26,'subobj','Riv','obj','Eta','subnum',1);
      !  riv=[25 26 27 30 31 32];
      !  for i=1:length(riv)
      !      c=riv(i);
      !      menu((i-1)*2+1)=struct('title',['r(',num2str(c),').Riv.E'],'input',0,'mainobj','r','gid',c,'subobj','Riv','obj','E','subnum',1);
      !      menu(i*2)=struct('title',['r(',num2str(c),').Riv.Eta'],'input',0,'mainobj','r','gid',c,'subobj','Riv','obj','Eta','subnum',1);
      !  end
      !  menu(length(menu)+1)=struct('title','w(1).Rec(1).tseries','input',0,'mainobj','w','gid',1,'subobj','Rec','obj','tseries','subnum',1);
      !  Prob.display = menu;
      !  %}
      !%% w stats pre-computed
      w%Stats%TA = g%DM%d(2)*g%DM%d(1)*g%DM%msize(2)*g%DM%msize(1)
      mask2D_dbl = 0.0D0
      where(w%DM%mask)
        mask2D_dbl = 1.0D0
      end where
      w%Stats%TAV= sum(mask2D_dbl)*g%DM%d(2)*g%DM%d(1)
      !w.tData.h1save = g.GW(1).h+0;
      !w.tData.h2save = g.GW(2).h+0;
      w%tData%h1save = g%GW(1)%h + 0.0D0
      w%tData%h2save = g%GW(2)%h + 0.0D0
      !%############### later add in lulc stats

      end subroutine adhoc_Set
