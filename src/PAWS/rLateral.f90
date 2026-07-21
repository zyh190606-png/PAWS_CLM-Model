
!
      subroutine rLateral(rid)
      ! (E) Compute Source Term and Boundary Conditions for the River Segments-
      ! Source Term: (m3/s/m)
      ! 1.Prcp; (m/s) --*w
      ! 2.Evaporation; (m/day)--*w/86400
      ! 3.Overland Inflow; 4.Tributaries (cms); --/dx
      ! 5.Groundwater Seepage/Contribution (m3/day/m) /86400
      ! 6.Upstream/Downstream BC
      ! Units of all source terms are m2/day (m3/hr/m). so Sprcp = prcp*w
      ! We are not expanding Sprcp when top width increases with increasing flow-
      ! with the thought that the prcp landing in the adjacent cell will-
      ! contribute that water anyway. In stead of doing thing on both sides,
      ! just ignore this step.
      ! Currently, overland exchange is directly applied in F_oc_dw, rather-
      ! than a source term. r.Riv.Qoc = 0 at all times. So this is now mainly-
      ! computing tributary inflows
      ! Input:
      ! tributary flows: r.Riv.Qx
      ! Output: 
      ! r(rid).Riv.S and upstream, downstream boundary conditions
      
      use Vdata
      use polyEval, Only: polyInterp1
      implicit none
      integer rid
      integer, dimension(:), pointer :: rL=>NULL()
      real*8 :: dt, dx, dt2, DTT, nn, frac,dtInflow
      integer i, l, idx, k
      integer, parameter :: n_poly = 2 ! order of accuracy will be n_poly - 1
      integer :: mode=1
      real*8,save :: MX(n_poly,n_poly)
      real*8 :: P(n_poly),ints(n_poly),interval(2),v1,res
      integer, save :: state = 0
      real*8,pointer::tt(:),vv(:)
      
      ! Compute Tributary Contributions (m3/s)
      dt = r(rid)%Riv%dt
      dx = r(rid)%DM%d(1)
      if (associated(r(rid)%Net%Trib)) then
        rL => r(rid)%Net%Trib
        do i = 1, size(rL)
          l = rL(i)
          idx = r(rid)%Net%TribC(i);
          ! Ab = r(rid).Riv.Ab(idx);
      
          dt2 = r(l)%Riv%dt
          DTT = max(dt, dt2)
          r(rid)%Riv%trib(idx) = r(rid)%Riv%trib(idx)+r(l)%Net%DSQ/(DTT/dt2)
          r(l)%Net%DSQ = 0
          !/Ab; % DSQ already upscaled. This is NOT normalized to Ab
        enddo
      endif

      ! Obtain Upstream and Downstream boundary conditions
      !if ~isempty(r(rid).Net.UPS) ! Upstream BC most likely Flux BC
      if(associated(r(rid)%Net%UPS)) then
        if (r(rid)%Net%UPS>size(r) .AND. r(rid)%Net%UPS > 1E6) then
        ! encoding of UPS: 10000+idx.mode (mode =0 for daily forward flow)
        ! currently requires usgs.t to be evenly spaced, meaning: integrated value for the day/interval
          ! if r(rid)%Net%UPS>1 million, then it means this is a recorded inflow condition
          ! r(rid)%Net%UPS - 1e6 gives us the id of the w.tData.usgs station number to retrieve the data from
          ! we put dt of inflow condition in the fraction. if it is 0, then dt is 1.
          ! note the data in the usgs time series must be consecutive
          ! r(rid)%Net%UPST = 2 r(rid)%Net%UPST is not valid
          r(rid)%Riv%UPST = 2D0;
          idx = int(r(rid)%Net%UPS - 1D6 + 1D-14) !r(rid)%Net%UPS is real*8
          vv=>w%tData%usgs(idx)%v
          tt=>w%tData%usgs(idx)%t
          mode = (r(rid)%Net%UPS - 1D6 - idx)*10 
          call polyInterp1(size(vv),tt,vv,(/ r(rid)%Riv%t, r(rid)%Riv%t+r(rid)%Riv%dt /),res,mode,n_poly,MX,state)
          r(rid)%Riv%UPSV = max(res,0D0) ! near 0, sometimes integral interpolation can give negative results
          
          !frac = (r(rid)%Riv%t - w%tData%usgs(idx)%v(k))/(w%tData%usgs(idx)%v(k+1)-w%tData%usgs(idx)%v(k)) ! linearly interpolate
          !r(rid)%Riv%UPSV = (1D0-frac) * w%tData%usgs(idx)%v(k) + frac * w%tData%usgs(idx)%v(k+1)
        else
          r(rid)%Net%UPST = 2
          l = r(rid)%Net%UPS
          nn= r(l)%DM%msize(1) - r(l)%DM%nGho(1)
          r(rid)%Riv%Ac(1) = r(l)%Riv%Ac(int(nn))
          r(rid)%Riv%UPSV = r(l)%Net%DSQ
        endif
      else
        r(rid)%Riv%UPST = 0.0D0
      end if

      ! now support fixed head Downstream BC
      !if ~isempty(r(rid).Net.DS)
      if(associated(r(rid)%Net%DS)) then
        r(rid)%Riv%DST = 4
        l = r(rid)%Net%DS
        idx = r(rid)%Net%DSC
        r(rid)%Riv%DSV(1) = max(r(l)%Riv%h(idx),0.0)
        !r(rid).Riv.DSV(1) = max(r(l).Riv.EstE(idx)-r(l).Riv.E(idx),0);
        r(rid)%Riv%DSV(2) = r(rid)%Riv%Qx(r(rid)%DM%msize(1)-1) ! Q from last time step. added w 01/30/2014
        r(rid)%Riv%w(size(r(rid)%Riv%E)) = r(l)%Riv%w(idx) ! added w 01/30/2014
        r(rid)%Riv%wF(size(r(rid)%Riv%E)) = r(l)%DM%d(1)      ! downstream dx. added w 01/30/2014. This is a lazy-add, wF(nx) means width at the interface, but it doesn't do anything in the code. To be formalized later
        r(rid)%Riv%E(size(r(rid)%Riv%E)) = r(l)%Riv%E(idx) 
        !######### NEED TO THINK ABOUT THE VALIDITY OF ONE EXTRA CELL HERE
        !######### SHOULD JUST MAKE IT AN OVERLAPPING CELL
      else
        r(rid)%Riv%DST = 3      ! Free flowing BC
      end if

      ! Sum up Source Terms
      r(rid)%Riv%S = 0.0D0
      r(rid)%Riv%S = r(rid)%Riv%prcp * r(rid)%Riv%w + (r(rid)%Riv%Qoc + r(rid)%Riv%trib)/dx +  &
     &   r(rid)%Riv%Qgc + r(rid)%Riv%Ql - r(rid)%Riv%evap       ! dA/dt
      ! Currently, only trib matters because Qoc is applied directly in F_oc_dw, prcp is 0 (rained on ground)
      ! Qgc is solved inside riv_cv, Ql and evap are 0s. trib has unit of m2/s because DSQ has m3/s
      
      !4;
      ! if ~isempty(find(r(rid).Riv.S<0))
      !     warning('5');
      ! end
      ! Clear Transitional Fluxes
      !em = zeros(size(r(rid).Riv.h));

      end subroutine rLateral
