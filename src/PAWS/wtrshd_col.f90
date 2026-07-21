#include <preproc.h> 
!

      subroutine wtrshd_col(T,JD,TT,sig)
      !% (E) The daily marching script that runs watershed model for a day.
      !% Basic time step is hourly (or the macro watershed time step w.dt).-
      !% Surface flow and river flow may run faster.
      !% Input:
      !% wid: watershed id, most likely 1, T: final time, most likely current-
      !% time + 1 (day), JD: julian day (of the year)
      !% Output:
      !% TT: CPU time collected for different components:
      !% [1vegation+vadoseZone, 2overlandFlow, 3riverNetwork, 4groundwater, 5output]
      !% sig: whether the model exited the day normally (0) or with error (1)-
      !% updated all relevant fields in g and r
      use Vdata
      use flow
      use surface_proc_mod
      use budgeter_mod, Only: budgeter
      implicit none
      !global g r w
      real*8 :: T, TT(5)
      integer :: JD, sig
      real*8, parameter :: PREC = 1D-6
      real*8 :: gid, nFT, diff, RDT
      integer, dimension(size(w%CR)) :: rL
      real*8, dimension(size(w%CR)) :: DTS
      integer, dimension(size(w%ER)) :: eL
      real*8, dimension(size(g%OVN%h,1),size(g%OVN%h,2)):: hf
      real*8, dimension(size(g%OVN%hc,1),size(g%OVN%hc,2)):: hc
      real*8, dimension(size(g%Veg%hI,1),size(g%Veg%hI,2)):: hI, h      
      integer :: i, j, SCODE, Hour, kk, m2(3), c
      character*32 :: date_char, time_char
      character*32 :: rMB, rSig, wSig, gSig, pondSig, colSig
      real*8 :: time_start, time_end, dtNew, dtOld
      real*8, dimension(:,:), pointer :: mask=>NULL()
      logical :: subp, ovnp, rivp
      integer :: dbflag, scalar(2),k,filter(size(g%OVN%h))
      
      dbflag = 0
      
      rMB = 'rMB'; rSig = 'r'; wSig = 'w'; gSig = 'g'; pondSig = 'pond'
      colSig = 'col'
      scalar = 1
      
      TT = 0.0D0
      gid = w%g      !% Should only contain root level grids
      !rL=w(wid).CR; % Control River List, should have been sorted from upstream to downstream
      rL = w%CR
      !eL=w(wid).ER; % Exchange River
      eL = w%ER
      !Sol = w(wid).Sol;
      !COMP= Sol.COMP;
      !Exch= w(wid).Exch;
      !for i=1:length(rL), DTS(i)=r(rL(i)).Riv.dt; end
      do i = 1, size(rL)
        DTS(i) = r(rL(i))%Riv%dt
      enddo
      SCODE = 1
      !% weather data for the day
      sig = 0
      mask => w%DM%mask2

      Hour = 0
      nFT = nint(w%dt/g%OVN%dt)      !% flow time step is smaller than wtrshd time step
   
      do while (w%t < (T-PREC))
        !write(*,*) 'g.OVN.h = ', sum(g%OVN%h)
        !write(*,*) 'g.Veg.h = ', sum(g%Veg%h)
        !write(*,*) 'g.VDZ.h = ', sum(g%VDZ%h)
        !write(*,*) 'g.GW(1).h = ', sum(g%GW(1)%h)
        !write(*,*) 'r(25).Riv.h = ', sum(r(25)%Riv%h)
        !write(*,*) '*****************'        !%disp([JD H])
        Hour = Hour + 1       !% hour of day
        if (JD>=316 .and. Hour>=6) then
            !%pp=maxALL((g.GW(1).h-g.Topo.E).*w.DM.mask);%TV;
            ![JD H];
            dbflag = 1
        endif
        !hf = g.OVN.h; hc  = g.OVN.hc; hI = g.Veg.hI; h = hf + hc + hI; w.Rec.hSave = h; % debug
        !hf = g%OVN%h
        !hc = g%OVN%hc
        !hI = g%Veg%hI
        !h = hf + hc + hI
        if(.not. associated(w%Rec%hSave)) allocate(w%Rec%hSave(size(h,1),size(h,2)))
        !w%Rec%hSave = h
        !Sol.wea(wid,Hour,JD)
        time_start = 0d0
        call timer(time_start)
        call distWea(Hour,JD)
        
        call surface_proc(Hour)
        g%OVN%hc = g%VDZ%h(:,:,1)
        
        !call budgeter(pondSig,filter,mask,Hour,k)
       
        !TT(1)=TT(1)+toc;
        time_end = 0d0  
        call timer(time_end)
        TT(2) = TT(2) + (time_end - time_start)
          
        w%t = w%t + w%dt
        
        !k=1; filter(1)=197
        
        !call budgeter(wSig,scalar(1),mask,Hour,1)
        call find1D(w%DM%MAP,w%m(1)*w%m(2),filter,k)
        call budgeter(colSig,filter,mask,Hour,k)

        call gClearTFlux()
        if ( any( isnan(g%GW(1)%h) ) .or. any( isnan(g%VDZ%h) ) ) then
            call ch_assert('find nan in g.GW(1).h or g.VDZ.h')
        endif
        
      enddo
      
      end subroutine wtrshd_col
