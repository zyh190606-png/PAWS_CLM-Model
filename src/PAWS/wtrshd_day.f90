#include <preproc.h>
!

      subroutine wtrshd_day(T,JD,TT,sig)
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
      use surface_proc_mod, Only: surface_proc
      use budgeter_mod, Only: budgeter
      use memoryChunk, Only:tempStorage
      use pawsControl
      
      implicit none
      !global g r w
      real*8 :: T, TT(5)
      integer :: JD, sig
      real*8, parameter :: PREC = 1D-6
      real*8 :: gid, nFT, diff, mxv, RDT
      real*8, save:: lastDt =0D0
      integer, dimension(size(w%CR)) :: rL
      real*8, dimension(size(w%CR)) :: DTS
      integer, dimension(size(w%ER)) :: eL
      real*8, dimension(size(g%OVN%h,1),size(g%OVN%h,2)):: hf
      integer, dimension(size(g%OVN%h,1)*size(g%OVN%h,2)):: filter
      real*8, dimension(size(g%OVN%hc,1),size(g%OVN%hc,2)):: hc
      real*8, dimension(size(g%Veg%hI,1),size(g%Veg%hI,2)):: hI, h      
      integer :: i, j, SCODE, Hour, kk, m2(3), c, nstep
      character*32 :: rMB, rSig, wSig, gSig, pondSig,colSig
      real*8 :: time_start=0d0, time_end=0d0, dtNew, dtOld
      real*8, dimension(:,:), pointer :: mask=>NULL(), GWBCflux=>null(), hcAfterLand=>null(),highGWFix=>null()
      logical*1, dimension(:,:), pointer :: nanMap2D=>NULL()
      integer :: scalar(2)
      logical*1 :: subp, ovnp, rivp,dbgsave = .false.
      integer :: dbflag,k,kMap,k1,k2,kRiv
      
      dbflag = 0
      
      rMB = 'rMB'; rSig = 'r'; wSig = 'w'; gSig = 'g'; pondSig = 'pond'
      colSig = 'col'
      scalar = (/1, 1/)
      
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

      Hour = 0; NSTEP = 0
      nFT = nint(w%dt/g%OVN%dt)      !% flow time step is smaller than wtrshd time step
   
      do while (w%t < (T-PREC))
        !write(*,*) 'g.OVN.h = ', sum(g%OVN%h)
        !write(*,*) 'g.Veg.h = ', sum(g%Veg%h)
        !write(*,*) 'g.VDZ.h = ', sum(g%VDZ%h)
        !write(*,*) 'g.GW(1).h = ', sum(g%GW(1)%h)
        !write(*,*) 'r(25).Riv.h = ', sum(r(25)%Riv%h)
        !write(*,*) '*****************'        !%disp([JD H])
        Hour = floor((w%t-floor(w%t))*24D0+1D-10)+1       !% hour of day
        NSTEP = NSTEP + 1
        
        if (JD .eq. 244 .and. Hour .eq. 1) then
            !%pp=maxALL((g.GW(1).h-g.Topo.E).*w.DM.mask);%TV;
            ![JD H];
            dbgsave = .false.
            if (dbgsave) call saveMat2('debugDay.mat')
        else 
            dbgsave = .false.
        endif
        
        !hf = g.OVN.h; hc  = g.OVN.hc; hI = g.Veg.hI; h = hf + hc + hI; w.Rec.hSave = h; % debug
        hf = g%OVN%h
        hc = g%OVN%hc
        hI = g%Veg%hI
        h = hf + hc + hI
        if(.not. associated(w%Rec%hSave)) allocate(w%Rec%hSave(size(h,1),size(h,2)))
        w%Rec%hSave = h
        !Sol.wea(wid,Hour,JD)
        call distWea(Hour,JD)
        
        !%releaseSchedule(wid,H,JD)
        !do i = 1, size(w%child)
        !    c = w(wid).child(i);
        !    w(c).watershed(c);
        !enddo

        
        !%GW_source(gid,1);
#if (defined CLM_Only)
#else
        if (SCODE == 0) then
          g%VDZ%h(:,:,1) = g%OVN%h        !% do not use F_contribute
        else
          g%VDZ%h(:,:,1) = g%OVN%hc
        endif
#endif
        !tic
        time_start = 0d0
        call timer(time_start)
        !Sol.veg(wid,H) % interception, ET, unconfined aquifer
        call surface_proc(Hour)
        
        if (dbgsave) call saveMat2('afterLand.mat')
        
#ifndef CLM_Only
        !TT(1)=TT(1)+toc;
        time_end = 0d0
        call timer(time_end)
        TT(1) = TT(1) + (time_end - time_start)
        
        if (SCODE == 0) then
            g%OVN%h = g%VDZ%h(:,:,1)
        else
            g%OVN%hc = g%VDZ%h(:,:,1)       ! WHAT IF hc is negative (from VDZ)
        endif
        
        call tempStorage('hcAfterLand',hcAfterLand,1)
        hcAfterLand = g%VDZ%h(:,:,1) ! for mapOutput2 calculation
        call tempStorage('highGWFix',highGWFix,1)
        
        
        !if (no_Surf) goto 200 !! THIS NEEDS TO BE CHANGED LATER!!!
        
        !budgeter('pond',wid,[],H)           ! debug        
        !call budgeter(pondSig,1,mask,Hour)
    
        !% Operator Splitting: First infiltrate then OVN flow
        !%disp([JD H g.GW(2).h(i,j) g.GW(1).h(i,j) g.VDZ.DF(i,j) g.VDZ.Dperc(i,j) g.VDZ.h(i,j,2) g.VDZ.SRC(i,j,2)])
        subp = .false.
        ovnp = .false.
        !if ~isempty(find(isnan(g(gid).GW(1).h), 1)) ||...
        !        any(isnan(g(gid).VDZ.h(:))) !%|| anyALL((g(gid).GW(1).h(:)>310))
        if ( any( isnan(g%GW(1)%h) ) .or. any( isnan(g%VDZ%h) ) ) then
            subp = .true.
            sig = 1
            call find1DL(isnan(g%GW(1)%h),w%m(1)*w%m(2),filter,k1)
            call display( 'JD: '//trim(num2str(JD))//' hour:'//trim(num2str(hour)))
            call display( 'Subsurface NaN! GW k1:'//trim(num2str(k1))// 'Location(1)'//trim(num2str(filter(1))))
            call find1DL(isnan(g%VDZ%h),w%m(1)*w%m(2)*w%m(3),filter,k2)
            call display( 'Subsurface NaN! VDZ k2:',k2,'Location(1)',filter(1))
            return
        endif
        !loc = (g(gid).GW(1).h>g(gid).GW(1).topBound) & w.DM.mask;
        !if anyALL(loc)
        !    diff = g(gid).GW(1).h(loc)-g(gid).GW(1).topBound(loc);
        !    g(gid).GW(1).h(loc) = g(gid).GW(1).topBound(loc);
        !    g(gid).OVN.hc(loc) = g(gid).OVN.hc(loc)+g(gid).GW(1).ST(loc).*diff;
        !end
        do j = 1, size(g%GW(1)%h,2)
          do i = 1, size(g%GW(1)%h,1)
            if(w%DM%mask(i,j) .AND. g%GW(1)%h(i,j) > g%GW(1)%topBound(i,j)) then              
              diff = g%GW(1)%h(i,j)-g%GW(1)%topBound(i,j)
              g%GW(1)%h(i,j) = g%GW(1)%topBound(i,j)
              highGWFix(i,j) = g%GW(1)%ST(i,j)*diff
              g%OVN%hc(i,j)  = g%OVN%hc(i,j)+highGWFix(i,j)
            else
              highGWFix(i,j) = 0D0
            endif
          enddo
        enddo
    
        do kk = 1, int(nFT)
          RDT = g%OVN%dt
          !tic
          time_start = 0d0
          call timer(time_start)
          if (SCODE == 0) then  !% not using F_contribute
          else
            !%F_contribute(gid)
          endif
          !Sol.ovn(gid);
          if (g%VDZ%dtP(3) > 0) then
            g%OVN%S = g%OVN%Rf * (1.0D0 / g%VDZ%dtP(3) / 86400.0D0)     !% The unit of Rf is m (come in g.VDZ.dtP(3) time)
          else
          endif
          
          !prism(6,g.OVN,g.Topo,g.DM.d)
          !call prism(6, g%OVN, g%Topo, g%DM%d)
          CALL ovn_cv(g%OVN%ms,g%OVN%h,g%OVN%U,g%OVN%V,g%OVN%S,g%DM%d,g%Topo%E,g%Topo%Ex,g%Topo%Ey    &
         &    , g%OVN%Mann, w%DM%Mask,g%OVN%dt,lastDt, size(g%OVN%Cond),g%OVN%Cond)
         
          g%OVN%t = g%OVN%t + g%OVN%dt
          
          !TT(2)=TT(2)+toc;
          time_end = 0d0  
          call timer(time_end)
          TT(2) = TT(2) + (time_end - time_start)
          !%[r(27).Riv.h' r(27).Riv.Qx']
          !%test=0;
          !%if min(g.OVN.h(:))<-0.02%g.OVN.t>=7.328091805558624e+005-2*g.OVN.dt-1e-5;%
          !%if g.OVN.t>=7.320321249999766e+005-g.VDZ.dt-1e-4;%datenum('31-Aug-2007 08:30:00')-1e-5
          !%    4;%TV;
          !%end
        
          subp = .false.
          !if ~isempty(find(isnan(g(gid).OVN.h), 1))
          if (any(isnan(g%OVN%h))) then
            call find1DL(isnan(g%OVN%h),w%m(1)*w%m(2),filter,k)
            call display( 'OVN NaN!!! k:',k,' location(1): ',filter(1))
            ovnp = .true.
            sig = 1
          endif
          !tic;
          time_start = 0d0          
          call timer(time_start)
          kRiv = 0
          do while (RDT > PREC)     !% River Network Module. This module needs to adaptively select time step.

            !write(777,'(F13.5,F20.16)') w%t, r(25)%Riv%umax
            !if(w%t>=731095.08D0) then
            !    write(*,*) 'err'
            !endif
            ![dtNew, dtOld] = rMinDtSet(rL,RDT);
            kRiv = kRiv+1
            call rMinDtSet(rL, RDT, dtNew, dtOld)
            
            !rAvgStages(wid);
            call rAvgStages()
            if(.not. associated(g%Riv%Qoc)) allocate(g%Riv%Qoc(size(rL)))
            do i = 1, size(rL)
                !Exch.OC(gid,i,dtNew);
                call F_oc_dw(i)
            enddo
            do i = 1, size(rL)
                c = rL(i)
                
                !%Exch.GC(gid)
                
                call rLateral(c)
                
                !Sol.rte(c);
                call r_MacCormack(c)
                !%disp([H kk c])
                !%if w.signals(1)==1
                !  % transport
                !  %runQ(r(c).q); % March
                !%end
                !%if c==30, 
                !    %budgeter('r',c); 
                !%end
                !recAcc('rMB',c)
                call recAcc(rMB,c)
                !rClearTFlux(c)
                call rClearTFlux(c)
            enddo
            !recAcc('r')
            call recAcc(rSig,0)
            !%budgeter('OVN',1,[],H); 
            rivp = .false.
            do j = 1, size(w%CR)
                !if (~isempty(find(isnan(r(w.CR(j)).Riv.h), 1)) || ~isempty(find(r(w.CR(j)).Riv.h<0, 1)))
                if(any(isnan(r(w%CR(j))%Riv%h)) .or. any(r(w%CR(j))%Riv%h<0)) then
                    rivp = .true.
                    sig = 1
                    call find1DL(isnan(r(w%CR(j))%Riv%h),size(r(w%CR(j))%Riv%h),filter,k)
                    call display( 'River NaN!!!, j:', j, ' k: ',k, ' location(1):',filter(1))
                    return
                endif
            enddo
            if (any(isnan(g%OVN%h))) then
                ovnp = .true.
                call find1DL(isnan(g%OVN%h),w%m(1)*w%m(2),filter,k)
                call display( 'OVN Exchange NaN!!! k:', k,' location(1): ', filter(1))
            endif
            
            !if rivp  || ~isempty(find(isnan(g(gid).OVN.h), 1)) || subp || ovnp
            if(rivp  .or. subp .or. ovnp) then
                !warning('4');
                call ch_assert('wtrshd_day :: NaN detected')
                return
            endif
            RDT = RDT - dtNew
          enddo
          do i = 1, size(rL)
            r(rL(i))%Riv%dt = DTS(i)
          enddo   !% TEMPORARY
          !TT(3)=TT(3)+toc;        
        enddo
          
        time_end = 0d0          
        call timer(time_end)
        TT(3) = TT(3) + (time_end - time_start)
        
        if (dbgsave) then
          call saveMat2('afterSurfaceFlow.mat')
        endif
        !%disp(['beyond surface ',num2str(H)])
        !% Saturated Zone. 
        !% ############################### GW(2)???
        !%for i=1:length(g(gid).GW)-1
        !tic
 200   time_start = 0d0        
        call timer(time_start)
        do i = 1, size(g%GW)
            if (abs(g%GW(i)%tt+g%GW(i)%dt-(w%t + w%dt)) < PREC) then
                !% Sum up recharge (or release)
                !GW_source(gid,i);
                call GW_source(i)
                !Sol.gw(gid,i);
                call GW_Sol(i)
            endif
        enddo
        
        call tempStorage('GWBCFlux',GWBCflux,1)
        call fixedHeadGWupdate(g%GW(1)%ST,int(g%DM%msize),g%GW(1)%h,w%tData%h1save,GWBCflux)
        if (dbgsave)then
          call saveMat2('afterGW.mat')
        endif
        
        !%ii=13;jj=14; [g.GW(1).W(ii,jj,:) g.GW(1).h(ii,jj,:) g.GW(2).h(ii,jj) g.GW(2).ST(ii,jj)]
        !TT(4)=TT(4)+toc;
        time_end = 0d0        
        call timer(time_end)
        TT(4) = TT(4) + (time_end - time_start)
        call columnOutput(0)
        call recAcc(gSig,0)

        !tic
        time_start = 0d0        
        call timer(time_start)
        !budgeter('w',wid,[],H); % BASIN Average mass balance output
        call writeRec()
        call budgeter(wSig,scalar(1),mask,Hour,1)
#endif        

        !k=1; filter(1)=483
        if (dbgsave) THEN
        call find1D(w%DM%MAP,w%m(1)*w%m(2),filter,k)
        call budgeter(colSig,filter,mask,Hour,k)
        endif
        !mapOutput(0); % Temporal Average spatial maps
        !call mapOutput(0, w%tData%maps)
        
        ! depending on many objects there are, evenly put them in 
        kMap = min(size(w%tData%maps), ceiling(JD*size(w%tData%maps)/365D0))
        call mapOutput2(0, w%tData%maps(kMap))
        
        !writeRec % recording
        !TT(5)=TT(5)+toc;
        time_end = 0d0        
        call timer(time_end)
        TT(5) = TT(5) + (time_end - time_start)
        w%t = w%t + w%dt

        !%movOutput(0);  %wtrshd_day, calendarMarch to do movie
        !% save for display::::
        if(.not. associated(w%tData%m3%R)) then
          allocate(w%tData%m3%R(size(g%VDZ%Dperc,1),size(g%VDZ%Dperc,2)))
        endif
        w%tData%m3%R = g%VDZ%Dperc
        if(.not. associated(w%tData%m3%E)) then
          allocate(w%tData%m3%E(size(g%Veg%ET,1),size(g%Veg%ET,2)))
        endif        
        w%tData%m3%E = g%Veg%ET
        !gClearTFlux(gid)
        call gClearTFlux()
        !% Transport
        !%checkEndDt(wid);
        
        !%figure(20); i=25; plot(r(i).DM.x,r(i).Riv.E+r(i).Riv.h); hold on; plot(r(i).DM.x,r(i).Riv.hg,'r'); plot(r(i).DM.x,r(i).Riv.E,'k'); hold off
        !%figure(30); ii=[19 15 14 17]; jj=[27 19 18 21]; cols='brgkmybrgkmy';
        !%plotColumn(ii,jj,cols);
        !%TV       
      enddo
      
      end subroutine wtrshd_day
