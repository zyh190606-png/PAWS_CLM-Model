
!
      subroutine calendarMarch(par,npar,init)
      !% (E) this function runs the model according to the calendar year/month/day
      !% it facilitates the scheduling of events
      !% it also applies essential operations like initialization, parameter
      !% change and calls adhoc_set
      !% outputs: budgeter:recfile(1)-Basin average fluxes (time series)
      !%          recAcc and writeRec: recfile(2)-Monitoring point (time series)
      !%          mapOutput: results saved in final mat file (time averaged Spatial maps)
      !%          movOutput: spatial maps saved over time
      !%          columnOutput: cell recording saved in final mat file w.tData.S1

      !global w g r Env Prob mov
      use Vdata
      USE paws_varctl, Only: WRITEOUT, WRITEOUTL, nyrs, ndays, nmonths, mode, nrepeat, nMatOutputFreq
      USE wMonthlyAverages_Mod, Only: wMonthlyAverages
      implicit none
      !  if (ispc && ~isempty(varargin)),
      !      handles=varargin{1}; timeP=1;
      !      hp=handles.togglebutton_pause;
      !      hc=handles.edit_current;
      !  else
      !      timeP = 0; % Non-GUI execution
      !  end
      integer :: npar,init
      type(par_type), dimension(npar) :: par
      real*8 :: Prec = 1e-8_8;
      real*8 :: PlotT, TT(5), r1
      integer :: datevec(3),Y1,M1,D1,Y2,M2,D2,gid,sig,Y,Mend,Ms,S,S2,JDS,M,Ds,Dend,D,JD,kkk, i
      character*32 :: date_char, time_char, monthf, yearlyf, temp, YY, MM
      logical,save:: columnMode = .false.,debug3=.false.
      integer :: doRuns = 1, doRunIter
#if (defined CLM_Only)
      columnMode = .true.
#endif
      !  %if (isfield(w(wid).Rec.tseries(1),'t'))
      !  %    tt = w(wid).Rec.tseries(1).t; QQ=w(wid).Rec.tseries(1).v;
      !  %else
      !  %    tt=[]; QQ=[]; S1=[];
      !  %end
      PlotT = w%plotT
      Env%KK = 0.0D0
      !S1 = [];
      !tic;
      !  %if length(w(wid).Rec.tseries)<3 || isempty(w(wid).Rec.tseries(3).v)
      !  %    w(wid).Rec.tseries(3).v=[];
      !  %end
      Env%t1 = 0.0D0
      Env%t2 = 0.0D0
      Env%t3 = 0.0D0
      datevec = datevec2(nint(w%t))
      Y1 = datevec(1)
      M1 = datevec(2)
      D1 = datevec(3)
      datevec = datevec2(nint(w%ModelEndTime))
      Y2 = datevec(1)
      M2 = datevec(2)
      D2 = datevec(3)

      gid = w%g
      sig = 1
      !if ~isfield(w.tData,'S1')
      !  w.tData.S1= [];
      !end if
      
      ! debug use:
      if (TRIM(mode) .eq. 'column')  columnMode = .true.
      doRuns =  nrepeat
      if (columnMode) then
        writeoutL = .true.
      endif
      !!!!!!!!!!TEMPORARY
      !! TESTING
      !g%GW(1)%K = 1D-12
      !g%GW(2)%K = 1D-12
      
      
      call allocVData
      call display( 'Base time step is: ',w%dt)
      if (w%t == w%ModelStartTime .and. init .eq. 1) then
        !call adhoc_Set
        !changepar(par)
        call changepar2(-1, par, npar)
        call display( 'Changing Parameters Done')
        !if ~isfield(g.VDZ,'GA') || isempty(g.VDZ.GA)
        !  g.VDZ.GA = GreenAmptMatrix(0,50,g.VDZ);
        !end if
        if(.not. associated(g%VDZ%GA%h)) then
          call GreenAmptMatrix(opt=0,m=50.0D0,VDZ=g%VDZ,GA=g%VDZ%GA)
        end if
        ![g.VDZ.GA.FM,g.VDZ.GA.KBar] = GreenAmptMatrix(2,g.VDZ);
        call GreenAmptMatrix(opt=2,VDZ=g%VDZ,FM=g%VDZ%GA%FM,KBar=g%VDZ%GA%KBar)
        call initValues()
        call display( 'InitValues done')
        !w.tData.par = par; % record the changes that are made

        !Prob.dem = []; Prob.ned = []; % remove two large fields to reduce saving data amount
        ! Don't do these things if it is a continuation of jobs
      end if
      
      ! OVERALL ADJUSTMENT SECTION
      ! THESE ARE FOR Hypothetical scenarios testing
      !g%Veg%RPT = 7 ! all of which are trees
      
      !openFiles % will re-open files even if not at the start time
      call openFiles()

      !kkk=find(Env.recfile{1}=='_');
      kkk = 0
      kkk = index(Env%recfile(1),'_')
      !if ~isempty(kkk), monthf =  ['month_',Env.recfile{1}(1:kkk(1)-1)];
      !else monthf =  ['month']; end
      if(kkk/=0) then
        temp = Env%recfile(1)
        monthf =  'month_'//temp(1:kkk-1)
      else
        monthf = 'month'
      endif
      call display( 'Entering time loop')
      call mapOutput2(10,w%tData%maps(1)) 

      do doRunIter = 1, doRuns
      !% doRuns loop%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      do Y = Y1, Y2
      !% year loop%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (Y == Y2) then
          Mend = M2
          S2 = 2
        else
          Mend = 12
          S2 = 0
        endif
        if (Y == Y1) then
          Ms = M1
          S = 1
        else
          Ms = 1
          S = 0
        endif
        !% At the beginning of year, schedule some variables
        !% Schedule Solar Radiation
        call scheduleSolar(Y)
        !% Schedule Plant Growth (Temporary)
        !g(gid).Veg.VDB=genVegDat(g(gid).Veg.VDB);
        !call genVegDat(g%Veg%VDB)
        !JDS = datenum(num2str(Y*10000+101),'yyyymmdd');
        JDS = datenum2((/Y*10000+101/))

        do M = Ms, Mend
        !% monthly loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          if ((S==1) .and. (M==M1)) then
            Ds = D1
          else
            Ds = 1
          endif
          if ((S2==2) .and. (M==M2)) then
            Dend = D2
          else
            Dend = eomday(Y,M)
          endif
          !g(gid).Wea.month=wMonthlyAverages(M);
          call wMonthlyAverages(M, g%Wea%month)

        !    %{
        !    seq = [11   12   1   2   3 4];
        !    fac = [0.9 0.7 0.3 0.2 0.5 0.8];
        !    if ismember(M,seq)
        !        km = find(seq == M);
        !        g.VDZ.KS = g.VDZ.KSSAVE*fac(km);
        !    else
        !        g.VDZ.KS = g.VDZ.KSSAVE+0;
        !    end
        !    %}


          do D = Ds, Dend     !% daily loop%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Env%KK = Env%KK + 1.0D0
            !%[dd,JD]=ymd2jd(round(w.t),2);
            !JD = ymd2jd([Y M D],3);
            call ymd2jd((/Y, M, D/), 3, JD)
#if (defined CLM)
#else            
            !if (dmod(Env%KK,2.0D0)==1.0D0) call dailyVegUpdate(JD)    !% update canopy state every 2 days
#endif
           call dailyVegUpdate2(JD)   ! This one only changes ho -- CURRENT RPT code is not right
        !
            !%tic;
            if (columnMode) then
              call wtrshd_col(w%t+1.0D0, JD, TT, sig); ! trying new code
            else
              call wtrshd_day(w%t+1.0D0, JD, TT, sig);
            endif
            
                !i = 35; j = 12;
                !%wtrshd_column(wid, w(wid).t+1,JD,i,j);
                !  %output=columnOutput(i,j,1); w.tData.S1 = [w.tData.S1; output];
                !%wtrshd_wea(wid, w(wid).t+1,JD);
                !%wtrshd_surface(wid, w(wid).t+1,JD);
                !%toc;
                !disp([Y JD TT, sum(TT)]) % Run a day
                !write(*,'(2I6,6F10.4)') Y, JD, TT, sum(TT)
                r1 = sum(TT)
                call display(num2str(Y),num2str(JD),num2str(TT,5),num2str(r1) )
                !%disp([Y JD])
                !% Set Current Time Flag
                !if (timeP==1), set(hc,'String',datestr(w(1).t,'yyyy-mm-dd HH:MM:SS')); end
                !% Read Pause Flag
                !if (timeP==1 && get(hp,'Value') == get(hp,'Max'))
                !    set(hp,'Value',get(hp,'Min'));
                !    break
                !end
                !%F = movOutput(1);
                !%if (D==15), cmd_Save('month.mat'); end
                !wSyncTime(wid,JDS+JD)   % sync time to remove small time roundoff error
                call wSyncTime(JDS+JD)
!                if (JD>=281) then
!                  if (debug3) call savemat2('daily.mat')
!                endif
                if (.not. columnMode) call timeSeriesMapControl()
          enddo
        !    !% Compute Monthly Average Values
        !    if sig==0
        !        cmd_Save(monthf)
        !    end
            !datestr(now)
            call date_and_time(date_char, time_char)
            call display( date_char(1:4)//'/'//date_char(5:6)//'/'//date_char(7:8), &
     &      '  ', time_char(1:2)//':'//time_char(3:4)//':'//time_char(5:6))
            if (sig .eq. 0) then
            ! the purpose of saving monthly mat file is for debugging (so that you are restart near where things went wrong
            ! so, we don't want to over-write the monthly save if things did go wrong. Only save it when s == 0
            if (.not. columnMode .and. nMatOutputFreq .eq. 1) call savemat2(trim(monthf)//'.mat')
            endif
            !%if (mod(M,4)==0)
            !%    cmd_Save('3monthly.mat')
            !%end
            write(YY,*) Y
            write(MM,*) M
            call savemat2(trim(YY)//'_'//trim(MM)//'.mat')
        enddo
        !% Calculate Yearly Average Values
        !call savemat2(trim(monthf)//'.mat')
        
        nyrs = nyrs + 1
      enddo
      
      
        if (columnMode) then
           datevec = datevec2(nint(w%ModelStartTime))
           
           call wSyncTime( datenum2((/(datevec(1)+1)*10000 + 101/) ))
            Y1 = datevec(1)+1
            M1 = 1
            D1 = 1
           if (doRunIter .eq. (doRuns - 1)) then
              WRITEOUT =  .true. ! for initial states
           else
              WRITEOUT = .false.
           endif
        endif
      
      enddo
      !% doRuns loop%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !%save rec S1
      !%movie2avi(mov.R,'Recharge.avi','compression','none')
      !%movie2avi(mov.E,'ET.avi','compression','none')

      end subroutine calendarMarch
