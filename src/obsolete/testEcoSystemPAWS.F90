      PROGRAM testEcoSystemPAWS
      USE VData
      USE clmtype
      USE clmtypeInitMod
      USE paws_varctl, Only: WRITEOUT
      type(par_type),dimension(1)::par
      CHARACTER*(30),file,file2,prj
      TYPE(GlobalG),pointer::gg
      integer status, nargin,m(3),init
      
      nargin = NARGS( )
      if (nargin > 1) then
      CALL GETARG (1, file2, status)
      file = file2
      else
      file = 'data_in.mat'
      endif
      
      call bridgeGlobal(file)
      
      gg => g
      
      
      w%DM%mask = .false.
      w%DM%mask(20,11:26)=.true. 
      
      !w%DM%mask(20,20) = .true.
      WHERE (w%DM%mask)
       g%VDZ%MAPP = 1D0
      ELSEWHERE
       g%VDZ%MAPP = 0D0
      ENDWHERE
      w%timezone = -5 ! Eastern Standard
      
      !g%Veg%RPT(20,20,1:3)=(/3,4,6/)  ! original code -- tested working
      DO I = 1,16
      g%Veg%RPT(20,10+i,1:3)=3
      ENDDO
      call genVegDat(g%Veg%VDB)
      call dailyVegUpdate(243) 
      
      DO I = 1,16
      g%Veg%RPT(20,10+i,1:3)=i
      g%Veg%Imp(20,i) = 0D0
      ENDDO
      
      !g%Veg%RPT(20,20,1:3)=(/21,7,13/) ! CLM code  ! original code -- tested working
      !g%Veg%Imp(20,20) = g%Veg%Frac(20,20,1)       ! original code -- tested working
      
      m = (/g%VDZ%m(1),g%VDZ%m(2),g%VDZ%m(3)/)
      allocate(g%VDZ%ICE(m(1),m(2),m(3)),g%VDZ%sand(m(1),m(2),m(3)),g%VDZ%clay(m(1),m(2),m(3)),  &
     &    g%VDZ%om(m(1),m(2),m(3)),g%VDZ%soic(m(1),m(2)),g%VDZ%gti(m(1),m(2)))
      
      g%VDZ%ICE = 0D0
      g%VDZ%sand = 50D0;
      g%VDZ%clay = 42D0
      g%VDZ%OM  = 0.25D0/100D0*2650D0 ! convert % to kg/m3, times the soil bulk density!
      g%VDZ%soic = 1
      g%VDZ%gti = 0.42D0;
      prj = 'testCLM'
      Env%recfile = (/trim(prj)//'.txt',trim(prj)//'_Rec.txt2'/)
      
      WRITEOUT = .FALSE.
      do i = 1,4
      WRITE(*,*) I
      IF (i .eq. 4) WRITEOUT = .TRUE.
      call calendarMarch(par,0,i);
      call wSyncTime(int(w%ModelStartTime))
      
      enddo
      
      END PROGRAM
      
      
      
      
      
      
      
      