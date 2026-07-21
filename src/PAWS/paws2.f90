      PROGRAM testEcoSystemPAWS
      USE VData
      USE clmtype
      USE clmtypeInitMod
      type(par_type),dimension(1)::par
      CHARACTER*(30),file,file2,prj
      TYPE(GlobalG),pointer::gg
      integer status, nargin,m(3)
      
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
      w%DM%mask(10,12)=.true. 
      WHERE (w%DM%mask)
       g%VDZ%MAPP = 1D0
      ELSEWHERE
       g%VDZ%MAPP = 0D0
      ENDWHERE
      w%timezone = -5 ! Eastern Standard
      
      !g%Veg%RPT(10,12,1:3)=(/3,4,6/)  ! original code
!      call genVegDat(g%Veg%VDB)
!      call dailyVegUpdate(243) 
      
      !g%Veg%RPT(20,20,1:3)=(/1,6,13/) ! CLM code
      
!      m = (/g%VDZ%m(1),g%VDZ%m(2),g%VDZ%m(3)/)
!      allocate(g%VDZ%ICE(m(1),m(2),m(3)),g%VDZ%sand(m(1),m(2),m(3)),g%VDZ%clay(m(1),m(2),m(3)),  &
!     &    g%VDZ%om(m(1),m(2),m(3)),g%VDZ%soic(m(1),m(2)),g%VDZ%gti(m(1),m(2)))
!      g%VDZ%ICE = 0D0
!      g%VDZ%sand = 50D0;
!      g%VDZ%clay = 42D0
!      g%VDZ%om  = 0.25D0/100D0*2650D0 ! convert % to kg/m3, times the soil bulk density!
!      g%VDZ%soic = 1
!      g%VDZ%gti = 0.42D0;
      prj = 'testCLM'
      Env%recfile = (/trim(prj)//'.txt',trim(prj)//'_Rec.txt'/)
      
      call calendarMarch(par,0,1);
      
      END PROGRAM
      
      
      
      
      
      
      
      