MODULE surface_proc_mod
CONTAINS
!====================================================================================
      subroutine surface_proc(H,varin)
      USE VData
      !% (E) land surface processes (vegetation, ET) and subsurface processes-
      !% (soil water flow).
      !% This function assembles some input in VDZ.Cond, and then calls the-
      !% fortran interface with code = 1
      !% Input:
      !% H: hour of day; varargin: irrelavant for watershed application
      !% Output:
      !% updated g.VDZ, g.GW.h and g.Veg
      implicit none
      !global w g
      integer,optional :: varin ! changed 06/20/13, optional statement for Unsat
      integer :: H, nr, nc, nz, s1, s2, i, j, mtd
      real*8 :: gid
      integer, dimension(size(g%VDZ%WTL,1),size(g%VDZ%WTL,2)) :: WTL

      real*8, dimension(size(g%GW(1)%h,1),size(g%GW(1)%h,2))::ES
      real*8, dimension(size(g%GW(1)%EB,1),size(g%GW(1)%EB,2))::EB
      real*8, dimension(size(g%OVN%SN,1),size(g%OVN%SN,2))::k

      integer, save :: state = 0, mode = 0;
      integer idx,idy

      !persistent state
      gid = w%g
      !%HB = g(gid).GW(1).h;
      IF (PRESENT(varin)) mode = varin
      !% ! [1 ST,2 zl (=DZ(uc)/2+D(aquitard)),3 DR, 4, Kl, 5, hl(pressure head), 6
      !% E(center) uc aquifer]
      !if isempty(state) || w(wid).t == w(wid).ModelStartTime
      !if (state==0 .or. w%t == w%ModelStartTime) then
      idy = int(g%VDZ%m(2)/2); idx = int(g%VDZ%m(1)/2); 
      if (g%VDZ%Cond(idy,idx,6) .eq. 0 .OR. w%t == w%ModelStartTime) then
      !%warning('quick access disabled for testing GGA')
        !% g(gid).VDZ.Cond(:,:,1)=g(gid).GW(1).ST; NEED TO BE UPDATED FOR NLC
        !%g(gid).VDZ.Cond(:,:,2)=g(gid).GW(1).h - g(gid).GW(1).EB;
        !%g(gid).VDZ.Cond(:,:,2)=g(gid).GW(1).DZ/2+g(gid).GW(1).botD;
        ![nr,nc,nz]=size(g.VDZ.h);
        nr = size(g%VDZ%h,1)
        nc = size(g%VDZ%h,2)
        nz = size(g%VDZ%h,3)
        if (size(g%GW) > 1) then
          g%VDZ%Cond(:,:,2) = g%GW(1)%E - g%GW(2)%E   !% dzl
          g%VDZ%Cond(:,:,4) = g%GW(1)%botK
        else
          g%VDZ%Cond(:,:,2) = g%GW(1)%E - g%GW(1)%EB
          g%VDZ%Cond(:,:,4) = 0D0
        endif
        !WTL = g(gid).VDZ.WTL; ES = zeros(size(g(gid).GW(1).h));
        WTL = g%VDZ%WTL
        ES = 0.0D0
        do j = 1, nc
            do i = 1, nr
                ES(i,j) = g%VDZ%EB(i,j,WTL(i,j)-1)
            enddo
        enddo
        !EB = g(gid).GW(1).EB;
        EB = g%GW(1)%EB
        g%VDZ%Cond(:,:,6) = (EB + ES) / 2.0D0
        if (size(g%GW) > 1) g%VDZ%Cond(:,:,2) = g%VDZ%Cond(:,:,6) - g%GW(2)%E
        !%g(gid).VDZ.Cond(:,:,7)=g(gid).GW(1).DZ; % OLD
        g%VDZ%Cond(:,:,7) = ES - EB   !%GW(1).DZ is the total thickness of Topo.E-GW(1).EB, D is the saturated thickness
        g%VDZ%Cond(:,:,8) = g%GW(1)%EB
        !state = 1
      endif
      g%VDZ%Cond(:,:,1) = g%GW(1)%ST
      g%VDZ%Cond(:,:,3) = g%GW(1)%DR
      if (size(g%GW) > 1) g%VDZ%Cond(:,:,5) = g%GW(2)%h - g%GW(2)%E
      ![nr,nc,nz]=size(g(gid).VDZ.R);
      nr = size(g%VDZ%R,1)
      nc = size(g%VDZ%R,2)
      nz = size(g%VDZ%R,3)
      if (nz<2) then
        !g(gid).VDZ.R = zeros(nr,nc,4);
        if(associated(g%VDZ%R)) deallocate(g%VDZ%R)
        allocate(g%VDZ%R(nr,nc,4))
        g%VDZ%R = 0.0D0
      endif
      !if ~isfield(g.OVN,'OVF')  % This is to deal with mat files saved before
      if (.not. associated(g%OVN%OVF)) then
        !g.OVN.OVF = zeros(nr,nc)+0.05;
        allocate(g%OVN%OVF(nr,nc))
        g%OVN%OVF = 0.05D0
      endif
      !if ~isfield(g.VDZ.th,'state')
      if (.not. associated(g%VDZ%th%state)) then
          !g.VDZ.th.state = zeros(1)+0;;
          allocate(g%VDZ%th%state)
          g%VDZ%th%state = 0.0D0
      elseif (g%VDZ%th%state > 0) then
          g%VDZ%th%state = 0.0D0
          g%VDZ%th%Soutv(:,:,[1, 7, 22]) = 0.0D0
      endif
      
      !% R(1):botdrain dt; R(2):saturation signal;
      ! INTERFACE CHANGED.
      !if length(varargin)>0 && isstruct(varargin{1})
      !if (present(varin)) then
        !VDZ = varargin{1};
        !warning('101 in surface_proc for testing GGA');
      !  call display( '101 in surface_proc for testing GGA')
      !  mtd = 101
      !    prism(101,g(gid).Veg.t+g(gid).Veg.dt,H,g(gid).Veg,VDZ,g(gid).Wea.day,g(gid).OVN,g(gid).GW(1).h);
      !  CALL UNSAT(g%VDZ%MAPP,g%VDZ,g%OVN%H,g%GW(1)%h) !!!!!!!! TEMP!
      !else
        mtd = 1
      !    prism(mtd,g(gid).Veg.t+g(gid).Veg.dt,H,g(gid).Veg,g(gid).VDZ,g(gid).Wea.day,g(gid).OVN,g(gid).GW(1).h);
      IF (mode .eq. 0) then
        CALL land(g%Veg%t+g%Veg%dt,int(H),g%GW(1)%h)
      ELSE
        CALL UNSAT(g%VDZ%MAPP,g%VDZ, k, g%GW(1)%H, mode) ! k will not be used, just passing around
      ENDIF
      !endif
      !% If depression storage is allowed, this will modify HB


      !% debug
      !s1 = sumALL(g.VDZ.K(:,:,1)==0 & g.VDZ.K(:,:,2)==1); % h1negative
      !s2 = sumALL(g.VDZ.K(:,:,1)==0 & g.VDZ.K(:,:,2)==2); % dt decrease
!      s1 = 0
!      s2 = 0
!      do j = 1, size(g%VDZ%K(:,:,1),2)
!        do i = 1, size(g%VDZ%K(:,:,1),1)
!          if (g%VDZ%K(i,j,1)==0 .and. g%VDZ%K(i,j,2)==1) then
!            s1 = s1 + 1
!          endif
!          if (g%VDZ%K(i,j,1)==0 .and. g%VDZ%K(i,j,2)==2) then
!            s2 = s2 + 1
!          endif
!        enddo
!      enddo
      !if ~isfield(w.Rec,'tss')
      if(.not. ASSOCIATED(w%Rec%tss)) then
          allocate(w%Rec%tss(3))
          w%Rec%tss = [0, 0, 0]
      endif
      !if (s1>0 | s2>0)
      !    4;
      !end
      w%Rec%tss = w%Rec%tss + dble([s1,  s2,  1])
      !% debug end
      !%g(gid).OVN.S = g(gid).OVN.Rf/(g.OVN.dt*86400); % Runoff becomes source term to surface flow
      !% also note the unit in ovn_cv is m/s
      !%g(gid).Veg.t = g(gid).Veg.t + g(gid).Veg.dt; %done inside the solver


      end subroutine surface_proc
END MODULE
