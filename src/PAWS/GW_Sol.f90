!

      subroutine GW_sol(iGW)
      !% (E) calls fortran program aquifer2D to solve the groundwater flow-
      !% equation in quasi-3D fashion
      !% Input: iGW (aquifer layer number)
      !% Output: 
      !% g.GW.h; g.GW.Qperc (Recharge to lower aquifer)
      !global g w
      use Vdata
      implicit none
      integer :: iGW, i
      real*8 :: gid  
      real*8, dimension(size(w%DM%mask,1),size(w%DM%mask,2)) :: MAPP
      real*8, dimension(size(g%GW(iGW)%h,1),size(g%GW(iGW)%h,2)) :: Dperc

      real*8, dimension(:), pointer :: Dperc1D=>NULL(), gGWh1D1=>NULL(), gGWh1D2=>NULL()
      integer,DIMENSION(:),POINTER::idx=>NULL()
      type(BC_type),dimension(:),pointer :: Cond=>NULL()
      real*8 h0, ww, h1
      
      gid = w%g

      !MAPP = double(w(wid).DM.mask);
      MAPP = 0.0D0
      where(W%DM%mask)
      	MAPP = 1.0D0
      end where
      !%aquifer2D(MAPP,g(gid).GW(iGW));

      !% print debug
      !% h0 = meanALL((g.GW(1).h-g.GW(1).EB).*g.GW(1).ST.*MAPP);
      !% ww  = meanALL(g.GW(1).W.*MAPP)*g.GW(1).dt;
!      h0 = sum(g%GW(iGW)%h*g.GW(iGW).ST*MAPP)/size(g%GW(iGW)%h)
!      ww = sum(g%GW(iGW)%W*MAPP)*g%GW(iGW)%dt/size(g%GW(iGW)%W)
      
      !%{
      !if length(g(gid).GW)>iGW && isfield(g(gid).GW(iGW),'NLC')
      !    g(gid).GW(iGW).NLC(:,:,9)=g(gid).GW(iGW+1).h;
      !    prism(5,MAPP,g.GW(iGW)) % code 5 is for nonlinear solver
      !else
      !    if ~isfield(g(gid).GW(iGW),'NLC') || isempty(g(gid).GW(iGW).NLC)
      !        % not to die at linkGW
      !        g(gid).GW(iGW).NLC = g(gid).GW(1).h*0;
      !        g(gid).GW(iGW).STO = g(gid).GW(1).h*0;
      !    end
      !    prism(3,MAPP,g.GW(iGW))
      !end
      !%}
      if (iGW==1) then
        do i = 1, size(g%GW(1)%Cond)
          g%GW(1)%Cond(i)%K = g%GW(1)%Cond(i)%K + 0.0D0
        enddo
      endif

      !prism(3,MAPP,g.GW(iGW))
!      CASE(3,5) ! Call Groundwater Solver
!        I = 2
!        Mapp = mxGetPr(prhs(I))
!
!        I = 3;
!        pm=mxGetPr(mxGetField(prhs(I),1,'STATE'))
!        GWM = prhs(I)
!        CALL mxCopyPtrToReal8(pm, STATE1, 1)
!        IF (abs(STATE1)<PREC) CALL bridgeGW(prhs(I),1)
!        CondM = mxGetField(prhs(I),1,'Cond')
!        CALL bridgeCond(BBC,CondM,NBC)
!        G => GW(1)
!        CALL Aquifer(G,%val(MAPP),BBC,NBC,mode)
      call Aquifer(g%GW(iGW),MAPP,g%GW(iGW)%Cond,size(g%GW(iGW)%Cond),3)
      
      g%GW(iGW)%tt = g%GW(iGW)%tt + g%GW(iGW)%dt

      if (iGW < size(g%GW)) then
          Cond => g%GW(iGW)%Cond
          idx => Cond(1)%idx
          !%g(gid).VDZ.t = g(gid).VDZ.t + g(gid).VDZ.dtP(3);
          !%g(gid).GW(1).tt = g(gid).GW(1).tt + g(gid).GW(1).dt;
          !Dperc = zeros(size(g(gid).GW(iGW).h));
          Dperc = 0.0D0
          Dperc1D => OneD_ptr(size(Dperc),Dperc)
          gGWh1D1 => OneD_ptr(size(g%GW(iGW)%h),g%GW(iGW)%h)
          gGWh1D2 => OneD_ptr(size(g%GW(iGW+1)%h),g%GW(iGW+1)%h)
          Dperc1D(idx) = (Cond(1)%K * (gGWh1D1(idx) - gGWh1D2(idx)) /Cond(1)%D) * g%GW(iGW)%dt
          g%GW(iGW)%Qperc = g%GW(iGW)%Qperc + Dperc     !% Unit m
      endif

      !% print debug
      !% Dp = meanALL(Dperc.*MAPP);
      !% hh = meanALL((g.GW(1).h-g.GW(1).EB).*g.GW(1).ST.*MAPP);
      !% test = hh - (h0 + ww - Dp);
      !% 4;
!      h1 = sum(g%GW(iGW)%h*g.GW(iGW).ST*MAPP)/size(g%GW(iGW)%h)
!      print*, 3, h1-h0, ww, h1
      end subroutine GW_sol
