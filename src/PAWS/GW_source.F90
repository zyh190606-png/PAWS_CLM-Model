!

      subroutine GW_source(iGW)
      !% (E) In this program source terms set elsewhere will be cleaned.-
      !% So if necessary they must be set elsewhere
      !% Input:
      !% iGW: aquifer layer number; g.Riv.RepLidx (The repetitive indices:-
      !% multiple segments in one land cell); r.Riv.fG (river exchange); -
      !% g.VDZ.Dperc (recharge)
      !% Output:
      !% g.GW.h, W, Cond.h
      !global g r
      use Vdata
      implicit none
      integer :: rid, k1, k2
      real*8 :: dt, Ar
      integer, dimension(:), pointer :: Cidx=>NULL(), idx=>NULL(), Rep=>NULL()
      real*8, dimension(:), pointer :: gGW1W1D=>NULL(), gGWh1D=>NULL()
      real*8, dimension(:,:), pointer :: TN=>NULL()
      real*8, dimension(:), allocatable :: Wgc
      real*8, dimension(size(g%VDZ%GWL,1),size(g%VDZ%GWL,2)) :: D
      real*8, dimension(size(g%GW(1)%h,1),size(g%GW(1)%h,2)) :: D2

      integer :: iGW, i, j, jj, k, kk
      integer :: pi,pj,pz,pk,indpump
      real*8 :: pumprate
      type(BC_type),dimension(:),pointer :: Cond=>NULL()
      
      g%GW(iGW)%W = 0.0D0
      gGW1W1D => OneD_ptr(size(g%GW(1)%W),g%GW(1)%W)
      dt = g%GW(iGW)%dt     
     
      
      if (iGW==1) then
          !% This is the layer that is interacting with surface. These sources include
          !% pumping and stream exchange, etc
          
          !% River Exchange
          !% Assume to go to first layer of GW
          Ar = g%DM%d(1) * g%DM%d(2)
          do i = 1, size(g%Riv%ID)
              rid = g%Riv%ID(i)
              Cidx => g%Riv%Cidx(i)%ptr
              idx => g%Riv%Lidx(i)%ptr
              TN => g%Riv%TN(i)%ptr
              !%g(gid).GW(1).W(idx) =g(gid).GW(1).W(idx)+(TN*(r(rid).Riv.fG(Cidx))')/(Ar*dt); BUG!SEE THE SIGN OF SOURCE TERM!!!!!
              if(allocated(Wgc)) deallocate(Wgc)
              allocate(Wgc(size(TN,1)))
              Wgc = matmul( TN, r(rid)%Riv%fG(Cidx) ) / (Ar*dt)
              !g%GW(1)%W(idx) = g%GW(1)%W(idx) - Wgc
              
              
              DO jj=1,size(Wgc)
                ! this loop avoids potential MBE with compilation-specific handling of this vector calculation below
                ! when there are repeated indices
                 gGW1W1D(idx(jj) ) = gGW1W1D(idx(jj) ) - Wgc(jj)
              ENDDO
              
              !gGW1W1D(idx) = gGW1W1D(idx) - Wgc
              !
              !Rep => g%Riv%RepLidx(i)%ptr
              !if (Rep(1)/=0) then
              !    do jj = 1, size(Rep)
              !        k1 = Rep(jj)
              !        k2 = idx(k1)
              !        !%g(gid).GW(1).W(idx) =g(gid).GW(1).W(k2)-Wgc(k1); % A BUG!!!!!!!!!!
              !        !g%GW(1)%W(k2) = g%GW(1)%W(k2) - Wgc(k1)
              !        gGW1W1D(k2) = gGW1W1D(k2) - Wgc(k1)
              !        !% g(gid).OVN.h(k2)=g(gid).OVN.h(k2)-Moc(k1)/Ar;
              !    enddo
              !endif
          enddo
          !%D = g(gid).GW(1).h - g(gid).GW(1).EB; D(D<0) = g(gid).GW(1).DZ(D<0)/2;
          D = g%VDZ%GWL - g%GW(1)%EB
          D2 = g%GW(1)%h - g%GW(1)%EB
          !D(g(gid).GW(1).h<g(gid).GW(1).ES)=D2(g(gid).GW(1).h<g(gid).GW(1).ES); % This is correct too much drainage
          where(g%GW(1)%h < g%GW(1)%ES .AND. g%GW(1)%h - g%VDZ%GWL<1D0 )
            D = D2
          endwhere
          !D(D<0) = g(gid).GW(1).DZ(D<0)/50;
          where(D < 0.0D0)
            !D = g%GW(1)%DZ / 50.0D0 !????
            D = 0D0
          endwhere
          g%GW(1)%D = D
          g%GW(1)%W = g%GW(1)%W + g%VDZ%Dperc / dt    !% THIS MUST BE ADDED IF NOT DONE INSIDE FORTRAN
          
          !kuai: add pumping source term
          pk=size(w%TDATA%PUMPING)
          do k = 1,pk    
              pj=w%tData%pumping(k)%IX
              pi=w%tData%pumping(k)%IY
              pz=w%tData%pumping(k)%layer
              
              indpump=0
              if (size(w%tData%pumping(k)%v)>1) then  
                do kk=1, size(w%tData%pumping(k)%v) - 1
                    if ((w%t>=w%tData%pumping(k)%datenums(kk)) .AND. (w%t<=w%tData%pumping(k)%datenums(kk+1))) then
                        indpump=kk
                    endif            
                enddo
              elseif (size(w%tData%pumping(k)%v) .eq. 1) then 
                if (w%t>=w%tData%pumping(k)%datenums(1)) then
                    indpump=1
                endif  
              else
                indpump=0
              endif
              if (indpump==0) then
                  if  (w%t>= w%tData%pumping(k)%datenums(size(w%tData%pumping(k)%v))) then
                      indpump=size(w%tData%pumping(k)%datenums)
                  endif                  
              endif
              if (pz==IGW) then
                  if (indpump==0) then
                      pumprate=0
                  else
                      pumprate=w%tData%pumping(k)%v(indpump) ! pumprate should be m3/day, postive value. 
                  endif                  
                  g%GW(1)%W(pi,pj) = g%GW(1)%W(pi,pj)  - pumprate/(g%DM%d(1)*g%DM%d(2))
              endif              
          enddo
      else
          !% percolation from upper layers
          g%GW(iGW)%W = g%GW(iGW-1)%Qperc / dt
          !% There may also be pumping terms:
          
          !kuai: add pumping source term
          pk=size(w%TDATA%PUMPING)
          do k = 1,pk
              pj=w%tData%pumping(k)%IX
              pi=w%tData%pumping(k)%IY
              pz=w%tData%pumping(k)%layer
              
              indpump=0
              if (size(w%tData%pumping(k)%v)>1) then
                do kk=1, size(w%tData%pumping(k)%v) - 1
                    if ((w%t>=w%tData%pumping(k)%datenums(kk)) .AND. (w%t<=w%tData%pumping(k)%datenums(kk+1))) then
                        indpump=kk
                    endif            
                enddo
              elseif (size(w%tData%pumping(k)%v) .eq. 1) then
                if (w%t>=w%tData%pumping(k)%datenums(1)) then
                    indpump=1
                endif
              else
                indpump=0
              endif

              if (indpump==0) then
                  if  (w%t>= w%tData%pumping(k)%datenums(size(w%tData%pumping(k)%v))) then
                      indpump=size (w%tData%pumping(k)%datenums)
                  endif                  
              endif
              
              if (pz==IGW) then
                  if (indpump==0) then
                      pumprate=0
                  else
                      pumprate=w%tData%pumping(k)%v(indpump)
                  endif                  
                  g%GW(1)%W(pi,pj) = g%GW(1)%W(pi,pj)  - pumprate/(g%DM%d(1)*g%DM%d(2))
              endif              
          enddo
      endif

      if (iGW < size(g%GW)) then
          Cond => g%GW(iGW)%Cond
          gGWh1D => OneD_ptr(size(g%GW(iGW+1)%h),g%GW(iGW+1)%h)
          do j = 1, size(Cond)
              if (Cond(j)%JType==2) then
                  !g%GW(iGW)%Cond(j)%h = g%GW(iGW+1)%h(Cond(j)%indices)
                  g%GW(iGW)%Cond(j)%h = gGWh1D(Cond(j)%idx)                  
              endif
          enddo          
      endif
      
      if(allocated(Wgc)) deallocate(Wgc)
      end subroutine GW_source
      
      
      subroutine GW_pack3D(nLayer)
      ! packing multiple layers into one 3D problem
      use memoryChunk
      implicit none
      integer nLayer
      real*8,pointer,dimension(:,:,:)::gwH3D
      call tempStorage('gwH3D',gwH3D, nLayer)
      
      
      
      
      end subroutine GW_pack3D
