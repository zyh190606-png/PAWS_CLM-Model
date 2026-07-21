
!

      subroutine dailyVegUpdate(D)
      use Vdata
      implicit none
      integer :: D
      !global g w
      !%[Y,M,D]=datevec(w(wid).t);
      !gid = w(wid).g;
      !% This function is a simplified vegetation growth module
      !% The vegetation are given prescribed LAI,Root,hc and Kc depending on the
      !% day of the year.
      character*32, dimension(4) :: ITEMS
      character*32 :: it
      type(VDB_type) VDB
      integer, DIMENSION(size(g%Veg%RPT,1),size(g%Veg%RPT,2),size(g%Veg%RPT,3)) :: RPT
      real*8, dimension(size(g%Veg%RPT,1),size(g%Veg%RPT,2),size(g%Veg%RPT,3)) :: ff

      real*8, dimension(size(g%Veg%RPT,1),size(g%Veg%RPT,2),size(g%Veg%RPT,3)) :: ho_temp

      real*8, dimension(size(g%Veg%VDB%ho(1:2,:),1),size(g%Veg%VDB%ho(1:2,:),2)) :: hoC
      real*8, dimension(size(hoC,2)) :: ho
      real*8, dimension(size(g%Veg%VDB%GrowthT(:,1))) :: GT
      real*8, dimension(size(g%Veg%VDB%LAICoeff,2)) :: tt
      real*8, dimension(size(g%Veg%VDB%LAICoeff(1:2,:),2)) :: tm
      real*8, dimension(size(tm)) :: tf


      real*8 :: temp
      integer :: i, j, k, nr, nc, ii, jj, kk

      ITEMS = (/'LAI ','Root','hc  ','Kc  '/)
      RPT = g%Veg%RPT
      VDB = g%Veg%VDB

      !%
      !hoC  = VDB.ho(1:2,:); [nr,nc]=size(hoC);
      hoC = VDB%ho(1:2,:)
      nr = size(hoC,1)
      nc = size(hoC,2)
      !ho = zeros(1,nc);
      ho = 0.0D0
      do i = 1, nc
        GT = VDB%GrowthT(:,i)
        if (D<GT(1) .or. D>GT(4)) then
            ho(i)=hoC(1,i)
        elseif (D>=GT(2) .and. D<=GT(3)) then
            ho(i)=hoC(2,i)
        elseif (D<GT(2)) then
            temp = (D-GT(1))/(GT(2)-GT(1))
            ho(i)=hoC(1,i)+temp*(hoC(2,i)-hoC(1,i))
        else
            temp = (D-GT(3))/(GT(4)-GT(3))
            ho(i)=hoC(2,i)-temp*(hoC(2,i)-hoC(1,i))
        endif
      enddo
      do i = 1, size(ITEMS)
        it = ITEMS(i)
        !g(gid).Veg.day.(it) = VDB.(it)(D,:);
        !g(gid).Veg.(it) = g(gid).Veg.day.(it)(RPT);
        select case (trim(it))
          case ('LAI')
            g%Veg%day%LAI = VDB%LAI(D,:)
            forall(ii=1:size(RPT,1),jj=1:size(RPT,2),kk=1:size(RPT,3))
              g%Veg%LAI(ii,jj,kk) = g%Veg%day%LAI(RPT(ii,jj,kk))
            endforall
          case ('Root')
            g%Veg%day%Root = VDB%Root(D,:)
            forall(ii=1:size(RPT,1),jj=1:size(RPT,2),kk=1:size(RPT,3))
              g%Veg%Root(ii,jj,kk) = g%Veg%day%Root(RPT(ii,jj,kk))
            endforall
          case ('hc')
            g%Veg%day%hc = VDB%hc(D,:)
            forall(ii=1:size(RPT,1),jj=1:size(RPT,2),kk=1:size(RPT,3))
              g%Veg%hc(ii,jj,kk) = g%Veg%day%hc(RPT(ii,jj,kk))
            endforall
          case ('Kc')
            g%Veg%day%Kc = VDB%Kc(D,:)
            forall(ii=1:size(RPT,1),jj=1:size(RPT,2),kk=1:size(RPT,3))
              g%Veg%Kc(ii,jj,kk) = g%Veg%day%Kc(RPT(ii,jj,kk))
            endforall
        end select
      enddo
      !% VDB stores a table, RPT on the column, day on the row, its length= 366

      !% Update Canopy Storage Capacity

      !% Should modify K also


      !% update ho--initial abstraction
      !tt = max(VDB.LAICoeff); tt(tt==0)=1;
      tt = maxval(VDB%LAICoeff, dim=1)
      where(tt==0.0D0)
        tt = 1.0D0
      endwhere
      !tm = min(VDB.LAICoeff(1:2,:)); tf = tt - tm;
      tm = minval(VDB%LAICoeff(1:2,:), dim=1)
      tf = tt - tm
      !tf(tf==0)=1;
      where(tf==0.0D0)
        tf = 1.0D0
      endwhere
      !ff = ((g%Veg%LAI-tm(RPT))/tf(RPT))**g%VDZ%VParm(6)
      do k = 1, size(RPT,3)
        do j = 1, size(RPT,2)
          do i = 1, size(RPT,1)
            ff(i,j,k) = ((g%Veg%LAI(i,j,k)-tm(RPT(i,j,k)))/tf(RPT(i,j,k)))**g%VDZ%VParm(6)
          enddo
        enddo
      enddo
      !ff(ff<0)=0;
      where(ff<0.0D0)
        ff = 0.0D0
      endwhere
      !%denom = 1-exp(-tt/2); denom(denom==0)=1;
      !% a simple linear function
      do k = 1, size(RPT,3)
        do j = 1, size(RPT,2)
          do i = 1, size(RPT,1)
            ho_temp(i,j,k) = ho(RPT(i,j,k))
          enddo
        enddo
      enddo
      g%OVN%ho = sum(ho_temp*ff*g%Veg%Frac, dim=3)
      g%OVN%hod = g%OVN%ho*2.0D0

      end subroutine dailyVegupdate



!

      subroutine dailyVegUpdate2(D)
      ! this is to work with CLM, since now LAI, htop, hbot, root are maintained by CLM, and kc is not used, we ONLY update ho!!!!!
      use Vdata
      implicit none
      integer :: D
      !global g w
      !%[Y,M,D]=datevec(w(wid).t);
      !gid = w(wid).g;
      !% This function is a simplified vegetation growth module
      !% The vegetation are given prescribed LAI,Root,hc and Kc depending on the
      !% day of the year.
      character*32, dimension(4) :: ITEMS
      character*32 :: it
      !type(VDB_type) VDB
      !integer, DIMENSION(size(g%Veg%RPT,1),size(g%Veg%RPT,2),size(g%Veg%RPT,3)) :: RPT
      type(VDB_type),pointer:: VDB
      integer, dimension(:,:,:), pointer:: RPT=>NULL()
      real*8, dimension(size(g%Veg%RPT,1),size(g%Veg%RPT,2),size(g%Veg%RPT,3)) :: ff

      real*8, dimension(size(g%Veg%RPT,1),size(g%Veg%RPT,2),size(g%Veg%RPT,3)) :: ho_temp

      real*8, dimension(size(g%Veg%VDB%ho(1:2,:),1),size(g%Veg%VDB%ho(1:2,:),2)) :: hoC
      !real*8, dimension(size(hoC,2)) :: ho
      real*8, dimension(size(g%Veg%VDB%GrowthT(:,1))) :: GT
      real*8, dimension(size(g%Veg%VDB%LAICoeff,2)) :: tt
      real*8, dimension(size(g%Veg%VDB%LAICoeff(1:2,:),2)) :: tm
      !real*8, dimension(size(tm)) :: tf
      real*8, dimension(0:23) :: tf,ho ! no time before I leave SF!! hard-coded~!!!!!~
      real*8 :: temp,ho_base
      integer :: i, j, k, nr, nc, ii, jj, kk, t
      tf = (/1,6,6,6,6,6,6,6,6,6,6,6,3,3,3,3,3,4,4,4,4,0,0,0/)
      ho = (/0.001,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.004,0.004,0.004,0.003,0.003,0.003,0.004,0.004,0.004,0.004,0.004,0.004,0.,0.,0.005/)
      ho_base = 0.001D0

      ITEMS = (/'LAI ','Root','hc  ','Kc  '/)
      RPT => g%Veg%RPT
      VDB => g%Veg%VDB

      !%
      !hoC  = VDB.ho(1:2,:); [nr,nc]=size(hoC);
      hoC = VDB%ho(1:2,:)
      nr = size(hoC,1)
      nc = size(hoC,2)
      !ho = zeros(1,nc);

      !% update ho--initial abstraction
      !tt = max(VDB.LAICoeff); tt(tt==0)=1;
!      tt = maxval(VDB%LAICoeff, dim=1)
!      where(tt==0.0D0)
!        tt = 1.0D0
!      endwhere
!      !tm = min(VDB.LAICoeff(1:2,:)); tf = tt - tm;
!      tm = minval(VDB%LAICoeff(1:2,:), dim=1)
!      tf = tt - tm
      !tf(tf==0)=1;
      !ff = ((g%Veg%LAI-tm(RPT))/tf(RPT))**g%VDZ%VParm(6)
      do j = 1, size(RPT,2)
        do i = 1, size(RPT,1)
           if (w%DM%mask(i,j)) then
           g%OVN%ho(i,j) =ho_base
           g%OVN%hod(i,j)=0D0
	
           do k = 1, g%Veg%nt
            t = g%Veg%RPT(i,j,k)
            !ff(i,j,k) = (g%Veg%LAI(i,j,k)/tf(t))**g%VDZ%VParm(6)
            if (tf(t) .eq. 0) then
		tf(t)=1D0
	    endif
            ff(i,j,k) = (g%Veg%LAI(i,j,k)/tf(t))**2D0  !now fixed!!
            ff(i,j,k) = max(min(ff(i,j,k),1D0),0D0)
            if (isnan(ff(i,j,k))) ff(i,j,k) = 0D0
            g%OVN%ho(i,j) = g%OVN%ho(i,j)+ho(t)*ff(i,j,k)*g%Veg%Frac(i,j,k)*g%VDZ%VParm(5)
            !ff(i,j,k) = ((g%Veg%LAI(i,j,k)-tm(RPT(i,j,k)))/tf(RPT(i,j,k)))**g%VDZ%VParm(6)
           enddo
           g%OVN%hod(i,j) = g%OVN%ho(i,j)*2D0
          endif
        enddo
      enddo
      end subroutine dailyVegupdate2
