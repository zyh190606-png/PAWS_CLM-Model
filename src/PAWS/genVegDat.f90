
!

      subroutine genVegDat(VDB)
      !% VDB is loaded
      !% a classification is also needed to do different things for different
      !% plants
      use Vdata
      implicit none
      type(VDB_type) VDB
      integer :: nr, nc, i, n
      real*8, dimension(366,size(VDB%LAICoeff, 2)) :: LAI, Kc, hc, Root

      real*8, dimension(size(VDB%LAICoeff(:,1))) :: LAIC
      real*8, dimension(size(VDB%KcCoeff(:,1))) :: KcC
      real*8, dimension(size(VDB%hcCoeff(:,1))) :: hcC
      real*8, dimension(size(VDB%RootCoeff(:,1))) :: RootC
      real*8, dimension(size(VDB%GrowthT(:,1))) :: GrowthT
      real*8, dimension(366) :: days
      
      real*8 :: temp
      logical, dimension(366) :: loc

      ![nr,nc]=size(VDB.LAICoeff);
      nr = size(VDB%LAICoeff, 1)
      nc = size(VDB%LAICoeff, 2)
      !LAI=zeros(366,nc); Kc=zeros(366,nc); hc=zeros(366,nc); Root=zeros(366,nc);
      LAI = 0.0D0
      Kc = 0.0D0
      hc = 0.0D0
      Root = 0.0D0
      days = (/(i, i = 1, 366)/)
      do i = 1, nc
        !LAIC= VDB.LAICoeff(:,i); % LAI coefficient [LAImin,LAImax,alpha]
        LAIC = VDB%LAICoeff(:,i)
        !KcC = VDB.KcCoeff(:,i); % [Kcmin,Kcmax]
        KcC = VDB%KcCoeff(:,i)
        !hcC = VDB.hcCoeff(:,i); % [hcmin,hcmax]
        hcC = VDB%hcCoeff(:,i)
        !RootC= VDB.RootCoeff(:,i); % [Rootmin,Rootmax]
        RootC = VDB%RootCoeff(:,i)
        !GrowthT = VDB.GrowthT(:,i); % [Dini,Dinf,Dw,Ddorm]
        GrowthT = VDB%GrowthT(:,i)

        !days =1:366;
        !% Before growth
        !loc=(days <= GrowthT(1));
        where (days <= GrowthT(1))
          LAI(:,i) = LAIC(1)
          Kc(:,i) = KcC(1)
          hc(:,i) = hcC(1)
          Root(:,i) = RootC(1)
        endwhere
        !% Growth Period
        !loc=(days>= GrowthT(1) & days<GrowthT(2));
        do n = 1, size(days)
          if((days(n) >= GrowthT(1)) .and. (days(n) < GrowthT(2))) then
            temp = (days(n)-GrowthT(1)) / (GrowthT(2)-GrowthT(1))
            LAI(n,i) = LAIC(1)+(LAIC(2)-LAIC(1))*temp     !% Growth Period, linear interp
            Kc(n,i) = KcC(1)+(KcC(2)-KcC(1))*(1.0D0-dexp(-0.7D0*LAI(n,i)))
            Root(n,i) = RootC(1)+(RootC(2)-RootC(1))*temp
            hc(n,i) = hcC(1)+(hcC(2)-hcC(1))*temp
          endif
        enddo
        !% Growth Peaked
        !loc=(days>= GrowthT(2) & days<GrowthT(3));
        where ((days >= GrowthT(2)) .and. (days < GrowthT(3)))
          LAI(:,i) = LAIC(2)
          Kc(:,i) = KcC(2)
          Root(:,i) = RootC(2)
          hc(:,i) = hcC(2)
        endwhere
        !% Plant wilting
        !loc=(days>= GrowthT(3) & days<GrowthT(4));
        do n = 1, size(days)
          if((days(n) >= GrowthT(3)) .and. (days(n) < GrowthT(4))) then
            temp = (days(n)-GrowthT(3)) / (GrowthT(4)-GrowthT(3))
            LAI(n,i) = LAIC(2)-(LAIC(2)-LAIC(1))*temp
            Kc(n,i) = KcC(1)+(KcC(2)-KcC(1))*(1.0D0-dexp(-0.7D0*LAI(n,i)))
            hc(n,i) = hcC(2)-(hcC(2)-hcC(1))*temp
            Root(n,i) = RootC(2)-(RootC(2)-RootC(1))*temp
          endif
        enddo
        !% Plant Dormant
        !loc=(days >= GrowthT(end));
        where(days >= GrowthT(size(GrowthT)))
          LAI(:,i) = LAIC(1)
          Kc(:,i) = KcC(1)
          Root(:,i) = RootC(1)
          hc(:,i) = hcC(1)
        endwhere

      enddo

      VDB%LAI = LAI
      VDB%Kc = Kc
      VDB%Root = Root
      VDB%hc = hc

      end subroutine genVegDat
