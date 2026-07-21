      subroutine allocVData
      ! This subroutine allocate data in VData module that has not been linked in from the Mat file
      ! This steop is done once at the beginning
      USE Vdata
      USE memoryChunk, Only: TEMP,TEMPZ,IT,ITZ,TEMPL,TEMPZL,ITL,ITZL,MM
      IMPLICIT NONE
      type(sta_type),pointer:: day
      integer ns, ny, I
      
      
      day => g%Wea%day
      ns = day%ns; ny = day%ny
      allocate(day%x(ns),day%y(ns),day%decl(ny,ns),day%radd(ny,ns),day%radi(ny,ns),day%CO2(ny,ns), &
     &     day%tau(ny,ns),day%vpd(ny,ns),day%snowf(ny,ns),day%pa(ny,ns),day%taumx(ny,ns,2),day%elev(ns),           &
     &     day%lat(ns),day%lon(ns))
!      
!      do i= 1, ns
!        allocate(w%Stations(i)%decl(366*24))
!      enddo
      
      allocate(w%DM%Map(w%m(1),w%m(2)))
      allocate(w%DM%MapL(w%m(1),w%m(2)))
      !allocate(g%VDZ%ICE(w%m(1),w%m(2),w%m(3)))
      allocate(g%Veg%DEW(w%m(1),w%m(2)))
      allocate(g%Veg%CSNEW(w%m(1),w%m(2)))
      !g%VDZ%ICE = 0D0
      g%Veg%DEW = 0D0
      g%Veg%CSNEW = 0D0
      
      ! These have to be in here for initializations:
      WHERE(w%DM%mask)
        w%DM%Map  = 1D0
        w%DM%MapL = 1
        g%VDZ%MAPP = 1D0
      ELSEWHERE
        w%DM%Map  = 0D0
        w%DM%MapL = 0
        g%VDZ%MAPP = 0D0
      ENDWHERE
      
      !CP added 10/23/10
      day => g%Wea%day
      ns = day%ns; ny = day%ny
      do i= 1, ns
        day%elev(i) = w%Stations(i)%XYELEV(3) !decimal degrees
        day%x(i)    = w%Stations(i)%XYELEV(1)
        day%y(i)    = w%Stations(i)%XYELEV(2)
        day%lat(i)  = w%Stations(i)%latlong(1)
        day%lon(i)  = w%Stations(i)%latlong(2)
        day%PA(:,i) = 101.3D0-0.01152D0*day%elev(i)+0.544D-6*day%elev(i)**2
        day%CO2(:,i)= 388D0  !!330D0  !! THIS IS CHANGING!! RIGHT NOW IT IS 388!! (Oct 2010! Chaopeng)
        day%vpd(:,i)= -99999999
        day%decl(:,i) = -9999999
      enddo
      
      
      ! TEMPORARY storage (to avoid stack overflow)
      ! USE memoryChunk in a program
      ! and link to some layers of these memory chunks
      allocate(TEMP(w%m(1),w%m(2),30))
      allocate(TEMPZ(w%m(1),w%m(2),w%m(3),10))
      ALLOCATE(IT(30),ITZ(10)); IT = 0; ITZ = 0
      allocate(TEMPL(w%m(1),w%m(2),30))
      allocate(TEMPZL(w%m(1),w%m(2),w%m(3),10))
      ALLOCATE(ITL(30),ITZL(10)); ITL = 0; ITZL = 0
      mm(1:4) = w%m(1:4)
      
      end subroutine
