      MODULE VDB
      ! only data types, procedures and constants are defined in here.
      PUBLIC :: pp,gd
      PUBLIC :: initVBC,readcrop
      
      type pp
        real*8, dimension(:), pointer :: p=>NULL()
        integer,dimension(:), pointer :: j=>NULL()
      end type pp
      
      type gd ! general data type
        ! it comes with a strs array, an integer array, a pointers array, and an integer nf, the number of valid fields
        integer,pointer::m(:)=>NULL()
        integer np,nj                      ! number of associated fields
        type(pp),dimension(:),pointer::p=>NULL() ! for real*8 pointers
        type(pp),dimension(:),pointer::j=>NULL() ! for integer pointers
        character*100,dimension(:),pointer::f=>NULL(),fj=>NULL()
      end type gd
      
      interface operator (.F.) ! pointer to field
      module procedure rptsearch
      end interface
      
      interface operator (.FI.) ! indexing
      module procedure rptsearchi
      end interface
      
      interface operator (.I.) ! pointer to field
      module procedure iptsearch
      end interface
      
      interface operator (.II.) ! indexing
      module procedure iptsearchi
      end interface
      
      type(gd) :: VBC ! This structure is used to store Crop database
      
      REAL*8 VDBCON_PPFD50,VDBCON_EPAR,VDBCON_VKC,VDBCON_CP
      REAL*8 VDBCON_SIGMA,VDBCON_RAD2PAR
      ! eco-physiological constants used:
      DATA VDBCON_EPAR /4.55D0/
      DATA VDBCON_PPFD50/75D0/
      DATA VDBCON_VKC /0.41D0/     ! von Karman
      DATA VDBCON_CP /1.013D0/
      DATA VDBCON_SIGMA /4.903D-9/
      DATA VDBCON_RAD2PAR /0.45D0/ ! (DIM) ratio PAR / SWtotal  */
      
      CONTAINS
      
      function rptsearch(v1,str)
      ! return pointer
      real*8,dimension(:),pointer :: rptsearch
      type(gd),intent(in)::v1
      character(*),intent(in)::str
      real*8,target:: z(1)
      integer I,J
        DO i = 1,v1%np
          IF (TRIM(str) == TRIM(v1%f(i))) THEN
            rptsearch => v1%p(i)%p
            RETURN
          ENDIF
        ENDDO
        z = -99D0
        rptsearch => z
      end function rptsearch
      
      function iptsearch(v1,str)
      ! return pointer (integer)
      integer,dimension(:),pointer :: iptsearch
      type(gd),intent(in)::v1
      character(*),intent(in)::str
      integer I,J
      integer,target:: z(1)
        DO i = 1,v1%nj
          IF (TRIM(str) == TRIM(v1%fj(i))) THEN
            iptsearch => v1%j(i)%j
            RETURN
          ENDIF
        ENDDO
        z = -99
        iptsearch => z
      end function iptsearch
      
      function rptsearchi(v1,str)
      ! return index
      integer :: rptsearchi
      type(gd),intent(in)::v1
      character(*),intent(in)::str
      integer I,J
        DO i = 1,v1%np
          IF (TRIM(str) == TRIM(v1%f(i))) THEN
            rptsearchi = i
            RETURN
          ENDIF
        ENDDO
        rptsearchi = -99
      end function rptsearchi
      
      function iptsearchi(v1,str)
      ! return index (integer)
      integer :: iptsearchi
      type(gd),intent(in)::v1
      character(*),intent(in)::str
      integer I,J
        DO i = 1,v1%nj
          IF (TRIM(str) == TRIM(v1%fj(i))) THEN
            iptsearchi = i
            RETURN
          ENDIF
        ENDDO
        iptsearchi = -99
      end function iptsearchi
      
      subroutine initVBC(mm,VBCD)
      integer::mm(2) ![#plants, #properties]
      real*8,target:: VBCD(mm(1),mm(2))
      CHARACTER*100::FF
      integer id,nf
      
      nf = mm(2)
      IF (ASSOCIATED(VBC%F)) DEALLOCATE(VBC%F)
      ALLOCATE(VBC%F(nf))
      ALLOCATE(VBC%p(nf))
      VBC%np = nf
      
      VBC%F(1:nf) = (/'ICNUM','CPNUM','IDC','BIO_E','HVSTI','BLAI',        &
     &  'FRGRW1','LAIMX1','FRGRW2','LAIMX2','DLAI','CHTMX','RDMX',        &
     &  'T_OPT','T_BASE','CNYLD','CPYLD','PLTNFR1','PLTNFR2','PLTNFR3','PLTPFR1','PLTPFR2','PLTPFR3',  &
     &  'WSYF','USLE_C','GSI','VPDFR','FRGMAX','WAVP','CO2HI','BIOEHI','RSDCO_PL','ALAI_MIN',          &
     &  'BIO_LEAF','MAT_YRS','BMX_TREES','EXT_COEF','vpd2'/) !vpd2 is added
      id = VBC.FI.'FRGMAX'
      DO I = 1, mm(2)
        VBC%p(i)%p => VBCD(:,i) ! points to memory chunk, as simple as this!!!!
      ENDDO
      end subroutine initVBC
      
      subroutine readCrop(file,mm,VBCD)
      real*8 rr(100)
      integer mcrdb,ic,j,id,mm(2),eof
      real*8,target:: VBCD(mm(1),mm(2))
      real*8,dimension(:),pointer::frgmax=>NULL(),vpdfr=>NULL()
      character*(*) file
      vpdfr  => VBCD(:,27)
      frgmax => VBCD(:,28)
      mcrdb  = 300
      open(unit=15,file=TRIM(file))
      do ic = 1, mcrdb
        read (15,*,iostat=eof) rr(1),rr(2),rr(3)
        if (eof < 0) exit
        id = 3
        do j=1,3
        read (15,*,iostat=eof) rr(id+1),rr(id+2),rr(id+3),rr(id+4),     &
     &   rr(id+5),rr(id+6),rr(id+7),rr(id+8),rr(id+9),rr(id+10)
        id = id + 10
        if (eof < 0) exit
        enddo
        read (15,*,iostat=eof) rr(id+1),rr(id+2),rr(id+3),rr(id+4)
        id =id+4
        if (eof < 0) exit 
        
        DO j=1,id
        VBCD(ic,j)=rr(j)
        ENDDO
        
        VBCD(ic,id+1)=max((1D0 - frgmax(ic)) / (vpdfr(ic) - 1D0),0D0) !vpdfr is kPa. this is for vpd2
      end do
      close(15)
      end subroutine readcrop
      
      END MODULE VDB
      
