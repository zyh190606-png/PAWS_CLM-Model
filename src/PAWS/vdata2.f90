!#include <preproc.h>
      ! Basic Library for the Watershed Model
      ! Chaopeng Shen, PhD. Environmental Engineering
      ! Michigan State University
      ! This module defines the data structure and holds the data
      ! It also contains some basic procedures used on the data
      ! Every Grid occupies one var object.
      ! Notes: When we send var(i) as an input argument, it will not be a pointer but data itself
      ! Notes: saved module actually still exist and associated. And we can't delete mex file without exiting matlab
      ! Notes: Cannot DEALLOCATE pointer that points to Matlab variables
      !Thus the memory management strategy should be association--deassociation, retain only one copy, sync with Matlab
      !This means all data of 'var' should have a copy inside Matlab and only use 'NULLIFY' to dissociate upon exiting Fortran
      ! Notes: If some argument is defined OPTIONAL then you can't use it if .NOT. PRESENT()
      ! Notes: A dummy variable is automatically assigned SAVE attribute (it only stores an address)
      ! Notes: If changes have been made in this module, if must be compiled for other programs who uses this module to know it's changed
      !This means it should be opened by the same project in MVSDEV
      ! Notes: Local Variables need to be initialized. Otherwise problem may occur
      ! Notes: REAL*8, DIMENSION(:,:),POINTER :: for pointer, it must have deferred shape
      ! Notes: v(:,I) doesn't work. vector subscript cannot be used for subsection.
      ! Notes: for fixed head BC. we need to modify not only d, but also b!
      ! NOTES: the mxGetField arguments are case sensitive!!!!
      ! Notes: should not modify a passed in constant (like a typed numeric value)
      ! NOTES: a declared pointer shouldn't do ptr = someVariable
      ! THESE RULES MAKES YOUR CODE MORE COMPATIBLE ON DIFFERENT PLATFORMS:
      ! use % instead of . for struct field connector
      ! use ! instead of C for comments, starting from the 7th column
      ! continuation line symbol should be on the 6th column, and use $ in stead of any other thing
      ! the line to be continued should have an & after 72th column
      ! always use mwpointer for mwpointers, including pm, on 64 bit machine this point is of type *8. so if defined as integer it will die!
      ! if the incoming matlab data is 2D, then when you do mxGetDimensions, just get 2 elements. also declare your accepting integer the correct size
      ! when use CALL mxCopyPtrToInteger4(mxGetDimensions(pm),m,n), initialize m to the number of dimensions and only extraction n dimensions
      ! if no problem during -g compiled mode but there is problem with normal mode, suspect the link program forgot an input argument (had this happen before)
      
      ! NEW NOTES:
      ! If A module uses B, by default, a module C that uses A also has access to B items. To avoid potential incorrect access, following protocol needs to be followed:
      ! The USE statement in module general section (not in subroutines) should have the format 'USE B, Only: ' unless using very basic modules and you know what you are doing
      ! We should have PRIVATE statement after implicit none so that anything in A not declared PUBLIC (including, especially, those made accessible in B by USE)-
      !    will not be shared by C
      ! subroutines that need to use a lot of items from another module should USE in its own scope, not

      MODULE locFunc
      PUBLIC :: snrm,sumdz
      CONTAINS

      DOUBLE PRECISION FUNCTION snrm(n,sx,itol)
      IMPLICIT NONE
      integer n
      REAL*8, DIMENSION(n), INTENT(IN) :: sx
      integer, INTENT(IN) :: itol
      !REAL*8 :: snrm !Compute one of two norms for a vector sx, as signaled by itol. Used by linbcg.
      if (itol <= 3) then
      snrm=sqrt(dot_product(sx,sx)) !Vector magnitude norm.
      else
      snrm=maxval(abs(sx)) !Largest component norm.
      end if
      END FUNCTION snrm

      FUNCTION SUMDZ(IZ,E,DZ)
      REAL*8 SUMDZ
      integer iz
      real*8,DIMENSION(iz)::E,DZ
      SUMDZ = (E(1)-E(IZ))+(DZ(IZ)+DZ(1))/2D0 ! E(1) is highest
      END FUNCTION SUMDZ
      End Module LocFunc

      MODULE Vdata
      !SAVE
      PRIVATE
      PUBLIC :: var, NSTRUCT
      PUBLIC :: vars, params
      PUBLIC :: alloVdata,deAlloVdata
      PUBLIC :: reAlloV, alloV, deAlloV
      PUBLIC :: BC_type, ACCUM !, PREC
      PUBLIC :: THOMAS,linbcg,pcg,insertLoc,First,THOMASBC,interfaceH
      PUBLIC :: rtnewt, rtnobd, rtsafe,SPDXU,elem,snrm
      
      PUBLIC :: oneD_Ptr,oneD_PtrI,twoD_Ptr,twoD_ptrI,twoD_ptrF
      PUBLIC :: threeD_ptr,threeD_ptrI,num2str,str2num,timer
      PUBLIC :: oneD_PtrO,oneD_PtrIO,twoD_PtrO,twoD_ptrIO
      PUBLIC :: threeD_ptrO,threeD_ptrIO
      
      PUBLIC :: datevec2, datenum2, GreenAmptMatrix, vanGenuchten
      PUBLIC :: eomday, ymd2jd
      PUBLIC :: gVegDay_type, acc_type, msg_type, Stats_type, Wea_type
      PUBLIC :: SOIL_type,GW2D_type,VEGE_type,STA_type,THERMAL_type, par_type
      PUBLIC :: OV_type,GAmod_type,Topomod_type
      PUBLIC :: Ptr_Array_1D, Ptr_Array_2D
      !PUBLIC :: VD,SOILS,GW,VE,DAY,THERM,OV,GA,Topo
      integer, PARAMETER :: NSTRUCT = 1 ! number of structure initiated
      !Niu Jie added
      public :: VDZ_type, DM_TYPE, Stations_type, robj_type, Rec_type, usgs_type,usgs_gw_type
      public :: VDB_type, tData_type, rRiv_type, Net_type, month_type, gRiv_type, maps_type
      public :: GlobalG, GlobalW, GlobalR, Env_type, sim_type, comp_type
      public :: g, w, r, Env
      LOGICAL, PUBLIC :: DEBUGG = .FALSE.
      integer, parameter, public :: length_of_character = 100
      
      !NNNNNNNNNNNNNNNNNNNNNNNNN
      
      ! TYPES 11111111111111111111111111111111111111111111111111111111111111111111111     
	  type sim_type
	    real*8, dimension(:), allocatable :: t, v
	    real*8, pointer :: sd=>NULL()
	  end type sim_type

      type comp_type
        real*8, dimension(:), allocatable :: obs, sim
        type(sim_type) :: ts(2)
      end type comp_type
		      
      TYPE vars
      ! All fields saved as pointer so they can be dynamically created.
      ! Just remember to deallocate after usage to prevent memory leak
        integer :: id, NY, NX, NZ
        integer,DIMENSION(:,:,:),POINTER::isat=>NULL(),jt=>NULL()
        REAL*8,DIMENSION(:,:,:),POINTER::S=>NULL(),h=>NULL(),phi=>NULL()
        REAL*8,DIMENSION(:,:,:),POINTER::phiS=>NULL(),K=>NULL(),KS=>NULL(),Z=>NULL()
      END TYPE vars

      TYPE params
        REAL*8:: the,thre,he,lam,Ke,eta
        REAL*8:: KSe,phie,phiSe
        REAL*8:: a ! a is lam*eta
      END TYPE params
      
      TYPE par_type
        character*32,pointer :: field=>NULL(), mapfile=>NULL()
        integer, dimension(:), pointer :: mapp=>NULL()
        integer,pointer :: id=>NULL(),imet=>NULL(), code=>NULL()
        real*8,pointer :: v=>NULL(), ul=>NULL() ,ll=>NULL()
      end TYPE par_type

      TYPE BC_type
        integer:: JType, typ, bnd
        integer,DIMENSION(:),POINTER::idx=>NULL(),Cidx=>NULL(),Sidx=>NULL()
        integer,DIMENSION(:),POINTER::Nidx=>NULL(),Eidx=>NULL(),Widx=>NULL()
        REAL*8,DIMENSION(:),POINTER::h=>NULL(),K=>NULL(),D=>NULL(),Q=>NULL()
        REAL*8 :: JIndex
        REAL*8 :: BData,BData2
        REAL*8,DIMENSION(:),POINTER:: BV=>NULL()
      END TYPE BC_type

      TYPE ACCUM_type
        REAL*8,DIMENSION(:,:),POINTER::infil=>NULL(),drn=>NULL(),evap=>NULL()
        REAL*8,DIMENSION(:,:,:),POINTER::wex
      END TYPE ACCUM_type

      TYPE TopoMod_type
        integer ms(2)
        REAL*8,DIMENSION(:,:),POINTER::E=>NULL(),Ex=>NULL(),Ey=>NULL()
        REAL*8,DIMENSION(:,:),POINTER::S0x=>NULL(),S0y=>NULL(),S0=>NULL()
        real*8 :: S0p
      END TYPE TopoMod_type

      TYPE GAMod_type
        integer :: ms(3)
        integer, dimension(:), pointer :: m
        integer, pointer :: nspan
        REAL*8,DIMENSION(:,:,:),POINTER::Kbar=>NULL(),FM=>NULL(),SW=>NULL()
        REAL*8,DIMENSION(:,:),POINTER::WFP=>NULL(),Code=>NULL(),MM=>NULL()
        REAL*8,DIMENSION(:,:),POINTER::hspan=>NULL(),THEW=>NULL()
        REAL*8,DIMENSION(:),POINTER:: h=>NULL(),l=>NULL(),n=>NULL()
        REAL*8,POINTER::dn=>NULL(),dl=>NULL()
        REAL*8,DIMENSION(:),POINTER::GCode=>NULL()
      END TYPE GAmod_type
      !ms,m,nspan,KBar,WFP,FM,BC,C,SW,M,                         &
      !&      h,n,l,dn,dl,hspan

      TYPE SOIL_type
        integer :: nl,nz,nt
        real*8, pointer :: STATE
        REAL*8,DIMENSION(:,:),POINTER::KS=>NULL(),THER=>NULL(),THES=>NULL()
        REAL*8,DIMENSION(:,:),POINTER::alpha=>NULL(),lambda=>NULL()
        REAL*8,DIMENSION(:,:),POINTER::nn=>NULL(),DP=>NULL(),FC=>NULL()
        REAL*8,DIMENSION(:),POINTER::DPUL=>NULL(),DPUU=>NULL()
        integer,DIMENSION(:),POINTER::DPN=>NULL()
      END TYPE SOIL_type

      type msg_type
        real*8 :: TF
      end type msg_type

      TYPE GW2D_type
        integer ny,nx,nG,AType,nc
        REAL*8,DIMENSION(:,:),POINTER::h=>NULL(),EB=>NULL(),ES=>NULL(),E=>NULL()
        REAL*8,DIMENSION(:,:),POINTER::D=>NULL(),K=>NULL(),W=>NULL(),W2=>NULL()
        REAL*8,DIMENSION(:,:),POINTER::ST=>NULL(),T=>NULL(),STO=>NULL(),Qperc=>NULL()
        REAL*8,DIMENSION(:,:),POINTER::hNEW=>NULL(),DZ=>NULL(),Seep=>NULL(),DR=>NULL()
        REAL*8,DIMENSION(:,:),POINTER::A1B=>NULL(),A2=>NULL(),A3=>NULL(),A4=>NULL(),A5=>NULL()
        REAL*8,DIMENSION(:,:),POINTER::topH=>NULL(),topK=>NULL(),TopD=>NULL(),botH=>NULL()
        real*8,DIMENSION(:,:),POINTER::topBound=>NULL(),DLDR=>NULL(),botK=>NULL(),botD=>NULL()
        logical*1,dimension(:,:),pointer::MAPP=>NULL()
        REAL*8,POINTER::tt=>NULL(),dt=>NULL()
        REAL*8,DIMENSION(:),POINTER::dd=>NULL()
        REAL*8,DIMENSION(:,:,:),POINTER:: NLC=>NULL()
        REAL*8,pointer:: STATE
        type(msg_type) :: msg
        type(BC_type),dimension(:),pointer :: Cond
      END TYPE GW2d_type

      type VDB_type
        real*8,dimension(:,:),pointer::LAICoeff=>NULL(),hcCoeff=>NULL(),KcCoeff=>NULL()
        REAL*8,DIMENSION(:,:),POINTER::RootCoeff=>NULL(),GrowthT=>NULL(),ho=>NULL()
        real*8,dimension(:,:),POINTER::LAI=>NULL(),Kc=>NULL(),Root=>NULL(),hc=>NULL()
      end type VDB_type

      type gVegDay_type
        real*8, dimension(:),POINTER::LAI=>NULL(),Kc=>NULL(),Root=>NULL(),hc=>NULL()
      end type gVegDay_type

      TYPE Vege_type
        integer ny,nx,nt
        integer,DIMENSION(:,:,:),POINTER:: RPT=>NULL(),RDFN=>NULL(),flag=>NULL()
        integer,DIMENSION(:,:),POINTER::SMAP=>NULL(),EDFN=>NULL()
        REAL*8,POINTER:: dt=>NULL(),t=>NULL()
        REAL*8,DIMENSION(:,:),POINTER::CS=>NULL(),Cpc=>NULL(),CSNEW=>NULL(),EvapCS=>NULL(),RPET=>NULL(),BSoil=>NULL(),DEW=>NULL()
        REAL*8,DIMENSION(:,:,:),POINTER::DELTA=>NULL(),LAI=>NULL(),FRAC=>NULL(),Root=>NULL(),Kc=>NULL(),hc=>NULL()
        REAL*8,DIMENSION(:,:,:),POINTER::phu=>NULL(),fLAI=>NULL(),biomass=>NULL()
        REAL*8,DIMENSION(:,:),POINTER::PET=>NULL(),ET=>NULL(),TP=>NULL(),h=>NULL(),Evap=>NULL(),EvapG=>NULL(),Gprcp=>NULL()
        REAL*8,DIMENSION(:,:),POINTER::ALBEDO=>NULL(),ASP=>NULL(),SWE=>NULL(),Salb=>NULL(),SNOCOV=>NULL()
        REAL*8,DIMENSION(:,:),POINTER::SNOM=>NULL(),TEMP1=>NULL(),SNOW=>NULL(),NewS=>NULL(),Imp=>NULL(),hI=>NULL()
        REAL*8,DIMENSION(:),POINTER::SNOCoef=>NULL()
        type(VDB_type) VDB
        type(gVegDay_type) day
      END TYPE Vege_type

      TYPE OV_type
        integer nt,ms(2)
        REAL*8,DIMENSION(:),pointer :: FPARM
        REAL*8,DIMENSION(:,:),POINTER::h=>NULL(),SN=>NULL(),SNSAVE=>NULL(),Mann=>NULL(),Qoc=>NULL(),hback=>NULL()
        REAL*8,DIMENSION(:,:),POINTER::ho=>NULL(),hoD=>NULL(),hc=>NULL(),hf=>NULL(),Rf=>NULL(),SATH=>NULL(),OVF=>NULL()
        REAL*8,DIMENSION(:,:),POINTER::ffrac=>NULL(),U=>NULL(),V=>NULL(),S=>NULL(),hd=>NULL(),Ex=>NULL(),prcp=>NULL(),dist=>NULL()
        REAL*8,DIMENSION(:,:,:),POINTER :: dP=>NULL()
        REAL*8,pointer::dt,t
      END TYPE OV_type


      TYPE Thermal_type
        integer ny,nx,nz,nd,np
        REAL*8,DIMENSION(:,:,:),POINTER::a=>NULL(),c=>NULL(),bbase=>NULL(),Temp=>NULL(),KT=>NULL(),CT=>NULL()
        REAL*8,DIMENSION(:,:,:),POINTER::Sstatev=>NULL(),Ssitev=>NULL(),Soutv=>NULL(),tss=>NULL(),tsa=>NULL()
        REAL*8,DIMENSION(:),POINTER::SParm=>NULL()
        REAL*8,DIMENSION(:,:),POINTER::dfc=>NULL(),SMR=>NULL(),SNOP=>NULL(),SNOE=>NULL(),Salb=>NULL(),BSub=>NULL()
        REAL*8,POINTER:: state=>NULL()
      END TYPE Thermal_type

      TYPE VDZ_type
        integer m(3)
        integer,DIMENSION(:,:),POINTER::JType=>NULL(),NZC=>NULL(),WTL=>NULL()
        REAL*8 SAT
        real*8,pointer::STATE
        REAL*8,POINTER::t=>NULL()
        REAL*8,DIMENSION(:),POINTER::dtP=>NULL()
        REAL*8,DIMENSION(:),POINTER::VParm=>NULL()
        REAL*8,DIMENSION(:,:),POINTER::DF=>NULL(),DPerc=>NULL(),dt=>NULL(),MAPP=>NULL(),PRatio=>NULL(),GWL=>NULL(),SMR=>NULL()
        REAL*8,DIMENSION(:,:,:),POINTER::h=>NULL(),THE=>NULL(),K=>NULL(),CI=>NULL(),SRC=>NULL(),DZ=>NULL(),DDZ=>NULL(),EB=>NULL()
        REAL*8,DIMENSION(:,:,:),POINTER::THER=>NULL(),THES=>NULL(),KS=>NULL(),N=>NULL(),LAMBDA=>NULL(),ALPHA=>NULL(),E=>NULL()
        REAL*8,DIMENSION(:,:,:),POINTER::FC=>NULL(),WP=>NULL(),ET=>NULL(),Cond=>NULL(),SAVES=>NULL(),R=>NULL(),Ko=>NULL()
        real*8,DIMENSION(:,:,:),POINTER::KSSAVE=>NULL(),REP=>NULL()
         ! newly added
        real*8,DIMENSION(:,:,:),POINTER::ICE,OM
        real*8,DIMENSION(:,:,:),POINTER::sand,clay
        integer,dimension(:,:),pointer:: soic
        real*8,DIMENSION(:,:),POINTER::gti
        
        type(GAmod_type) :: GA
        type(Thermal_type),pointer :: th=>NULL()
      END TYPE VDZ_type
      !CALL vandam1c(nz,nzl,DF(I,J),Dperc(I,J),h(I,J,:),                  &
      !&         THE(I,J,:),K(I,J,:),CI(I,J,:),S(I,J,:),dt,                   &
      !&         MAXIT,THER(TP),THES(TP),Lambda(TP),alpha(TP),                &
      !&         nn(TP),KS(TP),DZ(I,J,:),DDZ(I,J,:),BBCI,STATE)

      TYPE STA_type ! Hourly Distributed Weather Data (for the day)
        integer ns,ny
        REAL*8,DIMENSION(:,:),POINTER::prcp=>NULL(),tmin=>NULL(),tmax=>NULL(),temp=>NULL(),rad=>NULL(),ret=>NULL()
        REAL*8,DIMENSION(:,:),POINTER::rets=>NULL(),wnd=>NULL(),hmd=>NULL(),hb=>NULL(),hli=>NULL(),cosz=>NULL()
        REAL*8,DIMENSION(:,:),POINTER::pre=>NULL(),rn=>NULL(),dptp=>NULL()
        real*8,dimension(:,:),pointer::RADD,RADI,decl=>NULL(),SNOWF=>NULL()
        ! zenc is the cosine of zenith angle,pre is pressure
        REAL*8,DIMENSION(:,:),POINTER::CO2=>NULL(),TAU=>NULL(),vpd=>NULL(),Pa=>NULL() !CP newly added
        !CP newly added TAUMX. maximum air tranmittance of the site for different wavebands (last index), for Spokas scheme, normally 0.7
        REAL*8,DIMENSION(:,:,:),POINTER :: TAUMX=>NULL() 
        REAL*8,DIMENSION(:),POINTER::ELEV=>NULL(),lat=>NULL(),lon=>NULL(),X=>NULL(),Y=>NULL()
        LOGICAL,DIMENSION(1000):: SNOW ! number of stations must be less than 1000
        REAL*8,DIMENSION(1000)::Tavg
      END TYPE STA_type


      TYPE QTYPE_type
        integer ny,nz,nx,ks
        REAL*8,DIMENSION(:),POINTER::dat=>NULL(),Q=>NULL(),u=>NULL(),w=>NULL(),A=>NULL(),h=>NULL()
        REAL*8 :: dx
      END TYPE QType_type

      ! Niu Jie added
      type DM_TYPE
        real*8,dimension(:),pointer::Nxy=>NULL(),dim=>NULL(),Bindex=>NULL(),origin_idx=>NULL()
        real*8,dimension(:),pointer::origin=>NULL(),x=>NULL(),y=>NULL(),Pbound=>NULL(),Cbound=>NULL(),Ebound=>NULL()
        real*8,dimension(:,:),pointer::Map=>NULL()
        integer,dimension(:,:),pointer::MapL=>NULL()
        logical*1,dimension(:,:),pointer::mask=>NULL()
        real*8,dimension(:,:),pointer::mask2=>NULL()
        real*8,dimension(:),pointer::d=>NULL(),msize=>NULL(),nGho=>NULL()
      end type DM_TYPE

      type Stations_type
        integer,dimension(:),pointer::id=>NULL(),gIndices=>NULL()
        real*8,dimension(:),pointer::ptr=>NULL(),pct=>NULL()
        real*8,dimension(:),pointer::XYElev=>NULL(),LatLong=>NULL(),wPct=>NULL(),dates=>NULL(),datenums=>NULL()  ! lat and long in decimal degrees
        REAL*8,DIMENSION(:),POINTER::vratio=>NULL(),prcp=>NULL(),tmin=>NULL(),tmax=>NULL(),sol=>NULL(),wnd=>NULL(),hmd=>NULL()
        REAL*8,DIMENSION(:),POINTER::dptp=>NULL(),awnd=>NULL(),evap=>NULL(),Rad=>NULL(),tau=>NULL(),rrad=>NULL(),s1=>NULL() ! RAD is in WJ/m2/day
        REAL*8,DIMENSION(:),POINTER::s0=>NULL(),snow=>NULL(),mnrh=>NULL(),mntp=>NULL(),mxrh=>NULL(),ref=>NULL(),DData=>NULL()
        integer, dimension(:,:), pointer :: neighbor=>NULL()
        real*8, dimension(:), pointer :: dat_time=>NULL(), cosz=>NULL()
        !character*32, dimension(:), pointer :: Name
        !CP added 10/23/10
        real*8,dimension(:),pointer::RD=>NULL(),RI=>NULL() ! direct and diffuse radiation, unlike RAD, these are in W/m2 (to reduce conversion needed)
        real*8,dimension(:),pointer::decl=>NULL()  ! solar declination angle (in radians) (supposed to be the same for all stations)
      end type Stations_type

      type robj_type
        character, pointer :: type=>NULL()
        integer,pointer :: id,sidx,objnum
        real*8,pointer :: dist,ddist,accT,cc,weight
        real*8, dimension(:), pointer :: acc=>NULL(),xy=>NULL(),v=>NULL()
        integer,dimension(:),pointer :: ix=>NULL()
        character*32, dimension(:), pointer :: Field=>NULL()
      end type robj_type

      type acc_type
        real*8, dimension(:), POINTER :: Qoc=>NULL(), Qout=>NULL(), trib=>NULL()
      end type acc_type

      type Rec_type
        type(robj_type), dimension(:), POINTER:: robj=>NULL()
        type(acc_type) :: acc
        real*8, dimension(:,:), POINTER :: hG=>NULL(), hSave=>NULL()
        real*8, dimension(:), POINTER :: Mem=>NULL(), MemR=>NULL(), MemQgc=>NULL(), tss=>NULL()
      end type Rec_type

      type Stats_type
        real*8, pointer :: THES, THE2, THED, TA, TAV, THE_FULL
        real*8, dimension(:), pointer :: tGprcp=>NULL()
      end type Stats_type

      type usgs_type
        real*8, dimension(:), pointer :: dates=>NULL(),v=>NULL(),datenums=>NULL(),t=>NULL()
        real*8, pointer :: id=>NULL(),sd=>NULL(),ed=>NULL()
        character*32, dimension(:), pointer :: name=>NULL()
      end type usgs_type
      
      type usgs_gw_type
        real*8,dimension(:),pointer::dates=>NULL(),dtwt=>NULL(),datenums=>NULL(),v=>NULL(),t=>NULL(),XY=>NULL(),indices=>NULL()
        real*8,pointer::elev=>NULL(),sd=>NULL(),ed=>NULL()
      end type usgs_gw_type
      
      type m3_type
        real*8, dimension(:,:),pointer :: R=>NULL(),E=>NULL()
      end type m3_type
      
      type maps_type
        real*8,dimension(:,:),pointer::Prcp=>NULL(),Inf=>NULL(),PET=>NULL(),ET=>NULL(),EvapG=>NULL(),Dperc=>NULL(),Dperc2=>NULL()
        real*8,dimension(:,:),pointer::SatE=>NULL(),InfE=>NULL(),botT=>NULL(),Qoc=>NULL(),Qgc=>NULL(),SUBM=>NULL(),SMELT=>NULL()
        real*8,dimension(:,:),pointer::GW=>NULL(),ovnH=>NULL()
        real*8,pointer::n=>NULL()
      end type maps_type
      
      type tData_type
        type(usgs_type),dimension(:),POINTER :: usgs=>NULL()
        type(usgs_gw_type),dimension(:),pointer :: usgs_gw=>NULL()
        type(m3_type) :: m3
        type(par_type), dimension(:), pointer :: par=>NULL()
        type(maps_type) :: maps
        real*8, dimension(:,:), pointer :: dc=>NULL(), h1save=>NULL(), h2save=>NULL()
        real*8, pointer :: objfun=>NULL()
      end type tData_type

      type rRiv_type
        real*8,dimension(:),pointer::h=>NULL(),hx=>NULL(),Ac=>NULL(),U=>NULL(),w=>NULL(),wF=>NULL()
        REAL*8,DIMENSION(:),POINTER::R=>NULL(),E=>NULL(),Eta=>NULL(),Ex=>NULL(),ZBank=>NULL(),Qx=>NULL(),EstE=>NULL()
        real*8,dimension(:),pointer::Mann=>NULL(),S=>NULL(),SE=>NULL(),S0=>NULL(),SM=>NULL(),DSV=>NULL(),K=>NULL()
        real*8,dimension(:),pointer::hg=>NULL(),mr=>NULL(),fG=>NULL(),Qoc=>NULL(),Qgc=>NULL(),Ql=>NULL(),prcp=>NULL()
        REAL*8,DIMENSION(:),POINTER::evap=>NULL(),trib=>NULL(),Cond=>NULL(),Ab=>NULL()
        real*8,dimension(:),pointer::rtr=>NULL(),X=>NULL(),Y=>NULL(),DDist=>NULL()
        real*8,pointer::UPST=>NULL(),UPSV=>NULL(),DST=>NULL(),Dx=>NULL(),code=>NULL(),l=>NULL(),t=>NULL(),dt=>NULL(),umax=>NULL()
        integer,dimension(:),pointer:: CType=>NULL()
        integer,pointer::State=>NULL()
      end type rRiv_type

      type Ptr_Array_1D
        real*8, dimension(:), pointer :: ptr=>NULL()
      end type Ptr_Array_1D

      type Ptr_IntArray_1D
        integer, dimension(:), pointer :: ptr=>NULL()
      end type Ptr_IntArray_1D

      type Ptr_Array_2D
        real*8, dimension(:,:), pointer :: ptr=>NULL()
      end type Ptr_Array_2D

      type gRiv_type
        integer,dimension(:),pointer:: ID=>NULL()
        real*8, dimension(:), pointer :: sides=>NULL()
        type(Ptr_Array_1D),dimension(:),POINTER::ZBank=>NULL(),Len=>NULL()
        type(Ptr_IntArray_1D),dimension(:),pointer::Cidx=>NULL(),Lidx=>NULL(),RepLidx=>NULL()
        TYPE(PTR_ARRAY_1D),DIMENSION(:),POINTER::Stages=>NULL(),W=>NULL(),K=>NULL(),Qoc=>NULL()
        type(Ptr_Array_2D),dimension(:),POINTER::TM=>NULL(),TN=>NULL()
      end type gRiv_type

      type Net_type
        integer,dimension(:),pointer::level=>NULL(),UPSC=>NULL(),INSERT=>NULL(),Parent=>NULL(),Trib=>NULL(),TribC=>NULL()
        real*8,pointer::DSQ=>NULL(),UPST=>NULL(),UPS=>NULL(),DS=>NULL(),DSC=>NULL(),UPSQ=>NULL()
      end type Net_type

      type month_type
        real*8, pointer :: prcp, temp, tmx, tmn, rad, wnd, hmd, dew
      end type month_type

      type Wea_type
        type(month_type) month
        type(STA_type) :: day
        real*8,pointer :: dt, t
      end type Wea_type
      
      type Tol_type
        real*8, pointer :: COU=>NULL(), Qmax=>NULL()
      end type Tol_type
      
      type ET_type
        real*8,pointer :: dt, t
      end type ET_type
      
      type GlobalG
        type(DM_TYPE) :: DM
        type(OV_type) :: OVN
        type(VDZ_type) :: VDZ
        type(Vege_type) :: VEG
        type(SOIL_type):: SOILS
        type(GW2D_type), DIMENSION(:),POINTER :: GW=>NULL()
        type(Topomod_type) :: Topo
        type(Wea_type) Wea
        type(gRiv_type) :: Riv
        type(ET_type) :: ET
        real*8,pointer :: t, dt, status
      end type GlobalG

      type GlobalR
        type(rRiv_type) Riv
        type(Net_type) Net
        type(DM_TYPE) :: DM
        real*8, pointer :: status=>NULL()
      end type GlobalR

      type GlobalW
        type(Topomod_type) :: Topo
        type(Vege_type) :: VEG
        type(SOIL_type) :: SOILS
        type(Stations_type), DIMENSION(:),POINTER :: Stations=>NULL()
        type(Rec_type) :: Rec
        type(tData_type) :: tData
        type(Stats_type) :: Stats
        type(DM_TYPE) DM
        type(Tol_type) Tolerances
        integer, dimension(:), pointer :: CR=>NULL(), ER=>NULL()
        real*8,pointer :: t, dt, ModelStartTime, ModelEndTime, plotT, g
        real*8 timezone !CP newly added
        integer m(10)
        character*32 msg
      end type GlobalW

      type Env_type
        character*32, dimension(2) :: recfile
        integer, dimension(2) :: recfid
        real*8 :: KK, t1, t2, t3
      end type Env_type

      !NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
      ! Data 2222222222222222222222222222222222222222222222222222222222222222222222
      !REAL*8, PARAMETER :: PREC = 1.0E-15
      ! var is an assumed shape of TYPE vars. This is to provide possibility of running several grids in the
      ! same Fortran Session using the same module.
      ! when data comes in, it needs to occupy one var using 'alloV' and
      TYPE(vars), DIMENSION(:), POINTER :: var=>NULL()
      !TYPE(params), DIMENSION(:), POINTER :: par
      !TYPE(VDZ),DIMENSION(NSTRUCT), target:: VD
      !TYPE(SOIL),target :: SOILS
      !TYPE(GW2D),DIMENSION(NSTRUCT), target :: GW
      !TYPE(STA),target:: DAY
      !TYPE(Vege),DIMENSION(NSTRUCT), target::VE
      !TYPE(THERMAL),DIMENSION(NSTRUCT), target::THERM
      !TYPE(OV),DIMENSION(NSTRUCT), target::OVN
      !TYPE(GAmod), target::GA
      !TYPE(Topomod), target::Topo
      type(GlobalG),save,target :: g
      type(GlobalR),dimension(:),pointer :: r=>NULL()
      type(GlobalW), save, target :: w
      type(Env_type) :: Env



!====================================================================================
      interface vanGenuchten
        module procedure vanGenuchten3D1
        module procedure vanGenuchten3D2
      end interface
      
      
      interface num2str
        module procedure num2str
        module procedure num2strI
      end interface
!=================================================================================
      ! PROCEDURES 3333333333333333333333333333333333333333333333333333333333333333
      CONTAINS
        SUBROUTINE alloV(id) ! Allocate an var object
        integer :: id, n, i
        IF (ASSOCIATED(var)) THEN
          n=size(var)
          IF (id <=n) THEN ! free previous memory to avoid memory leak
          ! Will not DEALLOCATE Pointers Interfacing Matlab
            CALL deAlloV(id)
          ELSE
            var => reAlloV(var,id)
            DO i=id+1,n
              var(i)%id=i
              NULLIFY(var(i)%h,var(i)%phi,var(i)%phiS)
              NULLIFY(var(i)%isat,var(i)%K,var(i)%KS)
            ENDDO
          ENDIF
        ELSE
          ALLOCATE(var(id))
          DO i=1,id
            var(i)%id=i
            NULLIFY(var(i)%h,var(i)%phi,var(i)%phiS)
            NULLIFY(var(i)%isat,var(i)%K,var(i)%KS)
          ENDDO
        ENDIF
        END SUBROUTINE alloV

        ! Note the difference between alloV and alloVdata
        SUBROUTINE alloVdata(id,ny,nx,nz) ! Allocate data storage memory inside an var object
        integer id, ny, nx, nz
        CALL alloV(id) ! First construct the object
        ALLOCATE(var(id)%h(ny,nx,nz),var(id)%phi(ny,nx,nz))
        ALLOCATE(var(id)%phiS(ny,nx,nz),var(id)%K(ny,nx,nz))

        ALLOCATE(var(id)%KS(ny,nx,nz))
        ALLOCATE(var(id)%isat(ny,nx,nz))
        END SUBROUTINE alloVdata

        SUBROUTINE deAlloV(id)
        integer,OPTIONAL:: id
        integer :: n,i
        IF(PRESENT(id)) THEN
          NULLIFY(var(id)%isat,var(id)%phi,var(id)%h)
          NULLIFY(var(id)%KS,var(id)%K,var(id)%phiS)
        ELSE
          n = size(var)
          DO i = 1,n
          ! This is only to dissociate, not DEALLOCATE since you can't DEALLOCATE pointers for Matlab interfacing
            NULLIFY(var(i)%isat,var(i)%phi,var(i)%h)
            NULLIFY(var(i)%KS,var(i)%K,var(i)%phiS)
          ENDDO
          DEALLOCATE(var)
        ENDIF
        ! Deallocation order reversed to remove memory cleanly
        END SUBROUTINE deAlloV

        SUBROUTINE deAlloVdata(id)
        integer id
        DEALLOCATE(var(id)%isat,var(id)%phi,var(id)%h)
        DEALLOCATE(var(id)%KS,var(id)%K,var(id)%phiS)
        ! Deallocation order reversed to remove memory cleanly
        END SUBROUTINE deAlloVdata

        FUNCTION reAlloV(p,m)
        TYPE(vars), DIMENSION(:), POINTER:: p, reAlloV
        integer :: ierr, nold, m
        ALLOCATE(realloV(m),stat=ierr)
        !IF(ierr /= 0) !call nrerror('problem attempting to allocate memory for var(m)')
        IF(.not. associated(p)) RETURN
        nold = size(p)
        reAlloV(1:min(nold,m)) = p(1:min(nold,m))
        DEALLOCATE(p)
        END FUNCTION reAlloV


       SUBROUTINE THOMAS(N,N1,NM,a,b,c,d,x)
         ! This program does not care about a(1) and c(N)
         ! All BC should be embedded into coefficients
         ! this programs solves for tri-diagonal matrix in the form of
         ! aX(i-1)+bX(i)+cX(i+1) =d
         IMPLICIT NONE
         integer :: N, i, N1, NM
         REAL*8,DIMENSION(N):: a,b,c,d,x,cc,dd
         REAL*8 DEN
         cc(N1) = c(N1)/b(N1)
         dd(N1) = d(N1)/b(N1)

         DO i=N1+1,NM
           DEN = b(i) - cc(i-1) * a(i)
           cc(i) = c(i)/DEN
           dd(i) = (d(i) - a(i) * dd(i-1)) / DEN
         ENDDO

         IF (DEN > 1.D-15) THEN
         x(NM) = dd(NM)
         DO i=NM-1,N1,-1
           x(i)= dd(i) - cc(i)*x(i+1)
         ENDDO
         ENDIF ! else equation is over-determined
       END SUBROUTINE THOMAS

      SUBROUTINE elem(n,A,AO)
      integer n
      REAL*8,DIMENSION(n)::A,AO
      AO = A
      !A(26)=10000
      END SUBROUTINE


      !-------------------------------------------------------------------------
      !This function computes the first spatial derivatives using compact schemes:
      !Originally Written by Phanikumar Mantha
      !Modified and added Dirichlet BC by CP SHEN
      !--------------------------------------------------------------------------
	  SUBROUTINE FIRST(CC,N,DX,L1,A1,LM,AM,F)
	  ! L1 = 1, One sided BC from interior points
	  ! L1 = 0, Neumann, Prescribed First Order Derivative Value F
	  ! Same for LM
	  ! Only pass in interior points!!!!! (1 to N are all interior)
	  IMPLICIT REAL *8 (A-H, O-Z)
	  REAL*8 U,F,A,B,C,D,A1,AM,DX,SI
	  integer N
	  DIMENSION CC(N),F(N),A(N),B(N),C(N),D(N)
	  AI=1.0d0/4.0d0 ! Modified from original FIRST 'coz using the standard Thomas. CHAO PENG
	  BI=1.0d0
	  CI=AI
	  D1I=3.0d0/2.0d0
	  D2I=0.0d0
	  CONST1=D1I/(2.0d0*DX)
	  A = AI
	  B = BI
	  C = CI
	  DO I=2,N-1
	  D(I)=CONST1*(CC(I+1)-CC(I-1))
	  ENDDO
	  SI = 1D0/6D0
	  IF (L1 .eq. 1) THEN
	   A(1)=0.0d0
	   B(1)=1.0d0
	   C(1)=2.0d0
	   D(1)=(-2.5d0*CC(1)+2.0d0*CC(2)+0.5d0*CC(3))/DX
	  ELSEIF (L1 .eq. 2) THEN
	   D(1) = D(1) - .25D0*A1
	  ELSEIF (L1 .eq. 3) THEN ! 4-th order compact BC with Dirichlet BC
	   ![ -1/(6*dx), -3/(2*dx), 3/(2*dx), 1/(6*dx), -1]
       ![       u-1,        u0,       u1,       u2, w1]
       A(1)=0D0
       B(1)=1D0
       C(1)=1D0
       D(1)=(-SI*A1-1.5D0*CC(1)+1.5D0*CC(2)+SI*CC(3))/DX
      ELSEIF (L1 .eq. 4) THEN ! dCCdx(N) = 0
       A(1)=0.0
	   B(1)=1.0
	   C(1)=0.0
	   D(1)=0.0
	  ENDIF

	  IF (LM .eq. 1) THEN
	   A(N)=2D0
	   B(N)=1D0
	   C(N)=0D0
	   D(N)=(-2.5d0*CC(N)+2.0d0*CC(N-1)+0.5d0*CC(N-2))/DX
	  ELSEIF (LM .eq. 2) THEN
	   D(N)=D(N) - 0.25D0*AM
	  ELSEIF (LM .eq. 3) THEN
	   A(N)=1D0
	   B(N)=1D0
	   C(N)=0D0
	   D(N)=(-SI*AM-1.5D0*CC(N)+1.5D0*CC(N-1)+SI*CC(N-2))/DX
	  ELSEIF (LM .eq. 4) THEN ! dCCdx(N) = 0
       A(N)=0.0
	   B(N)=1.0
	   C(N)=0.0
	   D(N)=0.0
      ENDIF

      CALL THOMAS(N,1,N,A,B,C,D,F)
	  !CALL THOMASBC(A,B,C,D,N,F,L1,A1,LM,AM)
	  RETURN
	  END SUBROUTINE

      !--------------------------------------------------------------------------
      !	Thomas algorithm for tridiogonal matrix system with Boundary Conditions Included
      !--------------------------------------------------------------------------
      SUBROUTINE THOMASBC(A,B,C,D,MD,W,L1,A1,LM,AM)
      IMPLICIT REAL *8  (A-H,O-Z)
      DIMENSION A(MD),B(MD),C(MD),D(MD),E(MD),F(MD),W(MD)
      !	L1 =  1:   Dirichlet Boundary Condition; e.g., u=2
      !	L1 =  2:   Neumann Boundary Condition; e.g., du/dx = 0
      !	Notice that if L1=1 is specified, the input values of A1 are ignored.
	  OBC=-1.0D0/C(MD)
 	  IF(L1.EQ.1) THEN
	  E(1)=A(1)/B(1)
	  F(1)=D(1)/B(1)
	  ELSEIF (L1.EQ.2) THEN
	  E(1)=0.0
	  F(1)=A1
	  ENDIF
      MM=MD-1
      DO 1 M=2,MM
      DEN=B(M)-C(M)*E(M-1)
      E(M)=A(M)/DEN
1     F(M)=(D(M)+C(M)*F(M-1))/DEN
	  IF(LM.EQ.1)W(MD)=(OBC*D(MD)-F(MM))/(E(MM)+OBC)   !for SD or FD
      IF(LM.EQ.2) W(MD)=AM        ! i.e., F=0 is specified
      !	WRITE(*,*) 'W(MD)=', W(MD)
      DO 2 MK=1,MM
      M=MD-MK
      W(M)=E(M)*W(M+1)+F(M)
2     CONTINUE
      RETURN
      END SUBROUTINE


      SUBROUTINE linbcg(n,m,A,d,b,x,itol,tol,itmax,iter,err)
      ! Linear BiConjugate Gradient
      ! A is diagonally-stored
      ! Modified from numerical recipies
      ! CP Shen, MSU
      USE locFunc
      IMPLICIT NONE
      integer n,m
      integer d(m),dT(m)
      REAL*8 A(n,m)
      REAL*8, DIMENSION(n), INTENT(IN) :: b !Double precision is a good idea in this
      REAL*8, DIMENSION(n), INTENT(INOUT) :: x !routine.
      integer, INTENT(IN) :: itol,itmax
      REAL*8, INTENT(IN) :: tol
      integer :: iter
      REAL*8, INTENT(OUT) :: err
      REAL*8, PARAMETER :: EPS=1.0e-14
      !REAL*8 :: snrm
      !REAL*8 snrm
      !Solves Ax = b for x, giv en b of the same length, by the iterative biconjugate gradient
      !method. On input x should be set to an initial guess of the solution (or all zeros); itol is
      !1,2,3, or 4, specifying which convergence test is applied (see text); itmax is the maximum
      !number of allowed iterations; and tol is the desired convergence tolerance. On output, x
      !is reset to the improved solution, iter is the number of iterations actually taken, and err
      !is the estimated error. The matrix A is referenced only through the user-supplied routines
      !atimes, which computes the product of either A or its transpose on a vector; and asolve,
      integer :: I,J,K
      REAL*8 :: ak,akden,bk,bkden,bknum,bnrm,dxnrm,xnrm,zm1nrm,znrm
      REAL*8, DIMENSION(size(b)) :: p,pp,r,rr,z,zz


      n=size(b)
      iter=0
      r = 0.0D0; rr=0.0D0;z=0.0D0;zz =0.0D0;p= 0.0D0;pp=0.0D0
      call SPDXU(A,d,x,r,n,m,0) !atimes(x,r,0)
      r=b-r ! residual
      rr=r
      ! call atimes(r,rr,0) !Uncomment this line to get the minimum residual?variant of the algorithm.
      select case(itol)
      case(1)
         bnrm=snrm(n,b,itol)
         call asolve(A,d,r,z,n,m) !asolve(r,z,0)
      case(2)
         call asolve(A,d,b,z,n,m) !asolve(b,z,0) !
         bnrm=snrm(n,z,itol)
         call asolve(A,d,r,z,n,m) !asolve(r,z,0)
      case(3:4)
         call asolve(A,d,b,z,n,m) !asolve(b,z,0)
         bnrm=snrm(n,z,itol)
         call asolve(A,d,r,z,n,m) !asolve(r,z,0)
         znrm=snrm(n,z,itol)
      case default
         !call nrerror(illegal itol in linbcg?
      end select
      do !Main loop.
        if (iter > itmax) exit
        iter=iter+1
        call asolve(A,d,rr,zz,n,m) !asolve(rr,zz,1) !Final 1 indicates use of transpose matrix
        bknum=dot_product(z,rr) !Calculate coefficient bk and direction vectors
        if (iter == 1) then !p and pp.
        p=z
        pp=zz
        else
        bk=bknum/bkden
        p=bk*p+z
        pp=bk*pp+zz
        end if
        bkden=bknum !Calculate coefficient ak, new iterate x, and
        call SPDXU(A,d,p,z,n,m,0) ! atimes(p,z,0) !new residuals r and rr.
        akden=dot_product(z,pp)
        ak=bknum/akden
        call SPDXU(A,d,pp,zz,n,m,1) !atimes(pp,zz,1)
        x=x+ak*p
        r=r-ak*z
        rr=rr-ak*zz
        call asolve(A,d,r,z,n,m) !asolve(r,z,0) !Solve A  = r !and check stopping criterion.
        select case(itol)
        case(1)
        err=snrm(n,r,itol)/bnrm
        case(2)
        err=snrm(n,z,itol)/bnrm
        case(3:4)
        zm1nrm=znrm
        znrm=snrm(n,z,itol)
        if (abs(zm1nrm-znrm) > EPS*znrm) then
        dxnrm=abs(ak)*snrm(n,p,itol)
        err=znrm/abs(zm1nrm-znrm)*dxnrm
        else
        err=znrm/bnrm !Error may not be accurate, so loop again.
        cycle
        end if
        xnrm=snrm(n,x,itol)
        if (err <= 0.5D0*xnrm) then
        err=err/xnrm
        else
        err=znrm/bnrm !Error may not be accurate, so loop again.
        cycle
        end if
        end select
        !write (*,*) ?iter=?iter,?err=?err
        if (err <= tol) exit
      end do
      END SUBROUTINE linbcg

      SUBROUTINE pcg(n,m,A,d,b,x,itol,tol,itmax,iter,err)
      ! This is the pre-conditioned (using main diagonal, or called Jacobi preconditioner)
      !    CG algorithm given by Leveque FD book on page 84, but ak and bk is different
      ! A is a diagonally-stored large,sparse SPD matrix
      ! modified from Numerical Recipies, using Diagonal-Sparse Matrix Storage
      ! CHAOPENG SHEN, Environmental Engineering, Michigan State University
      ! This only works for SPD matrix, but in this case it's 2X faster than linbcg
      USE locFunc
      IMPLICIT NONE
      integer n,m
      integer d(m),dT(m)
      REAL*8 A(n,m)
      REAL*8, DIMENSION(n), INTENT(IN) :: b !Double precision is a good idea in this
      REAL*8, DIMENSION(n), INTENT(INOUT) :: x !routine.
      integer, INTENT(IN) :: itol,itmax
      REAL*8, INTENT(IN) :: tol
      integer :: iter
      REAL*8, INTENT(OUT) :: err
      REAL*8, PARAMETER :: EPS=1.0d-14
      !REAL*8 snrm
      !Solves Ax = b for x, giv en b of the same length, by the iterative biconjugate gradient
      !method. On input x should be set to an initial guess of the solution (or all zeros); itol is
      !1,2,3, or 4, specifying which convergence test is applied (see text); itmax is the maximum
      !number of allowed iterations; and tol is the desired convergence tolerance. On output, x
      !is reset to the improved solution, iter is the number of iterations actually taken, and err
      !is the estimated error. The matrix A is referenced only through the user-supplied routines
      !atimes, which computes the product of either A or its transpose on a vector; and asolve,
      integer :: I,J,K
      REAL*8 :: ak,akden,bk,bkden,bknum,bnrm,dxnrm,xnrm,zm1nrm,znrm
      REAL*8, DIMENSION(size(b)) :: p,r,z

      n=size(b)
      iter=0
      r = 0.0D0; z=0.0D0;p= 0.0D0;
      call SPDXU(A,d,x,r,n,m,0) !atimes(x,r,0)
      r=b-r ! residual
      if (maxval(abs(r))<1.D-12) THEN
            RETURN ! already the answer
      endif
      ! call atimes(r,rr,0) !Uncomment this line to get the minimum residual?variant of the algorithm.
      select case(itol)
      case(1)
         bnrm=snrm(n,b,itol)
         call asolve(A,d,r,z,n,m) !asolve(r,z,0)
      case(2)
         call asolve(A,d,b,z,n,m) !asolve(b,z,0) !
         bnrm=snrm(n,z,itol)
         call asolve(A,d,r,z,n,m) !asolve(r,z,0)
      case(3:4)
         call asolve(A,d,b,z,n,m) !asolve(b,z,0)
         bnrm=snrm(n,z,itol)
         call asolve(A,d,r,z,n,m) !asolve(r,z,0)
         znrm=snrm(n,z,itol)
      case default
         !call nrerror(illegal itol in linbcg?
      end select
      do !Main loop.
        if (iter > itmax) exit
        iter=iter+1
        bknum=dot_product(z,r) !Calculate coefficient bk and direction vectors
        if (iter == 1) then !p and pp.
          p=z
        else
          bk=bknum/bkden
          p=bk*p+z
        end if
        bkden=bknum !Calculate coefficient ak, new iterate x, and
        call SPDXU(A,d,p,z,n,m,0) ! atimes(p,z,0) !new residuals r and rr.
        akden=dot_product(z,p)
        ak=bknum/akden
        x=x+ak*p
        r=r-ak*z
        call asolve(A,d,r,z,n,m) !asolve(r,z,0) !Solve A  = r !and check stopping criterion.
        select case(itol)
        case(1)
          err=snrm(n,r,itol)/bnrm
          err=max(err,snrm(n,r,5)/bnrm)
        case(2)
          err=snrm(n,z,itol)/bnrm
        case(3:4)
          zm1nrm=znrm
          znrm=snrm(n,z,itol)
        if (abs(zm1nrm-znrm) > EPS*znrm) then
          dxnrm=abs(ak)*snrm(n,p,itol)
          err=znrm/abs(zm1nrm-znrm)*dxnrm
        else
        err=znrm/bnrm !Error may not be accurate, so loop again.
        cycle
        end if
        xnrm=snrm(n,x,itol)
        if (err <=  1D-12) THEN
        exit
        ENDIF

        if (err <= 0.5D0*xnrm) then
        err=err/xnrm
        else
        err=znrm/bnrm !Error may not be accurate, so loop again.
        cycle
        end if
        end select
        !write (*,*) ?iter=?iter,?err=?err
        if (err <= tol) exit
      end do
      END SUBROUTINE pcg


      SUBROUTINE asolve(A,d,b,x,n,m)
      ! This is used mainly for pre-conditioning
      ! To solve Main_Diagonal(A)*u=b
      ! Since it's only main diagonal, transpose operator is the same
      integer n,m
      integer d(m)
      REAL*8 A(n,m),x(n),b(n)
      integer I,J
      DO I = 1,m
        IF (d(I) .eq. 0) THEN
          DO J = 1,n
            if (abs(A(J,I)) > 1.d-10) THEN
            x(J)=b(J)/A(J,I)
            else
            x(J)=0.0D0
            endif
          ENDDO
          EXIT
        ENDIF
      ENDDO

      END SUBROUTINE
      SUBROUTINE SPDXU(A,d,u,b,n,m,flag)
      ! Compute b = full(A) * u
      ! A(:,i) is the d(i)-th diagonal of the full matrix
      ! if d(i)<0, then the first abs(d(i)) elements are zero
      ! if d(i)>0, then the last d(i) elements are zero
      ! it is expect u(K) can go out of bound. If no error is reported, let it be
      ! flag>=1, b = full(A)' *b
      integer m,n ! n is the length of u, m is the number of diagonals
      REAL*8 :: A(n,m),u(n),b(n)
      integer d(m),dt(m),tflag ! d stores the location of diagonals
      integer I,J,K(m),flag
      LOGICAL MK(m)

      REAL*8 :: x(0:n+1)
      integer MN,MX
      b = 0.0D0
      x = 0.0D0
      if (flag <1) then
        MN = abs(minval(d))
        MX = abs(maxval(d))
        DO I = 1,min(MN,n)
          K = d + I
#if (defined CBound)           
          ! when checkbound is desired for debugging purpose. This is a compilation-time preprocessor
          DO J = 1,m
           IF (K(J)>=1 .AND. K(J)<=n) THEN
            b(I) = b(I) + A(I,J)*u(K(J))
           ENDIF
          ENDDO
#else
          b(I)=sum(A(I,:)*u(K),mask=(K>=1 .AND. K<=n)) ! this may be faster
#endif
        ENDDO

        DO I = MN+1,n-MX
          K = d + I
          b(I)=sum(A(I,:)*u(K)) ! out of bind should not happen here
        ENDDO

        DO I = n-MX+1,n
          K = d + I
#if (defined CBound) 
          DO J = 1,m
           IF (K(J)>=1 .AND. K(J)<=n) THEN
            b(I) = b(I) + A(I,J)*u(K(J))
           ENDIF
          ENDDO
#else
          b(I)=sum(A(I,:)*u(K),mask=(K>=1 .AND. K<=n))
#endif
        ENDDO
        !DO I = 1,n
          !K = d + I
          !b(I)=sum(A(I,:)*u(K))
        !ENDDO
      else
        !MN = abs(minval(d))
        !MX = abs(maxval(d))
        MN = abs(maxval(d))
        MX = abs(minval(d))
        x(1:n) = u
        DO I = 1,min(MN,n)
          K = -d + I
          !b(I) = b(I) + sum(u(K)*A(K,J))
#if (defined CBound)  
          DO J = 1,m
            IF (K(J)>=1 .and. K(J)<=n) THEN
              b(I) = b(I) + x(K(J))*A(K(J),J)
            ENDIF
          ENDDO
#else
          WHERE (K < 0 .OR. K > n) K = 0
          DO J = 1,m
              b(I) = b(I) + x(K(J))*A(K(J),J)
          ENDDO
#endif
        ENDDO

        DO I = MN+1,n-MX
          K = -d + I
          !b(I) = b(I) + sum(u(K)*A(K,J)) ! wrong
          DO J = 1,m
            b(I) = b(I) + u(K(J))*A(K(J),J) ! Out of bound should not happen here
          ENDDO
        ENDDO

        DO I = n-MX+1,n
          K = -d + I
#if (defined CBound) 
          DO J = 1,m
            IF (K(J)>=1 .and. K(J)<=n) then
            b(I) = b(I) + x(K(J))*A(K(J),J) 
            ENDIF
          ENDDO
#else
          WHERE (K < 0 .OR. K > n) K = 0
          DO J = 1,m
            b(I) = b(I) + x(K(J))*A(K(J),J) 
          ENDDO
#endif        
          !b(I) = b(I) + sum(u(K)*A(K,J)) ! wrong
        ENDDO
      endif
      END subroutine

       FUNCTION rtnewt(funcd,x1,x2,xacc,INP1,INP2,INP3,INP4,INP5,INP6,E)
       ! newton iteration, modified from numerical recipies
       IMPLICIT NONE
       integer JMAX
       REAL*8 :: rtnewt, x1,x2,xacc,INP1,INP2,INP3,INP4,INP5,INP6
       PARAMETER (JMAX =20)
       REAL*8 :: PREC
       PARAMETER (PREC = 1e-12)
       integer :: j
       REAL*8 df,dx,f,E, rt0
       EXTERNAL :: funcd
       ! assuming first try will not generate NaN!!!! (to save some time)
       rtnewt = 0.5D0*(x1 + x2)
       rt0 = rtnewt
       DO j = 1,JMAX
         CALL funcd(rtnewt,f,df,INP1,INP2,INP3,INP4,INP5,INP6)

         if (isnan(f)) then
         ! jumped too much, assuming last rtnewt worked
          rtnewt = rt0 - dx/2D0
          dx = dx/2D0
         else
          dx = f/df
          rt0 = rtnewt
          rtnewt = rtnewt - dx
         endif
         if ((abs(dx)<xacc) .OR. (abs(f)<PREC)) THEN
           E = f
           return
         endif
       ENDDO
       E = f
       !pause 'rtnewt exceeded maximum iterations'
       END FUNCTION rtnewt

       FUNCTION rtnobd(funcd,x0,xacc,INP1,INP2,INP3,INP4,INP5,INP6)
       IMPLICIT NONE
       integer JMAX
       REAL*8 :: rtnobd,x0,xacc,INP1,INP2,INP3,INP4,INP5,INP6
       PARAMETER (JMAX =20)
       REAL*8 :: PREC
       PARAMETER (PREC = 1e-12)
       integer :: j
       REAL*8 df,dx,f
       EXTERNAL :: funcd
       rtnobd = x0
       DO j = 1,JMAX
         CALL funcd(rtnobd,f,df,INP1,INP2,INP3,INP4,INP5,INP6)
         dx = f/df
         rtnobd = rtnobd - dx
         if ((abs(dx)<xacc) .OR. (abs(f)<PREC)) return
       ENDDO
       !PAUSE 'rtnobd exceeded maximum iterations'
       END FUNCTION rtnobd

       FUNCTION rtsafe(funcd,x1,x2,xacc,INP1,INP2,INP3,INP4,INP5,INP6)
       ! newton iteration, modified from numerical recipies
       IMPLICIT NONE
       integer MAXIT
       PARAMETER (MAXIT = 50)
       REAL*8 :: rtsafe, x1,x2,xmax,INP1,INP2,INP3,INP4,INP5,INP6,xacc
       REAL*8 :: PREC
       PARAMETER (PREC = 1e-12)
       integer :: j
       EXTERNAL :: funcd
       REAL*8 df,dx,dxold,f,fh,fl,temp,xh,xl,dfl,dfh
       CALL funcd(x1,fl,dfl,INP1,INP2,INP3,INP4,INP5,INP6)
       CALL funcd(x2,fh,dfh,INP1,INP2,INP3,INP4,INP5,INP6)
       if ((fl > 0.0D0 .AND. fh > 0.0D0).OR.                                &
     &  (fl < 0.0D0 .AND. fh < 0.0D0)) THEN
         !IF (ABS(INP6) < PREC) THEN
           rtsafe = 1e20
           return
         !ELSE
           !PAUSE 'root out bound-rtsafe'
         !ENDIF
       ENDIF
       if (abs(fl)<PREC) then
         rtsafe = x1
         return
       else if (abs(fh)<PREC) then
         rtsafe = x2
         return
       else if (fl.lt.0) then
         xl=x1
         xh=x2
       else
         xl=x2
         xh=x1
       endif
       !rtsafe = 0.5*(x1+x2)
       rtsafe = 0.5D0*(x1-fl/dfl + x2 - fh/dfh)
       dxold=abs(x2-x1)
       dx=dxold
       CALL funcd(rtsafe,f,df,INP1,INP2,INP3,INP4,INP5,INP6)

       if (isnan(f)) then
         rtsafe = 0.5*(x1+x2)
         CALL funcd(rtsafe,f,df,INP1,INP2,INP3,INP4,INP5,INP6)
       endif

       do j = 1,MAXIT
         if (((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f).ge.0.                       &
     &         .or. abs(2.0D0*f).gt.abs(dxold*df) .or. isnan(f)) then
            dxold=dx
            dx=0.5*(xh-xl)
            rtsafe=xl+dx
            if(xl.eq.rtsafe) return
         else
            dxold=dx
            dx=f/df
            temp=rtsafe
            rtsafe=rtsafe-dx
            if(temp.eq.rtsafe) return
         endif
         if ((abs(dx).lt.xacc) .or.(abs(f)<PREC)) return
         CALL funcd(rtsafe,f,df,INP1,INP2,INP3,INP4,INP5,INP6)
         if (f.lt.0.) then
            xl=rtsafe
         else
            xh=rtsafe
         endif
       enddo
       !PAUSE 'rtsafe exceeding MAXIT'
       return
       END FUNCTION rtsafe

       FUNCTION insertLoc(nz,E,ins,I,mode)
       integer insertLoc,nz,I,mode
       REAL*8 E(nz),ins
       integer In,II
       ! find where ins is in E
       ! In is defined as ins between E(In) and E(In+1), may be equal to E(In)
       ! mode=1, E is in descending order
       IF (mode .eq. 1) THEN
         DO II = I,nz
           IF (E(II)<ins) EXIT
         ENDDO
         In = II
         DO II = In,1,-1
           IF (E(II)>=ins) EXIT
         ENDDO
         In = II
       ENDIF
       insertLoc = In

       END FUNCTION insertLoc

       FUNCTION num2strI(num)
       ! Converting an integer number into a string variable
       ! Limited to 16 digits, use TRIM outside to trim the trailing blanks
       integer I,N
       parameter (N = 16)
       INTEGER,intent(in):: num
       CHARACTER*(N) num2strI
       CHARACTER*(N) str1
       
       write(str1,'(I16)') num
       ! maybe:
       ! write(str1,'(I<ia+ib>)') num
       num2strI = ADJUSTL(str1)
       
       END FUNCTION
       
       FUNCTION num2str(num)
       ! Converting an real number into a string variable
       ! string is written in its default format in 16 digits
       ! Limited to F20.10 digits, use TRIM outside to trim the trailing blanks
       integer I,N
       parameter (N = 20)
       real*8,intent(in):: num
       CHARACTER*(N) num2str
       CHARACTER*(N) str1,str2
       
       write(str1,'(F20.10)') num
       ! maybe:
       ! write(str1,'(I<ia+ib>)') num
       num2str = ADJUSTL(str1)
       str2    = num2str
       ! make it more compact:
       DO I = N,1,-1
         IF (str2(I:I) .eq. ' ' .OR. str2(I:I) .eq. '0' .OR. str2(I:I) .eq. '.') THEN
           num2str(I:I) = ' '
         ELSE
           EXIT
         ENDIF
       ENDDO
       
       END FUNCTION

       FUNCTION str2num(str)
       ! Converting a string to a double variable
       REAL*8 str2num
       CHARACTER*(*) str
       
       read(str,*) str2num
       
       END FUNCTION
!=================================================================================
       ! ########################################################################
       ! ########  Pointer Shape Manipulating Functions #########################
       ! Transforming the rank/shape of data
       ! oneD_ptr:  return 1D pointer to a N-D real*8 matrix
       ! oneD_ptrO: return 1D pointer to a N-D real*8 matrix with a stated beginning index beg (so that the declaration is beg:nn)
       ! oneD_ptrI:  return 1D pointer to a N-D integer matrix
       ! oneD_ptrIO: return 1D pointer to a N-D integer matrix with a stated beginning index beg (so that the declaration is beg:nn)
       ! Similar things for twoD_ptr and threeD_ptr, differences is that the optional beg is for the last dimension
       ! Don't try to write explicit interface for these functions, it will make the compiler check for matching and thus these won't work anymore
       ! Also data type checks (performed by the compiler) is useful for debugging
       
       FUNCTION oneD_Ptr(nn,datND)
       ! returns a 1D pointer to a N-D matrices so that row-based index can be used-
       ! to query values
       INTEGER nn,nz,beg
       REAL*8,Target:: datND(nn)
       REAL*8,DIMENSION(:),POINTER::oneD_Ptr
       !print *, dat2D
       oneD_Ptr => datND
       END FUNCTION oneD_Ptr
       
       FUNCTION oneD_PtrO(nn,datND,beg)
       ! returns a 1D pointer to a N-D matrices so that row-based index can be used-
       ! to query values
       INTEGER nn,nz,beg
       REAL*8,Target:: datND(beg:nn)
       REAL*8,DIMENSION(:),POINTER::oneD_PtrO
       !print *, dat2D
       oneD_PtrO => datND
       END FUNCTION oneD_PtrO
       
       FUNCTION oneD_PtrI(nn,datND)
      ! returns a 1D pointer to a N-D matrices so that row-based index can be used-
      ! to query values
       INTEGER nn,nz,beg
       INTEGER,Target:: datND(nn)
       INTEGER,DIMENSION(:),POINTER::oneD_PtrI
       !print *, dat2D
       oneD_PtrI => datND
       END FUNCTION oneD_PtrI
       
       FUNCTION oneD_PtrIO(nn,datND,beg)
      ! returns a 1D pointer to a N-D matrices so that row-based index can be used-
      ! to query values
       INTEGER nn,nz,beg
       INTEGER,Target:: datND(beg:nn)
       INTEGER,DIMENSION(:),POINTER::oneD_PtrIO
       !print *, dat2D
       oneD_PtrIO => datND
       END FUNCTION oneD_PtrIO
       
      FUNCTION twoD_Ptr(nn,nz,dat3D)
      ! returns a 2D pointer to a 3D matrices to condense the first two dimensions into 1
      ! this works thanks to the way memory is stored in Fortran
       INTEGER nn,nz,beg
       REAL*8,Target:: dat3D(nn,nz)
       REAL*8,DIMENSION(:,:),POINTER::twoD_Ptr
       !print *, dat2D
       twoD_Ptr => dat3D
       END FUNCTION twoD_Ptr
       
       FUNCTION twoD_PtrO(nn,nz,dat3D,beg)
      ! returns a 2D pointer to a 3D matrices to condense the first two dimensions into 1
      ! this works thanks to the way memory is stored in Fortran
       INTEGER nn,nz,beg
       REAL*8,Target:: dat3D(nn,beg:nz)
       REAL*8,DIMENSION(:,:),POINTER::twoD_PtrO
       !print *, dat2D
       twoD_PtrO => dat3D
       END FUNCTION twoD_PtrO
       
       FUNCTION twoD_PtrF(nn,nz,dat3D)
      ! returns a 2D pointer to a 3D matrices to condense the first two dimensions into 1
      ! this works thanks to the way memory is stored in Fortran
      ! Full dimension descriptions
       INTEGER nn(2),nz(2),beg
       REAL*8,Target:: dat3D(nn(1):nn(2),nz(1):nz(2))
       REAL*8,DIMENSION(:,:),POINTER::twoD_PtrF
       !print *, dat2D
       twoD_PtrF => dat3D
       END FUNCTION twoD_PtrF
       
       FUNCTION twoD_PtrI(nn,nz,dat3D)
      ! returns a 2D pointer to a 3D matrices to condense the first two dimensions into 1
      ! this works thanks to the way memory is stored in Fortran
       INTEGER nn,nz,beg
       INTEGER,Target:: dat3D(nn,nz)
       INTEGER,DIMENSION(:,:),POINTER::twoD_PtrI
       !print *, dat2D
       twoD_PtrI => dat3D
       END FUNCTION twoD_PtrI
       
       FUNCTION twoD_PtrIO(nn,nz,dat3D,beg)
      ! returns a 2D pointer to a 3D matrices to condense the first two dimensions into 1
      ! this works thanks to the way memory is stored in Fortran
       INTEGER nn,nz,beg
       INTEGER,Target:: dat3D(nn,beg:nz)
       INTEGER,DIMENSION(:,:),POINTER::twoD_PtrIO
       !print *, dat2D
       twoD_PtrIO => dat3D
       END FUNCTION twoD_PtrIO
       
       FUNCTION threeD_Ptr(m,dat3D)
      ! returns a 3D pointer to a data chunk
       INTEGER m(3)
       REAL*8,Target:: dat3D(m(1),m(2),m(3))
       REAL*8,DIMENSION(:,:,:),POINTER::threeD_Ptr
       !print *, dat2D
       threeD_Ptr => dat3D
       END FUNCTION threeD_Ptr
       
       FUNCTION threeD_PtrO(m,dat3D,beg)
      ! returns a 3D pointer to a data chunk
       INTEGER m(3),beg
       REAL*8,Target:: dat3D(m(1),m(2),beg:m(3))
       REAL*8,DIMENSION(:,:,:),POINTER::threeD_PtrO
       !print *, dat2D
       threeD_PtrO => dat3D
       END FUNCTION threeD_PtrO
       
       FUNCTION threeD_PtrI(m,dat3D)
      ! returns a 3D pointer to a data chunk
       INTEGER m(3)
       INTEGER,Target:: dat3D(m(1),m(2),m(3))
       INTEGER,DIMENSION(:,:,:),POINTER::threeD_PtrI
       !print *, dat2D
       threeD_PtrI => dat3D
       END FUNCTION threeD_PtrI
       
       
       FUNCTION threeD_PtrIO(m,dat3D,beg)
      ! returns a 3D pointer to a data chunk
       INTEGER m(3),beg
       INTEGER,Target:: dat3D(m(1),m(2),beg:m(3))
       INTEGER,DIMENSION(:,:,:),POINTER::threeD_PtrIO
       !print *, dat2D
       threeD_PtrIO => dat3D
       END FUNCTION threeD_PtrIO

       ! END ARRAY MANIPULATION  #################################################

!=================================================================================
       Subroutine interfaceH(ny,nx,h,E,Ex,Ey,hx,hy,M)
       implicit none
       integer ny, nx
       real*8, DIMENSION(ny,nx):: h, E, Ex,Ey,hx,hy,t1,t2
       logical:: M(ny,nx)
      !!! THIS ONE DOES NOT LOOK AT Ex or Ey, DISREGARD EDGE EFFECT

       integer i,j
       real*8 ER(ny,nx)
       REAL*8 PREC,tt1,tt2,phi

       hx = 0.0D0
       hy = 0.0D0
       ER = 0.0D0
       PREC = 1.0D-10
      ! hx!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       IF (nx>1) THEN
       DO I = 1,ny
         DO J = 1,nx-1
          IF (M(I,J)) THEN
           if (E(I,J)>=E(I,J+1)) THEN ! Positive Wind
             if (h(I,J)<PREC) THEN
               hx(I,J)=0.0D0
             elseif (h(I,J+1)>PREC) THEN
               if (h(I,J)*4D0<h(I,J+1)) THEN ! changed to 10/10/2009
                 hx(I,J)=h(I,J)
               ELSE
                hx(I,J)=0.5D0*(h(I,J)+h(I,J+1)) ! modification 05/06/2009
               endif
             else ! half-wet
               hx(I,J)=h(I,J)
             endif

           else
             if (h(I,J+1)<PREC) THEN
               hx(I,J) = 0.0D0
             elseif (h(I,J)>PREC) THEN
              !hx(I,J) = 0.5D0*(h(I,J)+h(I,J+1))

               if (h(I,J+1)<h(I,J)) THEN
                 hx(I,J)=h(I,J+1)
               ELSE
                 hx(I,J)=0.5D0*(h(I,J)+h(I,J+1)) ! modification 05/06/2009
               endif
             else ! half-wet
               hx(I,J) = h(I,J+1)
             endif

           endif
         ENDIF ! mask if
         ENDDO
       ENDDO
       hx = max(hx,0.0D0)
       ENDIF

       ! hy!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       IF (ny>1) THEN
       DO I = 1,ny-1
         DO J = 1,nx
        !IF ((I .eq. 29) .AND. (J .eq. 45)) THEN
        !print *,3
        !ENDIF

           if (E(I,J)>=E(I+1,J)) THEN
             if (h(I,J)<PREC) THEN
               hy(I,J)=0.0D0
             elseif (h(I+1,J)>PREC) THEN
              !hy(I,J)=0.5D0*(h(I,J)+h(I+1,J))
              !tt1=h(I,J); tt2=h(I+1,J)
               if (h(I,J)<h(I+1,J)) THEN
                 hy(I,J)=h(I,J)
               ELSE
                 hy(I,J)=0.5D0*(h(I,J)+h(I+1,J)) ! modification 05/06/2009
               endif
             else ! half-wet
               hy(I,J)=h(I,J)
             endif
           else
             if (h(I+1,J)<PREC) THEN
               hy(I,J) = 0.0D0
             elseif (h(I,J)>PREC) THEN
              !hy(I,J) = 0.5D0*(h(I,J)+h(I+1,J))

               if (h(I+1,J)<h(I,J)) THEN
                 hy(I,J)=h(I+1,J)
               ELSE
                 hy(I,J)=0.5D0*(h(I,J)+h(I+1,J)) ! modification 05/06/2009
               endif
             else ! half-wet
               hy(I,J) = h(I+1,J)
             endif
           endif
         ENDDO
       ENDDO
       hy = max(hy,0.0D0)
       ENDIF

       END SUBROUTINE interfaceH

      !===============================================================================
      function datenum2(dd)
      ! datenum2([y m d]) or datenum2([y jd]) or datenum2(yyyymmdd)
      ! a faster version of datenum by Chaopeng
      ! datenums of first day of every year from 1950 to 2050:
      integer :: datenum2, datenums(101), d, y
      integer, dimension(:) :: dd
      datenums = [712224,712589,712954,713320,713685,714050,714415,714781,715146,715511, &
        715876,716242,716607,716972,717337,717703,718068,718433,718798,71916,  &
        719529,719894,720259,720625,720990,721355,721720,722086,722451,722816, &
        723181,723547,723912,724277,724642,725008,725373,725738,726103,726469, &
        726834,727199,727564,727930,728295,728660,729025,729391,729756,730121, &
        730486,730852,731217,731582,731947,732313,732678,733043,733408,733774, &
        734139,734504,734869,735235,735600,735965,736330,736696,737061,737426, &
        737791,738157,738522,738887,739252,739618,739983,740348,740713,741079, &
        741444,741809,742174,742540,742905,743270,743635,744001,744366,744731, &
        745096,745462,745827,746192,746557,746923,747288,747653,748018,748384,748749]

        if(size(dd)==1) then
            call ymd2jd(dd, 3, d, y);
        elseif(size(dd)==3) then
            y = dd(1)       !m =  dd(2); d = dd(3);
            call ymd2jd(dd, 3, d, y);
        elseif(size(dd)==2) then
            y = dd(1)       !m =  dd(2); d = dd(3);
            d = dd(2)
        end if
        if((y>=1950).and.(y<=2050)) then
                datenum2 = datenums(y-1949) + d - 1;
        else
            write(*,*) 'time exceeded 1950~2050 which is the range of datenum2'
            !pause
            !stop
        end if
      end function datenum2

      !=============================================================================
      function datevec2(dd)
      ! datevec2(datenum)
      ! a faster version of datevec by Chaopeng
      ! datenums of first day of every year from 1950 to 2050:
      implicit none
      integer :: datevec2(3), dd, datenums(101), k, y, jd, i, YJD
      dd = floor(dd+1e-12);
      datenums = [712224,712589,712954,713320,713685,714050,714415,714781,715146,715511, &
     &      715876,716242,716607,716972,717337,717703,718068,718433,718798,71916, &
     &      719529,719894,720259,720625,720990,721355,721720,722086,722451,722816, &
     &      723181,723547,723912,724277,724642,725008,725373,725738,726103,726469, &
     &      726834,727199,727564,727930,728295,728660,729025,729391,729756,730121, &
     &      730486,730852,731217,731582,731947,732313,732678,733043,733408,733774, &
     &      734139,734504,734869,735235,735600,735965,736330,736696,737061,737426, &
     &      737791,738157,738522,738887,739252,739618,739983,740348,740713,741079, &
     &      741444,741809,742174,742540,742905,743270,743635,744001,744366,744731, &
     &      745096,745462,745827,746192,746557,746923,747288,747653,748018,748384,748749]

      do i = 1, size(datenums)
        if(datenums(i) > dd) then
          k = i
          exit
        end if
      end do
      y = k + 1948
      jd = dd - datenums(k-1) + 1
      call ymd2jd([y, jd], 4, YJD, darray= datevec2)

      end function datevec2

      !-----------------------------------------------------------------------------
      subroutine ymd2jd(dd, option, d, YJD, darray)
      ! option=0, input yyyyjul,output yyyymmdd
      ! option=1, input yyyymmdd, output yyyyjul
      ! option=2, input MATLAB julian, output yyyyjul
      ! option=3, input dd=[y m d], output jd
      ! option=4, input dd=[y jd], output [y m d] and yyyymmdd
      ! option=5, input MATLAB datenum, output [y m d]
      implicit none
      integer, dimension(:) :: dd
      integer :: option, d, dd1, ceomd(12)
      integer, optional :: YJD
      integer, dimension(3), optional :: darray
      integer :: y, m, day, i, jd
      integer :: mm(12) = [(i, i = 1, 12)]
      select case(option)
        case(0)
            !year = floor(dd/1000);
            !start = datenum(num2str(year*10000+101),'yyyymmdd');
            !str=num2str(dd);
            !days=str2num(str(end-2:end));
            !day = start+days-1;
            !d=num2str(datestr(day,'yyyymmdd'));
        case(1)
            !year = floor(dd/10000);
            !start = datenum(num2str(year*10000+101),'yyyymmdd');
            !curr = datenum(num2str(dd),'yyyymmdd');
            !YJD = curr - start + 1;
            !d=year*1000+YJD;
        case(2)
            !year = str2num(datestr(dd,'yyyy'));
            !start = datenum(num2str(year*10000+101),'yyyymmdd');
            !YJD = dd - start + 1;
            !d=year*1000+YJD;
        case(3)
            if(size(dd)==3) then
                y = dd(1)
                m = dd(2)
                day = dd(3)
            elseif(size(dd)==1) then
                dd1 = dd(1)
                y = floor(real(dd1/10000.0D0))
                m = floor((dd1 - y*10000.0D0)/100.0D0)
                day = dd1 - y*10000 - m*100
            end if
            call ceomday(y, mm, ceomd)     ! use cumulative end of month day function. Chaopeng
            if(m > 1) then
              d = ceomd(m-1) + day
            else
              d = day
            end if
            if(present(YJD)) YJD = y
        case(4)
            y = dd(1)
            jd = dd(2)
            call ceomday(y, mm, ceomd)
            do i = 1, size(ceomd)
              if(ceomd(i) >= jd) then
                m = i
                exit
              end if
            end do
            if (m > 1) then
              day = jd - ceomd(m-1)
            else
              day = jd
            end if
            if(present(darray)) then
              darray = [y, m, day]
            end if
            if(present(YJD)) YJD = y * 1e4 + m * 1e2 + day
        case(5)
            !ymd = datestr(dd,'yyyymmdd'); d = str2num(ymd);
            !YJD = [str2num(ymd(1:4)), str2num(ymd(5:6)), str2num(ymd(7:8))];
        case default
      end select

      end subroutine ymd2jd

      !--------------------------------------------------------------------------------
      subroutine ceomday(y, m, d)
      ! m should be 1:12
      ! Cumulative EOMDAY End of month.
      ! Chaopeng: fast for determining jd
      implicit none
      integer :: y, m(12), d(12), loc(12)
      integer :: cdpm(12), n = 12, i
      cdpm = [31,59,90,120,151,181,212,243,273,304,334,365]
      !d = y - m
      d = cdpm(m)
      loc = 0
      do i = 1, n
        if((mod(y,4) == 0 .and. mod(y,100) /= 0) .and. m(i)>=2) then
            loc(i) = 1
        end if
      end do
      do i = 1, n
        if(mod(y,400)==0 .and. (m(i)>=2)) then
            loc(i) = 1
        end if
      end do
      do i = 1, n
        if(loc(i)==1) then
            d(i) = d(i) + 1
        end if
      end do


      ! a correction: 2000 is a leap year
      ! However, some exceptions to this rule are required since the duration of
      ! a solar year is slightly less than 365.25 days. Years that are evenly
      ! divisible by 100 are not leap years, unless they are also evenly
      ! divisible by 400, in which case they are leap years.[1][2]
      ! For example, 1600 and 2000 were leap years, but 1700, 1800 and 1900
      ! were not. Similarly, 2100, 2200, 2300, 2500, 2600, 2700, 2900, and
      ! 3000 will not be leap years, but 2400 and 2800 will be.

      end subroutine ceomday



!====================================================================================
      subroutine vanGenuchten3D1(dat,THES,THER,alpha,lambda,n,KS,option,THE,S,K,h)
      implicit none
      integer nr,nc,nz,ii,jj,kk
      real*8,dimension(:,:,:)::dat,THES,THER,alpha,lambda,n,KS,THE,S,K,h
      !real*8,dimension(size(dat,1),size(dat,2),size(dat,3))::m
      !logical,dimension(size(dat,1),size(dat,2),size(dat,3))::sat
      real*8,dimension(:,:,:),allocatable::m
      logical,dimension(:,:,:),allocatable::sat
      integer option

      if(allocated(m)) deallocate(m)
      allocate(m(size(dat,1),size(dat,2),size(dat,3)))
      if(allocated(sat)) deallocate(sat)
      allocate(sat(size(dat,1),size(dat,2),size(dat,3)))

      nr = size(dat,1)
      nc = size(dat,2)
      nz = size(dat,3)
      ![nr,nc]=size(THES);
      if (option==0) then
        !% head
        m = 1.0D0 - 1.0D0/n;
        h = dat
        THE = THER + (THES-THER)*(1.0D0/(1.0D0+abs(alpha*h)**n)**m)
        S = (THE-THER)/(THES-THER)
        K = KS*(S**lambda)*((1.0D0-(1.0D0-S**(1.0D0/m))**m)**2.0D0)
        !sat=(h>0);
        sat = .false.
        where(h>0.0D0)
          sat = .true.
        end where
        !if (nr>1 || nc>1)
        !    THE(sat)=THES(sat);
        !    S(sat)=1;
        !    K(sat)=KS(sat);
        !else
        !    THE(sat)=THES;
        !    S(sat)=1;
        !    K(sat)=KS;
        !end
        where(sat)
          THE = THES
          S = 1.0D0
          K = KS
        end where
      elseif (option==1) then
        !% theta
        THE = dat
        m = 1.0D0 - 1.0D0/n
        S = (THE-THER)/(THES-THER)
        h = ((S**(-1.0D0/m)-1.0D0)**(1.0D0/n))/(-abs(alpha))
        K = KS*(S**lambda)*((1.0D0-(1.0D0-S**(1.0D0/m))**m)**2.0D0)
      end if
      !  if any(K>KS)
      !      4;
      !  end

      if(allocated(m)) deallocate(m)
      if(allocated(sat)) deallocate(sat)

      end subroutine vanGenuchten3D1

!====================================================================================
      subroutine vanGenuchten3D2(dat,THES,THER,alpha,lambda,n,KS,option,THE)
      implicit none
      integer nr,nc,nz,ii,jj,kk
      real*8 :: dat,h
      real*8,dimension(:,:,:)::THES,THER,alpha,lambda,n,KS,THE
      !real*8,dimension(size(n,1),size(n,2),size(n,3))::m,S,K
      real*8,dimension(:,:,:),allocatable::m,S,K
      logical::sat=.false.
      integer option

      if(allocated(m)) deallocate(m)
      allocate(m(size(n,1),size(n,2),size(n,3)))
      if(allocated(S)) deallocate(S)
      allocate(S(size(n,1),size(n,2),size(n,3)))
      if(allocated(K)) deallocate(K)
      allocate(K(size(n,1),size(n,2),size(n,3)))

      nr = size(n,1)
      nc = size(n,2)
      nz = size(n,3)
      ![nr,nc]=size(THES);
      if (option==0) then
        !% head
        m = 1.0D0 - 1.0D0/n;
        h = dat
        THE = THER + (THES-THER)*(1.0D0/(1.0D0+abs(alpha*h)**n)**m)
        S = (THE-THER)/(THES-THER)
        K = KS*(S**lambda)*((1.0D0-(1.0D0-S**(1.0D0/m))**m)**2.0D0)
        !sat=(h>0);
        if(h>0) sat = .true.
        if(sat) then
          THE = THES
          S = 1
          K = KS
        end if
      elseif (option==1) then
        !% theta
        THE = dat
        m = 1.0D0 - 1.0D0/n
        S = (THE-THER)/(THES-THER)
        !h = ((S**(-1.0D0/m)-1.0D0)**(1.0D0/n))/(-abs(alpha))
        K = KS*(S**lambda)*((1.0D0-(1.0D0-S**(1.0D0/m))**m)**2.0D0)
      end if
      !  if any(K>KS)
      !      4;
      !  end
      if(allocated(m)) deallocate(m)
      if(allocated(S)) deallocate(S)
      if(allocated(K)) deallocate(K)

      end subroutine vanGenuchten3D2

!====================================================================================

      subroutine GreenAmptMatrix(opt,m,VDZ,GA,FM,KBar)
      implicit none
      !  % nargin==1 -- construct GA
      !  % nargin> 1 -- (GA,n,l,h)
      !  % hspan: [start,stepsize,end,endingpos]
      !  global g
      integer :: opt, i, j, k, nr, nc, nz
      real*8,optional :: m
      real*8 :: LAMBDA(2), N(2)
      real*8 :: pk(4), dd(4), nk, hspan(3,4)
      real*8, allocatable :: h(:), temp(:), KK(:,:,:), TK(:,:,:), MM(:,:), l(:), nn(:)
      real*8, allocatable :: DZ(:,:,:), sumDZ(:,:,:), KS(:,:,:), THE(:), S(:)
      real*8,optional :: FM(:,:,:), KBar(:,:,:)
      type(VDZ_type),optional :: VDZ
      type(GAMod_type),optional :: GA
      if (opt == 1) then
        !GA = varargin{1};
        !n = varargin{2};
        !l = varargin{3};
        !h = varargin{4};
        !i = round((n - GA.n(1))/GA.dn)+1; i = min(max(i,1),length(GA.n));
        !j = round((l - GA.l(1))/GA.dl)+1; j = min(max(j,1),length(GA.l));
        !if h<GA.hspan(1,1)
        !    gg = 1;
        !    k = 1;
        !elseif h > GA.hspan(end,3)
        !    gg = GA.nspan;
        !    k = GA.m(3);
        !else
        !    for ii = 1:GA.nspan
        !        if ii == 1
        !            sd = 1;
        !        else
        !            sd = GA.hspan(ii-1,4);
        !        end
        !        if h >= GA.hspan(ii,1) && h<GA.hspan(ii,3)
        !            gg = ii;
        !            k = round((h - GA.hspan(ii,1))/GA.hspan(ii,2))+sd;
        !            break
        !        end
        !    end
        !end
        !SW = GA.SW(i,j,end) - GA.SW(i,j,k);
        !varargout{1}=SW;
      elseif (opt == 0) then
        !m = varargin(1)     !if isempty(m),m=50; end
        if(.not. present(m)) then
          m = 50.0D0
        end if
        !VDZ = varargin(2)
        ! constructing the matrix
        ! n is on i and l is on j
        !LAMBDA=[minALL(VDZ.LAMBDA) maxALL(VDZ.LAMBDA)];
        !N=[minALL(VDZ.N) maxALL(VDZ.N)];
        !LAMBDA = [min(minALL(VDZ.LAMBDA),0.4) max(maxALL(VDZ.LAMBDA),0.7)];
        LAMBDA = [min(minval(VDZ%LAMBDA(:,:,:)),0.4), max(maxval(VDZ%LAMBDA(:,:,:)),0.7)]
        !N = [min(minALL(VDZ.N),1.2) max(maxALL(VDZ.N),2.4)];
        N = [min(minval(VDZ%N(:,:,:)),1.2), max(maxval(VDZ%N(:,:,:)),2.4)]
        !l = LAMBDA(1):(LAMBDA(2)-LAMBDA(1))/m:LAMBDA(2); if isempty(l), l = LAMBDA(1); end
        allocate(l(int(m)+1))
        do i = 1, int(m + 1)
          l(i) = LAMBDA(1)+(i-1.0D0)*(LAMBDA(2)-LAMBDA(1))/m
        end do
        if(LAMBDA(2)<LAMBDA(1)) then
          l = LAMBDA(1)
        end if
        !n = N(1):(N(2)-N(1))/m:N(2); if isempty(n), n = N(1); end
        allocate(nn(int(m)+1))
        do i = 1, int(m + 1)
          nn(i) = N(1)+(i-1)*(N(2)-N(1))/m
        end do
        if(N(2)<N(1)) then
          nn = N(1)
        end if
        pk = [-1000.0D0,-20.0D0,-6.0D0,0.0D0]
        dd = [20.0D0,2.0D0,0.1D0,0.0D0]
        !%hspan = [-1000 20 -20 0 ; -20 2 -5 0; -5 0.1 0]; [nr,nc]=size(hspan);
        !h = [];
        do i = 1, size(pk) - 1
            if (i == 1) then
                nk = 0;
            else
                nk = hspan(i-1, 4)
            end if
            !L = pk(i)+dd(i): dd(i): pk(i+1);
            do j = 1, int((pk(i+1)-(pk(i)+dd(i)))/dd(i)+1)
              L(j) = pk(i)+dd(i)+(j-1)*dd(i)
            end do
            hspan(i,:)=[pk(i),dd(i),pk(i+1),nk+size(L)]
            allocate(temp(size(h)))
            temp = h
            allocate(h(size(temp)+size(L)))
            h(1:size(temp)) = temp(:)
            h(size(temp)+1:size(h)) = L(:)
        end do
        !%h = [-1000:20:-20,-10:0.05:-0.1];
        !K = zeros(length(n),length(l),length(h));
        allocate(KK(size(n),size(l),size(h)))
        KK = 0.0D0
        !for i=1:length(n), for j=1:length(l), [h,THE,S,K(i,j,:)]=vanGenuchten(h,0.4,0.2,1,l(j),n(i),1,0); end; end
        do j = 1, size(l)
          do i = 1, size(n)
            !call vanGenuchten(h,0.4,0.2,1,l(j),n(i),1,0,h,THE,S,KK(i,j,:))
          end do
        end do
        !TK = zeros(size(K));
        allocate(TK(size(KK,1),size(KK,2),size(KK,3)))
        TK = 0.0D0
        !for i=1:length(n), for j=1:length(l), for k=2:length(h), TK(i,j,k)=trapz(h(1:k),K(i,j,1:k)); end; end; end
        do k = 2, size(h)
          do j = 1, size(l)
            do i = 1, size(n)
              ! call trapz
            end do
          end do
        end do
        !M = zeros(length(n),length(l));
        allocate(MM(size(n),size(l)))
        MM = 0.0D0
        !for i=1:length(n), for j=1:length(l), if any(K(i,j,2:end)-K(i,j,1:end-1)<0), M(i,j)=1; end; end; end
        do j = 1, size(l)
          do i = 1, size(n)
            if(any((KK(i,j,2:size(KK,3))-KK(i,j,1:(size(KK,3)-1)))<0)) then
              MM(i,j) = 1
            end if
          end do
        end do
        !%mesh(n,l,TK(:,:,30).*(1-M))
        GA%MM = MM
        GA%SW = TK
        GA%n =  nn
        GA%l = l
        GA%dn = (N(2)-N(1))/m
        GA%dl = (LAMBDA(2)-LAMBDA(1))/m
        GA%hspan = hspan
        GA%nspan = int(size(pk))-1
        GA%m = int([size(n),size(l),size(h)])
        GA%h =  h
        !GA.WFP = zeros(size(VDZ.NZC));
        if(.not. associated(GA%WFP)) allocate(GA%WFP(size(VDZ%NZC,1),size(VDZ%NZC,2)))
        GA%WFP = 0.0D0
        !GA.FM  = zeros(size(VDZ.h));
        if(.not. associated(GA%FM)) allocate(GA%FM(size(VDZ%h,1),size(VDZ%h,2),size(VDZ%h,3)))
        GA%FM = 0.0D0
        !GA.KBar = zeros(size(VDZ.h));
        if(.not. associated(GA%KBar)) allocate(GA%KBar(size(VDZ%h,1),size(VDZ%h,2),size(VDZ%h,3)))
        GA%KBar = 0.0D0
        !GA.Code  = zeros(size(VDZ.NZC)); % switch code
        if(.not. associated(GA%Code)) allocate(GA%Code(size(VDZ%NZC,1),size(VDZ%NZC,2)))
        GA%Code = 0.0D0
        !GA.THEW = zeros(size(VDZ.NZC)); % THE in the dry portion of the IWFP cell saved
        if(.not. associated(GA%THEW)) allocate(GA%THEW(size(VDZ%NZC,1),size(VDZ%NZC,2)))
        GA%THEW = 0.0D0
        GA%GCode = [1.0D0, 0.3D0, 0.12D0, 0.9D0]      !+zeros(1); % [GACode, TimesKthreshold, THE_threshold, Eff]
        !varargout{1} = GA;
      elseif (opt == 2) then        ! % pre-compute Green and Ampt coefficients
        !% FM, KBar
        !% FM = sum(DZ_i*(THES-THER)_i);
        !% Kbar = sum(DZ)/sum(DZ_i/KS_i)
        !VDZ = varargin{1};
        !VDZ = varargin
        allocate(DZ(size(VDZ%DZ,1),size(VDZ%DZ,2),size(VDZ%DZ,3)))
        DZ = VDZ%DZ
        !FM = zeros(size(VDZ.K));
        FM = 0.0D0
        !KBar = zeros(size(VDZ.K));
        KBar = 0.0D0
        !sumDZ = zeros(size(VDZ.K));
        allocate(sumDZ(size(VDZ%K,1),size(VDZ%K,2),size(VDZ%K,3)))
        sumDZ = 0.0D0
        ![nr,nc,nz]=size(VDZ.h);
        nr = size(VDZ%h,1)
        nc = size(VDZ%h,2)
        nz = size(VDZ%h,3)
        KBar(:,:,1) = 1.0D0     ! % doesn't matter
        sumDZ(:,:,1) = 0.0D0
        FM(:,:,1) = 0.0D0
        allocate(KS(size(VDZ%KS,1),size(VDZ%KS,2),size(VDZ%KS,3)))
        KS = VDZ%KS * VDZ%GA%GCode(4)
        do k = 2, nz
            do j = 1, nc
                do i = 1, nr
                    !if i==8 && j==38
                    !    4;
                    !end
                    sumDZ(i,j,k) = sumDZ(i,j,k-1) + DZ(i,j,k)
                    KBar(i,j,k) = sumDZ(i,j,k)/(sumDZ(i,j,k-1)/KBar(i,j,k-1)+DZ(i,j,k)/KS(i,j,k))
                    FM(i,j,k) = FM(i,j,k-1) + DZ(i,j,k) * (VDZ%THES(i,j,k) - VDZ%THER(i,j,k))
                end do
            end do
        end do
      end if

      if(allocated(h)) deallocate(h)
      if(allocated(temp)) deallocate(temp)
      if(allocated(KK)) deallocate(KK)
      if(allocated(TK)) deallocate(TK)
      if(allocated(MM)) deallocate(MM)
      if(allocated(l)) deallocate(l)
      if(allocated(nn)) deallocate(nn)
      if(allocated(DZ)) deallocate(DZ)
      if(allocated(sumDZ)) deallocate(sumDZ)
      if(allocated(KS)) deallocate(KS)
      if(allocated(THE)) deallocate(THE)
      if(allocated(S)) deallocate(S)

      end subroutine GreenAmptMatrix

!====================================================================================
      function eomday(y, m)
      !  %EOMDAY End of month.
      !  %   D = EOMDAY(Y,M) returns the last day of the month for the given
      !  %   year, Y, and month, M.
      !  %   Algorithm:
      !  %      "Thirty days hath September, ..."
      !  %
      !  %   See also WEEKDAY, DATENUM, DATEVEC.

      !  %   Copyright 1984-2002 The MathWorks, Inc.
      !  %   $Revision: 1.12 $  $Date: 2002/04/15 03:20:40 $
      implicit none
      integer :: d, y, m, dpm(12), eomday

      !% Number of days in the month.
      dpm = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

      ! Make result the right size and orientation.
      d = y - m

      d = dpm(m)
      if((m == 2) .and. ((mod(y,4) == 0 .and. mod(y,100) /= 0) .or. mod(y,400) == 0)) d = 29

      eomday = d

      end function eomday

!====================================================================================
      subroutine timer(time)
      integer, dimension(8) :: values
      real*8 :: year, month, day, hour, min, sec, milsec, time
      
      values = 0
      call date_and_time(values = values)
      year = 0d0
      year = values(1)
      month = 0d0
      month = values(2)
      day = 0d0
      day = values(3)
      hour = 0d0
      hour = values(5)
      min = 0d0
      min = values(6)
      sec = 0d0
      sec = values(7)
      milsec = 0d0
      milsec = values(8)
      
      time = 0D0
      time = hour * 3600.0D0 + min * 60.0D0 + sec + milsec / 1000.0D0
      
      
      end subroutine timer
!====================================================================================


      
      END MODULE VData
