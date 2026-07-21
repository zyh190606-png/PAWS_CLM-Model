      subroutine recAcc(obj,varargin)
      !subroutine recAcc(obj,varargin)
      ! (E) record accumulative fluxes before they are cleared for mass budgeting
      ! This is generally called for river flux accumulation
      ! Input/Output self-explanatory in the codes
      use Vdata
      use memoryChunk, Only: tempStorage
      implicit none
      character*32 obj
      integer:: varargin
      integer id, i, j
      type(robj_type) s
      integer rid, objnum, is
      integer ub, lb, nlayer
      real*8 :: temp, z,zz
      character*32 :: F1, F2
      real*8, dimension(:), pointer:: QinRiv=>null()
      logical::  recV= .true.
      
      call tempStorage('QinRiv',QinRiv, size(w%CR)) 
      select case (trim(obj))
        case ('rMB')
            ! This is for mass balance
            !id = varargin{1};
            id = varargin
            !j = find(w.CR==id,1);
            do i = 1, size(w%CR)
              if(w%CR(i)==id) then
                j = i
                exit
              end if
            end do
            ! all units are mass (m3). Qoc is accumulated inside F_oc_dw
            ! w.Rec.acc.Qoc(j) = w.Rec.acc.Qoc(j) + sum(r(id).Riv.Qoc)*(r(id).Riv.dt*86400); 
            ! cell1 = round((79976-8615)/r(25).DM.d);
            !if (id==w%CR(size(w%CR))) then
              w%Rec%acc%Qout(j) = w%Rec%acc%Qout(j) + r(id)%Riv%Qx(size(r(id)%Riv%Qx)-1)*(r(id)%Riv%dt*86400.0D0)
              QinRiv(j) = QinRiv(j) + r(id)%Riv%Qx(1)*(r(id)%Riv%dt*86400.0D0)
            !end if

            w%Rec%acc%trib(j) = w%Rec%acc%trib(j) + sum(r(id)%Riv%trib)*(r(id)%Riv%dt*86400.0D0)
            
            if (recV) then
            ! 10 is Qoc m3
            is = 16;
            r(id)%Riv%vrv(is,1) = r(id)%Riv%vrv(is,1) + r(id)%Riv%dt ! accumulated time -- day
            r(id)%Riv%vrv(is,2) = r(id)%Riv%vrv(is,2) + 1           ! # of accumulations -- day
            r(id)%Riv%vrv(is,3) = w%t                               ! time stamp
            is = is + 1 !17
            r(id)%Riv%vrv(is,:) = r(id)%Riv%vrv(is,:) + r(id)%Riv%Qx*(r(id)%Riv%dt*86400D0) ! discharge m3
            is = is + 1 !18
            r(id)%Riv%vrv(is,:) = r(id)%Riv%vrv(is,:) + r(id)%Riv%U ! velocity m/s
            is = is + 1 !19
            r(id)%Riv%vrv(is,:) = r(id)%Riv%vrv(is,:) + r(id)%Riv%h ! height m
            is = is + 1 !20
            r(id)%Riv%vrv(is,:) = r(id)%Riv%vrv(is,:) + r(id)%Riv%hg ! GW height m
            if (abs(w%t+w%dt - r(id)%Riv%t)<1D-8) then
               is = is + 1 !21
               r(id)%Riv%vrv(is,:) = r(id)%Riv%vrv(is,:) + r(id)%Riv%fG ! Qgc m3
            endif
            
            endif
    
        case ('r')
            do i = 1, size(w%Rec%robj)
                s = w%Rec%robj(i)
                if (associated(s%type)) then
                if (s%type == 'r') then
                    rid = s%id
                    w%Rec%robj(i)%acc = w%Rec%robj(i)%acc + r(rid)%Riv%Qx(s%ix(1))*(r(rid)%Riv%dt*86400.0D0)
                    w%Rec%robj(i)%accT = w%Rec%robj(i)%accT + r(rid)%Riv%dt
                end if
                end if
            end do
        case ('g')
            do i = 1, size(w%Rec%robj)
                s = w%Rec%robj(i)
                if (s%type == 'g') then
                    F1 = w%Rec%robj(i)%Field(1)
                    F2 = w%Rec%robj(i)%Field(2)
                    objnum = w%Rec%robj(i)%objnum
                    ! w%Rec%robj(i)%acc = g%(F1)(objnum)%(F2)(s%ix(1),s%ix(2)) ???
                    if (trim(F1)=='GW') then
                        if (trim(F2)=='h') then
                            if(.not. associated(w%Rec%robj(i)%acc)) allocate(w%Rec%robj(i)%acc(1))
                            w%Rec%robj(i)%acc = g%GW(objnum)%h(s%ix(1),s%ix(2))
                        endif
                    endif
                end if
                ! by Kuai, write acc of soil moisture
                if (s%type == 's') then                    
                    ub=s%ix(3)
                    lb=s%ix(4)
                    temp=0
                    zz=0
                    do j=ub,lb
                        z=g%VDZ%DZ(s%ix(1),s%ix(2),j)    
                        temp=temp+g%VDZ%THE(s%ix(1),s%ix(2),j)*z
                        zz=zz+z
                    enddo
                    if(.not. associated(w%Rec%robj(i)%acc)) allocate(w%Rec%robj(i)%acc(1))
                    w%Rec%robj(i)%acc = temp / zz                    
                endif
            end do
      end select
      ! output(rid,dist,ddist,xy,ix,sidx) --- objnum, Field{1:2}
      
      end subroutine recAcc
!====================================================================================
