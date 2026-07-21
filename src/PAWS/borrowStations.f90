
!

      subroutine borrowStations(Stations,nStations,i,field,startdate,enddate,dat,nk,bk)
      !% (E) borrow valid data records from neighboring stations for missing data.
      !% bk--# of bad data points in the end
      !% nk--# of borrowed data points in the end
      use Vdata
      implicit none
      integer :: nStations
      type(Stations_type) Stations(nStations)
      integer :: i, startdate, enddate, nk, bk, sp, ep, kkk, nr, nc, j, length, nbb, nkk
      character(len=*) field
      real*8 :: dat

      real*8, dimension(size(Stations)) :: shift

      real*8, dimension(size(Stations(i)%neighbor,1),size(Stations(i)%neighbor,2)) :: nb
      
      logical :: bor

      real*8 :: ub, dd, val

      ub = 200.0D0
      if (field=='rad' .or. field=='rrad') then
        ub = 20000.0D0
      endif
      if (startdate<19000000 .or. startdate>21000000 .or. startdate == enddate) then
        !% index themselves
        sp = startdate
        ep = enddate
        dat = 0.0D0
        bor = .true.
      else
        !sp = find(Stations(i).dates == startdate,1);
      !  do j = 1, size(Stations(i)%dates)
      !    if(Stations(i)%dates(i)==startdate) then
      !      sp = i
      !      exit
      !    endif
      !  enddo
        !ep = find(Stations(i).dates == enddate,1);
      !  do j = 1, size(Stations(i)%dates)
      !    if(Stations(i)%dates(j)==enddate) then
      !      ep = j
      !      exit
      !    endif
      !  enddo
        !dat = Stations(i).(field);
      !  select case (field)
      !    case ('tmin')
      !      dat = Stations(i)%tmin
      !    case ('tmax')
      !      dat = Stations(i)%tmax
      !  end select
        !bor = (dat<-90 | dat>ub | isnan(dat));
      !  where(dat<-90 .or. dat>ub .or. isnan(dat))
      !    bor = .true.
      !  endwhere
      endif
      !nb = Stations(i).neighbor; kkk = 1; [nr,nc]=size(nb);
      nb = Stations(i)%neighbor
      kkk = 1
      nr = size(nb,1)
      nc = size(nb,2)

      nk = 0
      bk = 0
      do j = 1, size(Stations)
        shift(j)=Stations(i)%datenums(1)-Stations(j)%datenums(1)
      enddo

      do j = sp, ep
        !dd = Stations(i).(field)(j)
        select case (field)
          case ('tmax')
            dd = Stations(i)%tmax(j)
          case ('tmin')
            dd = Stations(i)%tmin(j)
        end select
        kkk = 1
        !if (size(bor)==1 .or. bor(j)) then
        if(bor) then
            do while (kkk<=nr)
                nbb = nb(kkk,1)
                nkk = j + shift(nbb)
                !if nkk<1 || nkk>length(Stations(nbb).(field))
                select case (field)
                  case ('tmax')
                    length = size(Stations(nbb)%tmax)
                  case ('tmin')
                    length = size(Stations(nbb)%tmin)
                end select
                if (nkk < 1 .or. nkk > length) then
                    val = -99.0D0
                else
                  select case (field)
                    case ('tmax')
                      val = Stations(nbb)%tmax(nkk)
                    case ('tmin')
                      val = Stations(nbb)%tmin(nkk)
                  end select
                endif
                if (val<-90 .or. val>ub .or. isnan(val)) then
                  kkk=kkk+1;
                else
                  exit
                endif
            enddo
            !if (size(dat)>1) then
            !    if (kkk>nr) then
            !        dat(j) = -99.0D0
            !        bk = bk + 1
            !    else
            !        dat(j) = val
            !        nk = nk + 1
            !    endif
            !else
                if (kkk>nr) then
                    dat = -99.0D0
                    bk = bk + 1
                else
                    dat = val
                    nk = nk + 1
                endif
            !endif
        endif
      enddo
      !4;

      end subroutine borrowStations
