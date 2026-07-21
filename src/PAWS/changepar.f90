!

        subroutine changepar(opt, par, npar)
        !% this function is set to modify parameter values
        !% without an input, this will read 'pfile.dat'
        !% in the file pfile.dat:
        !% first line is reserved for comments,not used
        !% column 1 is the modification field to be modified. This is just the call
        !% name of the variable field, e.g. g.VDZ.th.SParm
        !% column 2 is the modification method, imet:
        !% 1:replace value by v; 2:add v to current; 3:multiply by 1+v
        !% column 3 is v
        !% column 4 is an upper limit; column 5 is an lower limit. Adjustment
        !% will be confined within these ranges. If no limit, put NaN. These are not
        !% the limits of parameter estimation, but to prevent against changing
        !% parameters to unreasonble values;
        !% column 6 is a reserved column currently not used, put 0
        !% rest of the columns, if not empty, then:
        !%    if it numeric, it is the index to be modifed--limit to one integer
        !%    if it is char, it is a file that contains the logical map to be
        !% modified (same format as w.DM.mask)
        !%           if empty, it will be applied to be entire field

        !% If this function has input
        !% it should be a structure array that contains field:
        !% field, imet, v, ul, ll, mask
        !% same as above
        !global g r w
        use Vdata
        implicit none
        integer opt, npar
        type(par_type), dimension(npar) :: par
        integer i, nr, nc, nz, k, n
        real*8, dimension(:,:), allocatable :: dat
        real*8, dimension(:), allocatable :: dat0
        real*8, dimension(:), pointer :: dat1D=>NULL()
        
        !if length(varargin)==0
        if(opt>0) then
            !par = readpar('changepar.dat');
            call readpar('changepar.dat', par, npar)
        elseif(opt==0) then
                return
            !par = varargin{1};
            !if isnumeric(par)
            !    return % no changepar, just a fake call
            !end
        endif

        do i = 1, size(par)
            if (i==1) then
                !disp(['Applying changepar to ',num2str(length(par)),' par set']);
                call display( 'Applying changepar to', size(par), ' par set')
            endif
            !str = ['[nr,nc,nz]=size(',par(i).field,');'];
            !eval(str)
            ! Since Fortran can not do the same thing as above
            ! have to do hard coding here
            select case (i)
              case(1)
                nr = 1
                nc = 1
                nz = 1
              case(2)
                nr = 1
                nc = 1
                nz = 1
              case(3)
                nr = size(g%VDZ%N,1)
                nc = size(g%VDZ%N,2)
                nz = size(g%VDZ%N,3)
              case(4)
                nr = size(g%VDZ%KS,1)
                nc = size(g%VDZ%KS,2)
                nz = size(g%VDZ%KS,3)
              case(5)
                nr = size(g%VDZ%ALPHA,1)
                nc = size(g%VDZ%ALPHA,2)
                nz = size(g%VDZ%ALPHA,3)              
              case(6)
                nr = size(g%GW(1)%K,1)
                nc = size(g%GW(1)%K,2)
                nz = 1              
              case(7)
                nr = size(g%OVN%dist,1)
                nc = size(g%OVN%dist,2)
                nz = 1              
              case(8)
                nr = 1
                nc = 1
                nz = 1
              case(9)
                nr = 1
                nc = 1
                nz = 1
              case(10)            
                nr = 1
                nc = 1
                nz = 1
              case(11)            
                nr = 1
                nc = 1
                nz = 1
              case(12)            
                nr = 1
                nc = 1
                nz = 1
            end select
            if(allocated(dat)) deallocate(dat)
            allocate(dat(nr,nc))
            do k = 1, nz
                !if (nz>1) then
                !    str = ['dat=',par(i).field,'(:,:,k);'];
                !else
                !    str = ['dat=',par(i).field,';'];
                !endif
                !eval(str)
                
                select case (i)
                  case(1)
                    dat = g%VDZ%VParm(2)
                  case(2)
                    dat = g%VDZ%VParm(1)
                  case(3)
                    dat = g%VDZ%N(:,:,k)
                  case(4)
                    dat = g%VDZ%KS(:,:,k)
                  case(5)
                    dat = g%VDZ%ALPHA(:,:,k)
                  case(6)
                    dat = g%GW(1)%K
                  case(7)
                    dat = g%OVN%dist
                  case(8)
                    dat = g%VDZ%VParm(5)
                  case(9)
                    dat = g%VDZ%th%SParm(15)
                  case(10)
                    dat = g%VDZ%VParm(9)
                  case(11)
                    dat = g%OVN%FParm(1)
                  case(12)
                    dat = g%OVN%FParm(3)
                end select
                
                !if ~isempty(par(i).id) && ~isnan(par(i).id) && par(i).id~=0
                dat1D => oneD_Ptr(size(dat),dat)
                if(allocated(dat0)) deallocate(dat0)
                if (associated(par(i)%id)) then
                  if((.not. isnan(dble(par(i)%id))) .and. (par(i)%id/=0)) then
                    allocate(dat0(1))
                    dat0 = dat1D(par(i)%id)
                  endif
                !elseif ~isempty(par(i).mapp)
                elseif (associated(par(i)%mapp)) then
                    allocate(dat0(size(par(i)%mapp)))
                    dat0 = dat1D(par(i)%mapp)
                else
                    allocate(dat0(size(dat1D)))
                    dat0 = dat1D
                endif
                
                !switch par(i).imet
                select case (par(i)%imet)
                    case (1)
                        dat0 = par(i)%v
                    case (2)
                        dat0 = par(i)%v + dat0
                    case (3)
                        dat0 = par(i)%v * dat0
                end select
                if (.not. isnan(par(i)%ll)) then
                    !dat0(dat0<par(i).ll)=par(i).ll;
                  where (dat0 < par(i)%ll) 
                    dat0 = par(i)%ll
                  endwhere
                endif
                if (.not. isnan(par(i)%ul)) then
                    !dat0(dat0>par(i).ul)=par(i).ul;
                  where (dat0 > par(i)%ul)
                    dat0 = par(i)%ul
                  endwhere
                endif
                !if ~isempty(par(i).id) && ~isnan(par(i).id) && par(i).id~=0
                if (associated(par(i)%id)) then
                  if((.not. isnan(dble(par(i)%id))) .and. (par(i)%id/=0)) then
                    dat1D([par(i)%id]) = oneD_ptr(size(dat0),dat0)
                  endif
                elseif (associated(par(i)%mapp)) then
                    dat1D(par(i)%mapp) = oneD_ptr(size(dat0),dat0)
                else
                    dat1D = dat0
                endif
                !if (nz>1) then
                !    str = [par(i).field,'(:,:,k)=dat;'];
                !else
                !    str = [par(i).field,'=dat;'];
                !endif
                !eval(str)
                select case (i)
                  case(1)
                    g%VDZ%VParm(2) = dat(1,1)
                  case(2)
                    g%VDZ%VParm(1) = dat(1,1)
                  case(3)
                    g%VDZ%N(:,:,k) = dat
                  case(4)
                    g%VDZ%KS(:,:,k) = dat
                  case(5)
                    g%VDZ%ALPHA(:,:,k) = dat
                  case(6)
                    g%GW(1)%K = dat
                  case(7)
                    g%OVN%dist = dat
                  case(8)
                    g%VDZ%VParm(5) = dat(1,1)
                  case(9)
                    g%VDZ%th%SParm(15) = dat(1,1)
                  case(10)
                    g%VDZ%VParm(9) = dat(1,1)
                  case(11)
                    g%OVN%FParm(1) = dat(1,1)
                  case(12)
                    g%OVN%FParm(3) = dat(1,1)
                end select
            enddo
        enddo
        
        if(allocated(dat)) deallocate(dat)
        if(allocated(dat0)) deallocate(dat0)
        end subroutine changepar


        function strmatch(nd,str_array,str)
        implicit none
        integer ns, nd, strmatch, i
        character(len=*),dimension(nd):: str_array
        character(len=*)::str
        strmatch = 0
        do i = 1,nd
          IF (TRIM(adjustl(str)) == TRIM( adjustl(str_array(i)) )) THEN
             strmatch = i
            RETURN
          ENDIF
        enddo
        end function strmatch

!

        subroutine changepar2(opt, par, npar)
        !cp: changepar2 is more flexible than changepar
        !% this function is set to modify parameter values
        !% without an input, this will read 'pfile.dat'
        !% in the file pfile.dat:
        !% first line is reserved for comments,not used
        !% column 1 is the modification field to be modified. This is just the call
        !% name of the variable field, e.g. g.VDZ.th.SParm
        !% column 2 is the modification method, imet:
        !% 1:replace value by v; 2:add v to current; 3:multiply by 1+v
        !% column 3 is v
        !% column 4 is an upper limit; column 5 is an lower limit. Adjustment
        !% will be confined within these ranges. If no limit, put NaN. These are not
        !% the limits of parameter estimation, but to prevent against changing
        !% parameters to unreasonble values;
        !% column 6 is a reserved column currently not used, put 0
        !% rest of the columns, if not empty, then:
        !%    if it numeric, it is the index to be modifed--limit to one integer
        !%    if it is char, it is a file that contains the logical map to be
        !% modified (same format as w.DM.mask)
        !%           if empty, it will be applied to be entire field

        !% If this function has input
        !% it should be a structure array that contains field:
        !% field, imet, v, ul, ll, mask
        !% same as above
        !global g r w
        use Vdata
        implicit none
        integer opt, npar
        type(par_type), dimension(npar) :: par
        integer i, nr, nc, nz, k, n, kk
        real*8, dimension(:), pointer :: dat=>NULL()
        real*8, dimension(:), allocatable :: dat0
        real*8, dimension(:), pointer :: dat1D=>NULL()
        real*8, dimension(:), pointer :: dat_temp=>NULL()
        integer,parameter:: npar_base = 9
        character(len=100),dimension(npar_base):: par_base
        integer,pointer,dimension(:) :: mask1D=>NULL()
        integer :: strmatch
        
        !if length(varargin)==0
        if(opt>0) then
            !par = readpar('changepar.dat');
            call readpar('changepar.dat', par, npar)
        elseif(opt==0) then
                return
            !par = varargin{1};
            !if isnumeric(par)
            !    return % no changepar, just a fake call
            !end
        endif
        
        par_base(1:npar_base) = (/                     &
                      &         'g.VDZ.VParm   ',         &
                      &         'g.VDZ.N       ',             &
                      &         'g.VDZ.KS      ',            &
                      &         'g.VDZ.ALPHA   ',         &
                      &         'g.GW(1).K     ',           &
                      &         'g.GW(2).K     ',           &
                      &         'g.OVN.dist    ',          &
                      &         'g.VDZ.th.SParm',      &
                      &         'g.OVN.FParm   '          &
                      &         /)
        mask1D  => oneD_ptrI(size(w%DM%MapL),w%DM%MapL)
        

        do i = 1, size(par)
            !str = ['[nr,nc,nz]=size(',par(i).field,');'];
            !eval(str)
            ! Since Fortran can not do the same thing as above
            ! have to do hard coding here
            
            kk = strmatch(npar_base,par_base,par(i)%field)
            IF (kk .ne. 0) THEN
            call display( 'Applying changepar to par',i,trim(par(i)%field)//' id: '//num2str(par(i)%id) )
            
                !if (nz>1) then
                !    str = ['dat=',par(i).field,'(:,:,k);'];
                !else
                !    str = ['dat=',par(i).field,';'];
                !endif
                !eval(str)
                
                select case (kk)
                  case(1)
                    dat => oneD_ptr(size( g%VDZ%VParm),g%VDZ%VParm)
                  case(2)
                    dat => oneD_ptr(size( g%VDZ%N)    ,g%VDZ%N)
                  case(3)
                    dat => oneD_ptr(size( g%VDZ%KS)   ,g%VDZ%KS)
                  case(4)
                    dat => oneD_ptr(size( g%VDZ%ALPHA),g%VDZ%ALPHA)
                  case(5)
                    dat => oneD_ptr(size( g%GW(1)%K)  ,g%GW(1)%K)
                  case(6)
                    dat => oneD_ptr(size( g%GW(2)%K)  ,g%GW(2)%K)
                  case(7)
                    dat => oneD_ptr(size( g%OVN%dist) ,g%OVN%dist)
                  case(8)
                    dat => oneD_ptr(size( g%VDZ%th%SParm),g%VDZ%th%SParm)
                  case(9)
                    dat => oneD_ptr(size( g%OVN%FParm),g%OVN%FParm)
                  case default
                    call display( 'invalid parameter field!')
                    stop
                end select
                
                !if ~isempty(par(i).id) && ~isnan(par(i).id) && par(i).id~=0
                dat1D => dat
                
                if(allocated(dat0)) deallocate(dat0)
                if (associated(par(i)%id)) then
                  if((.not. isnan(dble(par(i)%id))) .and. (par(i)%id/=0)) then
                    allocate(dat0(1))
                    dat0 = dat1D(par(i)%id)
                  elseif (par(i)%id .eq. 0) then
                    allocate(dat0( size(dat1D) ))
                    dat0 = dat1D
                  else
                    call display( 'invalid par!!!!')
                    stop
                  endif
                !elseif ~isempty(par(i).mapp)
                elseif (associated(par(i)%mapp)) then
                    allocate(dat0(size(par(i)%mapp)))
                    dat0 = dat1D(par(i)%mapp)
                else
                    allocate(dat0(size(dat1D)))
                    dat0 = dat1D
                endif
                
                !switch par(i).imet
                select case (par(i)%imet)
                    case (1)
                        dat0 = par(i)%v
                    case (2)
                        dat0 = par(i)%v + dat0
                    case (3)
                        dat0 = par(i)%v * dat0
                end select
                
                if (.not. isnan(par(i)%ll)) then
                    !dat0(dat0<par(i).ll)=par(i).ll;
                  where (dat0 < par(i)%ll) 
                    dat0 = par(i)%ll
                  endwhere
                endif
                if (.not. isnan(par(i)%ul)) then
                    !dat0(dat0>par(i).ul)=par(i).ul;
                  where (dat0 > par(i)%ul)
                    dat0 = par(i)%ul
                  endwhere
                endif
                !if ~isempty(par(i).id) && ~isnan(par(i).id) && par(i).id~=0
                if (associated(par(i)%id)) then
                  if((.not. isnan(dble(par(i)%id))) .and. (par(i)%id/=0)) then
                    dat1D([par(i)%id]) = oneD_ptr(size(dat0),dat0)
                  elseif (par(i)%id .eq. 0) then
                    dat1D = dat0
                  else
                    call display( 'invalid par!!!!')
                    stop
                  endif
                elseif (associated(par(i)%mapp)) then
                    dat1D(par(i)%mapp) = oneD_ptr(size(dat0),dat0)
                else
                    dat1D = dat0
                endif
            ELSE
                call display( 'invalid parameter field!')
                stop
            ENDIF
        enddo
        
        end subroutine changepar2
