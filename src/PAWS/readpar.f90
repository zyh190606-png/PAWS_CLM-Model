
      subroutine readpar(file, par, npar)
      use Vdata
      implicit none
      character(len=*) file
      integer npar
      type(par_type), dimension(100) :: par
      integer status, k
      character*200 line
      !% see changepar for the format of the file
      !if ischar(file)
      !  fid = fopen(file,'rt');
      !elseif isnumeric(file) && file>2
      !fid = file;
      !end
      !if fid <=2
      !  disp('no par to change');
      !  return
      !end
      open(unit=10,file=trim(file),status='old',action='read',iostat=status)
      if(status/=0) then
        call display( 'Error open file : '// file)
        npar = 0;
        return
      end if
      
      !line = fgetl(fid);
      call display( 'Reading parameter change file : '// file)
      read(10,'(A200)') line
      k = 0;
      do
        !line = fgetl(fid);
        read(10,'(A200)',iostat=status) line
        !if ~ischar(line) && line==-1, break, exit, end
        if(status/=0 .or. trim(line) == trim('   ')) exit
        k = k+1
        !A=textscan(line,'%s %f %f %f %f %f %s');
        !if isempty(A{end})
        !    id = 0; mapfile = '';
        !else
        !    try
        !        id = str2num(A{end}{1}); mapfile='';
        !    catch
        !        id = 0; mapfile = A{end};
        !    end
        !end
        !par(k).field = A{1}{1};
        !par(k).imet = A{2};
        !par(k).v = A{3};
        !par(k).ul = A{4};
        !par(k).ll = A{5};
        !par(k).code = A{6};
        allocate(par(k)%field, par(k)%imet, par(k)%v, par(k)%ul, par(k)%ll, par(k)%code, par(k)%id)
        read(line,'(A32,I6,e16.6,f9.2,e10.2,I7, I10)') par(k)%field, par(k)%imet, par(k)%v, par(k)%ul, par(k)%ll, par(k)%code, par(k)%id

        !par(k).id = id; par(k).mapfile = mapfile; par(k).mapp = [];
        !if ~isempty(mapfile)
        !    try
        !        par(k).mapp = logical(load(mapfile));
        !    catch
        !        error(['error while reading mapfile: ',mapfile,' of column ',num2str(k)]);
        !    end
        !end
      end do
      npar = k
      !if ischar(file)
      !  fclose(fid);
      !end      

      close(unit=10)
      
      end subroutine readpar
