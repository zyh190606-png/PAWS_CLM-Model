
    !
    subroutine run(file_name)
    ! run --without variable, it reads 'prism.in' as command file
    ! run('filename.txt') uses 'filename.txt' as command file
    ! run('matfile.mat') loads 'matfile.mat' and run
    ! run(0) just run current memory values
    ! run(arg1,par) will first modify parameters before run, par can be a
    ! filename sent to changepar, or a parameter change structure
    ! run(arg1,par,outfile) will save some output to a file for quick
    ! comparison purposes.
    ! the command file has project name in the first line
    ! mat file name (which to load) in the second line
    ! starttime in the second line and endtime in the third line (yyyymmdd)
    !#function matlabrc
    ! if isdeployed
    ! #function matlabrc
    ! else
    ! gpath
    ! end
    ! if 0
    ! matlabrc % to include this function is to avoid some compiler issue
    ! end

    ! gOBJ
    ! format short g
    use Vdata
    use paws_varctl, Only : mode, nrepeat, nMatOutputFreq
    use displayMod,  Only : printCode, printID
    implicit none
    character*100 :: file='', mat='', outfile='', pfile='', rfile='',fileEx=''
    character(*)::file_name
    real*8 :: cp = 0.0D0
    integer :: npar,eof
    type(par_type),dimension(100)::par
    character*100 :: prj = 'wBudgeter'
    character*100 :: date, dates, time, times, dummy, printFileName
    integer :: flag=0, status, printCodeR = 0, modeBC = 0
    real*8 :: ts=0D0, te=0D0, dt=0D0
    real*8, dimension(72) :: ET_sim
    ! if length(varargin)>=1
    ! if isnumeric(varargin{1})
    ! else
    ! if ischar(varargin{1}) && strcmp(varargin{1}(end-2:end),'.in')
    ! file = varargin{1};
    ! elseif ischar(varargin{1})
    ! mat = varargin{1};
    ! end
    ! end
    ! if length(varargin)>=2
    ! par = varargin{2};
    ! end
    ! if length(varargin)>=3
    ! prj = varargin{3};
    ! end
    ! else
    ! file = 'prism.in';
    ! end
    flag = index(file_name, '.in')
    if(flag > 0 .and. flag <= len_trim(file_name)) then
        file = file_name
    else
        file = ''
        !file = 'prism.in'
    endif
    flag = index(file_name, '.mat')
    if(flag > 0 .and. flag <= len_trim(file_name)) then
        mat = file_name
    endif
    flag = index(file_name, '.txt') ! this mode is reserved for testing only
    if(flag > 0 .and. flag <= len_trim(file_name)) then
        fileEx = file_name
    else
        fileEx = ''
        !file = 'prism.in'
    endif

    call display( 'Driver control file: '// file)
    call display( 'mat initalization  : '// mat)
    !  if ~isempty(file)
    if(len_trim(file) > 0) then
        !      clear global Prob w g r q
        !      fid = fopen(file,'rt');
        open(unit=15,file=file,status='old',action='read',iostat=status)
        if(status/=0) then
            call display( 'Error open file : '// file)
            return
        end if
        !      prj = fgetl(fid);
        read(15,'(A100)') prj
        !      mat = fgetl(fid);
        read(15,'(A100)') mat
        !      load(mat)
        call bridgeGlobal2(mat)
        call bridgeMChunk(mat) ! link CLM data or not

        !      global w
        !      line = fgetl(fid);
        !      if length(line)>1
        !          disp(['1',line])
        !          y = num2str(line(1:4)); m = num2str(line(5:6)); d = num2str(line(7:8));
        !          st = datenum2([y m d]);
        !          wSyncTime(1,st);
        !      end
        !      line = fgetl(fid);
        !      if length(line)>1
        !          disp('2',line)
        !          y = num2str(line(1:4)); m = num2str(line(5:6)); d = num2str(line(7:8));
        !          et = datenum2([y m d]);
        !          w(1).ModelEndTime = et;
        !      else
        !      end
        !      pfile = fgetl(fid);
        read(15,*)
        read(15,*)
        read(15,'(A100)') pfile
        !      if length(pfile)>1 && ((~isstruct(par) && ~ischar(par) && par == 999) || length(varargin)<2)
        !          par = readpar(pfile);
        npar = 100
        call readpar(pfile, par, npar)
        !          % only use the pfile specified if it is an existing file on the
        !          % path. Passed in valid par will over-ride this file, unless it is 999, which is a signal
        !      end
        !      fclose(fid);

        ! new control option: nMatOutputFreq:
        !0. not period mat control (but will output final _m file); 1. monthly output (original)
        nMatOutputFreq = 1
        read(15,*,iostat=eof) mode,nrepeat,nMatOutputFreq

        if (eof .ne. 0) then
            mode = 'normal'
            nrepeat = 1
            nMatOutputFreq = 1
        endif

        if (nMatOutputFreq .eq. 0) then
            call display('No month.mat will be saved')
        endif

        read(15,*,iostat=eof) printCodeR, printFileName
        if (eof .eq. 0) then
            if (printCodeR .ne. 0) call display('Now output will be piped to: '// trim(printFileName))
            printCode = printCodeR
            open(unit=printID, file=trim(printFileName), action='write')
        endif

        close(unit=15)
        !  elseif ~isempty(mat)
        !      clear global Prob w g r
        !      load(mat)
        !  end
    elseif(len_trim(mat) > 0) then
        call bridgeGlobal2(mat)
        call bridgeMChunk(mat) ! link CLM data or not
        npar = 0
    elseif(len_trim(fileEx) > 0) then
        open(unit=15,file=fileEx,status='old',action='read',iostat=status)
        if(status/=0) then
            call display( 'Error open file : '// file)
            return
        end if
        !      prj = fgetl(fid);
        read(15,'(A100)') prj
        !      mat = fgetl(fid);
        read(15,'(A100)') mat
        call bridgeGlobal2(mat)
        read(15,*) ts, te, dt
        read(15,*) modeBC
        !      load(mat)
        if (dt .eq. 0) dt = w%dt
        call wSyncTime(ts); w%ModelEndTime = te
        call solversSetDT(dT)
        close(unit=15)
    endif
    ! global Env
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Env%recfile(1) = trim(prj)//'.txt'
    Env%recfile(2) = trim(prj)//'_Rec.txt'
    ! %%% this filename may be modified
    ! %%% related operations are:
    ! %%% removing previous sim files
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! % modify parameters
    ! %changepar(par); % if par is numeric, will not do anything.
    ! global w
    rfile = trim(prj)//'_result.txt'
    call display( 'Model Project Name : ', prj)
    call date_and_time(date, time)
    call display( 'Start Running Model at wallclock: '//date(1:4)//'/'//date(5:6)//'/'//date(7:8)// &
        &      '  '// time(1:2)//':'//time(3:4)//':'//time(5:6))
    dates = date; times = time
    call display( 'File Output to: '// trim(Env%recfile(1))// ' & '// trim(Env%recfile(2))//' & '// trim(rfile))
    call display( 'Model Start Time: '// num2str(datevec2(int(w%t))) )
    call display( 'Model End Time: '// num2str(datevec2(int(w%ModelEndTime))) )

    w%timezone = -5
    call display( 'Using U.S. Eastern Standard Time')

    if (w.t < w.ModelEndTime) then

        if (len_trim(fileEx) > 0) then
            ! testcases entrance
            call timeStepping(modeBC, te)
        else
            ! main model entrance
            call calendarMarch(par,npar,1,ET_sim);
        endif

        call date_and_time(date, time)
        call display( 'Model started at wallclock: ', dates(1:4)//'/'//dates(5:6)//'/'//dates(7:8), &
            &      '  ', times(1:2)//':'//times(3:4)//':'//times(5:6))
        call display( 'Simulation Successfully Completed at wallclock: '// date(1:4)//'/'//date(5:6)//'/'//date(7:8)// &
            &      '  '// time(1:2)//':'//time(3:4)//':'//time(5:6))
        call display( 'Model End Time: '// num2str(datevec2(int(w%t))) )
        !  if isfield(w,'tData') && isfield(w.tData,'usgs') && isfield(w.tData.usgs(1),'sd') && ~isempty(w.tData.usgs(1).sd)
        close(UNIT=Env%recfid(1))
        close(UNIT=Env%recfid(2))

        !call mapOutput(1, w%tData%maps) ! using mapOutput2 in wtrshd_day
        call saveMat2([trim(prj)//'_m.mat'])
        call display( 'saveMat2 completed saving mapOutput data to mat file')
    else
        ! otherwise it is a result collection run
    endif

    if (associated(w%tData%usgs) .and. associated(w%tData%usgs(1)%sd)) then
        !      % evaluate usgs comparison
        !      %w.tData.stats=simResults(Env.recfile,rfile,w.tData.ns);
        !      %w.tData.stats=simResults2(Env.recfile,rfile,w.tData.ns);
        !      [w.tData.objfun,w.tData.comp,w.tData.stats]=simResults3(Env.recfile{2},rfile);
        call simResults3(Env%recfile(2),rfile,w%tData%objfun,ET_sim)
        !      %%% recfile, rfile,          input file, pfile, matfile
        !      %%% [prj.txt,prj_result.txt, (input) jobname.in, changepar.dat, (input)
        !      % stats = [NASH,RMSE,NRNASH,R2,F,p,errVar];
        !  end
        call display( 'objective function calculated')
    endif
    !  fclose('all');
    !  if ~isempty(prj)
    !      %w.tData.S1 = S1;
    !      w.tData.maps = mapOutput(1);
    !      cmd_Save([prj,'_m']);
    !  end
    call saveMat2([trim(prj)//'_m.mat'])
    if (printCodeR .ne. 0) close(printID)
    !  % if ~isempty(prj)
    !  %     save([prj,'_m'],'S1')
    !  % end


    end subroutine run
