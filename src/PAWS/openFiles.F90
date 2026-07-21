!

!===================================================================================      
      
      subroutine printheader(vars,n,skip,comm,line)
      !! Junk characters may be generated when character function length is not the same as declared in caller
      !% print a header that aligns well with data
      implicit none
      integer n, skip, comm
      character*32, dimension(n) :: vars
      character*32 :: var, item
      character*(*) :: line
      integer :: i, pad
       
      line = ''
      do i = skip+1, size(vars)
          var = vars(i)
          pad = max(comm-len_trim(var)-1, 0)
          item = repeat(' ',pad)//trim(var)
          line = trim(line)//item
      enddo
      
      end subroutine printheader

      subroutine openFiles()
      !% This function opens output files for writing
      !% Files to be opened: 
      !% Env.recfile (1)Basin average budget record file (printed each time step)
      !%             (2)point recording files (use Toolsgui to load)
      !% This function also writes headers/basic formatting guide for these files

      !% Other forms of output not included in here:
      !% movOutput: time series of fluxes (will be done in a better way)
      !% mat saved. This saves the final states of variables which has the
      !% mapOutput: Temporal average of spatial maps
      use Vdata
      use budgeter_mod, Only: nitemw
      implicit none
      integer,parameter:: len = 25
      integer,parameter:: nitem = 54
      integer :: i, fid(2)
      character*32 :: varsm(100)
      character*(len*nitem) :: line
      character*(20) fmtstr
      logical of
      
      fid(1) = 30
      fid(2) = 40
      nitemw = nitem+1
      
      fmtstr = '(A9,A9,A6,A'//TRIM(num2str(len*nitem))//')'
      
      if (abs(w%t - w%ModelStartTime) < 1D-4) then
          !fid=fopen(Env.recfile{1},'wt');
          INQUIRE (FILE = trim(Env%recfile(1)), OPENED = of)
          if ( .NOT. Of) then
          open(UNIT=fid(1),FILE=trim(Env%recfile(1)),STATUS='REPLACE',ACTION='WRITE')
          !output='[id,w.dt,t,TA,TAVoTA,hf,hc,h,test2,rh,test3,THE,test4,GW,test5,CS,SWE,Mass,test1,0,Qout,Dperc2,prcp,Temp,SatE,InfE,botT,Rad,Dperc,NewS,SMELT,SUBM,Gprcp,ET,Inf,Qoc,Qgc,0,PET,Tp,EvapCS,EvapG]';
          !vars=cmd2var(output);s
          varsm(1:nitem) =(/'w.dt','t','TA','TAVoTA','hf','hc','h','test2','rh','test3',     &
         & 'THE','test4','GW','GWfluxOut','test5','CS','SWE','Mass','test1','Qout',         &
         & 'Dperc2','prcp','Temp','SatE','InfE','Ex','ExC','botT','Rad','Dperc','NewS','SMELT', &
         & 'SUBM','Gprcp','ET','Inf','Qoc','Qgc','PET','Tp','EvapCS',           &
         & 'EvapG','meanTHE','var2d','npp','gpp','ar','hr','SatFrac','DepH','Irrigationm','FrozenInfm','Soil2AirDirectm','Mass'/)
          call printheader(varsm,nitem,0,len,line)
          !fprintf(fid,['day,     time, id,',line,'\n']);
          write(UNIT=fid(1),FMT=fmtstr) 'day','time','id',line
          !fclose(fid);
          close(UNIT=fid(1))
          !% point recording file will just be opened. No header will be written
          open(UNIT=fid(2),FILE=trim(Env%recfile(2)),STATUS='REPLACE',ACTION='WRITE')
          close(UNIT=fid(2))
          endif
      endif
      
      Env%recfid(1) = fid(1)
      Env%recfid(2) = fid(2)
      !do i = 1, size(Env%recfile)
      INQUIRE (FILE = trim(Env%recfile(1)), OPENED = of)
      if ( .NOT. Of) then
        !Env.recfid(i) = fopen(Env.recfile{i},'at');
        open(UNIT=Env%recfid(1),FILE=trim(Env%recfile(1)),ACCESS='APPEND',ACTION='WRITE')
		open(UNIT=Env%recfid(2),FILE=trim(Env%recfile(2)),ACCESS='APPEND',ACTION='WRITE')
      !enddo
      endif
      end subroutine openFiles
