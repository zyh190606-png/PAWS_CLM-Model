!
        subroutine simResults3(outfile,rfile,objfun,ET_sim)
        !function [objfun,comp,stats] = simResults3(outfile,rfile)
        !% established 04/30/2010
        !% evaluate both streamflow and groundwater
        !% using the observation and record files
        !global w
        use Vdata
        use pawsControl, Only: CaliStartTime
        implicit none
		character(len=*) :: outfile, rfile
		real*8 :: objfun, weight, res(3),ET_obsbar,ET_NASH
		real*8, dimension(:,:), allocatable :: P
		character*300 :: fm, line
		integer :: status, fid, k, I, n, m, sd, ed, eid, ncell
		integer :: ii, jj, n1, n2, nn, ns,nk, msd
		type(robj_type), pointer :: rec=>NULL()
		real*8, pointer :: cc=>NULL()
		type(sim_type) :: sim, obs
		real*8, dimension(4) :: s
		real*8, dimension(:,:),pointer :: map1=>NULL(),map2=>NULL(),tempdata=>NULL()
		real*8, dimension(72) :: ET_sim,ET_obs
        !% w.rec.robj.comp -- integer is the usgs (or usgs_gw) number to compare to,
        !% the following digits stands for m, the method for calculating the
        !% objective function
        !% w is the weight, how much weight is assigned to this single objective
        !% function

        !P = load(outfile); objfun = 0;
        objfun = 0.0D0
		if(allocated(P)) deallocate(P)
		allocate(P(int(24*(w%ModelEndTime-w%ModelStartTime+1)), size(w%Rec%robj)+1))
		fid = Env%recfid(2)
		open(UNIT=fid,FILE=outfile,STATUS='OLD',ACTION='READ',IOSTAT=status)
	    if(status/=0) then
		  call display( 'Error open file : '// outfile)
		  return
	    end if
		!fm = repeat('1F13.6,',size(w%Rec%robj)-1)
        !fm = '(1F13.3,'//trim(fm)//'1F13.6)'
		k = 0
		do
		    read(fid,'(A300)',iostat=status) line
		    !if ~ischar(line) && line==-1, break, exit, end
		    if(status/=0 .or. k>=size(P,1)) exit
		    k = k+1
		    !read(line,fmt=trim(fm)) (P(k,I),I=1,size(w%Rec%robj)+1)
		    read(line,*) (P(k,I),I=1,size(w%Rec%robj)+1)
		enddo
		close(unit=fid)

        do i = 1, size(w%Rec%robj)
            rec => w%Rec%robj(i)
            cc => rec%cc
            if(associated(cc)) then
              if (cc==0) cc=>NULL()
            endif
            msd = int(CaliStartTime)
	    select case (rec%type)
                case ('r')
                    !sim.t = P(:,1); sim.v = P(:,i+1); sim = wDaily(sim);
                    if(allocated(sim%t)) deallocate(sim%t)
                    allocate(sim%t(k))
                    if(allocated(sim%v)) deallocate(sim%v)
                    allocate(sim%v(k))
                    sim%t = P(1:k,1); sim%v = P(1:k,i+1)
                    call wDaily(sim, 0)
                    !if ~isempty(cc) && cc~=0
                    if(associated(cc) .and. cc/=0) then
                        n= floor(cc); m = nint((cc - DBLE(n))*10)
                        !obs = w.tData.usgs(n); weight = rec.weight;
                        if(allocated(obs%t)) deallocate(obs%t)
                        allocate(obs%t(size(w%tData%usgs(n)%t)))
                        obs%t = w%tData%usgs(n)%t
                        if(allocated(obs%v)) deallocate(obs%v)
                        allocate(obs%v(size(w%tData%usgs(n)%v)))
                        obs%v = w%tData%usgs(n)%v
                        obs%sd => w%tData%usgs(n)%sd
                        weight = rec%weight
                    endif
                case ('g')
                    !sim.t = P(:,1); sim.v = P(:,i+1); sim = wDaily(sim);
                    if(allocated(sim%t)) deallocate(sim%t)
                    allocate(sim%t(k))
                    if(allocated(sim%v)) deallocate(sim%v)
                    allocate(sim%v(k))
                    sim%t = P(1:k,1); sim%v = P(1:k,i+1)
                    call wDaily(sim, 0)
                    cc => rec%cc
                    !if ~isempty(cc) && cc~=0
                    if(associated(cc) .and. cc/=0) then
                        n= floor(cc); m = nint((cc - DBLE(n))*10)
                        !obs = w.tData.usgs_gw(n); weight = rec.weight;
                        if(allocated(obs%t)) deallocate(obs%t)
                        allocate(obs%t(size(w%tData%usgs_gw(n)%t)))
                        obs%t = w%tData%usgs_gw(n)%t
                        if(allocated(obs%v)) deallocate(obs%v)
                        allocate(obs%v(size(w%tData%usgs_gw(n)%v)))
                        obs%v = w%tData%usgs_gw(n)%v
                        obs%sd => w%tData%usgs_gw(n)%sd
                        weight = rec%weight
                    endif
                case ('m')  !grace
                    if(allocated(sim%t)) deallocate(sim%t)
                    allocate(sim%t(k))
                    if(allocated(sim%v)) deallocate(sim%v)
                    allocate(sim%v(k))
                    sim%t = P(1:k,1); sim%v = P(1:k,i+1)
                    call wDaily(sim, 0)
                    if(associated(cc) .and. cc/=0) then
                        n= floor(cc); m = nint((cc - DBLE(n))*10)
                        !obs = w.tData.usgs_gw(n); weight = rec.weight;
                        if(allocated(obs%t)) deallocate(obs%t)
                        allocate(obs%t(size(w%tData%usgs(n)%t)))
                        obs%t = w%tData%usgs(n)%t
                        if(allocated(obs%v)) deallocate(obs%v)
                        allocate(obs%v(size(w%tData%usgs(n)%v)))
                        obs%v = w%tData%usgs(n)%v
                        obs%sd => w%tData%usgs(n)%sd
                        weight = rec%weight
                    endif 
                case ('s')  !soil moisture
                    if(allocated(sim%t)) deallocate(sim%t)
                    allocate(sim%t(k))
                    if(allocated(sim%v)) deallocate(sim%v)
                    allocate(sim%v(k))
                    sim%t = P(1:k,1); sim%v = P(1:k,i+1)
                    call wDaily(sim, 0)
                    if(associated(cc) .and. cc/=0) then
                        n= floor(cc); m = nint((cc - DBLE(n))*10)
                        !obs = w.tData.usgs_gw(n); weight = rec.weight;
                        if(allocated(obs%t)) deallocate(obs%t)
                        allocate(obs%t(size(w%tData%usgs(n)%t)))
                        obs%t = w%tData%usgs(n)%t
                        if(allocated(obs%v)) deallocate(obs%v)
                        allocate(obs%v(size(w%tData%usgs(n)%v)))
                        obs%v = w%tData%usgs(n)%v
                        obs%sd => w%tData%usgs(n)%sd
                        weight = rec%weight
                    endif 
                    
                case ('h')
                   !compare spatial map of groundwater head
                   nn = size(w%DM%mask)
                   map1 => w%tData%h1save
                   tempdata => w%tData%maps(1)%GW
                   allocate(map2(size(map1,1), size(map1,2)) )
                   map2 = tempdata/w%tData%maps(1)%n
                   
                   ncell = 0
                   do ii = 1,size(w%DM%mask,1)
                   do jj = 1,size(w%DM%mask,2)
                     if (w%DM%mask(ii,jj)) then
                        ncell = ncell + 1
                     endif
                   enddo
                   enddo
                   if(allocated(obs%t)) deallocate(obs%t)
                   allocate(obs%t(ncell))
                   if(allocated(obs%v)) deallocate(obs%v)
                   allocate(obs%v(ncell))
                   if(allocated(sim%t)) deallocate(sim%t)
                   allocate(sim%t(ncell))
                   if(allocated(sim%v)) deallocate(sim%v)
                   allocate(sim%v(ncell))
                   
                   nk = 0
                   do jj = 1,size(w%DM%mask,2)
                   do ii = 1,size(w%DM%mask,1)
                     if (w%DM%mask(ii,jj)) then
                        nk = nk + 1
                        obs%t(nk)=nk
                        sim%t(nk)=nk
                        obs%v(nk)=g%Topo%E(ii,jj)-map1(ii,jj)
                        sim%v(nk)=g%Topo%E(ii,jj)-map2(ii,jj)
                     endif
                   enddo
                   enddo
                   if(associated(cc) .and. cc/=0) then
                     weight = rec%weight
                     if(associated(obs%sd)) nullify(obs%sd)
                     n= floor(cc); m = nint((cc - DBLE(n))*10)
                     msd = 0
                   endif
            end select
            
            !% data range
            !if ~isempty(cc)
            if(associated(cc)) then
                !if isfield(obs,'sd') && ~isempty(obs.sd)
                if(associated(obs%sd)) then
                    sd = int(obs%sd)     !% evaluation period, datenum
                    ed = int(w%ModelEndTime)         !% evaluation period end
                else
                    sd = nint(max(obs%t(1),sim%t(1)))
                    ed = nint(min(obs%t(size(obs%t)),sim%t(size(sim%t))))
                endif
                
                sd = max(msd,sd)
                !% method
                select case (m)
                    case (0,1,2)    !% {nash,rmse,rnash}
                        ![s,comp1]=compTS(obs,sim,sd+1,ed,1);
                        call compTS(obs,sim,sd+1,ed,1,s)
                        !comp(i) = comp1;
                        res(1) = s(m+1)
                        !if m==0 | m == 2
                        if (m==0 .or. m==2) then
                            res(1) = 1- res(1)    !% nash or rnash
                        endif
                        !stats(i,:)=s;
                    case (4)
                        ![s,comp1]=compTS(obs,sim,sd+1,ed,1);
                        call compTS(obs,sim,sd+1,ed,1,s)
                        !comp(i) = comp1;
                        res(1) = (s(1)+s(3)) / 2
                        res(1) = 1 - res(1)    !% (nash+rnash)/2
                        res(1) = min(res(1),1.5D0)
                        !stats(i,:)=s;
                    case (5)        !% nash with mean removed
                        ![s,comp1]=compTS(obs,sim,sd+1,ed,2);
                        call compTS(obs,sim,sd+1,ed,2,s)
                        !comp(i) = comp1;
                        res(2) = 1 - s(1)
                        res(2) = min(res(2),1.5D0)    !% nash don't go too big and interfere
                        !stats(i,:)=s;
                    case (6)    !soil moisture
                        call compTS(obs,sim,sd+1,ed,3,s)  
                        res(2) = 1-s(1)
                    case (7)    !grace
                        call compTS(obs,sim,sd+1,ed,1,s)  
                        res(2) = s(5)
                end select
                !if isnan(res) || isempty(res) || (sim.t(end)-sim.t(1))/(ed-sd)<0.8 
                
                ns = size(sim%t)
                !if(ns>0) then
                !  if((sim%t(ns)-sim%t(1))/(ed-sd)<0.8) then
                !    res = 2000.0D0
                !  endif
                !else
                !  res = 2000.0D0
                !endif
                !objfun = objfun + weight * res
                !call display( 'objfun',i,' w',weight,' res',res)
            endif
        enddo
        
!        ET_obs(:) = [17.43D0, 18.23D0, 39.68D0,	36.72D0, 59.59D0, 81.58D0, 89.91D0, 83.84D0, 57.58D0, 31.51D0, 19.78D0, 16.42D0, &
!                 &   16.76D0, 18.07D0, 37.61D0, 45.04D0, 54.72D0, 86.01D0, 107.94D0, 84.17D0, 52.06D0, 25.23D0, 22.09D0, 17.38D0, &
!                 &   13.51D0, 22.16D0, 37.88D0, 41.94D0, 56.65D0, 78.40D0, 87.88D0,  88.51D0, 59.46D0, 29.87D0, 22.45D0, 16.96D0, &
!                 &   16.82D0, 20.46D0, 32.83D0,	36.06D0, 62.63D0, 92.36D0, 105.19D0, 91.08D0, 42.43D0, 23.34D0,	19.65D0, 13.48D0, &
!                 &   16.27D0, 17.92D0, 34.07D0,	46.23D0, 64.42D0, 81.34D0, 91.61D0,	95.08D0, 55.95D0, 27.61D0,	22.11D0, 19.56D0, &
!                 &   20.59D0, 25.93D0, 39.92D0,	36.69D0, 61.38D0, 71.21D0, 86.86D0,	82.21D0, 46.40D0, 29.01D0,	23.22D0, 18.67D0]
        
        ET_obs(:) = [16.76D0, 18.07D0, 37.61D0, 45.04D0, 54.72D0, 86.01D0, 107.94D0, 84.17D0, 52.06D0, 25.23D0, 22.09D0, 17.38D0, &
                 &   13.51D0, 22.16D0, 37.88D0, 41.94D0, 56.65D0, 78.40D0, 87.88D0,  88.51D0, 59.46D0, 29.87D0, 22.45D0, 16.96D0, &
                 &   16.82D0, 20.46D0, 32.83D0,	36.06D0, 62.63D0, 92.36D0, 105.19D0, 91.08D0, 42.43D0, 23.34D0,	19.65D0, 13.48D0, &
                 &   16.27D0, 17.92D0, 34.07D0,	46.23D0, 64.42D0, 81.34D0, 91.61D0,	95.08D0, 55.95D0, 27.61D0,	22.11D0, 19.56D0, &
                 &   20.59D0, 25.93D0, 39.92D0,	36.69D0, 61.38D0, 71.21D0, 86.86D0,	82.21D0, 46.40D0, 29.01D0,	23.22D0, 18.67D0, &
                 &   17.43D0, 18.23D0, 39.68D0,	36.72D0, 59.59D0, 81.58D0, 89.91D0, 83.84D0, 57.58D0, 31.51D0, 19.78D0, 16.42D0]
        
        
        
        ET_obsbar = sum(ET_obs(6:72))/(size(ET_obs)-5)
        ET_NASH = 1- (1 - sum((ET_obs(6:72)-ET_sim(6:72))**2)/sum((ET_obs(6:72)-ET_obsbar)**2));
        res(3) =  ET_NASH  ! 0.2 is the wight for ET
        res(3) = min(res(3),1.5D0)
        !
        !eid = fopen(rfile,'wt');
        eid = 50
        open(unit=eid,file=rfile,status='UNKNOWN',action='WRITE')
        !fprintf(eid,'%10g',objfun); % 1-Nash
        write(eid,*) res(1)
        write(eid,*) res(2)
        write(eid,*) res(3)
        !fclose(eid);
        close(unit=eid)
        
		if(allocated(p)) deallocate(p)

        end subroutine simResults3
