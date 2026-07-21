!
        subroutine compTS(tobs, tsim, sd, ed, opt, stats)
        !function [stats,comp]=compTS(tobs,tsim,sd,ed,varargin)
        use Vdata
        implicit none
        
        interface
            subroutine intersect(A,B,c,ia,ib)
            real*8, dimension(:) :: A, B
            integer, dimension(:), allocatable :: ia, ib
            real*8, dimension(:), allocatable :: c             
            end subroutine intersect    
        end interface
            
        type(sim_type) :: tobs, tsim
        type(comp_type) :: comp
        real*8, dimension(4) :: stats
        integer :: opt, sd, ed, k1, k2, N, i
        logical, dimension(:), allocatable :: loc
        real*8, dimension(:), allocatable :: tempt, tempv
        integer, dimension(:), allocatable :: ia, ib
        real*8, dimension(:), allocatable :: c             
        !if length(varargin)>0
            !opt = varargin{1};
        !else
            !opt = 0;
        !end
        !% assuming ts1 and ts2 are both consecutive
        !% missing data should be marked as NaN or Inf
        !%if opt ==0 || opt ==2
        if (opt ==0) then
            !k1 = find(round(tobs.t) == round(sd),1);
            k1 = 0
            do i = 1, size(tobs%t)
              if(nint(tobs%t(i))==sd) then
                k1 = i
                exit
              endif
            enddo
            !k2 = find(round(tobs.t) == round(ed),1);
            k2 = 0
            do i = 1, size(tobs%t)
              if(nint(tobs%t(i))==ed) then
                k2 = i
                exit
              endif
            enddo
            N = k2 - k1
            !comp.obs = tobs.v(k1:k2);
            if(allocated(comp%obs)) deallocate(comp%obs)
            allocate(comp%obs(size(tobs%v(k1:k2))))
            comp%obs = tobs%v(k1:k2)
            !k1 = find(round(tsim.t) == round(sd),1);
            k1 = 0
            do i = 1, size(tsim%t)
              if(nint(tsim%t(i))==sd) then
                k1 = i
                exit
              endif
            enddo
            !k2 = find(round(tsim.t) == round(ed),1);
            k2 = 0
            do i = 1, size(tsim%t)
              if(nint(tsim%t(i))==ed) then
                k2 = i
                exit
              endif
            enddo
            !if isempty(k1) || isempty(k2)
            if(k1==0 .or. k2==0) then
                !comp.sim = zeros(size(comp.obs));
                if(allocated(comp%sim)) deallocate(comp%sim)
                allocate(comp%sim(size(comp%obs)))
                comp%sim = 0.0D0
            else
                !comp.sim = tsim.v(k1:k1+N);
                if(allocated(comp%sim)) deallocate(comp%sim)
                allocate(comp%sim(size(tsim%v(k1:k1+N))))
                comp%sim = tsim%v(k1:k1+N)
            endif
        else
            !loc1=(tsim.t<sd | tsim.t>ed);
            if(allocated(loc)) deallocate(loc)
            allocate(loc(size(tsim%t)))
            loc = .true.
            where (tsim%t < sd .or. tsim%t > ed)
              loc = .false.
            end where
            !tsim.t(loc1)=[]; tsim.v(loc1)=[];
            N = 0
            do i = 1, size(loc)
              if(loc(i)) N = N + 1
            enddo
            if(allocated(tempt)) deallocate(tempt)
            allocate(tempt(N))
            if(allocated(tempv)) deallocate(tempv)
            allocate(tempv(N))
            N = 0
            do i = 1, size(loc)
              if(loc(i)) then
                N = N + 1
                tempt(N) = tsim%t(i)
                tempv(N) = tsim%v(i)
              endif
            enddo
            if(allocated(tsim%t)) deallocate(tsim%t)
            allocate(tsim%t(N))
            tsim%t = tempt
            if(allocated(tsim%v)) deallocate(tsim%v)
            allocate(tsim%v(N))
            tsim%v = tempv
            
            !loc2=(tobs.t<sd | tobs.t>ed);            
            if(allocated(loc)) deallocate(loc)
            allocate(loc(size(tobs%t)))
            loc = .true.
            where (tobs%t < sd .or. tobs%t > ed)
              loc = .false.
            end where
            !tobs.t(loc2)=[]; tobs.v(loc2)=[];
            N = 0
            do i = 1, size(loc)
              if(loc(i)) N = N + 1
            enddo
            if(allocated(tempt)) deallocate(tempt)
            allocate(tempt(N))
            if(allocated(tempv)) deallocate(tempv)
            allocate(tempv(N))
            N = 0
            do i = 1, size(loc)
              if(loc(i)) then
                N = N + 1
                tempt(N) = tobs%t(i)
                tempv(N) = tobs%v(i)
              endif
            enddo
            if(allocated(tobs%t)) deallocate(tobs%t)
            allocate(tobs%t(N))
            tobs%t = tempt
            if(allocated(tobs%v)) deallocate(tobs%v)
            allocate(tobs%v(N))
            tobs%v = tempv
            
            if (size(tobs%t)==0) then
                stats = [1,1,1,1] !no result fit in record period. Didn't finished. Just for debug. 
                return
            endif
                        
            ![c,ia,ib]=intersect(tsim.t,tobs.t);
            call intersect(tsim%t, tobs%t, c, ia, ib)
            k1 = maxval(ib)
            k2 = minval(ib)
            if (k2<0 .OR. k1>1D8) THEN
              call display("Warning: compTS detects simulated values may not be consecutive. CompTS may die")
            endif
            
            !comp.sim = tsim.v(ia);
            if(allocated(comp%sim)) deallocate(comp%sim)
            allocate(comp%sim(size(tsim%v(ia))))
            comp%sim = tsim%v(ia)
            !comp.obs = tobs.v(ib);
            if(allocated(comp%obs)) deallocate(comp%obs)
            n = size(tobs%v(ib))
            allocate(comp%obs(n))
            comp%obs = tobs%v(ib)
            !comp.ts(1)=tobs;
            if(allocated(comp%ts(1)%t)) deallocate(comp%ts(1)%t)
            allocate(comp%ts(1)%t(size(tobs%t)))
            comp%ts(1)%t = tobs%t
            if(allocated(comp%ts(1)%v)) deallocate(comp%ts(1)%v)
            allocate(comp%ts(1)%v(size(tobs%v)))
            comp%ts(1)%v = tobs%v
            if(associated(comp%ts(1)%sd)) deallocate(comp%ts(1)%sd)
            allocate(comp%ts(1)%sd)
            if (associated(tobs%sd) ) comp%ts(1)%sd = tobs%sd
            !comp.ts(2)=copyStruct(tsim,tobs,0);
            if(allocated(comp%ts(2)%t)) deallocate(comp%ts(2)%t)
            allocate(comp%ts(2)%t(size(tsim%t)))
            comp%ts(2)%t = tsim%t
            if(allocated(comp%ts(2)%v)) deallocate(comp%ts(2)%v)
            allocate(comp%ts(2)%v(size(tsim%v)))
            comp%ts(2)%v = tsim%v
            !if(associated(comp%ts(2)%sd)) deallocate(comp%ts(2)%sd)
            !allocate(comp%ts(2)%sd)
            !comp%ts(2)%sd = tsim%sd
            
            if(size(c)>0) then
                k1 = 1; k2 = 1
            else
                k1 = 0; k2 = 0
            endif
        endif
        if (opt==2) then        !% remove mean
            !%comp.obs = comp.obs - min(comp.obs);
            !%comp.sim = comp.sim - min(comp.sim);
            !comp.sim = comp.sim - (mean(comp.sim)-mean(comp.obs));
            comp%sim = comp%sim - (sum(comp%sim)/size(comp%sim) - sum(comp%obs)/size(comp%obs))
        endif
        if (opt==3) then !grace
            comp%sim = comp%sim - (sum(comp%sim)/size(comp%sim))        
            comp%obs = comp%obs - (sum(comp%obs)/size(comp%obs))
        end if

        ![stats]=compStats(comp);
        stats = 0.0D0
        call compStats(comp, stats)
        !if isempty(k1) || isempty(k2)
        if (k1==0 .or. k2==0) then
            stats(:) = 2000.0D0
            stats([1,3]) = -2000.0D0
        endif

        if(allocated(loc)) deallocate(loc)
        if(allocated(tempt)) deallocate(tempt)
        if(allocated(tempv)) deallocate(tempv)
        if(allocated(ia)) deallocate(ia)
        if(allocated(ib)) deallocate(ib)
        if(allocated(c)) deallocate(c)
        
        end subroutine compTS

!======================================================================================
        subroutine intersect(A,B,c,ia,ib)
        use displayMod
        implicit none
        real*8, dimension(:) :: A, B
        integer, dimension(:), allocatable :: ia, ib
        real*8, dimension(:), allocatable :: c
        integer :: i1, i2, nA, nB, n
        real*8, parameter :: eps = 1D-15
        real*8 dt
        logical, dimension(:), allocatable :: lA, lB
        
        nA = size(A)
        nB = size(B)
        
        if(.not. allocated(lA)) allocate(lA(nA))
        if(.not. allocated(lB)) allocate(lB(nB))
        lA = .false.
        lB = .false.
        
        dt = B(size(B))-B(size(B)-1)
        dt = 1.1*dt;
        n = 0
        
        do i1 = 1, nA
          do i2 = n, nB
            if(abs(A(i1)-B(i2)) <= eps) then
            if (la(i1) .or. lb(i2)) then
               call display("intersect:: still find repetitive values?")
            endif
              lA(i1) = .true.
              lB(i2) = .true.
              exit
            endif
          enddo
        enddo
        
        n = 0
        do i1 = 1, nA
          if(lA(i1)) then
            n = n + 1
          endif
        enddo
        if(.not. allocated(ia)) allocate(ia(n))
        if(.not. allocated(ib)) allocate(ib(n))
        if(.not. allocated(c)) allocate(c(n))
        
        n = 0
        do i1 = 1, nA
          if(lA(i1)) then
            n = n + 1
            ia(n) = i1
            c(n) = A(i1)
          endif
        enddo

        n = 0
        do i2 = 1, nB
          if(lB(i2)) then
            n = n + 1
            ib(n) = i2
          endif
        enddo    
        
        end subroutine intersect