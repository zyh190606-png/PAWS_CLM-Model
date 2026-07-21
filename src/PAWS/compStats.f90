!
        subroutine compStats(comp, stats)
        !function [stats]=compStats(comp)
        use Vdata
        implicit none
        
        type(comp_type) :: comp
        real*8, dimension(4) :: stats
        real*8, dimension(:), allocatable :: obs, sim, temp
        real*8 :: obsbar, NASH, RMSE, ME, N, NRNASH
        logical, dimension(:), allocatable :: loc
        integer :: i, k
        
        !obs = comp.obs;
        if(allocated(obs)) deallocate(obs)
        allocate(obs(size(comp%obs)))
        obs = comp%obs
        
        !sim = comp.sim;
        if(allocated(sim)) deallocate(sim)
        allocate(sim(size(comp%sim)))
        sim = comp%sim
        
        !loc = find(isnan(obs) | isnan(sim) | sim>100000 | obs<-90);
        if(allocated(loc)) deallocate(loc)
        allocate(loc(size(obs)))
        loc = .true.
        where (isnan(obs) .or. isnan(sim) .or. sim>100000 .or. obs<-90)
          loc = .false.
        end where
        k = 0
        do i = 1, size(loc)
          if(loc(i)) k = k + 1
        enddo
        if(allocated(temp)) deallocate(temp)
        allocate(temp(k))
        !obs(loc)=[]; 
        k = 0
        do i = 1, size(loc)
          if(loc(i)) then
            k = k + 1
            temp(k) = obs(i)
          endif
        enddo
        if(allocated(obs)) deallocate(obs)
        allocate(obs(k))
        obs = temp
        
        !sim(loc)=[];
        k = 0
        do i = 1, size(loc)
          if(loc(i)) then
            k = k + 1
            temp(k) = sim(i)
          endif
        enddo
        if(allocated(sim)) deallocate(sim)
        allocate(sim(k))
        sim = temp
                
        !obsbar=mean(obs); 
        obsbar = sum(obs)/size(obs)
        NASH = 1 - sum((obs-sim)**2)/sum((obs-obsbar)**2);
        !RMSE=sqrt(mean((obs-sim).^2));
        RMSE = dsqrt( sum((obs-sim)**2) / size(obs) )
        !X = [ones(size(sim)) sim];
        ![b,bint,r,rint,stats] = regress(obs,X); % stats:[R2,F,p,errVar]

        !ME = mean(sim -  obs);
        ME = sum(sim-obs) / size(obs)
        !% now indices: Nroot NASH:
        !N = 1/2;
        N = 1.0D0 / 2.0D0
        !obs = (comp.obs).^N;
        if(allocated(obs)) deallocate(obs)
        allocate(obs(size(comp%obs)))
        WHERE(comp%obs<0) comp%obs=0;
        obs = comp%obs ** N
        
        !sim = (comp.sim).^N;
        if(allocated(sim)) deallocate(sim)
        allocate(sim(size(comp%sim)))
        WHERE(comp%sim<0) comp%sim=0;
        sim = comp%sim ** N

        !loc = find(isnan(obs) | isnan(sim) | sim>100000);
        if(allocated(loc)) deallocate(loc)
        allocate(loc(size(obs)))
        loc = .true.
        where (isnan(obs) .or. isnan(sim) .or. sim>100000)
          loc = .false.
        end where
        k = 0
        do i = 1, size(loc)
          if(loc(i)) k = k + 1
        enddo
        if(allocated(temp)) deallocate(temp)
        allocate(temp(k))
                        
        !obs(loc)=[];
        k = 0
        do i = 1, size(loc)
          if(loc(i)) then
            k = k + 1
            temp(k) = obs(i)
          endif
        enddo
        if(allocated(obs)) deallocate(obs)
        allocate(obs(k))
        obs = temp
        
        !sim(loc)=[];
        k = 0
        do i = 1, size(loc)
          if(loc(i)) then
            k = k + 1
            temp(k) = sim(i)
          endif
        enddo
        if(allocated(sim)) deallocate(sim)
        allocate(sim(k))
        sim = temp
        
        !obsbar=mean(obs); 
        obsbar = sum(obs)/size(obs)      
        NRNASH = 1 - sum((obs-sim)**2)/sum((obs-obsbar)**2);


        !stats = [NASH,RMSE,NRNASH,ME, stats];
        stats = [NASH, RMSE, NRNASH, ME]

        if(allocated(obs)) deallocate(obs)
        if(allocated(sim)) deallocate(sim)
        if(allocated(temp)) deallocate(temp)
        if(allocated(loc)) deallocate(loc)
        
        end subroutine compStats