!
		subroutine wDaily(ts, opt)
		!function tsn = wDaily(ts,varargin)
		!% convert a regularly spaced time series (ts.t, ts.v) into daily ts
		!% time stamps are at the end of the time step
		!% dt = ts.t(2)-ts.t(1);
		!% if dt<1 % subdaily
		!%     if abs(dt - 1/24)<1e-5 % hourly
		!%         H = ts.t(1)-floor(ts.t(1)); 
		!%         
		!%         
		!%         
		!%     end
		!%     
		!% end
		use Vdata
		implicit none
		type(sim_type) :: ts
		integer :: opt, n, m, i, is, ie, nB, i1
		character*32 :: oper
		real*8 dt
		real*8, dimension(:), allocatable :: tn, vn

		!if length(varargin)==0
		if(opt==0) then
			!oper = @mean;
			oper = 'mean'
		else
			!oper = str2func(varargin{1});
		endif
		!tsn = ts;
		!n = length(ts.v); m = ceil(n/24);
		
		nB = size(ts%t)
		dt = ts%t(size(ts%t))- ts%t(size(ts%t)-1)
		dt = dt*1.25
        do i1 = nB,2,-1
          if (ts%t(i1)-ts%t(i1-1)<=0) then
              call display("wDaily:: Found anachronistic time stamps, Job renewed? Only use last portion")
              n = i1
              EXIT
          endif
          n = i1-1
        enddo
        
		i1 = n-1
		n = size(ts%v(i1:size(ts%v))); 
		m = ceiling(dble(n)/24.0D0)
		!tn = zeros(m,1); vn = zeros(m,1);
		if(allocated(tn)) deallocate(tn)
		allocate(tn(m))
		if(allocated(vn)) deallocate(vn)
		allocate(vn(m))
		tn = 0.0D0; vn = 0.0D0
		do i = 1, m
			is = (i-1)*24+i1+1; ie = min(i*24+i1,n)
			tn(i)=floor(ts%t(is))
			select case (trim(oper))
			  case ('mean')
			    vn(i) = sum(ts%v(is:ie))/size(ts%v(is:ie))
			end select
		enddo
		if(allocated(ts%t)) deallocate(ts%t)
		allocate(ts%t(size(tn)))
		ts%t = tn
		if(allocated(ts%v)) deallocate(ts%v)
		allocate(ts%v(size(vn)))
		ts%v = vn
		
		if(allocated(tn)) deallocate(tn)
		if(allocated(vn)) deallocate(vn)
		
		end subroutine wDaily