!

      subroutine columnOutput(cmd)
      !% (E) function to grab output from a land column for a time step.-
      !% This is for output purposes, working with wtrshd_column, also wtrshd_day.
      !% Include fluxes and states
      !global g w
      use Vdata
      implicit none
      integer :: cmd, k
      real*8 :: THE, SWE, h, ET, PET, SUBLM, Inf, Rch, Prcp, GW, GW2, num(12)
      integer :: i, j
      real*8, dimension(:), allocatable :: acc, output
      
      if (cmd == 0) then    !% record data
          do k = 1, size(w%Rec%robj)
              !if isempty(w.Rec.robj(k).acc)
              if(size(w%Rec%robj(k)%acc)>1) then
                w%Rec%robj(k)%acc = 0.0D0
              endif
              
              if (w%Rec%robj(k)%type == 'c') then
                  i = w%Rec%robj(k)%IX(1)
                  j = w%Rec%robj(k)%IX(2)
                  THE = sum(g%VDZ%THE(i,j,:) * g%VDZ%DZ(i,j,:))
                  SWE = g%Veg%SWE(i,j)
                  h = g%OVN%hc(i,j)
                  ET =  g%Veg%ET(i,j)
                  PET = g%Veg%PET(i,j)
                  SUBLM = g%Veg%SNOM(i,j)
                  Inf = g%VDZ%DF(i,j)
                  Rch = g%VDZ%Dperc(i,j)
                  Prcp = g%OVN%prcp(i,j)
                  GW = g%GW(1)%h(i,j)
                  GW2 = g%GW(2)%h(i,j)
                  acc = w%Rec%robj(k)%acc
                  !num = [w.t h,THE,SWE,num(4)+ET,num(5)+PET,num(6)+SUBLM,num(7)+Inf,num(8)+Rch,num(9)+Prcp GW GW2];
                  num = (/w%t,h,THE,SWE,acc(4)+ET,acc(5)+PET,acc(6)+SUBLM,acc(7)+Inf,acc(8)+Rch,acc(9)+Prcp,GW, GW2/)
                  if(allocated(output)) deallocate(output)
                  allocate(output(12))
                  output = num
                  w%Rec%robj(k)%acc = num
              endif
          enddo
      else    !% output and clear data
          !error('need more work here')
          call display( 'error! (need more work here)' )
          if(allocated(output)) deallocate(output)
          output = w%Rec%robj(1)%acc
          !w.Rec.robj(k).acc = [];
      endif
      
      if(allocated(acc)) deallocate(acc)
      if(allocated(output)) deallocate(output)
      end subroutine columnOutput
