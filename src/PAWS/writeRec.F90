!
      subroutine writeRec
      ! (E) write accumulated fluxes at monitoring points to file, and clear value
      use Vdata
      implicit none
      real*8, dimension(size(w%Rec%robj)) :: output, den
      integer i, fid
      type(robj_type), pointer :: s
      character*500 :: fm
      
      !if isfield(w.Rec,'robj') && ~isempty(w.Rec.robj)
      if(associated(w%Rec%robj) .and. size(w.CR)>0) then
        output = 0.0d0
        den = output
        
        do i = 1, size(output)
            s => w%Rec%robj(i)
            !write(*,*) 'S%Type',s%type
            if(TRIM(s%type) == TRIM('r')) then
            	!write(*,*) 1, i, s%acc(1)
                output(i) = s%acc(1)/(s%accT*86400)        ! River discharge
            elseif (TRIM(s%type) == TRIM('g')) then
                output(i) = s%acc(1)
            elseif (TRIM(s%type) == TRIM('s')) then
                output(i) = s%acc(1)
            else
                output(i) = 0D0 
                !write(*,*) 2, i, s%acc(1)
            end if
            nullify(s)
            w%Rec%robj(i)%acc = 0.0D0
            if(associated(w%Rec%robj(i)%accT)) w%Rec%robj(i)%accT = 0.0D0
        end do
        fid = Env%recfid(2)
        ! fm = ['%13.9g,',repmat('%8.6g,',1,length(output)-1),'%8.6g\n'];
        if (size(output)>0) then
        fm = repeat('1F13.6,',size(output)-1)
        fm = '(1F13.3,'//trim(fm)//'1F13.6)'
        ! fprintf(fid,fm,w.t,output);
        write(UNIT=fid,FMT=trim(fm)) w%t,output
        endif
      end if
      
      end subroutine writeRec    
