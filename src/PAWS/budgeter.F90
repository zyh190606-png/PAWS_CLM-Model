#include <preproc.h> 
    module budgeter_mod
    integer nitemw

    CONTAINS
    subroutine budgeter(obj,id,mask,Hin,m)
    !subroutine budgeter(obj,id,mask,Hin,varargin)
    ! (E) This function summarizes the hydrologic budget of the watershed in
    ! the past time step and writes relevant mass balance output
    ! Input:
    ! obj: a string indicating the type of budgeting being done. w for-
    ! watershed, r for river, OVN for surface runoff, pond for surface water
    ! Output:
    ! text file output
    use Vdata
    use memoryChunk, Only: IT,FINDEMPTY,tempStorage, TEMPZ => TEMP
#ifndef MEX      
    USE clmtype ! to kick out clm variables
#endif
    implicit none
    character*32 obj
    integer Hin
    real*8, dimension(w%m(1),w%m(2)) :: mask
    integer:: varargin=1,m
    integer mode, i, j, id(m), fid, nn,L, ii, jj, k
    real*8 TA, TAV, N, Qgc, Qoc, rh, prcp, SMELT, Qout, Dperc, Dperc2, Inf
    real*8 ET, PET, EvapCS, SUBM, Tp, EvapG, NewS, Ex, Gprcp, driftOut,ExC,GWBCOut,highGW
    real*8 hd, h, hf, hc, hI, THE, Rf, evap, ICE, DEW, CSNEW
    real*8 GW, CS, SWE, Tmp, Rad, Mass, test1, test2, test3, test4, test5, botT
    real*8 dt, S, Area, Qin, trib, tt, kkk, test, SatE, InfE, TL, SNOM
    real*8 var2D, meanTHE, npp, gpp, ar, hr, SatFrac, depH
    real*8, dimension(:), allocatable :: T, output
    real*8, dimension(size(w%CR)) :: rL, dat
    integer c
    real*8, dimension(size(g%OVN%h,1),size(g%OVN%h,2)) :: hf2D,the2d
    real*8, dimension(size(g%OVN%hc,1),size(g%OVN%hc,2)) :: hc2D
    real*8, dimension(size(g%Veg%hI,1),size(g%Veg%hI,2)) :: hI2D
    real*8, dimension(size(g%OVN%Rf,1),size(g%OVN%Rf,2)) :: h2D, test22D
    real*8, dimension(size(g%Veg%Gprcp,1),size(g%Veg%Gprcp,2)) :: Gprcp2D
    real*8, dimension(size(g%Veg%EvapG,1),size(g%Veg%EvapG,2)) :: EvapG2D
    real*8, dimension(size(g%VDZ%DF,1),size(g%VDZ%DF,2)) :: Inf2D
    integer,dimension(:),pointer::idx=>NULL(),Rep=>NULL(),Cidx=>NULL()
    Real*8,dimension(:),pointer::temp1d=>NULL(),QinRiv=>null()
    Real*8,dimension(:,:),pointer::GWBCflux=>NULL(),highGWFix=>null(),prcpI=>null(),evapI=>null(),FrozenInf=>null(),Combo=>null()
    Real*8,dimension(:,:),pointer::RfF=>NULL(),RfI=>null(),hfSave=>null(),Soil2AirDirect=>null(),iceMove=>null(),hfd0=>null()
    Real*8,dimension(:,:),pointer::EvapSoil=>null(),the_sumSave=>null(),irrigation=>null()

    real*8, dimension(:,:), allocatable :: sum_temp
    real*8, dimension(size(g%VDZ%th%Soutv,1),size(g%VDZ%th%Soutv,2),size(g%VDZ%th%Soutv,3))::Soutv
    real*8, dimension(size(g%VDZ%THE,1),size(g%VDZ%THE,2),size(g%VDZ%THE,3)) :: temp

    real*8, dimension(:,:), allocatable :: mask1, mask2, nmask2
    integer :: date(3), yyyy, mm, dd, time
    real*8 err,m_err
    real*8, save :: col_mass, THE0,CS0,h0,SWE0,m_err2
    real*8, save :: hc0=0D0,hf0=0D0,hd0=0D0,hI0=0D0
    real*8, save :: rf1,rf2,evapIm,prcpIm,RfIm,RfFm,Rfm,FrozenInfm,Combom,Soil2AirDirectm,iceMovem,EvapSoilm,irrigationm
    real*8, save :: t1,t2,t3,t4,t5
    real*8, pointer, save:: col_mass_s(:,:)
    integer, save:: TNSTATEI = 0
    integer,dimension(size(g%OVN%h,1)*size(g%OVN%h,2)) :: filter
    integer,save:: nmap = 0
    integer msize(4), mpoint(2)
    real*8 test2I, test2c, test2fd, test2fdo

    msize(1) = size(g%VDZ%h,1)
    msize(2) = size(g%VDZ%h,2)
    msize(3) = size(g%VDZ%h,3)
    msize(4) = size(g%Veg%Frac,3)


    !if length(varargin)>0
    if(varargin>0) then
        mode = varargin
        ! mode == 1, normal run; mode==2, surface only
    else
        mode = 1
    end if

    select case (trim(obj))
    case ('w')
        ! state variables: h,rh,GWh,THE,SWE,CS, (all normalize to valid
        ! area)
        ! fluxes: external[prcp, Qout,Dperc2], internal(Inf, ET, Qoc, Qgc, SNOM)
        TA = w%Stats%TA
        TAV= w%Stats%TAV
        ! here the meaning of H is the fractional step in an hour
        N = dnint(g%Veg%dt/g%OVN%dt)
        rL = w%CR
        Qgc = 0.0D0
        Qoc = 0.0D0
        Qin = 0D0
        rh = 0.0D0
        !if isfield(w.DM,'mask2') && ~isempty(w.DM.mask2)
        if(associated(w%DM%mask2)) then
            if(.not. allocated(mask2)) then
                allocate(mask2(msize(1),msize(2)) )
            end if
            mask2 = w%DM%mask2
        else
            if(.not. allocated(mask2)) then
                allocate(mask2(msize(1),msize(2)) )
            end if
            mask2 = 0.0D0
            where(w%DM%mask)
                mask2 = 1.0D0
            end where
        end if

        if(.not. allocated(mask1)) then
            allocate(mask1(msize(1),msize(2)) )
        endif
        mask1 = 0.0D0
        where(w%DM%mask)
            mask1 = 1.0D0
        end where

        if (nmap .eq. 0) then
            do i = 1, msize(1)
                do j = 1, msize(2)
                    if (w%DM%mask(i,j)) nmap = nmap + 1
                enddo
            enddo
        endif

        !if abs(r(rL(end)).Riv.t - w.t)<1e-5
        ! fluxes
        ExC = 0
        temp1d => oneD_ptr(msize(1)*msize(2),g%OVN%Ex)
        call tempStorage('QinRiv',QinRiv, size(w%CR))
        do i = 1, size(rL)
            c = rL(i)
            !Qgc = Qgc + sum(r(c).Riv.fG(2:end-1)); % Qgc is cleared after each unconfined groundwater step
            Qgc = Qgc + sum(r(c)%Riv%fG(2:size(r(c)%Riv%fG)-1))/TA
            rh = rh + sum(r(c)%Riv%Ac(2:size(r(c)%Riv%Ac)-1))*r(c)%DM%d(1)
            dat(i) = Qgc

            idx => g%Riv%Lidx(i)%ptr
            ExC = ExC + sum(temp1d(idx))/(msize(1)*msize(2))
            if (associated(r(c)%Net%UPS)) then
                !if (r(c)%Net%UPS > 1D6) Qin = Qin + r(c)%Riv%Qx(1) * w%dt * 86400D0/TA
                if (r(c)%Net%UPS > 1D6) Qin = Qin + QinRiv(i)/TA
            endif
        end do
        !Qgc = Qgc/TA;
        rh = rh/TA
        Qoc = sum(w%Rec%acc%Qoc)/TA
        w%Rec%acc%Qoc(:) = 0.0
        if(.not. associated(g%VDZ%SMR)) then
            allocate(g%VDZ%SMR(size(g%OVN%prcp,1),size(g%OVN%prcp,2)))
        endif
        g%VDZ%SMR = g%OVN%prcp * mask1 * w%dt
        !prcp = mean(g.VDZ.SMR(:))
        prcp = sum(g%VDZ%SMR)/size(g%VDZ%SMR)
        !T=g.Wea.day.prcp(H,:); prcp2 = mean(T(T>-90)); % comment this line out if uses actually applied prcp#########

        !SatE = meanALL(g.OVN.Rf.*g.VDZ.R(:,:,2));
        SatE = sum(g%OVN%Rf * g%VDZ%R(:,:,2))/size(g%OVN%Rf)
        SatFrac = sum(g%VDZ%R(:,:,2))/size(g%OVN%Rf)
        depH = sum(g%OVN%hd)/size(g%OVN%hd)
        !InfE = meanALL(g.OVN.Rf.*(1-g.VDZ.R(:,:,2)));
        InfE = sum(g%OVN%Rf * (1-g%VDZ%R(:,:,2)))/size(g%OVN%Rf)
        !Soutv= g.VDZ.th.Soutv;

        Soutv = g%VDZ%th%Soutv
        !SMELT= meanALL(max(Soutv(:,:,7)-Soutv(:,:,1),0));
        if(allocated(sum_temp)) deallocate(sum_temp)
        allocate(sum_temp(size(Soutv,1),size(Soutv,2)))
        sum_temp = Soutv(:,:,7)-Soutv(:,:,1)
        where (sum_temp<=0)
            sum_temp = 0
        endwhere
        SMELT = sum(sum_temp)/size(Soutv(:,:,1))
        !SMELT = meanALL(Soutv(:,:,22));
        !            if (size(w.CR)>0) then
        !              Qout = w%Rec%acc%Qout(size(w%Rec%acc%Qout))/TA
        !            else
        !              Qout = 0D0
        !            endif

        Qout = 0D0
        do i = 1, size(rL)
            c= rL(i)
            if (.not. associated(r(c)%Net%Parent) ) THEN
                ! This MB check does not consider the case if a river flows into a overland pond/lake,etc
                Qout = Qout + w%Rec%acc%Qout(i)/TA
            endif
        enddo

        w%Rec%acc%Qout(:) = 0D0
        w%Rec%acc%trib(:) = 0D0
        QinRiv = 0D0
        !Dperc2= mean(g.GW(1).Qperc(:));
        Dperc2 = sum(g%GW(1)%Qperc)/size(g%GW(1)%Qperc)
        !Inf = mean(g.VDZ.DF(:))     ! the variables that are zero
        Inf = sum(g%VDZ%DF)/size(g%VDZ%DF)
        !Dperc = mean(g.VDZ.Dperc(:));
        Dperc = sum(g%VDZ%Dperc)/size(g%VDZ%Dperc)
        !ET  = mean(g.Veg.ET(:));
        ET = sum(g%Veg%ET)/size(g%Veg%ET)
        !PET  = mean(g.Veg.PET(:));
        PET = sum(g%Veg%PET)/size(g%Veg%PET)
        !EvapCS= mean(g.Veg.EvapCS(:)); % This is part of ET
        EvapCS = sum(g%Veg%EvapCS)/size(g%Veg%EvapCS)
        !SUBM=mean(g.Veg.SNOM(:));
        SUBM = sum(g%Veg%SNOM)/size(g%Veg%SNOM)
        !Tp = mean(g.Veg.Tp(:));
        Tp = sum(g%Veg%Tp)/size(g%Veg%Tp)
        !EvapG = mean(g.Veg.EvapG(:));
        EvapG = sum(g%Veg%EvapG)/size(g%Veg%EvapG)
        !NewS = mean(g.Veg.NewS(:));
        NewS = sum(g%Veg%NewS)/size(g%Veg%NewS)
        !Ex = meanALL(g.OVN.Ex);
        Ex = sum(g%OVN%Ex)/size(g%OVN%Ex)
        !%SNOW = mean(g.Veg.SNOW(:));
        if (mode == 1) then
            !Gprcp = mean(g.Veg.Gprcp(:));
            Gprcp = sum(g%Veg%Gprcp)/size(g%Veg%Gprcp)
        elseif (mode ==2) then
            !Gprcp = w.Stats.tGprcp;
            Gprcp = sum(w%Stats%tGprcp)/size(w%Stats%tGprcp)
        end if
        !% state variables
        !driftOut = meanALL(g.OVN.h.*(~mask2));
        if(allocated(nmask2)) deallocate(nmask2)
        allocate(nmask2(size(mask2,1),size(mask2,2)))
        nmask2 = 0.0D0
        where(mask2<0.1)
            nmask2 = 1.0D0
        endwhere
        driftOut = sum(g%OVN%h * nmask2)/size(g%OVN%h)
        !g.OVN.h(~mask2)=0;
        ice = 0
        do i = 1, size(g%OVN%h,1)
            do j = 1, size(g%OVN%h,2)
                !if(.not. mask2(i,j)) then
                if(mask2(i,j)<0.1) then
                    g%OVN%h(i,j) = 0
                end if
                ice = ice + sum(g%VDZ%ICE(I,J,:))/1000D0
            end do
        end do
        !if isnan(SatE)
        !    4;
        !end
        !hf = mean(g.OVN.h(:)); hc  = mean(g.OVN.hc(:)); hI = meanALL(g.Veg.hI);
        dew = sum(g%Veg%DEW)/size(g%Veg%DEW)
        hf = sum(g%OVN%h)/size(g%OVN%h)
        hc = sum(g%OVN%hc)/size(g%OVN%hc)
        hI = sum(g%Veg%hI)/size(g%Veg%hI)
        !hd = meanALL(g.OVN.hd.*g.OVN.dP(:,:,4));
        hd = sum(g%OVN%hd * g%OVN%dP(:,:,4))/size(g%OVn%hd)
        h = hf + hc + hI + hd
        !temp = g.VDZ.THE.*g.VDZ.DZ; temp = sum(temp,3);
        temp = g%VDZ%THE * g%VDZ%DZ
        if(allocated(sum_temp)) deallocate(sum_temp)
        allocate(sum_temp(size(temp,1),size(temp,2)))
        sum_temp = sum(temp, 3) ! vertically integrated THE
        !THE = mean(temp(:));

        ! THIS portion is for statistics!!!!!!!!
        the2d = sum_temp/(g%Topo%E - g%GW(1)%EB)
        meanTHE = sum(the2d,mask = w%DM%mask)/nmap
        var2d   = cvMap(the2d,w%DM%mask,size(w%DM%mask),meanTHE)

        THE = sum(sum_temp)/size(sum_temp) + ice/size(w%DM%mask)
        !if isfield(g.GW(1),'STO') && g.GW(1).STO(20,20)>0
        if(allocated(sum_temp)) then
            deallocate(sum_temp)
        end if
        allocate(sum_temp(size(w%DM%mask,1),size(w%DM%mask,2)))
        mpoint = (/(msize(1)+1)/2 ,(msize(2)+1)/2 /)
        if ((associated(g%GW(1)%STO)) .and. (g%GW(1)%STO(mpoint(1),mpoint(2))>0)) then
            !temp = (g.GW(1).STO).*w.DM.mask;
            sum_temp = g%GW(1)%STO * mask1
        else
            !temp = g.GW(1).ST.*(g.GW(1).h-g.GW(1).EB).*w.DM.mask;
            sum_temp = g%GW(1)%ST * (g%GW(1)%h - g%GW(1)%EB) * mask1
        end if
        !GW = mean(temp(:));
        GW = sum(sum_temp)/size(sum_temp)
        call tempStorage('GWBCFlux',GWBCflux,1)
        GWBCOut = sum(GWBCflux)/size(GWBCflux)
        call tempStorage('highGWFix',highGWFix,1)
        highGW  = sum(highGWFix)/size(highGWFix)
        call tempStorage('FrozenInf',FrozenInf,1)
        FrozenInfm  = sum(FrozenInf)/size(FrozenInf)
        call tempStorage('Soil2AirDirect',Soil2AirDirect,1)
        Soil2AirDirectm  = sum(Soil2AirDirect)/size(Soil2AirDirect)
        call tempStorage('Combo',Combo,1)
        Combom  = sum(Combo)/size(Combo)
        call tempStorage('iceMove',iceMove,1)
        iceMovem  = sum(iceMove)/size(iceMove)
        call tempStorage('EvapSoil',EvapSoil,1)
        EvapSoilm = sum(EvapSoil)/size(EvapSoil)
        call tempStorage('hfSave',hfSave,1)
        call tempStorage('irrigation',irrigation,1)
        irrigationm = sum(irrigation)/size(irrigation)

        !CS = mean(g.Veg.CS(:));
        CS = sum(g%Veg%CS)/size(g%Veg%CS)
        CSNEW = sum(g%Veg%CSNEW)/size(g%Veg%CSNEW)
        !SWE = mean(g.Veg.SWE(:));
        SWE = sum(g%Veg%SWE)/size(g%Veg%SWE)
        !T=g.Wea.day.temp(H,:); Temp = mean(T(T>-90));
        if(.not. allocated(T)) then
            allocate(T(size(g%Wea%day%temp(Hin,:))))
        end if
        T = g%Wea%day%temp(Hin,:)
        Tmp = 0.0D0
        nn = 0
        do i = 1, size(T)
            if(T(i)>-90.0D0) then
                Tmp = Tmp + T(i)
                nn = nn + 1
            endif
        enddo
        Tmp = Tmp / nn
        !T=g.Wea.day.rad(H,:); Rad = mean(T(T>-90));
        deallocate(T)
        allocate(T(size(g%Wea%day%rad(Hin,:))))
        T = g%Wea%day%rad(Hin,:)
        Rad = 0.0D0
        nn = 0
        do i = 1, size(T)
            if(T(i)>-90.0D0) then
                Rad = Rad + T(i)
                nn = nn + 1
            endif
        enddo
        Rad = Rad / nn


        !%[id,r(c).Riv.l,dt,Area,prcp,Qgc,Qoc,trib,Qin,Qout,S,test,test1];
        Mass = h + THE + GW + CS + SWE + rh
        test1 = Mass - (prcp +irrigationm-SUBM + Qin - Qout - Dperc2 - ET - Soil2AirDirectm - driftOut - GWBCOut) - w%Rec%Mem(1)

        ! TEMPORARY!!!! codes above are correct for whole system!!
        !test1 = Mass - (prcp -SUBM - Qout - Dperc2 - ET - driftOut) - w%Rec%Mem(1) + Dperc

        !%w.Rec.Mem(1) = Mass + Dperc2; % Memory are to remember states from last time step for mass balance check
        test2 = h - (Gprcp + Combom - Qoc - EvapG + EvapSoilm - Inf - driftOut + Ex + highGW - FrozenInfm) - w%Rec%Mem(2)
        !%w.Rec.Mem(2) = h;
        test3 = rh - (Qoc + Qgc + Qin  - Qout) - w%Rec%Mem(3)
        !%test= h-(Gprcp-EvapG-Inf-Qoc);
        !%w.Rec.Mem(3) = rh;
        !% test= Area - (prcp + Qoc + trib + Qin - Qout + Qgc);
        !% test1=Area - Qgc;
        test4 = THE - (Inf - EvapSoilm  + FrozenInfm + iceMovem - Soil2AirDirectm - Dperc - (ET - EvapG - EvapCS)) -w%Rec%Mem(4)
        !%w.Rec.Mem(4)= THE;
        test5 = GW - (Dperc - Dperc2 - Qgc - Ex - GWBCOut - highGW) -w%Rec%Mem(5)      !%############### PUMPING SOURCE TERMS

        ! analyze test2:
        IF (debugg) then
            call tempStorage('PrcpI',PrcpI,1)
            call tempStorage('the_sumSave',the_sumSave,1)
            call tempStorage('EvapI',EvapI,1)
            call tempStorage('RfI',RfI,1)
            evapIm = sum(EvapI)/size(EvapI)
            prcpIm = sum(prcpI)/size(prcpI)
            RfIm   = sum(RfI)/size(RfI)
            Rfm   = sum(g%OVN%Rf)/size(g%OVN%Rf)
            test2I = (hI - hI0) - (prcpIm - (RfIm + evapIm) )
            test2fd= (hf+hd-hf0-hd0) - (Rfm+Ex - (Qoc+driftOut))
            test2c = (hc-hc0)-(Gprcp+ Combom-prcpIm+highGW - (Inf + FrozenInfm + evapG -EvapSoilm - evapIm + Rfm - RfIm) )

            call tempStorage('hfd0',hfd0,1)
            t5 = sum(hfd0)/size(hfd0)
            test2fdo = (hf+hd-t5) -(Rfm - Qoc)
            If (abs(test2I)>1D-11 .or. abs(test2fd)>1D-12 .or. abs(test2c)>1D-12) then
                call display("MBE with test2I")
            Endif
            t5 = 0D0
            DO I=1,msize(1)
                DO J=1,msize(2)

                    IF (w%DM%mask(I,J)) THEN
                        t1 = g%OVN%hc(I,J) - hfSave(I,J)
                        t2 = g%Veg%Gprcp(I,J)+Combo(I,J)-prcpI(I,J)
                        t3 = g%VDZ%DF(I,J)+FrozenInf(I,J)+g%Veg%EvapG(I,J) -EvapSoil(I,J)-EvapI(I,J)+g%OVN%Rf(I,J)-RfI(I,J)
                        t4 = t1 - t2 + t3
                        if (abs(t4)>1D-11) THEN
                            call display("MBE with test2I, I:",I,"J: ",J)
                        endif
                        !

                        sum_temp(I,J) = 0D0
                        do k=1,g%VDZ%m(3)
                            sum_temp(I,J) = sum_temp(I,J)+g%VDZ%THE(I,J,K)*g%VDZ%DZ(I,J,K)+g%VDZ%ICE(I,J,K)/1000D0
                        enddo
                        !test4 = THE - (Inf - EvapSoilm  + FrozenInfm + iceMovem - Soil2AirDirectm - Dperc - (ET - EvapG - EvapCS)) -w%Rec%Mem(4)
                        t1 = (sum_temp(I,J)-the_sumSave(I,J) )
                        t2 = g%VDZ%DF(I,J)-EvapSoil(I,J)+FrozenInf(I,J)+iceMove(I,J)-Soil2AirDirect(I,J)-g%VDZ%Dperc(I,J)
                        t3 = g%Veg%ET(I,J)-g%Veg%EvapG(I,J)-g%Veg%EvapCS(I,J)
                        t4 = t1 -(t2-t3)
                        if (abs(t4)>1D-6) then
                            !   call display("MBE with test4 column, IJ:"//num2str(I)//num2str(J))
                        endif
                        t5 = t5 + t4
                        !

                    ENDIF
                ENDDO
            ENDDO


            !endif
            the_sumSave = sum_temp
            hfSave = g%OVN%hc
        ENDIF
        !%w.Rec.Mem(5) = GW + Dperc2; % Accumulative fluxes like Dperc2 need to be added to Mem t(n-1)
        !%w.Rec.Mem(6) = Dperc2;
        !botT = meanALL(g.VDZ.R(:,:,1));
        botT = sum(g%VDZ%R(:,:,1))/size(g%VDZ%R(:,:,1))
#ifndef MEX             
        npp = sum(clm3%g%l%c%ccf%pcf_a%npp)/size(clm3%g%l%c%ccf%pcf_a%npp)
        gpp = sum(clm3%g%l%c%ccf%pcf_a%gpp)/size(clm3%g%l%c%ccf%pcf_a%gpp)
        ar = sum(clm3%g%l%c%ccf%pcf_a%ar)/size(clm3%g%l%c%ccf%pcf_a%ar)
        hr = sum(clm3%g%l%c%ccf%hr)/size(clm3%g%l%c%ccf%hr)
#endif            
        !% temporary: test2 used to output Qout2, test3 to output
        !% farmlane
        if(allocated(output)) deallocate(output)
        allocate(output(nitemw))
        output=(/dble(id),w%dt,w%t,TA,TAV/TA,hf,hc,h,test2,rh,test3,THE,test4,GW,GWBCOut,test5,CS,SWE,Mass,test1,&
            &  Qout,Dperc2,prcp,Tmp,SatE,InfE,Ex,ExC,botT,Rad,Dperc,NewS,SMELT,SUBM,Gprcp,ET,Inf,Qoc,Qgc,PET,Tp,EvapCS,EvapG,&
            &     meanTHE,var2d,npp,gpp,ar,hr,SatFrac,depH,Irrigationm,FrozenInfm,Soil2AirDirectm,Mass/)
        ! line of statistics

        fid = Env%recfid(1)
        !fm = ['%s,%d,',repmat('%16.8g,',1,length(output)-2),'%16.8g\n'];
        !fprintf(fid,fm,datestr(r(c).Riv.t,'yyyymmdd,HHMMSS'),output);
        !date = datevec2(int(r(c)%Riv%t))
        date = datevec2(int(g%OVN%t))
        yyyy = date(1)
        mm = date(2)
        dd = date(3)
        !time = (r(c)%Riv%t - int(r(c)%Riv%t))*24
        time = (g%OVN%t - int(g%OVN%t))*24
        WHERE(abs(output)<1D-98) output = 0D0 ! otherwise we may have 0.4504139765909859-103 which may be incorrectly read into readBudget

        if (fid .ne. 0) then
            write(UNIT=fid,FMT='(I9,I9,1F6.0,60E24.16,E24.16)') yyyy*10000+mm*100+dd,time*10000,output
        endif
        !%fclose(fid);
        !%if ((abs(test1)>1e-7 .or. abs(test2)>1e-7) .and. test1<1)
        if ((abs(test1)>1e-6 .or. abs(test2)>1e-6 .or. abs(test3)>1e-8) .and. test1<1) then
            call display( 'budgeter Tests1-5 :: '// num2str(test1)//num2str(test2)//num2str(test3)//num2str(test4)//num2str(test5))
        endif
        if (abs(test4)>5e-7) then
            !write(*, *),3
        endif
        w%Rec%Mem(1:6) = (/Mass+Dperc2,h,rh,THE,GW+Dperc2,Dperc2/)
        CS0 = CS
        hc0 = hc;
        hf0 = hf
        hd0 = hd
        hI0 = hI

        !print debug
        if (Qgc > 0.01D0 .or. abs(PET)>1) then
            do i=1,g%VDZ%m(1)
                do j=1,g%VDZ%m(2)
                    if (abs(g%Veg%PET(i,j))>1 .or. g%Veg%SWE(i,j)<-0.1) then
                        call display( 'budgeter: suspicious')
                    endif
                enddo
            enddo
        endif
        !end if
    case ('r')
        !% mass balance equation:
        !% dA/dt = Qin-Qout+trib+Qoc+prcp-evap
        !% output=[id TA t dt(s) Qgc Qoc Qin Qout trib] all terms in unit of
        !% m2, mass/Length of River
        c = id(1)
        dt = r(c)%Riv%dt*86400
        TL = r(c)%Riv%l
        !S = mean(r(c).Riv.S(2:end-1)) * dt;
        S = sum(r(c)%Riv%S(2:(size(r(c)%Riv%S)-1)))/size(r(c)%Riv%S(2:(size(r(c)%Riv%S)-1))) * dt
        !Area = mean(r(c).Riv.Ac(2:end-1)); %[L2]
        Area = sum(r(c)%Riv%Ac(2:(size(r(c)%Riv%Ac)-1)))/size(r(c)%Riv%Ac(2:(size(r(c)%Riv%Ac)-1)))
        Qoc = sum(r(c)%Riv%Qoc(2:(size(r(c)%Riv%Qoc)-1)))/TL*dt
        Qgc = sum(r(c)%Riv%fG(2:(size(r(c)%Riv%fG)-1)))/TL
        Qin = r(c)%Net%UPSQ/TL*dt
        Qout= r(c)%Riv%Qx(size(r(c)%Riv%Qx)-1)/TL*dt
        !prcp= mean(r(c).Riv.prcp(2:end-1).*r(c).Riv.w(2:end-1))*dt;
        prcp = sum(r(c)%Riv%prcp(2:(size(r(c)%Riv%prcp)-1)) * r(c)%Riv%w(2:(size(r(c)%Riv%w)-1)))/ &
            &				 size(r(c)%Riv%prcp(2:(size(r(c)%Riv%prcp)-1)))*dt
        trib= sum(r(c)%Riv%trib(2:(size(r(c)%Riv%trib)-1)))/TL*dt
        tt=r(c)%Riv%t - floor(r(c)%Riv%t)
        if (abs(mod(tt-r(c)%Riv%dt,1.0D0/24.0D0))<1e-4) then     ! Qgc is fresh
            kkk = 0
            test = Area - (prcp + Qoc + trib + Qin - Qout + Qgc) - (w%Rec%MemR(c)+w%Rec%MemQgc(c))
        else
            test = Area - (prcp + Qoc + trib + Qin - Qout + Qgc) - w%Rec%MemR(c)
            kkk = 1
        end if
        w%Rec%MemR(c) = Area - Qgc
        !%         if (mod(tt,1/24)<1e-4) then
        !%             w%Rec%MemR(c) = Area
        !%         end if
        w%Rec%MemQgc(c) = Qgc
        !%S = r(c).Riv.S * dt;
        !%sumS = Qoc + prcp + trib;
        output=(/dble(id),r(c)%Riv%l,dt,Area,prcp,Qgc,Qoc,trib,Qin,Qout,S,test/)
        !fm = ['%s,%d,',repmat('%16.8g,',1,length(output)-2),'%16.8g\n'];
        !if abs(test)>1e-7 && kkk==1
        !  fid = fopen('rBudgeter.txt','at');
        !  fprintf(fid,fm,datestr(r(c).Riv.t,'yyyymmddHHMMSS'),output);
        !  fclose(fid);
        !end if
    case ('VDZ')

    case ('Veg')

    case ('GW')

    case ('OVN')
        !fid = fopen('gBudgeter.txt','at');
        !Inf = mean(g.VDZ.DF(:));
        Inf = sum(g%VDZ%DF)/size(g%VDZ%DF)
        !Gprcp = mean(g.Veg.Gprcp(:));
        Gprcp = sum(g%Veg%Gprcp)/size(g%Veg%Gprcp)
        !EvapG= mean(g.Veg.EvapG(:));
        EvapG = sum(g%Veg%EvapG)/size(g%Veg%EvapG)
        !h = mean(g.OVN.h(:));
        h = sum(g%OVN%h)/size(g%OVN%h)
        !SNOM = mean(g.Veg.SNOM(:));
        SNOM = sum(g%Veg%SNOM)/size(g%Veg%SNOM)
        TA = w%Stats%TA
        !T=g.Wea.day.temp(H,:); Temp = mean(T(T>-90));
        T = g%Wea%day%temp(Hin,:)
        !Tmp = sum(T(T>-90))/size(T(T>-90))
        nn = 0
        Tmp = 0.0D0
        do i = 1, size(T)
            if(T(i)>-90) then
                Tmp = Tmp + T(i)
                nn = nn + 1
            endif
        enddo
        Tmp = Tmp / nn
        !T=g.Wea.day.rad(H,:); Rad = mean(T(T>-90));
        T = g%Wea%day%rad(Hin,:)
        !Rad = sum(T(T>-90))/size(T(T>-90))
        nn = 0
        Rad = 0.0D0
        do i = 1, size(T)
            if(T(i)>-90) then
                Rad = Rad + T(i)
                nn = nn + 1
            endif
        enddo
        Rad = Rad / n
        Qoc = 0
        do j = 1, size(g%Riv%Qoc)
            Qoc = Qoc + sum(g%Riv%Qoc(j)%ptr) * (g%OVN%dt/TA)
        end do
        test = h-(Gprcp-EvapG-Inf-Qoc)
        output = (/dble(id),g%OVN%dt,TA,Gprcp,Temp,Rad,SNOM,EvapG,Inf,Qoc,h,test/)
        !fm = ['%s,%d,',repmat('%12g,',1,length(output)-2),'%12g\n'];
        !fprintf(fid,fm,datestr(g.OVN.t),output);
        !fclose(fid);
    case ('pond')
        !% matrix test (before flow and exchange with riv)
        !% activated only for debugging
        hf2D = g%OVN%h
        hc2D = g%OVN%hc
        hI2D = g%Veg%hI
        h2D = hf2D + hc2D + hI2D + g%OVN%Rf
        Gprcp2D = g%Veg%Gprcp
        EvapG2D = g%Veg%EvapG
        Inf2D = g%VDZ%DF
        test22D = h2D - (Gprcp2D - EvapG2D - Inf2D) - w%Rec%hSave
        m_err2  = maxval(abs(test22D))
        if (m_err2>1D-8) then
            !write(*,*), 2
        endif

        !if anyALL(abs(test2)>1e-7)
        !    4;
        !end
    case ('col')
        ! a column total mass conservation check
        ! produce the total column mass in a column, excluding groundwater
        ! including: hc:surface ponding layer, hf:flow domain, hI:impervious layer, hd:depression storage(the intercepted portion), CS:canopy storage, SWE:snow water, THE:soil water
        ! fluxes: prcp (forcing), evap, snow_out, Dperc2, Ex: exchange between groundwater and hd
        ! This is correct when Qoc is not considered, particularly, when wtrshd_col is used

        !!!!!!!!!!!!!!!! IN CASE YOU ARE WONDERING
        !!! IF THERE IS OVERLAND FLOW, THIS ACCOUNTING IS NOT COMPLETE!

        IF (TNSTATEI .eq. 0) THEN
            TNSTATEI = FINDEMPTY(IT,1,size(IT));
        ENDIF
        col_mass_s => TEMPZ(:,:,TNSTATEI);

        do ii = 1,size(id)
            j = int(id(ii)/w%m(1))+1
            i = id(ii) - (j-1)*w%m(1)
            L = g%Veg%SMAP(i,j)
            if (j .eq. 0) then
                j = j - 1
                i = w%m(1)
            endif
            ! Rf is runoff from ponding domain to flow domain, since flow domain has inter-cell fluxes, its mass balance cannot be examined on a column level
            ! as a result, we look at Rf
            Rf = g%OVN%Rf(i,j)
            hf = g%OVN%h(i,j)
            hI = g%Veg%hI(i,j)
            hd = g%OVN%hd(i,j) * g%OVN%dP(i,j,4)
            CS = g%Veg%CS(i,j)
            SWE= g%Veg%SWE(i,j)
            SUBM = g%Veg%SNOM(i,j)
            Qoc = g%OVN%Qoc(i,j)
            THE = sum(g%VDZ%THE(i,j,2:)*g%VDZ%DZ(i,j,2:))+sum(g%VDZ%ICE(i,j,2:))/1000D0
            !col_mass = hf+hc+hI+hd+CS+SWE+THE
            col_mass = g%VDZ%h(i,j,1)+hI+CS+SWE+THE  ! not considering anything after runoff
            prcp = g%Wea%Day%prcp(Hin,L)/24D0
            dperc= g%VDZ%Dperc(i,j)
            evap = g%Veg%Evap(i,j)
            Ex = g%OVN%Ex(i,j)
            Rf = g%OVN%Rf(i,j) ! That has not been applied to hf
            Tp = g%Veg%TP(i,j)
            test22D(i,j) = col_mass - (col_mass_s(i,j) + prcp - evap - dperc - Rf - Tp - SUBM)
            if (abs(test22D(i,j))>2D-6 .and. test22D(i,j)<100) then
                call display('budgeter col:'// num2str(test22D(i,j)) )
            endif

            col_mass_s(i,j) = col_mass
            THE0 = THE
            h0 = hc+hI
            CS0 = CS
            SWE0 = SWE
        enddo
    end select

    if(allocated(T)) deallocate(T)
    if(allocated(output)) deallocate(output)
    if(allocated(sum_temp)) deallocate(sum_temp)
    if(allocated(mask2)) deallocate(mask2)
    if(allocated(nmask2)) deallocate(nmask2)

    end subroutine budgeter


    function cvMap(dat,mask,n,mean)
    ! coefficient of variation
    implicit none
    integer n
    logical*1 mask(n)
    real*8 dat(n)
    real*8 mean,varMap,cvmap
    integer i,nn
    varMap = 0D0
    nn = 0
    do i = 1,n
        if (mask(i)) then
            varMap = varMap + (dat(i)-mean)**2D0
            nn = nn + 1
        endif
    enddo
    cvMap = sqrt(varMap/nn)/mean
    end function
    end module budgeter_mod
