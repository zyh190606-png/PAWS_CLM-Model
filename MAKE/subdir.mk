###############################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
F90_SRCS += \
../src/CLM/main/abortutils.f90 \
../src/CLM/main/ch_assert.f90 \
../src/CLM/main/surfrdMod.f90 \
../src/CLM/main/gd.f90 \
../src/CLM/main/matTable.f90 \
../src/PAWS/adhoc_Set.f90 \
../src/PAWS/allocVData.f90 \
../src/PAWS/linkVdataFromGd.f90 \
../src/CLM/biogeophys/BalanceCheckMod.f90 \
../src/CLM/biogeophys/BareGroundFluxesMod.f90 \
../src/CLM/biogeophys/Biogeophysics1Mod.f90 \
../src/CLM/biogeophys/Biogeophysics2Mod.f90 \
../src/CLM/biogeophys/BiogeophysicsLakeMod.f90 \
../src/PAWS/borrowStations.f90 \
../src/PAWS/budgeter.f90 \
../src/CLM/biogeochem/C13SummaryMod.f90 \
../src/PAWS/calendarMarch.f90 \
../src/CLM/biogeophys/CanopyFluxesMod.f90 \
../src/CLM/biogeochem/CASAPhenologyMod.f90 \
../src/PAWS/changepar.f90 \
../src/CLM/main/clm_atmlnd.f90 \
../src/CLM/main/clm_driverInitMod.f90 \
../src/CLM/main/clm_time_manager.f90 \
../src/CLM/main/clm_varcon.f90 \
../src/CLM/main/clm_varctl.f90 \
../src/CLM/main/clm_varorb.f90 \
../src/CLM/main/clm_varpar.f90 \
../src/CLM/main/clm_varsur.f90 \
../src/CLM/main/clmtype.f90 \
../src/CLM/main/clmtypeInitMod.f90 \
../src/CLM/biogeochem/CNAllocationMod.f90 \
../src/CLM/biogeochem/CNWHarvestMod.f90 \
../src/CLM/biogeochem/CNAnnualUpdateMod.f90 \
../src/CLM/biogeochem/CNBalanceCheckMod.f90 \
../src/CLM/biogeochem/CNC13FluxMod.f90 \
../src/CLM/biogeochem/CNC13StateUpdate1Mod.f90 \
../src/CLM/biogeochem/CNC13StateUpdate2Mod.f90 \
../src/CLM/biogeochem/CNC13StateUpdate3Mod.f90 \
../src/CLM/biogeochem/CNCStateUpdate1Mod.f90 \
../src/CLM/biogeochem/CNCStateUpdate2Mod.f90 \
../src/CLM/biogeochem/CNCStateUpdate3Mod.f90 \
../src/CLM/biogeochem/CNDecompMod.f90 \
../src/CLM/biogeochem/CNDVEcosystemDynIniMod.f90 \
../src/CLM/biogeochem/CNDVEstablishmentMod.f90 \
../src/CLM/biogeochem/CNDVLightMod.f90 \
../src/CLM/biogeochem/CNEcosystemDynMod.f90 \
../src/CLM/biogeochem/CNFireMod.f90 \
../src/CLM/biogeochem/CNGapMortalityMod.f90 \
../src/CLM/biogeochem/CNGRespMod.f90 \
../src/CLM/main/CNiniSpecial.f90 \
../src/CLM/main/CNiniTimeVar.f90 \
../src/CLM/biogeochem/CNMRespMod.f90 \
../src/CLM/biogeochem/CNNDynamicsMod.f90 \
../src/CLM/biogeochem/CNNStateUpdate1Mod.f90 \
../src/CLM/biogeochem/CNNStateUpdate2Mod.f90 \
../src/CLM/biogeochem/CNNStateUpdate3Mod.f90 \
../src/CLM/biogeochem/CNPhenologyMod.f90 \
../src/CLM/biogeochem/CNPrecisionControlMod.f90 \
../src/CLM/biogeochem/CNSetValueMod.f90 \
../src/CLM/biogeochem/CNSummaryMod.f90 \
../src/CLM/biogeochem/CNVegStructUpdateMod.f90 \
../src/CLM/biogeochem/CNWoodProductsMod.f90 \
../src/CLM/main/accFldsMod.f90 \
../src/CLM/main/accumulMod.f90 \
../src/PAWS/columnOutput.f90 \
../src/PAWS/comm.f90 \
../src/PAWS/compStats.f90 \
../src/PAWS/compTS.f90 \
../src/PAWS/dailyVegUpdate.f90 \
../src/CLM/main/decompMod.f90 \
../src/PAWS/distWea.f90 \
../src/CLM/main/domainMod.f90 \
../src/CLM/biogeochem/DUSTMod.f90 \
../src/CLM/main/dynlandMod.f90 \
../src/PAWS/F_oc.f90 \
../src/PAWS/F_oc_dw.f90 \
../src/CLM/main/fileutils.f90 \
../src/CLM/main/filterMod.f90 \
../src/PAWS/Flow.f90 \
../src/CLM/biogeophys/FracWetMod.f90 \
../src/CLM/biogeophys/FrictionVelocityMod.f90 \
../src/PAWS/gClearTFlux.f90 \
../src/PAWS/genVegDat.f90 \
../src/CLM/main/getdatetime.f90 \
../src/PAWS/GW_Sol.f90 \
../src/PAWS/GW_source.f90 \
../src/CLM/biogeophys/Hydrology1Mod.f90 \
../src/CLM/biogeophys/Hydrology2Mod.f90 \
../src/CLM/biogeophys/HydrologyLakeMod.f90 \
../src/CLM/main/iniTimeConst.f90 \
../src/CLM/main/initSurfAlbMod.f90 \
../src/PAWS/initValues.f90 \
../src/PAWS/makeDay.f90 \
../src/PAWS/mapOutput.f90 \
../src/PAWS/matRW.f90 \
../src/PAWS/memoryChunk.f90 \
../src/CLM/main/mkarbinitMod.f90 \
../src/CLM/main/nanMod.f90 \
../src/PAWS/openFiles.f90 \
../src/PAWS/paws.f90 \
../src/PAWS/PAWS_CLM_Mod.f90 \
../src/PAWS/paws_varctl.f90 \
../src/PAWS/PenmanMonteith.f90 \
../src/CLM/main/pft2colMod.f90 \
../src/CLM/main/pftvarcon.f90 \
../src/PAWS/prism.f90 \
../src/CLM/biogeophys/QSatMod.f90 \
../src/PAWS/r_MacCormack.f90 \
../src/PAWS/rAvgStages.f90 \
../src/PAWS/rClearTFlux.f90 \
../src/PAWS/readpar.f90 \
../src/PAWS/recAcc.f90 \
../src/PAWS/rLateral.f90 \
../src/PAWS/rMinDtSet.f90 \
../src/PAWS/run.f90 \
../src/PAWS/SaveDebug.f90 \
../src/PAWS/saveMat.f90 \
../src/PAWS/scheduleSolar.f90 \
../src/CLM/main/shr_const_mod.f90 \
../src/CLM/main/shr_file_mod.f90 \
../src/CLM/main/shr_kind_mod.f90 \
../src/CLM/main/shr_log_mod.f90 \
../src/CLM/main/shr_orb_mod.f90 \
../src/CLM/main/shr_sys_mod.f90 \
../src/PAWS/simResults3.f90 \
../src/CLM/biogeophys/SNICARMod.f90 \
../src/CLM/biogeophys/SnowHydrologyMod.f90 \
../src/PAWS/SnowLSub.f90 \
../src/PAWS/SNOWUEB_MOD.f90 \
../src/CLM/biogeophys/SoilHydrologyMod.f90 \
../src/CLM/biogeophys/SoilTemperatureMod.f90 \
../src/PAWS/Solar_Radiation_Mod.f90 \
../src/CLM/main/spmdMod.f90 \
../src/CLM/main/subgridAveMod.f90 \
../src/CLM/main/subgridMod.f90 \
../src/PAWS/surface_proc_mod.f90 \
../src/CLM/biogeophys/SurfaceAlbedoMod.f90 \
../src/CLM/biogeophys/SurfaceRadiationMod.f90 \
../src/CLM/main/system_messages.f90 \
../src/CLM/biogeophys/TridiagonalMod.f90 \
../src/PAWS/vdata.f90 \
../src/PAWS/vdb.f90 \
../src/CLM/biogeochem/VOCEmissionMod.f90 \
../src/PAWS/wDaily.f90 \
../src/PAWS/wMonthlyAverages.f90 \
../src/PAWS/writeRec.f90 \
../src/PAWS/wSyncTime.f90 \
../src/PAWS/wtrshd_col.f90 \
../src/PAWS/utilities.f90 \
../src/PAWS/wtrshd_day.f90

OBJS += \
../obj/abortutils.o \
../obj/ch_assert.o \
../obj/surfrdMod.o \
../obj/gd.o \
../obj/adhoc_Set.o \
../obj/matTable.o \
../obj/allocVData.o \
../obj/linkVdataFromGd.o \
../obj/BalanceCheckMod.o \
../obj/BareGroundFluxesMod.o \
../obj/Biogeophysics1Mod.o \
../obj/Biogeophysics2Mod.o \
../obj/BiogeophysicsLakeMod.o \
../obj/borrowStations.o \
../obj/accFldsMod.o \
../obj/accumulMod.o \
../obj/budgeter.o \
../obj/C13SummaryMod.o \
../obj/calendarMarch.o \
../obj/CanopyFluxesMod.o \
../obj/CASAPhenologyMod.o \
../obj/changepar.o \
../obj/clm_atmlnd.o \
../obj/clm_driverInitMod.o \
../obj/clm_time_manager.o \
../obj/clm_varcon.o \
../obj/clm_varctl.o \
../obj/clm_varorb.o \
../obj/clm_varpar.o \
../obj/clm_varsur.o \
../obj/clmtype.o \
../obj/clmtypeInitMod.o \
../obj/CNAllocationMod.o \
../obj/CNAnnualUpdateMod.o \
../obj/CNBalanceCheckMod.o \
../obj/CNC13FluxMod.o \
../obj/CNC13StateUpdate1Mod.o \
../obj/CNC13StateUpdate2Mod.o \
../obj/CNC13StateUpdate3Mod.o \
../obj/CNCStateUpdate1Mod.o \
../obj/CNCStateUpdate2Mod.o \
../obj/CNCStateUpdate3Mod.o \
../obj/CNDecompMod.o \
../obj/CNDVEcosystemDynIniMod.o \
../obj/CNHarvestMod.o \
../obj/CNDVEstablishmentMod.o \
../obj/CNDVLightMod.o \
../obj/CNEcosystemDynMod.o \
../obj/CNFireMod.o \
../obj/CNGapMortalityMod.o \
../obj/CNGRespMod.o \
../obj/CNiniSpecial.o \
../obj/CNiniTimeVar.o \
../obj/CNMRespMod.o \
../obj/CNNDynamicsMod.o \
../obj/CNNStateUpdate1Mod.o \
../obj/CNNStateUpdate2Mod.o \
../obj/CNNStateUpdate3Mod.o \
../obj/CNPhenologyMod.o \
../obj/CNPrecisionControlMod.o \
../obj/CNSetValueMod.o \
../obj/CNSummaryMod.o \
../obj/CNVegStructUpdateMod.o \
../obj/CNWoodProductsMod.o \
../obj/columnOutput.o \
../obj/comm.o \
../obj/compStats.o \
../obj/compTS.o \
../obj/dailyVegUpdate.o \
../obj/decompMod.o \
../obj/distWea.o \
../obj/domainMod.o \
../obj/DUSTMod.o \
../obj/dynlandMod.o \
../obj/F_oc.o \
../obj/F_oc_dw.o \
../obj/fileutils.o \
../obj/filterMod.o \
../obj/Flow.o \
../obj/FracWetMod.o \
../obj/FrictionVelocityMod.o \
../obj/gClearTFlux.o \
../obj/genVegDat.o \
../obj/getdatetime.o \
../obj/GW_Sol.o \
../obj/GW_source.o \
../obj/Hydrology1Mod.o \
../obj/Hydrology2Mod.o \
../obj/HydrologyLakeMod.o \
../obj/iniTimeConst.o \
../obj/initSurfAlbMod.o \
../obj/initValues.o \
../obj/makeDay.o \
../obj/mapOutput.o \
../obj/matRW.o \
../obj/memoryChunk.o \
../obj/mkarbinitMod.o \
../obj/nanMod.o \
../obj/openFiles.o \
../obj/paws.o \
../obj/PAWS_CLM_Mod.o \
../obj/paws_varctl.o \
../obj/PenmanMonteith.o \
../obj/pft2colMod.o \
../obj/pftvarcon.o \
../obj/prism.o \
../obj/QSatMod.o \
../obj/r_MacCormack.o \
../obj/rAvgStages.o \
../obj/rClearTFlux.o \
../obj/readpar.o \
../obj/recAcc.o \
../obj/rLateral.o \
../obj/rMinDtSet.o \
../obj/run.o \
../obj/SaveDebug.o \
../obj/saveMat.o \
../obj/scheduleSolar.o \
../obj/shr_const_mod.o \
../obj/shr_file_mod.o \
../obj/shr_kind_mod.o \
../obj/shr_log_mod.o \
../obj/shr_orb_mod.o \
../obj/shr_sys_mod.o \
../obj/simResults3.o \
../obj/SNICARMod.o \
../obj/SnowHydrologyMod.o \
../obj/SnowLSub.o \
../obj/SNOWUEB_MOD.o \
../obj/SoilHydrologyMod.o \
../obj/SoilTemperatureMod.o \
../obj/Solar_Radiation_Mod.o \
../obj/spmdMod.o \
../obj/subgridAveMod.o \
../obj/subgridMod.o \
../obj/surface_proc_mod.o \
../obj/SurfaceAlbedoMod.o \
../obj/SurfaceRadiationMod.o \
../obj/system_messages.o \
../obj/TridiagonalMod.o \
../obj/vdata.o \
../obj/vdb.o \
../obj/VOCEmissionMod.o \
../obj/wDaily.o \
../obj/wMonthlyAverages.o \
../obj/writeRec.o \
../obj/wSyncTime.o \
../obj/wtrshd_col.o \
../obj/utilities.o \
../obj/wtrshd_day.o 

VPATH = ../src/CLM/main:../src/CLM/biogeochem:../src/CLM/biogeophys:../src/PAWS

# Each subdirectory must supply rules for building sources it contributes
../obj/%.o: %.f90
	@echo 'Building file: $<'
	@echo 'Invoking: Intel(R) Intel(R) 64 Fortran Compiler'
	ifort -O2 -fpp -I"../src/include" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

../obj/abortutils.o:	../src/CLM/main/abortutils.f90 \
				../obj/clm_varctl.o \
				../obj/shr_sys_mod.o

../obj/adhoc_Set.o: 	../src/PAWS/adhoc_Set.f90 \
				../obj/vdata.o

../obj/surfrdMod.o: 	../src/CLM/main/surfrdMod.f90

../obj/allocVData.o:	../src/PAWS/allocVData.f90 \
				../obj/vdata.o \
				../obj/memoryChunk.o

../obj/accumulMod.o:	../src/CLM/main/accumulMod.f90 \
				../obj/shr_kind_mod.o \
				../obj/shr_const_mod.o \
				../obj/abortutils.o \
				../obj/clm_varctl.o \
				../obj/clm_time_manager.o \
				../obj/decompMod.o \
				../obj/spmdMod.o \
				../obj/nanMod.o \
				../obj/clm_varcon.o 

../obj/accFldsMod.o:	../src/CLM/main/accFldsMod.f90 \
				../obj/accumulMod.o \
				../obj/clm_atmlnd.o \
				../obj/pftvarcon.o \
				../obj/clmtype.o 

../obj/BalanceCheckMod.o:	../src/CLM/biogeophys/BalanceCheckMod.f90 \
				../obj/shr_kind_mod.o \
				../obj/abortutils.o \
				../obj/clm_varctl.o \
				../obj/clmtype.o \
				../obj/clm_varpar.o \
				../obj/subgridAveMod.o \
				../obj/clm_varcon.o \
				../obj/clm_atmlnd.o \
				../obj/clm_time_manager.o 

../obj/BareGroundFluxesMod.o:	../src/CLM/biogeophys/BareGroundFluxesMod.f90 \
					../obj/shr_kind_mod.o \
					../obj/clmtype.o \
					../obj/clm_atmlnd.o \
					../obj/clm_varpar.o \
					../obj/clm_varcon.o \
					../obj/shr_const_mod.o \
					../obj/FrictionVelocityMod.o \
					../obj/QSatMod.o 

../obj/Biogeophysics1Mod.o: 	../src/CLM/biogeophys/Biogeophysics1Mod.f90 \
					../obj/shr_kind_mod.o \
					../obj/clmtype.o \
					../obj/clm_atmlnd.o \
					../obj/clm_varcon.o \
					../obj/clm_varpar.o \
					../obj/QSatMod.o \
					../obj/shr_const_mod.o 

../obj/Biogeophysics2Mod.o: 	../src/CLM/biogeophys/Biogeophysics2Mod.f90 \
					../obj/shr_kind_mod.o \
					../obj/clmtype.o \
					../obj/clm_atmlnd.o \
					../obj/clm_time_manager.o \
					../obj/clm_varcon.o \
					../obj/clm_varpar.o \
					../obj/SoilTemperatureMod.o \
					../obj/subgridAveMod.o 

../obj/BiogeophysicsLakeMod.o: 	../src/CLM/biogeophys/BiogeophysicsLakeMod.f90 \
					../obj/shr_kind_mod.o \
					../obj/clmtype.o \
					../obj/clm_atmlnd.o \
					../obj/clm_time_manager.o \
					../obj/clm_varpar.o \
					../obj/clm_varcon.o \
					../obj/QSatMod.o \
					../obj/FrictionVelocityMod.o \
					../obj/TridiagonalMod.o 

../obj/borrowStaton.o: 		../src/PAWS/borrowStations.f90 \
					../obj/vdata.o 

../obj/budgeter.o: 		../src/PAWS/budgeter.f90 \
					../obj/vdata.o \
					../obj/memoryChunk.o 

../obj/C13SummaryMod.o: 		../src/CLM/biogeochem/C13SummaryMod.f90 \
					../obj/shr_kind_mod.o \
					../obj/clmtype.o \
					../obj/pft2colMod.o \
					../obj/clm_varctl.o \
					../obj/shr_sys_mod.o 

../obj/calendarMarch.o: 		../src/PAWS/calendarMarch.f90 \
					../obj/paws_varctl.o\
					../obj/wMonthlyAverages.o\
					../obj/vdata.o 

../obj/CanopyFluxesMod.o: 		../src/CLM/biogeophys/CanopyFluxesMod.f90 \
					../obj/abortutils.o \
					../obj/clm_varctl.o \
					../obj/shr_sys_mod.o \
					../obj/shr_kind_mod.o \
					../obj/clmtype.o \
					../obj/clm_atmlnd.o \
					../obj/clm_time_manager.o \
					../obj/clm_varpar.o \
					../obj/clm_varcon.o \
					../obj/QSatMod.o \
					../obj/FrictionVelocityMod.o \
					../obj/spmdMod.o \
					../obj/shr_const_mod.o \
					../obj/spmdMod.o \
					../obj/pftvarcon.o

../obj/CASAPhenologyMod.o: 	../src/CLM/biogeochem/CASAPhenologyMod.f90 \
					../obj/shr_kind_mod.o \
					../obj/clm_varpar.o \
					../obj/clmtype.o \
					../obj/decompMod.o \
					../obj/clm_varctl.o \
					../obj/clm_varcon.o \
					../obj/clm_time_manager.o

../obj/changepar.o: 		../src/PAWS/changepar.f90 \
					../obj/vdata.o 

../obj/clm_atmlnd.o: 		../src/CLM/main/clm_atmlnd.f90 \
					../obj/clm_varpar.o \
					../obj/clm_varcon.o \
					../obj/clm_varctl.o \
					../obj/decompMod.o \
					../obj/shr_kind_mod.o \
					../obj/nanMod.o \
					../obj/spmdMod.o \
					../obj/abortutils.o

../obj/clm_driverInitMod.o: 	../src/CLM/main/clm_driverInitMod.f90 \
					../obj/shr_kind_mod.o \
					../obj/clmtype.o \
					../obj/clm_varpar.o \
					../obj/subgridAveMod.o \
					../obj/clm_varcon.o \
					../obj/clm_atmlnd.o \
					../obj/QSatMod.o \
					../obj/filterMod.o 

../obj/clm_time_manager.o: 	../src/CLM/main/clm_time_manager.f90 \
                                        ../obj/paws_varctl.o

../obj/clm_varcon.o: 		../src/CLM/main/clm_varcon.f90 \
					../obj/shr_kind_mod.o \
					../obj/shr_const_mod.o \
					../obj/clm_varpar.o 

../obj/clm_varctl.o: 		../src/CLM/main/clm_varctl.f90 \
					../obj/shr_kind_mod.o \
					../obj/clm_varpar.o \
					../obj/shr_sys_mod.o 

../obj/clm_varorb.o: 		../src/CLM/main/clm_varorb.f90 \
					../obj/shr_kind_mod.o

../obj/clm_varpar.o: 		../src/CLM/main/clm_varpar.f90 \
					../obj/shr_kind_mod.o

../obj/clm_varsur.o: 		../src/CLM/main/clm_varsur.f90 \
					../obj/shr_kind_mod.o

../obj/clmtype.o: 			../src/CLM/main/clmtype.f90 \
					../obj/shr_kind_mod.o 

../obj/clmtypeInitMod.o: 		../src/CLM/main/clmtypeInitMod.f90 \
					../obj/shr_kind_mod.o \
					../obj/surfrdMod.o \
					../obj/nanMod.o \
					../obj/clmtype.o \
					../obj/clm_varpar.o \
					../obj/memoryChunk.o \
					../obj/vdata.o \
					../obj/clm_atmlnd.o \
					../obj/decompMod.o \
					../obj/clm_varcon.o

../obj/CNAllocationMod.o: 		../src/CLM/biogeochem/CNAllocationMod.f90 \
					../obj/shr_kind_mod.o \
					../obj/clmtype.o \
					../obj/clm_varctl.o \
					../obj/shr_sys_mod.o \
					../obj/clm_time_manager.o \
					../obj/surfrdMod.o \
					../obj/pft2colMod.o 

../obj/CNAnnualUpdateMod.o: 	../src/CLM/biogeochem/CNAnnualUpdateMod.f90 \
					../obj/shr_kind_mod.o \
					../obj/clmtype.o \
					../obj/clm_time_manager.o \
					../obj/clm_varcon.o \
					../obj/pft2colMod.o 

../obj/CNBalanceCheckMod.o: 	../src/CLM/biogeochem/CNBalanceCheckMod.f90 \
					../obj/abortutils.o \
					../obj/shr_kind_mod.o \
					../obj/clm_varctl.o \
					../obj/clmtype.o \
					../obj/clm_time_manager.o 

../obj/CNC13FluxMod.o: 		../src/CLM/biogeochem/CNC13FluxMod.f90 \
					../obj/shr_kind_mod.o \
					../obj/clmtype.o \
					../obj/pft2colMod.o \
					../obj/clm_varpar.o 

../obj/CNC13StateUpdate1Mod.o: 	../src/CLM/biogeochem/CNC13StateUpdate1Mod.f90 \
					../obj/shr_kind_mod.o \
					../obj/clmtype.o \
					../obj/clm_time_manager.o 

../obj/CNC13StateUpdate2Mod.o: 	../src/CLM/biogeochem/CNC13StateUpdate2Mod.f90 \
					../obj/shr_kind_mod.o \
					../obj/clmtype.o \
					../obj/clm_time_manager.o

../obj/CNC13StateUpdate3Mod.o: 	../src/CLM/biogeochem/CNC13StateUpdate3Mod.f90 \
					../obj/shr_kind_mod.o \
					../obj/clmtype.o \
					../obj/clm_time_manager.o

../obj/CNCStateUpdate1Mod.o: 	../src/CLM/biogeochem/CNCStateUpdate1Mod.f90 \
					../obj/shr_kind_mod.o \
					../obj/pftvarcon.o \
					../obj/clmtype.o \
					../obj/clm_time_manager.o

../obj/CNCStateUpdate2Mod.o: 	../src/CLM/biogeochem/CNCStateUpdate2Mod.f90 \
					../obj/shr_kind_mod.o \
					../obj/clmtype.o \
					../obj/clm_time_manager.o

../obj/CNCStateUpdate3Mod.o: 	../src/CLM/biogeochem/CNCStateUpdate3Mod.f90 \
					../obj/shr_kind_mod.o \
					../obj/clmtype.o \
					../obj/clm_time_manager.o

../obj/CNDecompMod.o: 		../src/CLM/biogeochem/CNDecompMod.f90 \
					../obj/shr_kind_mod.o \
					../obj/shr_const_mod.o \
					../obj/clmtype.o \
					../obj/CNAllocationMod.o \
					../obj/clm_time_manager.o \
					../obj/pft2colMod.o

../obj/CNDVEcosystemDynIniMod.o: 	../src/CLM/biogeochem/CNDVEcosystemDynIniMod.f90 \
					../obj/shr_kind_mod.o \
					../obj/clmtype.o \
					../obj/decompMod.o \
					../obj/shr_const_mod.o 

../obj/CNDVEstablishmentMod.o: 	../src/CLM/biogeochem/CNDVEstablishmentMod.f90 \
					../obj/shr_kind_mod.o \
					../obj/abortutils.o \
					../obj/clmtype.o \
					../obj/clm_varpar.o \
					../obj/clm_varcon.o \
					../obj/clm_varctl.o \
					../obj/pftvarcon.o \
					../obj/shr_const_mod.o 

../obj/CNDVLightMod.o: 		../src/CLM/biogeochem/CNDVLightMod.f90 \
					../obj/shr_kind_mod.o \
					../obj/shr_const_mod.o \
					../obj/clmtype.o 

../obj/CNEcosystemDynMod.o: 		../src/CLM/biogeochem/CNEcosystemDynMod.f90 \
					../obj/shr_kind_mod.o \
					../obj/clm_varctl.o \
					../obj/clmtype.o \
					../obj/CNSetValueMod.o \
					../obj/CNNDynamicsMod.o \
					../obj/CNMRespMod.o \
					../obj/CNDecompMod.o \
					../obj/CNPhenologyMod.o \
					../obj/CNGRespMod.o \
					../obj/CNCStateUpdate1Mod.o \
					../obj/CNNStateUpdate1Mod.o \
					../obj/CNGapMortalityMod.o \
					../obj/CNCStateUpdate2Mod.o \
					../obj/CNNStateUpdate2Mod.o \
					../obj/CNCStateUpdate3Mod.o \
					../obj/CNNStateUpdate3Mod.o \
					../obj/CNBalanceCheckMod.o \
					../obj/CNPrecisionControlMod.o \
					../obj/CNVegStructUpdateMod.o \
					../obj/CNAnnualUpdateMod.o \
					../obj/CNSummaryMod.o \
					../obj/CNC13StateUpdate1Mod.o \
					../obj/CNC13StateUpdate2Mod.o \
					../obj/CNC13StateUpdate3Mod.o \
					../obj/CNC13FluxMod.o \
					../obj/C13SummaryMod.o \
					../obj/CNHarvestMod.o \
					../obj/CNWoodProductsMod.o \
	                ../obj/CNFireMod.o

../obj/CNHarvestMod.o: 		../src/CLM/biogeochem/CNHarvestMod.f90 \
					../obj/spmdMod.o \
					../obj/clmtype.o \
					../obj/clm_varpar.o \
					../obj/decompMod.o \
					../obj/clm_varsur.o \
					../obj/shr_sys_mod.o \
					../obj/abortutils.o \
					../obj/pftvarcon.o \
					../obj/clm_varcon.o \
					../obj/clm_time_manager.o \
					../obj/shr_kind_mod.o 

../obj/CNFireMod.o: 		../src/CLM/biogeochem/CNFireMod.f90 \
					../obj/shr_kind_mod.o \
					../obj/shr_const_mod.o \
					../obj/pft2colMod.o \
					../obj/clm_varctl.o \
					../obj/clmtype.o \
					../obj/clm_time_manager.o \
					../obj/clm_varpar.o 

../obj/CNGapMortalityMod.o: 	../src/CLM/biogeochem/CNGapMortalityMod.f90 \
					../obj/shr_kind_mod.o \
					../obj/clmtype.o \
					../obj/clm_varcon.o \
					../obj/clm_time_manager.o \
					../obj/clm_varpar.o 

../obj/CNGRespMod.o: 		../src/CLM/biogeochem/CNGRespMod.f90 \
					../obj/pftvarcon.o \
					../obj/shr_kind_mod.o \
					../obj/clmtype.o 

../obj/CNiniSpecial.o: 		../src/CLM/main/CNiniSpecial.f90 \
					../obj/shr_kind_mod.o \
					../obj/pftvarcon.o \
					../obj/decompMod.o \
					../obj/clm_varcon.o \
					../obj/clm_varctl.o \
					../obj/clmtype.o \
					../obj/CNSetValueMod.o

../obj/CNiniTimeVar.o: 		../src/CLM/main/CNiniTimeVar.f90 \
					../obj/clmtype.o \
					../obj/clm_atmlnd.o \
					../obj/shr_kind_mod.o \
					../obj/clm_varcon.o \
					../obj/pftvarcon.o \
					../obj/decompMod.o   

../obj/CNMRespMod.o: 		../src/CLM/biogeochem/CNMRespMod.f90 \
					../obj/shr_kind_mod.o \
					../obj/clm_varpar.o \
					../obj/pftvarcon.o \
					../obj/shr_const_mod.o \
					../obj/clmtype.o 

../obj/CNNDynamicsMod.o: 		../src/CLM/biogeochem/CNNDynamicsMod.f90 \
					../obj/shr_kind_mod.o \
					../obj/clmtype.o \
					../obj/clm_atmlnd.o \
					../obj/clm_varctl.o \
					../obj/shr_sys_mod.o \
					../obj/clm_varpar.o \
					../obj/clm_time_manager.o   

../obj/CNNStateUpdate1Mod.o: 	../src/CLM/biogeochem/CNNStateUpdate1Mod.f90 \
					../obj/shr_kind_mod.o \
					../obj/pftvarcon.o \
					../obj/clmtype.o \
					../obj/clm_time_manager.o   

../obj/CNNStateUpdate2Mod.o: 	../src/CLM/biogeochem/CNNStateUpdate2Mod.f90 \
					../obj/shr_kind_mod.o \
					../obj/pftvarcon.o \
					../obj/clmtype.o \
					../obj/clm_time_manager.o   

../obj/CNNStateUpdate3Mod.o: 	../src/CLM/biogeochem/CNNStateUpdate3Mod.f90 \
					../obj/shr_kind_mod.o \
					../obj/pftvarcon.o \
					../obj/clmtype.o \
					../obj/clm_time_manager.o   

../obj/CNPhenologyMod.o: 		../src/CLM/biogeochem/CNPhenologyMod.f90 \
					../obj/clmtype.o \
					../obj/pftvarcon.o \
					../obj/clm_varcon.o \
					../obj/shr_kind_mod.o \
					../obj/clm_time_manager.o \
					../obj/shr_const_mod.o \
					../obj/clm_varpar.o 

../obj/CNPrecisionControlMod.o: 	../src/CLM/biogeochem/CNPrecisionControlMod.f90 \
					../obj/shr_kind_mod.o \
					../obj/clmtype.o \
					../obj/pftvarcon.o \
					../obj/abortutils.o \
					../obj/clm_varctl.o

../obj/CNSetValueMod.o: 		../src/CLM/biogeochem/CNSetValueMod.f90 \
					../obj/shr_kind_mod.o \
					../obj/clm_varpar.o \
					../obj/clm_varctl.o \
					../obj/clmtype.o \
					../obj/decompMod.o 

../obj/CNSummaryMod.o: 		../src/CLM/biogeochem/CNSummaryMod.f90 \
					../obj/shr_kind_mod.o \
					../obj/clmtype.o \
					../obj/pft2colMod.o \
					../obj/clm_varctl.o \
					../obj/shr_sys_mod.o 

../obj/CNVegStructUpdateMod.o: 	../src/CLM/biogeochem/CNVegStructUpdateMod.f90 \
					../obj/shr_kind_mod.o \
					../obj/clmtype.o \
					../obj/clm_atmlnd.o \
					../obj/pftvarcon.o \
					../obj/shr_const_mod.o \
					../obj/clm_time_manager.o

../obj/CNWoodProductsMod.o: 	../src/CLM/biogeochem/CNWoodProductsMod.f90 \
					../obj/decompMod.o \
					../obj/shr_kind_mod.o \
					../obj/clm_varcon.o \
					../obj/spmdMod.o \
					../obj/clmtype.o \
					../obj/clm_time_manager.o

../obj/columnOutput.o: 		../src/PAWS/columnOutput.f90 \
					../obj/vdata.o 

../obj/comm.o: 			../src/PAWS/comm.f90 \

../obj/compStats.o: 		../src/PAWS/compStats.f90 \
					../obj/vdata.o 

../obj/compTS.o: 			../src/PAWS/compTS.f90 \
					../obj/vdata.o 

../obj/dailyVegUpdate.o: 		../src/PAWS/dailyVegUpdate.f90 \
					../obj/vdata.o 

../obj/decompMod.o: 		../src/CLM/main/decompMod.f90 \
					../obj/memoryChunk.o \
					../obj/nanMod.o 

../obj/distWea.o: 			../src/PAWS/distWea.f90 \
					../obj/vdata.o \
					../obj/utilities.o

../obj/domainMod.o: 		../src/CLM/main/domainMod.f90 \
					../obj/shr_kind_mod.o \
					../obj/nanMod.o \
					../obj/spmdMod.o \
					../obj/abortutils.o \
					../obj/clm_varctl.o 

../obj/DUSTMod.o: 			../src/CLM/biogeochem/DUSTMod.f90 \
					../obj/shr_kind_mod.o \
					../obj/clmtype.o \
					../obj/clm_varpar.o \
					../obj/clm_varcon.o \
					../obj/clm_varctl.o \
					../obj/abortutils.o \
					../obj/subgridAveMod.o \
					../obj/clm_atmlnd.o \
					../obj/shr_const_mod.o \
					../obj/decompMod.o

../obj/dynlandMod.o: 		../src/CLM/main/dynlandMod.f90 \
					../obj/spmdMod.o \
					../obj/clmtype.o \
					../obj/decompMod.o \
					../obj/clm_varctl.o \
					../obj/shr_kind_mod.o \
					../obj/abortutils.o \
					../obj/clm_varcon.o \
					../obj/clm_varpar.o     

../obj/F_oc.o: 			../src/PAWS/F_oc.f90

../obj/F_oc_dw.o: 			../src/PAWS/F_oc_dw.f90 \
					../obj/vdata.o 

../obj/fileutils.o: 		../src/CLM/main/fileutils.f90 \
					../obj/abortutils.o \
					../obj/clm_varctl.o \
					../obj/shr_file_mod.o \
					../obj/clm_varctl.o  

../obj/filterMod.o: 		../src/CLM/main/filterMod.f90 \
					../obj/shr_kind_mod.o \
					../obj/abortutils.o \
					../obj/clm_varctl.o \
					../obj/pftvarcon.o \
					../obj/clmtype.o \
					../obj/decompMod.o \
					../obj/clm_varcon.o \
					../obj/clm_varpar.o   

../obj/Flow.o: 			../src/PAWS/Flow.f90 \
					../obj/clmtype.o \
					../obj/clm_varcon.o \
					../obj/vdata.o 

../obj/linkVdataFromGd.o: 	../src/PAWS/linkVdataFromGd.f90 \
					../obj/vdata.o \
					../obj/gd.o 

../obj/FracWetMod.o: 		../src/CLM/biogeophys/FracWetMod.f90 \
					../obj/shr_kind_mod.o \
					../obj/clmtype.o

../obj/FrictionVelocityMod.o: 	../src/CLM/biogeophys/FrictionVelocityMod.f90 \
					../obj/shr_kind_mod.o \
					../obj/clmtype.o \
					../obj/clm_atmlnd.o \
					../obj/clm_varcon.o \
					../obj/clm_varctl.o \
					../obj/shr_const_mod.o 

../obj/gClearTFlux.o: 		../src/PAWS/gClearTFlux.f90 \
					../obj/vdata.o 

../obj/genVegDat.o: 		../src/PAWS/genVegDat.f90 \
					../obj/vdata.o 

../obj/getdatetime.o: 		../src/CLM/main/getdatetime.f90

../obj/Gw_Sol.o: 			../src/PAWS/GW_Sol.f90 \
					../obj/vdata.o 

../obj/GW_source.o: 		../src/PAWS/GW_source.f90 \
					../obj/vdata.o  \
					../obj/memoryChunk.o 

../obj/Hydrology1Mod.o: 		../src/CLM/biogeophys/Hydrology1Mod.f90 \
					../obj/shr_kind_mod.o \
					../obj/clmtype.o \
					../obj/clm_atmlnd.o \
					../obj/clm_varcon.o \
					../obj/FracWetMod.o \
					../obj/clm_time_manager.o \
					../obj/subgridAveMod.o \
					../obj/SNICARMod.o 

../obj/Hydrology2Mod.o: 		../src/CLM/biogeophys/Hydrology2Mod.f90 \
					../obj/shr_kind_mod.o \
					../obj/clmtype.o \
					../obj/clm_atmlnd.o \
					../obj/clm_varcon.o \
					../obj/clm_varctl.o \
					../obj/clm_varpar.o \
					../obj/SnowHydrologyMod.o \
					../obj/SoilHydrologyMod.o \
					../obj/clm_time_manager.o 

../obj/HydrologyLakeMod.o: 	../src/CLM/biogeophys/HydrologyLakeMod.f90 \
					../obj/shr_kind_mod.o \
					../obj/clmtype.o \
					../obj/clm_atmlnd.o \
					../obj/clm_time_manager.o \
					../obj/clm_varcon.o 

../obj/iniTimeConst.o: 		../src/CLM/main/iniTimeConst.f90 \
					../obj/shr_kind_mod.o \
					../obj/nanMod.o \
					../obj/clmtype.o \
					../obj/decompMod.o \
					../obj/clm_atmlnd.o \
					../obj/clm_varpar.o \
					../obj/clm_varcon.o \
					../obj/clm_varctl.o \
					../obj/pftvarcon.o \
					../obj/clm_time_manager.o \
					../obj/abortutils.o \
					../obj/fileutils.o \
					../obj/spmdMod.o \
					../obj/SNICARMod.o \
					../obj/vdata.o \
					../obj/memoryChunk.o \
					../obj/filterMod.o 

../obj/initSurfAlbMod.o: 		../src/CLM/main/initSurfAlbMod.f90 \
					../obj/shr_kind_mod.o \
					../obj/abortutils.o \
					../obj/clm_varctl.o \
					../obj/shr_orb_mod.o \
					../obj/shr_const_mod.o \
					../obj/clmtype.o \
					../obj/spmdMod.o \
					../obj/decompMod.o \
					../obj/filterMod.o \
					../obj/clm_varpar.o \
					../obj/clm_varcon.o \
					../obj/clm_time_manager.o \
					../obj/FracWetMod.o \
					../obj/SurfaceAlbedoMod.o \
					../obj/CNEcosystemDynMod.o \
					../obj/CNVegStructUpdateMod.o \
					../obj/CNSetValueMod.o 

../obj/initValue.o: 		../src/PAWS/initValues.f90 \
					../obj/vdata.o 

../obj/makeDay.o: 			../src/PAWS/makeDay.f90 \
					../obj/vdata.o 

../obj/mapOutput.o: 		../src/PAWS/mapOutput.f90 \
					../obj/vdata.o \
					../obj/filterMod.o \
					../obj/matRW.o

../obj/matRW.o: 			../src/PAWS/matRW.f90 \
					../obj/vdata.o \
					../obj/ch_assert.o \
					../obj/memoryChunk.o

../obj/memoryChunk.o: 		../src/PAWS/memoryChunk.f90 \
					../obj/gd.o \
					../obj/ch_assert.o \
					../obj/matTable.o 
 

../obj/mkarbinitMod.o: 		../src/CLM/main/mkarbinitMod.f90 \
					../obj/shr_kind_mod.o \
					../obj/shr_const_mod.o \
					../obj/clmtype.o \
					../obj/clm_varpar.o \
					../obj/clm_varcon.o \
					../obj/clm_varctl.o \
					../obj/spmdMod.o \
					../obj/decompMod.o \
					../obj/shr_sys_mod.o \
					../obj/SNICARMod.o

../obj/nanMod.o: 			../src/CLM/main/nanMod.f90 \
					../obj/shr_kind_mod.o

../obj/openFiles.o: 		../src/PAWS/openFiles.f90 \
					../obj/budgeter.o \
					../obj/vdata.o 

../obj/paws.o: 			../src/PAWS/paws.f90 \
					../obj/vdata.o \
					../obj/matRW.o

../obj/PAWS_CLM_Mod.o: 		../src/PAWS/PAWS_CLM_Mod.f90 \
					../obj/ch_assert.o \
					../obj/matTable.o \
					../obj/clmtype.o \
					../obj/clmtypeInitMod.o \
					../obj/clm_varpar.o \
					../obj/vdata.o \
					../obj/memoryChunk.o \
					../obj/filterMod.o \
					../obj/abortutils.o \
					../obj/clm_varctl.o \
					../obj/mkarbinitMod.o \
					../obj/pftvarcon.o \
					../obj/initSurfAlbMod.o \
					../obj/clm_time_manager.o \
					../obj/shr_orb_mod.o \
					../obj/clm_varorb.o \
					../obj/shr_kind_mod.o \
					../obj/spmdMod.o \
					../obj/decompMod.o \
					../obj/clm_varcon.o \
					../obj/clm_driverInitMod.o \
					../obj/BalanceCheckMod.o \
					../obj/SurfaceRadiationMod.o \
					../obj/Hydrology1Mod.o \
					../obj/Hydrology2Mod.o \
					../obj/HydrologyLakeMod.o \
					../obj/Biogeophysics1Mod.o \
					../obj/BareGroundFluxesMod.o \
					../obj/CanopyFluxesMod.o \
					../obj/Biogeophysics2Mod.o \
					../obj/BiogeophysicsLakeMod.o \
					../obj/SurfaceAlbedoMod.o \
					../obj/pft2colMod.o \
					../obj/CNSetValueMod.o \
					../obj/CNEcosystemDynMod.o \
					../obj/CNAnnualUpdateMod.o \
					../obj/CNBalanceCheckMod.o \
					../obj/DUSTMod.o \
					../obj/SNICARMod.o \
					../obj/shr_const_mod.o \
					../obj/domainMod.o \
					../obj/clm_atmlnd.o \
					../obj/vdb.o \
					../obj/QSatMod.o \
					../obj/SNOWUEB_MOD.o \
					../obj/paws_varctl.o 

../obj/paws_varctl.o: 		../src/PAWS/paws_varctl.f90 

../obj/PenmanMonteith.o: 		../src/PAWS/PenmanMonteith.f90 

../obj/pft2colMod.o: 		../src/CLM/main/pft2colMod.f90 \
					../obj/shr_kind_mod.o \
					../obj/subgridMod.o \
					../obj/clmtype.o 

../obj/pftvarcon.o: 		../src/CLM/main/pftvarcon.f90 \
					../obj/shr_kind_mod.o \
					../obj/abortutils.o \
					../obj/clm_varpar.o \
					../obj/clm_varctl.o \
					../obj/fileutils.o \
					../obj/spmdMod.o \
					../obj/nanMod.o

../obj/prism.o: 			../src/PAWS/prism.f90 \
					../obj/surface_proc_mod.o \
					../obj/vdata.o \
					../obj/Flow.o \
					../obj/PAWS_CLM_Mod.o \
					../obj/comm.o \
					../obj/memoryChunk.o

../obj/QSatMod.o: 			../src/CLM/biogeophys/QSatMod.f90 \
					../obj/shr_kind_mod.o \
					../obj/shr_const_mod.o 

../obj/r_MacCormack.o: 		../src/PAWS/r_MacCormack.f90 \
					../obj/vdata.o \
					../obj/Flow.o

../obj/rAvgStages.o: 		../src/PAWS/rAvgStages.f90 \
					../obj/vdata.o 

../obj/rClearTFlux.o: 		../src/PAWS/rClearTFlux.f90 \
					../obj/vdata.o 

../obj/readpar.o: 			../src/PAWS/readpar.f90 \
					../obj/vdata.o 

../obj/recAcc.o: 			../src/PAWS/recAcc.f90 \
					../obj/vdata.o 

../obj/rLateral.o: 		../src/PAWS/rLateral.f90 \
					../obj/utilities.o \
					../obj/vdata.o 

../obj/rMinDtSet.o: 		../src/PAWS/rMinDtSet.f90 \
					../obj/vdata.o 

../obj/run.o: 			../src/PAWS/run.f90 \
					../obj/vdata.o 

../obj/SaveDebug.o: 		../src/PAWS/SaveDebug.f90 \
					../obj/matRW.o \
					../obj/vdata.o 

../obj/saveMat.o: 			../src/PAWS/saveMat.f90 \
					../obj/matRW.o \
					../obj/memoryChunk.o 

../obj/scheduleSolar.o: 		../src/PAWS/scheduleSolar.f90 \
					../obj/vdata.o \
					../obj/Solar_Radiation_Mod.o \
					../obj/shr_orb_mod.o \
					../obj/clm_varorb.o

../obj/shr_const_mod.o: 		../src/CLM/main/shr_const_mod.f90 \
					../obj/shr_kind_mod.o

../obj/shr_file_mod.o: 		../src/CLM/main/shr_file_mod.f90 \
					../obj/shr_kind_mod.o \
					../obj/shr_sys_mod.o \
					../obj/shr_log_mod.o

../obj/shr_kind_mod.o: 		../src/CLM/main/shr_kind_mod.f90

../obj/shr_log_mod.o: 		../src/CLM/main/shr_log_mod.f90 \
					../obj/shr_kind_mod.o

../obj/shr_orb_mod.o: 		../src/CLM/main/shr_orb_mod.f90 \
					../obj/shr_kind_mod.o \
					../obj/shr_sys_mod.o \
					../obj/shr_const_mod.o

../obj/shr_sys_mod.o: 		../src/CLM/main/shr_sys_mod.f90 \
					../obj/shr_kind_mod.o

../obj/simResults3.o: 		../src/PAWS/simResults3.f90 \
					../obj/vdata.o \
					../obj/matRW.o

../obj/SNICARMod.o: 		../src/CLM/biogeophys/SNICARMod.f90 \
					../obj/shr_kind_mod.o \
					../obj/shr_sys_mod.o \
					../obj/clm_varctl.o \
					../obj/shr_const_mod.o \
					../obj/clmtype.o \
					../obj/clm_varpar.o \
					../obj/clm_time_manager.o \
					../obj/abortutils.o \
					../obj/clm_varcon.o \
					../obj/fileutils.o \
					../obj/spmdMod.o \
					../obj/matTable.o

../obj/SnowHydrologyMod.o: 		../src/CLM/biogeophys/SnowHydrologyMod.f90 \
					../obj/shr_kind_mod.o \
					../obj/clm_varpar.o \
					../obj/clmtype.o \
					../obj/clm_varcon.o \
					../obj/clm_time_manager.o \
					../obj/clm_atmlnd.o \
					../obj/SNICARMod.o

../obj/SnowLSub.o: 		../src/PAWS/SnowLSub.f90 \
					../obj/comm.o \
					../obj/SNOWUEB_MOD.o

../obj/SNOWUEB_MOD.o: 		../src/PAWS/SNOWUEB_MOD.f90 \
					../obj/comm.o

../obj/SoilHydrologyMod.o: 	../src/CLM/biogeophys/SoilHydrologyMod.f90 \
					../obj/shr_kind_mod.o \
					../obj/clmtype.o \
					../obj/clm_varcon.o \
					../obj/clm_varpar.o \
					../obj/clm_time_manager.o \
					../obj/clm_varctl.o \
					../obj/shr_const_mod.o \
					../obj/TridiagonalMod.o

../obj/SoilTemperatureMod.o: 	../src/CLM/biogeophys/SoilTemperatureMod.f90 \
					../obj/shr_kind_mod.o \
					../obj/clmtype.o \
					../obj/clm_atmlnd.o \
					../obj/clm_time_manager.o \
					../obj/clm_varctl.o \
					../obj/clm_varcon.o \
					../obj/clm_varpar.o \
					../obj/TridiagonalMod.o

../obj/Solar_Radiation_Mod.o: 	../src/PAWS/Solar_Radiation_Mod.f90 \
					../obj/vdata.o \
					../obj/shr_orb_mod.o

../obj/spmdMod.o: 			../src/CLM/main/spmdMod.f90

../obj/subgridAveMod.o: 		../src/CLM/main/subgridAveMod.f90 \
					../obj/shr_kind_mod.o \
					../obj/clmtype.o \
					../obj/clm_varcon.o \
					../obj/clm_varctl.o \
					../obj/abortutils.o \
					../obj/clm_varpar.o

../obj/subgridMod.o: 		../src/CLM/main/subgridMod.f90 \
					../obj/shr_kind_mod.o \
					../obj/spmdMod.o \
					../obj/nanMod.o \
					../obj/abortutils.o \
					../obj/clm_varpar.o \
					../obj/clm_varctl.o \
					../obj/clm_varsur.o

../obj/surface_proc_mod.o: 	../src/PAWS/surface_proc_mod.f90 \
					../obj/vdata.o

../obj/SurfaceAlbedoMod.o: 	../src/CLM/biogeophys/SurfaceAlbedoMod.f90 \
					../obj/clm_varcon.o \
					../obj/clm_varpar.o \
					../obj/shr_kind_mod.o \
					../obj/SNICARMod.o \
					../obj/clmtype.o \
					../obj/shr_orb_mod.o \
					../obj/clm_time_manager.o

../obj/SurfaceRadiationMod.o: 	../src/CLM/biogeophys/SurfaceRadiationMod.f90 \
					../obj/clmtype.o \
					../obj/clm_atmlnd.o \
					../obj/clm_varpar.o \
					../obj/clm_varcon.o \
					../obj/clm_time_manager.o \
					../obj/SNICARMod.o \
					../obj/abortutils.o \
					../obj/shr_kind_mod.o \
					../obj/clm_varctl.o

../obj/system_messages.o: 		../src/CLM/main/system_messages.f90 \
					../obj/abortutils.o \
					../obj/clm_varctl.o

../obj/TridiagonalMod.o: 		../src/CLM/biogeophys/TridiagonalMod.f90 \
					../obj/shr_kind_mod.o

../obj/vdata.o: 			../src/PAWS/vdata.f90 \
					../obj/ch_assert.o

../obj/ch_assert.o: 			../src/PAWS/ch_assert.f90

../obj/gd.o: 			        ../src/PAWS/gd.f90 \
					../obj/ch_assert.o


../obj/matTable.o:			../src/PAWS/matTable.f90 \
					../obj/gd.o 


../obj/vdb.o: 			../src/PAWS/vdb.f90

../obj/VOCEmissionMod.o: 		../src/CLM/biogeochem/VOCEmissionMod.f90 \
					../obj/clm_varctl.o \
					../obj/abortutils.o \
					../obj/shr_kind_mod.o \
					../obj/clm_atmlnd.o \
					../obj/clmtype.o \
					../obj/clm_varpar.o \
					../obj/shr_const_mod.o \
					../obj/clm_varcon.o \
					../obj/pftvarcon.o

../obj/wDaily.o: 			../src/PAWS/wDaily.f90 \
					../obj/vdata.o

../obj/wMonthlyAverages.o: 	../src/PAWS/wMonthlyAverages.f90 \
					../obj/vdata.o

../obj/writeRec.o: 		../src/PAWS/writeRec.f90 \
					../obj/vdata.o

../obj/wSyncTime.o: 		../src/PAWS/wSyncTime.f90 \
					../obj/vdata.o

../obj/wtrshd_col.o: 		../src/PAWS/wtrshd_col.f90 \
					../obj/vdata.o \
					../obj/Flow.o \
					../obj/surface_proc_mod.o \
					../obj/budgeter.o

../obj/wtrshd_day.o: 		../src/PAWS/wtrshd_day.f90 \
					../obj/vdata.o \
					../obj/Flow.o \
					../obj/surface_proc_mod.o \
					../obj/budgeter.o
