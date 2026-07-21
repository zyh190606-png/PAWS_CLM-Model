################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include ../makefile.init

RM := rm -rf

# All of the sources participating in the build are defined here
-include sources.mk
-include subdir.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
endif

-include ../makefile.defs

# Add inputs and outputs from these tool invocations to the build variables 
suffix:
	@echo 'change all .F90 to .f90'
	./suffix.sh

# All Target
all: suffix PAWS_CLM

# Tool invocations
PAWS_CLM: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: Intel(R) Fortran Linker'
	ifort -O2 $(LIBS) -o "../bin/PAWS_CLM" $(OBJS) $(USER_OBJS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(OBJS)$(C_DEPS)$(EXECUTABLES) ../bin/PAWS_CLM ./*.mod ./*.i90
	-@echo ' '

.PHONY: all clean dependents suffix
.SECONDARY:

-include ../makefile.targets
