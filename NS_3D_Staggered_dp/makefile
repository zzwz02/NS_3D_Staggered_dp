#===============================================================================
# Copyright 2005-2018 Intel Corporation.
#
# This software and the related documents are Intel copyrighted  materials,  and
# your use of  them is  governed by the  express license  under which  they were
# provided to you (License).  Unless the License provides otherwise, you may not
# use, modify, copy, publish, distribute,  disclose or transmit this software or
# the related documents without Intel's prior written permission.
#
# This software and the related documents  are provided as  is,  with no express
# or implied  warranties,  other  than those  that are  expressly stated  in the
# License.
#===============================================================================

##  Content:
##      Intel(R) Math Kernel Library LAPACK95 examples creation and run
##******************************************************************************

#help:
#	@echo "Usage: make {libia32|soia32|libintel64|sointel64} [function=name]"
#	@echo "[compiler=compiler_name] [interface=interface_name] [threading=threading_name]"
#	@echo "[parallel=parallel_name] [omp=omp_name] [gnu_path=gnu_lib_path]"
#	@echo "[F95ROOT=<f95_root_path>] [MKLROOT=<MKL_root_path>]"
#	@echo "name           - function name. Please see lapack95.lst file."
#	@echo "compiler_name  - can be gnu or intel. Default value is intel."
#	@echo "                 Intel (R) Fortran Compiler as default."
#	@echo "                 If compiler=gnu then GNU gfortran compiler will be used."
#	@echo "                 Note: GNU g77 compiler is not supported."
#	@echo "interface_name - can be lp64 or ilp64 for intel64. Default value is lp64."
#	@echo "threading_name - can be parallel or sequential. Default value is parallel."
#	@echo "parallel_name  - can be intel or gnu. Default value is intel."
#	@echo "omp_name       - can be iomp5 if parallel=intel or"
#	@echo "                 can be iomp5 or gomp if parallel=gnu."
#	@echo "                 Default value is iomp5."
#	@echo "gnu_lib_path   - If you are using gnu threading layer, specify path to GNU libraries,"
#	@echo "                 such as libgomp and libgfortran, with gnu_path variable."
#	@echo "                 Default value is /usr/lib."
#	@echo "f95_root_path  - Path to lapack95 library and .mod files that are built on user side."
#	@echo "                 Mandatory if compiler is not intel. To build library and .mod files,"
#	@echo "                 go to the directory <MKLROOT>/interfaces/lapack95 and run:"
#	@echo "                 make <target> [options] FC=<f95compiler> INSTALL_DIR=<f95_root_path>"
#	@echo "                 See description of the FC and INSTALL_DIR in the interface makefile."
#	@echo "MKL_root_path  - Locate Intel(R) MKL libraries relative to <MKL_root_path>."
#	@echo "                 Default value: ../.., unless defined in environment."

##------------------------------------------------------------------------------
## examples of using:
##
## make libia32 function=geevx
##              - build  by  Intel(R) Fortran Compiler  (as default)
##                and run DGEBRD example for 32-bit applications,
##                static linking
##
## make libintel64 compiler=gnu F95ROOT=f95_root_path
##              - build by gfortran GNU  Fortran  compiler  and  run  all
##                examples of Intel(R) MKL for Intel(R) 64 processor family applications,
##                static linking
##
## make  sointel64
##              - build by Intel(R) Fortran Compiler (as default) and run all
##                examples of Intel(R) MKL for Intel(R) 64 processor family applications,
##                dynamic linking
##------------------------------------------------------------------------------

#include lapack95.lst

##  libia32 lib32:
##	$(MAKE) $(RES) EXT=a _IA=ia32 PLATF=ia32 RES_EXT=lib
##  soia32 so32:
##	$(MAKE) $(RES) EXT=so _IA=ia32 PLATF=ia32 RES_EXT=so
##  libintel64 libem64t:
##	$(MAKE) $(RES) EXT=a _IA=intel64 PLATF=intel64 RES_EXT=lib
##  sointel64 soem64t:
##	$(MAKE) $(RES) EXT=so _IA=intel64 PLATF=intel64 RES_EXT=so

EXT=a
_IA=intel64
PLATF=intel64
RES_EXT=lib

compiler=intel

interface=lp64
##threading=sequential

ifneq ($(threading),sequential)
   override threading=omp
endif

parallel=intel
omp=iomp5

## RES = $(addsuffix .res ,$(function))
ifndef MKLROOT
MKLROOT = ../..
endif
MKL_PATH = "$(MKLROOT)/lib/$(_IA)"
CMPLR_PATH = "$(MKLROOT)/../compiler/lib/$(_IA)"

FOPTS = -fpp -O2 -debug full -g -traceback -check all -fp-stack-check -qopt-matmul -init=snan -init=arrays -fpconstant -qopenmp -mkl=parallel ### compiler flag
LOPTS =          ### linking flag

FC=ifort
IFACE_COMP_PART=intel

IFACE_LIB = mkl_$(IFACE_COMP_PART)_lp64

## GFORTRAN_LIB =

ifeq ($(threading),sequential)
   THREAD_LIB = mkl_sequential
   LL_THREAD_LIB =
else
   ifeq ($(threading),omp)
      THREAD_LIB = mkl_$(parallel)_thread
      LL_THREAD_LIB = -L$(CMPLR_PATH)
      LL_THREAD_LIB += -l$(omp)
   endif
endif

CORE_LIB = mkl_core

ifeq ($(EXT),a)
    MKL_LIBS = -Wl,--start-group $(MKL_PATH)/$(RES_EXT)$(IFACE_LIB).$(EXT) $(MKL_PATH)/$(RES_EXT)$(THREAD_LIB).$(EXT) $(MKL_PATH)/$(RES_EXT)$(CORE_LIB).$(EXT) -Wl,--end-group
else
    MKL_LIBS = -L$(MKL_PATH) -l$(IFACE_LIB) -l$(THREAD_LIB) -l$(CORE_LIB)
endif

OUTPUT_DIR=.

ifndef F95ROOT
   F95ROOT = $(MKLROOT)
   ifneq ($(compiler),intel)
     $(error F95ROOT is mandatory parameter for compiler==$(compiler))
   endif
endif


MODPATH=$(F95ROOT)/include/$(_IA)/$(interface)
LIB95=$(F95ROOT)/lib/$(_IA)/libmkl_lapack95_$(interface).a

#-------------------------------------------------------------------------------

NAME :=	NS_3D_Staggered_dp
SRCS := FD_functions.f90 csr_mod.f90 coo_mod.f90 other_utility.f90 NS_functions.f90 ogpf.f90 NS_3D_Staggered_dp.f90 sparse_utility.f90 main.f90
OBJS := $(SRCS:.f90=.o)
EXE  := $(NAME).out
vpath %.f90 source

all: $(EXE)

$(EXE) : $(OBJS)
	$(FC) $^ $(FOPTS) -I $(MODPATH) $(LIB95) $(LOPTS) $(MKL_LIBS) $(LL_THREAD_LIB) -lpthread -lm -ldl -o $(OUTPUT_DIR)/$@

%.o : %.f90 
	$(FC) $(FOPTS) -I $(MODPATH) -c $^ -o $@
	
clean :	
	rm -fr $(EXE) $(OBJS) *.mod
	
#-------------------------------------------------------------------------------
