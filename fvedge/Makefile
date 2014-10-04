#------------------------------------------
# Author:       Andrew Kercher
# Program:      Makefile
# Last Updated: 06/17/2014 
# Description:  Makefile for 2D MHD Eqs.
#------------------------------------------

#-------------------------------------
# Source File
#-------------------------------------
SRC1 = src/main.cu
SRC2 = src/main.cpp

#-------------------------------------
# Target
#-------------------------------------
EXE1 = gpu.run
EXE2 = omp.run
TST = test.run

#-------------------------------------
# Definitions
#-------------------------------------
REAL = float
THR = THRUST_DEVICE_BACKEND_OMP

#----------------------------------------
# LIBRARIES AND COMPILER TO USE ON LINUX
#----------------------------------------
CMP1 = nvcc
CMP2 = g++ 
LOC = /home/akercher/
LIB = gomp

#----------------------------------------
# Compile and all targets 
#----------------------------------------
gpu_hd : $(EXE1)
	./$? bin/config.prob
omp_hd : $(EXE2)
	./$? bin/config.prob

gpu_mhd : $(EXE3)
	./$? bin/config.prob
omp_mhd : $(EXE4)
	./$? bin/config.prob

test : $(TST)
	./$? #bin/config.prob


#---------------------------------------------
# Dependencies
#---------------------------------------------
$(EXE1): $(SRC1)
	$(CMP1) $? -o $@ -O3 -DReal=$(REAL) -I$(LOC)

$(EXE2): $(SRC2)
	$(CMP2) $? -o $@ -O3 -fopenmp -DReal=$(REAL) -DTHRUST_DEVICE_SYSTEM=$(THR) -I$(LOC)

$(EXE3): $(SRC1)
	$(CMP1) $? -o $@ -O3 -DMHD -DReal=$(REAL) -I$(LOC)

$(EXE4): $(SRC2)
	$(CMP2) $? -o $@ -O3 -fopenmp -DMHD -DReal=$(REAL) -DTHRUST_DEVICE_SYSTEM=$(THR) -I$(LOC)

$(TST): src/test.cpp
	$(CMP2) $? -o $@ -O3 -fopenmp -DReal=$(REAL) -DTHRUST_DEVICE_SYSTEM=$(THR) -l$(LIB) -I$(LOC)

#-------------------------------
# CLEAN EXECUTABLES AND OBJECTS
#-------------------------------
clean:
	rm -f $(EXE)
	rm -f $(TST)
	rm -f *.run
