.SUFFIXES : .cu 

CC	:= gcc
CXX 	:= g++
NVCC	:= $(NVCC_PATH)

COMMONFLAGS := -c 
# the "--no-align-double" is necessary. Without it the nvcc compiler aligns structure sizes every 8 bytes which is 
#incompatible with some NCBI_BLAST data structures, i.e. the "BlastContextInfo"
NVCCFLAGS :=  --compiler-options -fno-strict-aliasing --compiler-options -fno-inline -DUNIX -DCUDA55 -O3 --no-align-double -gencode=arch=compute_30,code=sm_30 -gencode=arch=compute_35,code=sm_35 -maxrregcount 56  -Xptxas -dlcm=cg -Xptxas -dscm=cs  


all: libhblast.a

libhblast.a : h_blast.o h_blast.cu.o h_blast_kernel.h 
	ar cru libhblast.a h_blast.o h_blast.cu.o
	ranlib libhblast.a	
	cp libhblast.a                          $(PROJ_DIR)/lib/
	cp $(CUDA_LIB)/libcudart.so               $(PROJ_DIR)/lib/

h_blast.o : h_blast.c h_blast.h h_blast_kernel.cu 
	$(CXX) $(COMMONFLAGS) -O3 h_blast.c -o h_blast.o -I$(CXX_DIR)/include/ -I$(PROJ_DIR)/inc/

h_blast.cu.o : h_blast.cu h_blast.h h_blast_kernel.cu h_blast_kernel.h 
	$(NVCC) $(COMMONFLAGS) $(NVCCFLAGS) h_blast.cu -o h_blast.cu.o -I$(CXX_DIR)/include/ -I$(PROJ_DIR)/inc/


