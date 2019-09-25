ARP_LIB=/home/pcs155251/ARPACK/libarpack.a
UNI10INC := $(UNI10_ROOT_CONDA)/include
UNI10LIB := $(UNI10_ROOT_CONDA)/lib

CC := g++
FLAGS:= -m64 -DUNI_CPU -DUNI_LAPACK -std=c++11
#CC := icpc
#FLAGS:= -m64 -O2 -std=c++11 -DUNI_CPU -DUNI_LAPACK -DUNI_MKL 

OBJO := ../nsy_opt/operator.o ../nsy_ham/hamiltonian.o ../tools/general_tools.o ../arpack_wrapper/arpack_wrapper.o ../oneSiteiMPS/oneSiteiMPS.o channel.o
OBJE := envir.o gradient.o cg.o load.o test.o
EXU  := envir.e gradient.e cg.e load.e test.e

all: $(OBJO) $(OBJE) $(EXU)
	
$(OBJO): %.o: %.cpp %.h
	$(CC) -c $< -o $@ -I$(UNI10INC) $(FLAGS)

$(OBJE): %.o: %.cpp
	$(CC) -c $< -o $@ -I$(UNI10INC) $(FLAGS)

$(EXU): %.e: %.o $(OBJO)
	$(CC) -o $@ $< $(OBJO) $(FLAGS) -L$(UNI10LIB) -I$(UNI10INC) -luni10_lapack_cpu $(ARP_LIB) -lgfortran -Wl,-rpath,$(UNI10LIB)

clean:
	rm -f *.o *.e
