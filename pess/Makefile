UNI10_ROOT := /home/pcs155251/newLib
UNI10INC := $(UNI10_ROOT)/include
UNI10LIB := $(UNI10_ROOT)/lib

CC := g++
FLAGS:= -m64 -DUNI_CPU -DUNI_LAPACK -std=c++11
#CC := icpc
#FLAGS:= -m64 -g -std=c++11 -DUNI_CPU -DUNI_LAPACK -DUNI_MKL 

OBJO := pess.o ../ctmBase/ctmBase.o enPess.o enRedPess.o paraPess.o quant.o ../tools/general_tools.o ../nsy_ham/hamiltonian.o ../nsy_opt/operator.o ../parser/Parser.o
OBJE := main.o 
EXU  := main.e

all: $(OBJO) $(OBJE) $(EXU)

$(OBJO): %.o: %.cpp %.h
	$(CC) -c $< -o $@ -I$(UNI10INC) $(FLAGS) 

$(OBJE): %.o: %.cpp 
	$(CC) -c $< -o $@ -I$(UNI10INC) $(FLAGS) 

$(EXU) : %.e: %.o $(OBJO)
	$(CC) -o $@ $< $(OBJO) $(FLAGS) -L$(UNI10LIB) -luni10_lapack_cpu -Wl,-rpath,$(UNI10LIB)

clean:
	rm -f *.o *.e 
