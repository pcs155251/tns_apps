UNI10ROOT := /home/pcs155251/newLib/
UNI10INC := $(UNI10ROOT)/include
UNI10LIB := $(UNI10ROOT)/lib
UNI10STAT := $(UNI10LIB)/libuni10.a

CC := g++
FLAGS:= -g -m64 -DUNI_CPU -DUNI_LAPACK -std=c++11
#CC := icpc
#FLAGS:= -m64 -g -std=c++11 -DUNI_CPU -DUNI_LAPACK -DUNI_MKL 

OBJO := oneSiteiMPS.o ../nsy_opt/operator.o ../nsy_ham/hamiltonian.o ../tools/general_tools.o ../arpack_wrapper/arpack_wrapper.o 
OBJE := main.o
EXU  := main.e

all: $(OBJO) $(OBJE) $(EXU)

$(OBJO): %.o: %.cpp %.h
	$(CC) -c $< -o $@ -I$(UNI10INC) $(FLAGS) 

$(OBJE): %.o: %.cpp 
	$(CC) -c $< -o $@ -I$(UNI10INC) $(FLAGS) 

$(EXU): %.e: %.o $(OBJO)
	$(CC) -o $@ $< $(OBJO) $(FLAGS) -L$(UNI10LIB) $(UNI10STAT) -llapack -lblas $(ARP_LIB) -lgfortran

clean:
	rm -f *.o *.e 
