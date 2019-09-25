INC := $(UNI10_ROOT_CPU)/include
LIB := $(UNI10_ROOT_CPU)/lib

CC := g++
FLAGS:= -g -m64 -DUNI_CPU -DUNI_LAPACK -std=c++11
#CC := icpc
#FLAGS:= -m64 -g -std=c++11 -DUNI_CPU -DUNI_LAPACK -DUNI_MKL 

EXU := main.e
OBJ := oneSiteiMPS.o ../nsy_opt/operator.o ../nsy_ham/hamiltonian.o ../tools/general_tools.o ../arpack_wrapper/arpack_wrapper.o 

all: $(OBJ) $(EXU)

$(OBJ): %.o: %.cpp
	$(CC) -c $< -o $@ -I$(INC) $(FLAGS) 

$(EXU): %.e: %.cpp ../tools/general_tools.h $(OBJ)
	#icpc scopion
	#$(CC) -I$(INC) -L$(LIB) -L$(ARP) $(FLAGS) -o $@ $< $(OBJ) -luni10_lapack_cpu $(ARP_LIB) $(MKL_LIB) $(FOR_LIB)
	#icpc boromir & beehive
	#$(CC) -I$(INC) -L$(LIB) -L$(ARP) $(FLAGS) -o $@ $< $(OBJ) -luni10_lapack_cpu $(ARP_LIB) $(MKL_LIB) -lgfortran
	#g++ 
	#$(CC) -o $@ $< -I$(INC) $(OBJ) -L$(LIB) $(ARP_LIB) $(FLAGS) -luni10_lapack_cpu -lgfortran 
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ $< $(OBJ) -luni10_lapack_cpu $(ARP_LIB) -lgfortran
	#g++  beehive
	#$(CC) -o $@ $< -I$(INC) $(OBJ) -L$(LIB) $(ARP_LIB) $(FLAGS) -luni10_lapack_cpu -lgfortran /home/pcs155251/atlas/liblapack.so.3

leak_check:
	valgrind --leak-check=yes ./$(EXU)

clean:
	rm -f *.o *.e $(OBJ)
