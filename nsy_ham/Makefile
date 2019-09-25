INC := $(TUNI10_ROOT_CPU)/include
LIB := $(TUNI10_ROOT_CPU)/lib

CC := icpc
FLAGS:= -m64 -O2 -std=c++11 -DUNI_CPU -DUNI_LAPACK -DUNI_MKL 
#CC := g++
#FLAGS:= -m64 -std=c++11 -DUNI_CPU -DUNI_LAPACK
STATIC_LIB := $(LIB)/libuni10.a

EXU := main
OBJ := ../nsy_opt/operator.o hamiltonian.o 

all: $(OBJ) $(EXU)

$(OBJ): %.o: %.cpp
	$(CC) -c -I$(INC) $(FLAGS) $< -o $@

$(EXU): $(EXU).cpp
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ $< $(OBJ) -luni10_lapack_cpu

leak_check:
	valgrind --leak-check=yes ./$(EXU)

clean:
	rm -f $(EXU) *.o
