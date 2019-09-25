CC := g++
FLAGS:= -m64 --std=c++11 -g

EXU := example.e 
OBJ := Parser.o

all: $(EXU)

$(OBJ): %.o: %.cpp
	$(CC) -c $(FLAGS) $< -o $@
	
$(EXU): %.e: %.cpp $(OBJ)
	$(CC) $(FLAGS) -o $@ $< $(OBJ)

leak_check:
	valgrind --leak-check=yes ./$(EXU)

clean:
	rm -f *.o *.e
