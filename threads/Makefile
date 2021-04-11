# CFLAGS = -g -Wshadow -Wuninitialized -Warray-bounds -Wdiv-by-zero -std=c++14
# CFLAGS += -fsanitize=address -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith 
# CFLAGS += -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor
# CFLAGS += -Wcast-qual -Wno-suggest-attribute=format
# LDFLAGS= -fsanitize=address

CFLAGS = -O3
LDFLAGS= 

LIBS = -lm -lpthread

all: build

build: main.o matrix_operations.o solveLES.o
	g++ $(CFLAGS) main.o matrix_operations.o solveLES.o -o a.out $(LIBS) $(LDFLAGS)

main.o: main.cpp ./matrix_operations/matrix_operations.h
	g++ $(CFLAGS) main.cpp -c

matrix_operations.o: ./matrix_operations/matrix_operations.cpp ./matrix_operations/matrix_operations.h
	g++ $(CFLAGS) ./matrix_operations/matrix_operations.cpp -c

solveLES.o: ./matrix_operations/solveLES.cpp ./matrix_operations/matrix_operations.h
	g++ $(CFLAGS) ./matrix_operations/solveLES.cpp -c

clean:
	rm -rf *.o *.out build

  
