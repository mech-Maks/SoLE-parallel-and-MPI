# CFLAGS = -g -Wshadow -Wuninitialized -Warray-bounds -Wdiv-by-zero
# CFLAGS += -fsanitize=address -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith 
# CFLAGS += -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor
# CFLAGS += -Wno-suggest-attribute=format
# LDFLAGS= -fsanitize=address

CFLAGS = -O3
LDFLAGS =

# CFLAGS = -g
# LDFLAGS = 

COMPILER = mpicxx

all: build

build: main.o env_description.o
	$(COMPILER) main.o env_description.o -o a.out $(LDFLAGS)

main.o: main.cpp enviroment.h
	$(COMPILER) $(CFLAGS) main.cpp -c

env_description.o: env_description.cpp enviroment.h
	$(COMPILER) $(CFLAGS) env_description.cpp -c

clean:
	rm -rf *.o *.out build
