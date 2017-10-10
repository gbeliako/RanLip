#############################################################################
#									    #
#CLASS LIBRARY RANLIP FOR MULTIVARIATE NONUNIFORM RANDOM VARIATE GENERATION #
#	MAPLE version								    #
#	This makefile gives targets that show how to compile and link       #
#	user code to the ranlip shared library and statatic library.        #
#									    #
############################################################################

# compiler
CC = c++
CCC = cc
LL = cc

# Some options probably not needed: -g (which enables the debugger options).
FLAGS =  -O -Wno-deprecated -fPIC -g
LDFLAGS= -shared -lc -lmaplec -lrt

# Object file fo the example
OBJ = discrete.o
OBJ2 = ranlip1.o
OBJ3 = ranlipproc.o
OBJ4 = mwrap_ranlip.o

# LIB_PATH used to store the path in which the library files were installed.
# The commented out assignment is for when the library is installed into the
# users home directory. NOTE: $(HOME) referes to env varialble HOME.

#LIB_PATH = $(HOME)/ranlip/lib/
LIB_PATH = /opt/maple2017/bin.X86_64_LINUX/

# INCLUDE_PATH used to store the path in which the *.h files have been
# placed. The commented out assignment is for when the *.h files are placed
# in the users home directory.

#INCLUDE_PATH = $(HOME)/ranlip/include/
INCLUDE_PATH = /opt/maple2017/extern/include/

all:	mapleranlip.so



# shared_example target links ranliptest.o to the lip shared library. To make
# up shared_example executable.

mapleranlip.so:	$(OBJ) $(OBJ4) $(OBJ2) $(OBJ3)
		$(LL) $(FLAGS) $(OBJ4) $(OBJ) $(OBJ2) $(OBJ3) -o $@ $(LDFLAGS)  -L$(LIB_PATH) -v



# compiling examples to objectfiles.

discrete.o:			discrete.c
				$(CC) -c discrete.c $(FLAGS) -I$(INCLUDE_PATH)

ranlip1.o:			ranlip1.cpp
				$(CC) -c ranlip1.cpp $(FLAGS) -I$(INCLUDE_PATH)

ranlipproc.o:			ranlipproc.cpp
				$(CC) -c ranlipproc.cpp $(FLAGS) -I$(INCLUDE_PATH)

mwrap_ranlip.o:			 mwrap_ranlip.c
				$(CCC) -c  mwrap_ranlip.c $(FLAGS) -I$(INCLUDE_PATH)


.PHONY: clean all

clean:
		rm -f $(OBJ) 
		rm -f $(OBJ2)
		rm -f $(OBJ3)
		rm -f $(OBJ4)
		
