#############################################################################
#									    #
#CLASS LIBRARY RANLIP FOR MULTIVARIATE NONUNIFORM RANDOM VARIATE GENERATION #
#	MAPLE 2017 version	on OSX							    #
#									    #
############################################################################

# compiler

CC = c++
CCC = cc
LL = c++

# Some options probably not needed: -g (which enables the debugger options).
FLAGS =  -O -Wno-deprecated -fPIC -g -m64 -mmacosx-version-min=10.8
LDFLAGS=  -lc -lmaplec  -lhf -lmaple -arch x86_64 -shared

# Object file fo the example
OBJ = discrete.o
OBJ2 = ranlip1.o
OBJ3 = ranlipproc.o
OBJ4 = mwrap_ranlip.o

# LIB_PATH used to store the path in which the library files were installed.
# The commented out assignment is for when the library is installed into the
# users home directory. NOTE: $(HOME) refers to env variable HOME.

#LIB_PATH = $(HOME)/ranlip/lib/
LIB_PATH = /Library/Frameworks/Maple.framework/Versions/2017/bin.APPLE_UNIVERSAL_OSX/

# INCLUDE_PATH used to store the path in which the *.h files have been
# placed. The commented out assignment is for when the *.h files are placed
# in the users home directory.

#INCLUDE_PATH = $(HOME)/ranlip/include/
INCLUDE_PATH = /Library/Frameworks/Maple.framework/Versions/2017/extern/include/

all:	mapleranlip.so



# shared_example target links ranliptest.o to the lip shared library. To make
# up shared_example executable.

mapleranlip.so:	$(OBJ) $(OBJ4) $(OBJ2) $(OBJ3)
		$(LL)  $(OBJ) $(OBJ4) $(OBJ2) $(OBJ3) -o $@ $(LDFLAGS) -L$(LIB_PATH) -v



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
		
