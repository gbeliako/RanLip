#############################################################################
#									    #
#CLASS LIBRARY RANLIP FOR MULTIVARIATE NONUNIFORM RANDOM VARIATE GENERATION #
#									    #
#	This makefile gives targets that show how to compile and link       #
#	user code to the ranlip shared library and statatic library.        #
#									    #
#############################################################################
#
# This make file show how to compile and link examples included with this 
# distribution of RANLIP assuming different installations of the library. this
# include the following examples for both static and shared linking.
#
# shared_example and staic_example:	
#
#			shows how to compile and link library when install
#			in the library search path used to load libraries
#
# static_example2:	
#			shows how to compile and link by implicitly telling
#			the linker where to look for the library
#
# static_example3: s	
#			hows how to compile and link procedural C conde.
#
############################################################################

# compiler
CC = g++

# Some options probably not needed: -g (which enables the debugger options).
FLAGS = -g -O -Wno-deprecated

# Object file fo the example
OBJ = ranliptest.o
OBJ2 = ranliptest.o
OBJ3 = ranliptestproc.o

# LIB_PATH used to store the path in which the library files were installed.
# The commented out assignment is for when the library is installed into the
# users home directory. NOTE: $(HOME) referes to env varialble HOME.

#LIB_PATH = $(HOME)/ranlip/lib/
LIB_PATH = /usr/local/lib/

# INCLUDE_PATH used to store the path in which the *.h files have been
# placed. The commented out assignment is for when the *.h files are placed
# in the users home directory.

#INCLUDE_PATH = $(HOME)/ranlip/include/
INCLUDE_PATH = /usr/local/include/

all:	static_example2 static_example3 shared_example static_example 

#################################################################################
# linking ranliptest.o . If you have succesfully installed lip library and have
# LIB_PATH to /etc/ld.soconf Or you have added LIB_PATH TO LD_LIBRARY_PATH
# then compiling is as eassy as this.  

# shared_example target links ranliptest.o to the lip shared library. To make
# up shared_example executable.

shared_example:	$(OBJ)
		$(CC) -o shared_example $(OBJ) $(FLAGS) -lranlip -lm

# static_example target links ranliptest.o to the lip static library. To make
# up static_example executable.

static_example: $(OBJ)
		$(CC) -o static_example -static  $(OBJ) $(FLAGS)  -lranlip -lm


#################################################################################
# linking ranliptest.o . 

# static_example target links ranliptest.o to the lip static library. To make
# up static_example executable.

static_example2:$(OBJ2)
		$(CC) -o static_example2  $(OBJ2) $(FLAGS) $(LIB_PATH)libranlip.a -lm


#################################################################################
# linking ranlipporc.o

# static_example target links ranliptestproc.o to the lip static library. To make
# up static_example executable.

static_example3:$(OBJ3)
		$(CC) -o static_example3 -static $(OBJ3) $(FLAGS) $(LIB_PATH)libranlip.a -lm


#################################################################################
# compiling examples to objectfiles.

ranliptest.o:			ranliptest.cpp
				$(CC) -c ranliptest.cpp $(FLAGS) -I$(INCLUDE_PATH)

# compiling proccedual example using C compiler
#
ranliptestproc.o:		ranliptestproc.cpp
				gcc -c ranliptestproc.cpp $(FLAGS) -I$(INCLUDE_PATH)

.PHONY: clean all

clean:
		rm -f $(OBJ) shared_example static_example
		rm -f $(OBJ2) static_example2
		rm -f $(OBJ3) static_example3
		rm -f partition.txt
