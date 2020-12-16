SRC = ./linalg.cpp ./spin_chain.cpp

HDR = $(SRC:.cpp=.hpp)
HDR += ./timing.hpp  

LIB = -lm -lopenblas -lm

CC = g++ -O2 -fopenmp -I./ -I../clue/include/ -L./
CC += -DHAVE_LAPACK_CONFIG_H -DLAPACK_COMPLEX_CPP

#For the Gaia system in Giessen
ifeq ($(USER),pbuividovich)
CC += -std=c++14 -I../ -L../OpenBLAS/ -Wl,-R../OpenBLAS/
LIB += -lboost_program_options
endif

ifeq ($(USER),LocalAdmin)
LIB += -std=c++14 -lboost_program_options-mt
endif

RHOST = pbuividovich@134.176.18.144:/home/pbuividovich/ExDiag/

trans: transition.cpp $(SRC) $(HDR)
	$(CC) ./$< $(SRC) $(LIB) -o ./$@
	
trotter: trotter.cpp $(SRC) $(HDR)
	$(CC) $(SRC) ./$< $(LIB) -o ./$@
	
clean:
	rm -f -v ./trans ./trotter
	
RHOST = pbuividovich@134.176.18.144:~/ExDiag/
#RHOST = pbuivid@barkla5.liv.ac.uk:~/ExDiag/

FLIST =  ./Makefile $(SRC) $(HDR)
FLIST += ./transition.cpp ./trotter.cpp
FLIST += ./scripts/*.sh 
	
upload:
	dos2unix -k ./Makefile ./scripts/*.sh
	rsync -v -R -z -t -r $(FLIST) $(RHOST)
	

