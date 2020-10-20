
SRC = ./spin_chain.cpp
HDR = ./spin_chain.hpp

LIB = -lm -lboost_program_options-mt -lopenblas -lm

CC = g++ -std=c++11 -O2 -fopenmp -I./ -L./

RHOST = pbuividovich@134.176.18.144:/home/pbuividovich/ExDiag/

xxz: main.cpp $(SRC) $(HDR)
	$(CC) ./$< $(SRC) $(LIB) -o ./$@.exe
	
clean:
	rm -f -v ./xxz.exe
	
upload:
	dos2unix -k ./Makefile ./*.sh
	rsync -v -R -z -t -r ./Makefile ./*.sh ./main.cpp $(SRC) $(HDR) $(RHOST)
