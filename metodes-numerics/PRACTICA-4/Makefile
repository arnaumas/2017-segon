OPT=-g -Wall -pedantic
EXECS=arrels

all : $(EXECS)

# Executables
arrels : arrels.o polinomis.o
	gcc -o arrels arrels.o polinomis.o -lm

# Objectes
polinomis.o : polinomis.c
	gcc -c $(OPT) polinomis.c 

arrels.o : arrels.c
	gcc -c $(OPT) arrels.c 

# Neteja
clean :
	rm -f *.o

realclean : clean
	rm -f $(EXECS)
