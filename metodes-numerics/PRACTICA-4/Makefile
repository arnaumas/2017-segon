OPT=-g -Wall -pedantic
EXECS=arrels integralExp

all : $(EXECS)

# Executables
integralExp : integralExp.o polinomis.o
	gcc -o integralExp integralExp.o polinomis.o -lm

arrels : arrels.o polinomis.o
	gcc -o arrels arrels.o polinomis.o -lm

# Objectes
polinomis.o : polinomis.c
	gcc -c $(OPT) polinomis.c 

integralCheb.o : integralCheb.c
	gcc -c $(OPT) integralCheb.c 

arrels.o : arrels.c
	gcc -c $(OPT) arrels.c 

# Neteja
clean :
	rm -f *.o

realclean : clean
	rm -f $(EXECS)
