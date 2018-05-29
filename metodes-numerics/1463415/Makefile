OPT= -g
EXECS=integralExp integralAbs integralElipse

all : $(EXECS)

# Executables
arrels : arrels.o polinomis.o
	gcc -o arrels arrels.o polinomis.o -lm

integralExp : integralExp.o polinomis.o
	gcc -o integralExp integralExp.o polinomis.o -lm

integralAbs : integralAbs.o polinomis.o
	gcc -o integralAbs integralAbs.o polinomis.o -lm

integralElipse : integralElipse.o polinomis.o
	gcc -o integralElipse integralElipse.o polinomis.o -lm
	
# Objectes
polinomis.o : polinomis.c
	gcc -c $(OPT) polinomis.c 

integralExp.o : integralExp.c
	gcc -c $(OPT) integralExp.c 

integralAbs.o : integralAbs.c
	gcc -c $(OPT) integralAbs.c 

integralElipse.o : integralElipse.c
	gcc -c $(OPT) integralElipse.c 

arrels.o : arrels.c
	gcc -c $(OPT) arrels.c 

# Neteja
clean :
	rm -f *.o

realclean : clean
	rm -f $(EXECS)
