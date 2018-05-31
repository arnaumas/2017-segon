OPT= -g
EXECS=integralExp integralAbs integralElipse

all : $(EXECS)

# Executables
integralExp : integralExp.o polinomis.o matrius.o
	gcc -o integralExp integralExp.o polinomis.o matrius.o -lm

integralAbs : integralAbs.o polinomis.o matrius.o
	gcc -o integralAbs integralAbs.o polinomis.o matrius.o -lm

integralElipse : integralElipse.o polinomis.o matrius.o
	gcc -o integralElipse integralElipse.o polinomis.o matrius.o -lm

# Objectes
polinomis.o : polinomis.c
	gcc -c $(OPT) polinomis.c 

matrius.o : matrius.c
	gcc -c $(OPT) matrius.c 

integralExp.o : integralExp.c
	gcc -c $(OPT) integralExp.c 

integralAbs.o : integralAbs.c
	gcc -c $(OPT) integralAbs.c 

integralElipse.o : integralElipse.c
	gcc -c $(OPT) integralElipse.c 

# Neteja
clean :
	rm -f *.o

realclean : clean
	rm -f $(EXECS)
