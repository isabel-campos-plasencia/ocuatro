#CFLAGS= -g -Wall                             #para debug  
CFLAGS=  -O2 -funroll-loops -finline-functions #optimizacion maxima

OBJ   = o4ini.o o4.o o4med.o o4upd.o 

o4: $(OBJ)
	gcc $(CFLAGS) $(OBJ) -lm  -o o4
.c.o:
	gcc -DSUN4 -c $(CFLAGS) $<

clean: 
	/bin/rm -f $(OBJ) o4

#		*Individual File Dependencies*
o4.o: o4.c o4.h o4aritm.h 

o4ini.o: o4ini.c  o4.h o4aritm.h

o4med.o: o4med.c o4.h o4aritm.h

o4upd.o: o4upd.c o4.h o4aritm.h



