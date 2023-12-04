#MAC
COMMON=-O2 -I../../include -L../../bin -mavx -pthread
LIBS = -w -lmujoco200 -lnlopt -lglfw.3
CC = gcc

#LINUX
#COMMON=-O2 -I../../include -L../../bin -mavx -pthread -Wl,-rpath,'$$ORIGIN'
#LIBS = -lmujoco200 -lGL -lm -lnlopt -lglew ../../bin/libglfw.so.3
#CC = gcc

#Code to work with Windows is not ready yet

ROOT = projectile_opt

all:
	$(CC) $(COMMON) main.c $(LIBS) -o ../../bin/$(ROOT)

main.o:
	$(CC) $(COMMON) -c main.c

clean:
	rm *.o ../../bin/$(ROOT)
