CFLAGS = -O -DFOUR -m64

SHELL = /bin/sh
RM = /bin/rm
BIN = ../../bin


MakeLcircuit:	MakeLcircuit.o
	$(CC) -o MakeLcircuit $(CFLAGS) MakeLcircuit.o -lm
	mv MakeLcircuit $(BIN)/MakeLcircuit

ReadOutput: ReadOutput.o
	$(CC) -o ReadOutput $(CFLAGS) ReadOutput.o
	mv ReadOutput $(BIN)/ReadOutput

misc: MakeLcircuit ReadOutput

clean:
	$(RM) -f *.o

MakeLcircuit.o ReadOutput.o:	cmplx.h


