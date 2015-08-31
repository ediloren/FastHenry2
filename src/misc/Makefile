all: MakeLcircuit ReadOutput

MakeLcircuit: MakeLcircuit.c
	$(CC) -O -o MakeLcircuit MakeLcircuit.c -lm

ReadOutput: ReadOutput.c
	$(CC) -O -o ReadOutput ReadOutput.c -lm

clean::
	-@rm -f *.o MakeLcircuit ReadOutput
