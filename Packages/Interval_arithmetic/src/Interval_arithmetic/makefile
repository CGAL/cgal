
CFLAGS= -O2 -g

# all: workaround_4_irix.o workaround_4_irix_static.o workaround_4_irix_dynamic.o
all: workaround_4_irix.s

clean:
	rm -f *.o *.obj *.s

workaround_4_ms.obj: workaround_4_ms.c
	cl -c -O2 -nologo workaround_4_ms.c

workaround_4_irix.s: workaround_4_irix.c
	gcc -S workaround_4_irix.c

workaround_4_irix.o: workaround_4_irix_dynamic.o
	cp workaround_4_irix_dynamic.o workaround_4_irix.o

workaround_4_irix_static.o: workaround_4_irix.c
	gcc $(CFLAGS)       -c workaround_4_irix.c -o workaround_4_irix_static.o

workaround_4_irix_dynamic.o: workaround_4_irix.c
	gcc $(CFLAGS) -fpic -c workaround_4_irix.c -o workaround_4_irix_dynamic.o
