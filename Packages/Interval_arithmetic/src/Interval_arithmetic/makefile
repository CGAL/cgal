
CFLAGS= -O2 -g

all: workaround_4_irix.o workaround_4_irix_static.o workaround_4_irix_dynamic.o

clean:
	rm -f *.o *.s

workaround_4_irix.o: workaround_4_irix_dynamic.o
	cp $< $@

workaround_4_irix_static.o: workaround_4_irix.c
	gcc $(CFLAGS)       -c $< -o $@

workaround_4_irix_dynamic.o: workaround_4_irix.c
	gcc $(CFLAGS) -fpic -c $< -o $@
