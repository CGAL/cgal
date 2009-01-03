rounding_modes1: rounding_modes1.o
	g++ $^ -o $@

rounding_modes2: rounding_modes2.o
	g++ $^ -o $@ -lCGAL -L../../shared/lib

rounding_modes1.o: rounding_modes1.cpp
	g++ -c $< -o $@ -frounding-math

rounding_modes2.o: rounding_modes2.cpp
	g++ -c $< -o $@ -frounding-math -I../../include -I../../shared/include

clean:
	rm -f *.o rounding_modes1 rounding_modes2
