CXX = g++
INCLUDE = -I./include
CXXFLAGS = $(INCLUDE) -O1

random_set_test: random_set_test.cc 
	$(CXX) $(CXXFLAGS) random_set_test.cc -o random_set_test

automated_test: automated_test.cc 
	$(CXX) $(CXXFLAGS) automated_test.cc -o automated_test




        #	$(CXX) -O3 -fomit-frame-pointer -fexpensive-optimizations -falign-functions=4 -malign-double -fprefetch-loop-arrays -march=pentium3 -msse automated_test.cc -o segment_tree
# -funroll-loops brings 10 percent slowdown on duron750 (64kb cache)
