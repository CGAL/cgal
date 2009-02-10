include $(CGAL_MAKEFILE)

CXXFLAGS = -g \
           -I../include \
           -DCGAL_KINETIC_CHECK_EXPENSIVE\
           -DQT_THREAD_SUPPORT	\
           $(TESTSUITE_CXXFLAGS) \
           $(EXTRA_FLAGS) \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS)

LIBPATH = \
          $(TESTSUITE_LIBPATH) \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(TESTSUITE_LDFLAGS) \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_LDFLAGS)


DOTOS= file_0.o file_1.o file_2.o file_3.o file_4.o

lib: dotos
	ar cr ''libAOS3.a $(DOTOS)

dotos: $(DOTOS)

clean:
	rm -f *.o libAOS3.a
