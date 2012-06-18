user-config.jam
        Copy to tools/build/v2/

compile-boost
        Shell script to launch the compilation of Boost libraries with
        several configurations (different ABI).


The configurations are:
  - in stage/ for all g++ compilers with default options, gcc>=4.5,
  - in stage-intel/ for the Intel compilers,
  - in stage-cxxdebug/ for the g++ compilers with the STL debug mode
    (different ABI), for gcc>=4.6.

In addition, there is a configuration in stage-4.1/ for the g++-4.1
compiler (the one used on the Linux distribution RHEL 5).
