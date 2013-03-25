// ============================================================================
//
// Copyright (c) 1997-2002 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: $
// release_date  : $CGAL_Date: $
//
// file          : test/Nef_3/nef_3.h
// package       : Nef_3
// chapter       : 3D-Nef Polyhedra
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Ralf Osbild        <osbild@mpi-sb.mpg.de>
// maintainer    : Peter Hachenberger <hachenberger@mpi-sb.mpg.de>
// coordinator   : MPI Saarbruecken
//
// ============================================================================


#define CGAL_NEF3_SORT_OUTPUT 1


#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/OFF_to_nef_3.h>
#include <CGAL/Timer.h>
#include <fstream>
#include <cassert>


typedef CGAL::Arithmetic_kernel::Integer NT;
typedef CGAL::Arithmetic_kernel::Rational FNT;


template<typename Kernel>
class test {

  typedef CGAL::Nef_polyhedron_3<Kernel>                    Nef_polyhedron;

private:
  static const char* datadir;  

  bool are_files_equal(const char* name1, const char* name2) {
    std::ifstream in1(name1);
    std::ifstream in2(name2);
    std::string s1;
    std::string s2;
    bool OK = true;
    while(in1 && OK) {
      in1 >> s1;
      in2 >> s2;
      if(s1 != s2) {
        std::cerr << s1 << std::endl;
	std::cerr << s2 << std::endl;
	OK = false;
      }
    }
    if(in2)
      OK = false;

    if(!OK) {
      char c;
      std::ifstream err1(name1);
      while(err1.get(c))
	std::cerr << c;
      std::ifstream err2(name2);
      while(err2.get(c))
	std::cerr << c;
    }
	
    return OK;
  }

  bool does_nef3_equals_file(Nef_polyhedron& N, const char* name, const char* suffix) {
    char* fullname = new char[std::strlen(datadir)+std::strlen(name)+std::strlen(suffix)+1];
    std::strcpy(fullname, datadir);
    std::strcat(fullname, name);
    std::strcat(fullname, suffix);
    std::ofstream out("data/temp.nef3");
    out << N;
    bool b = are_files_equal("data/temp.nef3",fullname);
    delete [] fullname;
    return b;
  }

  Nef_polyhedron built_nef_from_off(const char *name) {

     Nef_polyhedron N;

     char* fullname = new char[std::strlen(datadir)+std::strlen(name)+1];
     std::strcpy(fullname, datadir);
     std::strcat(fullname, name);

     std::ifstream off_file (fullname);
     assert(off_file != NULL);

     std::size_t discarded = CGAL::OFF_to_nef_3 (off_file, N, true);
     assert(discarded == 0);
     return N;
  }

public:
  void run_test(bool compare,const char* suffix) {
    Nef_polyhedron N = built_nef_from_off( "nine_planes.off");
    if(compare)
      assert(does_nef3_equals_file(N,"nine_planes.nef3",suffix));
  }
};

template<typename Kernel>
const char* test<Kernel>::datadir="data/";

int main() {

  CGAL::Timer t;
  t.start();

#if defined( CGAL_USE_LEDA ) || defined ( CGAL_USE_GMP )
  typedef CGAL::Homogeneous<NT>              H_kernel;
  typedef CGAL::Cartesian<FNT>               C_kernel;
  
  test<H_kernel>  test_H;
  test<C_kernel>  test_C;

  test_H.run_test(true,".H");
# ifdef CGAL_USE_GMP
  test_C.run_test(false,".C");
# else
  test_C.run_test(false,".LC");
# endif
#endif
  t.stop();
  std::cout << "Time " << t.time() << std::endl;
  return 0;
}

