// Author(s)     : Michael Hemmer <mhemmer@uni-mainz.de>

#include <CGAL/basic.h>
#include <cassert>
#include <CGAL/Residue.h>
#include <CGAL/primes.h>
#include <CGAL/Test/_test_algebraic_structure.h>

#ifdef _OPENMP
#include <omp.h>
// This file needs Open MP.
// Use ,e.g. , gcc-4.3.1 with -fopenmp and -lgomp
int main ()  {
  int tid;
  
  //Beginning of parallel section. Fork a team of threads.
  //Specify variable scoping   
#pragma omp parallel private(tid)
  {
    // Enforce IEEE double precision and rounding mode to nearest
    CGAL::Protect_FPU_rounding<true> pfr(CGAL_FE_TONEAREST);
    tid = omp_get_thread_num();
    int old_prime = CGAL::internal::primes[0];
    int new_prime = CGAL::internal::primes[tid+1];
    assert(CGAL::Residue::get_current_prime() == old_prime);
    CGAL::Residue::set_current_prime(new_prime);
    assert(CGAL::Residue::get_current_prime() == new_prime);
    
    typedef CGAL::Residue NT;
    typedef CGAL::Field_tag Tag;
    typedef CGAL::Tag_true Is_exact;
    CGAL::test_algebraic_structure<NT,Tag, Is_exact>();    
  }  
}
#else
int main (){}
#endif
