#include <CGAL/Gmpfr.h>

int main(){

  CGAL::Gmpfr a(2);
  CGAL::Gmpfr b(3);

  // clear inexactitude flags and perform an exact operation
  CGAL::Gmpfr::clear_flags();
  std::cout<<a<<" + "<<b<<" = "<<a+b<<std::endl;
  std::cout<<"inexact flag: "<<CGAL::Gmpfr::inex_flag()<<std::endl;

  // clear flags again and perform an inexact operation
  CGAL::Gmpfr::clear_flags();
  std::cout<<a<<" / "<<b<<" = "<<a/b<<std::endl;
  std::cout<<"inexact flag: "<<CGAL::Gmpfr::inex_flag()<<std::endl;

  return 0;
}
