#include <CGAL/Gmpfr.h>

std::ostream& print_with_prec(std::ostream &o,const CGAL::Gmpfr &x){
  return (o<<x<<" prec "<<x.get_prec());
}

int main(){

  CGAL::Gmpfr a(200,5);
  CGAL::Gmpfr b(300,10);
  CGAL::Gmpfr c;  // constructed with default precision

  std::cout<<"a=";
  print_with_prec(std::cout,a);
  std::cout<<"\nb=";
  print_with_prec(std::cout,b);

  c=a+b;  // c has now the precision of a
  std::cout<<"\nc=a+b=";
  print_with_prec(std::cout,c);

  CGAL::Gmpfr d(a,7);  // d stores the value of a with precision 7
  d+=b;
  std::cout<<"\n\nd=a\n(d+=b)=";
  print_with_prec(std::cout,d);
  std::cout<<std::endl;

  return 0;
}
