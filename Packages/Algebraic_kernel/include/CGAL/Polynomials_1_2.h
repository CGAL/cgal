#ifndef CGAL_ALGEBRAIC_KERNEL_POLYNOMIALS_1_2_H
#define CGAL_ALGEBRAIC_KERNEL_POLYNOMIALS_1_2_H

#include <CGAL/enum.h>

CGAL_BEGIN_NAMESPACE

template < typename RT_ >
class Polynomial_1_2
{
  RT_ rep[3]; // stores a, b, c for line ax+by+c=0
  
public:
  
  typedef RT_ RT;
  
  Polynomial_1_2(){}
  
  Polynomial_1_2(const RT & a, const RT & b, const RT & c)
  { 
    rep[0]=a;
    rep[1]=b;
    rep[2]=c;
  }

  const RT & a() const
  { return rep[0]; }

  const RT & b() const
  { return rep[1]; }

  const RT & c() const
  { return rep[2]; }
};

template < typename RT >
bool 
operator == ( const Polynomial_1_2<RT> & p1,
	      const Polynomial_1_2<RT> & p2 )
{
  return( (p1.a() == p2.a()) && 
              (p1.b() == p2.b()) &&
              (p1.c() == p2.c()) );
}
    
CGAL_END_NAMESPACE

#endif //CGAL_ALGEBRAIC_KERNEL_POLYNOMIALS_1_2_H
