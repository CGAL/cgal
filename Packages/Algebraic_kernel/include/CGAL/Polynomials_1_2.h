#ifndef CGAL_ALGEBRAIC_KERNEL_POLYNOMIALS_1_2_H
#define CGAL_ALGEBRAIC_KERNEL_POLYNOMIALS_1_2_H

//////////// FIXME - pb RT (cas general Polynomial_1_2) ou FT (ici)

#include <CGAL/enum.h>

namespace CGAL {

  template < typename FT_ >
    class Polynomial_1_2
    {
      FT_ rep[3]; // stores a, b, c for line ax+by+c=0

    public:

      typedef FT_ FT;

      Polynomial_1_2(){}

      Polynomial_1_2(const FT & a,
		     const FT & b,
		     const FT & c)
	{ 
	  rep[0]=a;
	  rep[1]=b;
	  rep[2]=c;
	}

      const FT & a() const
	{ return rep[0]; }
      const FT & b() const
	{ return rep[1]; }
      const FT & c() const
	{ return rep[2]; }

    };

  template < typename FT >
  bool 
  operator == ( const Polynomial_1_2<FT> & p1,
		const Polynomial_1_2<FT> & p2 )
  {
    return( (p1.a() == p2.a()) && 
	    (p1.b() == p2.b()) &&
	    (p1.c() == p2.c()) );
  }
    

} // namespace CGAL

#endif //CGAL_ALGEBRAIC_KERNEL_POLYNOMIALS_1_2_H
