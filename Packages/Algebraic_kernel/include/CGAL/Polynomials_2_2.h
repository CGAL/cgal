#ifndef CGAL_ALGEBRAIC_KERNEL_POLYNOMIALS_2_2_H
#define CGAL_ALGEBRAIC_KERNEL_POLYNOMIALS_2_2_H

//////////// FIXME - pb RT (cas general Polynomial_2_2) ou FT (ici)

#include <CGAL/enum.h>

namespace CGAL {

  // polynomials of the form (X-a)^2 + (Y-b)^2 - R^2
  template < typename FT_ >
    class Polynomial_for_circles_2_2
    {
      FT_ rep[3]; // stores a, b, R^2

    public:

      typedef FT_ FT;

      Polynomial_for_circles_2_2(){}

      Polynomial_for_circles_2_2(const FT & a,
				 const FT & b,
				 const FT & rsq)
	{ 
	  rep[0]=a;
	  rep[1]=b;
	  rep[2]=rsq;
	}

      const FT & a() const
	{ return rep[0]; }
      const FT & b() const
	{ return rep[1]; }
      const FT & r_sq() const
	{ return rep[2]; }

    };

  template < typename FT >
  bool 
  operator == ( const Polynomial_for_circles_2_2<FT> & p1,
		const Polynomial_for_circles_2_2<FT> & p2 )
  {
    return( (p1.a() == p2.a()) && 
	    (p1.b() == p2.b()) &&
	    (p1.r_sq() == p2.r_sq()) );
  }
    

} // namespace CGAL

#endif //CGAL_ALGEBRAIC_KERNEL_POLYNOMIALS_2_2_H
