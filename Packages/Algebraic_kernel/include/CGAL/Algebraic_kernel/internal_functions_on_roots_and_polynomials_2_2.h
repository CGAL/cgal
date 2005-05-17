#ifndef CGAL_ALGEBRAIC_KERNEL_FUNCTIONS_ON_ROOTS_AND_POLYNOMIALS_2_H
#define CGAL_ALGEBRAIC_KERNEL_FUNCTIONS_ON_ROOTS_AND_POLYNOMIALS_2_H

namespace CGAL {
  namespace AlgebraicFunctors {

  template < class AK, class OutputIterator >
  inline 
  OutputIterator
  solve( const typename AK::Polynomial_for_circles_2_2 & e1,
	 const typename AK::Polynomial_for_circles_2_2 & e2,
	 OutputIterator res )
  {
    assert( ! (e1 == e2) ); // polynomials of this type cannot be multiple 
    // of one another if they are not equal

    typedef typename AK::FT FT;

    FT dx = e2.a() - e1.a();
    FT dy = e2.b() - e1.b();

    FT dx2 = CGAL::square(dx);
    FT dy2 = CGAL::square(dy);
    FT dist2 = dx2 + dy2; // squared distance between centers

    FT cond = 4*e1.r_sq()*e2.r_sq() - 
      CGAL::square( e1.r_sq() + e2.r_sq() - dist2 );

    if (cond < 0) return res;

    FT px = e2.a() + e1.a();
    FT py = e2.b() + e1.b();
    FT rx1 = e1.r_sq() - CGAL::square(e1.a());
    FT ry1 = e1.r_sq() - CGAL::square(e1.b());
    FT rx2 = e2.r_sq() - CGAL::square(e2.a());
    FT ry2 = e2.r_sq() - CGAL::square(e2.b());
    
    FT drx = rx2 - rx1;
    FT dry = ry2 - ry1;
    
    if (cond == 0) {
      // one double root, 
      // no need to care about the boolean of the Root_of
      *res++ = std::make_pair
	( std::make_pair
	  (make_root_of_2(4*dist2,
			  4*(dx*drx-px*dy2),
			  CGAL::square(drx) - dy2*(2*(rx1+rx2)-dy2),
			  true),
	   make_root_of_2(4*dist2,
			  4*(dy*dry-py*dx2),
			  CGAL::square(dry) - dx2*(2*(ry1+ry2)-dx2),
			  true)),
	  2 ); // multiplicity = 2
      return res;
    }

    // else, 2 distinct roots

    bool slope = CGAL::sign(dx) * CGAL::sign(dy) <= 0 ;

    // fixme : use more clever manipulation of booleans...
    bool low_y = slope ? true : false;
    * res++ = std::make_pair
	( std::make_pair
	  (make_root_of_2(4*dist2,
			  4*(dx*drx-px*dy2),
			  CGAL::square(drx) - dy2*(2*(rx1+rx2)-dy2),
			  true),
	   make_root_of_2(4*dist2,
			  4*(dy*dry-py*dx2),
			  CGAL::square(dry) - dx2*(2*(ry1+ry2)-dx2),
			  low_y)),
	  1 );
    
    low_y = slope ? false : true;
    * res++ = std::make_pair
	( std::make_pair
	  (make_root_of_2(4*dist2,
			  4*(dx*drx-px*dy2),
			  CGAL::square(drx) - dy2*(2*(rx1+rx2)-dy2),
			  false),
	   make_root_of_2(4*dist2,
			  4*(dy*dry-py*dx2),
			  CGAL::square(dry) - dx2*(2*(ry1+ry2)-dx2),
			  low_y)),
	  1 );


   return res;
  }

  } // namespace AlgebraicFunctors
} // namespace CGAL

#endif //  CGAL_ALGEBRAIC_KERNEL_FUNCTIONS_ON_ROOTS_AND_POLYNOMIALS_2_H
