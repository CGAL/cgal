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
    typedef typename AK::Root_for_circles_2_2 Root_for_circles_2_2;
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
	( Root_for_circles_2_2
	  (make_root_of_2(4*dist2,
			  4*(dx*drx-px*dy2),
			  CGAL::square(drx) - dy2*(2*(rx1+rx2)-dy2),
			  true),
	   make_root_of_2(4*dist2,
			  4*(dy*dry-py*dx2),
			  CGAL::square(dry) - dx2*(2*(ry1+ry2)-dx2),
			  true)),
	  static_cast<unsigned>(2) ); // multiplicity = 2
      return res;
    }

    // else, 2 distinct roots

    bool slope = CGAL::sign(dx) * CGAL::sign(dy) <= 0 ;

    // fixme : use more clever manipulation of booleans...
    bool low_y = slope ? true : false;
    * res++ = std::make_pair
	( Root_for_circles_2_2
	  (make_root_of_2(4*dist2,
			  4*(dx*drx-px*dy2),
			  CGAL::square(drx) - dy2*(2*(rx1+rx2)-dy2),
			  true),
	   make_root_of_2(4*dist2,
			  4*(dy*dry-py*dx2),
			  CGAL::square(dry) - dx2*(2*(ry1+ry2)-dx2),
			  low_y)),
	  static_cast<unsigned>(1) );
    
    low_y = slope ? false : true;
    * res++ = std::make_pair
	( Root_for_circles_2_2
	  (make_root_of_2(4*dist2,
			  4*(dx*drx-px*dy2),
			  CGAL::square(drx) - dy2*(2*(rx1+rx2)-dy2),
			  false),
	   make_root_of_2(4*dist2,
			  4*(dy*dry-py*dx2),
			  CGAL::square(dry) - dx2*(2*(ry1+ry2)-dx2),
			  low_y)),
	  static_cast<unsigned>(1) );


   return res;
  }

   template < class AK >
    inline 
    Sign sign_at( const typename AK::Polynomial_for_circles_2_2 & equation,
		  const typename AK::Root_for_circles_2_2 r){
    typedef typename AK::Root_of_2 Root_of_2;
    Root_of_2 part_left = square(r.x() - equation.a());
    Root_of_2 part_right = equation.r_sq() - square(r.y() - equation.b());
    if (part_left == part_right)
      return ZERO;
    return (part_left < part_right) ? NEGATIVE : POSITIVE; 
  }


  template <class AK>
  typename AK::Root_for_circles_2_2
  x_critical_points(const typename AK::Polynomial_for_circles_2_2 & c, bool i)
  {
            typedef typename AK::Root_of_2 Root_of_2;
            typedef typename AK::FT        FT;
            typedef typename AK::Root_for_circles_2_2
                                           Root_for_circles_2_2;

	    Root_of_2 a1= c.a() + make_root_of_2(FT(1),FT(0),-c.r_sq(),i);

            return Root_for_circles_2_2(a1, make_root_of_2(FT(1), FT(0),
                                   -square(c.b()), c.b()<0));
  }

  template <class AK, class OutputIterator>
  OutputIterator
  x_critical_points(const typename AK::Polynomial_for_circles_2_2 & c, OutputIterator res)
  {
            typedef typename AK::Root_of_2 Root_of_2;
            typedef typename AK::FT        FT;
            typedef typename AK::Root_for_circles_2_2
                                           Root_for_circles_2_2;

	    Root_of_2 a1= c.a() + make_root_of_2(FT(1),FT(0),-c.r_sq(),true);
	    Root_of_2 a2= c.a() + make_root_of_2(FT(1),FT(0),-c.r_sq(),false);

            *res++ =  Root_for_circles_2_2(a1, make_root_of_2(FT(1), FT(0),
							      -square(c.b()), c.b()<0));
	    *res++ =  Root_for_circles_2_2(a2, make_root_of_2(FT(1), FT(0),
							      -square(c.b()), c.b()<0));
	    return res;
  }



  template <class AK>
  typename AK::Root_for_circles_2_2
  y_critical_points(const typename AK::Polynomial_for_circles_2_2 &c, bool i)
  {
            typedef typename AK::Root_of_2 Root_of_2;
            typedef typename AK::FT        FT;
            typedef typename AK::Root_for_circles_2_2
                                           Root_for_circles_2_2;

            Root_of_2 b1= c.b()+make_root_of_2(FT(1),FT(0),-c.r_sq(),i);

            return Root_for_circles_2_2(make_root_of_2(FT(1),FT(0),
                                   -square(c.a()),c.a()<0),b1);
  }
  
  template <class AK, class OutputIterator>
    OutputIterator
    y_critical_points(const typename AK::Polynomial_for_circles_2_2 & c, OutputIterator res)
    {
      typedef typename AK::Root_of_2 Root_of_2;
      typedef typename AK::FT        FT;
      typedef typename AK::Root_for_circles_2_2
	Root_for_circles_2_2;
      
      Root_of_2 b1= c.b()+make_root_of_2(FT(1),FT(0),-c.r_sq(),true);
      Root_of_2 b2= c.b()+make_root_of_2(FT(1),FT(0),-c.r_sq(),false);
      
      *res++ = Root_for_circles_2_2(make_root_of_2(FT(1),FT(0),
						   -square(c.a()),c.a()<0),b1);
      *res++ = Root_for_circles_2_2(make_root_of_2(FT(1),FT(0),
						   -square(c.a()),c.a()<0),b2);
      return res;
  }



  } // namespace AlgebraicFunctors
} // namespace CGAL

#endif //  CGAL_ALGEBRAIC_KERNEL_FUNCTIONS_ON_ROOTS_AND_POLYNOMIALS_2_H
