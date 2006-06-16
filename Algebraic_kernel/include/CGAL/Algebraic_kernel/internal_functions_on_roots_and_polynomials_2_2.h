// Copyright (c) 2003-2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s)     : Monique Teillaud, Sylvain Pion

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

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

    typedef typename AK::RT RT;
    typedef typename AK::Root_for_circles_2_2 Root_for_circles_2_2;
    RT dx = e2.a() - e1.a();
    RT dy = e2.b() - e1.b();

    RT dx2 = CGAL::square(dx);
    RT dy2 = CGAL::square(dy);
    RT dist2 = dx2 + dy2; // squared distance between centers

    RT cond = 4*e1.r_sq()*e2.r_sq() - 
      CGAL::square( e1.r_sq() + e2.r_sq() - dist2 );

    if (cond < 0) return res;

    RT px = e2.a() + e1.a();
    RT py = e2.b() + e1.b();
    RT rx1 = e1.r_sq() - CGAL::square(e1.a());
    RT ry1 = e1.r_sq() - CGAL::square(e1.b());
    RT rx2 = e2.r_sq() - CGAL::square(e2.a());
    RT ry2 = e2.r_sq() - CGAL::square(e2.b());
    
    RT drx = rx2 - rx1;
    RT dry = ry2 - ry1;
    
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
	 const typename AK::Root_for_circles_2_2 & r)
  {
    typedef typename AK::Root_of_2 Root_of_2;
    Root_of_2 part_left = square(r.x() - equation.a());
    Root_of_2 part_right = equation.r_sq() - square(r.y() - equation.b());
    Comparison_result c = compare(part_left, part_right);
    if(c == EQUAL) return ZERO;
    if(c == LARGER) return POSITIVE;
    return NEGATIVE;
  }


  template <class AK>
  typename AK::Root_for_circles_2_2
  x_critical_point(const typename AK::Polynomial_for_circles_2_2 & c, 
		   bool i)
  {
    typedef typename AK::Root_of_2            Root_of_2;
    typedef typename AK::FT                   FT;
    typedef typename AK::Root_for_circles_2_2 Root_for_circles_2_2;

    Root_of_2 a1= c.a() + make_root_of_2(FT(1),FT(0),-c.r_sq(),i);
    
    return Root_for_circles_2_2(a1, c.b());
  }

  template <class AK, class OutputIterator>
  OutputIterator
  x_critical_points(const typename AK::Polynomial_for_circles_2_2 & c, 
		    OutputIterator res)
  {
    typedef typename AK::Root_of_2            Root_of_2;
    typedef typename AK::FT                   FT;
    typedef typename AK::Root_for_circles_2_2 Root_for_circles_2_2;

    Root_of_2 a1= c.a() + make_root_of_2(FT(1),FT(0),-c.r_sq(),true);
    Root_of_2 a2= c.a() + make_root_of_2(FT(1),FT(0),-c.r_sq(),false);
    
    *res++ =  Root_for_circles_2_2(a1, c.b());
    *res++ =  Root_for_circles_2_2(a2, c.b());
    
    return res;
  }

  template <class AK>
  typename AK::Root_for_circles_2_2
  y_critical_point(const typename AK::Polynomial_for_circles_2_2 &c, 
		   bool i)
  {
    typedef typename AK::Root_of_2            Root_of_2;
    typedef typename AK::FT                   FT;
    typedef typename AK::Root_for_circles_2_2 Root_for_circles_2_2;

    Root_of_2 b1= c.b()+make_root_of_2(FT(1),FT(0),-c.r_sq(),i);

    return Root_for_circles_2_2(c.a(),b1);
  }
  
  template <class AK, class OutputIterator>
  OutputIterator
  y_critical_points(const typename AK::Polynomial_for_circles_2_2 & c, 
		    OutputIterator res)
  {
    typedef typename AK::Root_of_2 Root_of_2;
    typedef typename AK::FT        FT;
    typedef typename AK::Root_for_circles_2_2 Root_for_circles_2_2;
      
    Root_of_2 b1= c.b()+make_root_of_2(FT(1),FT(0),-c.r_sq(),true);
    Root_of_2 b2= c.b()+make_root_of_2(FT(1),FT(0),-c.r_sq(),false);
    
    *res++ = Root_for_circles_2_2(c.a(), b1);
    *res++ = Root_for_circles_2_2(c.a(), b2);

    return res;
  }



  } // namespace AlgebraicFunctors
} // namespace CGAL

#endif //  CGAL_ALGEBRAIC_KERNEL_FUNCTIONS_ON_ROOTS_AND_POLYNOMIALS_2_H

