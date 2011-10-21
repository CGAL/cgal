// Copyright (c) 2010 CNRS, LIRIS, http://liris.cnrs.fr/, All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_LINEAR_CELL_COMPLEX_TRAITS_H
#define CGAL_LINEAR_CELL_COMPLEX_TRAITS_H 1

#include <CGAL/Cartesian.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

namespace CGAL {

  /** Trait class for Linear_cell_complex class.
   *  dD version (for the moment there is only one dD kernel in CGAL).
   */
  template <unsigned int d_, class Kernel>
  struct Linear_cell_complex_traits : public Kernel
  {
    static const unsigned int ambient_dimension = d_;
    
    typedef typename Kernel::FT          FT;
    typedef typename Kernel::Point_d     Point;
    typedef typename Kernel::Vector_d    Vector;
    typedef typename Kernel::Direction_d Direction;
    
    // Constructions
    struct Construct_translated_point
    {
      Point operator() (const Point&p, const Vector& v)
      { return p+v; }
    };

    // TODO THE Construct_vector
    struct Construct_vector : public Kernel::Construct_vector_d
    {
      using Kernel::Construct_vector_d::operator(); 
      Vector operator() (typename Kernel::FT x1)
      {
        Vector v(d_, NULL_VECTOR); v[0]=x1;
        return v;
      }
      Vector operator() (typename Kernel::FT x1, typename Kernel::FT x2)
      {
        Vector v(d_, NULL_VECTOR); v[0]=x1; v[1]=x2;
        return v;
      }
      Vector operator() (typename Kernel::FT x1, 
			 typename Kernel::FT x2, 
			 typename Kernel::FT x3)
      {
        Vector v(d_, NULL_VECTOR); v[0]=x1; v[1]=x2; v[2]=x3;
        return v;
      }
      Vector operator() (const Origin&, const Point& p)
      { return typename Kernel::Point_to_vector_d()(p); }
    };

    struct Construct_sum_of_vectors
    {
      Vector operator() (const Vector&v1, const Vector& v2)
      { return v1+v2; }
    };

    struct Construct_scaled_vector
    {
      Vector operator() (const Vector& v, 
			 typename Kernel::FT scale)
      { return scale*v; }
    };

    struct Construct_midpoint
    {
      Point operator() (const Point&p1, const Point& p2)
      { return typename Kernel::Midpoint_d()(p1, p2); }
    };

    // TODO Make the Construct_direction

    // Predicates
    struct Collinear
    {
      bool operator() (const Point&p1, const Point&p2, const Point&p3)
      { return ((p2-p1)*(p3-p2))==0; }
    };
  };

  /** Trait class for Linear_cell_complex class.
   *  2D version specialization.
   */
  template <class Kernel>
  struct Linear_cell_complex_traits<2,Kernel> : public Kernel
  {
    static const unsigned int ambient_dimension = 2;
    
    typedef typename Kernel::FT          FT;
    typedef typename Kernel::Point_2     Point;
    typedef typename Kernel::Vector_2    Vector;
    typedef typename Kernel::Direction_2 Direction;

    // Constructions
    typedef typename Kernel::Construct_translated_point_2
    Construct_translated_point;

    // TODO THE Construct_vector functor with two operators () (verify the kernel doc)
    struct Construct_vector : public Kernel::Construct_vector_2
    {
      using Kernel::Construct_vector_2::operator();      
      Vector operator() (typename Kernel::FT x1)
      {	return Kernel::Construct_vector_2()(x1, 0); }
    };

    typedef typename Kernel::Construct_sum_of_vectors_2
    Construct_sum_of_vectors;
    
    typedef typename Kernel::Construct_scaled_vector_2
    Construct_scaled_vector;

    typedef typename Kernel::Construct_midpoint_2
    Construct_midpoint;

    typedef typename Kernel::Construct_direction_2
    Construct_direction;

    /*    struct Vector_to_point TO REMOVE ?
    {
      Point operator() (const Vector&v)
      { return Kernel::Construct_translated_point(ORIGIN, v); }
      };*/
    
    // Predicates
    typedef typename Kernel::Collinear_2 Collinear;
  };

  /** Trait class for Linear_cell_complex class.
   *  3D version specialization.
   */
  template <class Kernel>
  struct Linear_cell_complex_traits<3,Kernel> : public Kernel
  {
    static const unsigned int ambient_dimension = 3;

    typedef typename Kernel::FT          FT;
    typedef typename Kernel::Point_3     Point;
    typedef typename Kernel::Vector_3    Vector;
    typedef typename Kernel::Direction_3 Direction;
    
    // Constructions
    typedef typename Kernel::Construct_translated_point_3 
    Construct_translated_point;

    // TODO the Construct_vector
    struct Construct_vector : public Kernel::Construct_vector_3
    {
      using Kernel::Construct_vector_3::operator();      
      Vector operator() (typename Kernel::FT x1)
      {	return Kernel::Construct_vector_3()(x1, 0, 0); }
      Vector operator() (typename Kernel::FT x1, typename Kernel::FT x2)
      {	return Kernel::Construct_vector_3()(x1, x2, 0); }
    };

    typedef typename Kernel::Construct_sum_of_vectors_3
    Construct_sum_of_vectors;

    typedef typename Kernel::Construct_scaled_vector_3 
    Construct_scaled_vector;
    
    typedef typename Kernel::Construct_midpoint_3
    Construct_midpoint;
    
    typedef typename Kernel::Construct_direction_3
    Construct_direction;

    typedef typename Kernel::Construct_normal_3
    Construct_normal;

    // Predicates
    typedef typename Kernel::Collinear_3 Collinear;
  };

} // namespace CGAL

#endif // CGAL_LINEAR_CELL_COMPLEX_TRAITS_H //
// EOF //
