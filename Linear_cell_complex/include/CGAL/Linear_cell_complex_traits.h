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
    typedef typename Kernel::FT       FT;
    typedef typename Kernel::Point_d  Point;
    typedef typename Kernel::Vector_d Vector;
    
    struct Collinear
    {
      bool operator() (const Point&p1, const Point&p2, const Point&p3)
      { return ((p2-p1)*(p3-p2))==0; }
    };

    struct Construct_translated_point
    {
      Point operator() (const Point&p, const Vector& v)
      { return p+v; }
    };
    struct Construct_midpoint
    {
      Point operator() (const Point&p1, const Point& p2)
      { return typename Kernel::Midpoint_d()(p1, p2); }
    };

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
    typedef typename Kernel::Vector_to_point_d 
    Vector_to_point;      
    struct Construct_scaled_vector
    {
      Vector operator() (const Vector& v, 
			 typename Kernel::FT scale)
      { return scale*v; }
    };
    struct Construct_sum_of_vectors
    {
      Vector operator() (const Vector&v1, const Vector& v2)
      { return v1+v2; }
    };
    struct Iso_rectangle
    {
      Iso_rectangle(const Point&p1, const Point& p2)
      {
	Point pmin,pmax;
	if ( compare_lexicographically(p1,p2)==SMALLER )
	  { pmin=p1; pmax=p2; }
	else
	  { pmin=p2; pmax=p1; }

	Vector v[2];
	unsigned int d=0;
	
 	for (unsigned int i=0; i<d_; ++i)
	  {
	    if ( p1[i]!=p2[i] )
	      {
		CGAL_assertion(d<2);
		v[d]=Vector(d_,typename Vector::Base_vector(),i);
		v[d] *=(pmax[i]-pmin[i]);
		++d;
	      }
	  }
	CGAL_assertion(d==2);
	
	p[0]=pmin;
	p[1]=Construct_translated_point()(pmin,v[0]);
	p[2]=pmax;
	p[3]=Construct_translated_point()(pmin,v[1]);
      }
      Iso_rectangle(const Iso_rectangle& air)
      {
	for (unsigned int i=0; i<4; ++i)
	  p[i]=air.p[i];	
      }
      Iso_rectangle& operator=(const Iso_rectangle& air) const
      {
	if ( this!=*air )
	  {
	    for (unsigned int i=0; i<4; ++i)
	      p[i]=air.p[i];	
	  }
	return *this;
      }
      Point& operator[] (unsigned int i)
      {
	CGAL_assertion(i<4);
	return p[i];
      }
      const Point& operator[] (unsigned int i) const
      {
	CGAL_assertion(i<4);
	return p[i];
      }
    private:
      Point p[4];
    };
    struct Iso_cuboid
    {
      Iso_cuboid(const Point&p1, const Point& p2)
      {
	Point pmin,pmax;
	if ( compare_lexicographically(p1,p2)==SMALLER )
	  { pmin=p1; pmax=p2; }
	else
	  { pmin=p2; pmax=p1; }

	Vector v[3];
	unsigned int d=0;
	
	for (unsigned int i=0; i<d_; ++i)
	  {
	    if ( p1[i]!=p2[i] )
	      {
		CGAL_assertion(d<3);
		v[d]=Vector(d_,typename Vector::Base_vector(),i);
		v[d] *=(pmax[i]-pmin[i]);
		++d;
	      }
	  }
	CGAL_assertion(d==3);
	
	p[0]=pmin;
	p[7]=pmax;

	p[1]=Construct_translated_point()(pmin,v[0]);
	p[2]=Construct_translated_point()(pmin,v[0]+v[1]);
	p[3]=Construct_translated_point()(pmin,v[1]);
	p[4]=Construct_translated_point()(pmin,v[1]+v[2]);
	p[5]=Construct_translated_point()(pmin,v[2]);
	p[6]=Construct_translated_point()(pmin,v[0]+v[2]);
      }
      Iso_cuboid(const Iso_cuboid& aic)
      {
	for (unsigned int i=0; i<8; ++i)
	  p[i]=aic.p[i];	
      }
      Iso_cuboid& operator=(const Iso_cuboid& aic) const
      {
	if ( this!=*aic )
	  {
	    for (unsigned int i=0; i<8; ++i)
	      p[i]=aic.p[i];	
	  }
	return *this;
      }
      Point& operator[] (unsigned int i)
      {
	CGAL_assertion(i<8);
	return p[i];
      }
      const Point& operator[] (unsigned int i) const
      {
	CGAL_assertion(i<8);
	return p[i];
      }
    private:
      Point p[8];
    };

    struct Construct_iso_cuboid
    {
      Iso_cuboid operator() (const Point&p1, const Point& p2)
      { return Iso_cuboid(p1,p2); }
    };
  };

  /** Trait class for Linear_cell_complex class.
   *  2D version specialization.
   */
  template <class Kernel>
  struct Linear_cell_complex_traits<2,Kernel> : public Kernel
  {
    typedef typename Kernel::FT       FT;
    typedef typename Kernel::Point_2  Point;
    typedef typename Kernel::Vector_2 Vector;

    typedef typename Kernel::Collinear_2 Collinear;

    typedef typename Kernel::Construct_translated_point_2
    Construct_translated_point;
    typedef typename Kernel::Construct_midpoint_2
    Construct_midpoint;

    struct Vector_to_point
    {
      Point operator() (const Vector&v)
      { return Kernel::Construct_translated_point(ORIGIN, v); }
    };
    struct Construct_vector : public Kernel::Construct_vector_2
    {
      using Kernel::Construct_vector_2::operator();      
      Vector operator() (typename Kernel::FT x1)
      {	return Kernel::Construct_vector_2()(x1, 0); }
    };
    typedef typename Kernel::Construct_scaled_vector_2
    Construct_scaled_vector;
    typedef typename Kernel::Construct_sum_of_vectors_2
    Construct_sum_of_vectors;
    
    typedef typename Kernel::Construct_direction_2
    Construct_direction;
    typedef typename CGAL::Direction_2<Kernel>
    Direction;

    typedef typename Kernel::Iso_rectangle_2
    Iso_rectangle;
  };

  /** Trait class for Linear_cell_complex class.
   *  3D version specialization.
   */
  template <class Kernel>
  struct Linear_cell_complex_traits<3,Kernel> : public Kernel
  {
    typedef typename Kernel::FT       FT;
    typedef typename Kernel::Point_3  Point;
    typedef typename Kernel::Vector_3 Vector;

    typedef typename Kernel::Collinear_3 Collinear;
    
    typedef typename Kernel::Construct_translated_point_3 
    Construct_translated_point;
    typedef typename Kernel::Construct_midpoint_3
    Construct_midpoint;
    
    struct Vector_to_point
    {
      Point operator() (const Vector&v)
      { return typename Kernel::Construct_translated_point_3()(ORIGIN, v); }
    };
    struct Construct_vector : public Kernel::Construct_vector_3
    {
      using Kernel::Construct_vector_3::operator();      
      Vector operator() (typename Kernel::FT x1)
      {	return Kernel::Construct_vector_3()(x1, 0, 0); }
      Vector operator() (typename Kernel::FT x1, typename Kernel::FT x2)
      {	return Kernel::Construct_vector_3()(x1, x2, 0); }
    };
    typedef typename Kernel::Construct_scaled_vector_3 
    Construct_scaled_vector;
    typedef typename Kernel::Construct_sum_of_vectors_3
    Construct_sum_of_vectors;
    
    typedef typename Kernel::Construct_normal_3
    Construct_normal;

    typedef typename Kernel::Iso_cuboid_3
    Iso_cuboid;

    typedef typename Kernel::Construct_iso_cuboid_3
    Construct_iso_cuboid;

    struct Iso_rectangle
    {
      Iso_rectangle(const Point&p1, const Point& p2)
      {
	Iso_cuboid ic(p1,p2);
	
	p[0]=p1;
	p[1]=ic[1];
	p[2]=p2;
	p[3]=ic[3];
      }
      Iso_rectangle(const Iso_rectangle& air)
      {
	for (unsigned int i=0; i<4; ++i)
	  p[i]=air.p[i];	
      }
      Iso_rectangle& operator=(const Iso_rectangle& air) const
      {
	if ( this!=*air )
	  {
	    for (unsigned int i=0; i<4; ++i)
	      p[i]=air.p[i];	
	  }
	return *this;
      }
      Point& operator[] (unsigned int i)
      {
	CGAL_assertion(i<4);
	return p[i];
      }
      const Point& operator[] (unsigned int i) const
      {
	CGAL_assertion(i<4);
	return p[i];
      }
    private:
      Point p[4];
    };
  };

} // namespace CGAL

#endif // CGAL_LINEAR_CELL_COMPLEX_TRAITS_H //
// EOF //
