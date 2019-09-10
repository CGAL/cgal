// Copyright (c) 2005  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Abdelkrim Mebarki <Abdelkrim.Mebarki@sophia.inria.fr>

#ifndef CGAL_EULER_INTEGRATOR_2_H_
#define CGAL_EULER_INTEGRATOR_2_H_

#include <CGAL/license/Stream_lines_2.h>


#include <CGAL/basic.h>

namespace CGAL {

// The class Euler_integrator_2 is a model of the concept Integrator
template <class VectorField_2>
class Euler_integrator_2
{
public:
  typedef Euler_integrator_2<VectorField_2> self;
  typedef typename VectorField_2::FT FT;
  typedef typename VectorField_2::Point_2 Point_2;
  typedef typename VectorField_2::Vector_2 Vector_2;
  typedef typename VectorField_2::Vector_field_2 Vector_field_2;
protected:
  FT default_integration_step;
public:
  Euler_integrator_2();

  Euler_integrator_2(const FT & integration_step);

  inline Point_2 operator()(const Point_2 & p, const Vector_field_2 &
			    vector_field_2, const bool & index) const;

  inline Point_2 operator()(const Point_2 & p, const Vector_field_2 &
			    vector_field_2, const FT &
			    integration_step, const bool & index)
    const;

  inline Point_2 operator()(const Point_2 & p, const Vector_field_2 &
			    vector_field_2, const FT &
			    integration_step, Vector_2 v, const bool &
			    index) const;

  inline FT get_default_integration_step()
    {
      return default_integration_step;
    }
  // just for debugging
  inline FT distance(const Point_2 & p, const Point_2 & q)
    {
      return sqrt(((p.x() - q.x())*(p.x() - q.x()))+((p.y() - q.y())*(p.y() - q.y())));
    }
};

template <class Vector_field>
Euler_integrator_2<Vector_field>::Euler_integrator_2()
  : default_integration_step(1.0)
{}

// An additional parameter in the constructor to specify the default integration step
template <class Vector_field>
Euler_integrator_2<Vector_field>::Euler_integrator_2(const FT & integration_step)
  : default_integration_step(integration_step)
{}

template <class Vector_field>
inline typename Euler_integrator_2<Vector_field>::Point_2 
Euler_integrator_2<Vector_field>::operator()
  (const Point_2 & p, const Vector_field_2 & , const FT & integration_step, Vector_2 v, const bool & index) const
{
  if (!index)
    {
      Vector_2 v_t(v.x()*(-1),v.y()*(-1));
      v = v_t;
    }
  Vector_2 Euler_step(v.x()*integration_step,
                      v.y()*integration_step);
  return Point_2(p.x() + Euler_step.x(), p.y() + Euler_step.y());
}

template <class Vector_field>
inline typename Euler_integrator_2<Vector_field>::Point_2 
Euler_integrator_2<Vector_field>::operator()
  (const Point_2 & p, const Vector_field_2 & vector_field_2, const FT & integration_step, const bool & index) const
{
  Vector_2 v;
  v = vector_field_2.get_field(p).first;
  return  this->operator()(p, vector_field_2, integration_step, v, index);
}

template <class Vector_field>
inline typename Euler_integrator_2<Vector_field>::Point_2 
Euler_integrator_2<Vector_field>::operator()
  (const Point_2 & p, const Vector_field_2 & vector_field_2, const bool & index) const
{
  return this->operator()(p, vector_field_2, default_integration_step,index);
}

} //namespace CGAL

#endif
