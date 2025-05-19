// Copyright (c) 2005  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Abdelkrim Mebarki <Abdelkrim.Mebarki@sophia.inria.fr>

#ifndef CGAL_RUNGE_KUTTA_INTEGRATOR_2_H_
#define CGAL_RUNGE_KUTTA_INTEGRATOR_2_H_

#include <CGAL/license/Stream_lines_2.h>


#include <CGAL/basic.h>
#include <CGAL/Euler_integrator_2.h>

namespace CGAL {

// The class Runge_kutta_integrator_2 is a model of the concept Integrator
template <class VectorField_2>
class Runge_kutta_integrator_2
{
public:
  typedef Runge_kutta_integrator_2<VectorField_2> self;
  typedef typename VectorField_2::FT FT;
  typedef typename VectorField_2::Point_2 Point_2;
  typedef typename VectorField_2::Vector_2 Vector_2;
  typedef typename VectorField_2::Vector_field_2 Vector_field_2;
protected:
  typedef CGAL::Euler_integrator_2<VectorField_2> Euler_integrator_2;

  Euler_integrator_2 euler_integrator_2;

  FT default_integration_step;

public:
  Runge_kutta_integrator_2();

  Runge_kutta_integrator_2(FT integration_step);

  ~Runge_kutta_integrator_2();

  Point_2 operator()(const Point_2 & p, const Vector_field_2 & vector_field_2, const bool & index) const;

  Point_2 operator()(const Point_2 & p, const Vector_field_2 & vector_field_2, const FT & integration_step, const bool & index) const;

  Point_2 operator()(const Point_2 & p, const Vector_field_2 & vector_field_2, const FT & integration_step, Vector_2 v, const bool & index) const;

  inline FT get_default_integration_step()
    {
      return default_integration_step;
    }
  Euler_integrator_2 * get_euler_integrator_2()
    {
      return &euler_integrator_2;
    }
  inline FT distance(const Point_2 & p, const Point_2 & q)
    {
      return sqrt(((p.x() - q.x())*(p.x() - q.x()))+((p.y()
                                                      -
                                                      q.y())*(p.y() - q.y())));
    }
};

template <class VectorField_2>
Runge_kutta_integrator_2<VectorField_2>::
Runge_kutta_integrator_2()
  : euler_integrator_2(),
    default_integration_step(euler_integrator_2.get_default_integration_step())
{}

template <class VectorField_2>
Runge_kutta_integrator_2<VectorField_2>::
Runge_kutta_integrator_2(FT integration_step)
  : euler_integrator_2(integration_step),
    default_integration_step(integration_step)
{}

template <class VectorField_2>
Runge_kutta_integrator_2<VectorField_2>::
~Runge_kutta_integrator_2()
{}

template <class VectorField_2>
inline
typename Runge_kutta_integrator_2<VectorField_2>::
Point_2 Runge_kutta_integrator_2<VectorField_2>::operator()
  (const Point_2 & p, const Vector_field_2& vector_field_2, const FT & integration_step, Vector_2 v, const bool & index) const
{
  Point_2 p1 = euler_integrator_2(p, vector_field_2, 0.5*integration_step, v, index);
  if(!vector_field_2.is_in_domain(p1))
    return p1;
  v = vector_field_2.get_field(p1).first;
  Point_2 p2 = euler_integrator_2(p, vector_field_2, integration_step,v, index);
  return p2;
}

template <class VectorField_2>
inline
typename Runge_kutta_integrator_2<VectorField_2>::Point_2 Runge_kutta_integrator_2<VectorField_2>::operator()
(const Point_2 & p, const Vector_field_2& vector_field_2, const FT &
 integration_step, const bool & index) const
{
  Vector_2 v;
  v = vector_field_2.get_field(p).first;
  Runge_kutta_integrator_2<VectorField_2>
    runge_kutta_integrator_2(default_integration_step);
  return this->operator()(p, vector_field_2, integration_step, v, index);
}

template <class VectorField_2>
inline
typename Runge_kutta_integrator_2<VectorField_2>::Point_2 Runge_kutta_integrator_2<VectorField_2>::operator()
(const Point_2 & p, const Vector_field_2& vector_field_2, const bool & index) const
{
  Vector_2 v;
  v = vector_field_2.get_field(p).first;
  return this->operator()(p, vector_field_2, default_integration_step, v, index);
}

} //namespace CGAL

#endif
