// Copyright (c) 1997  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
//
// Author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>

#ifndef CGAL_POINT_GENERATORS_3_H
#define CGAL_POINT_GENERATORS_3_H 1
#include <CGAL/generators.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/number_type_basic.h>

namespace CGAL {

template < class P, class Creator = 
                   Creator_uniform_3<typename Kernel_traits<P>::Kernel::RT,P> >
class Random_points_in_sphere_3 : public Random_generator_base<P> {
    void generate_point();
public:
    typedef Random_points_in_sphere_3<P,Creator> This;
    Random_points_in_sphere_3( double r = 1, Random& rnd = default_random)
        // g is an input iterator creating points of type `P' uniformly
        // distributed in the open sphere with radius r, i.e. |`*g'| < r .
        // Three random numbers are needed from `rnd' for each point
    : Random_generator_base<P>( r, rnd) { generate_point(); }
    This& operator++()    {
        generate_point();
        return *this;
    }
    This  operator++(int) {
        This tmp = *this;
        ++(*this);
        return tmp;
    }
};

template < class P, class Creator >
void
Random_points_in_sphere_3<P,Creator>::
generate_point() {
  // A strip between z and z+dz has an area independant of z
    typedef typename Creator::argument_type T;
    double alpha = this->_rnd.get_double() * 2.0 * CGAL_PI;
    double z     = 2 * this->_rnd.get_double() - 1.0;
    double r     = std::sqrt( 1 - z * z);
    r *= std::pow( this->_rnd.get_double() , 1.0/3.0 );  
    Creator creator;
    this->d_item = creator( T(this->d_range * r * std::cos(alpha)),
                            T(this->d_range * r * std::sin(alpha)),
                            T(this->d_range * z));
}


template < class P, class Creator = 
                   Creator_uniform_3<typename Kernel_traits<P>::Kernel::RT,P> >
class Random_points_on_sphere_3 : public Random_generator_base<P> {
    void generate_point();
public:
    typedef Random_points_on_sphere_3<P,Creator> This;
    Random_points_on_sphere_3( double r = 1, Random& rnd = default_random)
        // g is an input iterator creating points of type `P' uniformly
        // distributed on the sphere with radius r, i.e. |`*g'| == r . A
        // two random numbers are needed from `rnd' for each point.
    : Random_generator_base<P>( r, rnd) { generate_point(); }
    This& operator++()    {
        generate_point();
        return *this;
    }
    This  operator++(int) {
        This tmp = *this;
        ++(*this);
        return tmp;
    }
};

template < class P, class Creator >
void
Random_points_on_sphere_3<P,Creator>::
generate_point() {
  // A strip between z and z+dz has an area independant of z
    typedef typename Creator::argument_type T;
    double alpha = this->_rnd.get_double() * 2.0 * CGAL_PI;
    double z     = 2 * this->_rnd.get_double() - 1.0;
    double r     = std::sqrt( 1 - z * z);
    Creator creator;
    this->d_item = creator( T(this->d_range * r * std::cos(alpha)),
                            T(this->d_range * r * std::sin(alpha)),
                            T(this->d_range * z));
}


template < class P, class Creator = 
                   Creator_uniform_3<typename Kernel_traits<P>::Kernel::RT,P> >
class Random_points_in_cube_3 : public Random_generator_base<P>{
    void generate_point();
public:
    typedef Random_points_in_cube_3<P,Creator> This;
    Random_points_in_cube_3( double a = 1, Random& rnd = default_random)
    : Random_generator_base<P>( a, rnd) { generate_point(); }
    This& operator++()    {
        generate_point();
        return *this;
    }
    This  operator++(int) {
        This tmp = *this;
        ++(*this);
        return tmp;
    }
};

template < class P, class Creator >
void
Random_points_in_cube_3<P,Creator>::
generate_point() {
    typedef typename Creator::argument_type T;
    Creator creator;
    this->d_item =
	     creator( T(this->d_range * ( 2 * this->_rnd.get_double() - 1.0)),
                      T(this->d_range * ( 2 * this->_rnd.get_double() - 1.0)),
                      T(this->d_range * ( 2 * this->_rnd.get_double() - 1.0)));
}


template <class OutputIterator, class Creator>
OutputIterator
points_on_cube_grid_3( double a, std::size_t n, 
                         OutputIterator o, Creator creator)
{
    if  (n == 0)
        return o;

    int m = int(std::ceil(
                  std::sqrt(std::sqrt(static_cast<double>(n)))));

    while (m*m*m < int(n)) m++;

    double base = -a;  // Left and bottom boundary.
    double step = 2*a/(m-1);
    int j = 0;
    int k = 0;
    double px = base;
    double py = base;
    double pz = base;
    *o++ = creator( px, py, pz);
    for (std::size_t i = 1; i < n; i++) {
        j++;
        if ( j == m) {
           k++;
           if ( k == m) {
              py = base;
              px = base;
              pz = pz + step;
              k = 0;
           }
           else {
              px = base;
              py = py + step;
           }
           j = 0;
        } else {
           px = px + step;
        }
        *o++ = creator( px, py, pz);
    }
    return o;
}

template <class OutputIterator>
OutputIterator
points_on_cube_grid_3( double a, std::size_t n, OutputIterator o)
{
    typedef std::iterator_traits<OutputIterator> ITraits;
    typedef typename ITraits::value_type         P;
    return points_on_square_grid_3(a, n, o, 
                 Creator_uniform_3<typename Kernel_traits<P>::Kernel::RT,P>());
}


} //namespace CGAL

#endif // CGAL_POINT_GENERATORS_3_H //
// EOF //
