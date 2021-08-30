// Copyright (c) 2010  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Olivier Devillers



#ifndef CGAL_POINT_GENERATORS_D_H
#define CGAL_POINT_GENERATORS_D_H 1

#include <CGAL/disable_warnings.h>

#include <CGAL/generators.h>
#include <CGAL/number_type_basic.h>
#include <CGAL/number_type_config.h>

#include <cmath>
#include <vector>

namespace CGAL {


template < class P >
class Random_points_in_ball_d : public Random_generator_base<P>{
    void generate_point();
    int dimension;
public:
    typedef Random_points_in_ball_d<P> This;
 Random_points_in_ball_d( int dim, double a = 1,
                               Random& rnd = get_default_random())
        // g is an input iterator creating points of type `P' uniformly
        // distributed in the open sphere with radius r, i.e. |`*g'| < r .
   : Random_generator_base<P>( a, rnd), dimension(dim) { generate_point(); }
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

template < class P >
void
Random_points_in_ball_d<P>::
generate_point() {
  double norm = 0;
  std::vector< double > coord(dimension);

  for(int i=0; i<dimension; ++i) {
    // normal distribution
    //( a product of normal distib is a normal distrib in higher dim)
    double a=this->_rnd.get_double();
    a = std::sqrt( -2* std::log(1-a) );
    double b=this->_rnd.get_double();
    b = std::cos(2*CGAL_PI*b);
    coord[i]= a*b;
    norm += coord[i]*coord[i];
  }
  norm = this->d_range  * std::pow(this->_rnd.get_double(),1.0/dimension)
           /std::sqrt(norm);
  for(int i=0; i<dimension; ++i) coord[i] *= norm;
  this->d_item = P(dimension, coord.begin(), coord.end() );
}




template < class P >
class Random_points_on_sphere_d : public Random_generator_base<P>{
    void generate_point();
    int dimension;
public:
    typedef Random_points_on_sphere_d<P> This;
 Random_points_on_sphere_d( int dim, double a = 1,
                               Random& rnd = get_default_random())
        // g is an input iterator creating points of type `P' uniformly
        // distributed on the sphere with radius r, i.e. |`*g'| == r .
   : Random_generator_base<P>( a, rnd), dimension(dim) { generate_point(); }
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

template < class P >
void
Random_points_on_sphere_d<P>::
generate_point() {
  double norm = 0;
  std::vector< double > coord(dimension);

  for(int i=0; i<dimension; ++i) {
    // normal distribution
    double a=this->_rnd.get_double();
    a = std::sqrt( -2* std::log(1-a) );
    double b=this->_rnd.get_double();
    b = std::cos(2*CGAL_PI*b);
    coord[i]= a*b;
    norm += coord[i]*coord[i];
  }
  norm = this->d_range /std::sqrt(norm);
  for(int i=0; i<dimension; ++i) coord[i] *= norm;
  this->d_item = P(dimension, coord.begin(), coord.end() );
}


template < class P >
class Random_points_in_cube_d : public Random_generator_base<P>{
    void generate_point();
    int dimension;
public:
    typedef Random_points_in_cube_d<P> This;
 Random_points_in_cube_d( int dim, double a = 1,
                               Random& rnd = get_default_random())
   : Random_generator_base<P>( a, rnd), dimension(dim) { generate_point(); }
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

template < class P >
void
Random_points_in_cube_d<P>::
generate_point() {
    typedef typename Kernel_traits<P>::Kernel::RT RT;
    CGAL_assume(dimension>0);
    std::vector<RT> coord(dimension);
    for(int i=0; i<dimension; ++i)
      coord[i]=RT(this->d_range * ( 2 * this->_rnd.get_double() - 1.0));

    P p(dimension, coord.begin(), coord.end() );
    this->d_item = p;
}


template <class OutputIterator, class Creator>
OutputIterator
points_on_cube_grid_d( int dimension, double a,
                            std::size_t n, OutputIterator o, Creator creator)
{
  //  typedef typename OutputIterator::container_type::value_type Point;

    if  (n == 0)
        return o;

    // take m smallest such that m^dimension > n
    int m=int(std::floor(std::pow(static_cast<double>(n),1/(double)dimension)));
    while(true) {
      int nn=1;
      for (int i=0; i<dimension; ++i) nn*=m;
      if (std::size_t(nn)>=n) break;
      ++m;
    }

    double base = -a;  // Left and bottom boundary.
    double step = 2*a/(m-1);
    std::vector<int> indices(dimension);
    std::vector<double> coord(dimension);

    //initialize indices
    int j;
    for(j=0; j< dimension; ++j) { indices[j]=0; }

    std::size_t i=0;
    while (true) {
      //compute point
      for(j=0; j< dimension; ++j) {
        coord[j]=base+step*indices[j];
      }
      *o = creator(coord.begin(), coord.end() );
      ++i;
      if (i==n) break;
      ++o;
      // find index to increment
      for ( j=0; j < dimension; ++j) if ( indices[j]<m-1 ) break;
      // increment and reset smaller indices
      CGAL_assertion(j<dimension);
      indices[j]++;
      while(j>0) { --j; indices[j]=0; }
    }
    return o;
}


} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_POINT_GENERATORS_D_H //
// EOF //
