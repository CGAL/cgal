 // Copyright (c) 2009  GeometryFactory (France)
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
//
// Author(s)     : Laurent Rineau


#ifndef CGAL_TRIANGULATION_2_FILTERED_PROJECTION_TRAITS_3_H
#define CGAL_TRIANGULATION_2_FILTERED_PROJECTION_TRAITS_3_H

#include <CGAL/Triangulation_2_projection_traits_3.h>
#include <CGAL/Filtered_predicate.h>

namespace CGAL {

template < class Filtered_kernel >
class Triangulation_2_filtered_projection_traits_3
  : public Triangulation_2_projection_traits_3<Filtered_kernel>
{
  typedef Filtered_kernel K;
  typedef Triangulation_2_filtered_projection_traits_3<K> Self;
  typedef Triangulation_2_projection_traits_3<K> Base;

  typedef typename K::Exact_kernel Exact_kernel;
  typedef typename K::Approximate_kernel Approximate_kernel;
  typedef typename K::C2E C2E;
  typedef typename K::C2F C2F;

public:
  typedef Triangulation_2_projection_traits_3<Exact_kernel> Exact_traits;
  typedef Triangulation_2_projection_traits_3<Approximate_kernel> Filtering_traits;

private:
  // private data
  Exact_traits exact_traits_;
  Filtering_traits filtering_traits_;
  // NOTE: The traits are precomputed, here. My small bench showed that it
  // is at least twice slower if the predicates are constructed from a
  // vector, instead of from a traits reference (with a precomputed
  // normal), because the Filtered_predicate converts the object that is
  // passed as argument its constructor.
public:
  Triangulation_2_filtered_projection_traits_3(const typename K::Vector_3& n)
    : Base(n),
      exact_traits_(C2E()(n)),
      filtering_traits_(C2F()(n))
  {
  }

  Self& operator=(const Self& other)
  {
    if(this != &other) {
      Base::operator=(other);
      exact_traits_ = other.exact_traits_;
      filtering_traits_ = other.filtering_traits_;
    }
    return *this;
  }

  const Exact_traits& exact_traits() const { return exact_traits_; }
  const Filtering_traits& filtering_traits() const { return filtering_traits_; }

  struct MyC2E : public C2E {
#ifndef CGAL_CFG_MATCHING_BUG_6
    using C2E::operator();
#else
    typedef typename C2E Converter;
    typedef typename Converter::Source_kernel Source_kernel;
    typedef typename Converter::Target_kernel Target_kernel;

    CGAL::Point_3<Target_kernel >
    operator()(const CGAL::Point_3<Kernel> & p) const
    {
      return Converter::operator()(p);
    }
#endif

    Exact_traits operator()(const Self& traits) const
    {
      return traits.exact_traits();
    }
  };

  struct MyC2F : public C2F {
#ifndef CGAL_CFG_MATCHING_BUG_6
    using C2F::operator();
#else
    typedef typename C2F Converter;
    typedef typename Converter::Source_kernel Source_kernel;
    typedef typename Converter::Target_kernel Target_kernel;

    CGAL::Point_3<Target_kernel >
    operator()(const CGAL::Point_3<Kernel> & p) const
    {
      return Converter::operator()(p);
    }
#endif

    Filtering_traits operator()(const Self& traits) const
    {
      return traits.filtering_traits();
    }
  }; // end class MyC2F

#define CGAL_TRIANGULATION_2_PROJ_TRAITS_FILTER_PRED(P, Pf, obj)    \
  typedef  Filtered_predicate< \
    typename Exact_traits::P, \
    typename Filtering_traits::P, \
    MyC2E, \
    MyC2F > P; \
  const P& Pf() const { return P(*this); }

  CGAL_TRIANGULATION_2_PROJ_TRAITS_FILTER_PRED(Orientation_2,
                                               orientation_2_object,
                                               orientation)
  CGAL_TRIANGULATION_2_PROJ_TRAITS_FILTER_PRED(Side_of_oriented_circle_2,
                                               side_of_oriented_circle_2_object,
                                               side_of_oriented_circle)
}; // end class Triangulation_2_projection_traits_3<Filtered_kernel>

} // end namespace CGAL


#endif // CGAL_TRIANGULATION_2_FILTERED_PROJECTION_TRAITS_3_H
