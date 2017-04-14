// Copyright (c) 1999-2004,2006-2009,2014-2015   INRIA Sophia-Antipolis (France).
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
//
//
// Author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//                 Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//                 Andreas Fabri <Andreas.Fabri@sophia.inria.fr>
//                 Nico Kruithof <Nico.Kruithof@sophia.inria.fr>
//                 Manuel Caroli <Manuel.Caroli@sophia.inria.fr>
//                 Aymeric Pellé <Aymeric.Pelle@sophia.inria.fr>
#ifndef CGAL_REGULAR_TRAITS_ADAPTOR_H
#define CGAL_REGULAR_TRAITS_ADAPTOR_H

#include <CGAL/license/Triangulation_3.h>

#include <boost/mpl/if.hpp>

namespace CGAL {

  template < class RTT, class ConstructPoint, class Functor_>
class Regular_traits_adaptor
{
  const ConstructPoint& cp;
  const Functor_& f;

  typedef RTT                                              RTraits;
  typedef Functor_                                         Functor;

  typedef typename RTraits::FT                             FT;

  typedef typename RTraits::Tetrahedron_3                  Tetrahedron_3;
  typedef typename RTraits::Plane_3                        Plane_3;
  typedef typename RTraits::Sphere_3                       Sphere_3;

  typedef typename RTT::Point_3                            Point_3;
  typedef typename RTT::Weighted_point_3                   Weighted_point_3;

  template <class T>
  struct Conv_wp_to_p
  {
    typedef typename boost::remove_reference<T>::type T_no_ref;

    typedef typename boost::mpl::if_<
      typename boost::is_const<T_no_ref>::type,
      typename boost::mpl::if_<typename boost::is_reference<T>::type,
                          const Point_3&,
                          const Point_3>::type,
      typename boost::mpl::if_<typename boost::is_reference<T>::type,
                          Point_3&,
                          Point_3>::type >::type type;
  };

  template <typename> struct result {};

  template <typename F, typename A0> struct result<F(A0)> {
    typedef typename Conv_wp_to_p<A0>::type A0p;
    typedef typename cpp11::result_of<F(A0p)>::type type;
  };

  template <typename F, typename A0, typename A1> struct result<F(A0,A1)> {
    typedef typename Conv_wp_to_p<A0>::type A0p;
    typedef typename Conv_wp_to_p<A1>::type A1p;
    typedef typename cpp11::result_of<F(A0p, A1p)>::type type;
  };

  template <typename F, typename A0, typename A1, typename A2> struct result<F(A0,A1,A2)> {
    typedef typename Conv_wp_to_p<A0>::type A0p;
    typedef typename Conv_wp_to_p<A1>::type A1p;
    typedef typename Conv_wp_to_p<A2>::type A2p;
    typedef typename cpp11::result_of<F(A0p, A1p, A2p)>::type type;
  };

  template <typename F, typename A0, typename A1, typename A2, typename A3>
  struct result<F(A0,A1,A2,A3)> {
    typedef typename Conv_wp_to_p<A0>::type A0p;
    typedef typename Conv_wp_to_p<A1>::type A1p;
    typedef typename Conv_wp_to_p<A2>::type A2p;
    typedef typename Conv_wp_to_p<A3>::type A3p;
    typedef typename cpp11::result_of<F(A0p, A1p, A2p, A3p)>::type type;
  };

public:
  Regular_traits_adaptor (const ConstructPoint& cp, const Functor& f)
    : cp(cp), f(f)
  { }
  


  typename cpp11::result_of<Functor(Tetrahedron_3)>::type operator() (const Tetrahedron_3& t) const
  {
    return f(t);
  }

 typename cpp11::result_of< Functor(Point_3,Point_3) >::type operator() (const Point_3& p0, const Point_3& p1) const
  {
    return f(p0, p1);
  }

  typename cpp11::result_of< Functor(Point_3,Point_3, Point_3) >::type operator() (const Point_3& p0, const Point_3& p1, const Point_3& p2) const
  {
    return f(p0, p1, p2);
  }

  typename cpp11::result_of< Functor(Point_3,Point_3, Point_3,Point_3) >::type operator() (const Point_3& p0, const Point_3& p1, const Point_3& p2, const Point_3& p3) const
  {
    return f(p0, p1, p2, p3);
  }

 typename cpp11::result_of< Functor(Point_3,Origin) >::type operator() (const Point_3& p0, const Origin& o) const
  {
    return f(p0, o);
  }

  typename cpp11::result_of< Functor(Point_3,Origin) >::type operator() (const Origin& o, const Point_3& p0) const
  {
    return f(o, p0);
  }

  typename cpp11::result_of< Functor(Plane_3,Point_3) >::type operator() (const Plane_3& pl, const Point_3& p) const
  {
    return f(pl, p);
  }

   typename cpp11::result_of< Functor(Point_3,Point_3) >::type  operator() (const Weighted_point_3& p0, const Weighted_point_3& p1) const
  {
    return f(cp(p0), cp(p1));
  }

  typename cpp11::result_of<Functor(Plane_3,Point_3)>::type operator() (const Plane_3& p0, const Weighted_point_3& p1) const
  {
    return f(p0, cp(p1));
  }

  typename cpp11::result_of<Functor(Point_3,Point_3,Point_3)>::type operator() (const Weighted_point_3& p0, const Weighted_point_3& p1,
                          const Weighted_point_3& p2) const
  {
    return f(cp(p0), cp(p1), cp(p2));
  }

  typename cpp11::result_of<Functor(Point_3,Point_3,Point_3,Point_3)>::type operator() (const Weighted_point_3& p0, const Weighted_point_3& p1,
                          const Weighted_point_3& p2, const Weighted_point_3& p3) const
  {
    return f(cp(p0), cp(p1), cp(p2), cp(p3));
  }

  typename cpp11::result_of<Functor(Point_3,Point_3,Point_3,FT)>::type operator() (const Weighted_point_3& p0, const Weighted_point_3& p1,
                          const Weighted_point_3& p2, const Weighted_point_3& p3,
                          const FT w) const
  {
    return f(cp(p0), cp(p1), cp(p2), cp(p3), w);
  }
  
};
 
}  // namespace CGAL

#endif /* CGAL_REGULAR_TRAITS_ADAPTOR_H */
