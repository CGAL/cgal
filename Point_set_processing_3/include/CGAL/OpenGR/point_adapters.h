// Copyright (c) 2018  GeometryFactory(France).
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
// Author(s) : Necip Fazil Yildiran

#ifndef CGAL_OPENGR_POINT_ADAPTERS_H
#define CGAL_OPENGR_POINT_ADAPTERS_H

#include <CGAL/Simple_cartesian.h>

#include <boost/tuple/tuple.hpp>
//typedef CGAL::Simple_cartesian<double> K;

namespace CGAL {

namespace OpenGR {

namespace internal {

// typedef CGAL::Simple_cartesian<double> K;
template<typename Kernel, typename Range1, typename Range2, typename PointMap1, typename PointMap2, typename VectorMap1, typename VectorMap2>
struct PointAdapter {
  public:
    enum {Dim = 3};
    typedef typename Kernel::FT Scalar;
    typedef Eigen::Matrix<Scalar, Dim, 1> VectorType;
  
  private:
// TODO    Eigen::Map<const VectorType> m_pos, m_normal, m_color;
    VectorType col; // TODO =0
    
  public:
    PointAdapter(const boost::tuple<const typename Range1::value_type, const PointMap1*, const VectorMap1*> &t)
// TODO      : m_pos    (Eigen::Map<const VectorType >( get( *t.get<1>()/*pmap*/, t.get<0>()/*it*/ ).cartesian_begin() )), 
// TODO        m_normal (Eigen::Map<const VectorType >( get( *t.get<2>()/*vmap*/, t.get<0>()/*it*/ ).cartesian_begin() ))
    {
      col[0] = rand();
      col[1] = rand();
      col[2] = rand();
    }

// TODO: What if they are of different types?
//    PointAdapter(const boost::tuple<const typename Range2::value_type, const PointMap2*, const VectorMap2*> &t)
//      : m_pos    (Eigen::Map<const VectorType >( get( t.get<1>()/*pmap*/, t.get<0>()/*it*/ ).cartesian_begin() )), 
//        m_normal (Eigen::Map<const VectorType >( get( t.get<2>()/*vmap*/, t.get<0>()/*it*/ ).cartesian_begin() ))
//    { }

//    inline PointAdapter(const extlib1::PointType1& p)
//      : m_pos   (Eigen::Map<const VectorType >( p.pos )), 
//        m_normal(Eigen::Map<const VectorType >( p.n ))
//        /*, m_color (Eigen::Map<const VectorType >( p.color ))*/
//    { }
//    /*  get(point_map1, *it); */ 
//    // An instance of point_map1 is needed?
//    inline PointAdapter(const extlib1::PointType1& p)
//      : m_pos   (Eigen::Map<const VectorType >( p.pos )), 
//        m_normal(Eigen::Map<const VectorType >( p.n ))
//        /*, m_color (Eigen::Map<const VectorType >( p.color ))*/
//    { }

    inline /*const Eigen::Map<*/ const VectorType /*>*/& pos()    const { return col; /* TODO return m_pos;*/ }  
    inline /*const Eigen::Map<*/ const VectorType /*>*/& normal() const { return col; /* TODO return m_normal;*/ }
    inline /*const Eigen::Map<*/ const VectorType /*>*/& color()  const { return col; /* m_color;*/ }
    inline /*const Eigen::Map<*/ const VectorType /*>*/& rgb()    const { return color(); /*m_color; */}
    
};

}

} } // end of namespace CGAL::OpenGR

#endif // CGAL_OPENGR_POINT_ADAPTERS_H