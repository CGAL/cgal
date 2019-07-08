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

#ifndef CGAL_OPENGR_POINT_ADAPTER_H
#define CGAL_OPENGR_POINT_ADAPTER_H

#include <CGAL/Simple_cartesian.h>

#include <boost/tuple/tuple.hpp>
#include <iostream>

namespace CGAL {

namespace OpenGR {

namespace internal {

template<typename Kernel>
struct PointAdapter {
  public:
    enum {Dim = 3};
    using Scalar = typename Kernel::FT;
    using VectorType = typename Eigen::Matrix<Scalar, Dim, 1>;
  
  private:
    using Base = typename gr::Point3D<Scalar>;
    Base base;
    
  public:
    PointAdapter(const PointAdapter&) = default;

    //template<typename RangeVal, typename PointMap, typename VectorMap>
    //PointAdapter(const boost::tuple<RangeVal, PointMap, VectorMap> &t)
    template<typename RangePointVectorTuple>
    inline PointAdapter(const RangePointVectorTuple &t)
          : base( get(t.get<1>(), t.get<0>()).x(),
                  get(t.get<1>(), t.get<0>()).y(),
                  get(t.get<1>(), t.get<0>()).z() ) // pos
    {
      using BaseVectorType = typename Base::VectorType;
      
      // normal
      base.set_normal(
        BaseVectorType( get(t.get<2>(), t.get<0>()).x(),
                        get(t.get<2>(), t.get<0>()).y(),
                        get(t.get<2>(), t.get<0>()).z() ) 
                     );
    }

    inline const VectorType  pos()    const { return base.pos();    }  
    inline const VectorType  normal() const { return base.normal(); }
    inline const VectorType  rgb()    const { return base.rgb();    }
    
};

}

} }

#endif // CGAL_OPENGR_POINT_ADAPTER_H