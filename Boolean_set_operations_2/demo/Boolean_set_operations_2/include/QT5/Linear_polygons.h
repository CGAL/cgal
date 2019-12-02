// Copyright (c) 2012  Tel-Aviv University (Israel).
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
// Author(s) : Apurva Bhatt <response2apurva@gmail.com>
//             Ronnie Gandhi <ronniegandhi19999@gmail.com>
//             Efi Fogel <efifogel@gmain.com>

#ifndef CGAL_QT_LINEAR_POLYGONS_H
#define CGAL_QT_LINEAR_POLYGONS_H

#include <CGAL/Qt/Converter.h>

#include "Typedefs.h"
#include "QT5/Piecewise_set_graphics_item.h"
#include "QT5/Boundary_pieces_graphics_item.h"

namespace CGAL {
namespace Qt {

struct Linear_X_monotone_bbox {
  template <typename X_monotone_linear_segment_2>
  CGAL::Bbox_2 operator()(X_monotone_linear_segment_2 const& aC) const
  { return aC.bbox(); }
};

struct Linear_bbox {
  template<class Linear_segment_2>
  CGAL::Bbox_2 operator()(Linear_segment_2 const& aC) const
  {
    double x_min = to_double(aC.source().x());
    double x_max = to_double(aC.target().x());
    double y_min = to_double(aC.source().y());
    double y_max = to_double(aC.target().y());
    if(x_min > x_max) {
      std::swap(x_min, x_max);
      std::swap(y_min, y_max);
    }

    if(y_min > y_max) std::swap(y_min, y_max);

    return Bbox_2(x_min, y_min, x_max, y_max);
  }
};

struct Draw_linear_X_monotone_curve {
  template <typename X_monotone_linear_segment_2, typename Path>
  void operator()(X_monotone_linear_segment_2 const& curve, Path& aPath,
                  int aIdx) const
  {
    //typedef Simple_cartesian<double> Linear_kernel;

    //commenting it gives errors
    typedef Point_2<Kernel>             Linear_point;
    typedef CGAL::Qt::Converter<Kernel> Converter;
    Converter convert;

    Linear_point ps(CGAL::to_double(curve.source().x()),
                    CGAL::to_double(curve.source().y()));
    Linear_point pt(CGAL::to_double(curve.target().x()),
                    CGAL::to_double(curve.target().y()));

    if(aIdx == 0) aPath.moveTo(convert(ps));
    else aPath.lineTo(convert(ps));
    aPath.lineTo(convert(pt));
  }
};

struct Draw_linear_curve {
  template <typename Linear_segment_2, class Path>
  void operator()(Linear_segment_2 const& curve, Path& aPath, int aIdx) const
  {
    //typedef Simple_cartesian<double> Linear_kernel;

    //commenting it gives errors
    typedef Point_2<Kernel>             Linear_point;
    typedef Qt::Converter<Kernel>       Converter;
    Converter convert;

    Linear_point ps(CGAL::to_double(curve.source().x()),
                    CGAL::to_double(curve.source().y()));
    Linear_point pt(CGAL::to_double(curve.target().x()),
                    CGAL::to_double(curve.target().y()));

    if (aIdx == 0) aPath.moveTo(convert(ps));
    else aPath.lineTo(convert(ps));
    aPath.lineTo(convert(pt));
  }
};

template <typename Linear_boundary_pieces>
class Linear_boundary_pieces_graphics_item :
    public Boundary_pieces_graphics_item<Linear_boundary_pieces,
                                         Draw_linear_curve, Linear_bbox>
{
  typedef Boundary_pieces_graphics_item<Linear_boundary_pieces,
                                        Draw_linear_curve, Linear_bbox> Base;

public:
  Linear_boundary_pieces_graphics_item(Linear_boundary_pieces* aPieces) :
    Base(aPieces)
  {}
};
/*
template<class Linear_boundary>
class Linear_boundary_graphics_item : public Piecewise_boundary_graphics_item<Linear_boundary,Draw_linear_X_monotone_curve,Linear_X_monotone_bbox>
{
  typedef Piecewise_boundary_graphics_item<Linear_boundary,Draw_linear_X_monotone_curve,Linear_X_monotone_bbox> Base;

public :

  Linear_boundary_graphics_item(Linear_boundary* aBoundary) : Base(aBoundary) {}
};
*/
/*
template<class Linear_region>
class Linear_region_graphics_item : public Piecewise_region_graphics_item<Linear_region,Draw_linear_X_monotone_curve,Linear_X_monotone_bbox>
{

  typedef Piecewise_region_graphics_item<Linear_region,Draw_linear_X_monotone_curve,Linear_X_monotone_bbox> Base;

public:

  Linear_region_graphics_item(Linear_region* aRegion) : Base(aRegion) {}
};
*/
template<class Linear_set,class Gps_traits>
class Linear_set_graphics_item :
    public Piecewise_set_graphics_item<Linear_set, Gps_traits,
                                       Draw_linear_X_monotone_curve,
                                       Linear_X_monotone_bbox>
{

  typedef Piecewise_set_graphics_item<Linear_set, Gps_traits,
                                      Draw_linear_X_monotone_curve,
                                      Linear_X_monotone_bbox> Base;
public:
  Linear_set_graphics_item(Linear_set* aSet, Gps_traits Linear_traits) :
    Base(aSet,Linear_traits)
  {}
};

} // namespace Qt
} // namespace CGAL

#endif
