// Copyright (c) 2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Ophir Setter           <ophir.setter@cs.tau.ac.il>
//

/*!\file CGAL/IO/Fig_stream_Curve_renderer_2.h
 * \brief
  * provides \c CGAL::Fig_stream interface for the curve renderer
 */

#ifndef CGAL_FIG_STREAM_CURVE_RENDERER_2_H
#define CGAL_FIG_STREAM_CURVE_RENDERER_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>


#include <CGAL/IO/Fig_stream.h>
#include <CGAL/Curved_kernel_via_analysis_2/Curve_renderer_facade.h>

namespace CGAL {

/*! \brief
 * outputs a curve arc to \c Fig_stream
 */
template <class CKvA_2, class RatKernel>
CGAL::Fig_stream<RatKernel>&
operator <<
(CGAL::Fig_stream<RatKernel>& ws, const internal::Arc_2< CKvA_2 >& arc)
{
    typedef typename RatKernel::FT            NT;
    typedef typename RatKernel::Point_2       Point_2;
    typedef typename RatKernel::Segment_2     Segment_2;

    typedef Curve_renderer_facade<CKvA_2> Facade;

    typedef std::pair< int, int > Coord_2;
    typedef std::vector< Coord_2 > Coord_vec_2;

    boost::optional < Coord_2 > p1, p2;
    std::list<Coord_vec_2> points;

    Bbox_2 bbox (CGAL::to_double(ws.bounding_rect().xmin()),
                 CGAL::to_double(ws.bounding_rect().ymin()),
                 CGAL::to_double(ws.bounding_rect().xmax()),
                 CGAL::to_double(ws.bounding_rect().ymax()));

    // The maximum resolution of the renderer is 2048.
    // FIG units are 1200 per inch so we take 10 per inch. This should be
    // enough.
    const int width = (std::min)(ws.width() / 12, 2048);
    const int height = (std::min)(ws.height() / 12, 2048);
    Facade::setup(bbox, width, height);

    Facade::instance().draw(arc, points, &p1, &p2);
    if(points.empty())
        return ws;


    // When we create the polyline, we need to scale it back to "regular"
    // coordinates as the Fig_stream also scales them.
    const NT x_min = ws.bounding_rect().xmin();
    const NT y_min = ws.bounding_rect().ymin();
    const NT x_scale = NT(width) / (ws.bounding_rect().xmax() - x_min);
    const NT y_scale = NT(height) / (ws.bounding_rect().ymax() - y_min);
//    const NT scale = (std::min)(x_scale, y_scale);


    typename std::list<Coord_vec_2>::const_iterator lit = points.begin();
    while(lit != points.end())
    {
        const Coord_vec_2& vec = *lit;
        typename Coord_vec_2::const_iterator vit = vec.begin();

        std::vector< Point_2 > polyline;
        for (; vit != vec.end(); ++vit)
        {
            polyline.push_back(Point_2(vit->first / x_scale + x_min,
                                       vit->second / y_scale + y_min));
        }
        lit++;
        ws.write_polyline(polyline.begin(), polyline.end());
    }

    return ws;
}

/*! \brief
 *  outputs a curve point to \c Fig_stream
 */
template <class CKvA_2, class RatKernel>
CGAL::Fig_stream<RatKernel>&
operator <<
(CGAL::Fig_stream<RatKernel>& ws, const internal::Point_2< CKvA_2 >& pt)
{
#if 0 // TODO check why no point is drawn!
    typedef typename RatKernel::Point_2       Point_2;
    typedef typename RatKernel::Segment_2     Segment_2;

    typedef Curve_renderer_facade<CKvA_2> Facade;

    Bbox_2 bbox(CGAL::to_double(ws.bounding_rect().xmin()),
                CGAL::to_double(ws.bounding_rect().ymin()),
                CGAL::to_double(ws.bounding_rect().xmax()),
                CGAL::to_double(ws.bounding_rect().ymax()));

    Facade::setup(bbox, ws.width(), ws.height());

    std::pair< int, int > coord;
    if (!Facade::instance().draw(pt, coord)) {
        return ws;
    }

    int height = ws.height();

    Point_2 rat_pt(coord.first, height - coord.second);

    ws.set_fill_style(CGAL::FIG_FILLED);
    ws.set_fill_color(CGAL::FIG_BLACK);
    ws.set_point_style(CGAL::FIG_DISC);

    ws.write_point(rat_pt);
#endif
    return ws;
}

} //namespace CGAL

#endif // CGAL_FIG_STREAM_CURVE_RENDERER_2_H



