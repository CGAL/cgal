// Copyright (c) 2000  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
// Author(s)     : Oren Nechushtan

#ifndef CGAL_STRAIGHT_2_STREAM_H
#define CGAL_STRAIGHT_2_STREAM_H

#include <CGAL/Straight_2.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

#ifndef CGAL_NO_OSTREAM_INSERT_STRAIGHT_2
template <class R>
std::ostream& operator<<(std::ostream& os,const Straight_2_<R>& cv)
{
        typedef Pm_straight_exact_traits<R> Traits;
        typedef Straight_2_<R> Curve;
        switch(cv.current_state())
        {
        case Curve::SEGMENT:
                {
                        Segment_2<R> seg;
                        cv.current(seg);
                        return os << seg;
                }
        case Curve::RAY:
                {
                        Ray_2<R> ray;
                        cv.current(ray);
                        return os << ray;
                }
        case Curve::LINE:
                {
                        Line_2<R> line;
                        cv.current(line);
                        return os << line;
                }
        case Curve::POINT:
                {
                        Point_2<R> p;
                        cv.current(p);
                        return os << p;
                }
        case Curve::EMPTY:
          break;
        }
        CGAL_assertion_msg(
                cv.current_state()==Curve::SEGMENT||
                cv.current_state()==Curve::RAY||
                cv.current_state()==Curve::LINE||
                cv.current_state()==Curve::POINT||
                cv.current_state()==Curve::EMPTY,
                "\nUnknown type in  std:: ostream& operator<<( \
                std:: ostream& os,const Straight_2&)");
        return os;
}
#endif //CGAL_NO_OSTREAM_INSERT_STRAIGHT_2
#ifndef CGAL_NO_ISTREAM_EXTRACT_STRAIGHT_2
template <class R>
std:: istream& operator>>(std:: istream& is,Straight_2_<R>& cv)
{
        typedef Pm_straight_exact_traits<R> Traits;
        typedef Straight_2_<R> Curve;
        switch(cv.current_state())
        {
        case Curve::SEGMENT:
                {
                        Segment_2<R> seg;
                        cv.current(seg);
                        return os >> seg;
                }
        case Curve::RAY:
                {
                        Ray_2<R> ray;
                        cv.current(ray);
                        return os >> ray;
                }
        case Curve::LINE:
                {
                        Line_2<R> line;
                        cv.current(line);
                        return os >> line;
                }
        case Curve::POINT:
                {
                        Point_2<R> p;
                        cv.current(p);
                        return os >> p;
                }
        case Curve::EMPTY:
          break;
        }
        CGAL_assertion_msg(
                cv.current_state()==Curve::SEGMENT||
                cv.current_state()==Curve::RAY||
                cv.current_state()==Curve::LINE||
                cv.current_state()==Curve::POINT||
                cv.current_state()==Curve::EMPTY,
                "\nUnknown type in  std:: ostream& operator>>( \
                std:: ostream& os,Straight_2&)");
        return os;
}
#endif //CGAL_NO_ISTREAM_EXTRACT_STRAIGHT_2

} // CGALi

CGAL_END_NAMESPACE

#endif // CGAL_STRAIGHT_2_STREAM_H
