
// ============================================================================
//
// Copyright (c) 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision:  $
// release_date  : $CGAL_Date:  $
//
// file          : include/CGAL/Straight_2_stream.h
// source        : straight_2.fw
// author(s)     : Oren Nechushtan
//
// coordinator   : MPI, Saarbruecken
//
// ============================================================================


#ifndef CGAL_STRAIGHT_2_STREAM_H
#define CGAL_STRAIGHT_2_STREAM_H



CGAL_BEGIN_NAMESPACE

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

CGAL_END_NAMESPACE



#endif // CGAL_STRAIGHT_2_STREAM_H
