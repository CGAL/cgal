// Copyright (c) 2000
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Geert-Jan Giezeman


#ifndef CGAL_INTERSECTIONS_2_LINE_2_LINE_2_H
#define CGAL_INTERSECTIONS_2_LINE_2_LINE_2_H

#include <CGAL/config.h>
#include <CGAL/Line_2.h>
#include <CGAL/kernel_assertions.h>
#include <CGAL/number_utils.h>
#include <CGAL/Intersection_traits_2.h>

namespace CGAL {

class Cartesian_tag;
class Homogeneous_tag;

namespace Intersections {

namespace internal {

template <class K>
class Line_2_Line_2_pair {
public:
    enum Intersection_results {NOT_COMPUTED_YET, NO_INTERSECTION, POINT, LINE};
    Line_2_Line_2_pair(typename K::Line_2 const *line1,
                       typename K::Line_2 const *line2)
      : _line1(line1), _line2(line2), _result(NOT_COMPUTED_YET) {}

    Intersection_results intersection_type() const;

    typename K::Point_2 intersection_point() const;
    typename K::Line_2  intersection_line() const;
protected:
    typename K::Line_2 const*   _line1;
    typename K::Line_2 const *  _line2;
    mutable Intersection_results    _result;
    mutable typename K::Point_2         _intersection_point;
};

template <class K>
inline bool do_intersect(
    const typename K::Line_2 &p1,
    const typename K::Line_2 &p2,
    const K&)
{
    typedef Line_2_Line_2_pair<K> pair_t;
    pair_t pair(&p1, &p2);
    return pair.intersection_type() != pair_t::NO_INTERSECTION;
}



template <class K>
typename CGAL::Intersection_traits
<K, typename K::Line_2, typename K::Line_2>::result_type
intersection(const typename K::Line_2 &line1,
             const typename K::Line_2 &line2,
             const K&)
{
    typedef Line_2_Line_2_pair<K> is_t;
    is_t linepair(&line1, &line2);
    switch (linepair.intersection_type()) {
    case is_t::NO_INTERSECTION:
    default:
      return intersection_return<typename K::Intersect_2, typename K::Line_2, typename K::Line_2>();
    case is_t::POINT:
        return intersection_return<typename K::Intersect_2, typename K::Line_2, typename K::Line_2>(linepair.intersection_point());
    case is_t::LINE:
        return intersection_return<typename K::Intersect_2, typename K::Line_2, typename K::Line_2>(line1);
    }
}



template <class R, class POINT, class RT>
bool construct_if_finite(POINT &pt, RT x, RT y, RT w, R &, const Cartesian_tag &)
{
  typename R::Construct_point_2 construct_point;
    typedef typename R::FT FT;
    CGAL_kernel_precondition(CGAL_NTS is_finite(x)
                             && CGAL_NTS is_finite(y)
                             && w != RT(0));

    FT xw = FT(x)/FT(w);
    FT yw = FT(y)/FT(w);
    if (!CGAL_NTS is_finite(xw) || !CGAL_NTS is_finite(yw))
        return false;
    pt = construct_point(xw, yw);
    return true;
}



template <class R, class POINT, class RT>
bool construct_if_finite(POINT &pt, RT x, RT y, RT w, R &, const Homogeneous_tag&)
{
    typename R::Construct_point_2 construct_point;
    typedef typename R::FT FT;
    CGAL_kernel_precondition(CGAL_NTS is_finite(x)
                             && CGAL_NTS is_finite(y)
                             && w != RT(0));

    FT xw = FT(x)/FT(w);
    FT yw = FT(y)/FT(w);
    if (!CGAL_NTS is_finite(xw) || !CGAL_NTS is_finite(yw))
        return false;
    pt = construct_point(x, y, w);
    return true;
}


template <class R, class POINT, class RT>
inline
bool
construct_if_finite(POINT &pt, RT x, RT y, RT w, const R &r)
{
  typedef typename R::Kernel_tag Tag;
  Tag tag;
  return construct_if_finite(pt, x, y, w, r, tag);
}


template <class K>
typename Line_2_Line_2_pair<K>::Intersection_results
Line_2_Line_2_pair<K>::intersection_type() const
{
    typedef typename K::RT RT;
    if (_result != NOT_COMPUTED_YET)
        return _result;
    RT nom1, nom2, denom;
    // The non const this pointer is used to cast away const.
    denom = _line1->a()*_line2->b() - _line2->a()*_line1->b();
    if (denom == RT(0)) {
        if (RT(0) == (_line1->a()*_line2->c() - _line2->a()*_line1->c()) &&
            RT(0) == (_line1->b()*_line2->c() - _line2->b()*_line1->c()))
            _result = LINE;
        else
            _result = NO_INTERSECTION;
        return _result;
    }
    nom1 = (_line1->b()*_line2->c() - _line2->b()*_line1->c());
    if (!CGAL_NTS is_finite(nom1)) {
        _result = NO_INTERSECTION;
        return _result;
    }
    nom2 = (_line2->a()*_line1->c() - _line1->a()*_line2->c());
    if (!CGAL_NTS is_finite(nom2)) {
        _result = NO_INTERSECTION;
        return _result;
    }
    K dummyR;
    if (!construct_if_finite(_intersection_point,
                            nom1, nom2, denom, dummyR)){
        _result = NO_INTERSECTION;
        return _result;
    }
    _result = POINT;
    return _result;
}


template <class K>
typename K::Point_2
Line_2_Line_2_pair<K>::intersection_point() const
{
    if (_result == NOT_COMPUTED_YET)
        intersection_type();
    CGAL_kernel_assertion(_result == POINT);
    return _intersection_point;
}

template <class K>
typename K::Line_2
Line_2_Line_2_pair<K>::intersection_line() const
{
    if (_result == NOT_COMPUTED_YET)
        intersection_type();
    CGAL_kernel_assertion(_result == LINE);
    return *_line1;
}

} // namespace internal
} // namespace Intersections

CGAL_INTERSECTION_FUNCTION_SELF(Line_2, 2)
CGAL_DO_INTERSECT_FUNCTION_SELF(Line_2, 2)


} //namespace CGAL

#endif
