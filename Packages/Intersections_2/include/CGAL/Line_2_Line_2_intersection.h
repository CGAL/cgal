
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Geert-Jan Giezeman


#ifndef CGAL_LINE_2_LINE_2_INTERSECTION_H
#define CGAL_LINE_2_LINE_2_INTERSECTION_H

#include <CGAL/config.h>
#include <CGAL/Line_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>
#include <CGAL/Object.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

template <class K>
class Line_2_Line_2_pair {
public:
    enum Intersection_results {NO, POINT, LINE};
    Line_2_Line_2_pair() ;
    Line_2_Line_2_pair(typename K::Line_2 const *line1,
		       typename K::Line_2 const *line2);
    ~Line_2_Line_2_pair() {}

    Intersection_results intersection_type() const;

    bool                intersection(typename K::Point_2 &result) const;
    bool                intersection(typename K::Line_2 &result) const;
protected:
    typename K::Line_2 const*   _line1;
    typename K::Line_2 const *  _line2;
    mutable bool                    _known;
    mutable Intersection_results    _result;
    mutable typename K::Point_2         _intersection_point;
};

template <class K>
inline bool do_intersect(
    const typename CGAL_WRAP(K)::Line_2 &p1,
    const typename CGAL_WRAP(K)::Line_2 &p2,
    const K&)
{
    typedef Line_2_Line_2_pair<K> pair_t;
    pair_t pair(&p1, &p2);
    return pair.intersection_type() != pair_t::NO;
}



template <class K>
Object
intersection(const typename CGAL_WRAP(K)::Line_2 &line1, 
	     const typename CGAL_WRAP(K)::Line_2 &line2,
	     const K&)
{
    typedef Line_2_Line_2_pair<K> is_t;
    is_t linepair(&line1, &line2);
    switch (linepair.intersection_type()) {
    case is_t::NO:
    default:
        return Object();
    case is_t::POINT: {
        typename K::Point_2 pt;
        linepair.intersection(pt);
        return make_object(pt);
    }
    case is_t::LINE:
        return make_object(line1);
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
Line_2_Line_2_pair<K>::Line_2_Line_2_pair()
{
    _line1 = 0;
    _line2 = 0;
    _known = false;
}

template <class K>
Line_2_Line_2_pair<K>::Line_2_Line_2_pair(
    typename K::Line_2 const *line1, typename K::Line_2 const *line2)
{
    _line1 = line1;
    _line2 = line2;
    _known = false;
}

template <class K>
typename Line_2_Line_2_pair<K>::Intersection_results
Line_2_Line_2_pair<K>::intersection_type() const
{
    typedef typename K::RT RT;
    if (_known)
        return _result;
    RT nom1, nom2, denom;
    // The non const this pointer is used to cast away const.
    _known = true;
    denom = _line1->a()*_line2->b() - _line2->a()*_line1->b();
    if (denom == RT(0)) {
        if (RT(0) == (_line1->a()*_line2->c() - _line2->a()*_line1->c()) &&
            RT(0) == (_line1->b()*_line2->c() - _line2->b()*_line1->c()))
            _result = LINE;
        else
            _result = NO;
        return _result;
    }
    nom1 = (_line1->b()*_line2->c() - _line2->b()*_line1->c());
    if (!CGAL_NTS is_finite(nom1)) {
        _result = NO;
        return _result;
    }
    nom2 = (_line2->a()*_line1->c() - _line1->a()*_line2->c());
    if (!CGAL_NTS is_finite(nom2)) {
        _result = NO;
        return _result;
    }
    K dummyR;
    if (!construct_if_finite(_intersection_point,
                            nom1, nom2, denom, dummyR)){
        _result = NO;
        return _result;
    }
    _result = POINT;
    return _result;
}


template <class K>
bool
Line_2_Line_2_pair<K>::intersection(typename K::Point_2 &pt) const
{
    if (!_known)
        intersection_type();
    if (_result != POINT)
        return false;
    pt = _intersection_point;
    return true;
}

template <class K>
bool
Line_2_Line_2_pair<K>::intersection(typename K::Line_2 &l) const
{
    if (!_known)
        intersection_type();
    if (_result != LINE)
        return false;
    l = *_line1;
    return true;
}

} // namespace CGALi



template <class K>
inline bool do_intersect(
    const Line_2<K> &p1,
    const Line_2<K> &p2)
{
  typedef typename K::Do_intersect_2 Do_intersect;
  return Do_intersect()(p1, p2);
}

template <class K>
Object
intersection(const Line_2<K> &line1, const Line_2<K> &line2)
{
  typedef typename K::Intersect_2 Intersect;
  return Intersect()(line1, line2);
}


CGAL_END_NAMESPACE



#endif
