// Copyright (c) 2000  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Geert-Jan Giezeman

#include <CGAL/Segment_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/kernel_assertions.h>
#include <CGAL/number_utils.h>
#include <vector>

#include <CGAL/Line_2_Line_2_intersection.h>
#include <CGAL/Intersection_traits_2.h>

namespace CGAL {

namespace internal {

template <class K>
struct Pointlist_2_rec_ {
    Pointlist_2_rec_ *next;
    typename K::Point_2 point;
    Oriented_side side;
};

template <class K>
struct Pointlist_2_ {
    int size;
    Pointlist_2_rec_<K> *first;
    Pointlist_2_() ;
    ~Pointlist_2_() ;
};



template <class K>
Pointlist_2_<K>::Pointlist_2_()
{
    size = 0;
    first = 0;
}

template <class K>
Pointlist_2_<K>::~Pointlist_2_()
{
    Pointlist_2_rec_<K> *cur;
    for (int i=0; i<size; i++) {
        cur = first;
        first = cur->next;
        delete cur;
    }
}




template <class K>
void _init_list(Pointlist_2_<K> &list,
                const typename K::Triangle_2 &trian)
{
    // check on degeneracies of trian.
    if (!trian.is_degenerate()) {
        list.size = 3;
        list.first = 0;
        for (int i=0; i<3; i++) {
            Pointlist_2_rec_<K> *newrec =
                        new Pointlist_2_rec_<K>;
            newrec->next = list.first;
            list.first = newrec;
            newrec->point = trian[i];
        }
    } else {
        // _not_implemented();
        CGAL_kernel_assertion(false);
    }
}



template <class K>
void _cut_off(Pointlist_2_<K> &list,
                const typename K::Line_2 &cutter)
{
    int i;
    const int list_size=list.size;
    Pointlist_2_rec_<K> *cur, *last=0, *newrec;
    for (i=0, cur = list.first; i<list_size; i++, cur = cur->next) {
        cur->side = cutter.oriented_side(cur->point);
        last = cur;
    }
    for (cur = list.first, i=0; i<list_size; i++, cur = cur->next) {
        if ((cur->side == ON_POSITIVE_SIDE
             && last->side == ON_NEGATIVE_SIDE)
           || (cur->side == ON_NEGATIVE_SIDE
               && last->side == ON_POSITIVE_SIDE)) {
            // add a vertex after cur
            typename K::Line_2 l(cur->point, last->point);
            ++list.size;
            newrec = new Pointlist_2_rec_<K>;
            newrec->next = last->next;
            last->next = newrec;
            newrec->side = ON_ORIENTED_BOUNDARY;
            Line_2_Line_2_pair<K> linepair(&cutter,  &l);
            CGAL_kernel_assertion_code(typename Line_2_Line_2_pair<K>::Intersection_results isr =)
            linepair.intersection_type();
            CGAL_kernel_assertion(isr == Line_2_Line_2_pair<K>::POINT);
            newrec->point = linepair.intersection_point();
        }
        last = cur;
    }
    CGAL_kernel_assertion(list.size-list_size <= 2);
    Pointlist_2_rec_<K> **curpt;
    curpt = &list.first;
    while (*curpt != 0) {
        cur = *curpt;
        if (cur->side == ON_NEGATIVE_SIDE) {
            --list.size;
            *curpt = cur->next;
            delete cur;
        } else {
            curpt = &cur->next;
        }
    }
    if (list_size == 2 && list.size-list_size == 1) {
        --list.size;
        cur = list.first;
        if (cur->side == ON_ORIENTED_BOUNDARY) {
            list.first = cur->next;
            delete cur;
        } else {
            last = cur;
            cur = cur->next;
            last->next = cur->next;
            delete cur;
        }
    }
}


template <class K>
class Triangle_2_Triangle_2_pair {
public:
    enum Intersection_results {NO_INTERSECTION, POINT, SEGMENT, TRIANGLE, POLYGON};
    Triangle_2_Triangle_2_pair(typename K::Triangle_2 const *trian1,
                               typename K::Triangle_2 const *trian2)
	: _trian1(trian1), _trian2(trian2), _known(false) {}

    Intersection_results intersection_type() const;

    typename K::Point_2     intersection_point() const;
    typename K::Segment_2   intersection_segment() const;
    typename K::Triangle_2  intersection_triangle() const;
    bool                    intersection(/*Polygon_2<R> &result*/) const;
    int                     vertex_count() const;
    typename K::Point_2     vertex(int i) const;
protected:
    typename K::Triangle_2 const*   _trian1;
    typename K::Triangle_2 const *  _trian2;
    mutable bool                    _known;
    mutable Intersection_results    _result;
    mutable Pointlist_2_<K>    _pointlist;
};

template <class K>
typename Triangle_2_Triangle_2_pair<K>::Intersection_results
Triangle_2_Triangle_2_pair<K>::intersection_type() const
{
  typedef typename K::Line_2 Line_2;
    if (_known)
        return _result;
// The non const this pointer is used to cast away const.
    _known = true;
    if (!do_overlap(_trian1->bbox(), _trian2->bbox())) {
        _result = NO_INTERSECTION;
        return _result;
    }
    _init_list(_pointlist, *_trian1);
    if (_trian2->is_degenerate()) {
        // _not_implemented();
        CGAL_kernel_assertion(false);
    } else {
        Line_2 l(_trian2->vertex(0), _trian2->vertex(1));
        if (l.oriented_side(_trian2->vertex(2)) == ON_POSITIVE_SIDE) {
            // counterclockwise triangle
            _cut_off(_pointlist, l);
            l = Line_2(_trian2->vertex(1), _trian2->vertex(2));
            _cut_off(_pointlist, l);
            l = Line_2(_trian2->vertex(2), _trian2->vertex(0));
            _cut_off(_pointlist, l);
        } else {
            l = l.opposite();
            _cut_off(_pointlist, l);
            l = Line_2(_trian2->vertex(0), _trian2->vertex(2));
            _cut_off(_pointlist, l);
            l = Line_2(_trian2->vertex(2), _trian2->vertex(1));
            _cut_off(_pointlist, l);
        }
    }
    switch (_pointlist.size) {
    case 0:
        _result = NO_INTERSECTION;
        break;
    case 1:
        _result = POINT;
        break;
    case 2:
        _result = SEGMENT;
        break;
    case 3:
        _result = TRIANGLE;
        break;
    default:
        _result = POLYGON;
    }
    return _result;
}


template <class K>
bool
Triangle_2_Triangle_2_pair<K>::intersection(
        /* Polygon_2<R> &result */) const
{
    if (!_known)
        intersection_type();
    if (_result != TRIANGLE  &&  _result != POLYGON)
        return false;
    Pointlist_2_rec_<K> *cur;
    int i;
    for (i=0, cur = _pointlist.first;
         i<_pointlist.size;
         i++, cur = cur->next) {
      std::cout << to_double(cur->point.x()) << ' ';
      std::cout << to_double(cur->point.y()) << ' ';
    }
    std::cout << std::endl;
    return true;
}

template <class K>
int
Triangle_2_Triangle_2_pair<K>::vertex_count() const
{
    CGAL_kernel_assertion(_known);
    return _pointlist.size;
}

template <class K>
typename K::Point_2
Triangle_2_Triangle_2_pair<K>::vertex(int n) const
{
    CGAL_kernel_assertion(_known);
    CGAL_kernel_assertion(n >= 0 && n < _pointlist.size);
    Pointlist_2_rec_<K> *cur;
    int k;
    for (k=0, cur = _pointlist.first;
         k < n;
         k++, cur = cur->next) {
    }
    return cur->point;
}

template <class K>
typename K::Triangle_2
Triangle_2_Triangle_2_pair<K>::intersection_triangle() const
{
  typedef typename K::Triangle_2 Triangle_2;
    if (!_known)
        intersection_type();
    CGAL_kernel_assertion(_result == TRIANGLE);
    return Triangle_2(_pointlist.first->point,
		      _pointlist.first->next->point,
		      _pointlist.first->next->next->point);
}

template <class K>
typename K::Segment_2
Triangle_2_Triangle_2_pair<K>::intersection_segment() const
{
  typedef typename K::Segment_2 Segment_2;
    if (!_known)
        intersection_type();
    CGAL_kernel_assertion(_result == SEGMENT);
    return Segment_2(_pointlist.first->point,
		     _pointlist.first->next->point);
}

template <class K>
typename K::Point_2
Triangle_2_Triangle_2_pair<K>::intersection_point() const
{
    if (!_known)
        intersection_type();
    CGAL_kernel_assertion(_result == POINT);
    return _pointlist.first->point;
}



template <class K>
typename CGAL::Intersection_traits
<K, typename K::Triangle_2, typename K::Triangle_2>::result_type
intersection(const typename K::Triangle_2 &tr1, 
	     const typename K::Triangle_2 &tr2,
	     const K&)
{
    typedef Triangle_2_Triangle_2_pair<K> is_t;
    is_t ispair(&tr1, &tr2);
    switch (ispair.intersection_type()) {
    case is_t::NO_INTERSECTION:
    default:
        return intersection_return<typename K::Intersect_2, typename K::Triangle_2, typename K::Triangle_2>();
    case is_t::POINT:
        return intersection_return<typename K::Intersect_2, typename K::Triangle_2, typename K::Triangle_2>(ispair.intersection_point());
    case is_t::SEGMENT:
        return intersection_return<typename K::Intersect_2, typename K::Triangle_2, typename K::Triangle_2>(ispair.intersection_segment());
    case is_t::TRIANGLE:
        return intersection_return<typename K::Intersect_2, typename K::Triangle_2, typename K::Triangle_2>(ispair.intersection_triangle());
    case is_t::POLYGON: {
        typedef std::vector<typename K::Point_2> Container;
        Container points(ispair.vertex_count());
        for (int i =0; i < ispair.vertex_count(); i++) {
            points[i] = ispair.vertex(i);
        }
        return intersection_return<typename K::Intersect_2, typename K::Triangle_2, typename K::Triangle_2>(points);
    }
    }
}

} // namespace internal

CGAL_INTERSECTION_FUNCTION_SELF(Triangle_2, 2)

} //namespace CGAL
