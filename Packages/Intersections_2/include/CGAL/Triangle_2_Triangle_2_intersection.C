
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



#include <CGAL/Segment_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/kernel_assertions.h>
#include <CGAL/number_utils.h>
#include <vector>

#include <CGAL/Line_2_Line_2_intersection.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

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
class Triangle_2_Triangle_2_pair {
public:
    enum Intersection_results {NO, POINT, SEGMENT, TRIANGLE, POLYGON};
                        Triangle_2_Triangle_2_pair() ;
                        Triangle_2_Triangle_2_pair(
                                typename K::Triangle_2 const *trian1,
                                typename K::Triangle_2 const *trian2) ;
    ~Triangle_2_Triangle_2_pair() {}

    Intersection_results intersection_type() const;

    bool                intersection(typename K::Point_2 &result) const;
    bool                intersection(typename K::Segment_2 &result) const;
    bool                intersection(typename K::Triangle_2 &result) const;
    bool                intersection(/*Polygon_2<R> &result*/) const;
    int                 vertex_count() const;
    typename K::Point_2     vertex(int i) const;
protected:
    typename K::Triangle_2 const*   _trian1;
    typename K::Triangle_2 const *  _trian2;
    mutable bool                    _known;
    mutable Intersection_results    _result;
    mutable Pointlist_2_<K>    _pointlist;
};


// template <class R>
// inline bool do_intersect(
//     const Triangle_2<R> &p1,
//     const Triangle_2<R> &p2)
// {
//     typedef Triangle_2_Triangle_2_pair<R> pair_t;
//     pair_t pair(&p1, &p2);
//     return pair.intersection_type() != pair_t::NO;
// }


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
    int add = 0;
    Pointlist_2_rec_<K> *cur, *last=0, *newrec;
    for (i=0, cur = list.first; i<list.size; i++, cur = cur->next) {
        cur->side = cutter.oriented_side(cur->point);
        last = cur;
    }
    for (cur = list.first, i=0; i<list.size; i++, cur = cur->next) {
        if ((cur->side == ON_POSITIVE_SIDE
             && last->side == ON_NEGATIVE_SIDE)
           || (cur->side == ON_NEGATIVE_SIDE
               && last->side == ON_POSITIVE_SIDE)) {
            // add a vertex after cur
            add++;
            typename K::Line_2 l(cur->point, last->point);
            newrec = new Pointlist_2_rec_<K>;
            newrec->next = last->next;
            last->next = newrec;
            newrec->side = ON_ORIENTED_BOUNDARY;
            Line_2_Line_2_pair<K> linepair(&cutter,  &l);
            typename Line_2_Line_2_pair<K>::Intersection_results isr;
            isr = linepair.intersection_type();
            CGAL_kernel_assertion(isr == Line_2_Line_2_pair<K>::POINT);
            linepair.intersection(newrec->point);
        }
        last = cur;
    }
    CGAL_kernel_assertion(add <= 2);
    Pointlist_2_rec_<K> **curpt;
    curpt = &list.first;
    while (*curpt != 0) {
        cur = *curpt;
        if (cur->side == ON_NEGATIVE_SIDE) {
            add--;
            *curpt = cur->next;
            delete cur;
        } else {
            curpt = &cur->next;
        }
    }
    if (list.size == 2 && add == 1) {
        add = 0;
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
    list.size += add;
}

template <class K>
Triangle_2_Triangle_2_pair<K>::
Triangle_2_Triangle_2_pair()
{
    _trian1 = 0;
    _trian2 = 0;
    _known = false;
}

template <class K>
Triangle_2_Triangle_2_pair<K>::
Triangle_2_Triangle_2_pair(typename K::Triangle_2 const *trian1,
			   typename K::Triangle_2 const *trian2)
{
    _trian1 = trian1;
    _trian2 = trian2;
    _known = false;
}

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
        _result = NO;
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
        _result = NO;
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
bool
Triangle_2_Triangle_2_pair<K>::intersection(
        typename K::Triangle_2 &result) const
{
  typedef typename K::Triangle_2 Triangle_2;
    if (!_known)
        intersection_type();
    if (_result != TRIANGLE)
        return false;
    result = Triangle_2(_pointlist.first->point,
			_pointlist.first->next->point,
			_pointlist.first->next->next->point);
    return true;
}

template <class K>
bool
Triangle_2_Triangle_2_pair<K>::intersection(
        typename K::Segment_2 &seg) const
{
  typedef typename K::Segment_2 Segment_2;
    if (!_known)
        intersection_type();
    if (_result != SEGMENT)
        return false;
    seg = Segment_2(_pointlist.first->point,
		    _pointlist.first->next->point);
    return true;
}

template <class K>
bool
Triangle_2_Triangle_2_pair<K>::intersection(
        typename K::Point_2 &pt) const
{
    if (!_known)
        intersection_type();
    if (_result != POINT)
        return false;
    pt = _pointlist.first->point;
    return true;
}



template <class K>
Object
intersection(const typename CGAL_WRAP(K)::Triangle_2 &tr1, 
	     const typename CGAL_WRAP(K)::Triangle_2 &tr2,
	     const K& k)
{
    typedef Triangle_2_Triangle_2_pair<K> is_t;
    is_t ispair(&tr1, &tr2);
    switch (ispair.intersection_type()) {
    case is_t::NO:
    default:
        return Object();
    case is_t::POINT: {
        typename K::Point_2 pt;
        ispair.intersection(pt);
        return make_object(pt);
    }
    case is_t::SEGMENT: {
        typename K::Segment_2 iseg;
        ispair.intersection(iseg);
        return make_object(iseg);
    }
    case is_t::TRIANGLE: {
        typename K::Triangle_2 itr;
        ispair.intersection(itr);
        return make_object(itr);
    }
    case is_t::POLYGON: {
        typedef CGAL_STD::vector<typename K::Point_2> Container;
        Container points(ispair.vertex_count());
        for (int i =0; i < ispair.vertex_count(); i++) {
            points[i] = ispair.vertex(i);
        }
        return make_object(points);
    }
    }
}

} // namespace CGALi




CGAL_END_NAMESPACE

