
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
// file          : include/CGAL/Triangle_2_Triangle_2_intersection.C
// source        : intersection_2_2.fw
// author(s)     : Geert-Jan Giezeman
//
// coordinator   : Saarbruecken
//
// ============================================================================



#include <CGAL/Segment_2.h>
#include <CGAL/Triangle_2.h>

CGAL_BEGIN_NAMESPACE

template <class R>
struct Pointlist_2_rec_ {
    Pointlist_2_rec_ *next;
    Point_2<R> point;
    Oriented_side side;
};

template <class R>
struct Pointlist_2_ {
    int size;
    Pointlist_2_rec_<R> *first;
    Pointlist_2_() ;
    ~Pointlist_2_() ;
};

template <class R>
class Triangle_2_Triangle_2_pair {
public:
    enum Intersection_results {NO, POINT, SEGMENT, TRIANGLE, POLYGON};
                        Triangle_2_Triangle_2_pair() ;
                        Triangle_2_Triangle_2_pair(
                                Triangle_2<R> const *trian1,
                                Triangle_2<R> const *trian2) ;
    ~Triangle_2_Triangle_2_pair() {}
#ifdef CGAL_CFG_RETURN_TYPE_BUG_2
    Intersection_results intersection_type() const
    {
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
            Line_2<R> l(_trian2->vertex(0), _trian2->vertex(1));
            if (l.oriented_side(_trian2->vertex(2)) == ON_POSITIVE_SIDE) {
                // counterclockwise triangle
                _cut_off(_pointlist, l);
                l = Line_2<R>(_trian2->vertex(1), _trian2->vertex(2));
                _cut_off(_pointlist, l);
                l = Line_2<R>(_trian2->vertex(2), _trian2->vertex(0));
                _cut_off(_pointlist, l);
            } else {
                l = l.opposite();
                _cut_off(_pointlist, l);
                l = Line_2<R>(_trian2->vertex(0), _trian2->vertex(2));
                _cut_off(_pointlist, l);
                l = Line_2<R>(_trian2->vertex(2), _trian2->vertex(1));
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
    
#else
    Intersection_results intersection_type() const;
#endif // CGAL_CFG_RETURN_TYPE_BUG_2
    bool                intersection(Point_2<R> &result) const;
    bool                intersection(Segment_2<R> &result) const;
    bool                intersection(Triangle_2<R> &result) const;
    bool                intersection(/*Polygon_2<R> &result*/) const;
    int                 vertex_count() const;
    Point_2<R>     vertex(int i) const;
protected:
    Triangle_2<R> const*   _trian1;
    Triangle_2<R> const *  _trian2;
    mutable bool                    _known;
    mutable Intersection_results    _result;
    mutable Pointlist_2_<R>    _pointlist;
};

template <class R>
inline bool do_intersect(
    const Triangle_2<R> &p1,
    const Triangle_2<R> &p2)
{
    typedef Triangle_2_Triangle_2_pair<R> pair_t;
    pair_t pair(&p1, &p2);
    return pair.intersection_type() != pair_t::NO;
}

CGAL_END_NAMESPACE



#include <CGAL/Line_2.h>
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>
#include <vector>

CGAL_BEGIN_NAMESPACE

template <class R>
Pointlist_2_<R>::Pointlist_2_()
{
    size = 0;
    first = 0;
}

template <class R>
Pointlist_2_<R>::~Pointlist_2_()
{
    Pointlist_2_rec_<R> *cur;
    for (int i=0; i<size; i++) {
        cur = first;
        first = cur->next;
        delete cur;
    }
}




template <class R>
void _init_list(Pointlist_2_<R> &list,
                const Triangle_2<R> &trian)
{
    // check on degeneracies of trian.
    if (!trian.is_degenerate()) {
        list.size = 3;
        list.first = 0;
        for (int i=0; i<3; i++) {
            Pointlist_2_rec_<R> *newrec =
                        new Pointlist_2_rec_<R>;
            newrec->next = list.first;
            list.first = newrec;
            newrec->point = trian[i];
        }
    } else {
        // _not_implemented();
        CGAL_kernel_assertion(false);
    }
}

CGAL_END_NAMESPACE

#include <CGAL/Line_2_Line_2_intersection.h>

CGAL_BEGIN_NAMESPACE

template <class R>
void _cut_off(Pointlist_2_<R> &list,
                const Line_2<R> &cutter)
{
    int i;
    int add = 0;
    Pointlist_2_rec_<R> *cur, *last=0, *newrec;
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
            Line_2<R> l(cur->point, last->point);
            newrec = new Pointlist_2_rec_<R>;
            newrec->next = last->next;
            last->next = newrec;
            newrec->side = ON_ORIENTED_BOUNDARY;
            Line_2_Line_2_pair<R> linepair(&cutter,  &l);
            Line_2_Line_2_pair<R>::Intersection_results isr;
            isr = linepair.intersection_type();
            CGAL_kernel_assertion(isr == Line_2_Line_2_pair<R>::POINT);
            linepair.intersection(newrec->point);
        }
        last = cur;
    }
    CGAL_kernel_assertion(add <= 2);
    Pointlist_2_rec_<R> **curpt;
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

template <class R>
Triangle_2_Triangle_2_pair<R>::
Triangle_2_Triangle_2_pair()
{
    _trian1 = 0;
    _trian2 = 0;
    _known = false;
}

template <class R>
Triangle_2_Triangle_2_pair<R>::
Triangle_2_Triangle_2_pair(Triangle_2<R> const *trian1,
         Triangle_2<R> const *trian2)
{
    _trian1 = trian1;
    _trian2 = trian2;
    _known = false;
}

#ifndef CGAL_CFG_RETURN_TYPE_BUG_2
template <class R>
Triangle_2_Triangle_2_pair<R>::Intersection_results
Triangle_2_Triangle_2_pair<R>::intersection_type() const
{
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
        Line_2<R> l(_trian2->vertex(0), _trian2->vertex(1));
        if (l.oriented_side(_trian2->vertex(2)) == ON_POSITIVE_SIDE) {
            // counterclockwise triangle
            _cut_off(_pointlist, l);
            l = Line_2<R>(_trian2->vertex(1), _trian2->vertex(2));
            _cut_off(_pointlist, l);
            l = Line_2<R>(_trian2->vertex(2), _trian2->vertex(0));
            _cut_off(_pointlist, l);
        } else {
            l = l.opposite();
            _cut_off(_pointlist, l);
            l = Line_2<R>(_trian2->vertex(0), _trian2->vertex(2));
            _cut_off(_pointlist, l);
            l = Line_2<R>(_trian2->vertex(2), _trian2->vertex(1));
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

#endif


template <class R>
bool
Triangle_2_Triangle_2_pair<R>::intersection(
        /* Polygon_2<R> &result */) const
{
    if (!_known)
        intersection_type();
    if (_result != TRIANGLE  &&  _result != POLYGON)
        return false;
    Pointlist_2_rec_<R> *cur;
    int i;
    for (i=0, cur = _pointlist.first;
         i<_pointlist.size;
         i++, cur = cur->next) {
      cout << to_double(cur->point.x()) << ' ';
      cout << to_double(cur->point.y()) << ' ';
    }
    cout << endl;
    return true;
}

template <class R>
int
Triangle_2_Triangle_2_pair<R>::vertex_count() const
{
    CGAL_kernel_assertion(_known);
    return _pointlist.size;
}

template <class R>
Point_2<R>
Triangle_2_Triangle_2_pair<R>::vertex(int n) const
{
    CGAL_kernel_assertion(_known);
    CGAL_kernel_assertion(n >= 0 && n < _pointlist.size);
    Pointlist_2_rec_<R> *cur;
    int k;
    for (k=0, cur = _pointlist.first;
         k < n;
         k++, cur = cur->next) {
    }
    return cur->point;
}

template <class R>
bool
Triangle_2_Triangle_2_pair<R>::intersection(
        Triangle_2<R> &result) const
{
    if (!_known)
        intersection_type();
    if (_result != TRIANGLE)
        return false;
    result = Triangle_2<R>(_pointlist.first->point,
                    _pointlist.first->next->point,
                    _pointlist.first->next->next->point);
    return true;
}

template <class R>
bool
Triangle_2_Triangle_2_pair<R>::intersection(
        Segment_2<R> &seg) const
{
    if (!_known)
        intersection_type();
    if (_result != SEGMENT)
        return false;
    seg = Segment_2<R>(_pointlist.first->point,
                    _pointlist.first->next->point);
    return true;
}

template <class R>
bool
Triangle_2_Triangle_2_pair<R>::intersection(
        Point_2<R> &pt) const
{
    if (!_known)
        intersection_type();
    if (_result != POINT)
        return false;
    pt = _pointlist.first->point;
    return true;
}

CGAL_END_NAMESPACE



CGAL_BEGIN_NAMESPACE

template <class R>
Object
intersection(const Triangle_2<R> &tr1, const Triangle_2<R>&tr2)
{
    typedef Triangle_2_Triangle_2_pair<R> is_t;
    is_t ispair(&tr1, &tr2);
    switch (ispair.intersection_type()) {
    case is_t::NO:
    default:
        return Object();
    case is_t::POINT: {
        Point_2<R> pt;
        ispair.intersection(pt);
        return Object(new Wrapper< Point_2<R> >(pt));
    }
    case is_t::SEGMENT: {
        Segment_2<R> iseg;
        ispair.intersection(iseg);
        return Object(new Wrapper< Segment_2<R> >(iseg));
    }
    case is_t::TRIANGLE: {
        Triangle_2<R> itr;
        ispair.intersection(itr);
        return Object(new Wrapper< Triangle_2<R> >(itr));
    }
    case is_t::POLYGON: {
        typedef CGAL_STD::vector<Point_2<R> > Container;
        Container points(ispair.vertex_count());
        for (int i =0; i < ispair.vertex_count(); i++) {
            points[i] = ispair.vertex(i);
        }
        return Object(new Wrapper< Container >(points));
    }
    }
}

CGAL_END_NAMESPACE

