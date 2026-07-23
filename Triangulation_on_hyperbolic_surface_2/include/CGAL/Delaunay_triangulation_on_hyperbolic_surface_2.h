// Copyright (c) 2025
// INRIA Nancy (France), and Université de Lorraine (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Vincent Despré, Camille Lanuel, Marc Pouget, Monique Teillaud

#ifndef CGAL_DELAUNAY_TRIANGULATION_ON_HYPERBOLIC_SURFACE_2_H
#define CGAL_DELAUNAY_TRIANGULATION_ON_HYPERBOLIC_SURFACE_2_H

#include <CGAL/license/Triangulation_on_hyperbolic_surface_2.h>

#include <CGAL/Triangulation_on_hyperbolic_surface_2.h>
#include <CGAL/Root_of_traits.h>
#include <boost/numeric/interval.hpp>
#include <CGAL/Gmpq.h>
#include <CGAL/Algebraic_kernel_for_circles_2_2.h>
#include <CGAL/Circular_kernel_2.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_CK_traits_2.h>
#include <CGAL/Hyperbolic_surface_Delaunay_traits_2.h>

namespace CGAL{

template<class Traits>
struct Delaunay_triangulation_attributes {
    template<class CMap>
    struct Dart_wrapper{
        typedef Cell_attribute<CMap, Complex_number<typename Traits::FT>> Edge_attribute;
        typedef Cell_attribute<CMap, typename Triangulation_on_hyperbolic_surface_2<Traits, Delaunay_triangulation_attributes<Traits>>::Anchor> Face_attribute;
        typedef std::tuple<void,Edge_attribute,Face_attribute> Attributes;
    };
};

template<class Traits>
class Delaunay_triangulation_on_hyperbolic_surface_2: public Triangulation_on_hyperbolic_surface_2<Traits, Delaunay_triangulation_attributes<Traits>> {
public:
    typedef Triangulation_on_hyperbolic_surface_2<Traits,Delaunay_triangulation_attributes<Traits>>                     Base;
    typedef Combinatorial_map<2,Delaunay_triangulation_attributes<Traits>>                                              CMap;
    typedef typename Triangulation_on_hyperbolic_surface_2<Traits,Delaunay_triangulation_attributes<Traits>>::Anchor    Anchor;

    typedef typename Traits::FT                                         Number;
    typedef typename Traits::Complex                                    Complex_number;
    typedef typename Traits::Hyperbolic_point_2                         Point;
    typedef typename CMap::Dart_descriptor                              Dart_descriptor;
    typedef typename CMap::Dart_range                                   Dart_range;
    typedef typename CMap::template One_dart_per_cell_range<0>          Vertex_range;
    typedef typename CMap::template One_dart_per_cell_range<1>          Edge_range;
    typedef typename CMap::template One_dart_per_cell_range<2>          Face_range;
    typedef typename CMap::Dart_const_descriptor                        Dart_const_descriptor;
    typedef typename CMap::Dart_const_range                             Dart_const_range;
    typedef typename CMap::template One_dart_per_cell_const_range<1>    Edge_const_range;
    typedef typename CMap::template One_dart_per_cell_const_range<2>    Face_const_range;

    // undocumented types
    typedef typename Root_of_traits<Number>::Root_of_2                  Algebraic_number;
    typedef typename Traits::Hyperbolic_Voronoi_point_2                 Voronoi_point;
    typedef Hyperbolic_isometry_2<Traits>                               Isometry;
    typedef boost::numeric::interval<double>                            Interval;

    enum Locate_type {
        VERTEX = 0,
        EDGE,
        FACE,
        OUTSIDE
    };
    enum Locate_walk {
      STRAIGHT = 0,
      VISIBILITY,
    };

    //---------- CONSTRUCTORS
    Delaunay_triangulation_on_hyperbolic_surface_2(Traits & gt) : gt_(gt) {};
    Delaunay_triangulation_on_hyperbolic_surface_2(Traits & gt, CMap & cmap, Anchor & anch);
    Delaunay_triangulation_on_hyperbolic_surface_2(Traits & gt, Hyperbolic_fundamental_domain_2<Traits> const & domain);
    Delaunay_triangulation_on_hyperbolic_surface_2(Traits & gt, Base & triangulation);

    //---------- UTILITIES
    Anchor & anchor(Dart_descriptor const dart);
    Anchor const & anchor(Dart_const_descriptor const dart) const;
    Anchor & anchor();
    Anchor const & anchor() const;
    unsigned index_in_anchor(Dart_const_descriptor const dart) const;
    Dart_descriptor ith_dart(unsigned i, Anchor const & anch);
    bool is_valid() const;

    //---------- location and insertion
    void relative_locate(Point const & query, Locate_type& lt, unsigned & li, Anchor const & anch) const;
    Anchor locate(Point const & query, Locate_walk walk = STRAIGHT); // const ?
    Anchor locate(Point const & query, Locate_type & lt, unsigned & li, unsigned & ld, Anchor const & hint, Locate_walk walk = STRAIGHT); // const ?
    void insert(Point const & query, Anchor & hint);
    void insert(Point const & query);

    //---------- eps-net methods
    bool construct_epsilon_net(double const epsilon);
    bool is_epsilon_covering(const double epsilon) const;
    bool is_epsilon_packing(const double epsilon) const;
    bool is_epsilon_net(const double epsilon) const;
    double packing_value() const;
    double covering_value() const;
    void set_circumcenter_approximation_precision(unsigned p) {gt_.approximation_precision(p);}
    unsigned get_circumcenter_approximation_precision() {return gt_.approximation_precision();}
    // If any, its length gives an upper bound on the systole. If there
    //     is no loop edge eps<arcsinh1, then systole > 2*eps
    Dart_const_descriptor shortest_loop_edge() const;
    Number dart_cosh_length(Dart_const_descriptor dart) const;

    // undocumented
    void to_stream(std::ostream & s) const;
    void from_stream(std::istream & s);
    std::ostream& to_json(std::ostream& s) const;

 private:
    Traits gt_;
    double epsilon_;
    unsigned const NULL_INDEX = -1;
    unsigned const NB_SIDES = 3;
    unsigned const DOUBLE_PREC = 53;

    //---------- CONSTRUCTORS
    void set_anchors();

    //---------- UTILITIES
    unsigned ccw(unsigned i) const;
    unsigned cw(unsigned i) const;
    void set_attribute(Dart_descriptor dart, Anchor const & anchor);
    void set_attribute(Dart_descriptor dart, Complex_number const & cross_ratio);
    void update_infos(Dart_descriptor dart, Point const & r, Point const & s, Point const & t, Point const & query);

    //---------- location and insertion
    Anchor locate_visibility_walk(Point const & query, Locate_type & lt, unsigned & li, unsigned & ld, Anchor const & hint); // const ?
    Anchor locate_straight_walk(Point const & query, Locate_type & lt, unsigned & li, unsigned & ld, Anchor const & hint); // const ?
    std::vector<Dart_descriptor> insert_in_face(Point const & query, Anchor & anch);  // return dart_descriptor
    std::vector<Dart_descriptor> insert_in_edge(Point const & query, unsigned & li, Anchor & anch);
    std::vector<Dart_descriptor> split_insert(Point const & query, Anchor & anch, Locate_walk walk = STRAIGHT);

    //---------- Delaunay related methods
    void flip(Dart_descriptor dart);
    unsigned make_Delaunay();
    void push_flippable_edge(Dart_descriptor const dart, std::list<Dart_descriptor> & darts_to_flip);
    unsigned restore_Delaunay(std::list<Dart_descriptor> & darts_to_flip, std::vector<Dart_descriptor> & flipped_darts);

    //---------- eps-net methods
    Voronoi_point circumcenter(Anchor const & anch) const;
    Point approx_circumcenter_from_anchor(Anchor const & anch) const;
    //  Point approx_circumcenter_from_c(Voronoi_point const & c) const;
    void push_triangle(Dart_descriptor const dart, std::list<Dart_descriptor> & triangles, size_t & triangles_list_mark);
};


//---------- CONSTRUCTORS

template<class Traits>
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
Delaunay_triangulation_on_hyperbolic_surface_2( Traits & gt, CMap & cmap, Anchor & anch)
  : Base(cmap, anch), gt_(gt)
{
    Base::make_Delaunay();
    set_anchors();
}

template<class Traits>
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
Delaunay_triangulation_on_hyperbolic_surface_2( Traits & gt, Hyperbolic_fundamental_domain_2<Traits> const & domain)
  : Base(domain), gt_(gt)
{
    Base::make_Delaunay();
    set_anchors();
}

template<class Traits>
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
Delaunay_triangulation_on_hyperbolic_surface_2( Traits & gt, Base & triangulation)
: Base(triangulation.combinatorial_map(), triangulation.anchor()), gt_(gt)
{
    Base::make_Delaunay();
    set_anchors();
}

// only used in constructors
template<class Traits>
void
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
set_anchors()
{
    CMap & cmap = this->combinatorial_map_;

    std::queue<Anchor> bfs_queue;
    size_t in_queue = cmap.get_new_mark();
    cmap.unmark_all(in_queue);
    bfs_queue.push(this->anchor_);
    cmap.mark(this->anchor_.dart, in_queue);
    cmap.mark(Base::ccw(this->anchor_.dart), in_queue);
    cmap.mark(Base::cw(this->anchor_.dart), in_queue);

    while (!bfs_queue.empty()) {
        Anchor & current = bfs_queue.front();
        set_attribute(current.dart, current);
        Dart_descriptor invader = current.dart;
        for (int i = 0; i < 3; i++) {
            Dart_descriptor invaded = Base::opposite(invader);
            if (!cmap.is_marked(invaded, in_queue)) {
                Complex_number cross_ratio = Base::get_cross_ratio(invader);
                Point & c = current.vertices[i % 3];
                Point & a = current.vertices[(i + 1) % 3];
                Point & b = current.vertices[(i + 2) % 3];
                CGAL_assertion(gt_.is_Delaunay_hyperbolic_2_object()(a, b, c));
                Point d = Base::fourth_point_from_cross_ratio(a, b, c, cross_ratio);
                CGAL_assertion(norm(Complex_number(d.x(), d.y())) < Number(1));
                bfs_queue.push(Anchor(invaded, a, c, d));
                cmap.mark(invaded, in_queue);
                cmap.mark(Base::ccw(invaded), in_queue);
                cmap.mark(Base::cw(invaded), in_queue);
            }
            invader = Base::ccw(invader);
        }
        bfs_queue.pop();
    }
    cmap.free_mark(in_queue);

    this->anchor_ = Anchor();
    this->has_anchor_ = false;

    CGAL_assertion(is_valid());
}


//---------- UTILITIES

template<class Traits>
typename Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::Anchor&
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
anchor(Dart_descriptor const dart)
{
    return this->combinatorial_map_.template info<2>(dart);
}

template<class Traits>
typename Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::Anchor const &
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
anchor(Dart_const_descriptor const dart) const
{
    return this->combinatorial_map_.template info<2>(dart);
}

template<class Traits>
typename Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::Anchor &
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
anchor()
{
    return anchor(this->combinatorial_map_.darts().begin());
}

template<class Traits>
typename Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::Anchor const &
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
anchor() const
{
    return anchor(this->combinatorial_map_.darts().begin());
}

template<class Traits>
unsigned
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
index_in_anchor(Dart_const_descriptor const dart) const
{
    unsigned index = NULL_INDEX;
    Anchor const & anch = anchor(dart);
    Dart_const_descriptor current_dart = anch.dart;
    for (unsigned i = 0; i < NB_SIDES; ++i) {
        if (current_dart == dart) {
            index = i;
            break;
        }
        current_dart = Base::const_ccw(current_dart);
    }
    CGAL_assertion(index != NULL_INDEX);
    return index;
}

template<class Traits>
typename Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::Dart_descriptor
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
ith_dart(unsigned i, Anchor const & anch)
{
    Dart_descriptor dart = anch.dart;
    for (unsigned j = 0; j < i; ++j) {
        dart = Base::ccw(dart);
    }
    return dart;
}

template<class Traits>
unsigned
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
ccw(unsigned i) const
{
    CGAL_precondition( i <= NB_SIDES);
    return (i + 1) % NB_SIDES;
}

template<class Traits>
unsigned
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
cw(unsigned i) const
{
    CGAL_precondition( i <= NB_SIDES);
    return (i + 2) % NB_SIDES;
}

template<class Traits>
void
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
set_attribute(Dart_descriptor dart, Anchor const & anchor)
{
    this->combinatorial_map_.template set_attribute<2>(dart,
        this->combinatorial_map_.template create_attribute<2>(anchor));
}

template<class Traits>
void
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
set_attribute(Dart_descriptor dart, Complex_number const & cross_ratio)
{
    this->combinatorial_map_.template set_attribute<1>(dart,
        this->combinatorial_map_.template create_attribute<1>(cross_ratio));
}

template<class Traits>
bool
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
is_valid() const
{
    if (!Base::is_Delaunay()) {
        return false;
    }

    for (typename Face_const_range::const_iterator it = Base::faces_const_range().begin(); it != Base::faces_const_range().end(); ++it) {
        Anchor const & current = anchor(it);
        Dart_const_descriptor current_dart = current.dart;

        for (unsigned i = 0; i < NB_SIDES; ++i) {
            Dart_const_descriptor opposite_dart = Base::const_opposite(current_dart);
            Point const & c1 = current.vertices[i];
            Point const & a1 = current.vertices[ccw(i)];
            Point const & b1 = current.vertices[cw(i)];
            CGAL_assertion(norm(Complex_number(a1.x(), a1.y())) < Number(1));
            CGAL_assertion(norm(Complex_number(b1.x(), b1.y())) < Number(1));
            CGAL_assertion(norm(Complex_number(c1.x(), c1.y())) < Number(1));
            Complex_number cross_ratio = Base::get_cross_ratio(current_dart);
            Point d1 = Base::fourth_point_from_cross_ratio(a1, b1, c1, cross_ratio);

            unsigned j = index_in_anchor(opposite_dart);
            Anchor const & neighbor = anchor(opposite_dart);
            Point const & a2 = neighbor.vertices[j];
            Point const & c2 = neighbor.vertices[ccw(j)];
            Isometry pair_sides = isometry_pairing_the_sides<Traits>(a2, c2, a1, c1);
            CGAL_assertion(pair_sides.evaluate(a2) == a1);
            CGAL_assertion(pair_sides.evaluate(c2) == c1);

            Point d2 = pair_sides.evaluate(neighbor.vertices[cw(j)]);
            CGAL_assertion(d2 == d1);
            if (d2 != d1) {
                return false;
            }
            current_dart = Base::const_ccw(current_dart);
        }
    }
    return true;
}

//---------- location and insertion

//TODO DOC
// Computes the locate type lt of query relative to the anchor, and an index corresponding to:
// - if lt == FACE: NULL_INDEX,
// - if lt == EDGE: index of the edge on which query lies,
// - if lt == VERTEX: index of the vertex on which query lies,
// - if lt == OUTSIDE: index of the first edge such that query and the third point of the triangle lies on different sides.
template<class Traits>
void
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
  relative_locate(Point const & query, Locate_type& lt, unsigned & li, Anchor const & anch) const {
    lt = FACE;
    li = NULL_INDEX;
    //Traits gt;//MARC?? TO BE REMOVED
    typename Traits::Hyperbolic_orientation_2 ho2 = gt_.hyperbolic_orientation_2();
    for (unsigned i = 0; i < NB_SIDES; ++i) {
        Orientation ori_query = ho2(anch.vertices[i], anch.vertices[ccw(i)], query);
        if (ori_query == RIGHT_TURN) {
            lt = OUTSIDE;
            li = i;
            break;
        }
        if (ori_query == COLLINEAR) {
            lt = EDGE;
            li = i;
            if (ho2(anch.vertices[ccw(i)], anch.vertices[cw(i)], query) == COLLINEAR) {
                lt = VERTEX;
                li = ccw(i);
                break;
            }
        }
    }
}

// Output: an anchor of the triangle in which query lies,
// and an int corresponding to the number of traversed triangles to find it.
template<class Traits>
typename Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::Anchor
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
locate_visibility_walk(Point const & query, Locate_type & lt, unsigned & li, unsigned & ld, Anchor const & hint)
{
    CGAL_precondition(norm(Complex_number(query.x(), query.y())) < Number(1));
    CGAL_expensive_precondition(is_valid());

    // initialisation
    ld = 0;
    relative_locate(query, lt, li, hint);
    if (lt != OUTSIDE){
        return hint;
    }

    bool found = false;
    Dart_descriptor dart = ith_dart(li, hint);
    Point c = hint.vertices[li];
    Point a = hint.vertices[ccw(li)];
    Point b = hint.vertices[cw(li)];
    Point d;

    // visibility walk
    //Traits gt;
    typename Traits::Hyperbolic_orientation_2 ho2 = gt_.hyperbolic_orientation_2();
    while (!found) {
        Complex_number cross_ratio = Base::get_cross_ratio(dart);
        d = Base::fourth_point_from_cross_ratio(a, b, c, cross_ratio);
        dart = Base::opposite(dart);
        if (ho2(c, d, query) == RIGHT_TURN){
            b = a;
            a = d;
            dart = Base::ccw(dart);
        } else if (ho2(a, d, query) == LEFT_TURN) {
            b = c;
            c = d;
            dart = Base::cw(dart);
        } else {
            found = true;
        }
        ++ld;
    }
    Anchor res = Anchor(dart, a, c, d);
    relative_locate(query, lt, li, res);
    CGAL_assertion(lt != OUTSIDE);
    return res;
}

template<class Traits>
typename Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::Anchor
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
locate_straight_walk(Point const & query, Locate_type & lt, unsigned & li, unsigned & ld, Anchor const & hint)
{
    CGAL_precondition(norm(Complex_number(query.x(), query.y())) < Number(1));

    // initialisation
    Dart_descriptor dart = hint.dart;
    Point p0 = hint.vertices[0];
    Point r = hint.vertices[1];
    Point l = hint.vertices[2];
    ld = 0;
    //Traits gt;
    typename Traits::Hyperbolic_orientation_2 ho2 = gt_.hyperbolic_orientation_2();

    if (ho2(p0, query, r) == RIGHT_TURN) {
        while (ho2(p0, query, l) == RIGHT_TURN) {
            dart = Base::opposite(Base::cw(dart));
            Point old_l = l;
            l = Base::fourth_point_from_cross_ratio(p0, r, l, Base::get_cross_ratio(dart));
            r = old_l;
            ++ld;
        }
    } else {
        while (ho2(p0, query, r) == LEFT_TURN) {
            Point old_r = r;
            r = Base::fourth_point_from_cross_ratio(r, l, p0, Base::get_cross_ratio(dart));
            l = old_r;
            dart = Base::ccw(Base::opposite(dart));
            ++ld;
        }
    }

    // straight walk
    dart = Base::ccw(dart);
    Point p = p0;
    while (ho2(r, l, query) == RIGHT_TURN) {
        Complex_number cross_ratio = Base::get_cross_ratio(dart);
        Point s = Base::fourth_point_from_cross_ratio(l, p, r, cross_ratio);
        if (ho2(p0, query, s) == RIGHT_TURN) {
            p = r;
            r = s;
            dart = Base::cw(Base::opposite(dart));
        } else {
            p = l;
            l = s;
            dart = Base::ccw(Base::opposite(dart));
        }
        ++ld;
    }
    Anchor res = Anchor(dart, r, l, p);
    relative_locate(query, lt, li, res);
    CGAL_assertion(lt != OUTSIDE);
    return res;
}

template<class Traits>
typename Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::Anchor
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
locate(Point const & query, Locate_walk walk)
{
    Locate_type lt = OUTSIDE;
    unsigned li = 0;
    unsigned ld = 0;
    return locate(query, lt, li, ld, anchor(), walk);
}

template<class Traits>
typename Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::Anchor
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
locate(Point const & query, Locate_type & lt, unsigned & li, unsigned & ld, Anchor const & hint, Locate_walk walk)
{
    if (walk == VISIBILITY) {
        return locate_visibility_walk(query, lt, li, ld, hint);
    } else {
        return locate_straight_walk(query, lt, li, ld, hint);
    }
}

// Input: dart whose cross-ratio is those of the edge [t, r] and is computed with s and a fourth vertex.
// Output: modified anchor such that its vertices are {r, query, t} and its dart represents the edge between r and t.
// The cross-ratio of the edge [t, r] is updated with its value where query is replaced by s.
// The cross-ratio of the edge [r, q] is set.
// Warning: the other cross-ratios are not updated, they will be in insert_in_edge and insert_in_face.
template<class Traits>
void
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
update_infos(Dart_descriptor dart, Point const & r, Point const & s, Point const & t, Point const & query)
{
    // compute new cross-ratio
    Complex_number old_cr = Base::get_cross_ratio(dart);
    Point u = Base::fourth_point_from_cross_ratio(r, s, t, old_cr);
    Complex_number new_cr = Base::cross_ratio(r, query, t, u);

    // modify anch and cross-ratio
    std::shared_ptr<Anchor> anch = std::make_shared<Anchor>(anchor(dart));
    anch->dart = dart;
    this->combinatorial_map_.template info<1>(dart) = new_cr;
    anch->vertices[0] = t;
    anch->vertices[1] = r;
    anch->vertices[2] = query;
    this->combinatorial_map_.template info<2>(dart) = *anch;

    set_attribute(Base::ccw(dart), Base::cross_ratio(query, t, r, s));
}

// Output: the three anchors incident to query after its insertion inside anch.
// the darts of the new anchors correspond to the edges of the triangle the point was inserted in
// Note to self: No need to manage triangulation's anchor as it is done in update_infos.
template<class Traits>
std::vector<typename Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::Dart_descriptor>
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
insert_in_face(Point const & query, Anchor & anch)
{
    Dart_descriptor current_dart = anch.dart;
    this->combinatorial_map_.insert_cell_0_in_cell_2(anch.dart);
    std::vector<Dart_descriptor> darts_of_new_anchors;
    for (unsigned i = 0; i < NB_SIDES; ++i) {
        Point & c = anch.vertices[i];
        Point & a = anch.vertices[ccw(i)];
        Point & b = anch.vertices[cw(i)];
        update_infos(current_dart, a, b, c, query);
        darts_of_new_anchors.push_back(current_dart);
        current_dart = Base::ccw(Base::opposite(Base::ccw(current_dart)));
    }

    return darts_of_new_anchors;
}

// Output: the four anchors incident to query after its insertion on one edge of anch.
// the darts of the new anchors correspond to the edges of the triangle the point was inserted in
template<class Traits>
std::vector<typename Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::Dart_descriptor>
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
insert_in_edge(Point const & query, unsigned & li, Anchor & anch)
{
    // find dart on which we insert query
    Dart_descriptor insertion_dart = ith_dart(li, anch);

    // gather information
    std::vector<Dart_descriptor> darts_of_new_anchors;
    Point & c = anch.vertices[li];
    Point & a = anch.vertices[ccw(li)];
    Point & b = anch.vertices[cw(li)];
    Complex_number cross_ratio = Base::get_cross_ratio(insertion_dart);
    Point d = Base::fourth_point_from_cross_ratio(a, b, c, cross_ratio);

    Dart_descriptor dart_ab = Base::ccw(insertion_dart);
    Dart_descriptor dart_bc = Base::ccw(dart_ab);
    Dart_descriptor dart_cd = Base::ccw(Base::opposite(insertion_dart));
    Dart_descriptor dart_da = Base::ccw(dart_cd);

    // insert vertex on edge in the cmap and create new triangles
    this->combinatorial_map_.insert_cell_0_in_cell_1(insertion_dart);
    this->combinatorial_map_.insert_cell_1_in_cell_2(dart_bc, Base::cw(dart_ab));
    this->combinatorial_map_.insert_cell_1_in_cell_2(dart_da, Base::cw(dart_cd));
    CGAL_assertion(Base::ccw(dart_ab) == Base::opposite(Base::cw(dart_bc)));
    CGAL_assertion(Base::ccw(dart_bc) == Base::opposite(Base::cw(dart_cd)));
    CGAL_assertion(Base::ccw(dart_cd) == Base::opposite(Base::cw(dart_da)));
    CGAL_assertion(Base::ccw(dart_da) == Base::opposite(Base::cw(dart_ab)));

    // set new anchors and cross-ratios
    update_infos(dart_ab, b, c, a, query);
    update_infos(dart_bc, c, a, b, query);
    update_infos(dart_cd, d, a, c, query);
    update_infos(dart_da, a, c, d, query);

    darts_of_new_anchors.push_back(dart_ab);
    darts_of_new_anchors.push_back(dart_bc);
    darts_of_new_anchors.push_back(dart_cd);
    darts_of_new_anchors.push_back(dart_da);

    return darts_of_new_anchors;
}

// the darts of the new anchors correspond to the edges of the triangle the point was inserted in
template<class Traits>
std::vector<typename Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::Dart_descriptor>
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
split_insert(Point const & query, Anchor & anch, Locate_walk walk)
{
    Locate_type lt = OUTSIDE;
    unsigned li = 0;
    unsigned ld = 0;
    Anchor locate_anchor = locate(query, lt, li, ld, anch, walk);
    CGAL_precondition(lt != OUTSIDE);

    std::vector<Dart_descriptor> darts_of_new_anchors;
    if (lt == FACE) {
        darts_of_new_anchors = insert_in_face(query, locate_anchor);
    } else if (lt == EDGE) {
        darts_of_new_anchors = insert_in_edge(query, li, locate_anchor);
    }
    return darts_of_new_anchors;
}


//---------- Delaunay related methods

template<class Traits>
void
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
flip(Dart_descriptor dart)
{
    CGAL_expensive_precondition(is_valid());
    CGAL_precondition(!Base::has_anchor());

    // first gather all the information needed
    Dart_descriptor opposite = Base::opposite(dart);
    Anchor & anch = anchor(dart);

    unsigned index = index_in_anchor(dart);
    Point & C = anch.vertices[index];
    Point & A = anch.vertices[ccw(index)];
    Point & B = anch.vertices[cw(index)];
    Complex_number cross_ratio = Base::get_cross_ratio(dart);
    Point D = Base::fourth_point_from_cross_ratio(A, B, C, cross_ratio);

    // flip the edge in the cmap and update the cross-ratios
    Base::flip(dart);  // AC -> BD

    // update the two anchors
    // WARNING: do not update dart's anchor first,
    // else it will change vertices' values (due to refs) and cause errors
    this->combinatorial_map_.template info<2>(opposite) = Anchor(opposite, D, B, C);
    this->combinatorial_map_.template info<2>(dart) = Anchor(dart, B, D, A);

    CGAL_expensive_assertion(is_valid());
}

// Pushes dart in the list darts_to_flip if dart is flippable and is not already in the list
template<class Traits>
void
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
push_flippable_edge(Dart_descriptor const dart, std::list<Dart_descriptor>& darts_to_flip)
{
    if (Base::is_Delaunay_flippable(dart)) {
        bool already_there = false;
        for (Dart_descriptor const & dart_to_flip : darts_to_flip) {
            if (dart_to_flip == dart || dart_to_flip == Base::opposite(dart)) {
                already_there = true;
                break;
            }
        }
        if (!already_there) {
            darts_to_flip.push_back(dart);
        }
    }
}

// Output: number of flips done to make the triangulation Delaunay given a list of darts to flip
template<class Traits>
unsigned
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
restore_Delaunay(std::list<Dart_descriptor> & darts_to_flip, std::vector<Dart_descriptor> & flipped_darts)
{
    unsigned nb_of_flips_done = 0;
    while (!darts_to_flip.empty()) {
        Dart_descriptor current_dart = darts_to_flip.front();
        if (Base::is_Delaunay_flippable(current_dart)) {
            flip(current_dart);
            flipped_darts.push_back(current_dart);
            ++nb_of_flips_done;
            Dart_descriptor maybe_flippable[4] = {Base::ccw(current_dart), Base::cw(current_dart),
                                              Base::ccw(Base::opposite(current_dart)), Base::cw(Base::opposite(current_dart))};
            for (int i = 0; i < 4; ++i) {
                push_flippable_edge(maybe_flippable[i], darts_to_flip);
            }
        }
        darts_to_flip.pop_front();
    }

    CGAL_expensive_assertion(is_valid());
    CGAL_expensive_assertion(Base::is_Delaunay());
    return nb_of_flips_done;
}

template<class Traits>
unsigned
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
make_Delaunay()
{
    unsigned number_of_flips_done = 0;

    Dart_descriptor edge_to_flip = Base::pick_edge_to_flip();
    while (edge_to_flip != nullptr) {
        flip(edge_to_flip);
        edge_to_flip = Base::pick_edge_to_flip();
        number_of_flips_done++;
    }

    CGAL_expensive_assertion(is_valid());
    CGAL_expensive_assertion(Base::is_Delaunay());
    return number_of_flips_done;
}

// Inserts query in the triangulation, with a search starting from the given anchor, and restores the Delaunay property with flips
template<class Traits>
void
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
insert(Point const & query, Anchor & hint)
{
    CGAL_precondition(norm(Complex_number(query.x(), query.y())) < Number(1));
    CGAL_expensive_precondition(is_valid());

    std::vector<Dart_descriptor> darts_of_new_anchors = split_insert(query, hint);
    std::list<Dart_descriptor> darts_to_flip;
    for (Dart_descriptor dart : darts_of_new_anchors) {
        push_flippable_edge(darts_of_new_anchors, darts_to_flip);
        push_flippable_edge(Base::ccw(darts_of_new_anchors), darts_to_flip);
    }
    std::vector<Dart_descriptor> flipped_darts;
    restore_Delaunay(darts_to_flip, flipped_darts);
}

// Same but the search starts from main anchor
template<class Traits>
void
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
insert(Point const & query)
{
    insert(query, anchor());
}

//---------- eps-net methods

template<class Traits>
typename Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::Voronoi_point
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
circumcenter(Anchor const & anch) const
{ CGAL_precondition(gt_.is_Delaunay_hyperbolic_2_object()(anch.vertices[0], anch.vertices[1], anch.vertices[2]));
  typename Traits::Construct_hyperbolic_circumcenter_2 chc = gt_.construct_hyperbolic_circumcenter_2_object();
  return chc(anch.vertices[0], anch.vertices[1], anch.vertices[2]);
}

template<class Traits>
typename Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::Point
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
approx_circumcenter_from_anchor(Anchor const & anch) const
{
  typename Traits::Construct_approximate_hyperbolic_circumcenter_2 cahc = gt_.construct_approximate_hyperbolic_circumcenter_2_object();
  return  cahc(anch.vertices[0], anch.vertices[1], anch.vertices[2]);
}

/*  TO BE REMOVED
template<class Traits>
typename Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::Point
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
approx_circumcenter_from_c(Voronoi_point const & c) const
{
  //TEST
  typename Traits::Construct_approximate_hyperbolic_circumcenter_2_fct_style obj_approx_cc = gt_.construct_approximate_hyperbolic_circumcenter_2_object();
  obj_approx_cc(c);
  //END TEST
  return gt_.Construct_approximate_hyperbolic_circumcenter_2_non_fct(c);
}
*/

template<class Traits>
void
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
push_triangle(Dart_descriptor const dart, std::list<Dart_descriptor> & triangles, size_t & triangles_list_mark)
{
    this->combinatorial_map_.unmark(Base::ccw(dart), triangles_list_mark);
    this->combinatorial_map_.unmark(Base::cw(dart), triangles_list_mark);
    this->combinatorial_map_.mark(dart, triangles_list_mark);
    triangles.push_back(dart);
}

template<class Traits>
bool
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
construct_epsilon_net(double const epsilon)
{
    epsilon_ = epsilon;
    CGAL_precondition(is_epsilon_packing(epsilon));
    Interval const cosh_epsilon_interval = std::cosh(epsilon);
    Number const BOUND = upper(cosh_epsilon_interval);
    size_t triangles_list_mark = this->combinatorial_map_.get_new_mark();

    std::list<Dart_descriptor> triangles;
    for (typename Face_range::iterator it = this->combinatorial_map_.template one_dart_per_cell<2>().begin();
        it != this->combinatorial_map_.template one_dart_per_cell<2>().end(); ++it) {
        push_triangle(it, triangles, triangles_list_mark);
    }

    while (!triangles.empty()) {
        Dart_descriptor current_dart = triangles.front();
        Anchor & current_anchor = anchor(current_dart);
        triangles.pop_front();
        if(this->combinatorial_map_.is_marked(current_dart, triangles_list_mark)){
	  this->combinatorial_map_.unmark(current_dart, triangles_list_mark);
	  Voronoi_point c = circumcenter(current_anchor);
	  if (gt_.template cosh_hd<Algebraic_number,Voronoi_point>(c, current_anchor.vertices[0]) > BOUND) {
	    // Point approx_c = approx_circumcenter_from_c(c); TO BE REMOVED
	    Point approx_c = approx_circumcenter_from_anchor(current_anchor);
	    if(norm(Complex_number(approx_c.x(), approx_c.y())) >= Number(1)) {
	      break;  // avoid undefined behavior in case of bad approx outside of Poincaré
	    }
	    std::vector<Dart_descriptor> darts_of_new_anchors = split_insert(approx_c, current_anchor, VISIBILITY);
	    std::list<Dart_descriptor> darts_to_flip;
	    for (Dart_descriptor const & dart : darts_of_new_anchors) {
	      push_triangle(dart, triangles, triangles_list_mark);
	      push_flippable_edge(dart, darts_to_flip);
	      push_flippable_edge(Base::ccw(dart), darts_to_flip);
	    }

	    std::vector<Dart_descriptor> flipped_darts;
	    restore_Delaunay(darts_to_flip, flipped_darts);
	    for (Dart_descriptor const & dart : flipped_darts) {
	      push_triangle(dart, triangles, triangles_list_mark);
	      push_triangle(Base::opposite(dart), triangles, triangles_list_mark);
	    }
	  }
        }
    }
    this->combinatorial_map_.free_mark(triangles_list_mark);

    bool is_eps_net = is_epsilon_net(epsilon);
    if (!is_eps_net) {std::cout<< "WARNING: THE CONSTRUCTION OF THE EPSILON NET FAILED, THE APPROXIMATION PRECISION SHOULD BE INCREASED\n";}
    return is_eps_net;
}

template<class Traits>
bool
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
is_epsilon_covering(const double epsilon) const
{
  Interval eps = epsilon;
  bool is_covering = (upper(eps) >= covering_value());
  return is_covering;
}

template<class Traits>
bool
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
is_epsilon_packing(const double epsilon) const
{
  Interval eps = epsilon;
  bool is_packing = (lower(eps) <= packing_value());
  return is_packing;
}

template<class Traits>
bool
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
is_epsilon_net(const double epsilon) const
{
  return is_epsilon_covering(epsilon) && is_epsilon_packing(epsilon);
}

//TODO DOC
template<class Traits>
typename Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::Dart_const_descriptor
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
shortest_loop_edge() const
{
  Number   length;
  Number   min_length;
  Dart_const_descriptor shortest_loop_dart = nullptr;

  for (typename Edge_const_range::const_iterator it = this->combinatorial_map_.template one_dart_per_cell<1>().begin();
       it != this->combinatorial_map_.template one_dart_per_cell<1>().end(); ++it)
    {
      // check if the edge is a loop
      Dart_const_descriptor next = Base::const_ccw(it);
      auto doc = this->combinatorial_map_.template darts_of_cell<0>(it);
      for (auto dart = doc.begin(); dart != doc.end(); ++dart)
	{
	  if (next == dart) { // the edge is a loop
	    Number length = dart_cosh_length(it);
	    if (shortest_loop_dart == nullptr) {min_length = length;}
	    if (length < min_length) {
	      min_length = min(min_length,length);
	      shortest_loop_dart = it;
	      /*
	      	//TEST alternative loop test? May not be valid for a Cmap without vertices explicitly defined?
	Dart_const_descriptor d =  it;
	// Next dart in the face
	Dart_const_descriptor d2 = this->combinatorial_map_.beta(d, 1);
	// Check whether both darts belong to the same vertex
	bool is_loop2 =
	  this->combinatorial_map_.template belong_to_same_cell<0>(d, d2);
	std::cout<< "loop2 value, loop " <<is_loop2  << std::endl;
	      */
	    }
	  }
	}
    }
  return shortest_loop_dart;
}

//TODO DOC
template<class Traits>
double
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
covering_value() const
{
  Algebraic_number radius;
  Algebraic_number max_radius = 0;
  double res_double = 0.;

  for (typename Face_const_range::const_iterator it = this->combinatorial_map_.template one_dart_per_cell<2>().begin();
       it != this->combinatorial_map_.template one_dart_per_cell<2>().end(); ++it) {
    Dart_const_descriptor current_dart = it;
    const Anchor & current_anchor = anchor(current_dart);
    Voronoi_point c = circumcenter(current_anchor);
    radius = gt_.template cosh_hd<Algebraic_number,Voronoi_point>(c, current_anchor.vertices[0]);
    if (radius > max_radius) {
      max_radius = radius;
    }
  }
  //upper bound on acosh(max_radius), max_radius= a0 + a1 * sqrt(root), with a0,a1,root FT type
  if constexpr(std::is_same_v<Number, Gmpq>){
      mpfr_t res, a0, a1, root, root_sqrt;
      mpfr_set_default_prec(53);
      mpfr_init (res);
      mpfr_init (a0);
      mpfr_init (a1);
      mpfr_init (root);
      mpfr_init (root_sqrt);
      mpfr_set_q (a0, max_radius.a0().mpq(), MPFR_RNDU);
      mpfr_set_q (a1, max_radius.a1().mpq(), MPFR_RNDU);
      if (max_radius.a1() > 0) {
	mpfr_set_q (root, max_radius.root().mpq(), MPFR_RNDU);
	mpfr_sqrt(root_sqrt, root, MPFR_RNDU);
      } else {
	mpfr_set_q (root, max_radius.root().mpq(), MPFR_RNDD);
	mpfr_sqrt(root_sqrt, root, MPFR_RNDD);
      }
      mpfr_mul(res, a1, root_sqrt, MPFR_RNDU);
      mpfr_add(res, a0, res, MPFR_RNDU);
      mpfr_acosh (res, res, MPFR_RNDU);
      res_double = mpfr_get_d (res, MPFR_RNDZ);
      mpfr_clear (res);
      mpfr_clear (a0);
      mpfr_clear (a1);
      mpfr_clear (root);
      mpfr_clear (root_sqrt);
    } else {
    res_double = std::acosh(to_double(max_radius));
  }
  return res_double;
}

template<class Traits>
double
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
packing_value() const
{
  double res_double = std::numeric_limits<double>::infinity();
  // If all edges are loops, that is there is only 1 vertex:
  std::vector<unsigned int> count_vertices = {0};
  if (this->combinatorial_map_.count_cells(count_vertices)[0] == 1)
    {return res_double;}

  Number min_length = 0;
  Number length;

  for (typename Edge_const_range::const_iterator it = this->combinatorial_map_.template one_dart_per_cell<1>().begin();
       it != this->combinatorial_map_.template one_dart_per_cell<1>().end(); ++it) {
    // check if the edge is a loop
    bool is_loop = false;
    Dart_const_descriptor next = Base::const_ccw(it);
    auto doc = this->combinatorial_map_.template darts_of_cell<0>(it);
    for (auto dart = doc.begin(); dart != doc.end(); ++dart) {
      if (next == dart) {
	is_loop = true;
      }
    }
    if(!is_loop) {
      length = dart_cosh_length(it);
      if (min_length == 0) {min_length = length;}
      min_length = min(min_length, length);
   }
  }

  if constexpr(std::is_same_v<Number, CGAL::Gmpq>){
    //   std::cout<< "GOOD: constexpr(std::is_same_v<Number, Gmpq>) is true" << std::endl; //TOBE REMOVED
    //mpfr version
    int p = 53;
    mpfr_t res;
    mpfr_init2 (res, p);
    mpfr_set_q (res, min_length.mpq(), MPFR_RNDZ);
    mpfr_acosh (res, res, MPFR_RNDZ);
    // Gmpfr version
    // Gmpfr res = Gmpfr(min_length.numerator(), p, std::round_toward_zero) / Gmpfr( min_length.denominator(), p, std::round_toward_pos_infinity);
    res_double = mpfr_get_d (res, MPFR_RNDZ);
    mpfr_clear (res);
    } else {
    res_double = std::acosh(to_double(min_length));
  }
  return res_double;
}

template<class Traits>
typename Traits::FT
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
dart_cosh_length(Dart_const_descriptor dart) const
{
  Anchor const & current_anchor = anchor(dart);
  unsigned index = index_in_anchor(dart);
  Point const & a = current_anchor.vertices[index];
  Point const & b = current_anchor.vertices[ccw(index)];
  return gt_.cosh_hd(a, b);
}

//////////////////////////////////////////////////////
//       TO/FROM streams
//////////////////////////////////////////////////////
template<class Traits>
void
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
to_stream(std::ostream & s) const
{
    CGAL_precondition(is_valid());

    // Give indices to the darts
    std::map<Dart_const_descriptor, unsigned> darts_ids;
    unsigned current_dart_id = 0;
    for (typename Dart_const_range::const_iterator it=this->combinatorial_map_.darts().begin(); it!=this->combinatorial_map_.darts().end(); ++it) {
        darts_ids[it] = current_dart_id;
        current_dart_id++;
    }

    // Store the number of darts
    s << current_dart_id << std::endl;

    // Store the triangles and their anchor
    for (typename Face_const_range::const_iterator it = Base::faces_const_range().begin(); it != Base::faces_const_range().end(); ++it) {
        s << darts_ids[it] << std::endl;
        s << darts_ids[Base::const_cw(it)] << std::endl;
        s << darts_ids[Base::const_ccw(it)] << std::endl;
        Anchor anch = anchor(it);
        s << darts_ids[anch.dart] << std::endl;
        for(int i = 0; i < NB_SIDES ; ++i) {
            s << anch.vertices[i].x() << std::endl;
            s << anch.vertices[i].y() << std::endl;
        }
    }

    // Store the edges and their cross-ratio
    for (typename Edge_const_range::const_iterator it = Base::edges_const_range().begin(); it != Base::edges_const_range().end(); ++it) {
        s << darts_ids[it] << std::endl;
        s << darts_ids[Base::const_opposite(it)] << std::endl;
        s << Base::get_cross_ratio(it);
    }
}

template<class Traits>
void
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
from_stream(std::istream & s)
{
    this->combinatorial_map_.clear();

    // Load the number of darts
    std::string line;
    s >> line;
    int nb_darts = std::stoi(line);

    // Load the triangles and the anchors
    Dart_descriptor darts_by_id[nb_darts];
    for (int k = 0; k < nb_darts / NB_SIDES; ++k) {
        Dart_descriptor triangle_dart = this->combinatorial_map_.make_combinatorial_polygon(3);
        s >> line;
        int id_0 = std::stoi(line);
        s >> line;
        int id_1 = std::stoi(line);
        s >> line;
        int id_2 = std::stoi(line);
        darts_by_id[id_0] = triangle_dart;
        darts_by_id[id_1] = Base::cw(triangle_dart);
        darts_by_id[id_2] = Base::ccw(triangle_dart);

        Anchor anch = Anchor();
        s >> line;
        anch.dart = darts_by_id[std::stoi(line)];
        for(int i = 0; i < NB_SIDES; ++i) {
            Number x;
            s >> x;
            Number y;
            s >> y;
            anch.vertices[i] = Point(x, y);
        }
        set_attribute(anch.dart, anch);
  }

    // Load the edges
    for (int k = 0; k < nb_darts / 2; ++k) {
        s >> line;
        int id_0 = std::stoi(line);
        s >> line;
        int id_1 = std::stoi(line);
        Dart_descriptor dart_0 = darts_by_id[id_0];
        Dart_descriptor dart_1 = darts_by_id[id_1];
        this->combinatorial_map_.template sew<2>(dart_0, dart_1);
        Complex_number cross_ratio;
        s >> cross_ratio;
        set_attribute(dart_0, cross_ratio);
  }

  CGAL_assertion(is_valid());
}

//////////////////////////////////////////////////////
//       TO JSON OUTPUT
//////////////////////////////////////////////////////
template<class Traits>
std::ostream&
Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::
to_json(std::ostream& s) const
{
    CGAL_precondition(is_valid());

    std::map<Dart_const_descriptor, unsigned> darts_ids;
    unsigned current_dart_id = 0;

    for (typename Dart_const_range::const_iterator it =
             this->combinatorial_map_.darts().begin();
         it != this->combinatorial_map_.darts().end();
         ++it)
    {
        darts_ids[it] = current_dart_id++;
    }

    s << "{\n";
    s << "  \"type\": " << "\"Delaunay_triangulation\"" << ",\n";
    s << "  \"num_darts\": " << current_dart_id << ",\n";

    // Faces
    s << "  \"faces\": [\n";

    bool first_face = true;

    for (typename Face_const_range::const_iterator it =
             Base::faces_const_range().begin();
         it != Base::faces_const_range().end();
         ++it)
    {
        if (!first_face) s << ",\n";
        first_face = false;

        Anchor anch = anchor(it);

        s << "    {\n";
        s << "      \"dart\": " << darts_ids[it] << ",\n";
        s << "      \"cw\": " << darts_ids[Base::const_cw(it)] << ",\n";
        s << "      \"ccw\": " << darts_ids[Base::const_ccw(it)] << ",\n";
        s << "      \"anchor_dart\": " << darts_ids[anch.dart] << ",\n";

        s << "      \"vertices\": [\n";

        for (int i = 0; i < NB_SIDES; ++i)
        {
            if (i > 0) s << ",\n";

            s << "        { "
              << "\"x\": \"" << anch.vertices[i].x()
              << "\", \"y\": \"" << anch.vertices[i].y()
              << "\" }";
        }

        s << "\n      ]\n";
        s << "    }";
    }

    s << "\n  ]\n";

    /*
    // Edges
    s << "  \"edges\": [\n";

    bool first_edge = true;

    for (typename Edge_const_range::const_iterator it =
             Base::edges_const_range().begin();
         it != Base::edges_const_range().end();
         ++it)
    {
        if (!first_edge) s << ",\n";
        first_edge = false;

        s << "    {\n";
        s << "      \"dart\": " << darts_ids[it] << ",\n";
        s << "      \"opposite\": "
          << darts_ids[Base::const_opposite(it)] << ",\n";
        s << "      \"cross_ratio\": \""
          << Base::get_cross_ratio(it) << "\"\n";
        s << "    }";
    }

    s << "\n  ]\n";
    */
    s << "}";


    return s;
}

}  // namespace CGAL

#endif  //CGAL_DELAUNAY_TRIANGULATION_ON_HYPERBOLIC_SURFACE_2
