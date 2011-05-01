// Copyright (c) 2009 INRIA Sophia-Antipolis (France),
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)    : Samuel Hornus

#ifndef CGAL_DELAUNAY_COMPLEX_H
#define CGAL_DELAUNAY_COMPLEX_H

#include <CGAL/Pure_complex.h>
#include <CGAL/Dimension.h>
#include <CGAL/Default.h>

namespace CGAL {

template< typename DCTraits, typename _PCDS = Default >
class Delaunay_complex
: public Pure_complex<DCTraits,
            typename Default::Get<_PCDS, Pure_complex_data_structure<
                             typename Ambient_dimension<typename DCTraits::Point_d>::type,
                             Pure_complex_vertex<DCTraits>,
                             Pure_complex_simplex<DCTraits> >
                    >::type >
{
    typedef typename Ambient_dimension<typename DCTraits::Point_d>::type
                                                    Ambient_dimension_;
    typedef typename Default::Get<_PCDS, Pure_complex_data_structure<
                         Ambient_dimension_,
                         Pure_complex_vertex<DCTraits>,
                         Pure_complex_simplex<DCTraits> >
                >::type                         PCDS;
    typedef Pure_complex<DCTraits, PCDS>        Base;
    typedef Delaunay_complex<DCTraits, _PCDS>   Self;

    typedef typename DCTraits::Side_of_oriented_subsphere_d
                                                    Side_of_oriented_subsphere_d;
    typedef typename DCTraits::Side_of_oriented_sphere_d
                                                    Side_of_oriented_sphere_d;
    typedef typename DCTraits::Coaffine_orientation_d
                                                    Coaffine_orientation_d;
    typedef typename DCTraits::Orientation_d        Orientation_d;

public: // PUBLIC NESTED TYPES

    typedef DCTraits                                Geom_traits;
    typedef typename Base::Pure_complex_ds          Pure_complex_ds;

    typedef typename Base::Vertex                   Vertex;
    typedef typename Base::Simplex                  Simplex;
    typedef typename Base::Facet                    Facet;
    typedef typename Base::Face                     Face;

    typedef typename Base::Ambient_dimension        Ambient_dimension;
    typedef typename DCTraits::Point_d              Point;
    typedef typename DCTraits::Point_d              Point_d;

    typedef typename Base::Vertex_handle            Vertex_handle;
    typedef typename Base::Vertex_iterator          Vertex_iterator;
    typedef typename Base::Vertex_const_handle      Vertex_const_handle;
    typedef typename Base::Vertex_const_iterator    Vertex_const_iterator;

    typedef typename Base::Simplex_handle           Simplex_handle;
    typedef typename Base::Simplex_iterator         Simplex_iterator;
    typedef typename Base::Simplex_const_handle     Simplex_const_handle;
    typedef typename Base::Simplex_const_iterator   Simplex_const_iterator;

    typedef typename Base::size_type                size_type;
    typedef typename Base::difference_type          difference_type;

    typedef typename Base::Locate_type              Locate_type;

protected: // DATA MEMBERS

    Side_of_oriented_subsphere_d side_of_oss_;

public:
    
    using Base::ambient_dimension;
    using Base::are_incident_simplices_valid;
    using Base::coaffine_orientation_predicate;
    using Base::current_dimension;
    using Base::gather_adjacent_simplices;
    using Base::gather_incident_simplices;
    using Base::geom_traits;
    using Base::get_visited;
    using Base::index_of_covertex;
    using Base::infinite_vertex;
    using Base::insert_in_hole;
    using Base::insert_outside_convex_hull_1;
    using Base::is_finite;
    using Base::is_infinite;
    using Base::is_valid;
    using Base::locate;
    using Base::make_empty_face;
    using Base::new_simplex;
    using Base::number_of_vertices;
    using Base::orientation;
    using Base::pcds;
    using Base::reorient_simplices;
    using Base::set_visited;
    using Base::simplices_begin;
    using Base::simplices_end;
    using Base::vertices_begin;
    using Base::vertices_end;
    // using Base::
    
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - CREATION

    Delaunay_complex(const int dim, const Geom_traits k = Geom_traits())
    : Base(dim, k), side_of_oss_()
    {
    }

    ~Delaunay_complex() {}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ACCESS

    // Not Documented
    Side_of_oriented_subsphere_d & side_of_oriented_subsphere_predicate()
    {
       return side_of_oss_;
    }

    // Not Documented
    const Side_of_oriented_subsphere_d & side_of_oriented_subsphere_predicate() const
    {
        return side_of_oss_;
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - REMOVALS

    Simplex_handle remove(Vertex_handle);
    Simplex_handle remove(const Point & p, Simplex_handle hint = Simplex_handle())
    {
        Vertex_handle v;
        if( is_vertex(p, v, hint) )
            return remove(v);
        return Simplex_handle();
    }

    template< typename ForwardIterator >
    size_type remove(ForwardIterator start, ForwardIterator end)
    {
        size_type n = number_of_vertices();
        while( start != end )
            remove(*start++);
        return ( n - number_of_vertices() );
    }

    // Not documented
    void remove_decrease_dimension(Vertex_handle);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - INSERTIONS

    template< typename ForwardIterator >
    size_type insert(ForwardIterator start, ForwardIterator end)
    {
        size_type n = number_of_vertices();
        std::vector<Point> points(start, end);
        std::random_shuffle(points.begin(), points.end());
        spatial_sort(points.begin(), points.end(), geom_traits());
        Simplex_handle hint;
        for( typename std::vector<Point>::const_iterator p = points.begin(); p != points.end(); ++p )
        {
            hint = insert(*p, hint)->simplex();
        }
        return number_of_vertices() - n;
    }
    Vertex_handle insert(const Point &, const Locate_type, const Face &, const Facet &, const Simplex_handle);
    Vertex_handle insert(const Point & p, const Simplex_handle start = Simplex_handle())
    {
        Locate_type lt;
        Face f(ambient_dimension());
        Facet ft;
        Simplex_handle s = locate(p, lt, f, ft, start);
        return insert(p, lt, f, ft, s);
    }
    Vertex_handle insert(const Point & p, const Vertex_handle hint)
    {
        CGAL_assertion( Vertex_handle() != hint );
        return insert(p, hint->simplex());
    }
    Vertex_handle insert_outside_affine_hull(const Point &);
    Vertex_handle insert_in_conflict_zone(const Point &, const Simplex_handle);

// - - - - - - - - - - - - - - - - - - - - - - - - - GATHERING CONFLICTING SIMPLICES

    bool conflict(const Point &, Simplex_const_handle) const;
    template< class OrientationPredicate >
    Oriented_side perturbed_side_of_positive_sphere(const Point &,
            Simplex_const_handle, const OrientationPredicate &) const;

    template< typename OutputIterator >
    Facet compute_conflict_zone(const Point &, const Simplex_handle, OutputIterator) const;

    template < typename OrientationPredicate, typename SideOfOrientedSpherePredicate >
    class Conflict_predicate
    {
        const Self & dc_;
        const Point & p_;
        const OrientationPredicate & ori_;
        const SideOfOrientedSpherePredicate & side_of_s_;
        int cur_dim_;
    public:
        Conflict_predicate(
                const Self & dc,
                const Point & p,
                const OrientationPredicate & ori,
                const SideOfOrientedSpherePredicate & side)
        : dc_(dc), p_(p), ori_(ori), side_of_s_(side), cur_dim_(dc.current_dimension()) {}
        inline
        bool operator()(Simplex_const_handle s) const
        {
            bool ok;
            if( dc_.is_finite(s) )
            {
                Oriented_side side = side_of_s_(s->points_begin(), s->points_begin() + cur_dim_ + 1, p_);
                if( ON_POSITIVE_SIDE == side )
                    ok = true;
                else if( ON_NEGATIVE_SIDE == side )
                    ok = false;
                else
                    ok = ON_POSITIVE_SIDE == dc_.perturbed_side_of_positive_sphere<OrientationPredicate>(p_, s, ori_);
            }
            else
            {
                s->vertex(0)->set_point(p_); // set position of point at infinity to |p_|
                Orientation o =  ori_(s->points_begin(), s->points_begin() + cur_dim_ + 1);
                if( POSITIVE == o )
                    ok = true;
                else if( o == NEGATIVE )
                    ok = false;
                else
                    ok = (*this)(s->neighbor(0));
            }
            return ok;
        }
    };

    template < typename ConflictPredicate >
    class Conflict_traversal_predicate
    {
        const Self & dc_;
        const ConflictPredicate & pred_;
    public:
        Conflict_traversal_predicate(const Self & dc, const ConflictPredicate & pred)
        : dc_(dc), pred_(pred)
        {}
        inline
        bool operator()(const Facet & f) const
        {
            return pred_(dc_.simplex_of(f)->neighbor(dc_.index_of_covertex(f)));
        }
    };

private:
    // Some internal types to shorten notation
    typedef Conflict_predicate<Coaffine_orientation_d, Side_of_oriented_subsphere_d>
            Conflict_pred_in_subspace;
    typedef Conflict_predicate<Orientation_d, Side_of_oriented_sphere_d>
            Conflict_pred_in_fullspace;
    typedef Conflict_traversal_predicate<Conflict_pred_in_subspace>
            Conflict_traversal_pred_in_subspace;
    typedef Conflict_traversal_predicate<Conflict_pred_in_fullspace>
            Conflict_traversal_pred_in_fullspace;
};

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
// FUNCTIONS THAT ARE MEMBER METHODS:

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - REMOVALS

template< typename DCTraits, typename PCDS >
typename Delaunay_complex<DCTraits, PCDS>::Simplex_handle
Delaunay_complex<DCTraits, PCDS>
::remove( Vertex_handle v )
{
    CGAL_precondition( is_finite(v) );
    CGAL_expensive_precondition( is_vertex(v) );

    // THE CASE cur_dim == 0
    if( 0 == current_dimension() )
    {
        remove_decrease_dimension(v);
        return Simplex_handle();
    }
    else if( 1 == current_dimension() )
    {   // THE CASE cur_dim == 1
        if( 2 == number_of_vertices() )
        {
            remove_decrease_dimension(v);
            return Simplex_handle();
        }
        Simplex_handle left = v->simplex();
        if( is_infinite(left) && left->neighbor(0)->index_of(left) == 0 ) // we are on the infinite right.
            left = left->neighbor(0);
        if( 0 == left->index_of(v) )
            left = left->neighbor(1);
        CGAL_assertion( 1 == left->index_of(v) );
        Simplex_handle right = left->neighbor(0);
        if( is_finite(right) )
        {
            pcds().associate_vertex_with_simplex(left, 1, right->vertex(1));
            set_neighbors(left, 0, right->neighbor(0), right->mirror_index(0));
        }
        else
        {
            pcds().associate_vertex_with_simplex(left, 1, left->vertex(0));
            pcds().associate_vertex_with_simplex(left, 0, infinite_vertex());
            set_neighbors(left, 0, left->neighbor(1), left->mirror_index(1));
            set_neighbors(left, 1, right->neighbor(1), right->mirror_index(1));
        }
        pcds().delete_vertex(v);
        pcds().delete_simplex(right);
        return left;
    }

    // THE CASE cur_dim >= 2
    // Gather the finite vertices sharing an edge with |v|
    typedef std::vector<Simplex_handle> Simplices;
    Simplices simps;
    std::back_insert_iterator<Simplices> out(simps);
    gather_incident_simplices(v, out);
    typedef std::set<Vertex_handle> Vertex_set;
    Vertex_set verts;
    Vertex_handle vh;
    for( typename Simplices::iterator it = simps.begin(); it != simps.end(); ++it )
        for( int i = 0; i <= current_dimension(); ++i )
        {
            vh = (*it)->vertex(i);
            if( is_infinite(vh) )
                continue;
            if( vh == v )
                continue;
            verts.insert(vh);
        }

    // OK, create a Dark Delaunay complex
    typedef Pure_complex_vertex<Geom_traits, Vertex_handle> Dark_vertex_base;
    typedef Pure_complex_simplex<Geom_traits,
        internal::Pure_complex::Dark_simplex_data<Self> > Dark_simplex_base;
    typedef Pure_complex_data_structure<Ambient_dimension, Dark_vertex_base, Dark_simplex_base> Dark_pcds;
    typedef Delaunay_complex<DCTraits, Dark_pcds>   Dark_complex;
    typedef typename Dark_complex::Face             Dark_face;
    typedef typename Dark_complex::Facet            Dark_facet;
    typedef typename Dark_complex::Vertex_handle    Dark_v_handle;
    typedef typename Dark_complex::Simplex_handle   Dark_s_handle;
    Dark_complex dark_side(ambient_dimension());
    Dark_s_handle dark_s;
    Dark_v_handle dark_v;
    typedef std::map<Vertex_handle, Dark_v_handle> Vertex_map;
    Vertex_map light_to_dark;
    typename Vertex_set::iterator vit = verts.begin();
    while( vit != verts.end() )
    {
        dark_v = dark_side.insert((*vit)->point(), dark_s);
        dark_s = dark_v->simplex();
        dark_v->data() = *vit;
        light_to_dark[*vit] = dark_v;
        ++vit;
    }

    if( dark_side.current_dimension() != current_dimension() )
    {
        CGAL_assertion( dark_side.current_dimension() + 1 == current_dimension() );
        if( (size_type)(verts.size() + 1) == number_of_vertices() )
        {
            remove_decrease_dimension(v);
            return Simplex_handle();
        }
        else
        {   // |v| is strictly outside the convex hull of the rest of the points. This is an
            // easy case: first, modify the finite simplices, then, delete the infinite ones.
            // We don't even need the Dark complex.
            // First, mark the infinite simplices
            for( typename Simplices::iterator it = simps.begin(); it != simps.end(); ++it )
            {
                if( is_infinite(*it) )
                    set_visited(*it, true); // mark it for deletion
            }
            // Then, modify the finite simplices:
            for( typename Simplices::iterator it = simps.begin(); it != simps.end(); ++it )
            {
                if( get_visited(*it) )
                        continue;
                int v_idx = (*it)->index_of(v);
                pcds().associate_vertex_with_simplex(*it, v_idx, infinite_vertex());
                if( v_idx != 0 )
                {
                    // we must put the infinite vertex at index 0
                    (*it)->swap_vertices(0, v_idx);
                    // FIXME: are we sure this preseves the positive orientation of the simplex ?
                    (*it)->swap_vertices(current_dimension() - 1, current_dimension());
                }
            }
            // Then, modify the neighboring relation
            for( typename Simplices::iterator it = simps.begin(); it != simps.end(); ++it )
            {
                if( get_visited(*it) )
                    continue;
                for( int i = 1; i <= current_dimension(); ++i )
                {
                    (*it)->vertex(i)->set_simplex(*it);
                    Simplex_handle n = (*it)->neighbor(i);
                    if( ! get_visited(n) )
                        continue;
                    int n_idx = n->index_of(v);
                    set_neighbors(*it, i, n->neighbor(n_idx), n->neighbor(n_idx)->index_of(n));
                }
            }
            Simplex_handle ret_s;
            // Then, we delete the infinite simplices
            for( typename Simplices::iterator it = simps.begin(); it != simps.end(); ++it )
            {
                if( get_visited(*it) )
                    pcds().delete_simplex(*it);
                else
                    ret_s = *it;
            }
            pcds().delete_vertex(v);
            return ret_s;
        }
    }
    else //  From here on, dark_side.current_dimension() == current_dimension()
    {
        dark_side.infinite_vertex()->data() = infinite_vertex();
        light_to_dark[infinite_vertex()] = dark_side.infinite_vertex();
    }

    // Now, compute the conflict zone of v->point() in
    // the dark side. This is precisely the set of simplices
    // that we have to glue back into the light side.
    Dark_face       dark_f = dark_side.make_empty_face();
    Dark_facet      dark_ft;
    typename Dark_complex::Locate_type     lt;
    dark_s = dark_side.locate(v->point(), lt, dark_f, dark_ft);
    CGAL_assertion( lt != Dark_complex::ON_VERTEX
        && lt != Dark_complex::OUTSIDE_AFFINE_HULL );

    // |ret_s| is the simplex that we return
    Dark_s_handle dark_ret_s = dark_s;
    Simplex_handle ret_s;

    typedef std::vector<Dark_s_handle> Dark_simplices;
    Dark_simplices conflict_zone;
    std::back_insert_iterator<Dark_simplices> dark_out(conflict_zone);
    
    dark_ft = dark_side.compute_conflict_zone(v->point(), dark_s, dark_out);

    // THE FOLLOWING SHOULD MAYBE GO IN PCDS.
    // IF SO: make sure to remove set/get_visited from Pure_complex.
    // Here is the plan:
    // 1. Pick any Facet from boundary of the light zone
    // 2. Find corresponding Facet on boundary of dark zone
    // 3. stitch.

    // 1. Build a facet on the boudary of the light zone:
    Simplex_handle light_s = *simps.begin();
    Facet light_ft(light_s, light_s->index_of(v));

    // 2. Find corresponding Dark_facet on boundary of the dark zone
    Dark_simplices dark_incident_s;
    for( int i = 0; i <= current_dimension(); ++i )
    {
        if( index_of_covertex(light_ft) == i )
            continue;
        Dark_v_handle dark_v = light_to_dark[simplex_of(light_ft)->vertex(i)];
        dark_incident_s.clear();
        dark_out = std::back_inserter(dark_incident_s);
        dark_side.gather_incident_simplices(dark_v, dark_out);
        for( typename Dark_simplices::iterator it = dark_incident_s.begin(); it != dark_incident_s.end(); ++it )
        {
            (*it)->data().count_ += 1;
        }
    }

    for( typename Simplices::iterator it = simps.begin(); it != simps.end(); ++it )
        set_visited(*it, true);
    CGAL_assertion( get_visited(simplex_of(light_ft)) );
    for( typename Dark_simplices::iterator it = conflict_zone.begin(); it != conflict_zone.end(); ++it )
        dark_side.set_visited(*it, true);

    bool dark_facet_found(false);
    for( typename Dark_simplices::iterator it = dark_incident_s.begin(); it != dark_incident_s.end(); ++it )
    {
        if( current_dimension() != (*it)->data().count_ )
            continue;
        if( ! dark_side.get_visited(*it) )
            continue;
        // We found a simplex incident to the dark facet corresponding to the light facet |light_ft|
        dark_facet_found = true;
        int ft_idx = 0;
        while( light_s->has_vertex((*it)->vertex(ft_idx)->data()) )
            ++ft_idx;
        dark_ft = Dark_facet(*it, ft_idx);
    }
    // Now, we are ready to traverse both boundary and do the stiching.

    // But first, we create the new simplices in the light complex,
    // with as much adjacency information as possible.

    // Create new simplices with vertices
    for( typename Dark_simplices::iterator it = conflict_zone.begin(); it != conflict_zone.end(); ++it )
    {
        Simplex_handle new_s = new_simplex();
        (*it)->data().light_copy_ = new_s;
        for( int i = 0; i <= current_dimension(); ++i )
            pcds().associate_vertex_with_simplex(new_s, i, (*it)->vertex(i)->data());
        if( dark_ret_s == *it )
            ret_s = new_s;
    }

    // Setup adjacencies inside the hole
    for( typename Dark_simplices::iterator it = conflict_zone.begin(); it != conflict_zone.end(); ++it )
    {
        Simplex_handle new_s = (*it)->data().light_copy_;
        for( int i = 0; i <= current_dimension(); ++i )
            if( dark_side.get_visited((*it)->neighbor(i)) )
                pcds().set_neighbors(new_s, i, (*it)->neighbor(i)->data().light_copy_, (*it)->mirror_index(i));
    }

    // 3. Stitch
    typedef std::queue<std::pair<Facet, Dark_facet> > Queue;
    Queue q;
    q.push(std::make_pair(light_ft, dark_ft));
    dark_s = dark_side.simplex_of(dark_ft);
    int dark_i = dark_side.index_of_covertex(dark_ft);
    // mark dark_ft as visited:
    // TODO try by marking with Dark_v_handle (vertex)
    dark_s->neighbor(dark_i)->set_neighbor(dark_s->mirror_index(dark_i), Dark_s_handle());
    while( ! q.empty() )
    {
        std::pair<Facet, Dark_facet> p = q.front();
        q.pop();
        light_ft = p.first;
        dark_ft = p.second;
        light_s = simplex_of(light_ft);
        int light_i = index_of_covertex(light_ft);
        dark_s = dark_side.simplex_of(dark_ft);
        int dark_i = dark_side.index_of_covertex(dark_ft);
        Simplex_handle light_n = light_s->neighbor(light_i);
        set_neighbors(dark_s->data().light_copy_, dark_i, light_n, light_s->mirror_index(light_i));
        for( int di = 0; di <= current_dimension(); ++di )
        {
            if( di == dark_i )
                continue;
            int li = light_s->index_of(dark_s->vertex(di)->data());
            typename Pure_complex_ds::Rotor light_r(light_s, li, light_i);
            typename Dark_complex::Pure_complex_ds::Rotor dark_r(dark_s, di, dark_i);
            while( ! pcds().is_boundary_facet(light_r) )
                light_r = pcds().rotate_rotor(light_r);
            while( ! dark_side.pcds().is_boundary_facet(dark_r) )
                dark_r = dark_side.pcds().rotate_rotor(dark_r);
            Dark_s_handle dark_ns = dark_side.simplex_of(dark_r);
            int dark_ni = dark_side.index_of_covertex(dark_r);
            Simplex_handle light_ns = simplex_of(light_r);
            int light_ni = index_of_covertex(light_r);
            // mark dark_r as visited:
            // TODO try by marking with Dark_v_handle (vertex)
            Dark_s_handle outside = dark_ns->neighbor(dark_ni);
            Dark_v_handle mirror = dark_ns->mirror_vertex(dark_ni, current_dimension());
            int dn = outside->index_of(mirror);
            if( Dark_s_handle() == outside->neighbor(dn) )
                continue;
            outside->set_neighbor(dn, Dark_s_handle());
            q.push(std::make_pair(Facet(light_ns, light_ni), Dark_facet(dark_ns, dark_ni)));
        }
    }
    pcds().delete_simplices(simps.begin(), simps.end());
    pcds().delete_vertex(v);
    return ret_s;
}

template< typename DCTraits, typename PCDS >
void
Delaunay_complex<DCTraits, PCDS>
::remove_decrease_dimension(Vertex_handle v)
{
    CGAL_precondition( current_dimension() >= 0 );
    pcds().remove_decrease_dimension(v, infinite_vertex());
    // reset the predicates:
    coaffine_orientation_predicate() = geom_traits().coaffine_orientation_d_object();
    side_of_oriented_subsphere_predicate() = geom_traits().side_of_oriented_subsphere_d_object();
    if( 1 <= current_dimension() )
    {
        Simplex_handle s = infinite_vertex()->simplex()->neighbor(0);
        Orientation o = orientation(s);
        CGAL_assertion( COPLANAR != o );
        if( NEGATIVE == o )
            reorient_simplices();
    }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - INSERTIONS

template< typename DCTraits, typename PCDS >
typename Delaunay_complex<DCTraits, PCDS>::Vertex_handle
Delaunay_complex<DCTraits, PCDS>
::insert(const Point & p, const Locate_type lt, const Face & f, const Facet & ft, const Simplex_handle s)
{
    switch( lt )
    {
        case Base::OUTSIDE_AFFINE_HULL:
            return insert_outside_affine_hull(p);
            break;
        case Base::ON_VERTEX:
        {
            Vertex_handle v = s->vertex(f.index(0));
            v->set_point(p);
            return v;
            break;
        }
        default:
            if( 1 == current_dimension() )
            {
                if( Base::OUTSIDE_CONVEX_HULL == lt )
                {
                    return insert_outside_convex_hull_1(p, s);
                }
                Vertex_handle v = pcds().insert_in_simplex(s);
                v->set_point(p);
                return v;
            }
            else
                return insert_in_conflict_zone(p, s);
            break;
    }
}

template< typename DCTraits, typename PCDS >
typename Delaunay_complex<DCTraits, PCDS>::Vertex_handle
Delaunay_complex<DCTraits, PCDS>
::insert_outside_affine_hull(const Point & p)
{
    // we don't use Base::insert_outside_affine_hull(...) because here, we
    // also need to reset the side_of_oriented_subsphere functor.
    CGAL_precondition( current_dimension() < ambient_dimension() );
    Vertex_handle v = pcds().insert_increase_dimension(infinite_vertex());
    // reset the predicates:
    coaffine_orientation_predicate() = geom_traits().coaffine_orientation_d_object();
    side_of_oriented_subsphere_predicate() = geom_traits().side_of_oriented_subsphere_d_object();
    v->set_point(p);
    if( current_dimension() >= 1 )
    {
        Simplex_handle s = infinite_vertex()->simplex()->neighbor(0);
        Orientation o = orientation(s);
        CGAL_assertion( COPLANAR != o );
            if( NEGATIVE == o )
                reorient_simplices();
    }
    return v;
}

template< typename DCTraits, typename PCDS >
typename Delaunay_complex<DCTraits, PCDS>::Vertex_handle
Delaunay_complex<DCTraits, PCDS>
::insert_in_conflict_zone(const Point & p, const Simplex_handle s)
{
    typedef std::vector<Simplex_handle> Simplex_h_vector;
    typedef typename Simplex_h_vector::iterator SHV_iterator;
    static Simplex_h_vector cs; // for storing conflicting simplices.
    cs.clear();
    // cs.reserve(64);
    std::back_insert_iterator<Simplex_h_vector> out(cs);
    Facet ft = compute_conflict_zone(p, s, out);
    return insert_in_hole(p, cs.begin(), cs.end(), ft);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - GATHERING CONFLICTING SIMPLICES

// NOT DOCUMENTED
template< typename DCTraits, typename PCDS >
template< typename OrientationPred >
Oriented_side
Delaunay_complex<DCTraits, PCDS>
::perturbed_side_of_positive_sphere(const Point & p, Simplex_const_handle s,
        const OrientationPred & ori) const
{
    CGAL_precondition_msg( is_finite(s), "simplex must be finite");
    CGAL_expensive_precondition( POSITIVE == orientation(s) );
    typedef std::vector<const Point *> Points;
    Points points(current_dimension() + 2);
    int i(0);
    for( ; i <= current_dimension(); ++i )
        points[i] = &(s->vertex(i)->point());
    points[i] = &p;
    std::sort(points.begin(), points.end(),
            internal::Pure_complex::Compare_points_for_perturbation<Self>(*this));
    typename Points::const_reverse_iterator cut_pt = points.rbegin();
    Points test_points;
    while( cut_pt != points.rend() )
    {
        if( &p == *cut_pt )
            // because the simplex "s" is assumed to be positively oriented
            return ON_NEGATIVE_SIDE; // we consider |p| to lie outside the sphere
        test_points.clear();
        typename Simplex::Point_const_iterator spit = s->points_begin();
        int adjust_sign = -1;
        for( i = 0; i < current_dimension(); ++i )
        {
            if( &(*spit) == *cut_pt )
            {
                ++spit;
                adjust_sign = (((current_dimension() + i) % 2) == 0) ? -1 : +1;
            }
            test_points.push_back(&(*spit));
            ++spit;
        }
        test_points.push_back(&p);

        typedef typename CGAL::Iterator_project<typename Points::iterator,
                        internal::Pure_complex::Point_from_pointer<Self>,
                const Point &, const Point *> Point_pointer_iterator;

        Orientation ori_value = ori(
                Point_pointer_iterator(test_points.begin()),
                Point_pointer_iterator(test_points.end()));

        if( COPLANAR != ori_value )
            return Oriented_side( - adjust_sign * ori_value );

        ++cut_pt;
    }
    CGAL_assertion(false);
    return ON_NEGATIVE_SIDE; // outside circumsphere
}

template< typename DCTraits, typename PCDS >
bool
Delaunay_complex<DCTraits, PCDS>
::conflict(const Point & p, Simplex_const_handle s) const
{
    CGAL_precondition( 2 <= current_dimension() );
    if( current_dimension() < ambient_dimension() )
    {
        Conflict_pred_in_subspace c(*this, p, coaffine_orientation_predicate(), side_of_oriented_subsphere_predicate());
        return c(s);
    }
    else
    {
        Orientation_d ori = geom_traits().orientation_d_object();
        Side_of_oriented_sphere_d side = geom_traits().side_of_oriented_sphere_d_object();
        Conflict_pred_in_fullspace c(*this, p, ori, side);
        return c(s);
    }
}

template< typename DCTraits, typename PCDS >
template< typename OutputIterator >
typename Delaunay_complex<DCTraits, PCDS>::Facet
Delaunay_complex<DCTraits, PCDS>
::compute_conflict_zone(const Point & p, const Simplex_handle s, OutputIterator out) const
{
    CGAL_precondition( 2 <= current_dimension() );
    if( current_dimension() < ambient_dimension() )
    {
        Conflict_pred_in_subspace c(*this, p, coaffine_orientation_predicate(), side_of_oriented_subsphere_predicate());
        Conflict_traversal_pred_in_subspace tp(*this, c);
        return pcds().gather_simplices(s, tp, out);
    }
    else
    {
        Orientation_d ori = geom_traits().orientation_d_object();
        Side_of_oriented_sphere_d side = geom_traits().side_of_oriented_sphere_d_object();
        Conflict_pred_in_fullspace c(*this, p, ori, side);
        Conflict_traversal_pred_in_fullspace tp(*this, c);
        return pcds().gather_simplices(s, tp, out);
    }
}

} //namespace CGAL

#endif // CGAL_DELAUNAY_COMPLEX_H
