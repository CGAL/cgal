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

#ifndef CGAL_TRIANGULATION_H
#define CGAL_TRIANGULATION_H

#include <boost/iterator/filter_iterator.hpp>
#include <CGAL/internal/Triangulation/utilities.h>
#include <CGAL/Pure_complex_data_structure.h>
#include <CGAL/Pure_complex_simplex.h>
#include <CGAL/Pure_complex_vertex.h>
#include <CGAL/spatial_sort.h>
#include <CGAL/Dimension.h>
#include <CGAL/iterator.h>
#include <CGAL/Default.h>

namespace CGAL {

template <  class PCTraits, class TDS_ = Default >
class Pure_complex
{
    typedef typename Ambient_dimension<typename PCTraits::Point_d>::type
                                                    Ambient_dimension_;
    typedef typename Default::Get<TDS_, Pure_complex_data_structure
                    <   Ambient_dimension_,
                        Pure_complex_vertex<PCTraits>,
                        Pure_complex_simplex<PCTraits> >
                        >::type                     TDS;
    typedef Pure_complex<PCTraits, TDS_>           Self;

    typedef typename PCTraits::Coaffine_orientation_d
                                                    Coaffine_orientation_d;
    typedef typename PCTraits::Orientation_d        Orientation_d;

public:

    typedef PCTraits                                Geom_traits;
    typedef TDS                                    Pure_complex_ds;

    typedef typename TDS::Vertex                   Vertex;
    typedef typename TDS::Simplex                  Simplex;
    typedef typename TDS::Facet                    Facet;
    typedef typename TDS::Face                     Face;

    typedef Ambient_dimension_                      Ambient_dimension;
    typedef typename Geom_traits::Point_d           Point;
    typedef typename Geom_traits::Point_d           Point_d;

    typedef typename TDS::Vertex_handle            Vertex_handle;
    typedef typename TDS::Vertex_iterator          Vertex_iterator;
    typedef typename TDS::Vertex_const_handle      Vertex_const_handle;
    typedef typename TDS::Vertex_const_iterator    Vertex_const_iterator;

    typedef typename TDS::Simplex_handle           Simplex_handle;
    typedef typename TDS::Simplex_iterator         Simplex_iterator;
    typedef typename TDS::Simplex_const_handle     Simplex_const_handle;
    typedef typename TDS::Simplex_const_iterator   Simplex_const_iterator;
    
    typedef typename TDS::Facet_iterator           Facet_iterator;

    typedef typename TDS::size_type                size_type;
    typedef typename TDS::difference_type          difference_type;

    /// The type of location a new point is found lying on
    enum  Locate_type
    {
          ON_VERTEX  = 0
        , IN_FACE    = 1 // simplex of dimension in [1, |current_dimension()| - 2]
        , IN_FACET   = 2 // simplex of dimension |current_dimension()| - 1
        , IN_SIMPLEX = 3
        , OUTSIDE_CONVEX_HULL = 4
        , OUTSIDE_AFFINE_HULL = 5
    };

    // Finite elements iterators

    class Finiteness_predicate;

    typedef boost::filter_iterator<Finiteness_predicate, Vertex_iterator>
        Finite_vertex_iterator;
    typedef boost::filter_iterator<Finiteness_predicate, Vertex_const_iterator>
        Finite_vertex_const_iterator;
    typedef boost::filter_iterator<Finiteness_predicate, Simplex_iterator>
        Finite_simplex_iterator;
    typedef boost::filter_iterator<Finiteness_predicate, Simplex_const_iterator>
        Finite_simplex_const_iterator;
    typedef boost::filter_iterator<Finiteness_predicate, Facet_iterator>
        Finite_facet_iterator;
    
protected: // DATA MEMBERS

    Pure_complex_ds                     pcds_;
    const Geom_traits                   kernel_;
    Vertex_handle                       infinity_;
    mutable std::vector<Oriented_side>  orientations_;
    Coaffine_orientation_d              coaffine_orientation_;
    // for stochastic walk in the locate() function:
    mutable Random                      rng_;
#ifdef CGAL_TRIANGULATION_STATISTICS
    mutable unsigned long walk_size_;
#endif

public:

    //       FACETS OPERATIONS

	Face make_empty_face() const
	{
		return pcds().make_empty_face();
	}

    // works for Face_ = Facet and Face_ = Rotor.
    // NOT DOCUMENTED for the Rotor case...
    template< typename Face_ >
    Simplex_handle simplex_of(const Face_ & f) const
    {
        return pcds().simplex_of(f);
    }

    // works for Face_ = Facet and Face_ = Rotor.
    // NOT DOCUMENTED for the Rotor case...
    template< class Face_ >
    int index_of_covertex(const Face_ & f) const
    {
        return pcds().index_of_covertex<Face_>(f);
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - CREATION

    Pure_complex(const int dim, const Geom_traits k = Geom_traits())
        : pcds_(dim)
        , kernel_(k)
        , infinity_()
        , rng_((long)0)
#ifdef CGAL_TRIANGULATION_STATISTICS
        ,walk_size_(0)
#endif
    {
        clear();
    }

    ~Pure_complex() {}

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ACCESS FUNCTIONS

    // NOT DOCUMENTED -- TRY TO REMOVE IF POSSIBLE (CF Delaunay_complex::remove)
    bool get_visited(Simplex_handle s) const
    {
        return pcds().get_visited(s);
    }
    // NOT DOCUMENTED -- TRY TO REMOVE IF POSSIBLE (CF Delaunay_complex::remove)
    bool get_visited(Simplex_const_handle s) const
    {
        return pcds().get_visited(s);
    }

    // NOT DOCUMENTED -- TRY TO REMOVE IF POSSIBLE (CF Delaunay_complex::remove)
    void set_visited(Simplex_handle s, bool b) const
    {
        pcds().set_visited(s, b);
    }

    Coaffine_orientation_d & coaffine_orientation_predicate()
    {
        return coaffine_orientation_;
    }

    const Coaffine_orientation_d & coaffine_orientation_predicate() const
    {
        return coaffine_orientation_;
    }

    const Pure_complex_ds & pcds() const
    {
        return pcds_;
    }

    Pure_complex_ds & pcds()
    {
        return pcds_;
    }

    const Geom_traits & geom_traits() const
    {
        return kernel_;
    }

    int ambient_dimension() const { return pcds().ambient_dimension(); }
    int current_dimension() const { return pcds().current_dimension(); }

    bool empty() const
    {
        return current_dimension() == -1;
    }

    size_type number_of_vertices() const
    {
        return pcds().number_of_vertices() - 1;
    }

    size_type number_of_simplices() const
    {
        return pcds().number_of_simplices();
    }

    Vertex_handle infinite_vertex() const
    {
        return infinity_;
    }

    Simplex_handle infinite_simplex() const
    {
        CGAL_assertion(infinite_vertex()->simplex()->has_vertex(infinite_vertex()));
        return infinite_vertex()->simplex();
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - NON CONSTANT-TIME ACCESS FUNCTIONS

    size_type number_of_finite_simplices() const
    {
        Simplex_const_iterator s = simplices_begin();
        size_type result = number_of_simplices();
        for( ; s != simplices_end(); ++s )
        {
            if( is_infinite(s) )
                --result;
        }
        return result;
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - TRAVERSAL

    Vertex_iterator vertices_begin() { return pcds().vertices_begin(); }
    Vertex_iterator vertices_end()   { return pcds().vertices_end(); }

    Vertex_const_iterator vertices_begin() const { return pcds().vertices_begin(); }
    Vertex_const_iterator vertices_end()   const { return pcds().vertices_end(); }

    Finite_vertex_iterator finite_vertices_begin()
    { return Finite_vertex_iterator(Finiteness_predicate(*this), vertices_begin(), vertices_end()); }
    Finite_vertex_iterator finite_vertices_end()
    { return Finite_vertex_iterator(Finiteness_predicate(*this), vertices_end(), vertices_end()); }
    Finite_vertex_const_iterator finite_vertices_begin() const
    { return Finite_vertex_const_iterator(Finiteness_predicate(*this), vertices_begin(), vertices_end()); }
    Finite_vertex_const_iterator finite_vertices_end() const
    { return Finite_vertex_const_iterator(Finiteness_predicate(*this), vertices_end(), vertices_end()); }

    Simplex_iterator simplices_begin() { return pcds().simplices_begin(); }
    Simplex_iterator simplices_end()   { return pcds().simplices_end(); }

    Simplex_const_iterator simplices_begin() const { return pcds().simplices_begin(); }
    Simplex_const_iterator simplices_end()   const { return pcds().simplices_end(); }

    Finite_simplex_iterator finite_simplices_begin()
    { return Finite_simplex_iterator(Finiteness_predicate(*this), simplices_begin(), simplices_end()); }
    Finite_simplex_iterator finite_simplices_end()
    { return Finite_simplex_iterator(Finiteness_predicate(*this), simplices_end(), simplices_end()); }
    Finite_simplex_const_iterator finite_simplices_begin() const
    { return Finite_simplex_const_iterator(Finiteness_predicate(*this), simplices_begin(), simplices_end()); }
    Finite_simplex_const_iterator finite_simplices_end() const
    { return Finite_simplex_const_iterator(Finiteness_predicate(*this), simplices_end(), simplices_end()); }

    Facet_iterator facets_begin() { return pcds().facets_begin(); }
    Facet_iterator facets_end() { return pcds().facets_end(); }
    Facet_iterator finite_facets_begin()
    { return Finite_facet_iterator(Finiteness_predicate(*this), facets_begin(), facets_end()); }
    Facet_iterator finite_facets_end()
    { return Finite_facet_iterator(Finiteness_predicate(*this), facets_end(), facets_end()); }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - SOME PREDICATE FUNCTORS

    class Finiteness_predicate
    {
        const Self & pc_;
    public:
        Finiteness_predicate(const Self & pc) : pc_(pc) {}
        template < class T >
        bool operator()(const T & t)
        {
            return ( pc_.is_finite(t.vertex(0)) );
        }
    };

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - SIMPLE QUERIES

    bool is_vertex(const Point & p, Vertex_handle & v, Simplex_handle hint = Simplex_handle()) const
    {
        Locate_type lt;
        Face f(ambient_dimension());
        Facet ft;
        Simplex_handle s = locate(p, lt, f, ft, hint);
        if( ON_VERTEX == lt )
        {
            v = s->vertex(f.index(0));
            return true;
        }
        return false;
    }

    bool is_vertex(Vertex_const_handle v) const
    {
        return pcds().is_vertex(v);
    }

    bool is_simplex(Simplex_const_handle s) const
    {
        return pcds().is_simplex(s);
    }

    bool is_infinite(Vertex_const_handle v) const
    {
        CGAL_precondition(Vertex_const_handle() != v);
        return (infinite_vertex() == v);
    }

    bool is_infinite(Simplex_const_handle s) const
    {
        CGAL_precondition(Simplex_const_handle() != s);
        return is_infinite(s->vertex(0));
    }

    bool is_infinite(const Facet & ft) const
    {
        Simplex_const_handle s = simplex_of(ft);
        CGAL_precondition(s != Simplex_handle());
        if( is_infinite(s) )
            return (index_of_covertex(ft) != 0);
        return false;
    }

    bool is_infinite(const Face & f) const
    {
        Simplex_const_handle s = f.simplex();
        CGAL_precondition(s != Simplex_handle());
        if( is_infinite(s) )
        {
            int i(0);
            Vertex_handle v;
            while( v = f.vertex(i), Vertex_handle() != v )
            {
                if( is_infinite(v) )
                    return true;
                ++i;
            }
            return false;
        }
        return false;
    }

    bool is_finite(Vertex_const_handle v) const
    {
        return (! is_infinite(v));
    }

    bool is_finite(Simplex_const_handle s) const
    {
        return (! is_infinite(s));
    }

    bool is_finite(const Facet & ft) const
    {
        return (! is_infinite(ft));
    }

    bool is_finite(const Face & f) const
    {
        return (! is_infinite(f));
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ELEMENT GATHERING

    template< typename OutputIterator >
    OutputIterator gather_incident_simplices(const Face & f, OutputIterator out) const
    {
        return pcds().gather_incident_simplices(f, out);
    }
    template< typename OutputIterator >
    OutputIterator gather_incident_simplices(Vertex_const_handle v, OutputIterator out) const
    {
        return pcds().gather_incident_simplices(v, out);
    }
    template< typename OutputIterator >
    OutputIterator gather_adjacent_simplices(const Face & f, OutputIterator out) const
    {
        return pcds().gather_incident_simplices(f, out);
    }

    template< typename OutputIterator >
    OutputIterator gather_incident_faces(Vertex_const_handle v, const int d, OutputIterator out)
    {
        return pcds().gather_incident_faces(v, d, out);
    }

    template< typename OutputIterator, class Comparator >
    OutputIterator gather_incident_upper_faces( Vertex_const_handle v, const int d,
                                                OutputIterator out, Comparator cmp = Comparator())
    {
        return pcds().gather_incident_upper_faces(v, d, out, cmp);
    }
    template< typename OutputIterator >
    OutputIterator gather_incident_upper_faces( Vertex_const_handle v, const int d,
                                                OutputIterator out)
    {
        internal::Triangulation::Compare_vertices_for_upper_face<Self> cmp(*this);
        return pcds().gather_incident_upper_faces(v, d, out, cmp);
    }

    Orientation orientation(Simplex_const_handle s, bool in_is_valid = false) const
    {
        if( ! in_is_valid )
            CGAL_assertion( is_finite(s) );
        if( 0 == current_dimension() )
            return POSITIVE;
        if( current_dimension() == ambient_dimension() )
        {
            Orientation_d ori = geom_traits().orientation_d_object();
            return ori(s->points_begin(), s->points_begin() + 1 + current_dimension());
        }
        else
        {
            return coaffine_orientation_predicate()(s->points_begin(), s->points_begin() + 1 + current_dimension());
        }
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - UPDATE OPERATIONS

    void clear()
    {
        pcds_.clear();
        infinity_ = pcds().insert_increase_dimension();
        // A simplex has at most 1 + ambient_dimension() faces:
        orientations_.resize(1 + ambient_dimension());
        // Our coaffine orientation predicates HAS state member variables
        coaffine_orientation_predicate() = geom_traits().coaffine_orientation_d_object();
#ifdef CGAL_TRIANGULATION_STATISTICS
        walk_size_ = 0;
#endif
    }

    void set_current_dimension(const int d)
    {
        pcds().set_current_dimension(d);
    }

    Simplex_handle new_simplex()
    { 
        return pcds().new_simplex();
    }

    Vertex_handle  new_vertex(const Point & p) 
    {
        return pcds().new_vertex(p);
    }

    void set_neighbors(Simplex_handle s, int i, Simplex_handle s1, int j)
    {
        pcds().set_neighbors(s, i, s1, j);
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - VALIDITY

    bool is_valid(bool = true, int = 0) const;
    bool are_incident_simplices_valid(Vertex_const_handle, bool = true, int = 0) const;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - POINT LOCATION

protected:
    template< typename OrientationPredicate >
    Simplex_handle do_locate(   const Point &, Locate_type &, Face &, Facet &,
                                Simplex_handle start = Simplex_handle(),
                                OrientationPredicate & o = geom_traits().orientation_d_object()) const;
public:
    Simplex_handle locate(  const Point &, Locate_type &, Face &, Facet &,
                            Simplex_handle start = Simplex_handle()) const;
    Simplex_handle locate(  const Point &, Locate_type &, Face &, Facet &,
                            Vertex_handle) const;
    Simplex_handle locate(const Point & p, Simplex_handle s = Simplex_handle()) const;
    Simplex_handle locate(const Point & p, Vertex_handle v) const;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - REMOVALS

    Vertex_handle contract_face(const Point &, const Face &);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - POINT INSERTION

    template< typename ForwardIterator >
    size_type insert(ForwardIterator start, ForwardIterator end)
    {
        size_type n = number_of_vertices();
        std::vector<Point> points(start, end);
        std::random_shuffle(points.begin(), points.end());
        spatial_sort(points.begin(), points.end(), geom_traits());
        Simplex_handle hint = Simplex_handle();
        typename std::vector<Point>::const_iterator s = points.begin();
        while( s != points.end() )
        {
            hint = insert(*s++, hint)->simplex();
        }
        return number_of_vertices() - n;
    }
    Vertex_handle insert(const Point &, const Locate_type, const Face &, const Facet &, const Simplex_handle);
    Vertex_handle insert(const Point &, Simplex_handle start = Simplex_handle());
    Vertex_handle insert(const Point &, Vertex_handle);
    template< typename ForwardIterator >
    Vertex_handle insert_in_hole(const Point & p, ForwardIterator start, ForwardIterator end, const Facet & ft)
    {
        Emptyset_iterator out;
        return insert_in_hole(p, start, end, ft, out);
    }
    template< typename ForwardIterator, typename OutputIterator >
    Vertex_handle insert_in_hole(const Point & p, ForwardIterator start, ForwardIterator end, const Facet & ft, 
                                 OutputIterator out)
    {
        Vertex_handle v = pcds().insert_in_hole(start, end, ft, out);
        v->set_point(p);
        return v;
    }
    Vertex_handle insert_in_face(const Point &, const Face &);
    Vertex_handle insert_in_facet(const Point &, const Facet &);
    Vertex_handle insert_in_simplex(const Point &, Simplex_handle);
    Vertex_handle insert_outside_convex_hull_1(const Point &, Simplex_handle);
    Vertex_handle insert_outside_convex_hull(const Point &, Simplex_handle);
    Vertex_handle insert_outside_affine_hull(const Point &);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - FACET-TRAVERSAL PREDICATES

    template< typename OrientationPredicate >
    class Outside_convex_hull_traversal_predicate
    {
        Pure_complex & pc_;
        const Point & p_;
        OrientationPredicate & ori_;
        int cur_dim_;
    public:
        Outside_convex_hull_traversal_predicate(Pure_complex & pc, const Point & p,
                                    OrientationPredicate & ori)
        : pc_(pc), p_(p), ori_(ori), cur_dim_(pc.current_dimension()) {}
        // FUTURE change parameter to const reference
        bool operator()(Facet f) const
        {
            Simplex_handle s = pc_.simplex_of(f);
            const int i = pc_.index_of_covertex(f);
            Simplex_handle n = s->neighbor(i);
            if( pc_.is_finite(n) )
                return false;
            n->vertex(0)->set_point(p_);
            bool ok = (POSITIVE == ori_(n->points_begin(), n->points_begin() + cur_dim_ + 1));
            return ok;
        }
    };
    // make sure all simplices have positive orientation
    void reorient_simplices();

}; // Pure_complex<...>

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

// CLASS MEMBER FUNCTIONS

template < class PCT, class TDS >
void
Pure_complex<PCT, TDS>
::reorient_simplices()
{
    if( current_dimension() < 1 )
        return;
    Simplex_iterator sit = simplices_begin();
    Simplex_iterator send = simplices_end();
    while( sit != send )
    {
        if( is_infinite(sit) && (1 == current_dimension()) )
        {
            ++sit;
            continue;
        }
        sit->swap_vertices(current_dimension() - 1, current_dimension());
        ++sit;
    }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
// - - - - - - - - - - - - - - - - - - - - - - - - THE REMOVAL METHODS

template < class PCT, class TDS >
typename Pure_complex<PCT, TDS>::Vertex_handle
Pure_complex<PCT, TDS>
::contract_face(const Point & p, const Face & f)
{
    CGAL_precondition( is_finite(f) );
    Vertex_handle v = pcds().contract_face(f);
    v->set_point(p);
    CGAL_expensive_postcondition_msg(are_incident_simplices_valid(v), "new point is not where it should be");
    return v;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
// - - - - - - - - - - - - - - - - - - - - - - - - THE INSERTION METHODS

template < class PCT, class TDS >
typename Pure_complex<PCT, TDS>::Vertex_handle
Pure_complex<PCT, TDS>
::insert(const Point & p, const Locate_type lt, const Face & f, const Facet & ft, const Simplex_handle s)
{
    switch( lt )
    {
        case IN_SIMPLEX:
            //std::cerr << " IS";
            return insert_in_simplex(p, s);
            break;
        case OUTSIDE_CONVEX_HULL:
            //std::cerr << " OCH";
            return insert_outside_convex_hull(p, s);
            break;
        case OUTSIDE_AFFINE_HULL:
            //std::cerr << " OAF";
            return insert_outside_affine_hull(p);
            break;
        case IN_FACET:
        {
            //std::cerr << " IFT";
            return insert_in_facet(p, ft);
            break;
        }
        case IN_FACE:
            //std::cerr << " IF" << f.feature_dimension();
            return insert_in_face(p, f);
            break;
        case ON_VERTEX:
            //std::cerr << " OV";
            s->vertex(f.index(0))->set_point(p);
            return s->vertex(f.index(0));
            break;
    }
    CGAL_assertion(false);
    return Vertex_handle();
}

template < class PCT, class TDS >
typename Pure_complex<PCT, TDS>::Vertex_handle
Pure_complex<PCT, TDS>
::insert(const Point & p, Simplex_handle start)
{
    Locate_type lt;
    Face f(ambient_dimension());
    Facet ft;
    Simplex_handle s = locate(p, lt, f, ft, start);
    return insert(p, lt, f, ft, s);
}

template < class PCT, class TDS >
typename Pure_complex<PCT, TDS>::Vertex_handle
Pure_complex<PCT, TDS>
::insert(const Point & p, Vertex_handle v)
{
    if( Vertex_handle() == v )
        v = infinite_vertex();
    return insert(p, v->simplex());
}

template < class PCT, class TDS >
typename Pure_complex<PCT, TDS>::Vertex_handle
Pure_complex<PCT, TDS>
::insert_in_face(const Point & p, const Face & f)
{
    CGAL_precondition( ! is_infinite(f) );
    Vertex_handle v = pcds().insert_in_face(f);
    v->set_point(p);
    return v;
}

template < class PCT, class TDS >
typename Pure_complex<PCT, TDS>::Vertex_handle
Pure_complex<PCT, TDS>
::insert_in_facet(const Point & p, const Facet & ft)
{
    CGAL_precondition( ! is_infinite(ft) );
    Vertex_handle v = pcds().insert_in_facet(ft);
    v->set_point(p);
    return v;
}

template < class PCT, class TDS >
typename Pure_complex<PCT, TDS>::Vertex_handle
Pure_complex<PCT, TDS>
::insert_in_simplex(const Point & p, Simplex_handle s)
{
    CGAL_precondition( is_finite(s) );
    Vertex_handle v = pcds().insert_in_simplex(s);
    v->set_point(p);
    return v;
}

// NOT DOCUMENTED...
template < class PCT, class TDS >
typename Pure_complex<PCT, TDS>::Vertex_handle
Pure_complex<PCT, TDS>
::insert_outside_convex_hull_1(const Point & p, Simplex_handle s)
{
    // This is a special case for dimension 1, because in that case, the right
    // infinite simplex is not correctly oriented... (sice its first vertex is the
    // infinite one...
    CGAL_precondition( is_infinite(s) );
    CGAL_precondition( 1 == current_dimension() );
    bool swap = (0 == s->neighbor(0)->index_of(s));
    Vertex_handle v = pcds().insert_in_simplex(s);
    v->set_point(p);
    if( swap )
    {
        s->swap_vertices(0, 1);
    }
    return v;
}

template < class PCT, class TDS >
typename Pure_complex<PCT, TDS>::Vertex_handle
Pure_complex<PCT, TDS>
::insert_outside_convex_hull(const Point & p, Simplex_handle s)
{
    if( 1 == current_dimension() )
    {
        return insert_outside_convex_hull_1(p, s);
    }
    CGAL_precondition( is_infinite(s) );
    CGAL_assertion( current_dimension() >= 2 );
    std::vector<Simplex_handle> simps;
    simps.reserve(64);
    std::back_insert_iterator<std::vector<Simplex_handle> > out(simps);
    if( current_dimension() < ambient_dimension() )
    {
        Outside_convex_hull_traversal_predicate<Coaffine_orientation_d>
            ochtp(*this, p, coaffine_orientation_predicate());
        pcds().gather_simplices(s, ochtp, out);
    }
    else
    {
        Orientation_d ori = geom_traits().orientation_d_object();
        Outside_convex_hull_traversal_predicate<Orientation_d>
            ochtp(*this, p, ori);
        pcds().gather_simplices(s, ochtp, out);
    }
    Vertex_handle v = insert_in_hole(p, simps.begin(), simps.end(), Facet(s, 0));
    return v;
}

template < class PCT, class TDS >
typename Pure_complex<PCT, TDS>::Vertex_handle
Pure_complex<PCT, TDS>
::insert_outside_affine_hull(const Point & p)
{
    CGAL_precondition( current_dimension() < ambient_dimension() );
    Vertex_handle v = pcds().insert_increase_dimension(infinite_vertex());
    // reset the orientation predicate:
    coaffine_orientation_predicate() = geom_traits().coaffine_orientation_d_object();
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

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
// - - - - - - - - - - - - - - - - - - - - THE MAIN LOCATE(...) FUNCTION

template < class PCT, class TDS >
template< typename OrientationPredicate >
typename Pure_complex<PCT, TDS>::Simplex_handle
Pure_complex<PCT, TDS>
::do_locate(   const Point & p, // query point
            Locate_type & loc_type,// type of result (simplex, face, vertex)
            Face & face,// the face containing the query in its interior (when appropriate)
            Facet & facet,// the facet containing the query in its interior (when appropriate)
            const Simplex_handle start// starting simplex for the walk
            , OrientationPredicate & orientation_pred
        ) const
{
    const int cur_dim = current_dimension();

    if( cur_dim == -1 )
    {
        loc_type = OUTSIDE_AFFINE_HULL;
        return Simplex_handle();
    }
    else if( cur_dim == 0 )
    {
        Vertex_handle vit = infinite_simplex()->neighbor(0)->vertex(0);
        if( EQUAL != geom_traits().compare_lexicographically_d_object()(p, vit->point()) )
        {
            loc_type = OUTSIDE_AFFINE_HULL;
        }
        else
        {
            loc_type = ON_VERTEX;
            face.set_simplex(vit->simplex());
            face.set_index(0, 0);
        }
        return vit->simplex();
    }

    Simplex_handle s;

    // if we don't know where to start, we start from any bounded simplex
    if( Simplex_handle() == start )
        // THE HACK THAT NOBODY SHOULD DO... BUT DIFFICULT TO WORK AROUND
        // THIS... TODO: WORK AROUND IT
        s = const_cast<Self*>(this)->infinite_simplex()->neighbor(0);
    else
    {
        s = start;
        if( is_infinite(s) )
            s = s->neighbor(0);
    }

    // Check if query |p| is outside the affine hull
    if( cur_dim < ambient_dimension() )
    {
        if( ! geom_traits().contained_in_affine_hull_d_object()(
            s->points_begin(),
            s->points_begin() + current_dimension() + 1,
            p) )
        {
            loc_type = OUTSIDE_AFFINE_HULL;
            return Simplex_handle();
        }
    }

    // we remember the |previous|ly visited simplex to avoid the evaluation
    // of one |orientation| predicate
    Simplex_handle previous = Simplex_handle();
    bool simplex_not_found = true;
    while(simplex_not_found) // we walk until we locate the query point |p|
    {
    #ifdef CGAL_TRIANGULATION_STATISTICS
        ++walk_size_;
    #endif
        // For the remembering stochastic walk, we need to start trying
        // with a random index:
        int j, i = rng_.get_int(0, cur_dim);
        // we check |p| against all the simplex's hyperplanes in turn
       
        for(j = 0; j <= cur_dim; ++j, i = (i + 1) % (cur_dim + 1) )
        {
            Simplex_handle next = s->neighbor(i);
            if( previous == next )
            {   // no need to compute the orientation, we already know it
                orientations_[i] = POSITIVE;
                continue; // go to next simplex's facet
            }

            // we temporarily substitute |p| to the |i|-th point of the
            // simplex
            Point backup = s->vertex(i)->point();
            s->vertex(i)->set_point(p);

            orientations_[i] = orientation_pred(
                s->points_begin(),
                s->points_begin() + cur_dim + 1);

            // restore the correct point for vertex |i| of the simplex
            s->vertex(i)->set_point(backup);

            if( orientations_[i] != NEGATIVE )
            {
                // from this facet's point of view, we are inside the
                // simplex or on its boundary, so we continue to next facet
                continue;
            }

            // At this point, we know that we have to jump to the |next|
            // simplex because orientation_[i] == NEGATIVE
            previous = s;
            s = next;
            if( is_infinite(next) )
            {   // we have arrived OUTSIDE the convex hull of the triangulation,
                // so we stop the search
                simplex_not_found = false;
                loc_type = OUTSIDE_CONVEX_HULL;
                face.set_simplex(s);
            }
            break;
        } // end of the 'for' loop
        if( ( cur_dim + 1 ) == j ) // we found the simplex containing |p|
            simplex_not_found = false;
    }
    // Here, we know in which simplex |p| is in.
    // We now check more precisely where |p| landed:
    // vertex, facet, face or simplex.
    if( ! is_infinite(s) )
    {
        face.set_simplex(s);
        int num(0);
        int verts(0);
        for(int i = 0; i <= cur_dim; ++i)
        {
            if( orientations_[i] == COPLANAR )
            {
                ++num;
                facet = Facet(s, i);
            }
            else
                face.set_index(verts++, i);
        }
        if( 0 == num )
        {
            loc_type = IN_SIMPLEX;
            face.clear();
        }
        else if( cur_dim == num )
            loc_type = ON_VERTEX;
        else if( 1 == num )
            loc_type = IN_FACET;
        else
            loc_type = IN_FACE;
            
#if 0 // if ! defined( CGAL_NDEBUG )
        // informative test that can be removed:
        if( loc_type == IN_FACE || loc_type == IN_FACET )
        {
            std::cerr << "\n[Degenerate input : loc_type = " << loc_type << ':';
            for(int i = 0; i <= cur_dim; ++i)
                std::cerr << ' ' << orientations_[i];
            std::cerr << ']';
        }
#endif
    }
    return s;
}

template < class PCT, class TDS >
typename Pure_complex<PCT, TDS>::Simplex_handle
Pure_complex<PCT, TDS>
::locate(   const Point & p, // query point
            Locate_type & loc_type,// type of result (simplex, face, vertex)
            Face & face,// the face containing the query in its interior (when appropriate)
            Facet & facet,// the facet containing the query in its interior (when appropriate)
            Simplex_handle start// starting simplex for the walk
        ) const
{
    if( current_dimension() == ambient_dimension() )
    {
        Orientation_d ori = geom_traits().orientation_d_object();
        return do_locate(p, loc_type, face, facet, start, ori);
    }
    else
        return do_locate(p, loc_type, face, facet, start, coaffine_orientation_predicate());
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
// - - - - - - - - - - - - - - - - - - - - the locate(...) variants

template < class PCT, class TDS >
typename Pure_complex<PCT, TDS>::Simplex_handle
Pure_complex<PCT, TDS>
::locate(   const Point & p,
            Locate_type & loc_type,
            Face & face,
            Facet & facet,
            Vertex_handle start) const
{
    if( Vertex_handle() == start )
        start = infinite_vertex();
    return locate(p, loc_type, face, facet, start->simplex());
}

template < class PCT, class TDS >
typename Pure_complex<PCT, TDS>::Simplex_handle
Pure_complex<PCT, TDS>
::locate(const Point & p, Simplex_handle s) const
{
    Locate_type lt;
    Face face(ambient_dimension());
    Facet facet;
    return locate(p, lt, face, facet, s);
}

template < class PCT, class TDS >
typename Pure_complex<PCT, TDS>::Simplex_handle
Pure_complex<PCT, TDS>
::locate(const Point & p, Vertex_handle v) const
{
    if( Vertex_handle() != v )
        v = infinite_vertex();
    return this->locate(p, v->simplex());
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - VALIDITY

template < class PCT, class TDS >
bool
Pure_complex<PCT, TDS>
::is_valid(bool verbose, int level) const
{ 
    if( ! pcds().is_valid(verbose, level) )
        return false;

    Simplex_const_iterator s;
    for( s = simplices_begin(); s != simplices_end(); ++s )
    {
        if( s->has_vertex(infinite_vertex()) )
        {
            if( s->vertex(0) != infinite_vertex() )
            {
                if( verbose ) CGAL_warning_msg(false, "the vertex at infinity must have index 0");
                return false;
            }
        }
    }
    if( current_dimension() < 0 )
        return true;
    Orientation o;
    for( s = simplices_begin(); s != simplices_end(); ++s )
    {
        if( is_infinite(s) )
        {
            if( current_dimension() > 1 )
            {
                Simplex_handle fs = s->neighbor(0);
                infinite_vertex()->set_point(fs->vertex(s->mirror_index(0))->point());
                o = - orientation(s, true);
            }
            else
                o = POSITIVE;
        }
        else
            o = orientation(s, true);
        if( NEGATIVE == o )
        {
            if( verbose ) CGAL_warning_msg(false, "simplex is not correctly oriented");
            return false;
        }
        if( COPLANAR == o )
        {
            if( verbose ) CGAL_warning_msg(false, "simplex is flat");
            return false;
        }
    }
    return true;
}

template < class PCT, class TDS >
bool Pure_complex<PCT, TDS>::are_incident_simplices_valid(Vertex_const_handle v, bool verbose, int level) const
{
    if( current_dimension() <= 0 )
        return true;
    typedef std::vector<Simplex_const_handle> Simps;
    Simps simps;
    simps.reserve(64);
    std::back_insert_iterator<Simps> out(simps);
    gather_incident_simplices(v, out);
    typename Simps::const_iterator sit = simps.begin();
    for( ; sit != simps.end(); ++sit )
    {
        if( is_infinite(*sit) )
            continue;
        Orientation o = orientation(*sit);
        if( NEGATIVE == o )
        {
            if( verbose ) CGAL_warning_msg(false, "simplex is not correctly oriented");
            return false;
        }
        if( COPLANAR == o )
        {
            if( verbose ) CGAL_warning_msg(false, "simplex is flat");
            return false;
        }
    }
    return true;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// FUNCTIONS THAT ARE NOT MEMBER FUNCTIONS:

template < class PCT, class TDS >
std::istream & 
operator>>(std::istream & is, Pure_complex<PCT, TDS> & tr)
  // reads :
  // - the dimensions (ambient and current)
  // - the number of finite vertices
  // - the non combinatorial information on vertices (point, etc)
  // - the number of simplices
  // - the simplices by the indices of their vertices in the preceding list
  // of vertices, plus the non combinatorial information on each simplex
  // - the neighbors of each simplex by their index in the preceding list
{
    typedef Pure_complex<PCT, TDS> PC;
    typedef typename PC::Vertex_handle         Vertex_handle;
    typedef typename PC::Vertex_iterator       Vertex_iterator;
    typedef typename PC::Simplex_handle        Simplex_handle;
    typedef typename PC::Simplex_iterator      Simplex_iterator;

    // read current dimension and number of vertices
    size_t n;
    int cd;
    if( is_ascii(is) )
        is >> cd >> n;
    else
    {
        read(is, cd);
        read(is, n, io_Read_write());
    }

    CGAL_assertion_msg( cd <= tr.ambient_dimension(), "input Pure_complex has too high dimension");

    tr.clear();
    tr.set_current_dimension(cd);

    if( n == 0 )
        return is;

    std::vector<Vertex_handle> vertices;
    vertices.resize(n+1);
    vertices[0] = tr.infinite_vertex();
    is >> (*vertices[0]);

    // read the vertices:
    size_t i(1);
    while( i <= n )
    {
        vertices[i] = tr.new_vertex();
        is >> (*vertices[i]); // read a vertex
        ++i;
    }

    // now, read the combinatorial information
   return tr.pcds().read_simplices(is, vertices);
}

template < class PCT, class TDS >
std::ostream & 
operator<<(std::ostream & os, const Pure_complex<PCT, TDS> & tr)
  // writes :
  // - the dimensions (ambient and current)
  // - the number of finite vertices
  // - the non combinatorial information on vertices (point, etc)
  // - the number of simplices
  // - the simplices by the indices of their vertices in the preceding list
  // of vertices, plus the non combinatorial information on each simplex
  // - the neighbors of each simplex by their index in the preceding list
{
    typedef Pure_complex<PCT, TDS> PC;
    typedef typename PC::Vertex_const_handle         Vertex_handle;
    typedef typename PC::Vertex_const_iterator       Vertex_iterator;
    typedef typename PC::Simplex_const_handle        Simplex_handle;
    typedef typename PC::Simplex_const_iterator      Simplex_iterator;

    // outputs dimensions and number of vertices
    size_t n = tr.number_of_vertices();
    if( is_ascii(os) )
        os << tr.current_dimension() << std::endl << n << std::endl;
    else
    {
        write(os, tr.current_dimension());
        write(os, n, io_Read_write());
    }

    if( n == 0 )
        return os;

    size_t i(0);
    // write the vertices
    std::map<Vertex_handle, int> index_of_vertex;

    // infinite vertex has index 0
    index_of_vertex[tr.infinite_vertex()] = i++;
    os << *tr.infinite_vertex();
    for( Vertex_iterator it = tr.vertices_begin(); it != tr.vertices_end(); ++it )
    {
        if( tr.is_infinite(it) )
            continue;
        os << *it; // write the vertex
        index_of_vertex[it] = i++;
    }
    CGAL_assertion( i == n+1 );

    // output the combinatorial information
    return tr.pcds().write_simplices(os, index_of_vertex);
}

} //namespace CGAL

#endif // CGAL_TRIANGULATION_H
