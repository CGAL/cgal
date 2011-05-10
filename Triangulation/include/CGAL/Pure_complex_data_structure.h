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

#ifndef CGAL_TRIANGULATION_DATA_STRUCTURE_H
#define CGAL_TRIANGULATION_DATA_STRUCTURE_H

#include <CGAL/basic.h>
#include <CGAL/tuple.h>
#include <CGAL/Default.h>
#include <CGAL/iterator.h>
#include <CGAL/Compact_container.h>
#include <CGAL/Triangulation_face.h>
#include <CGAL/Pure_complex_ds_vertex.h>
#include <CGAL/Pure_complex_ds_simplex.h>
#include <CGAL/internal/Combination_enumerator.h>
#include <CGAL/internal/Triangulation/utilities.h>
#include <CGAL/internal/Triangulation/Triangulation_ds_iterators.h>

#include <algorithm>
#include <vector>
#include <queue>
#include <set>

namespace CGAL {

template<   class Dimen,
            class Vb = Default,
            class Sb = Default >
class Pure_complex_data_structure
{
    typedef Pure_complex_data_structure<Dimen, Vb, Sb>                  Self;
    typedef typename Default::Get<Vb, Pure_complex_ds_vertex<> >::type  V_base;
    typedef typename Default::Get<Sb, Pure_complex_ds_simplex<> >::type S_base;

public:
    typedef typename V_base::template Rebind_TDS<Self>::Other  Vertex;
    typedef typename S_base::template Rebind_TDS<Self>::Other  Simplex;

protected:
    typedef Compact_container<Vertex>   Vertex_container;
    typedef Compact_container<Simplex>  Simplex_container;

public:
    typedef Dimen                      Ambient_dimension;

    typedef typename Vertex_container::size_type        size_type;
    typedef typename Vertex_container::difference_type  difference_type;

    typedef typename Vertex_container::iterator         Vertex_handle;
    typedef typename Vertex_container::iterator         Vertex_iterator;
    typedef typename Vertex_container::const_iterator   Vertex_const_handle;
    typedef typename Vertex_container::const_iterator   Vertex_const_iterator;

    typedef typename Simplex_container::iterator        Simplex_handle;
    typedef typename Simplex_container::iterator        Simplex_iterator;
    typedef typename Simplex_container::const_iterator  Simplex_const_handle;
    typedef typename Simplex_container::const_iterator  Simplex_const_iterator;

    typedef internal::Triangulation::Triangulation_ds_facet_iterator<Self>
                                                        Facet_iterator;

    /* The 2 types defined below, |Facet| and |Rotor| are used when traversing
     the boundary `B' of the union of a set of simplices. |Rotor| makes it
     easy to rotate around itself, in the search of neighbors in `B' (see
     |rotate_rotor| and |insert_in_tagged_hole|) */

    // A co-dimension 1 sub-simplex.
    typedef cpp0x::tuple<Simplex_handle, int>          Facet;
    
    // A co-dimension 2 sub-simplex. called a Rotor because we can rotate
    // the two "covertices" around the sub-simplex. Useful for traversing the
    // boundary of a hole. NOT DOCUMENTED
    typedef cpp0x::tuple<Simplex_handle, int, int>    Rotor;

    typedef Triangulation_face<Self>                 Face;

protected: // DATA MEMBERS

    int dmax_, dcur_; // dimension of the current complex
    Vertex_container  vertices_;  // list of all vertices
    Simplex_container simplices_; // list of all simplices

private:

    void clean_dynamic_memory()
    { 
        vertices_.clear();
        simplices_.clear();
    }

    template < class Dim_tag >
    struct get_ambient_dimension
    {
        static int value(const int D) { return D; }
    };
    // specialization
    template < int D >
    struct get_ambient_dimension<Dimension_tag<D> >
    {
        static int value(const int) { return D; }
    };

public:

    Pure_complex_data_structure(const int dim) 
        : dmax_(get_ambient_dimension<Dimen>::value(dim)), dcur_(-2), vertices_(), simplices_()
    {
        CGAL_assertion_msg(dmax_ > 0, "ambient dimension must be positive.");
    }

    ~Pure_complex_data_structure()
    {
        clean_dynamic_memory();
    }

    // QUERIES

protected:

    bool check_range(const int i) const
    {
        if( current_dimension() < 0 )
        {
            return (0 == i);
        }
        return ( (0 <= i) && (i <= current_dimension()) );
    }

public:

    /* returns the current dimension of the simplices in the complex. */
    int ambient_dimension() const { return dmax_; }
    int current_dimension() const { return dcur_; }

    size_type number_of_vertices() const  { return this->vertices_.size();}
    size_type number_of_simplices() const  { return this->simplices_.size();}

    bool empty() const { return current_dimension() == -2; }

    Vertex_container & vertices() { return vertices_; }
    const Vertex_container & vertices() const { return vertices_; }
    Simplex_container & simplices() { return simplices_; }
    const Simplex_container & simplices() const { return simplices_; }

    Vertex_handle vertex(const Simplex_handle s, const int i) const 
    {
        CGAL_precondition(s != Simplex_handle() && check_range(i));
        return s->vertex(i);
    }

    Vertex_const_handle vertex(const Simplex_const_handle s, const int i) const 
    {
        CGAL_precondition(s != Simplex_handle() && check_range(i));
        return s->vertex(i);
    }

    bool is_vertex(const Vertex_const_handle & v) const
    {
        if( Vertex_const_handle() == v )
            return false;
        Vertex_const_iterator vit = vertices_begin();
        while( vit != vertices_end() && ( v != vit ) )
            ++vit;
        return v == vit;
    }

    bool is_simplex(const Simplex_const_handle & s) const
    {
        if( Simplex_const_handle() == s )
            return false;
        Simplex_const_iterator sit = simplices_begin();
        while( sit != simplices_end() && ( s != sit ) )
            ++sit;
        return s == sit;
    }

    Simplex_handle simplex(const Vertex_handle v) const 
    {
        CGAL_precondition(v != Vertex_handle());
        return v->simplex();
    } 

    Simplex_const_handle simplex(const Vertex_const_handle v) const 
    {
        CGAL_precondition(Vertex_const_handle() != v);
        return v->simplex();
    }

    Simplex_handle neighbor(const Simplex_handle s, const int i) const 
    {
        CGAL_precondition(Simplex_handle() != s && check_range(i));
        return s->neighbor(i);
    }

    Simplex_const_handle neighbor(const Simplex_const_handle s, const int i) const 
    {
        CGAL_precondition(Simplex_const_handle() != s && check_range(i));
        return s->neighbor(i);
    }

    int mirror_index(const Simplex_handle s, const int i) const
    {
        CGAL_precondition(Simplex_handle() != s && check_range(i));
        return s->mirror_index(i);
    }

    int mirror_index(const Simplex_const_handle s, const int i) const
    {
        CGAL_precondition(Simplex_const_handle() != s && check_range(i));
        return s->mirror_index(i);
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - FACETS OPERATIONS

	Face make_empty_face() const
	{
		return Face(ambient_dimension());
	}

    // works for Face_ = Facet and Face_ = Rotor.
    // NOT DOCUMENTED for the Rotor case...
    template< typename Face_ >
    Simplex_handle simplex_of(const Face_ & f) const
    {
        return cpp0x::get<0>(f);
    }

    // works for Face_ = Facet and Face_ = Rotor.
    // NOT DOCUMENTED for the Rotor case...
    template< class Face_ >
    int index_of_covertex(const Face_ & f) const
    {
        return cpp0x::get<1>(f);
    }

    // NOT DOCUMENTED
    // A Rotor has two covertices
    int index_of_second_covertex(const Rotor & f) const
    {
        return cpp0x::get<2>(f);
    }

    // works for Face_ = Facet and Face_ = Rotor.
    // NOT DOCUMENTED...
    template< class Face_ >
    bool is_boundary_facet(const Face_ & f) const
    {
        if( get_visited(neighbor(simplex_of(f), index_of_covertex(f))) )
            return false;
        if( ! get_visited(simplex_of(f)) )
            return false;
        return true;
    }

    // NOT DOCUMENTED...
    Rotor rotate_rotor(Rotor & f)
    {
        int opposite = mirror_index(simplex_of(f), index_of_covertex(f));
        Simplex_handle s = neighbor(simplex_of(f), index_of_covertex(f));
        int new_second = s->index_of(vertex(simplex_of(f), index_of_second_covertex(f)));
        return Rotor(s, new_second, opposite);
    }

    //       NICE UPDATE OPERATIONS

protected:
    void do_insert_increase_dimension(const Vertex_handle, const Vertex_handle);
public:
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - REMOVALS

    Vertex_handle contract_face(const Face &);
    void remove_decrease_dimension(Vertex_handle, Vertex_handle);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - INSERTIONS

    Vertex_handle insert_in_simplex(Simplex_handle);
    Vertex_handle insert_in_face(const Face &);
    Vertex_handle insert_in_facet(const Facet &);
    template< typename Forward_iterator >
    Vertex_handle insert_in_hole(Forward_iterator, const Forward_iterator, Facet);
    template< typename Forward_iterator, typename OutputIterator >
    Vertex_handle insert_in_hole(Forward_iterator, const Forward_iterator, Facet, OutputIterator);

    template< typename OutputIterator >
    Simplex_handle insert_in_tagged_hole(Vertex_handle, Facet, OutputIterator);

	Vertex_handle insert_increase_dimension(Vertex_handle);

    // NOT DOCUMENTED
    void clear_visited_marks(Simplex_handle) const;

    //  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  DANGEROUS UPDATE OPERATIONS

    // NOT DOCUMENTED
    bool get_visited(Simplex_handle s) const
    {
        CGAL_precondition(s != Simplex_handle());
        return static_cast<bool>(s->get_flags() & (unsigned int)1);
    }
    bool get_visited(Simplex_const_handle s) const
    {
        CGAL_precondition(s != Simplex_const_handle());
        return static_cast<bool>(s->get_flags() & (unsigned int)1);
    }

    // NOT DOCUMENTED
    void set_visited(Simplex_handle s, bool b) const
    {
        CGAL_precondition(s != Simplex_handle());
        unsigned int flags = s->get_flags();
        if( b )
            flags = (flags | (unsigned int)1);
        else
            flags = (flags & (~(unsigned int)1));
        s->set_flags(flags);
    }

    void clear()
    {
        clean_dynamic_memory();
        dcur_ = -2; 
    }

    void set_current_dimension(const int d)
    {
        CGAL_precondition(-1<=d && d<=ambient_dimension());
        dcur_ = d;
    }

    Simplex_handle new_simplex(const Simplex_handle s)
    { 
        return simplices_.emplace(*s);
    }

    Simplex_handle new_simplex()
    { 
        return simplices_.emplace(dmax_);
    }

    void delete_simplex(Simplex_handle s)
    {
        CGAL_precondition(Simplex_handle() != s);
        // CGAL_expensive_precondition(is_simplex(s));
        simplices_.erase(s);
    }

    template< typename Forward_iterator >
    void delete_simplices(Forward_iterator start, Forward_iterator end)
    {
        Forward_iterator s = start;
        while( s != end )
            simplices_.erase(*s++);
    }

    template< class T >
    Vertex_handle new_vertex( const T & t )
    { 
        return vertices_.emplace(t);
    }

    Vertex_handle new_vertex()
    { 
        return vertices_.emplace();
    }

    void delete_vertex(Vertex_handle v)
    {
        CGAL_precondition( Vertex_handle() != v );
        vertices_.erase(v);
    }

    void associate_vertex_with_simplex(Simplex_handle s, const int i, Vertex_handle v)
    {
        CGAL_precondition(check_range(i));
        CGAL_precondition(s != Simplex_handle());
        CGAL_precondition(v != Vertex_handle());
        s->set_vertex(i, v);
        v->set_simplex(s);
    }

    void set_neighbors(Simplex_handle s, int i, Simplex_handle s1, int j)
    {
        CGAL_precondition(check_range(i));
        CGAL_precondition(check_range(j));
        CGAL_precondition(s  != Simplex_handle());
        CGAL_precondition(s1 != Simplex_handle());
        s->set_neighbor(i, s1);
        s1->set_neighbor(j, s);
        s->set_mirror_index(i, j);
        s1->set_mirror_index(j, i);
    }

    // SANITY CHECKS

    bool is_valid(bool = true, int = 0) const;
    /*  op Partially checks whether |\Mvar| is an abstract simplicial
        complex. This function terminates without error if each vertex is a
        vertex of the simplex of which it claims to be a vertex, if the
        vertices of all simplices are pairwise distinct, if the neighbor
        relationship is symmetric, and if neighboring simplices share exactly
        |dcur_| vertices.  It returns an error message if one of these
        conditions is violated.  Note that it is not checked whether simplices
        that share |dcur_| vertices are neighbors in the data structure.
    */

    // NOT DOCUMENTED
    template< class OutStream> void write_graph(OutStream &);

    Vertex_iterator vertices_begin() { return vertices_.begin(); }
    Vertex_iterator vertices_end()   { return vertices_.end(); }
    Simplex_iterator simplices_begin() { return simplices_.begin(); }
    Simplex_iterator simplices_end()   { return simplices_.end(); }

    Vertex_const_iterator vertices_begin() const { return vertices_.begin(); }
    Vertex_const_iterator vertices_end()   const { return vertices_.end(); }
    Simplex_const_iterator simplices_begin() const { return simplices_.begin(); }
    Simplex_const_iterator simplices_end()   const { return simplices_.end(); }

    Facet_iterator facets_begin()
    {
        if( current_dimension() <= 0 )
            return facets_end();
        return Facet_iterator(*this);
    }
    Facet_iterator facets_end()
    {
        return Facet_iterator(*this, 0);
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - SIMPLEX GATHERING

    // a traversal predicate for gathering simplices incident to a given face
    // ``incident'' means that the given face is a subface of the simplex
    class Incident_simplex_traversal_predicate
    {
        const Face & f_;
        int dim_;
        const Pure_complex_data_structure & pcds_;
    public:
        Incident_simplex_traversal_predicate(const Pure_complex_data_structure & pcds,
                                            const Face & f)
        : f_(f), pcds_(pcds)
        {
            dim_ = f.feature_dimension();
        }
        bool operator()(const Facet & facet) const
        {
            Vertex_handle v = pcds_.simplex_of(facet)->vertex(pcds_.index_of_covertex(facet));
            for( int i = 0; i <= dim_; ++i )
            {
                if( v == f_.vertex(i) )
                    return false;
            }
            return true;
        }
    };

    // a traversal predicate for gathering simplices adjacent to a given face
    // ``adjacent'' means that the given face shares at least one vertex with the simplex
    class Adjacent_simplex_traversal_predicate
    {
        const Face & f_;
        int dim_;
        const Pure_complex_data_structure & pcds_;
    public:
        Adjacent_simplex_traversal_predicate(const Pure_complex_data_structure & pcds,
                                            const Face & f)
        : f_(f), pcds_(pcds)
        {
            dim_ = f.feature_dimension();
        }
        bool operator()(const Facet & facet) const
        {
            Simplex_handle s = pcds_.simplex_of(facet)->neighbor(pcds_.index_of_covertex(facet));
            for( int j = 0; j <= pcds_.current_dimension(); ++j )
            {
                for( int i = 0; i <= dim_; ++i )
                    if( s->vertex(j) == f_.vertex(i) )
                        return true;
            }
            return false;
        }
    };

    template< typename TraversalPredicate, typename OutputIterator >
    Facet gather_simplices(Simplex_handle, TraversalPredicate &, OutputIterator &) const;
    template< typename OutputIterator >
    OutputIterator gather_incident_simplices(const Face &, OutputIterator) const;
    template< typename OutputIterator >
    OutputIterator gather_incident_simplices(Vertex_const_handle, OutputIterator) const;
    template< typename OutputIterator >
    OutputIterator gather_adjacent_simplices(const Face &, OutputIterator) const;
#ifndef CGAL_CFG_NO_CPP0X_DEFAULT_TEMPLATE_ARGUMENTS_FOR_FUNCTION_TEMPLATES 
    template< typename OutputIterator, typename Comparator = std::less<Vertex_const_handle> >
    OutputIterator gather_incident_upper_faces(Vertex_const_handle v, const int d, OutputIterator out, Comparator cmp = Comparator())
    {
        return gather_incident_faces(v, d, out, cmp, true);
    }
    template< typename OutputIterator, typename Comparator = std::less<Vertex_const_handle> >
    OutputIterator gather_incident_faces(Vertex_const_handle, const int, OutputIterator, Comparator = Comparator(), bool = false);
#else
    template< typename OutputIterator, typename Comparator >
    OutputIterator gather_incident_upper_faces(Vertex_const_handle v, const int d, OutputIterator out, Comparator cmp = Comparator())
    {
        return gather_incident_faces(v, d, out, cmp, true);
    }
    template< typename OutputIterator >
    OutputIterator gather_incident_upper_faces(Vertex_const_handle v, const int d, OutputIterator out)
    {
        return gather_incident_faces(v, d, out, std::less<Vertex_const_handle>(), true);
    }
    template< typename OutputIterator, typename Comparator >
    OutputIterator gather_incident_faces(Vertex_const_handle, const int, OutputIterator, Comparator = Comparator(), bool = false);
    template< typename OutputIterator >
    OutputIterator gather_incident_faces(Vertex_const_handle, const int, OutputIterator,
        std::less<Vertex_const_handle> = std::less<Vertex_const_handle>(), bool = false);
#endif

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - INPUT / OUTPUT

    std::istream & read_simplices(std::istream &, const std::vector<Vertex_handle> &);
    std::ostream & write_simplices(std::ostream &, std::map<Vertex_const_handle, int> &) const;

}; // end of ``declaration/definition'' of Pure_complex_data_structure<...>

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

// FUNCTIONS THAT ARE MEMBER FUNCTIONS:

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
// - - - - - - - - - - - - - - - - - - - - - - - - THE GATHERING METHODS

template< class Dim, class Vb, class Sb >
template< typename OutputIterator >
OutputIterator
Pure_complex_data_structure<Dim, Vb, Sb>
::gather_incident_simplices(const Face & f, OutputIterator out) const
{
    // CGAL_expensive_precondition_msg(is_simplex(f.simplex()), "the facet does not belong to the Pure_complex");
    Incident_simplex_traversal_predicate tp(*this, f);
    gather_simplices(f.simplex(), tp, out);
    return out;
}

template< class Dim, class Vb, class Sb >
template< typename OutputIterator >
OutputIterator
Pure_complex_data_structure<Dim, Vb, Sb>
::gather_incident_simplices(Vertex_const_handle v, OutputIterator out) const
{
//    CGAL_expensive_precondition(is_vertex(v));
    CGAL_precondition(Vertex_handle() != v);
    Face f(v->simplex());
    f.set_index(0, v->simplex()->index_of(v));
    return gather_incident_simplices(f, out);
}

template< class Dim, class Vb, class Sb >
template< typename OutputIterator >
OutputIterator
Pure_complex_data_structure<Dim, Vb, Sb>
::gather_adjacent_simplices(const Face & f, OutputIterator out) const
{
    // CGAL_precondition_msg(is_simplex(f.simplex()), "the facet does not belong to the Pure_complex");
    Adjacent_simplex_traversal_predicate tp(*this, f);
    gather_simplices(f.simplex(), tp, out);
    return out;
}

template< class Dim, class Vb, class Sb >
template< typename TraversalPredicate, typename OutputIterator >
typename Pure_complex_data_structure<Dim, Vb, Sb>::Facet
Pure_complex_data_structure<Dim, Vb, Sb>
::gather_simplices(   Simplex_handle start,
                    TraversalPredicate & tp,
                    OutputIterator & out) const
{
    std::queue<Simplex_handle> queue;
    set_visited(start, true);
    queue.push(start);
    const int cur_dim = current_dimension();
    Facet ft;
    while( ! queue.empty() )
    {
        Simplex_handle s = queue.front();
        queue.pop();
        *out = s;
        ++out;
        for( int i = 0; i <= cur_dim; ++i )
        {
            Simplex_handle n = s->neighbor(i);
            if( ! get_visited(n) )
            {
                set_visited(n, true);
                if( tp(Facet(s, i)) )
                    queue.push(n);
                else
                    ft = Facet(s, i);
            }
        }
    }
    clear_visited_marks(start);
    return ft;
}

#ifdef CGAL_CFG_NO_CPP0X_DEFAULT_TEMPLATE_ARGUMENTS_FOR_FUNCTION_TEMPLATES 
template< class Dim, class Vb, class Sb >
template< typename OutputIterator >
OutputIterator
Pure_complex_data_structure<Dim, Vb, Sb>
::gather_incident_faces(Vertex_const_handle v, const int d, OutputIterator out,
    std::less<Vertex_const_handle> cmp, bool upper_faces)
{
    return gather_incident_faces<OutputIterator, std::less<Vertex_const_handle> >(v, d, out, cmp, upper_faces);
}
#endif

template< class Dim, class Vb, class Sb >
template< typename OutputIterator, typename Comparator >
OutputIterator
Pure_complex_data_structure<Dim, Vb, Sb>
::gather_incident_faces(Vertex_const_handle v, const int d, OutputIterator out, Comparator cmp, bool upper_faces)
{
    CGAL_precondition( 0 < d );
    if( d >= current_dimension() )
        return out;
    typedef std::vector<Simplex_handle> Simplices;
    Simplices simps;
    simps.reserve(64);
    // gather incident simplices
    std::back_insert_iterator<Simplices> sout(simps);
    gather_incident_simplices(v, sout);
    // for storing the handles to the vertices of a simplex
    typedef std::vector<Vertex_const_handle> Vertices;
    typedef std::vector<int> Indices;
    Vertices vertices(1 + current_dimension());
    Indices sorted_idx(1 + current_dimension());
    // setup Face comparator and Face_set
    typedef internal::Triangulation::Compare_faces_with_common_first_vertex<Self>
        Upper_face_comparator;
    Upper_face_comparator ufc(d);
    typedef std::set<Face, Upper_face_comparator> Face_set;
    Face_set face_set(ufc);
    for( typename Simplices::const_iterator s = simps.begin(); s != simps.end(); ++s )
    {
        int v_idx(0); // the index of |v| in the sorted simplex
        // get the vertices of the simplex and sort them
        for( int i = 0; i <= current_dimension(); ++i )
            vertices[i] = (*s)->vertex(i);
        if( upper_faces )
        {
            std::sort(vertices.begin(), vertices.end(), cmp);
            while( vertices[v_idx] != v )
                ++v_idx;
        }
        else
        {
            while( vertices[v_idx] != v )
                ++v_idx;
            if( 0 != v_idx )
                std::swap(vertices[0], vertices[v_idx]);
            v_idx = 0;
            typename Vertices::iterator vbegin(vertices.begin());
            ++vbegin;
            std::sort(vbegin, vertices.end(), cmp);
        }
        if( v_idx + d > current_dimension() )
            continue; // |v| is too far to the right
        // stores the index of the vertices of s in the same order
        // as in |vertices|:
        for( int i = 0; i <= current_dimension(); ++i )
            sorted_idx[i] = (*s)->index_of(vertices[i]);
        // init state for enumerating all candidate faces:
        internal::Combination_enumerator f_idx(d, v_idx + 1, current_dimension());
        Face f(*s);
        f.set_index(0, v_idx);
        while( ! f_idx.end() )
        {
            // check if face has already been found
            for( int i = 0; i < d; ++i )
                f.set_index(1 + i, sorted_idx[f_idx[i]]);
            face_set.insert(f);
            // compute next sorted face (lexicographic enumeration)
            ++f_idx;
        }
    }
    typename Face_set::iterator fit = face_set.begin();
    while( fit != face_set.end() )
        *out++ = *fit++;
    return out;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
// - - - - - - - - - - - - - - - - - - - - - - - - THE REMOVAL METHODS

template <class Dim, class Vb, class Sb>
typename Pure_complex_data_structure<Dim, Vb, Sb>::Vertex_handle
Pure_complex_data_structure<Dim, Vb, Sb>
::contract_face(const Face & f)
{
    const int fd = f.feature_dimension();
    CGAL_precondition( (1 <= fd ) && (fd < current_dimension()));
    std::vector<Simplex_handle> simps;
    // save the Face's vertices:
    Simplex s;
    for( int i = 0; i <= fd; ++i )
        s.set_vertex(i, f.vertex(i));
    // compute adjacent simplices
    simps.reserve(64);
    std::back_insert_iterator<std::vector<Simplex_handle> > out(simps);
    gather_adjacent_simplices(f, out);
    Vertex_handle v = insert_in_hole(simps.begin(), simps.end(), Facet(f.simplex(), f.index(0)));
    for( int i = 0; i <= fd; ++i )
        delete_vertex(s.vertex(i));
    return v;
}

template <class Dim, class Vb, class Sb>
void
Pure_complex_data_structure<Dim, Vb, Sb>
::remove_decrease_dimension(Vertex_handle v, Vertex_handle star)
{
    CGAL_assertion( current_dimension() >= -1 );
    if( -1 == current_dimension() )
    {
        clear();
        return;
    }
    else if( 0 == current_dimension() )
    {
        delete_simplex(v->simplex());
        delete_vertex(v);
        star->simplex()->set_neighbor(0, Simplex_handle());
        set_current_dimension(-1);
        return;
    }
    else if( 1 == current_dimension() )
    {
        Simplex_handle s = v->simplex();
        int star_index;
        if( s->has_vertex(star, star_index) )
            s = s->neighbor(star_index);
        // Here, |s| is not adjacent to |star|, so it's the only finite
        // simplex
        Simplex_handle inf1 = s->neighbor(0);
        Simplex_handle inf2 = s->neighbor(1);
        Vertex_handle v2 = s->vertex(1 - s->index_of(v));
        delete_vertex(v);
        delete_simplex(s);
        inf1->set_vertex(1, Vertex_handle());
        inf1->set_vertex(1, Vertex_handle());
        inf2->set_neighbor(1, Simplex_handle());
        inf2->set_neighbor(1, Simplex_handle());
        associate_vertex_with_simplex(inf1, 0, star);
        associate_vertex_with_simplex(inf2, 0, v2);
        set_neighbors(inf1, 0, inf2, 0);
        set_current_dimension(0);
        return;
    }
    typedef std::vector<Simplex_handle> Simplices;
    Simplices simps;
    gather_incident_simplices(v, std::back_inserter(simps));
    for( typename Simplices::iterator it = simps.begin(); it != simps.end(); ++it )
    {
        int v_idx = (*it)->index_of(v);
        if( ! (*it)->has_vertex(star) )
        {
            delete_simplex((*it)->neighbor(v_idx));
            for( int i = 0; i <= current_dimension(); ++i )
                (*it)->vertex(i)->set_simplex(*it);
        }
        else
            star->set_simplex(*it);
        if( v_idx != current_dimension() )
        {
            (*it)->swap_vertices(v_idx, current_dimension());
            if( ( ! (*it)->has_vertex(star) ) || (current_dimension() > 2) )
                (*it)->swap_vertices(current_dimension() - 2, current_dimension() - 1);
        }
        (*it)->set_vertex(current_dimension(), Vertex_handle());
        (*it)->set_neighbor(current_dimension(), Simplex_handle());
    }
    set_current_dimension(current_dimension()-1);
    delete_vertex(v);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
// - - - - - - - - - - - - - - - - - - - - - - - - THE INSERTION METHODS

template <class Dim, class Vb, class Sb>
typename Pure_complex_data_structure<Dim, Vb, Sb>::Vertex_handle
Pure_complex_data_structure<Dim, Vb, Sb>
::insert_in_simplex(Simplex_handle s)
{
    CGAL_precondition(0 < current_dimension());
    CGAL_precondition(Simplex_handle() != s);
    // CGAL_expensive_precondition(is_simplex(s));

    const int cur_dim = current_dimension();
    Vertex_handle v = new_vertex();
    // the simplex simps is just used to store the handle to all the new simplices.
    Simplex simps(ambient_dimension());
    for( int i = 1; i <= cur_dim; ++i )
    {
        Simplex_handle new_s = new_simplex(s);
        simps.set_neighbor(i, new_s);
        associate_vertex_with_simplex(new_s, i, v);
        s->vertex(i-1)->set_simplex(new_s);
        set_neighbors(new_s, i, neighbor(s, i), mirror_index(s, i));
    }
    simps.set_neighbor(0, s);
    associate_vertex_with_simplex(s, 0, v);
    for( int i = 0; i <= cur_dim; ++i )
        for( int j = 0; j <= cur_dim; ++j )
        {
            if( j == i ) continue;
            set_neighbors(simps.neighbor(i), j, simps.neighbor(j), i);
        }
    return v;
}

template <class Dim, class Vb, class Sb >
typename Pure_complex_data_structure<Dim, Vb, Sb>::Vertex_handle
Pure_complex_data_structure<Dim, Vb, Sb>
::insert_in_face(const Face & f)
{
    std::vector<Simplex_handle> simps;
    simps.reserve(64);
    std::back_insert_iterator<std::vector<Simplex_handle> > out(simps);
    gather_incident_simplices(f, out);
    return insert_in_hole(simps.begin(), simps.end(), Facet(f.simplex(), f.index(0)));
}
template <class Dim, class Vb, class Sb >
typename Pure_complex_data_structure<Dim, Vb, Sb>::Vertex_handle
Pure_complex_data_structure<Dim, Vb, Sb>
::insert_in_facet(const Facet & ft)
{
    Simplex_handle s[2];
    s[0] = simplex_of(ft);
    int i = index_of_covertex(ft);
    s[1] = s[0]->neighbor(i);
    i = ( i + 1 ) % current_dimension();
    return insert_in_hole(s, s+2, Facet(s[0], i));
}

template <class Dim, class Vb, class Sb >
template < typename OutputIterator >
typename Pure_complex_data_structure<Dim, Vb, Sb>::Simplex_handle
Pure_complex_data_structure<Dim, Vb, Sb>
::insert_in_tagged_hole(Vertex_handle v, Facet f,
                        OutputIterator new_simplices)
{
    CGAL_assertion_msg(is_boundary_facet(f), "starting facet should be on the hole boundary");

    const int cur_dim = current_dimension();

    Simplex_handle old_s = simplex_of(f);
    Simplex_handle new_s = new_simplex();
    const int facet_index = index_of_covertex(f);

    int i(0);
    for( ; i < facet_index; ++i )
        associate_vertex_with_simplex(new_s, i, old_s->vertex(i));
    ++i; // skip facet_index
    for( ; i <= cur_dim; ++i )
        associate_vertex_with_simplex(new_s, i, old_s->vertex(i));
    associate_vertex_with_simplex(new_s, facet_index, v);
    set_neighbors(  new_s,
                    facet_index,
                    neighbor(old_s, facet_index),
                    mirror_index(old_s, facet_index));

    // add the new simplex to the list of new simplices
    *new_simplices++ = new_s;

    // check all of |Facet f|'s neighbors
    for( i = 0; i <= cur_dim; ++i )
    {
        if( facet_index == i )
            continue;
        // we define a |Rotor| because it makes it easy to rotate around
        // in a self contained fashion. The corresponding potential
        // boundary facet is Facet(simplex_of(rot), index_of_covertex(rot))
        Rotor rot(old_s, i, facet_index);
        // |rot| on line above, stands for Candidate Facet
        while( ! is_boundary_facet(rot) )
            rot = rotate_rotor(rot);

        // we did find the |i|-th neighbor of Facet(old_s, facet_index)...
        // has it already been extruded to center point |v| ?
        Simplex_handle outside = neighbor(simplex_of(rot), index_of_covertex(rot));
        Simplex_handle inside  = simplex_of(rot);
        Vertex_handle m = inside->mirror_vertex(index_of_covertex(rot), current_dimension());
        int index = outside->index_of(m);
        Simplex_handle new_neighbor = outside->neighbor(index);

        if( new_neighbor == inside )
        {   // not already extruded... we do it recursively
            new_neighbor = insert_in_tagged_hole(   v,
                                                    Facet(simplex_of(rot), index_of_covertex(rot)),
                                                    new_simplices);
        }
        // now the new neighboring simplex exists, we link both
        set_neighbors(new_s, i, new_neighbor, index_of_second_covertex(rot));
    }
    return new_s;
}

template< class Dim, class Vb, class Sb >
template< typename Forward_iterator, typename OutputIterator >
typename Pure_complex_data_structure<Dim, Vb, Sb>::Vertex_handle
Pure_complex_data_structure<Dim, Vb, Sb>
::insert_in_hole(   Forward_iterator start, Forward_iterator end, Facet f,
                    OutputIterator out)
{
    CGAL_expensive_precondition(
            ( std::distance(start, end) == 1 )
         || ( current_dimension() > 1 ) );
    Forward_iterator sit = start;
    while( end != sit )
        set_visited(*sit++, true);
    Vertex_handle v = new_vertex();
    insert_in_tagged_hole(v, f, out);
    delete_simplices(start, end);
    return v;
}

template< class Dim, class Vb, class Sb >
template< typename Forward_iterator >
typename Pure_complex_data_structure<Dim, Vb, Sb>::Vertex_handle
Pure_complex_data_structure<Dim, Vb, Sb>
::insert_in_hole(Forward_iterator start, Forward_iterator end, Facet f)
{
    Emptyset_iterator out;
    return insert_in_hole(start, end, f, out);
}

template <class Dim, class Vb, class Sb>
void
Pure_complex_data_structure<Dim, Vb, Sb>
::clear_visited_marks(Simplex_handle start) const
{
    CGAL_precondition(start != Simplex_handle());

    std::queue<Simplex_handle> queue;
    set_visited(start, false);
    queue.push(start);
    const int cur_dim = current_dimension();
    while( ! queue.empty() )
    {
        Simplex_handle s = queue.front();
        queue.pop();
        for( int i = 0; i <= cur_dim; ++i )
        {
            if( get_visited(s->neighbor(i)) )
            {
                set_visited(s->neighbor(i), false);
                queue.push(s->neighbor(i));
            }
        }
    }
}

template <class Dim, class Vb, class Sb>
void Pure_complex_data_structure<Dim, Vb, Sb>
::do_insert_increase_dimension(const Vertex_handle x, const Vertex_handle star)
{
    Simplex_handle start = simplices_begin();
    Simplex_handle swap_me;
    const int cur_dim = current_dimension();
    for( Simplex_iterator S = simplices_begin(); S != simplices_end(); ++S )
    {
        if( Vertex_handle() != S->vertex(cur_dim) )
            continue;
        set_visited(S, true);
        // extends simplex |S| to include the new vertex as the
        // current_dimension()-th vertex
        associate_vertex_with_simplex(S, cur_dim, x);
        if( ! S->has_vertex(star) )
        {   // S is bounded, we create its unbounded "twin" simplex
            Simplex_handle S_new = new_simplex();
            set_neighbors(S, cur_dim, S_new, 0);
            associate_vertex_with_simplex(S_new, 0, star);
            // here, we could be clever so as to get consistent orientation
            for( int k = 1; k <= cur_dim; ++k )
                associate_vertex_with_simplex(S_new, k, vertex(S, k - 1));
        }
        else if( cur_dim == 2 )
        {   // if cur. dim. is 2, we must take care of the 'rightmost' infinite vertex.
            if( S->mirror_index(S->index_of(star)) == 0 )
                swap_me = S;
        }
    }
    // now we setup the neighbors
    set_visited(start, false);
    std::queue<Simplex_handle> queue;
    queue.push(start);
    while( ! queue.empty() )
    {
        Simplex_handle S = queue.front();
        queue.pop();
        // here, the first visit above ensured that all neighbors exist now.
        // Now we need to connect them with adjacency relation
        int star_index;
        if( S->has_vertex(star, star_index) )
        {
            set_neighbors(  S, cur_dim, neighbor(neighbor(S, star_index), cur_dim),
                            // this is tricky :-)  :
                            mirror_index(S, star_index) + 1);
        }
        else
        {
            Simplex_handle S_new = neighbor(S, cur_dim);
            for( int k = 0 ; k < cur_dim ; ++k )
            {
                Simplex_handle S_opp = neighbor(S, k);
                if( ! S_opp->has_vertex(star) )   
                    set_neighbors(S_new, k + 1, neighbor(S_opp, cur_dim), mirror_index(S, k) + 1);
                    // neighbor of S_new opposite to v is S_new'
                    // the vertex opposite to v remains the same but ...
                    // remember the shifting of the vertices one step to the right
            }
        }
        for( int k = 0 ; k < cur_dim ; ++k )
            if( get_visited(neighbor(S, k)) )
            {
                set_visited(neighbor(S, k), false);
                queue.push(neighbor(S, k));
            }
    }
    if( ( ( cur_dim % 2 ) == 0 ) && ( cur_dim > 1 ) )
    {
        for( Simplex_iterator S = simplices_begin(); S != simplices_end(); ++S )
        {
            if( x != S->vertex(cur_dim) )
                S->swap_vertices(cur_dim - 1, cur_dim);
        }
    }
    if( Simplex_handle() != swap_me )
        swap_me->swap_vertices(1, 2);
}

template <class Dim, class Vb, class Sb>
typename Pure_complex_data_structure<Dim, Vb, Sb>::Vertex_handle
Pure_complex_data_structure<Dim, Vb, Sb>
::insert_increase_dimension(Vertex_handle star = Vertex_handle())
{
    const int prev_cur_dim = current_dimension();
    CGAL_precondition(prev_cur_dim < ambient_dimension());
    if( -2 != current_dimension() )
    {
        CGAL_precondition( Vertex_handle() != star );
        CGAL_expensive_precondition(is_vertex(star));
    }
    else
    {
        CGAL_precondition( Vertex_handle() == star );
    }

    set_current_dimension(prev_cur_dim + 1);
    Vertex_handle v = new_vertex();
    switch( prev_cur_dim )
    {
        case -2:
        {   // insertion of the first vertex
            // ( geometrically : infinite vertex )
            Simplex_handle s = new_simplex();
            associate_vertex_with_simplex(s, 0, v);
            break;
        }
        case -1:
        {   // insertion of the second vertex
            // ( geometrically : first finite vertex )
            //we create a triangulation of the 0-sphere, with
            // vertices |star| and |v|
            Simplex_handle infinite_simplex = star->simplex();
            Simplex_handle finite_simplex = new_simplex();
            associate_vertex_with_simplex(finite_simplex, 0, v);
            set_neighbors(infinite_simplex, 0, finite_simplex, 0);
            break;
        }
        default:
            do_insert_increase_dimension(v, star);
            break;
    }
    return v;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
// - - - - - - - - - - - - - - - - - - - - - - - - VALIDITY CHECKS

template <class Dimen, class Vb, class Sb>
bool Pure_complex_data_structure<Dimen, Vb, Sb>
::is_valid(bool verbose, int /* level */) const
{ 
    Simplex_const_handle s, t;
    Vertex_const_handle v;
    int i, j, k;

    if( dcur_ == -2 )
    {
        if( ! vertices_.empty() || ! simplices_.empty() )
        {
            if( verbose ) CGAL_warning_msg(false, "current dimension is -2 but there are vertices or simplices");
            return false;
        }
    }

    if( dcur_ == -1 )
    {
        if ( (number_of_vertices() != 1) || (number_of_simplices() != 1) )
        {
            if( verbose ) CGAL_warning_msg(false, "current dimension is -1 but there isn't one vertex and one simplex");
            return false;
        }
    }

    int fake_dcur = (dcur_ > 0) ? dcur_ : 0;
    for( v = vertices_begin(); v != vertices_end(); ++v )
    {
        if( ! v->is_valid(verbose) )
            return false;
        bool ok(false);
        // check that |v|'s simplex actually contains |v|
        for( i = 0; i <= fake_dcur; ++i )
        {
            if( v->simplex()->vertex(i) == v )
            {
                ok = true;
                break;
            }
        }
        if( ! ok )
        {
            if( verbose ) CGAL_warning_msg(false, "the simplex incident to some vertex does not contain that vertex.");
            return false;
        }
    }
    // FUTURE: for each vertex v, gather incident simplices. then, check that
    // any simplex containing v is among those gathered simplices...

    if( dcur_ < 0 )
        return true;

    for( s = simplices_begin(); s != simplices_end(); ++s )
    {
        if( ! s->is_valid(verbose) )
            return false;
        for( i = 0; i <= dcur_; ++i )
            for( j = i + 1; j <= dcur_; ++j )
                if( vertex(s,i) == vertex(s,j) )
                {
                    CGAL_warning_msg(false, "a simplex has two equal vertices");
                    return false;
                }
    }

    for( s = simplices_begin(); s != simplices_end(); ++s )
    {
        for( i = 0; i <= dcur_; ++i )
            if( (t = neighbor(s,i)) != Simplex_const_handle() )
            {
                int l = mirror_index(s,i);
                if( s != neighbor(t,l) || i != mirror_index(t,l) )
                {
                    if( verbose ) CGAL_warning_msg(false, "neighbor relation is not symmetric");
                    return false;
                }
                for( j = 0; j <= dcur_; ++j )
                    if( j != i )
                    {
                        // j must also occur as a vertex of t
                        for( k = 0; k <= dcur_ && ( vertex(s,j) != vertex(t,k) || k == l); ++k )
                            ;
                        if( k > dcur_ )
                        {
                            if( verbose ) CGAL_warning_msg(false, "too few shared vertices between neighbors simplices.");
                            return false;
                        }
                    }
            }
            else
            {
                if( verbose ) CGAL_warning_msg(false, "simplex has a NULL neighbor");
                return false;
            }
    }
    return true;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
// - - - - - - - - - - - - - - - - - - - - - - - - INPUT / OUTPUT

// NOT DOCUMENTED
template <class Dim, class Vb, class Sb>
template <class OutStream>
void Pure_complex_data_structure<Dim, Vb, Sb>
::write_graph(OutStream & os)
{
    std::vector<std::set<int> > edges;
    os << number_of_vertices() + 1; // add the vertex at infinity
    int count(1);
    for( Vertex_iterator vit = vertices_begin(); vit != vertices_end(); ++vit )
        vit->idx_ = count++;
    edges.resize(number_of_vertices()+1);
    for( Simplex_iterator sit = simplices_begin(); sit != simplices_end(); ++sit )
    {
        int v1 = 0;
        while( v1 < current_dimension() )
        {
            int v2 = v1 + 1;
            while( v2 <= current_dimension() )
            {
                int i1, i2;
                if( Vertex_handle() != sit-> vertex(v1) )
                    i1 = sit->vertex(v1)->idx_;
                else
                    i1 = 0;
                if( Vertex_handle() != sit-> vertex(v2) )
                    i2 = sit->vertex(v2)->idx_;
                else
                    i2 = 0;
                edges[i1].insert(i2);
                edges[i2].insert(i1);
                ++v2;
            }
            ++v1;
        }
    }
    for( int i = 0; i < edges.size(); ++i )
    {
        os << std::endl << edges[i].size();
        for( std::set<int>::const_iterator nit = edges[i].begin();
        nit !=  edges[i].end(); ++nit )
        {
            os << ' ' << (*nit);
        }
    }
}

// NOT DOCUMENTED...
template<class Dimen, class Vb, class Sb>
std::istream &
Pure_complex_data_structure<Dimen, Vb, Sb>
::read_simplices(std::istream & is, const std::vector<Vertex_handle> & vertices)
{
    size_t m; // number of simplices
    int index;
    const int cd = current_dimension();
    if( is_ascii(is) )
        is >> m;
    else
        read(is, m, io_Read_write());

    std::vector<Simplex_handle> simplices;
    simplices.reserve(m);
    // read the vertices of each simplex
    size_t i = 0;
    while( i < m )
    {
        Simplex_handle s = new_simplex();
        simplices.push_back(s);
        for( int j = 0; j <= cd; ++j )
        {
            if( is_ascii(is) )
                is >> index;
            else
                read(is, index);
            s->set_vertex(j, vertices[index]);
        }
        // read other non-combinatorial information for the simplices
        is >> (*s);
        ++i;
    }

    // read the neighbors of each simplex
    i = 0;
    if( is_ascii(is) )
        while( i < m )
    {
        for( int j = 0; j <= cd; ++j )
        {
            is >> index;
            simplices[i]->set_neighbor(j, simplices[index]);
        }
        ++i;
    }
    else
        while( i < m )
    {
        for( int j = 0; j <= cd; ++j )
        {
            read(is, index);
            simplices[i]->set_neighbor(j, simplices[index]);
        }
        ++i;
    }

    // compute the mirror indices
    for( i = 0; i < m; ++i )
    {
        Simplex_handle s = simplices[i];
        for( int j = 0; j <= cd; ++j )
        {
            if( -1 != s->mirror_index(j) )
                continue;
            Simplex_handle n = s->neighbor(j);
            int k = 0;
            Simplex_handle nn = n->neighbor(k);
            while( s != nn )
                nn = n->neighbor(++k);
            s->set_mirror_index(j,k);
            n->set_mirror_index(k,j);
        }
    }
    return is;
}

// NOT DOCUMENTED...
template<class Dimen, class Vb, class Sb>
std::ostream &
Pure_complex_data_structure<Dimen, Vb, Sb>
::write_simplices(std::ostream & os, std::map<Vertex_const_handle, int> & index_of_vertex) const
{
    std::map<Simplex_const_handle, int> index_of_simplex;

    size_t m = number_of_simplices();

    if( is_ascii(os) )
        os << std::endl << m;
    else
        write(os, m, io_Read_write());

    const int cur_dim = current_dimension();
    // write the vertex indices of each simplex
    size_t i = 0;
    for( Simplex_const_iterator it = simplices_begin(); it != simplices_end(); ++it )
    {
        index_of_simplex[it] = i++;
        if( is_ascii(os) )
            os << std::endl;
        for( int j = 0; j <= cur_dim; ++j )
        {
            if( is_ascii(os) )
                os << ' ' << index_of_vertex[it->vertex(j)];
            else
                write(os, index_of_vertex[it->vertex(j)]);
        }
        // write other non-combinatorial information for the simplices
        os << (*it);
    }

    CGAL_assertion( i == m );

    // write the neighbors of each simplex
    if( is_ascii(os) )
        for( Simplex_const_iterator it = simplices_begin(); it != simplices_end(); ++it )
        {
            os << std::endl;
            for( int j = 0; j <= cur_dim; ++j )
                os << ' ' << index_of_simplex[it->neighbor(j)];
        }
    else
        for( Simplex_const_iterator it = simplices_begin(); it != simplices_end(); ++it )
        {
            for( int j = 0; j <= cur_dim; ++j )
                write(os, index_of_simplex[it->neighbor(j)]);
        }

    return os;
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

// FUNCTIONS THAT ARE NOT MEMBER FUNCTIONS:

template<class Dimen, class Vb, class Sb>
std::istream & 
operator>>(std::istream & is, Pure_complex_data_structure<Dimen, Vb, Sb> & tr)
  // reads :
  // - the dimensions (ambient and current)
  // - the number of finite vertices
  // - the non combinatorial information on vertices (point, etc)
  // - the number of simplices
  // - the simplices by the indices of their vertices in the preceding list
  // of vertices, plus the non combinatorial information on each simplex
  // - the neighbors of each simplex by their index in the preceding list
{
    typedef Pure_complex_data_structure<Dimen, Vb, Sb> PC;
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

    CGAL_assertion_msg( cd <= tr.ambient_dimension(), "input Pure_complex_data_structure has too high dimension");

    tr.clear();
    tr.set_current_dimension(cd);

    if( n == 0 )
        return is;

    std::vector<Vertex_handle> vertices;
    vertices.resize(n);

    // read the vertices:
    size_t i(0);
    while( i < n )
    {
        vertices[i] = tr.new_vertex();
        is >> (*vertices[i]); // read a vertex
        ++i;
    }

    // now, read the combinatorial information
    return tr.read_simplices(is, vertices);
}

template<class Dimen, class Vb, class Sb>
std::ostream & 
operator<<(std::ostream & os, const Pure_complex_data_structure<Dimen, Vb, Sb> & tr)
  // writes :
  // - the dimensions (ambient and current)
  // - the number of finite vertices
  // - the non combinatorial information on vertices (point, etc)
  // - the number of simplices
  // - the simplices by the indices of their vertices in the preceding list
  // of vertices, plus the non combinatorial information on each simplex
  // - the neighbors of each simplex by their index in the preceding list
{
    typedef Pure_complex_data_structure<Dimen, Vb, Sb> PC;
    typedef typename PC::Vertex_const_handle         Vertex_handle;
    typedef typename PC::Vertex_const_iterator       Vertex_iterator;
    typedef typename PC::Simplex_const_handle        Simplex_handle;
    typedef typename PC::Simplex_const_iterator      Simplex_iterator;

    // outputs dimension and number of vertices
    size_t n = tr.number_of_vertices();
    if( is_ascii(os) )
        os << tr.current_dimension() << std::endl << n;
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

    for( Vertex_iterator it = tr.vertices_begin(); it != tr.vertices_end(); ++it, ++i )
    {
        os << *it; // write the vertex
        index_of_vertex[it] = i;
    }
    CGAL_assertion( i == n );

    // output the combinatorial information
    return tr.write_simplices(os, index_of_vertex);
}

} //namespace CGAL

#endif // CGAL_TRIANGULATION_DATA_STRUCTURE_H
