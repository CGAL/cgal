// Copyright (c) 2009-2014 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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

#include <CGAL/license/Triangulation.h>


#include <CGAL/basic.h>
#include <CGAL/Default.h>
#include <CGAL/iterator.h>
#include <CGAL/Compact_container.h>
#include <CGAL/Triangulation_face.h>
#include <CGAL/Triangulation_ds_vertex.h>
#include <CGAL/Triangulation_ds_full_cell.h>
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
            class Fcb = Default >
class Triangulation_data_structure
{
    typedef Triangulation_data_structure<Dimen, Vb, Fcb>                     Self;
    typedef typename Default::Get<Vb, Triangulation_ds_vertex<> >::type     V_base;
    typedef typename Default::Get<Fcb, Triangulation_ds_full_cell<> >::type  FC_base;

public:
    typedef typename V_base::template Rebind_TDS<Self>::Other   Vertex; /* Concept */
    typedef typename FC_base::template Rebind_TDS<Self>::Other  Full_cell; /* Concept */

  // Tools to change the Vertex and Cell types of the TDS.
  template < typename Vb2 >
  struct Rebind_vertex {
    typedef Triangulation_data_structure<Dimen, Vb2, Fcb>  Other;
  };

  template < typename Fcb2 >
  struct Rebind_full_cell {
    typedef Triangulation_data_structure<Dimen, Vb, Fcb2>  Other;
  };



    // we want to store an object of this class in every Full_cell:
    class Full_cell_data
    {
        unsigned char bits_;
        public:
        Full_cell_data() : bits_(0) {}
        Full_cell_data(const Full_cell_data & fcd) : bits_(fcd.bits_) {}

        void clear()         { bits_ = 0; }
        void mark_visited()  { bits_ = 1; }
        void clear_visited() { bits_ = 0; }

        bool is_clear()   const { return bits_ == 0; }
        bool is_visited() const { return bits_ == 1; }
        // WARNING: if we use more bits and several bits can be set at once,
        // then make sure to use bitwise operation above, instead of direct
        // affectation.
    };

protected:
    typedef Compact_container<Vertex>   Vertex_container;
    typedef Compact_container<Full_cell>  Full_cell_container;

public:
    typedef Dimen                      Maximal_dimension;

    typedef typename Vertex_container::size_type        size_type; /* Concept */
    typedef typename Vertex_container::difference_type  difference_type; /* Concept */

    typedef typename Vertex_container::iterator         Vertex_handle; /* Concept */
    typedef typename Vertex_container::iterator         Vertex_iterator; /* Concept */
    typedef typename Vertex_container::const_iterator   Vertex_const_handle;
    typedef typename Vertex_container::const_iterator   Vertex_const_iterator;

    typedef typename Full_cell_container::iterator        Full_cell_handle; /* Concept */
    typedef typename Full_cell_container::iterator        Full_cell_iterator; /* Concept */
    typedef typename Full_cell_container::const_iterator  Full_cell_const_handle;
    typedef typename Full_cell_container::const_iterator  Full_cell_const_iterator;

    typedef internal::Triangulation::
            Triangulation_ds_facet_iterator<Self>         Facet_iterator; /* Concept */

    /* The 2 types defined below, |Facet| and |Rotor| are used when traversing
     the boundary `B' of the union of a set of full cells. |Rotor| makes it
     easy to rotate around itself, in the search of neighbors in `B' (see
     |rotate_rotor| and |insert_in_tagged_hole|) */

    // A co-dimension 1 sub-simplex.
    class Facet /* Concept */
    {
        Full_cell_handle full_cell_;
        int index_of_covertex_;
    public:
        Facet() : full_cell_(), index_of_covertex_(0) {}
        Facet(Full_cell_handle f, int i) : full_cell_(f), index_of_covertex_(i) {}
        Full_cell_handle full_cell() const { return full_cell_; }
        int index_of_covertex() const { return index_of_covertex_; }
    };

    // A co-dimension 2 sub-simplex. called a Rotor because we can rotate
    // the two "covertices" around the sub-simplex. Useful for traversing the
    // boundary of a hole. NOT DOCUMENTED
    class Rotor : public Facet
    {
        int index_of_second_covertex_;
    public:
        Rotor() : Facet(), index_of_second_covertex_(0) {}
        Rotor(Full_cell_handle f, int first, int second) : Facet(f, first), index_of_second_covertex_(second) {}
        int index_of_second_covertex() const { return index_of_second_covertex_; }
    };

    typedef Triangulation_face<Self>                    Face; /* Concept */

protected: // DATA MEMBERS

    int dmax_, dcur_; // dimension of the current triangulation
    Vertex_container  vertices_;  // list of all vertices
    Full_cell_container full_cells_; // list of all full cells

private:

    void clean_dynamic_memory()
    {
        vertices_.clear();
        full_cells_.clear();
    }

    template < class Dim_tag >
    struct get_maximal_dimension
    {
        static int value(int D) { return D; }
    };
    // specialization
    template < int D >
    struct get_maximal_dimension<Dimension_tag<D> >
    {
        static int value(int) { return D; }
    };

public:
    Triangulation_data_structure( int dim=0)  /* Concept */
        : dmax_(get_maximal_dimension<Dimen>::value(dim)), dcur_(-2), 
          vertices_(), full_cells_()
    {
        CGAL_assertion_msg(dmax_ > 0, "maximal dimension must be positive.");
    }
  
    ~Triangulation_data_structure()
    {
        clean_dynamic_memory();
    }

    Triangulation_data_structure(const Triangulation_data_structure & tds)
        : dmax_(tds.dmax_), dcur_(tds.dcur_),
        vertices_(tds.vertices_), full_cells_(tds.full_cells_)
    {
        typedef std::map<Vertex_const_handle, Vertex_handle> V_map;
        typedef std::map<Full_cell_const_handle, Full_cell_handle> C_map;
        V_map vmap;
        C_map cmap;
        Vertex_const_iterator vfrom = tds.vertices_begin();
        Vertex_iterator vto = vertices_begin();
        Full_cell_const_iterator cfrom = tds.full_cells_begin();
        Full_cell_iterator cto = full_cells_begin();
        while( vfrom != tds.vertices_end() )
            vmap[vfrom++] = vto++;
        while( cfrom != tds.full_cells_end() )
            cmap[cfrom++] = cto++;
        cto = full_cells_begin();
        while( cto != full_cells_end() )
        {
            for( int i = 0; i <= (std::max)(0, current_dimension()); ++i )
            {
                associate_vertex_with_full_cell(cto, i, vmap[cto->vertex(i)]);
                cto->set_neighbor(i, cmap[cto->neighbor(i)]);
            }
            ++cto;
        }
    }

    // QUERIES

protected:

    bool check_range(int i) const
    {
        if( current_dimension() < 0 )
        {
            return (0 == i);
        }
        return ( (0 <= i) && (i <= current_dimension()) );
    }

public:

    /* returns the current dimension of the full cells in the triangulation. */
    int maximal_dimension() const { return dmax_; } /* Concept */
    int current_dimension() const { return dcur_; } /* Concept */

    size_type number_of_vertices() const /* Concept */
    {
        return this->vertices_.size();
    }
    size_type number_of_full_cells() const /* Concept */
    {
        return this->full_cells_.size();
    }

    bool empty() const /* Concept */
    {
        return current_dimension() == -2;
    }

    Vertex_container & vertices() { return vertices_; }
    const Vertex_container & vertices() const { return vertices_; }
    Full_cell_container & full_cells() { return full_cells_; }
    const Full_cell_container & full_cells() const { return full_cells_; }

    Vertex_handle vertex(Full_cell_handle s, int i) const /* Concept */
    {
        CGAL_precondition(s != Full_cell_handle() && check_range(i));
        return s->vertex(i);
    }

    Vertex_const_handle vertex(Full_cell_const_handle s, int i) const /* Concept */
    {
        CGAL_precondition(s != Full_cell_handle() && check_range(i));
        return s->vertex(i);
    }

    bool is_vertex(Vertex_const_handle v) const /* Concept */
    {
        if( Vertex_const_handle() == v )
            return false;
        Vertex_const_iterator vit = vertices_begin();
        while( vit != vertices_end() && ( v != vit ) )
            ++vit;
        return v == vit;
    }

    bool is_full_cell(Full_cell_const_handle s) const /* Concept */
    {
        if( Full_cell_const_handle() == s )
            return false;
        Full_cell_const_iterator sit = full_cells_begin();
        while( sit != full_cells_end() && ( s != sit ) )
            ++sit;
        return s == sit;
    }

    Full_cell_handle full_cell(Vertex_handle v) const /* Concept */
    {
        CGAL_precondition(v != Vertex_handle());
        return v->full_cell();
    }

    Full_cell_const_handle full_cell(Vertex_const_handle v) const /* Concept */
    {
        CGAL_precondition(Vertex_const_handle() != v);
        return v->full_cell();
    }

    Full_cell_handle neighbor(Full_cell_handle s, int i) const /* Concept */
    {
        CGAL_precondition(Full_cell_handle() != s && check_range(i));
        return s->neighbor(i);
    }

    Full_cell_const_handle neighbor(Full_cell_const_handle s, int i) const/* Concept */
    {
        CGAL_precondition(Full_cell_const_handle() != s && check_range(i));
        return s->neighbor(i);
    }

    int mirror_index(Full_cell_handle s, int i) const /* Concept */
    {
        CGAL_precondition(Full_cell_handle() != s && check_range(i));
        return s->mirror_index(i);
    }

    int mirror_index(Full_cell_const_handle s, int i) const
    {
        CGAL_precondition(Full_cell_const_handle() != s && check_range(i)); /* Concept */
        return s->mirror_index(i);
    }

    int mirror_vertex(Full_cell_handle s, int i) const /* Concept */
    {
        CGAL_precondition(Full_cell_handle() != s && check_range(i));
        return s->mirror_vertex(i);
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - FACETS OPERATIONS

    // works for Face_ = Facet and Face_ = Rotor.
    // NOT DOCUMENTED for the Rotor case...
    template< typename Face_ >
    Full_cell_handle full_cell(const Face_ & f) const /* Concept */
    {
        return f.full_cell();
    }

    // works for Face_ = Facet and Face_ = Rotor.
    // NOT DOCUMENTED for the Rotor case...
    template< class Face_ >
    int index_of_covertex(const Face_ & f) const /* Concept */
    {
        return f.index_of_covertex();
    }

    // NOT DOCUMENTED
    // A Rotor has two covertices
    int index_of_second_covertex(const Rotor & f) const
    {
        return f.index_of_second_covertex();
    }

    // works for Face_ = Facet and Face_ = Rotor.
    // NOT DOCUMENTED...
    template< class Face_ >
    bool is_boundary_facet(const Face_ & f) const
    {
        if( get_visited(neighbor(full_cell(f), index_of_covertex(f))) )
            return false;
        if( ! get_visited(full_cell(f)) )
            return false;
        return true;
    }

    // NOT DOCUMENTED...
    Rotor rotate_rotor(Rotor & f)
    {
        int opposite = mirror_index(full_cell(f), index_of_covertex(f));
        Full_cell_handle s = neighbor(full_cell(f), index_of_covertex(f));
        int new_second = s->index(vertex(full_cell(f), index_of_second_covertex(f)));
        return Rotor(s, new_second, opposite);
    }

    //       NICE UPDATE OPERATIONS

protected:
    void do_insert_increase_dimension(Vertex_handle, Vertex_handle);
public:
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - REMOVALS

    Vertex_handle collapse_face(const Face &); /* Concept */
    void remove_decrease_dimension(Vertex_handle, Vertex_handle); /* Concept */

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - INSERTIONS

    Vertex_handle insert_in_full_cell(Full_cell_handle); /* Concept */
    Vertex_handle insert_in_face(const Face &); /* Concept */
    Vertex_handle insert_in_facet(const Facet &); /* Concept */
    template< typename Forward_iterator >
    Vertex_handle insert_in_hole(Forward_iterator, Forward_iterator, Facet); /* Concept */
    template< typename Forward_iterator, typename OutputIterator >
    Vertex_handle insert_in_hole(Forward_iterator, Forward_iterator, Facet, OutputIterator); /* Concept */

    template< typename OutputIterator >
    Full_cell_handle insert_in_tagged_hole(Vertex_handle, Facet, OutputIterator);

    Vertex_handle insert_increase_dimension(Vertex_handle=Vertex_handle()); /* Concept */

private:

  // Used by insert_in_tagged_hole
  struct IITH_task
  {
    IITH_task(
      Facet boundary_facet_,
      int index_of_inside_cell_in_outside_cell_,
      Full_cell_handle future_neighbor_ = Full_cell_handle(),
      int new_cell_index_in_future_neighbor_ = -1,
      int index_of_future_neighbor_in_new_cell_ = -1)
    : boundary_facet(boundary_facet_),
      index_of_inside_cell_in_outside_cell(index_of_inside_cell_in_outside_cell_),
      future_neighbor(future_neighbor_),
      new_cell_index_in_future_neighbor(new_cell_index_in_future_neighbor_),
      index_of_future_neighbor_in_new_cell(index_of_future_neighbor_in_new_cell_)
    {}

    // "new_cell" is the cell about to be created
    Facet boundary_facet;
    int index_of_inside_cell_in_outside_cell;
    Full_cell_handle future_neighbor;
    int new_cell_index_in_future_neighbor;
    int index_of_future_neighbor_in_new_cell;
  };

  // NOT DOCUMENTED
  void clear_visited_marks(Full_cell_handle) const;

  //  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  DANGEROUS UPDATE OPERATIONS

private:

    // NOT DOCUMENTED
    template< typename FCH > // FCH = Full_cell_[const_]handle
    bool get_visited(FCH c) const
    {
        return c->tds_data().is_visited();
    }

    // NOT DOCUMENTED
    template< typename FCH > // FCH = Full_cell_[const_]handle
    void set_visited(FCH c, bool m) const
    {
        if( m )
            c->tds_data().mark_visited();
        else
            c->tds_data().clear_visited();
    }

public:

    void clear() /* Concept */
    {
        clean_dynamic_memory();
        dcur_ = -2;
    }

    void set_current_dimension(int d) /* Concept */
    {
        CGAL_precondition(-2<=d && d<=maximal_dimension());
        dcur_ = d;
    }

    Full_cell_handle new_full_cell(Full_cell_handle s)
    {
        return full_cells_.emplace(*s);
    }

    Full_cell_handle new_full_cell() /* Concept */
    {
        return full_cells_.emplace(dmax_);
    }

    void delete_full_cell(Full_cell_handle s) /* Concept */
    {
        CGAL_precondition(Full_cell_handle() != s);
        // CGAL_expensive_precondition(is_full_cell(s));
        full_cells_.erase(s);
    }

    template< typename Forward_iterator >
    void delete_full_cells(Forward_iterator start, Forward_iterator end) /* Concept */
    {
        Forward_iterator s = start;
        while( s != end )
            full_cells_.erase(*s++);
    }

    template< class T >
    Vertex_handle new_vertex( const T & t )
    {
        return vertices_.emplace(t);
    }

    Vertex_handle new_vertex() /* Concept */
    {
        return vertices_.emplace();
    }

    void delete_vertex(Vertex_handle v) /* Concept */
    {
        CGAL_precondition( Vertex_handle() != v );
        vertices_.erase(v);
    }

    void associate_vertex_with_full_cell(Full_cell_handle s, int i, Vertex_handle v) /* Concept */
    {
        CGAL_precondition(check_range(i));
        CGAL_precondition(s != Full_cell_handle());
        CGAL_precondition(v != Vertex_handle());
        s->set_vertex(i, v);
        v->set_full_cell(s);
    }

    void set_neighbors(Full_cell_handle s, int i, Full_cell_handle s1, int j) /* Concept */
    {
        CGAL_precondition(check_range(i));
        CGAL_precondition(check_range(j));
        CGAL_precondition(s  != Full_cell_handle());
        CGAL_precondition(s1 != Full_cell_handle());
        s->set_neighbor(i, s1);
        s1->set_neighbor(j, s);
        s->set_mirror_index(i, j);
        s1->set_mirror_index(j, i);
    }

    // SANITY CHECKS

    bool is_valid(bool = true, int = 0) const; /* Concept */

    // NOT DOCUMENTED
    template< class OutStream> void write_graph(OutStream &);

    Vertex_iterator vertices_begin() { return vertices_.begin(); } /* Concept */
    Vertex_iterator vertices_end()   { return vertices_.end();   } /* Concept */
    Full_cell_iterator full_cells_begin() { return full_cells_.begin(); } /* Concept */
    Full_cell_iterator full_cells_end()   { return full_cells_.end();   } /* Concept */

    Vertex_const_iterator vertices_begin() const { return vertices_.begin(); } /* Concept */
    Vertex_const_iterator vertices_end()   const { return vertices_.end();   } /* Concept */
    Full_cell_const_iterator full_cells_begin() const { return full_cells_.begin(); } /* Concept */
    Full_cell_const_iterator full_cells_end()   const { return full_cells_.end();   } /* Concept */

    Facet_iterator facets_begin() /* Concept */
    {
        if( current_dimension() <= 0 )
            return facets_end();
        return Facet_iterator(*this);
    }
    Facet_iterator facets_end() /* Concept */
    {
        return Facet_iterator(*this, 0);
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - FULL CELL GATHERING

    // a traversal predicate for gathering full_cells incident to a given face
    // ``incident'' means that the given face is a subface of the full_cell
    class Incident_full_cell_traversal_predicate
    {
        const Face & f_;
        int dim_;
        const Triangulation_data_structure & tds_;
    public:
        Incident_full_cell_traversal_predicate(const Triangulation_data_structure & tds,
                                               const Face & f)
        : f_(f), tds_(tds)
        {
            dim_ = f.face_dimension();
        }
        bool operator()(const Facet & facet) const
        {
            Vertex_handle v = tds_.full_cell(facet)->vertex(tds_.index_of_covertex(facet));
            for( int i = 0; i <= dim_; ++i )
            {
                if( v == f_.vertex(i) )
                    return false;
            }
            return true;
        }
    };

    // a traversal predicate for gathering full_cells having a given face as subface
    class Star_traversal_predicate
    {
        const Face & f_;
        int dim_;
        const Triangulation_data_structure & tds_;
    public:
        Star_traversal_predicate(const Triangulation_data_structure & tds,
                                 const Face & f)
        : f_(f), tds_(tds)
        {
            dim_ = f.face_dimension();
        }
        bool operator()(const Facet & facet) const
        {
            Full_cell_handle s = tds_.full_cell(facet)->neighbor(tds_.index_of_covertex(facet));
            for( int j = 0; j <= tds_.current_dimension(); ++j )
            {
                for( int i = 0; i <= dim_; ++i )
                    if( s->vertex(j) == f_.vertex(i) )
                        return true;
            }
            return false;
        }
    };

    template< typename TraversalPredicate, typename OutputIterator >
    Facet gather_full_cells(Full_cell_handle, TraversalPredicate &, OutputIterator &) const; /* Concept */
    template< typename OutputIterator >
    OutputIterator incident_full_cells(const Face &, OutputIterator) const; /* Concept */
    template< typename OutputIterator >
    OutputIterator incident_full_cells(Vertex_const_handle, OutputIterator) const; /* Concept */
    template< typename OutputIterator >
    OutputIterator star(const Face &, OutputIterator) const; /* Concept */
#ifndef CGAL_CFG_NO_CPP0X_DEFAULT_TEMPLATE_ARGUMENTS_FOR_FUNCTION_TEMPLATES
    template< typename OutputIterator, typename Comparator = std::less<Vertex_const_handle> >
    OutputIterator incident_upper_faces(Vertex_const_handle v, int dim, OutputIterator out, Comparator cmp = Comparator())
    {
        return incident_faces(v, dim, out, cmp, true);
    }
    template< typename OutputIterator, typename Comparator = std::less<Vertex_const_handle> >
    OutputIterator incident_faces(Vertex_const_handle, int, OutputIterator, Comparator = Comparator(), bool = false) const;
#else
    template< typename OutputIterator, typename Comparator >
    OutputIterator incident_upper_faces(Vertex_const_handle v, int dim, OutputIterator out, Comparator cmp = Comparator())
    {
        return incident_faces(v, dim, out, cmp, true);
    }
    template< typename OutputIterator >
    OutputIterator incident_upper_faces(Vertex_const_handle v, int dim, OutputIterator out)
    {
        return incident_faces(v, dim, out, std::less<Vertex_const_handle>(), true);
    }
    template< typename OutputIterator, typename Comparator >
    OutputIterator incident_faces(Vertex_const_handle, int, OutputIterator, Comparator = Comparator(), bool = false) const;
    template< typename OutputIterator >
    OutputIterator incident_faces(Vertex_const_handle, int, OutputIterator,
        std::less<Vertex_const_handle> = std::less<Vertex_const_handle>(), bool = false) const;
#endif

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - INPUT / OUTPUT

    std::istream & read_full_cells(std::istream &, const std::vector<Vertex_handle> &);
    std::ostream & write_full_cells(std::ostream &, std::map<Vertex_const_handle, int> &) const;

}; // end of ``declaration/definition'' of Triangulation_data_structure<...>

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

// FUNCTIONS THAT ARE MEMBER FUNCTIONS:

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - THE GATHERING METHODS

template< class Dim, class Vb, class Fcb >
template< typename OutputIterator >
OutputIterator
Triangulation_data_structure<Dim, Vb, Fcb>
::incident_full_cells(const Face & f, OutputIterator out) const /* Concept */
{
    // CGAL_expensive_precondition_msg(is_full_cell(f.full_cell()), "the facet does not belong to the Triangulation");
    Incident_full_cell_traversal_predicate tp(*this, f);
    gather_full_cells(f.full_cell(), tp, out);
    return out;
}

template< class Dim, class Vb, class Fcb >
template< typename OutputIterator >
OutputIterator
Triangulation_data_structure<Dim, Vb, Fcb>
::incident_full_cells(Vertex_const_handle v, OutputIterator out) const /* Concept */
{
//    CGAL_expensive_precondition(is_vertex(v));
    CGAL_precondition(Vertex_handle() != v);
    Face f(v->full_cell());
    f.set_index(0, v->full_cell()->index(v));
    return incident_full_cells(f, out);
}

template< class Dim, class Vb, class Fcb >
template< typename OutputIterator >
OutputIterator
Triangulation_data_structure<Dim, Vb, Fcb>
::star(const Face & f, OutputIterator out) const /* Concept */
{
    // CGAL_precondition_msg(is_full_cell(f.full_cell()), "the facet does not belong to the Triangulation");
    Star_traversal_predicate tp(*this, f);
    gather_full_cells(f.full_cell(), tp, out);
    return out;
}

template< class Dim, class Vb, class Fcb >
template< typename TraversalPredicate, typename OutputIterator >
typename Triangulation_data_structure<Dim, Vb, Fcb>::Facet
Triangulation_data_structure<Dim, Vb, Fcb>
::gather_full_cells(Full_cell_handle start,
                    TraversalPredicate & tp,
                    OutputIterator & out) const /* Concept */
{
    std::queue<Full_cell_handle> queue;
    set_visited(start, true);
    queue.push(start);
    const int cur_dim = current_dimension();
    Facet ft;
    while( ! queue.empty() )
    {
        Full_cell_handle s = queue.front();
        queue.pop();
        *out = s;
        ++out;
        for( int i = 0; i <= cur_dim; ++i )
        {
            Full_cell_handle n = s->neighbor(i);
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
template< class Dim, class Vb, class Fcb >
template< typename OutputIterator >
OutputIterator
Triangulation_data_structure<Dim, Vb, Fcb>
::incident_faces(Vertex_const_handle v, int dim, OutputIterator out,
    std::less<Vertex_const_handle> cmp, bool upper_faces) const
{
    return incident_faces<OutputIterator, std::less<Vertex_const_handle> >(v, dim, out, cmp, upper_faces);
}
#endif

template< class Dim, class Vb, class Fcb >
template< typename OutputIterator, typename Comparator >
OutputIterator
Triangulation_data_structure<Dim, Vb, Fcb>
::incident_faces(Vertex_const_handle v, int dim, OutputIterator out, Comparator cmp, bool upper_faces) const
{
    CGAL_precondition( 0 < dim );
    if( dim >= current_dimension() )
        return out;
    typedef std::vector<Full_cell_handle> Simplices;
    Simplices simps;
    simps.reserve(64);
    // gather incident full_cells
    std::back_insert_iterator<Simplices> sout(simps);
    incident_full_cells(v, sout);
    // for storing the handles to the vertices of a full_cell
    typedef std::vector<Vertex_const_handle> Vertices;
    typedef std::vector<int> Indices;
    Vertices vertices(1 + current_dimension());
    Indices sorted_idx(1 + current_dimension());
    // setup Face comparator and Face_set
    typedef internal::Triangulation::Compare_faces_with_common_first_vertex<Self>
        Upper_face_comparator;
    Upper_face_comparator ufc(dim);
    typedef std::set<Face, Upper_face_comparator> Face_set;
    Face_set face_set(ufc);
    for( typename Simplices::const_iterator s = simps.begin(); s != simps.end(); ++s )
    {
        int v_idx(0); // the index of |v| in the sorted full_cell
        // get the vertices of the full_cell and sort them
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
        if( v_idx + dim > current_dimension() )
            continue; // |v| is too far to the right
        // stores the index of the vertices of s in the same order
        // as in |vertices|:
        for( int i = 0; i <= current_dimension(); ++i )
            sorted_idx[i] = (*s)->index(vertices[i]);
        // init state for enumerating all candidate faces:
        internal::Combination_enumerator f_idx(dim, v_idx + 1, current_dimension());
        Face f(*s);
        f.set_index(0, sorted_idx[v_idx]);
        while( ! f_idx.end() )
        {
            for( int i = 0; i < dim; ++i )
                f.set_index(1 + i, sorted_idx[f_idx[i]]);
            face_set.insert(f); // checks if face has already been found

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

template <class Dim, class Vb, class Fcb>
typename Triangulation_data_structure<Dim, Vb, Fcb>::Vertex_handle
Triangulation_data_structure<Dim, Vb, Fcb>
::collapse_face(const Face & f) /* Concept */
{
    const int fd = f.face_dimension();
    CGAL_precondition( (1 <= fd ) && (fd < current_dimension()));
    std::vector<Full_cell_handle> simps;
    // save the Face's vertices:
    Full_cell s;
    for( int i = 0; i <= fd; ++i )
        s.set_vertex(i, f.vertex(i));
    // compute the star of f
    simps.reserve(64);
    std::back_insert_iterator<std::vector<Full_cell_handle> > out(simps);
    star(f, out);
    Vertex_handle v = insert_in_hole(simps.begin(), simps.end(), Facet(f.full_cell(), f.index(0)));
    for( int i = 0; i <= fd; ++i )
        delete_vertex(s.vertex(i));
    return v;
}

template <class Dim, class Vb, class Fcb>
void
Triangulation_data_structure<Dim, Vb, Fcb>
::remove_decrease_dimension(Vertex_handle v, Vertex_handle star) /* Concept */
{
    CGAL_assertion( current_dimension() >= -1 );
    if( -1 == current_dimension() )
    {
        clear();
        return;
    }
    else if( 0 == current_dimension() )
    {
        delete_full_cell(v->full_cell());
        delete_vertex(v);
        star->full_cell()->set_neighbor(0, Full_cell_handle());
        set_current_dimension(-1);
        return;
    }
    else if( 1 == current_dimension() )
    {
        Full_cell_handle s = v->full_cell();
        int star_index;
        if( s->has_vertex(star, star_index) )
            s = s->neighbor(star_index);
        // Here, |star| is not a vertex of |s|, so it's the only finite
        // full_cell
        Full_cell_handle inf1 = s->neighbor(0);
        Full_cell_handle inf2 = s->neighbor(1);
        Vertex_handle v2 = s->vertex(1 - s->index(v));
        delete_vertex(v);
        delete_full_cell(s);
        inf1->set_vertex(1, Vertex_handle());
        inf1->set_vertex(1, Vertex_handle());
        inf2->set_neighbor(1, Full_cell_handle());
        inf2->set_neighbor(1, Full_cell_handle());
        associate_vertex_with_full_cell(inf1, 0, star);
        associate_vertex_with_full_cell(inf2, 0, v2);
        set_neighbors(inf1, 0, inf2, 0);
        set_current_dimension(0);
        return;
    }
    typedef std::vector<Full_cell_handle> Simplices;
    Simplices simps;
    incident_full_cells(v, std::back_inserter(simps));
    for( typename Simplices::iterator it = simps.begin(); it != simps.end(); ++it )
    {
        int v_idx = (*it)->index(v);
        if( ! (*it)->has_vertex(star) )
        {
            delete_full_cell((*it)->neighbor(v_idx));
            for( int i = 0; i <= current_dimension(); ++i )
                (*it)->vertex(i)->set_full_cell(*it);
        }
        else
            star->set_full_cell(*it);
        if( v_idx != current_dimension() )
        {
            (*it)->swap_vertices(v_idx, current_dimension());
            (*it)->swap_vertices(current_dimension() - 2, current_dimension() - 1);
        }
        (*it)->set_vertex(current_dimension(), Vertex_handle());
        (*it)->set_neighbor(current_dimension(), Full_cell_handle());
    }
    set_current_dimension(current_dimension()-1);
    delete_vertex(v);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - THE INSERTION METHODS

template <class Dim, class Vb, class Fcb>
typename Triangulation_data_structure<Dim, Vb, Fcb>::Vertex_handle
Triangulation_data_structure<Dim, Vb, Fcb>
::insert_in_full_cell(Full_cell_handle s) /* Concept */
{
    CGAL_precondition(0 < current_dimension());
    CGAL_precondition(Full_cell_handle() != s);
    // CGAL_expensive_precondition(is_full_cell(s));

    const int cur_dim = current_dimension();
    Vertex_handle v = new_vertex();
    // the full_cell 'fc' is just used to store the handle to all the new full_cells.
    Full_cell fc(maximal_dimension());
    for( int i = 1; i <= cur_dim; ++i )
    {
        Full_cell_handle new_s = new_full_cell(s);
        fc.set_neighbor(i, new_s);
        associate_vertex_with_full_cell(new_s, i, v);
        s->vertex(i-1)->set_full_cell(new_s);
        set_neighbors(new_s, i, neighbor(s, i), mirror_index(s, i));
    }
    fc.set_neighbor(0, s);
    associate_vertex_with_full_cell(s, 0, v);
    for( int i = 0; i <= cur_dim; ++i )
        for( int j = 0; j <= cur_dim; ++j )
        {
            if( j == i ) continue;
            set_neighbors(fc.neighbor(i), j, fc.neighbor(j), i);
        }
    return v;
}

template <class Dim, class Vb, class Fcb >
typename Triangulation_data_structure<Dim, Vb, Fcb>::Vertex_handle
Triangulation_data_structure<Dim, Vb, Fcb>
::insert_in_face(const Face & f) /* Concept */
{
    std::vector<Full_cell_handle> simps;
    simps.reserve(64);
    std::back_insert_iterator<std::vector<Full_cell_handle> > out(simps);
    incident_full_cells(f, out);
    return insert_in_hole(simps.begin(), simps.end(), Facet(f.full_cell(), f.index(0)));
}
template <class Dim, class Vb, class Fcb >
typename Triangulation_data_structure<Dim, Vb, Fcb>::Vertex_handle
Triangulation_data_structure<Dim, Vb, Fcb>
::insert_in_facet(const Facet & ft) /* Concept */
{
    Full_cell_handle s[2];
    s[0] = full_cell(ft);
    int i = index_of_covertex(ft);
    s[1] = s[0]->neighbor(i);
    i = ( i + 1 ) % current_dimension();
    return insert_in_hole(s, s+2, Facet(s[0], i));
}

template <class Dim, class Vb, class Fcb >
template < typename OutputIterator >
typename Triangulation_data_structure<Dim, Vb, Fcb>::Full_cell_handle
Triangulation_data_structure<Dim, Vb, Fcb>
::insert_in_tagged_hole(Vertex_handle v, Facet f,
                        OutputIterator new_full_cells)
{
    CGAL_assertion_msg(is_boundary_facet(f), "starting facet should be on the hole boundary");

    const int cur_dim = current_dimension();
    Full_cell_handle new_s;

    std::queue<IITH_task> task_queue;
    task_queue.push(
    IITH_task(f, mirror_index(full_cell(f), index_of_covertex(f))) );

  while (!task_queue.empty())
  {
    IITH_task task = task_queue.front();
    task_queue.pop();
    
    Full_cell_handle old_s = full_cell(task.boundary_facet);
    const int facet_index = index_of_covertex(task.boundary_facet);
    
    Full_cell_handle outside_neighbor = neighbor(old_s, facet_index);
    // Here, "new_s" might actually be a new cell, but it might also be "old_s"
    // if it has not been treated already in the meantime
    new_s = neighbor(outside_neighbor, task.index_of_inside_cell_in_outside_cell);
    // If the cell has not been treated yet
    if (old_s == new_s)
    {
      new_s = new_full_cell();

      int i(0);
      for ( ; i < facet_index ; ++i)
        associate_vertex_with_full_cell(new_s, i, old_s->vertex(i));
      ++i; // skip facet_index
      for ( ; i <= cur_dim ; ++i)
        associate_vertex_with_full_cell(new_s, i, old_s->vertex(i));
      associate_vertex_with_full_cell(new_s, facet_index, v);
      set_neighbors(new_s,
                    facet_index,
                    outside_neighbor,
                    mirror_index(old_s, facet_index));

      // add the new full_cell to the list of new full_cells
      *new_full_cells++ = new_s;
  
      // check all of |Facet f|'s neighbors
      for (i = 0 ; i <= cur_dim ; ++i)
      {
        if (facet_index == i)
          continue;
        // we define a |Rotor| because it makes it easy to rotate around
        // in a self contained fashion. The corresponding potential
        // boundary facet is Facet(full_cell(rot), index_of_covertex(rot))
        Rotor rot(old_s, i, facet_index);
        // |rot| on line above, stands for Candidate Facet
        while (!is_boundary_facet(rot))
          rot = rotate_rotor(rot);

        // we did find the |i|-th neighbor of Facet(old_s, facet_index)...
        // has it already been extruded to center point |v| ?
        Full_cell_handle inside  = full_cell(rot);
        Full_cell_handle outside = neighbor(inside, index_of_covertex(rot));
        // "m" is the vertex of outside which is not on the boundary
        Vertex_handle m = inside->mirror_vertex(index_of_covertex(rot), current_dimension()); // CJTODO: use mirror_index?
        // "index" is the index of m in "outside"
        int index = outside->index(m);
        // new_neighbor is the inside cell which is registered as the neighbor
        // of the outside cell => it's either a newly created inside cell or an
        // old inside cell which we are about to delete
        Full_cell_handle new_neighbor = outside->neighbor(index);

        // Is new_neighbor still the old neighbor?
        if (new_neighbor == inside)
        {
          task_queue.push(IITH_task(
            Facet(inside, index_of_covertex(rot)), // boundary facet
            index,                        // index_of_inside_cell_in_outside_cell
            new_s,                        // future_neighbor
            i,                            // new_cell_index_in_future_neighbor
            index_of_second_covertex(rot) // index_of_future_neighbor_in_new_cell 
          ));
        }
      }
    }

    // If there is some neighbor stories to fix
    if (task.future_neighbor != Full_cell_handle())
    {
      // now the new neighboring full_cell exists, we link both
      set_neighbors(new_s, 
                    task.index_of_future_neighbor_in_new_cell, 
                    task.future_neighbor, 
                    task.new_cell_index_in_future_neighbor);
    }
  }

  return new_s;
}

template< class Dim, class Vb, class Fcb >
template< typename Forward_iterator, typename OutputIterator >
typename Triangulation_data_structure<Dim, Vb, Fcb>::Vertex_handle
Triangulation_data_structure<Dim, Vb, Fcb>
::insert_in_hole(Forward_iterator start, Forward_iterator end, Facet f,
                 OutputIterator out) /* Concept */
{
    CGAL_expensive_precondition(
            ( std::distance(start, end) == 1 )
         || ( current_dimension() > 1 ) );
    Forward_iterator sit = start;
    while( end != sit )
        set_visited(*sit++, true);
    Vertex_handle v = new_vertex();
    insert_in_tagged_hole(v, f, out);
    delete_full_cells(start, end);
    return v;
}

template< class Dim, class Vb, class Fcb >
template< typename Forward_iterator >
typename Triangulation_data_structure<Dim, Vb, Fcb>::Vertex_handle
Triangulation_data_structure<Dim, Vb, Fcb>
::insert_in_hole(Forward_iterator start, Forward_iterator end, Facet f) /* Concept */
{
    Emptyset_iterator out;
    return insert_in_hole(start, end, f, out);
}

template <class Dim, class Vb, class Fcb>
void
Triangulation_data_structure<Dim, Vb, Fcb>
::clear_visited_marks(Full_cell_handle start) const // NOT DOCUMENTED
{
    CGAL_precondition(start != Full_cell_handle());

    std::queue<Full_cell_handle> queue;
    set_visited(start, false);
    queue.push(start);
    const int cur_dim = current_dimension();
    while( ! queue.empty() )
    {
        Full_cell_handle s = queue.front();
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

template <class Dim, class Vb, class Fcb>
void Triangulation_data_structure<Dim, Vb, Fcb>
::do_insert_increase_dimension(Vertex_handle x, Vertex_handle star)
{
    Full_cell_handle start = full_cells_begin();
    Full_cell_handle swap_me;
    const int cur_dim = current_dimension();
    for( Full_cell_iterator S = full_cells_begin(); S != full_cells_end(); ++S )
    {
        if( Vertex_handle() != S->vertex(cur_dim) )
            continue;
        set_visited(S, true);
        // extends full_cell |S| to include the new vertex as the
        // current_dimension()-th vertex
        associate_vertex_with_full_cell(S, cur_dim, x);
        if( ! S->has_vertex(star) )
        {   // S is bounded, we create its unbounded "twin" full_cell
            Full_cell_handle S_new = new_full_cell();
            set_neighbors(S, cur_dim, S_new, 0);
            associate_vertex_with_full_cell(S_new, 0, star);
            // here, we could be clever so as to get consistent orientation
            for( int k = 1; k <= cur_dim; ++k )
                associate_vertex_with_full_cell(S_new, k, vertex(S, k - 1));
        }
    }
    // now we setup the neighbors
    set_visited(start, false);
    std::queue<Full_cell_handle> queue;
    queue.push(start);
    while( ! queue.empty() )
    {
        Full_cell_handle S = queue.front();
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
            Full_cell_handle S_new = neighbor(S, cur_dim);
            for( int k = 0 ; k < cur_dim ; ++k )
            {
                Full_cell_handle S_opp = neighbor(S, k);
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
        for( Full_cell_iterator S = full_cells_begin(); S != full_cells_end(); ++S )
        {
            if( x != S->vertex(cur_dim) )
                S->swap_vertices(cur_dim - 1, cur_dim);
        }
    }
    if( Full_cell_handle() != swap_me )
        swap_me->swap_vertices(1, 2);
}

template <class Dim, class Vb, class Fcb>
typename Triangulation_data_structure<Dim, Vb, Fcb>::Vertex_handle
Triangulation_data_structure<Dim, Vb, Fcb>
::insert_increase_dimension(Vertex_handle star) /* Concept */
{
    const int prev_cur_dim = current_dimension();
    CGAL_precondition(prev_cur_dim < maximal_dimension());
    if( -2 != current_dimension() )
    {
        CGAL_precondition( Vertex_handle() != star );
        CGAL_expensive_precondition(is_vertex(star));
    }

    set_current_dimension(prev_cur_dim + 1);
    Vertex_handle v = new_vertex();
    switch( prev_cur_dim )
    {
        case -2:
        {   // insertion of the first vertex
            // ( geometrically : infinite vertex )
            Full_cell_handle s = new_full_cell();
            associate_vertex_with_full_cell(s, 0, v);
            break;
        }
        case -1:
        {   // insertion of the second vertex
            // ( geometrically : first finite vertex )
            //we create a triangulation of the 0-sphere, with
            // vertices |star| and |v|
            Full_cell_handle infinite_full_cell = star->full_cell();
            Full_cell_handle finite_full_cell = new_full_cell();
            associate_vertex_with_full_cell(finite_full_cell, 0, v);
            set_neighbors(infinite_full_cell, 0, finite_full_cell, 0);
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

template <class Dimen, class Vb, class Fcb>
bool Triangulation_data_structure<Dimen, Vb, Fcb>
::is_valid(bool verbose, int /* level */) const /* Concept */
{
    Full_cell_const_handle s, t;
    Vertex_const_handle v;
    int i, j, k;

    if( current_dimension() == -2 )
    {
        if( ! vertices_.empty() || ! full_cells_.empty() )
        {
            if( verbose ) CGAL_warning_msg(false, "current dimension is -2 but there are vertices or full_cells");
            return false;
        }
    }

    if( current_dimension() == -1 )
    {
        if ( (number_of_vertices() != 1) || (number_of_full_cells() != 1) )
        {
            if( verbose ) CGAL_warning_msg(false, "current dimension is -1 but there isn't one vertex and one full_cell");
            return false;
        }
    }

    for( v = vertices_begin(); v != vertices_end(); ++v )
    {
        if( ! v->is_valid(verbose) )
            return false;
    }
    
    // FUTURE: for each vertex v, gather incident full_cells. then, check that
    // any full_cell containing v is among those gathered full_cells...

    if( current_dimension() < 0 )
        return true;

    for( s = full_cells_begin(); s != full_cells_end(); ++s )
    {
        if( ! s->is_valid(verbose) )
            return false;
        // check that the full cell has no duplicate vertices
        for( i = 0; i <= current_dimension(); ++i )
            for( j = i + 1; j <= current_dimension(); ++j )
                if( vertex(s,i) == vertex(s,j) )
                {
                    CGAL_warning_msg(false, "a full_cell has two equal vertices");
                    return false;
                }
    }

    for( s = full_cells_begin(); s != full_cells_end(); ++s )
    {
        for( i = 0; i <= current_dimension(); ++i )
            if( (t = neighbor(s,i)) != Full_cell_const_handle() )
            {
                int l = mirror_index(s,i);
                if( s != neighbor(t,l) || i != mirror_index(t,l) )
                {
                    if( verbose ) CGAL_warning_msg(false, "neighbor relation is not symmetric");
                    return false;
                }
                for( j = 0; j <= current_dimension(); ++j )
                    if( j != i )
                    {
                        // j must also occur as a vertex of t
                        for( k = 0; k <= current_dimension() && ( vertex(s,j) != vertex(t,k) || k == l); ++k )
                            ;
                        if( k > current_dimension() )
                        {
                            if( verbose ) CGAL_warning_msg(false, "too few shared vertices between neighbors full_cells.");
                            return false;
                        }
                    }
            }
            else
            {
                if( verbose ) CGAL_warning_msg(false, "full_cell has a NULL neighbor");
                return false;
            }
    }
    return true;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - INPUT / OUTPUT

// NOT DOCUMENTED
template <class Dim, class Vb, class Fcb>
template <class OutStream>
void Triangulation_data_structure<Dim, Vb, Fcb>
::write_graph(OutStream & os)
{
    std::vector<std::set<int> > edges;
    os << number_of_vertices() + 1; // add the vertex at infinity
    int count(1);
    for( Vertex_iterator vit = vertices_begin(); vit != vertices_end(); ++vit )
        vit->idx_ = count++;
    edges.resize(number_of_vertices()+1);
    for( Full_cell_iterator sit = full_cells_begin(); sit != full_cells_end(); ++sit )
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
    for( std::size_t i = 0; i < edges.size(); ++i )
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
template<class Dimen, class Vb, class Fcb>
std::istream &
Triangulation_data_structure<Dimen, Vb, Fcb>
::read_full_cells(std::istream & is, const std::vector<Vertex_handle> & vertices)
{
    std::size_t m; // number of full_cells
    int index;
    const int cd = current_dimension();
    if( is_ascii(is) )
        is >> m;
    else
        read(is, m, io_Read_write());

    std::vector<Full_cell_handle> full_cells;
    full_cells.reserve(m);
    // read the vertices of each full_cell
    std::size_t i = 0;
    while( i < m )
    {
        Full_cell_handle s = new_full_cell();
        full_cells.push_back(s);
        for( int j = 0; j <= cd; ++j )
        {
            if( is_ascii(is) )
                is >> index;
            else
                read(is, index);
            s->set_vertex(j, vertices[index]);
        }
        // read other non-combinatorial information for the full_cells
        is >> (*s);
        ++i;
    }

    // read the neighbors of each full_cell
    i = 0;
    if( is_ascii(is) )
        while( i < m )
    {
        for( int j = 0; j <= cd; ++j )
        {
            is >> index;
            full_cells[i]->set_neighbor(j, full_cells[index]);
        }
        ++i;
    }
    else
        while( i < m )
    {
        for( int j = 0; j <= cd; ++j )
        {
            read(is, index);
            full_cells[i]->set_neighbor(j, full_cells[index]);
        }
        ++i;
    }

    // compute the mirror indices
    for( i = 0; i < m; ++i )
    {
        Full_cell_handle s = full_cells[i];
        for( int j = 0; j <= cd; ++j )
        {
            if( -1 != s->mirror_index(j) )
                continue;
            Full_cell_handle n = s->neighbor(j);
            int k = 0;
            Full_cell_handle nn = n->neighbor(k);
            while( s != nn )
                nn = n->neighbor(++k);
            s->set_mirror_index(j,k);
            n->set_mirror_index(k,j);
        }
    }
    return is;
}

// NOT DOCUMENTED...
template<class Dimen, class Vb, class Fcb>
std::ostream &
Triangulation_data_structure<Dimen, Vb, Fcb>
::write_full_cells(std::ostream & os, std::map<Vertex_const_handle, int> & index_of_vertex) const
{
    std::map<Full_cell_const_handle, int> index_of_full_cell;

    std::size_t m = number_of_full_cells();

    if( is_ascii(os) )
        os << std::endl << m;
    else
        write(os, m, io_Read_write());

    const int cur_dim = current_dimension();
    // write the vertex indices of each full_cell
    int i = 0;
    for( Full_cell_const_iterator it = full_cells_begin(); it != full_cells_end(); ++it )
    {
        index_of_full_cell[it] = i++;
        if( is_ascii(os) )
            os << std::endl;
        for( int j = 0; j <= cur_dim; ++j )
        {
            if( is_ascii(os) )
                os << ' ' << index_of_vertex[it->vertex(j)];
            else
                write(os, index_of_vertex[it->vertex(j)]);
        }
        // write other non-combinatorial information for the full_cells
        os << (*it);
    }

    CGAL_assertion( (std::size_t) i == m );

    // write the neighbors of each full_cell
    if( is_ascii(os) )
        for( Full_cell_const_iterator it = full_cells_begin(); it != full_cells_end(); ++it )
        {
            os << std::endl;
            for( int j = 0; j <= cur_dim; ++j )
                os << ' ' << index_of_full_cell[it->neighbor(j)];
        }
    else
        for( Full_cell_const_iterator it = full_cells_begin(); it != full_cells_end(); ++it )
        {
            for( int j = 0; j <= cur_dim; ++j )
                write(os, index_of_full_cell[it->neighbor(j)]);
        }

    return os;
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

// FUNCTIONS THAT ARE NOT MEMBER FUNCTIONS:

template<class Dimen, class Vb, class Fcb>
std::istream &
operator>>(std::istream & is, Triangulation_data_structure<Dimen, Vb, Fcb> & tr)
  // reads :
  // - the dimensions (maximal and current)
  // - the number of finite vertices
  // - the non combinatorial information on vertices (point, etc)
  // - the number of full_cells
  // - the full_cells by the indices of their vertices in the preceding list
  // of vertices, plus the non combinatorial information on each full_cell
  // - the neighbors of each full_cell by their index in the preceding list
{
    typedef Triangulation_data_structure<Dimen, Vb, Fcb> TDS;
    typedef typename TDS::Vertex_handle         Vertex_handle;

    // read current dimension and number of vertices
    std::size_t n;
    int cd;
    if( is_ascii(is) )
        is >> cd >> n;
    else
    {
        read(is, cd);
        read(is, n, io_Read_write());
    }

    CGAL_assertion_msg( cd <= tr.maximal_dimension(), "input Triangulation_data_structure has too high dimension");

    tr.clear();
    tr.set_current_dimension(cd);

    if( n == 0 )
        return is;

    std::vector<Vertex_handle> vertices;
    vertices.resize(n);

    // read the vertices:
    std::size_t i(0);
    while( i < n )
    {
        vertices[i] = tr.new_vertex();
        is >> (*vertices[i]); // read a vertex
        ++i;
    }

    // now, read the combinatorial information
    return tr.read_full_cells(is, vertices);
}

template<class Dimen, class Vb, class Fcb>
std::ostream &
operator<<(std::ostream & os, const Triangulation_data_structure<Dimen, Vb, Fcb> & tr)
  // writes :
  // - the dimensions (maximal and current)
  // - the number of finite vertices
  // - the non combinatorial information on vertices (point, etc)
  // - the number of full cells
  // - the full cells by the indices of their vertices in the preceding list
  // of vertices, plus the non combinatorial information on each full_cell
  // - the neighbors of each full_cell by their index in the preceding list
{
    typedef Triangulation_data_structure<Dimen, Vb, Fcb> TDS;
    typedef typename TDS::Vertex_const_handle         Vertex_handle;
    typedef typename TDS::Vertex_const_iterator       Vertex_iterator;

    // outputs dimension and number of vertices
    std::size_t n = tr.number_of_vertices();
    if( is_ascii(os) )
        os << tr.current_dimension() << std::endl << n;
    else
    {
        write(os, tr.current_dimension());
        write(os, n, io_Read_write());
    }

    if( n == 0 )
        return os;

    // write the vertices
    std::map<Vertex_handle, int> index_of_vertex;
    int i = 0;
    for( Vertex_iterator it = tr.vertices_begin(); it != tr.vertices_end(); ++it, ++i )
    {
        os << *it; // write the vertex
        if (is_ascii(os))
            os << std::endl;
        index_of_vertex[it] = i;
    }
    CGAL_assertion( (std::size_t) i == n );

    // output the combinatorial information
    return tr.write_full_cells(os, index_of_vertex);
}

} //namespace CGAL

#endif // CGAL_TRIANGULATION_DATA_STRUCTURE_H
