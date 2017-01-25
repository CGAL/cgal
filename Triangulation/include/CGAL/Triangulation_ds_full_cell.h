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

#ifndef CGAL_TRIANGULATION_DS_FULL_CELL_H
#define CGAL_TRIANGULATION_DS_FULL_CELL_H

#include <CGAL/license/Triangulation.h>


#include <CGAL/TDS_full_cell_default_storage_policy.h>
#include <CGAL/TDS_full_cell_mirror_storage_policy.h>
#include <CGAL/internal/Triangulation/Dummy_TDS.h>
#include <CGAL/Dimension.h>
#include <CGAL/Default.h>
#include <CGAL/array.h>

namespace CGAL {

template< class TDS = void, typename FullCellStoragePolicy = Default >
class Triangulation_ds_full_cell
{
    typedef typename Default::Get<FullCellStoragePolicy, TDS_full_cell_default_storage_policy>::type
                                                   Storage_policy;
    typedef Triangulation_ds_full_cell<TDS>        Self;
    typedef typename TDS::Maximal_dimension        Maximal_dimension;

public:
    typedef TDS                                    Triangulation_data_structure;
    typedef typename TDS::Face                     Face;
    typedef typename TDS::Vertex_handle            Vertex_handle; /* Concept */
    typedef typename TDS::Vertex_const_handle      Vertex_const_handle;
    typedef typename TDS::Full_cell_handle         Full_cell_handle; /* Concept */
    typedef typename TDS::Full_cell_const_handle   Full_cell_const_handle;
    typedef typename TDS::Full_cell_data           TDS_data; /* data that the TDS wants to be stored here */
    template< typename TDS2 >
    struct Rebind_TDS /* Concept */
    {
        typedef Triangulation_ds_full_cell<TDS2, FullCellStoragePolicy> Other;
    };

private: // STORAGE
    typedef TFC_data< Vertex_handle, Full_cell_handle,
                      Maximal_dimension, Storage_policy >   Combinatorics;
    friend struct TFC_data< Vertex_handle, Full_cell_handle,
                      Maximal_dimension, Storage_policy >;
    // array of vertices
    typedef typename Combinatorics::Vertex_handle_array     Vertex_handle_array;
    // neighbor simplices
    typedef typename Combinatorics::Full_cell_handle_array    Full_cell_handle_array;

    // NOT DOCUMENTED...
    typename Combinatorics::Xor_type xor_of_vertices(const int cur_dim) const
    {
        return combinatorics_.xor_of_vertices(cur_dim);
    }

public:
    typedef typename Vertex_handle_array::const_iterator    Vertex_handle_const_iterator;
    typedef Vertex_handle_const_iterator    Vertex_handle_iterator; /* Concept */

    Triangulation_ds_full_cell(const int dmax) /* Concept */
    : combinatorics_(dmax), tds_data_()
    {
        CGAL_assertion( dmax > 0 );
        for( int i = 0; i <= dmax; ++i )
        {
            set_neighbor(i, Full_cell_handle());
            set_vertex(i, Vertex_handle());
            set_mirror_index(i, -1);
        }
    }

    Triangulation_ds_full_cell(const Triangulation_ds_full_cell & s) /* Concept */
    : combinatorics_(s.combinatorics_), tds_data_(s.tds_data_)
    {}

    ~Triangulation_ds_full_cell() {}

    int maximal_dimension() const /* Concept */
    {
        return static_cast<int>(vertices().size() - 1);
    }

    Vertex_handle_const_iterator vertices_begin() const /* Concept */
    {
        return vertices().begin();
    }

    Vertex_handle_const_iterator vertices_end() const /* Concept */
    {
        return vertices().end();
    }

    Vertex_handle vertex(const int i) const /* Concept */
    {
        CGAL_precondition(0<=i && i<=maximal_dimension());
        return vertices()[i];
    }

    Full_cell_handle neighbor(const int i) const /* Concept */
    {
        CGAL_precondition(0<=i && i<=maximal_dimension());
        return neighbors()[i];
    }

    int mirror_index(const int i) const /* Concept */
    {
        CGAL_precondition(0<=i && i<=maximal_dimension());
        return combinatorics_.mirror_index(i);
    }

    // Advanced...
    Vertex_handle mirror_vertex(const int i, const int cur_dim) const /* Concept */
    {
        CGAL_precondition(0<=i && i<=maximal_dimension());
        return combinatorics_.mirror_vertex(i, cur_dim);
    }

    int index(Full_cell_const_handle s) const /* Concept */
    {
        // WE ASSUME THE FULL CELL WE ARE LOOKING FOR INDEED EXISTS !
        CGAL_precondition(has_neighbor(s));
        int index(0);
        while( neighbor(index) != s )
            ++index;
        return index;
    }

    int index(Vertex_const_handle v) const /* Concept */
    {
        // WE ASSUME THE VERTEX WE ARE LOOKING FOR INDEED EXISTS !
        CGAL_precondition(has_vertex(v));
        int index(0);
        while( vertex(index) != v )
            ++index;
        return index;
    }

    void set_vertex(const int i, Vertex_handle v) /* Concept */
    {
        CGAL_precondition(0<=i && i<=maximal_dimension());
        vertices()[i] = v;
    }

    void set_neighbor(const int i, Full_cell_handle s) /* Concept */
    {
        CGAL_precondition(0<=i && i<=maximal_dimension());
        neighbors()[i] = s;
    }

    void set_mirror_index(const int i, const int index) /* Concept */
    {
        CGAL_precondition(0<=i && i<=maximal_dimension());
        combinatorics_.set_mirror_index(i, index);
    }

    bool has_vertex(Vertex_const_handle v) const /* Concept */
    {
        int index;
        return has_vertex(v, index);
    }

    bool has_vertex(Vertex_const_handle v, int & index) const /* Concept */
    {
        const int d = maximal_dimension();
        index = 0;
        while( (index <= d) && (vertex(index) != v) )
            ++index;
        return (index <= d);
    }

    bool has_neighbor(Full_cell_const_handle s) const /* Concept */
    {
        int index;
        return has_neighbor(s, index);
    }

    bool has_neighbor(Full_cell_const_handle s, int & index) const /* Concept */
    {
        const int d = maximal_dimension();
        index = 0;
        while( (index <= d) && (neighbor(index) != s) )
            ++index;
        return (index <= d);
    }

    void swap_vertices(const int d1, const int d2) /* Concept */
    {
        CGAL_precondition(0 <= d1 && d1<=maximal_dimension());
        CGAL_precondition(0 <= d2 && d2<=maximal_dimension());
        combinatorics_.swap_vertices(d1, d2);
    }

    const TDS_data & tds_data() const { return tds_data_; } /* Concept */
    TDS_data & tds_data() { return tds_data_; } /* Concept */

    void*   for_compact_container() const { return combinatorics_.for_compact_container(); }
    void* & for_compact_container() { return combinatorics_.for_compact_container(); }

    bool is_valid(bool verbose = false, int = 0) const /* Concept */
    {
        const int d = maximal_dimension();
        int i(0);
        // test that the non-null Vertex_handles come first, before all null ones
        while( i <= d && vertex(i) != Vertex_handle() ) ++i;
        while( i <= d && vertex(i) == Vertex_handle() ) ++i;
        if( i <= d )
        {
            if( verbose ) CGAL_warning_msg(false, "full cell has garbage handles to vertices.");
            return false;
        }
        for( i = 0; i <= d; ++i )
        {
            if( Vertex_handle() == vertex(i) )
                break; // there are no more vertices
            Full_cell_handle n(neighbor(i));
            if( Full_cell_handle() != n )
            {
                int mirror_idx(mirror_index(i));
                if( n->neighbor(mirror_idx) == Full_cell_handle() )
                {
                    if( verbose ) CGAL_warning_msg(false, "neighbor has no back-neighbor.");
                    return false;
                }
                if( &(*(n->neighbor(mirror_idx))) != this )
                {
                    if( verbose ) CGAL_warning_msg(false, "neighbor does not point back to correct full cell.");
                    return false;
                }
            }
        }
        return true;
    }

private:
    // access to data members:
    Full_cell_handle_array & neighbors() {return combinatorics_.neighbors_; }
    const Full_cell_handle_array & neighbors() const {return combinatorics_.neighbors_; }
    Vertex_handle_array & vertices() {return combinatorics_.vertices_; }
    const Vertex_handle_array & vertices() const {return combinatorics_.vertices_; }

    // DATA MEMBERS
    Combinatorics       combinatorics_;
    mutable TDS_data    tds_data_;
};

// FUNCTIONS THAT ARE NOT MEMBER FUNCTIONS:

template < typename TDS, typename SSP >
std::ostream &
operator<<(std::ostream & O, const Triangulation_ds_full_cell<TDS,SSP> &) /* Concept */
{
    /*if( is_ascii(O) )
    {
        // os << '\n';
    }
    else {}*/
    return O;
}

template < typename TDS, typename SSP >
std::istream &
operator>>(std::istream & I, Triangulation_ds_full_cell<TDS,SSP> &) /* Concept */
{
    /*if( is_ascii(I) )
    {}
    else {}*/
    return I;
}

// Special case: specialization when template parameter is void.

// we must declare it for each possible full_cell storage policy because :
// (GCC error:) default template arguments may not be used in partial specializations
template< typename StoragePolicy >
class Triangulation_ds_full_cell<void, StoragePolicy>
{
public:
    typedef internal::Triangulation::Dummy_TDS  TDS;
    typedef TDS                                 Triangulation_data_structure;
    typedef TDS::Vertex_handle                  Vertex_handle;
    typedef TDS::Vertex_const_handle            Vertex_const_handle;
    typedef TDS::Full_cell_handle               Full_cell_handle;
    typedef TDS::Full_cell_const_handle         Full_cell_const_handle;
    typedef TDS::Vertex_handle_const_iterator   Vertex_handle_const_iterator;
    typedef TDS::Full_cell_data                 TDS_data;
    template <typename TDS2>
    struct Rebind_TDS
    {
        typedef Triangulation_ds_full_cell<TDS2, StoragePolicy> Other;
    };
    Vertex_handle_const_iterator vertices_begin();
    Vertex_handle_const_iterator vertices_end();
};

} //namespace CGAL

#endif // CGAL_TRIANGULATION_DS_FULL_CELL_H
