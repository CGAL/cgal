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

#ifndef CGAL_TRIANGULATION_DS_FULL_CELL_H
#define CGAL_TRIANGULATION_DS_FULL_CELL_H

#include <CGAL/TDS_full_cell_default_storage_policy.h>
#include <CGAL/TDS_full_cell_mirror_storage_policy.h>
#include <CGAL/internal/Triangulation/Dummy_TDS.h>
#include <CGAL/Dimension.h>
#include <CGAL/Default.h>
#include <CGAL/array.h>
#include <vector>

namespace CGAL {

template< class TDS = void, typename FullCellStoragePolicy = Default >
class Triangulation_ds_full_cell
{
    typedef typename Default::Get<FullCellStoragePolicy, TDS_full_cell_default_storage_policy>::type
                                                    Storage_policy;
    typedef Triangulation_ds_full_cell<TDS>           Self;
    typedef typename TDS::Ambient_dimension        Ambient_dimension;

public:
    typedef typename TDS::Face                     Face;
    typedef typename TDS::Vertex_handle            Vertex_handle;
    typedef typename TDS::Full_cell_handle           Full_cell_handle;
    typedef typename TDS::Vertex_const_handle      Vertex_const_handle;
    typedef typename TDS::Full_cell_const_handle     Full_cell_const_handle;
    template< typename TDS2 >
    struct Rebind_TDS
    {
        typedef Triangulation_ds_full_cell<TDS2, FullCellStoragePolicy> Other;
    };

private: // STORAGE
    typedef TFC_data< Vertex_handle, Full_cell_handle,
                      Ambient_dimension, Storage_policy >   Combinatorics;
    friend class TFC_data< Vertex_handle, Full_cell_handle,
                      Ambient_dimension, Storage_policy >;
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

    Triangulation_ds_full_cell(const int dmax)
    : flags_(0), combinatorics_(dmax)
    {
		CGAL_assertion( dmax > 0 );
        for( int i = 0; i <= dmax; ++i )
        {
            set_neighbor(i, Full_cell_handle());
            set_vertex(i, Vertex_handle());
            set_mirror_index(i, -1);
        }
    }

    Triangulation_ds_full_cell(const Triangulation_ds_full_cell & s)
    : flags_(s.flags_), combinatorics_(s.combinatorics_)
    {}

    ~Triangulation_ds_full_cell() {}

    int ambient_dimension() const
    {
        return vertices().size() - 1;
    }

    Vertex_handle_const_iterator vertices_begin() const
    {
        return vertices().begin();
    }

    Vertex_handle_const_iterator vertices_end() const
    {
        return vertices().end();
    }

    Vertex_handle vertex(const int i) const
    {
        CGAL_precondition(0<=i && i<=ambient_dimension());
        return vertices()[i];
    }

    Full_cell_handle neighbor(const int i) const
    {
        CGAL_precondition(0<=i && i<=ambient_dimension());
        return neighbors()[i];
    }

    int mirror_index(const int i) const
    {
        CGAL_precondition(0<=i && i<=ambient_dimension());
        return combinatorics_.mirror_index(i);
    }

	// Advanced...
    Vertex_handle mirror_vertex(const int i, const int cur_dim) const
    {
        CGAL_precondition(0<=i && i<=ambient_dimension());
        return combinatorics_.mirror_vertex(i, cur_dim);
    }

    int index_of(Full_cell_const_handle s) const
    {
        // WE ASSUME THE SIMPLEX WE ARE LOOKING FOR INDEED EXISTS !
        CGAL_precondition(has_neighbor(s));
        int index(0);
        while( neighbor(index) != s )
            ++index;
        return index;
    }

    int index_of(Vertex_const_handle v) const
    {
        // WE ASSUME THE VERTEX WE ARE LOOKING FOR INDEED EXISTS !
        CGAL_precondition(has_vertex(v));
        int index(0);
        while( vertex(index) != v )
            ++index;
        return index;
    }

	void set_vertex(const int i, Vertex_handle v)
	{
		CGAL_precondition(0<=i && i<=ambient_dimension());
		vertices()[i] = v;
	}

	void set_neighbor(const int i, Full_cell_handle s)
	{
		CGAL_precondition(0<=i && i<=ambient_dimension());
		neighbors()[i] = s;
	}

	void set_mirror_index(const int i, const int index)
	{
		CGAL_precondition(0<=i && i<=ambient_dimension());
		combinatorics_.set_mirror_index(i, index);
	}

    bool has_vertex(Vertex_const_handle v) const
    {
        int index;
        return has_vertex(v, index);
    }

    bool has_vertex(Vertex_const_handle v, int & index) const
    {
        const int d = ambient_dimension();
        index = 0;
        while( (index <= d) && (vertex(index) != v) )
            ++index;
        return (index <= d);
    }

    bool has_neighbor(Full_cell_const_handle s) const
    {
        int index;
        return has_neighbor(s, index);
    }

    bool has_neighbor(Full_cell_const_handle s, int & index) const
    {
        const int d = ambient_dimension();
        index = 0;
        while( (index <= d) && (neighbor(index) != s) )
            ++index;
        return (index <= d);
    }

    void swap_vertices(const int d1, const int d2)
    {
        CGAL_precondition(0 <= d1 && d1<=ambient_dimension());
        CGAL_precondition(0 <= d2 && d2<=ambient_dimension());
        combinatorics_.swap_vertices(d1, d2);
    }

    unsigned int get_flags() const { return flags_; }
    // Don't forget that member variable flags_ is mutable...
    void set_flags(unsigned int f) const { flags_ = f; }

    void*   for_compact_container() const { return combinatorics_.for_compact_container(); }
    void* & for_compact_container() { return combinatorics_.for_compact_container(); }

    bool is_valid(bool verbose = true, int /* level */ = 0) const
    {
        const int d = ambient_dimension();
        for( int i = 0; i <= d; ++i )
        {
            if( Vertex_handle() != vertex(i) )
            {
                if( Full_cell_handle() == neighbor(i) )
                {
                    if( verbose ) CGAL_warning_msg(false, "vertex has no opposite full cell.");
                    return false;
                }
                // Here, we can't check if neighbor(i) counts *this as a neighbor
                // because we can't construct a Full_cell_handle to *this...
                // So we have to do this check in the `parent' class (TDS)
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
    // |flags_| is the 'visited' mark when traversing. |flag_| is also used
    // in delaunay_triangulation_d for finding simplices in conflict with a
    // newly inserted point
    mutable unsigned int    flags_;
    Combinatorics           combinatorics_;
};

// FUNCTIONS THAT ARE NOT MEMBER FUNCTIONS:

template < typename TDS, typename SSP >
std::ostream &
operator<<(std::ostream & O, const Triangulation_ds_full_cell<TDS,SSP> & s)
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
operator>>(std::istream & I, Triangulation_ds_full_cell<TDS,SSP> & s)
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
    typedef internal::Triangulation::Dummy_TDS          TDS;
    typedef TDS::Vertex_handle   Vertex_handle;
    typedef TDS::Vertex_const_handle   Vertex_const_handle;
    typedef TDS::Full_cell_handle  Full_cell_handle;
    typedef TDS::Full_cell_const_handle  Full_cell_const_handle;
    typedef TDS::Vertex_handle_const_iterator Vertex_handle_const_iterator;
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
