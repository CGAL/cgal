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

#ifndef CGAL_TRIANGULATION_DS_VERTEX_H
#define CGAL_TRIANGULATION_DS_VERTEX_H

#include <CGAL/Compact_container.h>
#include <CGAL/internal/Triangulation/Dummy_TDS.h>

namespace CGAL {

/* The template parameter TDS must be a model of the concept
 * 'TriangulationDataStructure' that stores vertices of type
 * 'Triangulation_ds_vertex<TDS>'
 */
template< class TDS = void >
class Triangulation_ds_vertex 
{
    typedef Triangulation_ds_vertex<TDS> Self;

public:
    typedef typename TDS::Simplex_handle       Simplex_handle;

    template <typename TDS2>
    struct Rebind_TDS
    {
        typedef Triangulation_ds_vertex<TDS2> Other;
    };

protected: // DATA MEMBERS
    Simplex_handle simplex_; // A handle to an adjacent simplex

public:	
    // Constructs a vertex with adjacent simplex 's'
    Triangulation_ds_vertex(Simplex_handle s) : simplex_(s)
    {
        CGAL_assertion( Simplex_handle() != s );
    }
    // Constructs a vertex with no adjacent simplex
    Triangulation_ds_vertex() : simplex_() {}

    ~Triangulation_ds_vertex() {}

    /// Set 's' as an adjacent simplex
    void set_simplex(Simplex_handle s)
    {
		CGAL_precondition( Simplex_handle() != s );
        simplex_ = s;
    }

    /// Returns a simplex adjacent to the vertex
    Simplex_handle simplex() const
    {
        return simplex_;
    }

    bool is_valid(bool verbose = true, int /* level */ = 0) const
    {
        if( Simplex_handle() == simplex() )
        {
            if( verbose ) CGAL_warning_msg(false, "vertex has no incident simplex.");
            return false;
        }
        return true;
    }

public: // FOR MEMORY MANAGEMENT

    void*   for_compact_container() const { return simplex_.for_compact_container(); }
    void* & for_compact_container()       { return simplex_.for_compact_container(); }

};  // end of Triangulation_ds_vertex

// FUNCTIONS THAT ARE NOT MEMBER FUNCTIONS:

template < class TDS >
std::istream &
operator>>(std::istream & is, Triangulation_ds_vertex<TDS> & v)
{
    /*if( is_ascii(is) )
    {}
    else {}*/
    return is;
}

template< class TDS >
std::ostream &
operator<<(std::ostream & os, const Triangulation_ds_vertex<TDS> & v)
{
    /*if( is_ascii(os) )
    {
        os << '\n';
    }
    else {}*/
    return os;
}

// Special case: specialization when template parameter is void.

template<>
class Triangulation_ds_vertex<void>
{
public:
    typedef internal::Triangulation::Dummy_TDS  Triangulation_ds;
    typedef Triangulation_ds::Simplex_handle     Simplex_handle;
    template <typename TDS2>
    struct Rebind_TDS
    {
        typedef Triangulation_ds_vertex<TDS2> Other;
    };
};

} //namespace CGAL

#endif // CGAL_TRIANGULATION_DS_VERTEX_H
