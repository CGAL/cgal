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

//#include <CGAL/basic.h>
//#include <CGAL/Iterator_project.h>
#include <CGAL/Compact_container.h>
#include <CGAL/internal/Pure_complex/Dummy_PCDS.h>

namespace CGAL {

/* TDS must be a model of the concept 'PureComplexDataStructure' that 
    stores vertices of type 'Pure_complex_ds_vertex<Pure_complex>'
 */

template< class TDS = void >
class Pure_complex_ds_vertex 
{
    typedef Pure_complex_ds_vertex<TDS> Self;

public:
    typedef typename TDS::Simplex_handle       Simplex_handle;

    template <typename PC2>
    struct Rebind_TDS
    {
        typedef Pure_complex_ds_vertex<PC2> Other;
    };

protected: // DATA MEMBERS
    Simplex_handle simplex_; // A handle to an adjacent simplex

public:	
    // Constructs a vertex with adjacent simplex 's'
    Pure_complex_ds_vertex(Simplex_handle s) : simplex_(s)
    {
        CGAL_assertion( Simplex_handle() != s );
    }
    // Constructs a vertex with no adjacent simplex
    Pure_complex_ds_vertex() : simplex_() {}

    ~Pure_complex_ds_vertex() {}

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

};  // end of Pure_complex_ds_vertex

// FUNCTIONS THAT ARE NOT MEMBER FUNCTIONS:

template < class TDS >
std::istream &
operator>>(std::istream & is, Pure_complex_ds_vertex<TDS> & v)
{
    /*if( is_ascii(is) )
    {}
    else {}*/
    return is;
}

template< class TDS >
std::ostream &
operator<<(std::ostream & os, const Pure_complex_ds_vertex<TDS> & v)
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
class Pure_complex_ds_vertex<void>
{
public:
    typedef internal::Triangulation::Dummy_TDS  Pure_complex_ds;
    typedef Pure_complex_ds::Simplex_handle     Simplex_handle;
    template <typename PC2>
    struct Rebind_TDS
    {
        typedef Pure_complex_ds_vertex<PC2> Other;
    };
};

} //namespace CGAL

#endif // CGAL_TRIANGULATION_DS_VERTEX_H
