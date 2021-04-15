// Copyright (c) 2009-2014 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)    : Samuel Hornus

#ifndef CGAL_TRIANGULATION_DS_VERTEX_H
#define CGAL_TRIANGULATION_DS_VERTEX_H

#include <CGAL/license/Triangulation.h>

#include <CGAL/disable_warnings.h>

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
    typedef Triangulation_ds_vertex<TDS>    Self;

public:
    typedef TDS                             Triangulation_data_structure;
    typedef typename TDS::Full_cell_handle  Full_cell_handle; /* Concept */

    template <typename TDS2>
    struct Rebind_TDS /* Concept */
    {
        typedef Triangulation_ds_vertex<TDS2> Other;
    };

protected: // DATA MEMBERS
    Full_cell_handle full_cell_; // A handle to an incident full_cell

public:
    // Constructs a vertex with incident full_cell 's'
    Triangulation_ds_vertex(Full_cell_handle s) : full_cell_(s) /* Concept */
    {
        CGAL_assertion( Full_cell_handle() != s );
    }
    // Constructs a vertex with no incident full_cell
    Triangulation_ds_vertex() : full_cell_() {} /* Concept */

    ~Triangulation_ds_vertex() {}

    /// Set 's' as an incident full_cell
    void set_full_cell(Full_cell_handle s) /* Concept */
    {
        full_cell_ = s;
    }

    /// Returns a full_cell incident to the vertex
    Full_cell_handle full_cell() const /* Concept */
    {
        return full_cell_;
    }

    bool is_valid(bool verbose = false, int /* level */ = 0) const /* Concept */
    {
        if( Full_cell_handle() == full_cell() )
        {
            if( verbose )
                CGAL_warning_msg(false, "vertex has no incident full cell.");
            return false;
        }
        bool found(false);
        // These two typename below are OK because TDS fulfils the
        // TriangulationDataStructure concept.
        typename TDS::Full_cell::Vertex_handle_iterator vit(full_cell()->vertices_begin());
        typedef typename TDS::Vertex_handle Vertex_handle;
        while( vit != full_cell()->vertices_end() )
        {
            if( Vertex_handle() == *vit )
                break; // The full cell has no more vertices
            if( this == &(**vit) )
            {
                found = true;
                break;
            }
            ++vit;
        }
        if( ! found )
        {
            if( verbose )
                CGAL_warning_msg(false, "vertex's adjacent full cell does not contain that vertex.");
            return false;
        }
        return true;
    }

public: // FOR MEMORY MANAGEMENT

    void*   for_compact_container() const { return full_cell_.for_compact_container(); }
    void    for_compact_container(void *p){ full_cell_.for_compact_container(p); }

};  // end of Triangulation_ds_vertex

// FUNCTIONS THAT ARE NOT MEMBER FUNCTIONS:

template < class TDS >
std::istream &
operator>>(std::istream & is, Triangulation_ds_vertex<TDS> &) /* Concept */
{
    /*if( is_ascii(is) )
    {}
    else {}*/
    return is;
}

template< class TDS >
std::ostream &
operator<<(std::ostream & os, const Triangulation_ds_vertex<TDS> &) /* Concept */
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
    typedef internal::Triangulation::Dummy_TDS Triangulation_data_structure;
    typedef Triangulation_data_structure::Full_cell_handle Full_cell_handle; /* Concept */
    template <typename TDS2>
    struct Rebind_TDS /* Concept */
    {
        typedef Triangulation_ds_vertex<TDS2> Other;
    };
};

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_TRIANGULATION_DS_VERTEX_H
