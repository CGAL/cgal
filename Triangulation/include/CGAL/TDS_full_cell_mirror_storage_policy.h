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

#ifndef CGAL_TDS_FULL_CELL_MIRROR_STORAGE_POLICY_H
#define CGAL_TDS_FULL_CELL_MIRROR_STORAGE_POLICY_H

#include <CGAL/TDS_full_cell_default_storage_policy.h>

namespace CGAL {

// POLICY TAGS

struct TDS_full_cell_mirror_storage_policy {}; // Stores the mirror index of all vertices.

template< typename Vertex_handle, typename Full_cell_handle, typename Ambient_dimension >
struct   TFC_data< Vertex_handle, Full_cell_handle, Ambient_dimension, TDS_full_cell_mirror_storage_policy >
: public TFC_data< Vertex_handle, Full_cell_handle, Ambient_dimension, TDS_full_cell_default_storage_policy >
{
    typedef TFC_data< Vertex_handle, Full_cell_handle, Ambient_dimension, TDS_full_cell_default_storage_policy > Base;
    typedef typename Base::Vertex_handle_array          Vertex_handle_array;
    typedef typename Base::Full_cell_handle_array         Full_cell_handle_array;
    typedef typename internal::S_or_D_array< int, typename Base::Dimen_plus >   Int_array;

private:
    Int_array            mirror_vertices_;

public:
    TFC_data(const int dmax)
    : Base(dmax), mirror_vertices_(dmax+1)
    {}

    void set_mirror_index(const int i, const int index) 
    {
        mirror_vertices_[i] = index;
    }
    int mirror_index(const int i) const
    {
        return mirror_vertices_[i];
    }
    Vertex_handle mirror_vertex(const int i, const int) const
    {
        return Base::neighbors_[i]->vertex(mirror_index(i));
    }
    // FIXME: rename to switch_vertices
    void swap_vertices(const int d1, const int d2)
    {
        Base::swap_vertices(d1, d2);
        std::swap(mirror_vertices_[d1], mirror_vertices_[d2]);
        Base::neighbors_[d1]->set_mirror_index(mirror_vertices_[d1], d1);
        Base::neighbors_[d2]->set_mirror_index(mirror_vertices_[d2], d2);
    }
};

} //namespace CGAL

#endif // CGAL_TDS_FULL_CELL_MIRROR_STORAGE_POLICY_H
