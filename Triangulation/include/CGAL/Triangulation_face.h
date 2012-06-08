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

#ifndef CGAL_TRIANGULATION_FACE_H
#define CGAL_TRIANGULATION_FACE_H

#include <CGAL/basic.h>
#include <CGAL/internal/Static_or_dynamic_array.h>

namespace CGAL {

template< typename TDS >
class Triangulation_face
{
    typedef typename internal::Dimen_plus_one<typename TDS::Ambient_dimension>::type Dimen_plus;
public:
    typedef typename TDS::Full_cell_handle          Full_cell_handle; /* Concept */
    typedef typename TDS::Vertex_handle             Vertex_handle; /* Concept */
    typedef internal::S_or_D_array<int, Dimen_plus> Indices;

protected:
    Full_cell_handle        full_cell_;
    Indices                 indices_;

public:
    explicit Triangulation_face(Full_cell_handle s) /* Concept */
    : full_cell_(s), indices_(s->ambient_dimension()+2) 
    {
        CGAL_assertion( Full_cell_handle() != s );
        clear();
    }

    explicit Triangulation_face(const int ambient_dim) /* Concept */
    : full_cell_(), indices_(ambient_dim+2) 
    {
        clear();
    }

    Triangulation_face(const Triangulation_face & f) /* Concept */
    : full_cell_(f.full_cell_), indices_(f.indices_)
    {}

    int face_dimension() const /* Concept */
    {
        int i(0);
        while( -1 != indices_[i] ) ++i;
        return (i-1);
    }

    Full_cell_handle full_cell() const /* Concept */
    {
        return full_cell_;
    }

    int index(const int i) const /* Concept */
    {
        CGAL_precondition( (0 <= i) && (i <= face_dimension()) );
        return indices_[i];
    }

    Vertex_handle vertex(const int i) const /* Concept */
    {
        int j = index(i);
        if( j == -1 )
            return Vertex_handle();
        return full_cell()->vertex(j);
    }

// - - - - - - - - - - - - - - - - - -  UPDATE FUNCTIONS	

    void clear() /* Concept */
    {    
        const int d = indices_.size();
        for(int i = 0; i < d; ++i )
            indices_[i] = -1;
    }

    void set_full_cell(Full_cell_handle s) /* Concept */
    {
        CGAL_precondition( Full_cell_handle() != s );
        full_cell_ = s;
    }

    void set_index(const int i, const int idx) /* Concept */
    {
        CGAL_precondition( (0 <= i) && ((size_t)i+1 < indices_.size()) );
        CGAL_precondition( (0 <= idx) && ((size_t)idx < indices_.size()) );
        indices_[i] = idx;
    }
};

} //namespace CGAL

#endif // CGAL_TRIANGULATION_FACE_H
