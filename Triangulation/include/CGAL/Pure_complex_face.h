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

#ifndef CGAL_PURE_COMPLEX_FACE_H
#define CGAL_PURE_COMPLEX_FACE_H

#include <CGAL/basic.h>
#include <CGAL/internal/Static_or_dynamic_array.h>

namespace CGAL {

template< typename PCDS >
class Pure_complex_face
{
    typedef typename internal::Dimen_plus_one<typename PCDS::Ambient_dimension>::type Dimen_plus;
public:
    typedef typename PCDS::Simplex_handle           Simplex_handle;
    typedef typename PCDS::Vertex_handle            Vertex_handle;
    typedef internal::S_or_D_array<int, Dimen_plus> Indices;

protected:
    Simplex_handle              simplex_;
    Indices                     indices_;

public:
    explicit Pure_complex_face(Simplex_handle s)
    : simplex_(s), indices_(s->ambient_dimension()+1)
    {
        CGAL_assertion( Simplex_handle() != s );
        clear();
    }

    explicit Pure_complex_face(const int ambient_dim)
    : simplex_(), indices_(ambient_dim+1)
    {
        clear();
    }

    Pure_complex_face(const Pure_complex_face & f)
    : simplex_(f.simplex_), indices_(f.indices_)
    {}

    int feature_dimension() const
    {
        int i(0);
        while( -1 != indices_[i] ) ++i;
        return (i-1);
    }

    Simplex_handle simplex() const
    {
        return simplex_;
    }

    int index(const int i) const
    {
        CGAL_precondition( (0 <= i) && (i <= feature_dimension()) );
        return indices_[i];
    }

    Vertex_handle vertex(const int i) const
    {
        int j = index(i);
        if( j == -1 )
            return Vertex_handle();
        return simplex()->vertex(j);
    }

// - - - - - - - - - - - - - - - - - -  UPDATE FUNCTIONS	

    void clear()
    {    
        const int d = indices_.size();
        for(int i = 0; i < d; ++i )
            indices_[i] = -1;
    }

    void set_simplex(Simplex_handle s)
    {
        CGAL_precondition( Simplex_handle() != s );
        simplex_ = s;
    }

    void set_index(const int i, const int idx)
    {
        CGAL_precondition( (0 <= idx) && ((size_t)idx < indices_.size()) );
        indices_[i] = idx;
    }
};

} //namespace CGAL

#endif // CGAL_PURE_COMPLEX_FACE_H
