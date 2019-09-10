// Copyright (c) 2007  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Andreas Fabri, Fernando Cacciola

#ifndef CGAL_HALFEDGEDS_FACE_MAX_BASE_WITH_ID_H
#define CGAL_HALFEDGEDS_FACE_MAX_BASE_WITH_ID_H 1

#include <CGAL/HalfedgeDS_face_base.h>

namespace CGAL {

template < class Refs, class Pln, class ID>
class HalfedgeDS_face_max_base_with_id : public HalfedgeDS_face_base< Refs, Tag_true, Pln>
{
public:
    
    typedef HalfedgeDS_face_base< Refs, Tag_true, Pln> Base ;
    
    typedef ID size_type ;
    
private:

    size_type mID ;
    
public:

    HalfedgeDS_face_max_base_with_id() : mID ( size_type(-1) ) {}
    HalfedgeDS_face_max_base_with_id( Pln const& p) : Base(p), mID ( size_type(-1) ) {}
    HalfedgeDS_face_max_base_with_id( Pln const& p, size_type i ) : Base(p), mID (i) {}
    
    size_type&       id()       { return mID; }
    size_type const& id() const { return mID; }
};

} //namespace CGAL

#endif // CGAL_HALFEDGEDS_FACE_MAX_BASE_WITH_ID_H //
// EOF //
