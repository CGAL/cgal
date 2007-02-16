// Copyright (c) 2007  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
//
// Author(s)     : Andreas Fabri, Fernando Cacciola

#ifndef CGAL_HALFEDGEDS_HALFEDGE_MAX_BASE_WITH_ID_H
#define CGAL_HALFEDGEDS_HALFEDGE_MAX_BASE_WITH_ID_H 1

#include <CGAL/HalfedgeDS_halfedge_base.h>

CGAL_BEGIN_NAMESPACE

template < class Refs, class ID>
class HalfedgeDS_halfedge_max_base_with_id : public HalfedgeDS_halfedge_base< Refs, Tag_true, Tag_true, Tag_true >
{
public:
    typedef HalfedgeDS_halfedge_base< Refs, Tag_true, Tag_true, Tag_true> Base ;
    
    typedef typename Base::Base_base Base_base ;
    
    typedef ID size_type ;
    
private:

    size_type mID ;
    
public:

    HalfedgeDS_halfedge_max_base_with_id( size_type i = size_type(-1) ) : mID(i) {}
    
    size_type&       id()       { return mID; }
    size_type const& id() const { return mID; }
};

CGAL_END_NAMESPACE

#endif // CGAL_HALFEDGEDS_HALFEDGE_MAX_BASE_WITH_ID_H //
// EOF //
