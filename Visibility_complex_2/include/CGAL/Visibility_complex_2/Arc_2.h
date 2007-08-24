// Copyright (c) 2001-2004  ENS of Paris (France).
// All rights reserved.
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
// Author(s)     : Pierre Angelier, Michel Pocchiola

#ifndef CGAL_VISIBILITY_COMPLEX_2_ARC_2_H
#define CGAL_VISIBILITY_COMPLEX_2_ARC_2_H

#include <list>

CGAL_BEGIN_NAMESPACE
namespace Visibility_complex_2_details {
// -----------------------------------------------------------------------------
// Abstract unsigned convex arc

template < class D_ >
class Arc_base 
{
public:
    // -------------------------------------------------------------------------
    typedef D_                Disk;
    typedef const Disk*   Disk_handle;
    // -------------------------------------------------------------------------
public:
    // -------------------------------------------------------------------------
    Arc_base() : object_(0) { }
    Arc_base(Disk_handle P) : object_(P) { }
    // -------------------------------------------------------------------------
    Disk_handle   object()               const { return object_; }
    void             set_object(Disk_handle P) { object_ = P;    }
    // -------------------------------------------------------------------------
private:
    Disk_handle object_;
};

// -----------------------------------------------------------------------------
// --------------------- Arc_2 general definition -----------------------
// ----- This definition is sufficient if arcs have constant complexity --------
// -----------------------------------------------------------------------------

template < class Gtr_>
struct Arc_2 : public Arc_base<typename Gtr_::Disk>
{
    // -------------------------------------------------------------------------
    typedef Gtr_ Gt;
    typedef typename Gt::Disk                Disk;
    typedef typename Gt::R                   R;
    typedef typename R::FT                    FT;
    typedef Arc_base<Disk>                    Base;
    typedef typename Base::Disk_handle        Disk_handle;
    typedef typename Gt::Bitangent_2         Bitangent_2;
    // -------------------------------------------------------------------------
    Arc_2() : Base(0) { }
    Arc_2(Disk_handle P) : Base(P) { }
    Arc_2(Disk_handle P,const Bitangent_2& p, const Bitangent_2& q) 
	: Base(P) { }
    // -------------------------------------------------------------------------
  void split (Arc_2& /*tmp*/, const Bitangent_2& /*p*/) {  } 
  void split_cw(Arc_2& /*tmp*/, const Bitangent_2& /*p*/) { }
  void update_begin(const Bitangent_2& /*p*/) { }
  void update_end(const Bitangent_2& /*p*/) { }
  void join (Arc_2& /*y*/) { }
    // -------------------------------------------------------------------------
};

// -----------------------------------------------------------------------------

}
CGAL_END_NAMESPACE

#endif // CGAL_VISIBILITY_COMPLEX_2_ARC_2_H
