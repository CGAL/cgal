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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Pierre Angelier, Michel Pocchiola

#ifndef VISIBILITY_COMPLEX_SCENE_H
#define VISIBILITY_COMPLEX_SCENE_H

#include <list>

CGAL_BEGIN_NAMESPACE

// -----------------------------------------------------------------------------

template < class V , class I >
struct Dereference_iterator : public I
{
    typedef V                      value_type;
    typedef value_type&            reference;
    typedef typename I::value_type pointer;
    Dereference_iterator( I j ) : I(j) { }

    reference  operator*()  const { return *(I::operator*()); }
    pointer    operator->() const { return   I::operator*();  }
};

// -----------------------------------------------------------------------------

template < class Vc >
class Visibility_complex_scene
{
public:
    typedef typename Vc::Disk                            Disk;
    typedef typename Vc::Disk_handle                     Disk_handle;
    typedef std::list<Disk_handle>::const_iterator       Disk_handle_it;
    typedef Dereference_iterator<Disk,Disk_handle_it>    Disk_iterator;
    typedef typename Vc::Vertex                          Vertex;
    typedef std::list<Vertex>::const_iterator            Vertex_iterator;
public:
    void push_back(const Disk& d , bool is_convex = false);
    void push_back(const Vertex& v) { constraints.push_back(v); }
    
    Disk_iterator disks_begin() const { return Disk_iterator(disks.begin()); }
    Disk_iterator disks_end()   const { return Disk_iterator(disks.end());   }
    Vertex_iterator constraints_begin() const { return constraints.begin(); }
    Vertex_iterator constraints_end()   const { return constraints.end();   }
private:
    std::list<Disk_handle> disks;
    std::list<Vertex>        constraints;
};

// -----------------------------------------------------------------------------

template < class Gt >
void
Visibility_complex_scene<Gt>::push_back(const Disk& d , bool is_convex)
{
    if (is_convex) { disks.push_back(&d); return; }

}

// -----------------------------------------------------------------------------

CGAL_END_NAMESPACE

#endif // VISIBILITY_COMPLEX_SCENE_H
