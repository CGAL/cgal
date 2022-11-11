// Copyright (c) 1997-2000  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Peter Hachenberger <hachenberger@mpi-sb.mpg.de>

#ifndef CGAL_SNC_LIST_H
#define CGAL_SNC_LIST_H

#include <CGAL/license/Nef_3.h>


#include <CGAL/Nef_S2/SM_list.h>

namespace CGAL {

template < class Sphere_map>
class SNC_in_place_list_sm
    : public Sphere_map,
      public In_place_list_base<SNC_in_place_list_sm<Sphere_map> > {
public:
    typedef SNC_in_place_list_sm<Sphere_map> Self;

    SNC_in_place_list_sm() {}
    SNC_in_place_list_sm(const Self&)=default;
    SNC_in_place_list_sm(const Sphere_map& sm)   // down cast
        : Sphere_map(sm) {}
    Self& operator=( const Self& sm) {
        // This self written assignment avoids that assigning vertices will
        // overwrite the list linking of the target vertex.
        *((Sphere_map*)this) = ((const Sphere_map&)sm);
        return *this;
    }
};

template < class Halffacet>
class SNC_in_place_list_halffacet
    : public Halffacet,
      public In_place_list_base<SNC_in_place_list_halffacet<Halffacet> > {
public:
    typedef SNC_in_place_list_halffacet<Halffacet> Self;

    SNC_in_place_list_halffacet() {}
    SNC_in_place_list_halffacet(const Halffacet& v)   // down cast
        : Halffacet(v) {}
    SNC_in_place_list_halffacet(const Self&)=default;
    Self& operator=( const Self& v) {
        // This self written assignment avoids that assigning vertices will
        // overwrite the list linking of the target vertex.
        *((Halffacet*)this) = ((const Halffacet&)v);
        return *this;
    }
};

template < class Volume>
class SNC_in_place_list_volume
    : public Volume,
      public In_place_list_base<SNC_in_place_list_volume<Volume> > {
public:
    typedef SNC_in_place_list_volume<Volume> Self;

    SNC_in_place_list_volume() {}
    SNC_in_place_list_volume(const Volume& v)   // down cast
        : Volume(v) {}
    SNC_in_place_list_volume(const Self&)=default;
    Self& operator=( const Self& v) {
        // This self written assignment avoids that assigning vertices will
        // overwrite the list linking of the target vertex.
        *((Volume*)this) = ((const Volume&)v);
        return *this;
    }
};

template < class SHalfloop>
class SNC_in_place_list_shalfloop
    : public SHalfloop,
      public In_place_list_base<SNC_in_place_list_shalfloop<SHalfloop> > {
public:
    typedef SNC_in_place_list_shalfloop<SHalfloop> Self;

    SNC_in_place_list_shalfloop() {}
    SNC_in_place_list_shalfloop(const SHalfloop& v)   // down cast
        : SHalfloop(v) {}
    SNC_in_place_list_shalfloop(const Self&)=default;
    Self& operator=( const Self& v) {
        // This self written assignment avoids that assigning vertices will
        // overwrite the list linking of the target vertex.
        *((SHalfloop*)this) = ((const SHalfloop&)v);
        return *this;
    }
};

} //namespace CGAL

#endif // CGAL_SNC_LIST_H
