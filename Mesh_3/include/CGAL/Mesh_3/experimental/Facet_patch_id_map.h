// Copyright (c) 2010,2012  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau
//

#ifndef CGAL_MESH_3_FACET_PATCH_ID_MAP_H
#define CGAL_MESH_3_FACET_PATCH_ID_MAP_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/property_map.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/is_iterator.h>

namespace CGAL { namespace Mesh_3 {

/// A property map that, from a primitive of a AABB tree
/// retrieve the patch_id() of the facet.
template <typename MeshDomain,
          typename Primitive,
          bool id_is_iterator = CGAL::is_iterator<typename Primitive::Id>::value >
struct Facet_patch_id_map;

// Primitive::Id is an iterator type
template <typename MeshDomain,
          typename Primitive>
struct Facet_patch_id_map<MeshDomain, Primitive, true>
{
  typedef typename Primitive::Id Id;
  typedef typename std::iterator_traits<Id>::value_type Face;
  typedef typename Face::Patch_id value_type;
  typedef const value_type& reference;
  typedef typename Primitive::Id key_type;
  typedef boost::readable_property_map_tag category;

  friend reference get(Facet_patch_id_map<MeshDomain, Primitive, true>, key_type primitive_id)
  {
    return primitive_id->patch_id();
  }
};

// Primitive::Id is a std::pair
template <typename MeshDomain,
          typename Primitive>
struct Facet_patch_id_map<MeshDomain, Primitive, false>
{
  typedef typename MeshDomain::AABB_primitive::Id Id;
  typedef typename MeshDomain::Patch_id value_type;
  typedef value_type reference;
  typedef typename MeshDomain::AABB_primitive::Id key_type;
  typedef boost::readable_property_map_tag category;

  friend reference get(Facet_patch_id_map<MeshDomain, Primitive, false>, key_type primitive_id)
  {
    typedef typename boost::property_map<
         typename MeshDomain::Polyhedron,
         face_patch_id_t<typename MeshDomain::Patch_id> >::type Fpim;
    Fpim fpim = get(face_patch_id_t<typename MeshDomain::Patch_id>(),
                    *(primitive_id.second));
    typename MeshDomain::Patch_id patch_index = get(fpim,
                                                    primitive_id.first);
    return patch_index;
  }
};

}} // end namespace CGAL::Mesh_3

#endif // CGAL_MESH_3_FACET_PATCH_ID_MAP_H
