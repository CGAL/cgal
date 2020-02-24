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

#ifndef CGAL_MESH_3_GET_FACET_PATCH_ID_H
#define CGAL_MESH_3_GET_FACET_PATCH_ID_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/property_map.h>
#include <CGAL/boost/graph/properties.h>
#include <boost/mpl/has_xxx.hpp>

namespace CGAL { namespace Mesh_3 {

namespace internal {
BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Primitive_has_Id, Id, false)
} // end namespace CGAL::Mesh_3::internal

/// A property map that, from a primitive of a AABB tree of polyhedron
/// facets, retrieve the patch_id() of the facet.
template <typename Primitive>
struct Get_facet_patch_id{};


// generic version, that test if Primitive::Id exists
template <typename Primitive,
          bool = internal::Primitive_has_Id<Primitive>::value >
struct Get_facet_patch_id_property_traits {
};

// specialization when Primitive::Id exists
template <typename Primitive>
struct Get_facet_patch_id_property_traits<Primitive, true>
{
  typedef typename Primitive::Id Id;
  typedef typename std::iterator_traits<Id>::value_type Face;
  typedef typename Face::Patch_id value_type;
  typedef value_type& reference;
  typedef Primitive key_type;
  typedef boost::readable_property_map_tag category;
};

}} // end namespace CGAL::Mesh_3

namespace boost {
  // specialization for using pointers as property maps
  template <typename Primitive>
  struct property_traits<CGAL::Mesh_3::Get_facet_patch_id<Primitive> >
    : public CGAL::Mesh_3::Get_facet_patch_id_property_traits<Primitive> {};
}

namespace CGAL { namespace Mesh_3 {

template <typename Primitive>
typename boost::property_traits< Get_facet_patch_id<Primitive> >::value_type
get(const Get_facet_patch_id<Primitive>, const typename Primitive::Id& primitive_id) {
  return primitive_id->patch_id();
}



/// A property map that, from a primitive of a AABB tree of Gwdwg<Surface_mesh>
/// facets, retrieve the patch_id() of the facet.
template < typename MeshDomain>
struct Get_facet_patch_id_sm{};


// generic version, that test if Primitive::Id exists
template <typename MeshDomain,
          bool = internal::Primitive_has_Id<typename MeshDomain::AABB_primitive>::value>
struct Get_facet_patch_id_sm_property_traits {
};

// specialization when Primitive::Id exists
template <typename MeshDomain>
struct Get_facet_patch_id_sm_property_traits<MeshDomain, true>
{
  typedef typename MeshDomain::AABB_primitive::Id Id;
  typedef typename MeshDomain::Patch_id value_type;

  typedef value_type& reference;
  typedef typename MeshDomain::AABB_primitive key_type;
  typedef boost::readable_property_map_tag category;
};

}} // end namespace CGAL::Mesh_3

namespace boost {
  // specialization for using pointers as property maps
  template <typename MeshDomain>
  struct property_traits<CGAL::Mesh_3::Get_facet_patch_id_sm<MeshDomain> >
    : public CGAL::Mesh_3::Get_facet_patch_id_sm_property_traits<MeshDomain> {};
}

namespace CGAL { namespace Mesh_3 {

template <typename MeshDomain>
typename boost::property_traits< Get_facet_patch_id_sm<MeshDomain> >::value_type
get(const Get_facet_patch_id_sm<MeshDomain>,
    const typename MeshDomain::AABB_primitive::Id& primitive_id)
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
}} // end namespace CGAL::Mesh_3

#endif // CGAL_MESH_3_GET_FACET_PATCH_ID_H
