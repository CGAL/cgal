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

#ifndef CGAL_MESH_3_GET_CURVE_INDEX_H
#define CGAL_MESH_3_GET_CURVE_INDEX_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/property_map.h>

namespace CGAL { namespace Mesh_3 {

// template <typename Primitive, typename Mesh_domain>
// struct Get_curves_indices {
//   const Mesh_domain* mesh_domain;

//   Get_curves_indices(const Mesh_domain* mesh_domain)
//     : mesh_domain(mesh_domain) {}

// }; // end Get_curves_indices

// namespace boost {
//   // specialization for using pointers as property maps
//   template <typename Primitive, typename Mesh_domain>
//   struct property_traits<Get_curves_indices<Primitive, Mesh_domain> > {
//     typedef std::set<typename Mesh_domain::Surface_patch_index> value_type;
//     typedef Primitive key_type;
//     typedef readable_property_map_tag category;
//   };
// }

// template <typename Primitive, typename Mesh_domain>
// const std::set<typename Mesh_domain::Surface_patch_index>&
// get(Get_curves_indices<Primitive, Mesh_domain>& get_curves_indices, const Primitive& primitive) {
//   return get_curves_indices.mesh_domain->get_incidences(primitive.id().first->first);
// }

template <typename Primitive>
struct Get_curve_index {
}; // end Get_curve_index

template <typename Primitive>
typename boost::property_traits<Get_curve_index<Primitive> >::value_type
get(Get_curve_index<Primitive>, const typename Primitive::Id id) {
  return id.first->first;
}

}} // end namespace CGAL::Mesh_3

#include <boost/type_traits/remove_const.hpp>

namespace boost {
  // specialization for using pointers as property maps
  template <typename Primitive>
  struct property_traits<CGAL::Mesh_3::Get_curve_index<Primitive> > {
    typedef typename std::iterator_traits<typename Primitive::Id::first_type>::value_type ConstPair;
    typedef typename boost::remove_const<typename ConstPair::first_type>::type value_type;
    typedef value_type& reference;
    typedef typename Primitive::Id key_type;
    typedef readable_property_map_tag category;
  };
}


#endif // CGAL_MESH_3_GET_CURVE_INDEX_H
