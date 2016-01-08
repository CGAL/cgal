// Copyright (c) 2015 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Jane Tournois

#ifndef CGAL_POLYGON_MESH_PROCESSING_AABB_TREE_REMESHING_H
#define CGAL_POLYGON_MESH_PROCESSING_AABB_TREE_REMESHING_H

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/Polygon_mesh_processing/internal/Isotropic_remeshing/AABB_filtered_projection_traits.h>

#include <CGAL/Polygon_mesh_processing/connected_components.h>

namespace PMP = CGAL::Polygon_mesh_processing;

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

template<typename GeomTraits
       , typename TriangleIterator
       , typename IdType>
class AABB_triangle_with_id_primitive
  : public CGAL::AABB_triangle_primitive<GeomTraits, TriangleIterator>
{
  typedef CGAL::AABB_triangle_primitive<GeomTraits, TriangleIterator> Base;
  typedef typename GeomTraits::Triangle_3 Triangle_3;

public:
  AABB_triangle_with_id_primitive(TriangleIterator it)
    : Base(it)
    , patch_id_(-1)
  {}

  void set_patch_id(const IdType& patch_id) {
    patch_id_ = patch_id;
  }
  IdType patch_id() const {
    return patch_id_;
  }

private:
  std::size_t patch_id_;
};

template<typename GeomTraits,
         typename TriangleList,
         typename PatchIdList>
class AABB_tree_remeshing
  : public CGAL::AABB_tree<
      CGAL::AABB_traits<
        GeomTraits,
        AABB_triangle_with_id_primitive<GeomTraits,
          typename TriangleList::iterator,
          typename PatchIdList::value_type
    > > >
{
  typedef typename TriangleList::iterator     TriangleIterator;
  typedef typename PatchIdList::iterator      PatchIdIterator;
  typedef typename PatchIdList::value_type    IdType;

  typedef AABB_triangle_with_id_primitive<GeomTraits,
            TriangleIterator, IdType>         Primitive;

public:
  typedef CGAL::AABB_traits<GeomTraits, Primitive> AABB_traits;

public:
  AABB_tree_remeshing(TriangleIterator tb, TriangleIterator te,
                      PatchIdIterator pib, PatchIdIterator pie)
  {
    TriangleIterator tr_it = tb;
    PatchIdIterator pid_it = pib;
    for (; tr_it != te; ++tr_it, ++pid_it)
    {
      CGAL_assertion(pid_it != pie);

      Primitive prim(tr_it);
      prim.set_patch_id(*pid_it);
      this->insert(prim);
    }
  }
};

}//namespace internal
}//namespace PMP
}//namespace CGAL

#endif //CGAL_POLYGON_MESH_PROCESSING_AABB_TREE_REMESHING_H