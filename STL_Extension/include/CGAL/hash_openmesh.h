// Copyright (c) 2016  GeometryFactory (France).  All rights reserved.
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
// 
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_HASH_OPENMESH_H
#define CGAL_HASH_OPENMESH_H

#include <OpenMesh/Core/Mesh/Handles.hh>

#if OM_VERSION < 0x60200

namespace OpenMesh {

inline std::size_t hash_value(const BaseHandle& h) { return h.idx(); }

} // namespace OpenMesh
#endif


#ifndef OM_HAS_HASH

#include <functional>


namespace std {

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4099) // For VC10 it is class hash 
#endif

#ifndef CGAL_CFG_NO_STD_HASH

template <>
struct hash<OpenMesh::BaseHandle >
  : public std::unary_function<OpenMesh::BaseHandle, std::size_t>
{

  std::size_t operator()(const OpenMesh::BaseHandle& h) const
  {
    return h.idx();
  }
};

template <>
struct hash<OpenMesh::VertexHandle >
  : public std::unary_function<OpenMesh::VertexHandle, std::size_t>
{

  std::size_t operator()(const OpenMesh::VertexHandle& h) const
  {
    return h.idx();
  }
};

template <>
struct hash<OpenMesh::HalfedgeHandle >
  : public std::unary_function<OpenMesh::HalfedgeHandle, std::size_t>
{

  std::size_t operator()(const OpenMesh::HalfedgeHandle& h) const
  {
    return h.idx();
  }
};

template <>
struct hash<OpenMesh::EdgeHandle >
  : public std::unary_function<OpenMesh::EdgeHandle, std::size_t>
{

  std::size_t operator()(const OpenMesh::EdgeHandle& h) const
  {
    return h.idx();
  }
};

template <>
struct hash<OpenMesh::FaceHandle >
  : public std::unary_function<OpenMesh::FaceHandle, std::size_t>
{

  std::size_t operator()(const OpenMesh::FaceHandle& h) const
  {
    return h.idx();
  }
};

#endif // CGAL_CFG_NO_STD_HASH

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

} // namespace std


#endif  // OM_HAS_HASH

#endif // CGAL_HASH_OPENMESH_H
