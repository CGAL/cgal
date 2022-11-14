// Copyright (c) 2016  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_HASH_OPENMESH_H
#define CGAL_HASH_OPENMESH_H

#include <OpenMesh/Core/Mesh/Handles.hh>
#include <CGAL/algorithm.h>

namespace CGAL { namespace internal {


template<typename Halfedge_handle>
class OMesh_edge {
public:
  OMesh_edge() : halfedge_() {}
  explicit OMesh_edge(const Halfedge_handle& h) : halfedge_(h) {}
  Halfedge_handle halfedge() const { return halfedge_; }
  bool is_valid() const { return halfedge_.is_valid(); }

  bool
  operator==(const OMesh_edge& other) const {
    if(halfedge_ == other.halfedge_) {
      return true;
    } else if(halfedge_ != Halfedge_handle()) {
      return opposite() == other.halfedge_;
    } else {
      return false;
    }
  }

  bool operator<(const OMesh_edge& other) const
  {
    return this->idx() < other.idx();
  }

  bool
  operator!=(const OMesh_edge& other) const
  {
    return !(*this == other);
  }

  Halfedge_handle
  opposite() const { return Halfedge_handle((halfedge_.idx() & 1) ? halfedge_.idx()-1 : halfedge_.idx()+1); }

  OMesh_edge
  opposite_edge() const { return OMesh_edge(Halfedge_handle((halfedge_.idx() & 1) ? halfedge_.idx()-1 : halfedge_.idx()+1)); }

  unsigned int idx() const { return halfedge_.idx() / 2; }
private:
  Halfedge_handle halfedge_;
};

} } // CGAL::internal

#if OM_VERSION < 0x60200

namespace OpenMesh {

inline std::size_t hash_value(const BaseHandle& h) { return h.idx(); }

} // namespace OpenMesh
#endif

namespace OpenMesh {

inline std::size_t hash_value(const CGAL::internal::OMesh_edge<OpenMesh::HalfedgeHandle>& h) { return h.idx(); }

} // namespace OpenMesh



#include <functional>


namespace std {

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4099) // For VC10 it is class hash
#endif

#ifndef CGAL_CFG_NO_STD_HASH

#ifndef OM_HAS_HASH

template <>
struct hash<OpenMesh::BaseHandle >
  : public CGAL::cpp98::unary_function<OpenMesh::BaseHandle, std::size_t>
{

  std::size_t operator()(const OpenMesh::BaseHandle& h) const
  {
    return h.idx();
  }
};

template <>
struct hash<OpenMesh::VertexHandle >
  : public CGAL::cpp98::unary_function<OpenMesh::VertexHandle, std::size_t>
{

  std::size_t operator()(const OpenMesh::VertexHandle& h) const
  {
    return h.idx();
  }
};

template <>
struct hash<OpenMesh::HalfedgeHandle >
  : public CGAL::cpp98::unary_function<OpenMesh::HalfedgeHandle, std::size_t>
{

  std::size_t operator()(const OpenMesh::HalfedgeHandle& h) const
  {
    return h.idx();
  }
};

template <>
struct hash<OpenMesh::EdgeHandle >
  : public CGAL::cpp98::unary_function<OpenMesh::EdgeHandle, std::size_t>
{

  std::size_t operator()(const OpenMesh::EdgeHandle& h) const
  {
    return h.idx();
  }
};


template <>
struct hash<OpenMesh::FaceHandle >
  : public CGAL::cpp98::unary_function<OpenMesh::FaceHandle, std::size_t>
{

  std::size_t operator()(const OpenMesh::FaceHandle& h) const
  {
    return h.idx();
  }
};

#endif  // OM_HAS_HASH

template <typename H>
struct hash<CGAL::internal::OMesh_edge<H> >
  : public CGAL::cpp98::unary_function<CGAL::internal::OMesh_edge<H>, std::size_t>
{

  std::size_t operator()(const CGAL::internal::OMesh_edge<H>& h) const
  {
    return h.idx();
  }
};


#endif // CGAL_CFG_NO_STD_HASH

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

} // namespace std




#endif // CGAL_HASH_OPENMESH_H
