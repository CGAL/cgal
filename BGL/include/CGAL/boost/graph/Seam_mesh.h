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

#ifndef CGAL_SEAM_MESH_H
#define CGAL_SEAM_MESH_H

#include <boost/unordered_set.hpp>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/graph_traits_Seam_mesh.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>

namespace CGAL {

template <typename HD>
struct Seam_mesh_halfedge_descriptor
{
  typedef HD                                      TM_halfedge_descriptor;

  TM_halfedge_descriptor tmhd;
  bool seam;

  Seam_mesh_halfedge_descriptor()
    : tmhd(), seam(false)
  { }

  Seam_mesh_halfedge_descriptor(const Seam_mesh_halfedge_descriptor& other)
    : tmhd(other.tmhd), seam(other.seam)
  { }

  Seam_mesh_halfedge_descriptor(TM_halfedge_descriptor tmhd, bool seam=false)
    : tmhd(tmhd),seam(seam)
  { }

  bool operator==(const Seam_mesh_halfedge_descriptor& other) const
  {
    return (tmhd == other.tmhd) && (seam == other.seam);
  }

  bool operator!=(const Seam_mesh_halfedge_descriptor& other) const
  {
    return (tmhd != other.tmhd) || (seam != other.seam);
  }

  bool operator<(const Seam_mesh_halfedge_descriptor& other) const
  {
    return tmhd < other.tmhd;
  }

  operator TM_halfedge_descriptor() const
  {
    return tmhd;
  }

  friend
  std::ostream& operator<<(std::ostream& os, const Seam_mesh_halfedge_descriptor& hd)
  {
    os << hd.tmhd  << ((hd.seam)?" on seam":"");
    return os;
  }
};

/*!
\ingroup PkgBGLHelper

An adaptor for a triangle mesh which turns some marked edges into
boundary edges when exploring the seam mesh with the generic BGL style
functions.

\cgalModels `FaceGraph`

\tparam TM a model of `FaceGraph`
\tparam SEM a model of `ReadablePropertyMap` with `boost::graph_traits<TM>::%edge_descriptor` as key type and `bool` as value type.
\tparam SVM a model of `ReadablePropertyMap` with `boost::graph_traits<TM>::%vertex_descriptor` as key type and `bool` as value type.

\sa \link BGLSeam_meshGT `boost::graph_traits<Seam_mesh<TM> >` \endlink

*/

template <class TM, class SEM, class SVM>
class Seam_mesh
{
  typedef Seam_mesh<TM, SEM, SVM>                                 Self;

  typedef typename boost::graph_traits<TM>::halfedge_descriptor   TM_halfedge_descriptor;
  typedef typename boost::graph_traits<TM>::halfedge_iterator     TM_halfedge_iterator;
  typedef typename boost::graph_traits<TM>::edge_descriptor       TM_edge_descriptor;
  typedef typename boost::graph_traits<TM>::vertex_descriptor     TM_vertex_descriptor;

public:
/// @cond CGAL_DOCUMENT_INTERNALS
  typedef Seam_mesh_halfedge_descriptor<TM_halfedge_descriptor>   halfedge_descriptor;
  typedef typename boost::graph_traits<TM>::face_descriptor       face_descriptor;

  const TM& mesh()const
  {
    return tm;
  }

  ///  A vertex
  struct vertex_descriptor
  {
    halfedge_descriptor hd;

    vertex_descriptor() { }

    vertex_descriptor(const halfedge_descriptor& h)
      : hd(h)
    { }

    vertex_descriptor(const vertex_descriptor& other)
      : hd(other.hd)
    { }

    bool operator ==(const vertex_descriptor& other) const
    {
      return (hd == other.hd);
    }

    bool operator !=(const vertex_descriptor& other) const
    {
      return (hd != other.hd);
    }

    bool operator<(const vertex_descriptor& other) const
    {
      return hd < other.hd;
    }

    operator TM_halfedge_descriptor() const
    {
      return hd;
    }

    friend std::ostream& operator<<(std::ostream& os, const vertex_descriptor vd)
    {
      os << "seam mesh vertex: " <<  vd.hd;
      return os;
    }

    friend std::size_t hash_value(const vertex_descriptor&  vd)
    {
      return hash_value(vd.hd.tmhd);
    }
  };

  class vertex_iterator
    : public boost::iterator_facade<vertex_iterator,
                                    vertex_descriptor,
                                    std::forward_iterator_tag,
                                    vertex_descriptor>
  {
    typedef boost::iterator_facade<vertex_iterator,
                                    vertex_descriptor,
                                    std::forward_iterator_tag,
                                    vertex_descriptor>                  Facade;

  public:
    vertex_iterator() : hd(), end(), mesh_(NULL) { }

    vertex_iterator(const Iterator_range<TM_halfedge_iterator>& ir, const Self* m)
      : hd(ir.first), end(ir.second), mesh_(m)
    {
      //std::cerr << "vertex_iterator(..)\n";
      //std::cerr << *hd << std::endl;
      if(hd == end)
        return;

      TM_vertex_descriptor tvd = target(*hd, mesh_->mesh());
      if( (!mesh_->has_on_seam(tvd)) && (halfedge(tvd, mesh_->mesh()) == *hd))
        return;
      if(mesh_->has_on_seam(edge(*hd, mesh_->mesh())))
        return;
      if(mesh_->has_on_seam(tvd) && is_border(opposite(*hd, mesh_->mesh()), mesh_->mesh()))
        return;
      increment();
      //std::cerr << *hd << "  after increment" << std::endl;
      //std::cerr << "leave vertex_iterator(..)\n";
    }

    // constructor for the past the end iterator
    vertex_iterator(const TM_halfedge_iterator& hd, const Self* m)
      : hd(hd), end(hd), mesh_(m)
    { }

  private:
    friend class boost::iterator_core_access;

    void increment()
    {
      //std::cerr << "increment\n";
      if(hd == end)
        return;

      do {
        ++hd;
        //std::cerr << *hd << "  ++" << std::endl;
        if(hd == end) return;
        TM_vertex_descriptor tvd = target(*hd, mesh_->mesh());
        //std::cerr << "tvd = " << tvd << std::endl;
        if( (!mesh_->has_on_seam(tvd))&& (halfedge(tvd, mesh_->mesh()) == *hd)) {
          //std::cerr <<"return as not on seam and reverse incidence\n";
          return;
        }

        if(mesh_->has_on_seam(edge(*hd, mesh_->mesh()))) {
          //std::cerr <<"return as edge on seam\n";
          return;
        }

        if(mesh_->has_on_seam(tvd) && is_border(opposite(*hd, mesh_->mesh()), mesh_->mesh())) {
          //std::cerr <<"return as edge on border and target on seam\n";
          return;
        }
      } while(true);
    }

    bool equal(const vertex_iterator& other) const
    {
      return (this->mesh_ == other.mesh_) && (this->hd == other.hd);
    }

    vertex_descriptor dereference() const
    {
      return vertex_descriptor(*hd);
    }

    TM_halfedge_iterator hd, end;
    const Self* mesh_;
  };

  bool has_on_seam(TM_vertex_descriptor vd) const
  {
    return get(svm, vd);
  }

  bool has_on_seam(TM_edge_descriptor ed) const
  {
    return get(sem, ed);
  }

  bool has_on_seam(TM_halfedge_descriptor tmhd) const
  {
    return get(sem, edge(tmhd, tm));
  }

  /// returns `true` if the halfedge is on the seam.
  bool has_on_seam(const halfedge_descriptor& hd) const
  {
    return has_on_seam(edge(hd, tm));
  }

  Iterator_range<vertex_iterator> m_vertices() const
  {
    Iterator_range<TM_halfedge_iterator> ir = halfedges(tm);
    vertex_iterator beg(ir, this);
    vertex_iterator end(ir.second, this);
    return make_range(beg,end);
  }

  struct edge_descriptor
  {
    halfedge_descriptor hd;

    edge_descriptor(const halfedge_descriptor& hd)
      : hd(hd)
    { }
  };

  const TM& tm;
  SEM sem;
  SVM svm;
  int index;

  /// @endcond

public:

  /// Constructs a seam mesh for a triangle mesh and an edge and vertex property map
  /// \param tm the adapted mesh
  /// \param sem the edge property map with value `true` for seam edges
  /// \param svm the vertex property map with value `true` for seam vertices

  /// @note the vertices must be exactly the vertices on the seam edges. Maybe a bad design.
  Seam_mesh(const TM& tm, const SEM& sem, const SVM& svm)
    : tm(tm), sem(sem), svm(svm), index(0)
  { }

  /// Puts in `vim` indices of vertices in the connected component with the boundary on which lies `bhd`.
  /// Indices are consecutive and start at 0 or the last index+1 that was assigned in a previous call.
  /// That value is cached and returned when calling `num_vertices(sm)`.
  /// \tparam VertexIndexMap a model of `ReadWritePropertyMap` with `vertex_descriptor` as key and `int` as value-type
  template <typename VertexIndexMap>
  void initialize_vertex_index_map(halfedge_descriptor bhd, VertexIndexMap& vim)
  {
    Self& mesh = *this;

    index = 0;
    std::vector<face_descriptor> faces;
    typename boost::graph_traits<Seam_mesh>::halfedge_descriptor shd(opposite(bhd, *this));
    CGAL::Polygon_mesh_processing::connected_component(face(shd,*this),
                                                       *this,
                                                       std::back_inserter(faces));

    BOOST_FOREACH(face_descriptor fd, faces) {
      BOOST_FOREACH(TM_halfedge_descriptor tmhd , halfedges_around_face(halfedge(fd, tm), tm)) {
        halfedge_descriptor hd(tmhd);
        vertex_descriptor vd = target(hd, mesh);
        put(vim,vd,-1);
      }
    }

    BOOST_FOREACH(face_descriptor fd, faces) {
      BOOST_FOREACH(TM_halfedge_descriptor tmhd , halfedges_around_face(halfedge(fd, tm), tm)) {
        halfedge_descriptor hd(tmhd);
        vertex_descriptor vd = target(hd, mesh);
        if(get(vim,vd) == -1) {
          put(vim,vd,index);
          ++index;
        }
      }
    }
  }

/// @cond CGAL_DOCUMENT_INTERNALS
  // this is the number of different halfedge indices
  int m_num_vertices() const
  {
    return index;
  }

  halfedge_descriptor m_next(const halfedge_descriptor& hd) const
  {
    if((!hd.seam) && (!is_border(hd.tmhd, tm))) {
      return halfedge_descriptor(next(hd.tmhd, tm));
    }
    Halfedge_around_target_circulator<TM> hatc(hd.tmhd, tm);
    do {
      --hatc;
    } while((!has_on_seam(*hatc))&&(!is_border(opposite(*hatc, tm), tm)));

    return halfedge_descriptor(opposite(*hatc, tm), !is_border(opposite(*hatc, tm), tm));
  }

  halfedge_descriptor m_prev(const halfedge_descriptor& hd) const
  {
    if((!hd.seam)&& (!is_border(hd.tmhd, tm))) {
      return halfedge_descriptor(prev(hd.tmhd, tm));
    }
    Halfedge_around_source_circulator<TM> hatc(hd.tmhd, tm);
    do {
      ++hatc;
    } while((!has_on_seam(*hatc))&&(!is_border(opposite(*hatc, tm), tm)));

    return halfedge_descriptor(opposite(*hatc, tm), !is_border(opposite(*hatc, tm), tm));
  }

  halfedge_descriptor m_opposite(const halfedge_descriptor& hd) const
  {
    if(!hd.seam) {
      return halfedge_descriptor(opposite(hd.tmhd, tm), has_on_seam(hd));
    }

    return halfedge_descriptor(opposite(hd.tmhd, tm));
  }

  vertex_descriptor m_target(halfedge_descriptor hd) const
  {
    TM_halfedge_descriptor tmhd(hd);
    if(!has_on_seam(target(tmhd, tm))) {
      tmhd = halfedge(target(tmhd, tm), tm);
      return vertex_descriptor(halfedge_descriptor(tmhd));
    }

    if(hd.seam) {
      return m_target(halfedge_descriptor(prev(opposite(tmhd, tm), tm)));
    }

    while((!has_on_seam(tmhd)) && (!is_border(opposite(tmhd, tm), tm))) {
      tmhd = prev(opposite(tmhd, tm), tm);
    }

    return vertex_descriptor(halfedge_descriptor(tmhd));
  }

  vertex_descriptor m_source(const halfedge_descriptor& hd) const
  {
    return m_target(opposite(hd.tmhd, tm));
  }

  /// @endcond
};

} // namespace CGAL

#endif //CGAL_SEAM_MESH_H
