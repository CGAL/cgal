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
// Author(s)     : Andreas Fabri, Mael Rouxel-Labbé

#ifndef CGAL_SEAM_MESH_H
#define CGAL_SEAM_MESH_H

#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/graph_traits_Seam_mesh.h>

#include <CGAL/circulator.h>
#include <CGAL/Unique_hash_map.h>

#include <boost/foreach.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/unordered_set.hpp>

#include <fstream>
#include <iostream>
#include <iterator>
#include <utility>
#include <vector>

namespace CGAL {

#ifndef DOXYGEN_RUNNING
template <typename HD>
class Seam_mesh_halfedge_descriptor
{
public:
  typedef HD                                      TM_halfedge_descriptor;

  TM_halfedge_descriptor tmhd;
  bool seam;

  Seam_mesh_halfedge_descriptor()
    : tmhd(), seam(false)
  { }

  Seam_mesh_halfedge_descriptor(const Seam_mesh_halfedge_descriptor& other)
    : tmhd(other.tmhd), seam(other.seam)
  { }

  Seam_mesh_halfedge_descriptor(TM_halfedge_descriptor tmhd, bool seam = false)
    : tmhd(tmhd), seam(seam)
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
    if(tmhd == other.tmhd)
      return static_cast<int>(seam) < static_cast<int>(other.seam);

    return tmhd < other.tmhd;
  }

  operator TM_halfedge_descriptor() const
  {
    return tmhd;
  }

#ifdef CGAL_SEAM_MESH_INSERT_OPERATOR
  friend
  std::ostream& operator<<(std::ostream& os, const Seam_mesh_halfedge_descriptor& hd)
  {
    os << hd.tmhd << ((hd.seam)?" on seam":"");
    return os;
  }
#endif

  friend std::size_t hash_value(const Seam_mesh_halfedge_descriptor& hd)
  {
    return 2 * hash_value(hd.tmhd) + static_cast<std::size_t>(hd.seam);
  }
};
#endif

/// \ingroup PkgBGLAdaptors
///
/// This class is a data structure that takes a triangle mesh, further refered
/// to as `underlying mesh` and turns some marked edges of that mesh into
/// virtual boundary edges.
///
/// Note that a seam edge must be an edge of the underlying mesh that is not
/// a border edge (that is, such that either halfedge or the edge have null_face()
/// as incident face). Marking a border edge as seam will not do anything.
///
/// \cgalModels `FaceGraph` or `FaceListGraph`, depending on the underlying mesh `TM`.
///
/// \tparam TM a model of `FaceGraph` or `FaceListGraph`
/// \tparam SEM a model of `ReadablePropertyMap` with `boost::graph_traits<TM>::%edge_descriptor` as key type and `bool` as value type.
/// \tparam SVM a model of `ReadablePropertyMap` with `boost::graph_traits<TM>::%vertex_descriptor` as key type and `bool` as value type.
///
/// \sa \link BGLSeam_meshGT `boost::graph_traits<Seam_mesh<TM> >` \endlink
///
template <class TM, class SEM, class SVM>
class Seam_mesh
{
  typedef Seam_mesh<TM, SEM, SVM>                                 Self;

public:
  typedef TM                                                      TriangleMesh;

public:
  /// \name Underlying mesh's descriptor and iterator typedefs
  /// @{

  /// The type for the objects used to identify vertices in the underlying mesh.
  typedef typename boost::graph_traits<TM>::vertex_descriptor     TM_vertex_descriptor;

  /// The type for the objects used to identify halfedges in the underlying mesh.
  typedef typename boost::graph_traits<TM>::halfedge_descriptor   TM_halfedge_descriptor;

 /// The type for the iterators that traverse through the complete halfedge set of the underlying mesh.
 typedef typename boost::graph_traits<TM>::halfedge_iterator     TM_halfedge_iterator;

  /// The type for the objects used to identify edges in the underlying mesh.
  typedef typename boost::graph_traits<TM>::edge_descriptor       TM_edge_descriptor;
  /// @}

  /// \name Size types
  /// @{

  ///  The unsigned integer type used for representing the degree of vertices in the seam mesh.
  typedef typename boost::graph_traits<TM>::degree_size_type      degree_size_type;

  ///  The unsigned integer type used to represent the number of vertices in the seam mesh.
  typedef typename boost::graph_traits<TM>::vertices_size_type    vertices_size_type;

  ///  The unsigned integer type used for representing the number of edges in the seam mesh.
  typedef typename boost::graph_traits<TM>::edges_size_type       edges_size_type;

  ///  The unsigned integer type used for representing the number of halfedges in the seam mesh.
  typedef typename boost::graph_traits<TM>::halfedges_size_type   halfedges_size_type;

  ///  The unsigned integer type used for representing the number of faces in the seam mesh.
  typedef typename boost::graph_traits<TM>::faces_size_type       faces_size_type;
  /// @}

private:
  typedef CGAL::Unique_hash_map<TM_edge_descriptor, bool>         Seam_edge_uhm;
  typedef CGAL::Unique_hash_map<TM_vertex_descriptor, bool>       Seam_vertex_uhm;
  typedef boost::associative_property_map<Seam_edge_uhm>          Seam_edge_pmap;
  typedef boost::associative_property_map<Seam_vertex_uhm>        Seam_vertex_pmap;

private:
  const TM& tm;

  // seam edge property maps
  SEM sem;
  SVM svm;

  // combinatorics
  mutable edges_size_type number_of_seams;
  mutable vertices_size_type number_of_vertices;

public:
  /// Returns the underlying mesh.
  const TM& mesh() const
  {
    return tm;
  }

#ifdef DOXYGEN_RUNNING
  /// This class represents a halfedge of the seam mesh.
  ///
  /// Implementation note: a halfedge of the seam mesh is represented as a halfedge
  /// of the mesh and a boolean to indicate whether the halfedge is on a seam or not.
  ///
  /// \cgalModels `Descriptor`
  /// \cgalModels `LessThanComparable`
  /// \cgalModels `Hashable`
  ///
  class halfedge_descriptor
  {
public:
    TM_halfedge_descriptor tmhd;
    bool seam;

    /// %Default constructor
    halfedge_descriptor() : tmhd(), seam(false) { }

    halfedge_descriptor(TM_halfedge_descriptor tmhd, bool seam = false)
      : tmhd(tmhd), seam(seam)
    { }

#ifdef CGAL_SEAM_MESH_INSERT_OPERATOR
    /// Print the halfedge and if it is on a seam.
    friend std::ostream& operator<<(std::ostream& os, const halfedge_descriptor& hd)
    {
      os << hd.tmhd << ((hd.seam)?" on seam":"");
      return os;
    }
#endif
  };
#else
  typedef Seam_mesh_halfedge_descriptor<TM_halfedge_descriptor>   halfedge_descriptor;
#endif

#ifndef DOXYGEN_RUNNING
  class halfedge_iterator
    : public boost::iterator_facade<halfedge_iterator,
                                    halfedge_descriptor,
                                    std::forward_iterator_tag,
                                    halfedge_descriptor>
  {
    typedef boost::iterator_facade<halfedge_iterator,
                                   halfedge_descriptor,
                                   std::forward_iterator_tag,
                                   halfedge_descriptor>                  Facade;

  public:
    TM_halfedge_iterator hd, end;
    bool seam;
    const Self* mesh_;

    halfedge_iterator() : hd(), end(), seam(false), mesh_(NULL) { }

    halfedge_iterator(const Iterator_range<TM_halfedge_iterator>& ir, const Self* m)
      : hd(ir.first), end(ir.second), seam(false), mesh_(m)
    { }

    // constructor for the past the end iterator
    halfedge_iterator(const TM_halfedge_iterator& hd, const Self* m)
      : hd(hd), end(hd), seam(false), mesh_(m)
    { }

  private:
    friend class boost::iterator_core_access;

    void increment()
    {
      if(hd == end)
        return;

      if(mesh_->has_on_seam(*hd) && seam == false)
        seam = true;
      else {
        ++hd;
        seam = false;
      }
    }

    bool equal(const halfedge_iterator& other) const
    {
      return (this->mesh_ == other.mesh_) &&
             (this->hd == other.hd) && (this->seam == other.seam);
    }

    halfedge_descriptor dereference() const
    {
      return halfedge_descriptor(*hd, seam);
    }
  };
#endif

  /// This class represents a vertex of the seam mesh.
  ///
  /// Implementation note: to properly duplicate vertices that are on seams,
  /// a vertex_descriptor is in fact represented as a halfedge of the underlying
  /// mesh.
  ///
  /// \cgalModels `Descriptor`
  /// \cgalModels `LessThanComparable`
  /// \cgalModels `Hashable`
  ///
  class vertex_descriptor
  {
  public:
    halfedge_descriptor hd;

    /// %Default constructor
    vertex_descriptor() { }

    vertex_descriptor(const halfedge_descriptor& h)
      : hd(h)
    { }

    vertex_descriptor(const vertex_descriptor& other)
      : hd(other.hd)
    { }

    bool operator==(const vertex_descriptor& other) const
    {
      return (hd == other.hd);
    }

    bool operator!=(const vertex_descriptor& other) const
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

#ifdef CGAL_SEAM_MESH_INSERT_OPERATOR
    friend std::ostream& operator<<(std::ostream& os, const vertex_descriptor vd)
    {
      os << "seam mesh vertex: " << vd.hd;
      return os;
    }
#endif

    friend std::size_t hash_value(const vertex_descriptor& vd)
    {
      return hash_value(vd.hd.tmhd);
    }
  };

   // iterator
#ifndef DOXYGEN_RUNNING
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
    TM_halfedge_iterator hd, end;
    const Self* mesh_;

  public:
    /// Constructors
    vertex_iterator() : hd(), end(), mesh_(NULL) { }

    vertex_iterator(const Iterator_range<TM_halfedge_iterator>& ir, const Self* m)
      : hd(ir.first), end(ir.second), mesh_(m)
    {
      if(!is_current_vertex_valid())
        increment();
    }

    // constructor for the past the end iterator
    vertex_iterator(const TM_halfedge_iterator& hd, const Self* m)
      : hd(hd), end(hd), mesh_(m)
    { }

  private:
    friend class boost::iterator_core_access;

    bool is_current_vertex_valid() const
    {
      if(hd == end)
        return true;

      TM_vertex_descriptor tvd = CGAL::target(*hd, mesh_->mesh());
      if((!mesh_->has_on_seam(tvd))&& (CGAL::halfedge(tvd, mesh_->mesh()) == *hd)) {
        return true;
      }

      if(mesh_->has_on_seam(CGAL::edge(*hd, mesh_->mesh()))) {
        return true;
      }

      if(mesh_->has_on_seam(tvd) &&
         is_border(CGAL::opposite(*hd, mesh_->mesh()), mesh_->mesh())) {
        return true;
      }

      return false;
    }

    void increment()
    {
      if(hd == end)
        return;

      do {
        ++hd;
        if(is_current_vertex_valid())
          return;
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
  };
#endif

  /// This class represents an edge of the seam mesh.
  ///
  /// \cgalModels `Descriptor`
  ///
  class edge_descriptor
  {
  public:
    halfedge_descriptor hd;
    const Self* mesh_;

#ifdef CGAL_SEAM_MESH_INSERT_OPERATOR
    friend
    std::ostream& operator<<(std::ostream& os, const edge_descriptor& ed)
    {
      os << ed.hd;
      return os;
    }
#endif

    edge_descriptor()
      : mesh_(NULL)
    {}

    edge_descriptor(const halfedge_descriptor& hd, const Self* m)
      : hd(hd), mesh_(m)
    {}

    friend bool operator==(edge_descriptor e1, edge_descriptor e2)
    {
      return (e1.hd == e2.hd) || (e1.hd == e2.mesh_->opposite(e2.hd)); 
    }

    friend bool operator!=(edge_descriptor e1, edge_descriptor e2)
    {
      return ! (e1 == e2);
    }
  };

#ifndef DOXYGEN_RUNNING
   // iterator
  class edge_iterator
    : public boost::iterator_facade<edge_iterator,
                                    edge_descriptor,
                                    std::forward_iterator_tag,
                                    edge_descriptor>
  {
  public:
    typedef boost::iterator_facade<edge_iterator,
                                   edge_descriptor,
                                   std::forward_iterator_tag,
                                   edge_descriptor>                  Facade;

    TM_halfedge_iterator hd, end;
    bool seam;
    const Self* mesh_;

  public:
    edge_iterator() : hd(), end(), seam(false), mesh_(NULL) { }

    edge_iterator(const Iterator_range<TM_halfedge_iterator>& ir, const Self* m)
      : hd(ir.first), end(ir.second), seam(false), mesh_(m)
    {
      if(!is_current_edge_valid())
        increment();
    }

    // constructor for the past the end iterator
    edge_iterator(const TM_halfedge_iterator& hd, const Self* m)
      : hd(hd), end(hd), seam(false), mesh_(m)
    { }

  private:
    friend class boost::iterator_core_access;

    bool is_current_edge_valid() const
    {
      if(hd == end)
        return true;

      return CGAL::source(*hd, mesh_->mesh()) < CGAL::target(*hd, mesh_->mesh());
    }

    void increment()
    {
      if(hd == end)
        return;

      do {
        // a seam gives two edges
        if(mesh_->has_on_seam(*hd) && seam == false)
          seam = true;
        else {
          ++hd;
          seam = false;
        }

        if(is_current_edge_valid())
          return;

      } while(true);
    }

    bool equal(const edge_iterator& other) const
    {
      return (this->mesh_ == other.mesh_) &&
             (this->hd == other.hd) && (this->seam == other.seam);
    }

    edge_descriptor dereference() const
    {
      // if 'seam' is true, output an interior halfedge rather than the halfedge
      // on the opposite virtual border

      if(seam)
        return edge_descriptor(mesh_->opposite(halfedge_descriptor(*hd, true)), mesh_);
      else
        return edge_descriptor(halfedge_descriptor(*hd, false), mesh_);
    }
  };
#endif

  /// This class represents a face of the seam mesh.
  ///
  /// \cgalModels `Descriptor`
  ///
  typedef typename boost::graph_traits<TM>::face_descriptor     face_descriptor;

#ifndef DOXYGEN_RUNNING
  typedef typename boost::graph_traits<TM>::face_iterator       face_iterator;
#endif

public:
  /// \name Seam query functions
  /// @{

  /// Returns `true` if the vertex is on the seam.
  bool has_on_seam(TM_vertex_descriptor vd) const
  {
    return get(svm, vd);
  }

  /// Returns `true` if the edge is on the seam.
  bool has_on_seam(TM_edge_descriptor ed) const
  {
    return get(sem, ed);
  }

  /// Returns `true` if the halfedge is on the seam.
  bool has_on_seam(TM_halfedge_descriptor tmhd) const
  {
    return get(sem, CGAL::edge(tmhd, tm));
  }

  /// Returns `true` if the halfedge is on the seam.
  bool has_on_seam(const halfedge_descriptor& hd) const
  {
    return has_on_seam(CGAL::edge(hd, tm));
  }

  /// Return the number of seam edges in the seam mesh.
  edges_size_type number_of_seam_edges() const
  {
    return number_of_seams;
  }
  /// @}

  // Set the number of seam edges.
  void set_seam_edges_number(const edges_size_type sn) const
  {
    number_of_seams = sn;
  }

public:
  /// \name Range Types
  ///
  ///@{

  /// @cond CGAL_BEGIN_END
  /// Start iterator for vertices.
  vertex_iterator vertices_begin() const
  {
    Iterator_range<TM_halfedge_iterator> ir = CGAL::halfedges(tm);
    return vertex_iterator(ir, this);
  }

  /// End iterator for vertices.
  vertex_iterator vertices_end() const
  {
    Iterator_range<TM_halfedge_iterator> ir = CGAL::halfedges(tm);
    return vertex_iterator(ir.second, this);
  }
  /// @endcond

  /// Returns the iterator range of the vertices of the mesh.
  Iterator_range<vertex_iterator> vertices() const
  {
    Iterator_range<TM_halfedge_iterator> ir = CGAL::halfedges(tm);
    vertex_iterator beg(ir, this);
    vertex_iterator end(ir.second, this);
    return make_range(beg, end);
  }

  /// @cond CGAL_BEGIN_END
  /// Start iterator for halfedges.
  halfedge_iterator halfedges_begin() const
  {
    Iterator_range<TM_halfedge_iterator> ir = CGAL::halfedges(tm);
    return halfedge_iterator(ir, this);
  }


  /// End iterator for halfedges.
  halfedge_iterator halfedges_end() const
  {
    Iterator_range<TM_halfedge_iterator> ir = CGAL::halfedges(tm);
    return halfedge_iterator(ir.second, this);
  }
  /// @endcond

  /// Returns the iterator range of the halfedges of the mesh.
  Iterator_range<halfedge_iterator> halfedges() const
  {
    Iterator_range<TM_halfedge_iterator> ir = CGAL::halfedges(tm);
    halfedge_iterator beg(ir, this);
    halfedge_iterator end(ir.second, this);
    return make_range(beg, end);
  }

  /// @cond CGAL_BEGIN_END
  /// Start iterator for edges.
  edge_iterator edges_begin() const
  {
    Iterator_range<TM_halfedge_iterator> ir = CGAL::halfedges(tm);
    return edge_iterator(ir, this);
  }

  /// End iterator for edges.
  edge_iterator edges_end() const
  {
    Iterator_range<TM_halfedge_iterator> ir = CGAL::halfedges(tm);
    return edge_iterator(ir.second, this);
  }
  /// @endcond

  /// Returns the iterator range of the edges of the mesh.
  Iterator_range<edge_iterator> edges() const
  {
    Iterator_range<TM_halfedge_iterator> ir = CGAL::halfedges(tm);
    edge_iterator beg(ir, this);
    edge_iterator end(ir.second, this);
    return make_range(beg, end);
  }

  /// @cond CGAL_BEGIN_END
  /// Start iterator for faces.
  face_iterator faces_begin() const
  {
    return CGAL::faces(tm).begin();
  }

  /// End iterator for faces.
  face_iterator faces_end() const
  {
    return CGAL::faces(tm).end();
  }
  /// @endcond

  /// Returns the iterator range of the faces of the mesh.
  Iterator_range<face_iterator> faces() const
  {
    return CGAL::faces(tm);
  }
  ///@}

public:
  /// \name Memory Management
  /// @{

  /// Returns the number of vertices in the seam mesh.
  vertices_size_type num_vertices() const
  {
    if(number_of_vertices == static_cast<vertices_size_type>(-1)) {
      number_of_vertices = static_cast<vertices_size_type>(vertices().size());
    }

    return number_of_vertices;
  }

  /// Returns the number of halfedges in the seam mesh.
  halfedges_size_type num_halfedges() const
  {
    return CGAL::num_halfedges(tm) + 2 * number_of_seams;
  }

  /// Returns the number of edges in the seam mesh.
  halfedges_size_type num_edges() const
  {
    return CGAL::num_edges(tm) + number_of_seams;
  }

  /// Returns the number of faces in the seam mesh.
  faces_size_type num_faces() const
  {
    return CGAL::num_faces(tm);
  }
  /// @}

public:
  /// \name Degree Functions
  /// @{

  /// Returns the number of incident halfedges of vertex `v`.
  degree_size_type degree(vertex_descriptor v) const
  {
    degree_size_type count(0);

    if(v.hd == halfedge_descriptor())
      return count;

    Halfedge_around_target_circulator<Self> hatc(v.hd, *this), end = hatc;
    CGAL_For_all(hatc, end) {
      ++count;
    }

    return count;
  }
  /// @}

public:
  /// \name Low-Level Connectivity Functions
  ///

#ifndef DOXYGEN_RUNNING
  ///@{

  /// Returns the edge that contains halfedge `h` as one of its two halfedges.
  edge_descriptor edge(halfedge_descriptor h) const
  {
    return edge_descriptor(h,this);
  }

  /// Returns the halfedge corresponding to the edge `e`.
  halfedge_descriptor halfedge(edge_descriptor e) const
  {
    return e.hd;
  }
  /// @}

  ///@{

  /// Returns an incoming halfedge of vertex `v`.
  /// If `v` is a seam vertex, this will be the halfedge whose target is `v` and
  /// whose opposite is a virtual border halfedge.
  /// Otherwise, the rules of the underlying mesh are followed.
  /// \invariant `target(halfedge(v)) == v`
  halfedge_descriptor halfedge(vertex_descriptor v) const
  {
    TM_halfedge_descriptor h(v);
    return halfedge_descriptor(h, false /*not on seam*/);
  }

  /// Finds a halfedge between two vertices. Returns a default constructed
  /// `halfedge_descriptor`, if  `source` and  `target` are not connected.
  std::pair<halfedge_descriptor, bool> halfedge(vertex_descriptor u,
                                                vertex_descriptor v) const
  {
    halfedge_descriptor hdv(v);
    Halfedge_around_target_circulator<Self> hatcv(hdv, *this), endv = hatcv;
    CGAL_For_all(hatcv, endv) {
      halfedge_descriptor hd_around_v = *hatcv;
      TM_halfedge_descriptor tmhd_around_v = hd_around_v.tmhd;
      if(CGAL::source(tmhd_around_v, tm) == CGAL::target(u, tm)) {
        // found a u next to v in the base mesh 'tm',
        // but we must check that is also the case in the seam mesh 'this'
        // that means that the halfedge 'u' is incident to the source of hd_around_v

        halfedge_descriptor opp_hd_around_v = opposite(hd_around_v);
        Halfedge_around_target_circulator<Self> hatcu(opp_hd_around_v, *this),
                                                end2 = hatcu;
        CGAL_For_all(hatcu, end2) {
          halfedge_descriptor hd_around_u = *hatcu;
          if(hd_around_u == u.hd) {
            return std::make_pair(hd_around_v, true/*valid*/);
          }
        }
      }
    }

    // halfedge doesn't exist
    return std::make_pair(halfedge_descriptor(), false/*invalid*/);
  }

  /// Finds an edge between two vertices. Returns a default constructed
  /// `edge`, if  `source` and  `target` are not connected.
  std::pair<edge_descriptor, bool> edge(vertex_descriptor u, vertex_descriptor v) const
  {
    std::pair<halfedge_descriptor, bool> he = halfedge(u, v);
    return std::make_pair(edge_descriptor(he.first, this), he.second);
  }

  /// Returns a halfedge of face `f`.
  halfedge_descriptor halfedge(face_descriptor f) const
  {
    TM_halfedge_descriptor hd = CGAL::halfedge(f, tm);
    return halfedge_descriptor(hd, false/*not on seam*/);
  }

  /// Returns the face incident to halfedge `h`.
  face_descriptor face(halfedge_descriptor h) const
  {
    if(h.seam)
      return boost::graph_traits<CGAL::Seam_mesh<TM, SEM, SVM> >::null_face();

    return CGAL::face(h, tm);
  }

public:
  /// Returns the next halfedge within the incident face.
  halfedge_descriptor next(const halfedge_descriptor& hd) const
  {
    if((!hd.seam) && (!is_border(hd.tmhd, tm)))
      return halfedge_descriptor(CGAL::next(hd.tmhd, tm));

    Halfedge_around_target_circulator<TM> hatc(hd.tmhd, tm);
    do {
      --hatc;
    } while((!has_on_seam(*hatc)) && (!is_border(CGAL::opposite(*hatc, tm), tm)));

    return halfedge_descriptor(CGAL::opposite(*hatc, tm),
                               !is_border(CGAL::opposite(*hatc, tm), tm));
  }

  /// Returns the previous halfedge within the incident face.
  halfedge_descriptor prev(const halfedge_descriptor& hd) const
  {
    if((!hd.seam) && (!is_border(hd.tmhd, tm)))
      return halfedge_descriptor(CGAL::prev(hd.tmhd, tm));

    Halfedge_around_source_circulator<TM> hatc(hd.tmhd, tm);
    do {
      ++hatc;
    } while((!has_on_seam(*hatc)) && (!is_border(CGAL::opposite(*hatc, tm), tm)));

    return halfedge_descriptor(CGAL::opposite(*hatc, tm),
                               !is_border(CGAL::opposite(*hatc, tm), tm));
  }

  /// Returns the opposite halfedge of `hd`.
  halfedge_descriptor opposite(const halfedge_descriptor& hd) const
  {
    if(!hd.seam)
      return halfedge_descriptor(CGAL::opposite(hd.tmhd, tm), has_on_seam(hd));

    return halfedge_descriptor(CGAL::opposite(hd.tmhd, tm), false /*not on seam*/);
  }

  /// Returns the vertex the halfedge `h` emanates from.
  vertex_descriptor target(halfedge_descriptor hd) const
  {
    TM_halfedge_descriptor tmhd(hd);

    if(!has_on_seam(CGAL::target(tmhd, tm))) {
      tmhd = CGAL::halfedge(CGAL::target(tmhd, tm), tm);
      return vertex_descriptor(halfedge_descriptor(tmhd, false /*not on seam*/));
    }

    if(hd.seam)
      return target(halfedge_descriptor(CGAL::prev(CGAL::opposite(tmhd, tm), tm)));

    while((!has_on_seam(tmhd)) && (!is_border(CGAL::opposite(tmhd, tm), tm))) {
      tmhd = CGAL::prev(CGAL::opposite(tmhd, tm), tm);
    }

    return vertex_descriptor(halfedge_descriptor(tmhd));
  }

  /// Returns the vertex the halfedge `h` emanates from.
  vertex_descriptor source(const halfedge_descriptor& hd) const
  {
    return target(opposite(hd));
  }

  vertex_descriptor source(edge_descriptor e) const
  {
    return source(e.hd);
  }

  vertex_descriptor target(edge_descriptor e) const
  {
    return target(e.hd);
  }

  /// @}

  void build_TM_vertices_vector(std::vector<TM_vertex_descriptor>& tm_vds) const
  {
    assert(tm_vds.empty());

    // If the input is a list of integers, we need to build a correspondence
    // between vertices and integers.
    tm_vds.reserve(CGAL::num_vertices(tm));
    typedef typename boost::graph_traits<TM>::vertex_iterator  TM_vertex_iterator;
    TM_vertex_iterator tmvi = CGAL::vertices(tm).begin(),
                       tmvi_end = CGAL::vertices(tm).end();
    CGAL_For_all(tmvi, tmvi_end) {
      tm_vds.push_back(*tmvi);
    }
  }
#endif // ndef DOXYGEN_RUNNING

public:
  /// \name Seam selection
  /// @{

  /// Mark the edge of the underlying mesh that has extremities the vertices
  /// `tm_vd_s` and `tm_vd_s` as a seam edge.
  ///
  /// \return whether the edge was successfully marked or not.
  ///         Marking will fail if:
  ///         - No edge of the underlying mesh exist with extremities `tm_vd_s` and `tm_vd_s`,
  ///         - the edge of the underlying mesh with extremities `tm_vd_s` and `tm_vd_s` is a border edge, or
  ///         - the edge of the underlying mesh with extremities `tm_vd_s` and `tm_vd_s` is already a seam edge.
  ///
  bool add_seam(TM_vertex_descriptor tm_vd_s, TM_vertex_descriptor tm_vd_t)
  {
    std::pair<TM_edge_descriptor, bool> tmed = CGAL::edge(tm_vd_s, tm_vd_t, tm);
    if(!tmed.second) {
#ifdef SEAM_MESH_DEBUG
      std::cerr << "Warning: Ignored a constraint because it is not a valid edge of the mesh" << std::endl;
#endif
      return false;
    }

    if(!is_border(tmed.first, tm)) { // ignore seams that are also a border edge
      if(get(sem, tmed.first) == true) {
#ifdef SEAM_MESH_DEBUG
        std::cerr << "Warning: Ignored a constraint because it is already marked as a seam" << std::endl;
#endif
        return false;
      }

      put(sem, tmed.first, true);
      put(svm, tm_vd_s, true);
      put(svm, tm_vd_t, true);
      ++number_of_seams;
    } else {
#ifdef SEAM_MESH_DEBUG
      std::cerr << "Warning: Ignored a constraint because it is on the border of the mesh" << std::endl;
#endif
      return false;
    }

    return true;
  }

  /// Create new seams.
  ///
  /// The edges to be marked as seams are described by the range [first, last) of
  /// vertices of the underlying mesh. Each edge to be marked is described
  /// by two consecutive iterators.
  ///
  /// \returns one of the halfedges of the seam mesh that is on a seam.
  ///
  /// \pre InputIterator must be a model of `InputIterator`.
  /// \pre The value type of `InputIterator` must be `boost::graph_traits<TM>::%vertex_descriptor`.
  /// \pre There is an even number of vertices.
  template<class InputIterator>
  TM_halfedge_descriptor add_seams(InputIterator first, InputIterator last)
  {
    // must have an even number of input vertices
    assert(std::distance(first, last) % 2 == 0);

    TM_halfedge_descriptor tmhd = boost::graph_traits<TM>::null_halfedge();
    InputIterator it = first;

    while(it!=last) {
      TM_vertex_descriptor v1 = *(it++);
      TM_vertex_descriptor v2 = *(it++);

      if(!add_seam(v1, v2))
        continue;

      if(tmhd == boost::graph_traits<TM>::null_halfedge()) {
        tmhd = CGAL::halfedge(CGAL::edge(v1, v2, tm).first, tm);
      }
    }

    return tmhd;
  }

  /// Create new seams.
  ///
  /// A seam edge is described by a pair of integers. The integer index
  /// of a vertex of the underlying mesh is given by its position
  /// in the container `tm_vds`.
  ///
  /// \tparam VdContainer must be a model of <a href="http://en.cppreference.com/w/cpp/concept/SequenceContainer"><tt>SequenceContainer</tt></a> (that is, provide
  ///         the functions: `operator[]` and `at()`).
  ///
  /// \returns one of the halfedges of the seam mesh that is on a seam.
  ///
  /// \pre The stream must contain an even number of values.
  template<typename VdContainer>
  TM_halfedge_descriptor add_seams(std::ifstream& in,
                                   const VdContainer& tm_vds)
  {
    std::vector<TM_vertex_descriptor> seam_vertices;
    std::size_t s, t;
    while(in >> s >> t) {
      seam_vertices.push_back(tm_vds.at(s));
      seam_vertices.push_back(tm_vds.at(t));
    }

    return add_seams(seam_vertices.begin(), seam_vertices.end());
  }

  /// Create new seams.
  ///
  /// A seam edge is described by a pair of integers. The integer
  /// index of a vertex of the underlying mesh is defined as its position when
  /// iterating over the vertices of the underlying mesh with
  /// `boost::graph_traits<TM>::%vertices()`.
  ///
  /// \returns one of the halfedges of the seam mesh that is on a seam.
  ///
  /// \pre The stream must contain an even number of values.
  TM_halfedge_descriptor add_seams(std::ifstream& in)
  {
    std::vector<TM_vertex_descriptor> tm_vds;
    build_TM_vertices_vector(tm_vds);
    return add_seams(in, tm_vds);
  }

  /// Create new seams.
  ///
  /// A seam edge is described by a pair of integers. The integer index
  /// of a vertex of the underlying mesh is given by its position
  /// in the container `tm_vds`.
  ///
  /// \returns one of the halfedges of the seam mesh that is on a seam.
  ///
  /// \tparam VdContainer must be a model of <a href="http://en.cppreference.com/w/cpp/concept/SequenceContainer"><tt>SequenceContainer</tt></a> (that is, provide
  ///         the functions: `operator[]` and `at()`).
  ///
  /// \pre filename should be the name of a CGAL selection file: edges are
  ///      described by pairs of integers, on the third line of the file.
  template<typename VdContainer>
  TM_halfedge_descriptor add_seams(const char* filename,
                                   const VdContainer& tm_vds)
  {
    TM_halfedge_descriptor tmhd = boost::graph_traits<TM>::null_halfedge();

    // Check the file type
    std::string str = filename;
    if(str.substr(str.length() - 14) != ".selection.txt") {
      std::cerr << "Error: seams must be given by a *.selection.txt file" << std::endl;
      return tmhd;
    }

    // A '.selection.txt' file has:
    // -- the first line for selected vertices,
    // -- the second line for selected faces,
    // -- the third line for selected edges
    std::ifstream in(filename);
    std::string line;

    // skip two lines to get the istream to be at the beginning of the third line
    if(!std::getline(in, line) || !std::getline(in, line)) {
#ifdef SEAM_MESH_DEBUG
      std::cerr << "Error: no seams in input file: " << filename << std::endl;
#endif
      return tmhd;
    }

    return add_seams(in, tm_vds);
  }

  /// Create new seams.
  ///
  /// A seam edge is described by a pair of integers. The integer
  /// index of a vertex of the underlying mesh is defined as its position when
  /// iterating over the vertices of the underlying mesh with
  /// `boost::graph_traits<TM>::%vertices()`.
  ///
  /// \returns one of the halfedges of the seam mesh that is on a seam.
  ///
  /// \pre filename should be the name of a CGAL selection file: edges are
  ///      described by pairs of integers, on the third line of the file.
  TM_halfedge_descriptor add_seams(const char* filename)
  {
    std::vector<TM_vertex_descriptor> tm_vds;
    build_TM_vertices_vector(tm_vds);
    return add_seams(filename, tm_vds);  }

  /// @}

  /// Constructs a seam mesh for a triangle mesh and an edge and vertex property map
  ///
  /// \param tm the underlying mesh
  /// \param sem the edge property map with value `true` for seam edges
  /// \param svm the vertex property map with value `true` for seam vertices
  ///
  /// @note the vertices must be exactly the vertices on the seam edges.
  Seam_mesh(const TM& tm, const SEM& sem, const SVM svm)
    : tm(tm),
      sem(sem), svm(svm),
      number_of_seams(0), number_of_vertices(-1)
  { }
};

} // namespace CGAL

#endif //CGAL_SEAM_MESH_H
