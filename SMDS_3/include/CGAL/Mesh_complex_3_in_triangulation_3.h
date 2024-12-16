// Copyright (c) 2009-2014 INRIA Sophia-Antipolis (France).
// Copyright (c) 2010-2013 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stephane Tayeb, Clement Jamin
//
//******************************************************************************
// File Description :
//******************************************************************************

#ifndef CGAL_MESH_COMPLEX_3_IN_TRIANGULATION_3_H
#define CGAL_MESH_COMPLEX_3_IN_TRIANGULATION_3_H

#include <CGAL/license/SMDS_3.h>

#include <CGAL/SMDS_3/Mesh_complex_3_in_triangulation_3_fwd.h>
#include <CGAL/disable_warnings.h>
#include <CGAL/iterator.h>
#include <CGAL/SMDS_3/utilities.h>
#include <CGAL/SMDS_3/internal/Boundary_of_subdomain_of_complex_3_in_triangulation_3_to_off.h>
#include <CGAL/Time_stamper.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Union_find.h>
#include <CGAL/SMDS_3/io_signature.h>

#include <CGAL/IO/File_medit.h>
#include <CGAL/IO/File_maya.h>

#include <boost/bimap/bimap.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <CGAL/boost/iterator/transform_iterator.hpp>
#include <boost/iterator/iterator_adaptor.hpp>
#include <boost/unordered_map.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/functional/hash.hpp>

#include <iostream>
#include <fstream>


#ifdef CGAL_LINKED_WITH_TBB
#include <atomic>
#include <tbb/concurrent_hash_map.h>

namespace CGAL {
  template < class DSC, bool Const >
  std::size_t tbb_hasher(const CGAL::internal::CC_iterator<DSC, Const>& it)
  {
    return CGAL::internal::hash_value(it);
  }

  // As Marc Glisse pointed out the TBB hash of a std::pair is
  // simplistic and leads to the
  // TBB Warning: Performance is not optimal because the hash function
  //              produces bad randomness in lower bits in class
  //              tbb::interface5::concurrent_hash_map
  template < class DSC, bool Const >
  std::size_t tbb_hasher(const std::pair<CGAL::internal::CC_iterator<DSC, Const>,
                                         CGAL::internal::CC_iterator<DSC, Const> >& p)
  {
    return boost::hash<std::pair<CGAL::internal::CC_iterator<DSC, Const>,
                                 CGAL::internal::CC_iterator<DSC, Const> > >()(p);
  }

  struct Hash_compare_for_TBB {
    template < class DSC, bool Const >
    std::size_t hash(const std::pair<CGAL::internal::CC_iterator<DSC, Const>,
                                     CGAL::internal::CC_iterator<DSC, Const> >& p) const
    {
      return tbb_hasher(p);
    }
    template < class DSC, bool Const >
    std::size_t operator()(const CGAL::internal::CC_iterator<DSC, Const>& it)
    {
      return CGAL::internal::hash_value(it);
    }
    template <typename T>
    bool equal(const T& v1, const T& v2) const {
      return v1 == v2;
    }
  };
}//end namespace CGAL
#endif

namespace CGAL {

  namespace SMDS_3 {

    namespace details {

      template <typename Tr>
      class C3t3_helper_class
      {
      protected:
        typedef typename Tr::Vertex_handle Vertex_handle;
        typedef typename Tr::Cell_handle   Cell_handle;
        typedef typename Tr::Facet         Facet;
        typedef typename Tr::Edge          Edge;

        typedef std::pair<Vertex_handle, Vertex_handle> Pair_of_vertices;

        // computes and returns an ordered pair of Vertex
        Pair_of_vertices
        make_ordered_pair(const Vertex_handle vh1, const Vertex_handle vh2) const {
          if (vh1 < vh2) {
            return std::make_pair(vh1, vh2);
          }
          else {
            return std::make_pair(vh2, vh1);
          }
        }
        // same from an Edge
        Pair_of_vertices
        make_ordered_pair(const Edge e) const {
          return make_ordered_pair(e.first->vertex(e.second),
            e.first->vertex(e.third));
        }
        Facet canonical_facet(Cell_handle c, int i) const {
          Cell_handle c2 = c->neighbor(i);
          return (c2 < c) ? std::make_pair(c2, c2->index(c)) : std::make_pair(c, i);
        }
      }; // end class template C3t3_helper_class

    } // end namespace SMDS_3::details
  } //end namespace SMDS_3

/*!
  \ingroup PkgSMDS3Classes

  \brief A data structure to represent and maintain a 3D complex embedded
  in a 3D triangulation.

  The class `Mesh_complex_3_in_triangulation_3` implements a data structure
  to store the 3D restricted Delaunay triangulation used by a mesh
  generation process.

  This class is a model of the concept
  `MeshComplexWithFeatures_3InTriangulation_3`.

  \tparam Tr can be instantiated with any 3D
  triangulation of \cgal provided that its
  vertex and cell base class are models of the concepts
  `SimplicialMeshVertexBase_3` and `SimplicialMeshCellBase_3`, respectively.

  \tparam CornerIndex Type of indices for corners (i.e.\f$ 0\f$--dimensional features)
  of the discretized geometric domain.
  It must be a model of `CopyConstructible`, `Assignable`, `DefaultConstructible` and
  `LessThanComparable`.
  It must match the `Corner_index` of the model
  of the `MeshDomainWithFeatures_3` concept when used for mesh generation.

  \tparam CurveIndex Type of indices for curves (i.e. \f$ 1\f$-dimensional features)
  of the discretized geometric domain.
  It must be a model of `CopyConstructible`, `Assignable`, `DefaultConstructible` and
  `LessThanComparable`. The default constructed value must be the value for an edge which
  does not approximate a 1-dimensional feature of the geometric domain.
  It must match the `Curve_index` types of the model
  of the `MeshDomainWithFeatures_3` concept when used for mesh generation.

  Those two last template parameters default to `int`, so that they can be ignored
  if the domain used for mesh generation does not include 0 and 1-dimensionnal features (i.e
  is only a model of the concept `MeshDomain_3`).

  \cgalModels{MeshComplexWithFeatures_3InTriangulation_3}

  \sa \link make_mesh_3() `CGAL::make_mesh_3()`\endlink
  \sa \link refine_mesh_3() `CGAL::refine_mesh_3()`\endlink
  \sa `MeshComplex_3InTriangulation_3`
  \sa `MeshComplexWithFeatures_3InTriangulation_3`
  \sa `SimplicialMeshCellBase_3`,
  \sa `SimplicialMeshVertexBase_3`

*/
template <typename Tr,
          typename CornerIndex,
          typename CurveIndex>
class Mesh_complex_3_in_triangulation_3
#ifndef DOXYGEN_RUNNING
  : public CGAL::SMDS_3::details::C3t3_helper_class<Tr>
  , public CGAL::SMDS_3::internal::Debug_messages_tools
#endif
{
public:
  typedef typename Tr::Concurrency_tag                   Concurrency_tag;

private:
  typedef Mesh_complex_3_in_triangulation_3<
    Tr,CornerIndex,CurveIndex>                            Self;
  typedef SMDS_3::details::C3t3_helper_class<Tr>           Base;
  typedef CGAL::Hash_handles_with_or_without_timestamps   Hash_fct;

public:

#ifndef CGAL_NO_DEPRECATED_CODE
  typedef CurveIndex Curve_segment_index;
#endif

/// \name Types
/// @{
  typedef Tr                                            Triangulation;
  typedef typename Tr::size_type                        size_type;
  typedef typename Tr::Point                            Point;
  typedef typename Tr::Edge                             Edge;
  typedef typename Tr::Facet                            Facet;
  typedef typename Tr::Vertex_handle                    Vertex_handle;
  typedef typename Tr::Cell_handle                      Cell_handle;
  /*!
  Index type.
  */
  typedef typename Tr::Vertex::Index Index;
  /*!
  Surface index type.
  */
  typedef typename Tr::Cell::Surface_patch_index Surface_patch_index;
  /*!
  Subdomain index type.
  */
  typedef typename Tr::Cell::Subdomain_index Subdomain_index;
  /*!
  Corner index type.
  */
  typedef CornerIndex Corner_index;
  /*!
  Curve index type.
  */
  typedef CurveIndex Curve_index;
/// @}


private:
  // Type to store the edges:
  //  - a set of std::pair<Vertex_handle,Vertex_handle> (ordered at insertion)
  //  - which allows fast lookup from one Vertex_handle
  //  - each element of the set has an associated info (Curve_index) value
  typedef boost::bimaps::bimap<
    boost::bimaps::multiset_of<Vertex_handle>,
    boost::bimaps::multiset_of<Vertex_handle>,
    boost::bimaps::set_of_relation<>,
    boost::bimaps::with_info<Curve_index> >           Edge_map;

  typedef typename Edge_map::value_type               Internal_edge;

  // Type to store the corners
  typedef boost::unordered_map<Vertex_handle,
                               Corner_index,
                               Hash_fct>              Corner_map;

  // Type to store far vertices
  typedef std::vector<Vertex_handle>                  Far_vertices_vec;

public:
  enum Face_status {
    NOT_IN_COMPLEX = 0,
    ISOLATED = 1, // - An ISOLATED edge is a marked edge,
                  //   without any incident facets.
    BOUNDARY,     // - An edge is on BOUNDARY if it has only
                  //   one incident facet.
                  // - A vertex is on BOUNDARY if all its
                  //   incident edges are REGULAR or on
                  //   BOUNDARY, at least one is on
                  //   BOUNDARY, and the incident facets
                  //   form only one connected component.
    REGULAR,      // - A facet that is in the complex is
                  //   REGULAR.
                  // - An edge is REGULAR if it has
                  //   exactly two incident facets.
                  // - A vertex is REGULAR if all it
                  //   incident edges are REGULAR, and the
                  //   incident facets form only one
                  //   connected component.
    SINGULAR      // - SINGULAR is for all other cases.
  };

/// \name Creation
/// @{
   /**
   * @brief Constructor
   * builds an empty 3D complex.
   */
  Mesh_complex_3_in_triangulation_3();
  /**
   * Copy constructor
   */
  Mesh_complex_3_in_triangulation_3(const Self& rhs);

  /**
   * Move constructor
   */
  Mesh_complex_3_in_triangulation_3(Self&& rhs);

  /**
  * Assignment operator, also serves as move-assignment
  */
  Self& operator=(Self rhs)
  {
    swap(rhs);
    return *this;
  }

  /**
   * swaps `this` and `rhs`
   */
  void swap(Self& rhs)
  {
    Swap_elements<Concurrency_tag> swapper;
    swapper(rhs.number_of_facets_, number_of_facets_);
    tr_.swap(rhs.tr_);
    swapper(rhs.number_of_cells_, number_of_cells_);

    edges_.swap(rhs.edges_);
    corners_.swap(rhs.corners_);
    far_vertices_.swap(rhs.far_vertices_);
  }

/// @}

/// \name Access Functions
/// @{
  /// returns a const reference to the triangulation
  const Triangulation& triangulation() const { return tr_; }
/// @}

/// \name Non const access
/// @{
    /// returns a reference to the triangulation
  Triangulation& triangulation() { return tr_; }
/// @}


/// \name Modifiers
/// @{
  /**
   * clears data of the complex
   */
  void clear()
  {
    number_of_cells_ = 0;
    number_of_facets_ = 0;
    clear_manifold_info();
    tr_.clear();
    edges_.clear();
    corners_.clear();
    far_vertices_.clear();
  }

  /** adds cell \p cell to the 3D complex, with subdomain index \p index
  */
  void add_to_complex(const Cell_handle& cell, const Subdomain_index& index)
  {
    CGAL_precondition(!(index == Subdomain_index()));

    if (!is_in_complex(cell))
    {
      set_subdomain_index(cell, index);
      ++number_of_cells_;
    }
  }
  /** adds facet \p facet to the 2D complex, with surface index \p index
  */
  void add_to_complex(const Facet& facet, const Surface_patch_index& index)
  {
    add_to_complex(facet.first, facet.second, index);
  }
  /** adds facet(\p cell, \p i) to the 2D complex, with surface index \p index
  */
  void add_to_complex(const Cell_handle& cell,
                      const int i,
                      const Surface_patch_index& index);

  /** adds edge \p e to complex, with curve index \p index
   */
  void add_to_complex(const Edge& e,
                      const Curve_index& index)
  {
    add_to_complex(e.first->vertex(e.second),
                   e.first->vertex(e.third),
                   index);
  }

  /**
   * adds edge (\p v1, \p v2) to the 1D complex, with `Curve_index` \p index
   */
  void add_to_complex(const Vertex_handle& v1,
                      const Vertex_handle& v2,
                      const Curve_index& index)
  {
    add_to_complex(make_internal_edge(v1,v2), index);
  }

  /**
   * marks vertex \p v as a corner of the complex
   */
  void add_to_complex(const Vertex_handle& v, const Corner_index& index)
  {
    v->set_dimension(0);
    corners_.insert(std::make_pair(v,index));
  }

  /** removes cell \p cell from the 3D complex
  */
  void remove_from_complex(const Cell_handle& cell)
  {
    if (is_in_complex(cell))
    {
      set_subdomain_index(cell, Subdomain_index());
      --number_of_cells_;
    }
  }
  /** removes facet \p facet from the 2D complex
  */
  void remove_from_complex(const Facet& facet);

  /** removes facet(\p cell, \p i) from the 2D complex
  */
  void remove_from_complex(const Cell_handle& c, const int i) {
    remove_from_complex(Facet(c, i));
  }
  /**
   * removes edge \p e from the 1D complex
   */
  void remove_from_complex(const Edge& e)
  {
    remove_from_complex(e.first->vertex(e.second), e.first->vertex(e.third));
  }

  /**
   * removes edge (v1,v2) from the 1D complex
   */
  void remove_from_complex(const Vertex_handle& v1, const Vertex_handle& v2)
  {
    remove_from_complex(make_internal_edge(v1,v2));
  }

  /**
   * removes vertex \p v from the complex
   */
  void remove_from_complex(const Vertex_handle& v)
  {
    v->set_dimension(-1);
    corners_.erase(v);
  }

  /** sets the index of vertex \p vertex to \p index
  */
  void set_index(const Vertex_handle& vertex, const Index& index) const
  {
    vertex->set_index(index);
  }
  /** sets the surface index of facet \p facet to \p index
  */
  void set_surface_patch_index(const Facet& f, const Surface_patch_index& index)
  {
    set_surface_patch_index(f.first, f.second, index);
  }
  /** sets the surface index of facet(\p cell, \p i) to \p index
  */
  void set_surface_patch_index(const Cell_handle& cell,
    const int i,
    const Surface_patch_index& index) const
  {
    cell->set_surface_patch_index(i, index);
  }
  /** sets the subdomain index of cell \p cell to \p index
  */
  void set_subdomain_index(const Cell_handle& cell,
    const Subdomain_index& index) const
  {
    cell->set_subdomain_index(index);
  }
  /** sets the dimension of vertex \p vertex to \p dimension
  */
  void set_dimension(const Vertex_handle& vertex, int dimension) const
  {
    vertex->set_dimension(dimension);
  }
/// @}

/// \name Queries on the identifier of the face complex including triangulation cells, facets, and vertices.
/// @{
  /** returns the index of vertex \p v
  */
  Index index(const Vertex_handle& v) const { return v->index(); }

  /** returns the subdomain index of cell \p cell
  */
  Subdomain_index subdomain_index(const Cell_handle& cell) const
  {
    return cell->subdomain_index();
  }
  /** returns the surface index of facet \p f
  */
  Surface_patch_index surface_patch_index(const Facet& f) const
  {
    return surface_patch_index(f.first, f.second);
  }

  /** returns the surface index of facet(\p cell, \p i)
  */
  Surface_patch_index surface_patch_index(const Cell_handle& cell,
                                          const int i) const
  {
    return cell->surface_patch_index(i);
  }
  /** returns the dimension of the lowest dimensional face of the input 3D
  * complex that contains the vertex
  */
  int in_dimension(const Vertex_handle& v) const { return v->in_dimension(); }

  /**
  * returns the curve index of edge \p e
  */
  Curve_index curve_index(const Edge& e) const
  {
    return curve_index(e.first->vertex(e.second),
                       e.first->vertex(e.third));
  }

  /**
  * returns the curve index of the edge formed by \p v1 and \p v2
  */
  Curve_index curve_index(const Vertex_handle& v1,
                          const Vertex_handle& v2) const
  {
    return curve_index(make_internal_edge(v1, v2));
  }

  /**
   * returns the corner index of vertex \p v
   */
  Corner_index corner_index(const Vertex_handle& v) const
  {
    typename Corner_map::const_iterator it = corners_.find(v);
    if (corners_.end() != it) { return it->second; }
    return Corner_index();
  }
/// @}

#ifndef CGAL_NO_DEPRECATED_CODE
  CGAL_DEPRECATED
    Curve_index curve_segment_index(const Edge& e) const
  {
    return curve_index(e);
  }

  CGAL_DEPRECATED
    Curve_index curve_segment_index(const Vertex_handle& v1,
      const Vertex_handle& v2) const
  {
    return curve_index(v1, v2);
  }
#endif // CGAL_NO_DEPRECATED_CODE

  std::size_t number_of_far_points() const
  {
    return far_vertices_.size();
  }

  void add_far_point(const Point &p)
  {
    far_vertices_.push_back(triangulation().insert(p));
  }

  void add_far_point(Vertex_handle vh)
  {
    far_vertices_.push_back(vh);
  }


  void remove_isolated_vertex(Vertex_handle v)
  {
    Triangulation& tr = triangulation();

    std::vector<Cell_handle> new_cells;
    new_cells.reserve(32);
    tr.remove_and_give_new_cells(v, std::back_inserter(new_cells));

    typename std::vector<Cell_handle>::iterator nc_it = new_cells.begin();
    typename std::vector<Cell_handle>::iterator nc_it_end = new_cells.end();
    for (; nc_it != nc_it_end; ++nc_it)
    {
      Cell_handle c = *nc_it;
      for (int i = 0; i < 4; ++i)
      {
        Facet mirror_facet = tr.mirror_facet(std::make_pair(c, i));
        if (is_in_complex(mirror_facet))
        {
          set_surface_patch_index(c, i,
            surface_patch_index(mirror_facet));
          c->set_facet_surface_center(i,
            mirror_facet.first->get_facet_surface_center(mirror_facet.second));
        }
      }
      /*int i_inf;
      if (c->has_vertex(tr.infinite_vertex(), i_inf))
      {
        Facet mirror_facet = tr.mirror_facet(std::make_pair(c, i_inf));
        if (is_in_complex(mirror_facet))
        {
          set_surface_patch_index(c, i_inf,
                                  surface_patch_index(mirror_facet));
        }
      }*/
    }
  }

  void remove_far_points()
  {
    //triangulation().remove(far_vertices_.begin(), far_vertices_.end());
    typename Far_vertices_vec::const_iterator it = far_vertices_.begin();
    typename Far_vertices_vec::const_iterator it_end = far_vertices_.end();
    for (; it != it_end; ++it)
    {
      remove_isolated_vertex(*it);
    }
    far_vertices_.clear();
  }

  /*!
    The tetrahedral mesh generation algorithm implemented in
    \link make_mesh_3() `CGAL::make_mesh_3()`\endlink
    and \link refine_mesh_3() `CGAL::refine_mesh_3()`\endlink
    does not guarantee that all the points inserted
    by the algorithm are actually present in the final mesh.

    In most cases, all points are used, but if the geometry of the object
    has small features compared to the size of the simplices (triangles and tetrahedra),
    it might be that the Delaunay facets that are selected in the restricted Delaunay
    triangulation miss some vertices of the triangulation.
    The concurrent version of the tetrahedral mesh generation algorithm
    also inserts a small set of auxiliary vertices that belong to the triangulation
    but are isolated from the complex at the end of the meshing process.

    This function removes these so-called \em isolated vertices, that belong to the
    triangulation but not to any simplex of the `C3T3`, from the triangulation.
  */
  void remove_isolated_vertices()
  {
    Triangulation& tr = triangulation();
    for (Vertex_handle v : tr.finite_vertex_handles())
      v->set_meshing_info(0);

    for (Cell_handle c : this->cells_in_complex())
    {
      for (int i = 0; i < 4; ++i)
      {
        Vertex_handle vi = c->vertex(i);
        vi->set_meshing_info(vi->meshing_info() + 1);
      }
    }

    for (Facet f :this->facets_in_complex())
    {
          for (int i = 1; i < 4; ++i)
      {
        Vertex_handle vi = f.first->vertex((f.second + i) % 4);
        vi->set_meshing_info(vi->meshing_info() + 1);
      }
    }

    std::vector<Vertex_handle> isolated;
    for (Vertex_handle v : tr.finite_vertex_handles())
    {
      if (v->meshing_info() == 0.
        && (v->in_dimension() > 1 || v->in_dimension() < 0))
        isolated.push_back(v);
    }

#ifdef CGAL_MESH_3_VERBOSE
    std::cout << "Remove " << isolated.size() << " isolated vertices...";
    std::cout.flush();
#endif

    CGAL_assertion(far_vertices_.size() <= isolated.size());
    far_vertices_.clear();

    for (Vertex_handle v : isolated)
      remove_isolated_vertex(v);

#ifdef CGAL_MESH_3_VERBOSE
    std::cout << "\nRemove " << isolated.size() << " isolated vertices done." << std::endl;
#endif
  }

  void rescan_after_load_of_triangulation();

/// \name Queries on the faces of the embedded complex
/// @{
  /**
  * returns the number of cells which belong to the 3D complex
  */
  size_type number_of_cells_in_complex() const { return number_of_cells_; }
  /**
  * returns the number of cells which belong to the 3D complex
  */
  size_type number_of_cells() const
  {
    return number_of_cells_in_complex();
  }
  /**
  * returns the number of surface facets of the complex
  */
  size_type number_of_facets_in_complex() const { return number_of_facets_; }
  /**
  * returns the number of surface facets of the complex
  */
  size_type number_of_facets() const
  {
    return number_of_facets_in_complex();
  }
  /**
   * returns the number of edges of the complex
   */
  size_type number_of_edges_in_complex() const
  {
    return edges_.size();
  }
  /**
   * returns the number of edges of the complex
   */
  size_type number_of_edges() const
  {
    return edges_.size();
  }
  /**
   * returns the number of corners of the complex
   */
  size_type number_of_vertices_in_complex() const
  {
    return corners_.size();
  }
  /**
   * returns the number of corners of the complex
   */
  size_type number_of_corners() const
  {
    return corners_.size();
  }
  /**
  * returns \c true if cell \p cell belongs to the 3D complex
  */
  bool is_in_complex(const Cell_handle& cell) const
  {
    return !(subdomain_index(cell) == Subdomain_index());
  }
  /** returns true if facet \p facet belongs to the 2D complex
  */
  bool is_in_complex(const Facet& facet) const
  {
    return is_in_complex(facet.first, facet.second);
  }

  /** returns true if facet (\p cell, \p i) belongs to the 2D complex
  */
  bool is_in_complex(const Cell_handle& cell, const int i) const
  {
    return (cell->is_facet_on_surface(i));
  }
  /**
   * returns true if edge \p e belongs to the 1D complex
   */
  bool is_in_complex(const Edge& e) const
  {
    return is_in_complex(e.first->vertex(e.second), e.first->vertex(e.third));
  }

  /**
   * returns true if edge (v1,v2) belongs to the 1D complex
   */
  bool is_in_complex(const Vertex_handle& v1, const Vertex_handle& v2) const
  {
    return is_in_complex(make_internal_edge(v1,v2));
  }

  /**
   * returns true if \p v is a 0-dimensionnal feature in the complex
   */
  bool is_in_complex(const Vertex_handle& v) const
  {
    return (corners_.find(v) != corners_.end());
  }
/// @}


  /// \name I/O Functions
  /// @{
  /**
   * outputs the outer boundary of the entire domain, with facets oriented outward.
   */
  std::ostream& output_boundary_to_off(std::ostream& out) const
  {
    internal::output_boundary_of_c3t3_to_off(*this, 0, out, false);
    return out;
  }

  /**
   * outputs the outer boundary of the selected subdomain, with facets oriented outward.
   */
  std::ostream& output_boundary_to_off(std::ostream& out, Subdomain_index subdomain) const
  {
    output_boundary_of_c3t3_to_off(*this, subdomain, out);
    return out;
  }

  /**
   * outputs the surface facets, with a consistent orientation at the interface of two subdomains.
   */
  std::ostream& output_facets_in_complex_to_off(std::ostream& out) const
  {
    internal::output_facets_in_complex_to_off(*this, out);
    return out;
  }

  /*!
  outputs the mesh to `os`
  in Medit format.
  */
#ifdef DOXYGEN_RUNNING
  void output_to_medit(std::ostream& os) const
#else
  void output_to_medit(std::ostream& os,
                       bool rebind = true,
                       bool show_patches = false) const
#endif
  {
    // Call global function
    bool all_vertices = true;
    bool all_cells = false;
    CGAL::IO::output_to_medit(os, *this, rebind, show_patches,
      all_vertices, all_cells);
  }

  void output_to_maya(std::ostream& os, bool surfaceOnly = true) const
  {
    // Call global function
    CGAL::IO::output_to_maya(os, *this, surfaceOnly);
  }

  /// @}

  /**
   * fills \p out with incident edges (1-dimensional features of \p v).
   * OutputIterator value type is std::pair<Vertex_handle,Curve_index>
   * \pre v->in_dimension() < 2
   */
  template <typename OutputIterator>
  OutputIterator
  adjacent_vertices_in_complex(const Vertex_handle& v, OutputIterator out) const;

  // -----------------------------------
  // Undocumented
  // -----------------------------------

  bool is_valid(bool verbose = false) const;

  // -----------------------------------
  // Complex traversal
  // -----------------------------------
private:
  /**
 * @class Cell_not_in_complex
 * @brief A class to filter cells which do not belong to the complex
 */
  class Cell_not_in_complex
  {
    const Self* r_self_;
    Subdomain_index index_;//needed by SWIG, should be const Subdomain_index
  public:
    Cell_not_in_complex() {}//needed by SWIG
    Cell_not_in_complex(const Self& self,
      const Subdomain_index& index = Subdomain_index())
      : r_self_(&self)
      , index_(index) { }

    bool operator()(Cell_handle ch) const
    {
      if (index_ == Subdomain_index()) { return !r_self_->is_in_complex(ch); }
      else { return !(r_self_->subdomain_index(ch) == index_); }
    }
  }; // end class Cell_not_in_complex

  typedef SMDS_3::internal::Iterator_not_in_complex<Self> Iterator_not_in_complex;

  class Facet_iterator_not_in_complex
  {
    const Self* c3t3_;
    Surface_patch_index index_; //need by SWIG: should be const Surface_patch_index
  public:
    Facet_iterator_not_in_complex() {} //need by SWIG
    Facet_iterator_not_in_complex(const Self& c3t3,
      const Surface_patch_index& index = Surface_patch_index())
      : c3t3_(&c3t3)
      , index_(index) { }

    template <typename Iterator>
    bool operator()(Iterator it) const
    {
      if (index_ == Surface_patch_index()) { return !c3t3_->is_in_complex(*it); }
      else { return !(c3t3_->surface_patch_index(*it) == index_); }
    }
  };

  class Edge_iterator_not_in_complex
  {
    const Self& c3t3_;
    const Curve_index index_;
  public:
    Edge_iterator_not_in_complex(const Self& c3t3,
                                 const Curve_index& index = Curve_index())
    : c3t3_(c3t3)
    , index_(index) { }

    template <typename Iterator>
    bool operator()(Iterator it) const
    {
      if ( index_ == Curve_index() ) { return ! c3t3_.is_in_complex(*it); }
      else { return c3t3_.curve_index(*it) != index_;  }
    }
  };

  class Vertex_iterator_not_in_complex
  {
    const Self& c3t3_;
    const Corner_index index_;
  public:
    Vertex_iterator_not_in_complex(const Self& c3t3,
                                   const Corner_index& index = Corner_index())
    : c3t3_(c3t3)
    , index_(index) { }

    template <typename ItMap>
    bool operator()(const ItMap it) const
    {
      if ( index_ == Corner_index() ) { return false; }
      else { return it->second != index_;  }
    }
  };

  // Filtered iterator
  typedef Filter_iterator<
    typename Corner_map::const_iterator,
    Vertex_iterator_not_in_complex >            Vertex_map_filter_iterator;

  // Iterator type to get the first element of pair
  typedef boost::transform_iterator <
    SMDS_3::internal::First_of<typename Vertex_map_filter_iterator::value_type>,
    Vertex_map_filter_iterator >                Vertex_map_iterator_first;

  // Iterator type to remove a level of referencing
  class Vertex_map_iterator_first_dereference
    : public boost::iterator_adaptor <
        Vertex_map_iterator_first_dereference,
        Vertex_map_iterator_first,
        typename Vertex_map_iterator_first::value_type::value_type,
        boost::use_default,
        typename Vertex_map_iterator_first::value_type::reference >
  {
    typedef Vertex_map_iterator_first_dereference Self;
    typedef boost::iterator_adaptor <
        Vertex_map_iterator_first_dereference,
        Vertex_map_iterator_first,
        typename Vertex_map_iterator_first::value_type::value_type,
        boost::use_default,
        typename Vertex_map_iterator_first::value_type::reference > iterator_adaptor_;
  public:
    typedef typename  Vertex_map_iterator_first::reference  pointer;
    typedef typename iterator_adaptor_::reference           reference;

    Vertex_map_iterator_first_dereference() : Self::iterator_adaptor_() { }

    template < typename Iterator >
    Vertex_map_iterator_first_dereference(Iterator i)
      : Self::iterator_adaptor_(typename Self::iterator_adaptor_::base_type(i))
    { }

    pointer operator->() const { return *(this->base()); }
    reference operator*() const { return **(this->base()); }

    operator const Vertex_handle&() const { return *(this->base()); }
  };

public:
/// \name Traversal of the complex
/// @{
#ifdef DOXYGEN_RUNNING
  /// Iterator type to visit the cells of the 3D complex
  typedef unspecified_type Cells_in_complex_iterator;
  /// Iterator type to visit the facets of the 2D complex
  typedef unspecified_type Facets_in_complex_iterator;
  /// Iterator type to visit the edges of the 1D complex
  typedef unspecified_type Edges_in_complex_iterator;
  /// Iterator type to visit the vertices of the 0D complex
  typedef unspecified_type Vertices_in_complex_iterator;

  /// Range type for iterating over all cells of the 3D complex,
  /// with a nested type iterator that has as value type `Cell_handle`.
  typedef Iterator_range<unspecified_type> Cells_in_complex;
  /// Range type for iterating over all facets of the 2D complex,
  /// with a nested type iterator that has as value type `Facet`.
  typedef Iterator_range<unspecified_type> Facets_in_complex;
  /// Range type for iterating over all cells of the 1D complex,
  /// with a nested type iterator that has as value type `Edge`.
  typedef Iterator_range<unspecified_type> Edges_in_complex;
  /// Range type for iterating over all vertices of the 0D complex,
  /// with a nested type iterator that has as value type `Vertex_handle`.
  typedef Iterator_range<unspecified_type> Vertices_in_complex;

#else
  typedef Filter_iterator<
    typename Triangulation::Finite_edges_iterator,
    Edge_iterator_not_in_complex >              Edges_in_complex_iterator;
  typedef Vertex_map_iterator_first_dereference Vertices_in_complex_iterator;
  typedef Filter_iterator<
    typename Triangulation::Finite_facets_iterator,
    Facet_iterator_not_in_complex >               Facets_in_complex_iterator;

  /**
 * @class Cells_in_complex_iterator
 * @brief Iterator type to visit the cells of the triangulation that belong
 * to the 3D complex
 *
 * This class is useful to ensure that Cells_in_complex_iterator is convertible
 * to Cell_handle
 */
  class Cells_in_complex_iterator :
    public Filter_iterator<typename Triangulation::Finite_cells_iterator,
                           Cell_not_in_complex>
  {
  private:
    typedef typename Triangulation::Finite_cells_iterator Tr_iterator;
    typedef Filter_iterator<typename Triangulation::Finite_cells_iterator,
      Cell_not_in_complex> Base;
    typedef Cells_in_complex_iterator Self;

  public:
    Cells_in_complex_iterator() : Base() { }
    Cells_in_complex_iterator(Base i) : Base(i) { }

    Self& operator++() { Base::operator++(); return *this; }
    Self& operator--() { Base::operator--(); return *this; }
    Self operator++(int) { Self tmp(*this); ++(*this); return tmp; }
    Self operator--(int) { Self tmp(*this); --(*this); return tmp; }

    operator const Cell_handle&() const { return this->base(); }
  }; // end class Cells_in_complex_iterator

  typedef Iterator_range<Prevent_deref<Vertices_in_complex_iterator, const Vertex_handle&>> Vertices_in_complex;
  typedef Iterator_range<Edges_in_complex_iterator> Edges_in_complex;
  typedef Iterator_range<Facets_in_complex_iterator> Facets_in_complex;
  typedef Iterator_range<Prevent_deref<Cells_in_complex_iterator, const Cell_handle&>> Cells_in_complex;

#endif

/// \name Iterators
/// @{

  /// returns a \c Cells_in_complex_iterator to the first cell of the 3D complex
  Cells_in_complex_iterator cells_in_complex_begin() const
  {
    return CGAL::filter_iterator(tr_.finite_cells_end(),
                                 Cell_not_in_complex(*this),
                                 tr_.finite_cells_begin());
  }

  /// returns a \c Cells_in_complex_iterator to the first cell of the 3D complex
  Cells_in_complex_iterator cells_in_complex_begin(const Subdomain_index& index) const
  {
    return CGAL::filter_iterator(tr_.finite_cells_end(),
                                 Cell_not_in_complex(*this, index),
                                 tr_.finite_cells_begin());
  }

  /// returns the past-the-end iterator for the cells of the 3D complex
  Cells_in_complex_iterator cells_in_complex_end() const
  {
    return CGAL::filter_iterator(tr_.finite_cells_end(),
                                 Cell_not_in_complex(*this));
  }

  /// returns a `Facets_in_complex_iterator` to the first facet of the 2D complex
  Facets_in_complex_iterator facets_in_complex_begin() const
  {
    return CGAL::filter_iterator(tr_.finite_facets_end(),
                                 Facet_iterator_not_in_complex(*this),
                                 tr_.finite_facets_begin());
  }

  /// returns a `Facets_in_complex_iterator` to the first facet of the 2D complex
  Facets_in_complex_iterator
    facets_in_complex_begin(const Surface_patch_index& index) const
  {
    return CGAL::filter_iterator(tr_.finite_facets_end(),
      Facet_iterator_not_in_complex(*this, index),
      tr_.finite_facets_begin());
  }

  /// returns past-the-end iterator on facet of the 2D complex
  Facets_in_complex_iterator facets_in_complex_end(const Surface_patch_index = Surface_patch_index()) const
  {
    return CGAL::filter_iterator(tr_.finite_facets_end(),
      Facet_iterator_not_in_complex(*this));
  }

  /// returns a `Edges_in_complex_iterator` to the first edge of the 1D complex
  Edges_in_complex_iterator edges_in_complex_begin() const
  {
    return CGAL::filter_iterator(this->triangulation().finite_edges_end(),
                                 Edge_iterator_not_in_complex(*this),
                                 this->triangulation().finite_edges_begin());
  }

  /// returns a `Edges_in_complex_iterator` to the first edge of the 1D complex
  Edges_in_complex_iterator
  edges_in_complex_begin(const Curve_index& index) const
  {
    return CGAL::filter_iterator(this->triangulation().finite_edges_end(),
                                 Edge_iterator_not_in_complex(*this,index),
                                 this->triangulation().finite_edges_begin());
  }

  /// returns past-the-end iterator on edges of the 1D complex
  Edges_in_complex_iterator edges_in_complex_end(const Curve_index& = Curve_index()) const
  {
    return CGAL::filter_iterator(this->triangulation().finite_edges_end(),
                                 Edge_iterator_not_in_complex(*this));
  }

  /// returns a `Vertices_in_complex_iterator` to the first vertex of the 0D complex
  Vertices_in_complex_iterator vertices_in_complex_begin() const
  {
    return CGAL::filter_iterator(corners_.end(),
                                 Vertex_iterator_not_in_complex(*this),
                                 corners_.begin());
  }

  /// returns a `Vertices_in_complex_iterator` to the first vertex of the 0D complex
  Vertices_in_complex_iterator
  vertices_in_complex_begin(const Corner_index& index) const
  {
    return CGAL::filter_iterator(corners_.end(),
                                 Vertex_iterator_not_in_complex(*this,index),
                                 corners_.begin());
  }

  /// returns past-the-end iterator on vertices of the 0D complex
  Vertices_in_complex_iterator vertices_in_complex_end() const
  {
    return CGAL::filter_iterator(corners_.end(),
                                 Vertex_iterator_not_in_complex(*this));
  }

  /*!
    returns a range of iterators over vertices of the 0D complex
    \note The value type of `Vertices_in_complex::iterator` is `Vertex_handle`.
  */
  Vertices_in_complex vertices_in_complex() const
  {
      return { vertices_in_complex_begin(), vertices_in_complex_end() };
  }
  /*!
    returns a range of iterators over the edges of the 1D complex,
    starting at an arbitrary edge.
    Returns an empty range when `t.dimension() < 2`.
  */
  Edges_in_complex edges_in_complex() const
  {
    return Edges_in_complex(edges_in_complex_begin(),
                            edges_in_complex_end());
  }
  /*!
    returns a range of iterators over the facets of the 2D complex,
    starting at an arbitrary facet.
    Returns an empty range when `t.dimension() < 2`.
  */
  Facets_in_complex facets_in_complex() const
  {
    return Facets_in_complex(facets_in_complex_begin(),
                             facets_in_complex_end());
  }
  /*!
    returns a range of iterators over cells of the 3D complex.
    Returns an empty range when `triangulation().number_of_cells() == 0`
    or complex is empty.
    \note The value type of `Cells_in_complex::iterator` is `Cell_handle`.
  */
  Cells_in_complex cells_in_complex() const
  {
    return { cells_in_complex_begin(), cells_in_complex_end() };
  }
///  @}

public:
  template <typename Tr2, typename CoI, typename CuI>
  friend
  std::istream& operator>>(std::istream& is,
    Mesh_complex_3_in_triangulation_3<Tr2, CoI, CuI>& c3t3);

  static std::string io_signature()
  {
    return Get_io_signature<Tr>()();
  }

  /**
   * @cond SKIP_IN_MANUAL
   * creates an `Internal_edge` object (i.e a pair of ordered `Vertex_handle`)
   * @endcond
   */
  Internal_edge make_internal_edge(const Vertex_handle& v1,
                                   const Vertex_handle& v2) const
  {
    if ( v1 < v2 ) { return Internal_edge(v1,v2); }
    else { return Internal_edge(v2,v1); }
  }

  /**
   * @cond SKIP_IN_MANUAL
   * returns true if \p edge is in the 1D complex
   * @endcond
   */
  bool is_in_complex(const Internal_edge& edge) const
  {
    return (curve_index(edge) != Curve_index() );
  }

  /**
   * @cond SKIP_IN_MANUAL
   * adds edge \p edge to the 1D complex, with curve index `index`
   * @endcond
   */
  void add_to_complex(const Internal_edge& edge, const Curve_index& index)
  {
    CGAL_precondition(!is_in_complex(edge));
#if CGAL_MESH_3_PROTECTION_DEBUG & 1
    std::cerr << "Add edge ( " << disp_vert(edge.left)
              << " , " << disp_vert(edge.right) << " ), curve_index=" << index
              << " to c3t3.\n";
#endif // CGAL_MESH_3_PROTECTION_DEBUG
    std::pair<typename Edge_map::iterator, bool> it = edges_.insert(edge);
    it.first->info = index;
  }

  /**
   * @cond SKIP_IN_MANUAL
   * removes edge \p edge from the 1D complex
   * @endcond
   */
  void remove_from_complex(const Internal_edge& edge)
  {
    edges_.erase(edge);
  }

  /**
  * @cond SKIP_IN_MANUAL
  * returns the curve index of edge \p edge
  * @endcond
  */
  Curve_index curve_index(const Internal_edge& edge) const
  {
    typename Edge_map::const_iterator it = edges_.find(edge);
    if ( edges_.end() != it ) { return it->info; }
    return Curve_index();
  }

  /// @cond SKIP_IN_MANUAL
  /// Returns `NOT_IN_COMPLEX`, `BOUNDARY`, `REGULAR`, or `SINGULAR`,
  /// depending on the number of incident facets in the complex, and the
  /// number of connected components of its link
  /// @endcond
  Face_status face_status(const Vertex_handle v) const
  {
    if (!manifold_info_initialized_) init_manifold_info();
    const std::size_t n = v->cached_number_of_incident_facets();

    if (n == 0) return NOT_IN_COMPLEX;

    //test incident edges for REGULARITY and count BOUNDARY edges
    typename std::vector<Edge> edges;
    edges.reserve(64);
    if (tr_.is_parallel()) {
      tr_.incident_edges_threadsafe(v, std::back_inserter(edges));
    }
    else {
      tr_.incident_edges(v, std::back_inserter(edges));
    }
    int number_of_boundary_incident_edges = 0; // could be a bool
    for (typename std::vector<Edge>::iterator
      eit = edges.begin(), end = edges.end();
      eit != end; eit++)
    {
      switch (face_status(*eit))
      {
      case NOT_IN_COMPLEX: case REGULAR: break;
      case BOUNDARY: ++number_of_boundary_incident_edges; break;
      default:
#ifdef CGAL_MESHES_DEBUG_REFINEMENT_POINTS
        std::cerr << "singular edge...\n";
        std::cerr << tr_.point(v) << std::endl;
#endif // CGAL_MESHES_DEBUG_REFINEMENT_POINTS
        return SINGULAR;
      }
    }

    // From here all incident edges (in complex) are REGULAR or BOUNDARY.
    const std::size_t nb_components = union_find_of_incident_facets(v);
    if (nb_components > 1) {
#ifdef CGAL_MESHES_DEBUG_REFINEMENT_POINTS
      std::cerr << "singular vertex: nb_components=" << nb_components << std::endl;
      std::cerr << tr_.point(v) << std::endl;
#endif // CGAL_MESHES_DEBUG_REFINEMENT_POINTS
      return SINGULAR;
    }
    else { // REGULAR OR BOUNDARY
#ifdef CGAL_MESHES_DEBUG_REFINEMENT_POINTS
      std::cerr << "regular or boundary: " << tr_.point(v) << std::endl;
#endif // CGAL_MESHES_DEBUG_REFINEMENT_POINTS
      if (number_of_boundary_incident_edges != 0)
        return BOUNDARY;
      else
        return REGULAR;
    }
  }

  /// @cond SKIP_IN_MANUAL
  /// This function should be called only when incident edges
  /// are known to be `REGULAR` or `BOUNDARY`
  /// @endcond
  bool is_regular_or_boundary_for_vertices(Vertex_handle v) const {
    return union_find_of_incident_facets(v) == 1;
  }

  /// @cond SKIP_IN_MANUAL
  /// Returns `NOT_IN_COMPLEX`, `BOUNDARY`, `REGULAR`, or `SINGULAR`,
  /// depending on the number of incident facets in the complex
  /// @endcond
  Face_status face_status(const Edge& edge) const
  {
    if (!manifold_info_initialized_) init_manifold_info();

#ifdef CGAL_LINKED_WITH_TBB
    typename Edge_facet_counter::const_accessor accessor;
    if (!edge_facet_counter_.find(accessor,
      this->make_ordered_pair(edge)))
      return NOT_IN_COMPLEX;
    switch (accessor->second)
#else // not CGAL_LINKED_WITH_TBB
    switch (edge_facet_counter_[this->make_ordered_pair(edge)])
#endif // not CGAL_LINKED_WITH_TBB
    {
    case 0: return NOT_IN_COMPLEX;
    case 1: return BOUNDARY;
    case 2: return REGULAR;
    default: return SINGULAR;
    }
  }

  /// Returns true if the vertex \p v has is incident to at least a facet
  /// of the complex
  bool has_incident_facets_in_complex(const Vertex_handle& v) const
  {
    if (!manifold_info_initialized_) init_manifold_info();
    return v->cached_number_of_incident_facets() > 0;
  }

 /**
  * @cond SKIP_IN_MANUAL
  * @brief inserts \p [first,last[ in the triangulation (with dimension 2)
  * @param first the iterator on the first point to insert
  * @param last the iterator past the last point to insert
  *
  * InputIterator value type must be \c std::pair<Tr::Point,Index>
  * @endcond
  */
  template <typename InputIterator>
  void insert_surface_points(InputIterator first, InputIterator last)
  {
    typename Tr::Geom_traits::Construct_weighted_point_3 cwp =
      tr_.geom_traits().construct_weighted_point_3_object();

    while (first != last)
    {
      Vertex_handle vertex = tr_.insert(cwp((*first).first));
      vertex->set_index((*first).second);
      vertex->set_dimension(2);
      ++first;
    }
  }

  /**
  * @cond SKIP_IN_MANUAL
  * @brief inserts \p [first,last[ in the triangulation (with dimension 2 and
  * index \p default_index)
  * @param first the iterator on the first point to insert
  * @param last the iterator past the last point to insert
  * @param default_index the index to be used to insert points
  *
  * InputIterator value type must be \c Tr::Point
  * @endcond
  */
  template <typename InputIterator>
  void insert_surface_points(InputIterator first,
                             InputIterator last,
                             const Index& default_index)
  {
    typename Tr::Geom_traits::Construct_weighted_point_3 cwp =
      tr_.geom_traits().construct_weighted_point_3_object();

    while (first != last)
    {
      Vertex_handle vertex = tr_.insert(cwp(*first));
      vertex->set_index(default_index);
      vertex->set_dimension(2);
      ++first;
    }
  }

  void clear_cells_and_facets_from_c3t3()
  {
    //clear cells
    for (typename Tr::All_cells_iterator cit = this->triangulation().all_cells_begin();
         cit != this->triangulation().all_cells_end();
         ++cit)
    {
      set_subdomain_index(cit, Subdomain_index());
    }
    this->number_of_cells_ = 0;

    //clear facets
    for (typename Tr::All_facets_iterator fit = this->triangulation().all_facets_begin();
         fit != this->triangulation().all_facets_end();
         ++fit)
    {
      const auto& facet = *fit;
      set_surface_patch_index(facet.first, facet.second, Surface_patch_index());
      if (this->triangulation().dimension() > 2) {
        const Facet& mirror = tr_.mirror_facet(facet);
        set_surface_patch_index(mirror.first, mirror.second, Surface_patch_index());
      }
    }
    this->number_of_facets_ = 0;

    //clear manifold info
    clear_manifold_info();
  }

  void clear_manifold_info() {
    edge_facet_counter_.clear();
    manifold_info_initialized_ = false;
  }

  /** @cond SKIP_IN_MANUAL
  * Returns bbox
  * @endcond
  */
  Bbox_3 bbox() const;

private:
  // Sequential: non-atomic
  // "dummy" is here to allow the specialization (see below)
  // See https://groups.google.com/group/comp.lang.c++.moderated/browse_thread/thread/285ab1eec49e1cb6
  template<typename Concurrency_tag2, typename dummy = void>
  struct Number_of_elements
  {
    typedef size_type type;
  };

  template<typename Concurrency_tag2, typename dummy = void>
  struct Init_number_of_elements
  {
    template<typename T>
    void operator()(T& a, const T& b)
    {
      a = b;
    }
    template<typename T>
    void operator()(T& a)
    {
      a = 0;
    }
  };

  template<typename Concurrency_tag2, typename dummy = void>
  struct Swap_elements
  {
    template<typename T>
    void operator()(T& a, T& b)
    {
      std::swap(a, b);
    }
  };
#ifdef CGAL_LINKED_WITH_TBB
  // Parallel: atomic
  template<typename dummy>
  struct Number_of_elements<Parallel_tag, dummy>
  {
    typedef std::atomic<size_type> type;
  };

  template<typename dummy>
  struct Init_number_of_elements<Parallel_tag, dummy>
  {
    template<typename T>
    void operator()(T& a, const T& b)
    {
      a = b.load();
    }
    template<typename T>
    void operator()(T& a)
    {
      a = 0;
    }
  };

  template<typename dummy>
  struct Swap_elements<Parallel_tag, dummy>
  {
    template<typename T>
    void operator()(T& a, T& b)
    {
      T tmp;
      tmp.exchange(a);
      a.exchange(b);
      b.exchange(tmp);
    }
  };
#endif // CGAL_LINKED_WITH_TBB

private:
  void init_manifold_info() const
  {
    for (typename Tr::All_vertices_iterator
      vit = triangulation().finite_vertices_begin(),
      end = triangulation().finite_vertices_end();
      vit != end; ++vit)
    {
      vit->set_c2t3_cache(0, (std::numeric_limits<size_type>::max)());
    }

    edge_facet_counter_.clear();

    for (typename Tr::Finite_facets_iterator
      fit = triangulation().finite_facets_begin(),
      end = triangulation().finite_facets_end();
      fit != end; ++fit)
    {
      if (is_in_complex(*fit)) {
        const Cell_handle cell = fit->first;
        const int i = fit->second;
        for (int j = 0; j < 3; ++j)
        {
          const int edge_index_va = tr_.vertex_triple_index(i, j);
          const int edge_index_vb = tr_.vertex_triple_index(i, (j == 2) ? 0 : (j + 1));
          const Vertex_handle edge_va = cell->vertex(edge_index_va);
          const Vertex_handle edge_vb = cell->vertex(edge_index_vb);
#ifndef CGAL_LINKED_WITH_TBB
          ++edge_facet_counter_[this->make_ordered_pair(edge_va, edge_vb)];
#else // CGAL_LINKED_WITH_TBB
          {
            typename Edge_facet_counter::accessor accessor;
            edge_facet_counter_.insert(accessor,
              this->make_ordered_pair(edge_va, edge_vb));
            ++accessor->second;
          }
#endif // CGAL_LINKED_WITH_TBB

          const std::size_t n = edge_va->cached_number_of_incident_facets();
          edge_va->set_c2t3_cache(n + 1, (std::numeric_limits<size_type>::max)());
        }
      }
    }
    manifold_info_initialized_ = true;
  }

  /// Extract the subset `F` of facets of the complex incident to `v` and
  /// return the number of connected component of the adjacency graph of `F`.
  std::size_t union_find_of_incident_facets(const Vertex_handle v) const
  {
    if (v->is_c2t3_cache_valid())
    {
      const std::size_t n = v->cached_number_of_components();
      if (n != (std::numeric_limits<size_type>::max)()) return n;
    }

    Union_find<Facet> facets;
    { // fill the union find
      std::vector<Facet> non_filtered_facets;
      if (tr_.is_parallel()) {
        tr_.incident_facets_threadsafe(v, std::back_inserter(non_filtered_facets));
      }
      else {
        tr_.incident_facets(v, std::back_inserter(non_filtered_facets));
      }

      for (typename std::vector<Facet>::iterator
        fit = non_filtered_facets.begin(),
        end = non_filtered_facets.end();
        fit != end; ++fit)
      {
        if (is_in_complex(*fit)) facets.push_back(*fit);
      }
    }

    typedef boost::unordered_map<Vertex_handle,
                                 typename Union_find<Facet>::handle,
                                 Hash_fct>    Vertex_set_map;
    typedef typename Vertex_set_map::iterator Vertex_set_map_iterator;

    Vertex_set_map vsmap;

    for (typename Union_find<Facet>::iterator
         it = facets.begin(), end = facets.end();
         it != end; ++it)
    {
      const Cell_handle& ch = (*it).first;
      const int& i = (*it).second;
      for (int j = 0; j < 3; ++j) {
        const Vertex_handle w = ch->vertex(tr_.vertex_triple_index(i, j));
        if (w != v) {
          Vertex_set_map_iterator vsm_it = vsmap.find(w);
          if (vsm_it != vsmap.end()) {
            facets.unify_sets(vsm_it->second, it);
          }
          else {
            vsmap.insert(std::make_pair(w, it));
          }
        }
      }
    }
    const std::size_t nb_components = facets.number_of_sets();

    const std::size_t n = v->cached_number_of_incident_facets();
    v->set_c2t3_cache(n, nb_components);
    return nb_components;
  }

public:
  // -----------------------------------
  // Backward Compatibility
  // -----------------------------------
#ifndef DOXYGEN_RUNNING
#ifndef CGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX
  typedef Surface_patch_index   Surface_index;

  void set_surface_index(const Facet& f, const Surface_index& index)
  {
    set_surface_patch_index(f, index);
  }

  void set_surface_index(const Cell_handle& c, const int i, const Surface_index& index)
  {
    set_surface_patch_index(c, i, index);
  }

  Surface_index surface_index(const Facet& f) const
  {
    return surface_patch_index(f);
  }

  Surface_index surface_index(const Cell_handle& c, const int i) const
  {
    return surface_patch_index(c, i);
  }
#endif // CGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX

#ifndef CGAL_MESH_3_NO_DEPRECATED_C3T3_ITERATORS
  typedef Facets_in_complex_iterator  Facet_iterator;
  typedef Cells_in_complex_iterator   Cell_iterator;

  Facet_iterator facets_begin() const
  {
    return facets_in_complex_begin();
  }

  Facet_iterator facets_end() const
  {
    return facets_in_complex_end();
  }

  Cell_iterator cells_begin() const
  {
    return cells_in_complex_begin();
  }

  Cell_iterator cells_end() const
  {
    return cells_in_complex_end();
  }
#endif // CGAL_MESH_3_NO_DEPRECATED_C3T3_ITERATORS
#endif // DOXYGEN_RUNNING
  // -----------------------------------
  // End backward Compatibility
  // -----------------------------------


private:
  // Private data members
  Triangulation tr_;

  typedef typename Base::Pair_of_vertices Pair_of_vertices;
#ifdef CGAL_LINKED_WITH_TBB
  typedef tbb::concurrent_hash_map<Pair_of_vertices, int,
                                   Hash_compare_for_TBB> Edge_facet_counter;
#else // not CGAL_LINKED_WITH_TBB
  typedef std::map<Pair_of_vertices, int> Edge_facet_counter;
#endif // not CGAL_LINKED_WITH_TBB

  mutable Edge_facet_counter edge_facet_counter_;

  typename Number_of_elements<Concurrency_tag>::type number_of_facets_;
  typename Number_of_elements<Concurrency_tag>::type number_of_cells_;

  mutable bool manifold_info_initialized_;

  Edge_map edges_;
  Corner_map corners_;
  Far_vertices_vec far_vertices_;
};

template <typename Tr, typename CI_, typename CSI_>
Mesh_complex_3_in_triangulation_3<Tr, CI_, CSI_>::
Mesh_complex_3_in_triangulation_3()
  : Base()
  , tr_()
  , edge_facet_counter_() //TODO: parallel!
  , manifold_info_initialized_(false) //TODO: parallel!
{
  // We don't put it in the initialization list because
  // std::atomic has no constructors
  number_of_facets_ = 0;
  number_of_cells_ = 0;
}

template <typename Tr, typename CI_, typename CSI_>
Mesh_complex_3_in_triangulation_3<Tr,CI_,CSI_>::
Mesh_complex_3_in_triangulation_3(const Self& rhs)
  : Base(rhs)
  , tr_(rhs.tr_)
  , edge_facet_counter_(rhs.edge_facet_counter_)
  , manifold_info_initialized_(rhs.manifold_info_initialized_)
  , edges_()
  , corners_()
{
  Init_number_of_elements<Concurrency_tag> init;
  init(number_of_facets_, rhs.number_of_facets_);
  init(number_of_cells_, rhs.number_of_cells_);

  // Copy edges
  for ( typename Edge_map::const_iterator it = rhs.edges_.begin(),
       end = rhs.edges_.end() ; it != end ; ++it )
  {
    const Vertex_handle& va = it->right;
    const Vertex_handle& vb = it->left;

    Vertex_handle new_va;
    this->triangulation().is_vertex(rhs.triangulation().point(va), new_va);

    Vertex_handle new_vb;
    this->triangulation().is_vertex(rhs.triangulation().point(vb), new_vb);

    this->add_to_complex(make_internal_edge(new_va,new_vb), it->info);
  }

  // Copy corners
  for ( typename Corner_map::const_iterator it = rhs.corners_.begin(),
       end = rhs.corners_.end() ; it != end ; ++it )
  {
    Vertex_handle new_v;
    this->triangulation().is_vertex(rhs.triangulation().point(it->first), new_v);
    this->add_to_complex(new_v, it->second);
  }

  // Parse vertices to identify far vertices
  if (rhs.far_vertices_.size() > 0)
  {
    Triangulation &tr = triangulation();
    typename Tr::Finite_vertices_iterator vit = tr.finite_vertices_begin();
    for(typename Tr::Finite_vertices_iterator end = tr.finite_vertices_end();
        vit != end ; ++vit)
    {
      if (vit->in_dimension() == -1)
        far_vertices_.push_back(vit);
    }
    CGAL_assertion(far_vertices_.size() == rhs.far_vertices_.size());
  }
}

template <typename Tr, typename CI_, typename CSI_>
Mesh_complex_3_in_triangulation_3<Tr, CI_, CSI_>::
Mesh_complex_3_in_triangulation_3(Self&& rhs)
  : Base()
  , tr_(std::move(rhs.tr_))
  , edge_facet_counter_(std::move(rhs.edge_facet_counter_))
  , manifold_info_initialized_(std::exchange(rhs.manifold_info_initialized_, false))
  , edges_(std::move(rhs.edges_))
  , corners_(std::move(rhs.corners_))
  , far_vertices_(std::move(rhs.far_vertices_))
{
  Init_number_of_elements<Concurrency_tag> init;
  init(number_of_facets_, rhs.number_of_facets_);
  init(number_of_cells_, rhs.number_of_cells_);
  init(rhs.number_of_facets_); // set to 0
  init(rhs.number_of_cells_); // set to 0
}

template <typename Tr, typename CI_, typename CSI_>
template <typename OutputIterator>
OutputIterator
Mesh_complex_3_in_triangulation_3<Tr,CI_,CSI_>::
adjacent_vertices_in_complex(const Vertex_handle& v, OutputIterator out) const
{
  CGAL_precondition(v->in_dimension() < 2);

  typedef typename Edge_map::right_const_iterator Rcit;
  typedef typename Edge_map::left_const_iterator Lcit;

  // Add edges containing v is on the left
  std::pair<Rcit,Rcit> range_right = edges_.right.equal_range(v);
  for ( Rcit rit = range_right.first ; rit != range_right.second ; ++rit )
  {
    *out++ = std::make_pair(rit->second, rit->info);
  }

  // Add edges containing v on the right
  std::pair<Lcit,Lcit> range_left = edges_.left.equal_range(v);
  for ( Lcit lit = range_left.first ; lit != range_left.second ; ++lit )
  {
    *out++ = std::make_pair(lit->second, lit->info);
  }

  return out;
}


template <typename Tr, typename CI_, typename CSI_>
bool
Mesh_complex_3_in_triangulation_3<Tr,CI_,CSI_>::
is_valid(bool verbose) const
{
  typedef typename Tr::Weighted_point                         Weighted_point;
  typedef boost::unordered_map<Vertex_handle, int, Hash_fct>  Vertex_map;

  Vertex_map vertex_map;

  // Fill map counting neighbor number for each vertex of an edge
  for ( typename Edge_map::const_iterator it = edges_.begin(),
       end = edges_.end() ; it != end ; ++it )
  {
    const Vertex_handle& v1 = it->right;
    if ( vertex_map.find(v1) == vertex_map.end() ) { vertex_map[v1] = 1; }
    else { vertex_map[v1] += 1; }

    const Vertex_handle& v2 = it->left;
    if ( vertex_map.find(v2) == vertex_map.end() ) { vertex_map[v2] = 1; }
    else { vertex_map[v2] += 1; }
  }

  // Verify that each vertex has 2 neighbors if it's not a corner
  for ( typename Vertex_map::iterator vit = vertex_map.begin(),
       vend = vertex_map.end() ; vit != vend ; ++vit )
  {
    if ( vit->first->in_dimension() != 0 && vit->second != 2 )
    {
      if(verbose)
        std::cerr << "Validity error: vertex " << (void*)(&*vit->first)
                  << " (" << this->triangulation().point(vit->first) << ") "
                  << "is not a corner (dimension " << vit->first->in_dimension()
                  << ") but has " << vit->second << " neighbor(s)!\n";
      return false;
    }
  }

  // Verify that balls of each edge intersect
  for ( typename Edge_map::const_iterator it = edges_.begin(),
       end = edges_.end() ; it != end ; ++it )
  {
    typename Tr::Geom_traits::Compute_weight_3 cw =
      this->triangulation().geom_traits().compute_weight_3_object();
    typename Tr::Geom_traits::Construct_point_3 cp =
      this->triangulation().geom_traits().construct_point_3_object();
    typename Tr::Geom_traits::Construct_sphere_3 sphere =
      this->triangulation().geom_traits().construct_sphere_3_object();
    typename Tr::Geom_traits::Do_intersect_3 do_intersect =
      this->triangulation().geom_traits().do_intersect_3_object();

    const Weighted_point& itrwp = this->triangulation().point(it->right);
    const Weighted_point& itlwp = this->triangulation().point(it->left);

    if ( ! do_intersect(sphere(cp(itrwp), cw(itrwp)), sphere(cp(itlwp), cw(itlwp))) )
    {
      std::cerr << "Points p[" << disp_vert(it->right) << "], dim=" << it->right->in_dimension()
                << " and q[" << disp_vert(it->left) << "], dim=" << it->left->in_dimension()
                << " form an edge but do not intersect !\n";
      return false;
    }
  }

  return true;
}

template <typename Tr, typename CI_, typename CSI_>
void
Mesh_complex_3_in_triangulation_3<Tr, CI_, CSI_>::
add_to_complex(const Cell_handle& cell,
               const int i,
               const Surface_patch_index& index)
{
  CGAL_precondition(!(index == Surface_patch_index()));

  if (!is_in_complex(cell, i))
  {
    Facet mirror = tr_.mirror_facet(std::make_pair(cell, i));
    set_surface_patch_index(cell, i, index);
    set_surface_patch_index(mirror.first, mirror.second, index);
    ++number_of_facets_;
    if (manifold_info_initialized_) {
      for (int j = 0; j < 3; ++j)
      {
        int edge_index_va = tr_.vertex_triple_index(i, j);
        int edge_index_vb = tr_.vertex_triple_index(i, (j == 2) ? 0 : (j + 1));
        Vertex_handle edge_va = cell->vertex(edge_index_va);
        Vertex_handle edge_vb = cell->vertex(edge_index_vb);
#ifdef CGAL_LINKED_WITH_TBB
        {
          typename Edge_facet_counter::accessor accessor;
          edge_facet_counter_.insert(accessor,
            this->make_ordered_pair(edge_va, edge_vb));
          ++accessor->second;
        }
#else // not CGAL_LINKED_WITH_TBB
        ++edge_facet_counter_[this->make_ordered_pair(edge_va, edge_vb)];
#endif // not CGAL_LINKED_WITH_TBB

        const std::size_t n = edge_va->cached_number_of_incident_facets();
        const std::size_t m = edge_va->cached_number_of_components();
        edge_va->set_c2t3_cache(n + 1, m);
      }
      const int dimension_plus_1 = tr_.dimension() + 1;
      // update c2t3 for vertices of f
      for (int j = 0; j < dimension_plus_1; j++) {
        if (j != i) {
#ifdef CGAL_MESHES_DEBUG_REFINEMENT_POINTS
          if (cell->vertex(j)->is_c2t3_cache_valid())
            std::cerr << "(" << tr_.point(cell, j) << ")->invalidate_c2t3_cache()\n";
#endif // CGAL_MESHES_DEBUG_REFINEMENT_POINTS
          cell->vertex(j)->invalidate_c2t3_cache();
        }
      }
    }
  }
}

template <typename Tr, typename CI_, typename CSI_>
void
Mesh_complex_3_in_triangulation_3<Tr, CI_, CSI_>::
remove_from_complex(const Facet& facet)
{
  if (is_in_complex(facet))
  {
    Facet mirror = tr_.mirror_facet(facet);
    set_surface_patch_index(facet.first, facet.second, Surface_patch_index());
    set_surface_patch_index(mirror.first, mirror.second, Surface_patch_index());
    --number_of_facets_;
    if (manifold_info_initialized_) {
      const Cell_handle cell = facet.first;
      const int i = facet.second;
      for (int j = 0; j < 3; ++j)
      {
        const int edge_index_va = tr_.vertex_triple_index(i, j);
        const int edge_index_vb = tr_.vertex_triple_index(i, (j == 2) ? 0 : (j + 1));
        const Vertex_handle edge_va = cell->vertex(edge_index_va);
        const Vertex_handle edge_vb = cell->vertex(edge_index_vb);
#ifdef CGAL_LINKED_WITH_TBB
        {
          typename Edge_facet_counter::accessor accessor;
          edge_facet_counter_.insert(accessor,
            this->make_ordered_pair(edge_va, edge_vb));
          --accessor->second;
        }
#else // not CGAL_LINKED_WITH_TBB
        --edge_facet_counter_[this->make_ordered_pair(edge_va, edge_vb)];
#endif // not CGAL_LINKED_WITH_TBB

        const std::size_t n = edge_va->cached_number_of_incident_facets();
        CGAL_assertion(n > 0);
        const std::size_t m = edge_va->cached_number_of_components();
        edge_va->set_c2t3_cache(n - 1, m);
      }
      const int dimension_plus_1 = tr_.dimension() + 1;
      // update c2t3 for vertices of f
      for (int j = 0; j < dimension_plus_1; j++) {
        if (j != facet.second) {
#ifdef CGAL_MESHES_DEBUG_REFINEMENT_POINTS
          if (cell->vertex(j)->is_c2t3_cache_valid())
            std::cerr << "(" << tr_.point(cell, j) << ")->invalidate_c2t3_cache()\n";
#endif // CGAL_MESHES_DEBUG_REFINEMENT_POINTS
          cell->vertex(j)->invalidate_c2t3_cache();
        }
      }
    }
  }
}

template <typename Tr, typename CI_, typename CSI_>
Bbox_3
Mesh_complex_3_in_triangulation_3<Tr, CI_, CSI_>::
bbox() const
{
  if (0 == triangulation().number_of_vertices())
  {
    return Bbox_3();
  }

  typename Tr::Finite_vertices_iterator vit = tr_.finite_vertices_begin();
  Bbox_3 result = tr_.point(vit++).bbox();

  for (typename Tr::Finite_vertices_iterator end = tr_.finite_vertices_end();
    vit != end; ++vit)
  {
    result = result + tr_.point(vit).bbox();
  }

  return result;
}

template <typename Tr, typename CI_, typename CSI_>
void
Mesh_complex_3_in_triangulation_3<Tr,CI_,CSI_>::
rescan_after_load_of_triangulation()
{
  corners_.clear();
  for(typename Tr::Finite_vertices_iterator
        vit = this->triangulation().finite_vertices_begin(),
        end = this->triangulation().finite_vertices_end();
      vit != end; ++vit)
  {
    if ( vit->in_dimension() == 0 ) {
      add_to_complex(vit, Corner_index(1));
    }
  }

  this->number_of_facets_ = 0;
  for (typename Tr::Finite_facets_iterator
    fit = this->triangulation().finite_facets_begin(),
    end = this->triangulation().finite_facets_end();
    fit != end; ++fit)
  {
    if (this->is_in_complex(*fit)) {
      ++this->number_of_facets_;
    }
  }

  this->number_of_cells_ = 0;
  for (typename Tr::Finite_cells_iterator
    cit = this->triangulation().finite_cells_begin(),
    end = this->triangulation().finite_cells_end();
    cit != end; ++cit)
  {
    if (this->is_in_complex(cit)) {
      ++this->number_of_cells_;
    }
  }
}

template <typename Tr, typename CI_, typename CSI_>
std::ostream &
operator<< (std::ostream& os,
            const Mesh_complex_3_in_triangulation_3<Tr,CI_,CSI_> &c3t3)
{
  // TODO: implement edge saving
  return os << c3t3.triangulation();
}


template <typename Tr, typename CI_, typename CSI_>
std::istream &
operator>> (std::istream& is,
            Mesh_complex_3_in_triangulation_3<Tr,CI_,CSI_> &c3t3)
{
  // TODO: implement edge loading
  c3t3.clear();
  is >> c3t3.triangulation();

  if (!is) {
    c3t3.clear();
    return is;
  }

  c3t3.rescan_after_load_of_triangulation();
  return is;
}

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_MESH_COMPLEX_3_IN_TRIANGULATION_3_H
