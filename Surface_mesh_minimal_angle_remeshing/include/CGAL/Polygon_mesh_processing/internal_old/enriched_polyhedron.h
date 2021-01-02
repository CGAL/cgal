// Copyright (c) 2019  INRIA Sophia-Antipolis (France).
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
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Kaimo Hu

#ifndef ENRICHEDPOLYHEDRON_H_
#define ENRICHEDPOLYHEDRON_H_

// CGAL
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Triangulation_hierarchy_2.h>
#include <CGAL/Timer.h>
#include <fstream>
// local
#include "bvd.h"
#include "random.h"
#include "HalfedgeDS_decorator.h"
#include "Polyhedron_3.h"

#define DOUBLE_MAX 1000000.0
#define DOUBLE_MIN -1000000.0
#define MAX_VALUE 10000
#define MIN_VALUE 0.0001  // specified for numerical stability
#define SQUARED_MIN_VALUE 0.00000001

enum SampleNumberStrategy {
  // #samples per facet is roughtly fixed (with respect to the sample strategy)
  k_fixed = 0,
  // #samples per facet is variable with respect to size_of_facets()
  k_variable
};

enum SampleStrategy {
  k_uniform = 0,  // #samples per facet is proportional to its area
  k_adaptive      // #samples per facet is roughly the same
};

enum OptimizeType {
  k_none = 0,
  k_input_to_remesh,
  k_remesh_to_input,
  k_both
};

enum OptimizeStrategy {
  k_approximation = 0,
  k_Interpolation
};

enum EdgeFlipStrategy {
  k_improve_valence = 0,
  k_improve_angle
};

enum RelocateStrategy {
  k_barycenter = 0,
  k_cvt_barycenter
};

enum VertexType {
  k_feature_vertex = 0,
  k_crease_vertex,
  k_smooth_vertex
};

enum DrawType {
  k_polyhedron = 0,
  k_all_voronoi,
  k_vertex_voronoi,
  k_edge_voronoi,
  k_facet_voronoi
};

enum RenderType {
  k_plain_facets = 0, // facets with plain color
  k_ifi_facets,       // facets with interpolated feature intensity colors
  k_mr_facets,
  k_classifications,  // vertex type (feature, crease, smooth)
  k_gaussian_curvature,
  k_maximal_halfedge_dihedral,
  k_normal_dihedral,
  k_feature_intensity,
  k_capacity,
  k_weight
};

typedef CGAL::Bbox_3 Bbox;
typedef CGAL::Simple_cartesian<double> Kernel;            // Kernel
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef std::pair<Point, Point> Point_pair;               // for local links
typedef std::pair<FT, Point_pair> Link;
typedef std::list<Link> Link_list;                        // for out links
typedef Link_list::iterator Link_list_iter;
typedef Link_list::const_iterator Link_list_const_iter;
typedef std::list<Link_list_iter> Link_iter_list;         // for in links
typedef Link_iter_list::iterator Link_iter_list_iter;
typedef Link_iter_list::const_iterator Link_iter_list_const_iter;
typedef std::list<Link*> Link_pointer_list;
typedef std::list<Link*>::iterator Link_pointer_iter;
typedef std::list<Link*>::const_iterator Link_pointer_const_iter;

namespace CGAL {
  // a refined facet with a normal
  template<class Refs, class T, class P, class Vector>
  class Enriched_facet : public CGAL::HalfedgeDS_face_base<Refs, T> {
   public:
    typedef Vector Normal;

    // life cycle
    Enriched_facet() {
      m_tag = 0;
      m_normal = CGAL::NULL_VECTOR;
      m_max_squared_error = 0.0;
    }
    // tag
    int& tag() { return m_tag; }
    const int& tag() const { return m_tag; }
    // normal
    Normal& normal() { return m_normal; }
    const Normal& normal() const { return m_normal; }
    // max squared error
    FT& max_squared_error() { return m_max_squared_error; }
    const FT& max_squared_error() const { return m_max_squared_error; }
    // links
    Link_list& facet_out_links() { return m_facet_out_links; }
    const Link_list& facet_out_links() const { return m_facet_out_links; }
    Link_iter_list& facet_in_links() { return m_facet_in_links; }
    const Link_iter_list& facet_in_links() const { return m_facet_in_links; }
    Link_iter_list& edge_in_links() { return m_edge_in_links; }
    const Link_iter_list& edge_in_links() const { return m_edge_in_links; }
    Link_pointer_list& vertex_in_links() { return m_vertex_in_links; }
    const Link_pointer_list& vertex_in_links() const
        { return m_vertex_in_links; }

   private:
    int m_tag;                        // general-purpose tag
    Normal m_normal;
    FT m_max_squared_error;           // the maximum squared error
    // links
    Link_list m_facet_out_links;      // point links from this facet
    Link_iter_list m_facet_in_links;  // point links to this facet
    Link_iter_list m_edge_in_links;
    Link_pointer_list m_vertex_in_links;
  };

  // a refined halfedge with a general tag
  template <class Refs, class Tprev, class Tvertex, class Tface, class Norm>
  class Enriched_halfedge :
    public CGAL::HalfedgeDS_halfedge_base<Refs, Tprev, Tvertex, Tface> {
   public:
    // life cycle
    Enriched_halfedge() {
      m_tag = 0;
      // 0.0: flat, 1.0: the maximum; -1.0 means data is in opposite
      m_normal_dihedral = -1.0;
      m_is_crease = false;
      //m_is_feature = false;
    }
    // tag
    int& tag() { return m_tag; }
    const int& tag() const { return m_tag; }
    // normal_dihedral
    FT& normal_dihedral() { return m_normal_dihedral; }
    const FT& normal_dihedral() const { return m_normal_dihedral; }
    // is crease
    bool& is_crease() { return m_is_crease; }
    const bool& is_crease() const { return m_is_crease; }
    // feature
    //bool& is_feture() { return m_is_feature; }
    //const bool& is_feature() const { return m_is_feature; }
    // links
    Link_list& edge_out_links() { return m_edge_out_links; }
    const Link_list& edge_out_links() const { return m_edge_out_links; }

   private:
    int m_tag;              // general-purpose tag
    FT m_normal_dihedral;   // the normalized face normal dihedral
    bool m_is_crease;       // boundary is also regarded as crease
    //bool m_is_feature;    // mark whether it is a feature edge (not used now)
    // links
    Link_list m_edge_out_links;         // point link from this edge
  };

  // a refined vertex with a normal and a tag
  template <class Refs, class T, class P, class Vector>
  class Enriched_vertex : public CGAL::HalfedgeDS_vertex_base<Refs, T, P> {
   public:
    typedef Vector Normal;

    // life cycle
    Enriched_vertex() {
      m_tag = 0;
      m_max_halfedge_dihedral = 0.0;
      m_gaussian_curvature = 0.0;
      //m_vertex_type = VertexType::k_smooth_vertex;
      //m_is_feature = false;
    }
    Enriched_vertex(const P &pt) :
      CGAL::HalfedgeDS_vertex_base<Refs, T, P>(pt) {
      m_tag = 0;
      m_max_halfedge_dihedral = 0.0;
      m_gaussian_curvature = 0.0;
      //m_vertex_type = VertexType::k_smooth_vertex;
      //m_is_feature = false;
    }
    // tag
    int& tag() { return m_tag; }
    const int& tag() const { return m_tag; }
    // normal
    Normal& normal() { return m_normal; }
    const Normal& normal() const { return m_normal; }
    // max halfedge dihedral
    FT& max_halfedge_dihedral() { return m_max_halfedge_dihedral; }
    const FT& max_halfedge_dihedral() const { return m_max_halfedge_dihedral; }
    // gaussian curvature
    FT& gaussian_curvature() { return m_gaussian_curvature; }
    const FT& gaussian_curvature() const { return m_gaussian_curvature; }
    // feature
    FT feature_intensity() const {
      // range: [1, (PI + 1) ^ 2]
      return (m_gaussian_curvature + 1) * (m_max_halfedge_dihedral + 1);
      // range: [0, (PI + 1) ^ 2 - 1]
      //return (m_gaussian_curvature + 1) * (m_max_halfedge_dihedral + 1) - 1;
    }
    //VertexType& vertex_type() { return m_vertex_type; }
    //const VertexType& vertex_type() const { return m_vertex_type; }
    // feature
    //bool& is_feature() { return m_is_feature; }
    //const bool& is_feature() const { return m_is_feature; }
    // links
    Link& vertex_out_link() { return m_vertex_out_link; }
    const Link& vertex_out_link() const { return m_vertex_out_link; }

   private:
    int m_tag;
    Normal m_normal;
    FT m_max_halfedge_dihedral; // maximum halfedge weights around it
    FT m_gaussian_curvature;    // normalized gaussian curvature
    //VertexType m_vertex_type;
    //bool m_is_feature;
    // links
    Link m_vertex_out_link;                 // Links from this vertex
  };

  // the combination of the refined items
  struct Enriched_item : public CGAL::Polyhedron_items_3 {
    // wrap vertex
    template<class Refs, class Traits>
    struct Vertex_wrapper {
      typedef typename Traits::Point_3 Point;
      typedef typename Traits::Vector_3 Normal;
      typedef Enriched_vertex<Refs, CGAL::Tag_true, Point, Normal> Vertex;
    };
    // wrap face
    template<class Refs, class Traits>
    struct Face_wrapper {
      typedef typename Traits::Point_3 Point;
      typedef typename Traits::Vector_3 Normal;
      typedef Enriched_facet<Refs, CGAL::Tag_true, Point, Normal> Face;
    };
    // wrap halfedge
    template<class Refs, class Traits>
    struct Halfedge_wrapper {
      typedef typename Traits::Vector_3 Normal;
      typedef Enriched_halfedge<Refs, CGAL::Tag_true, CGAL::Tag_true,
        CGAL::Tag_true, Normal> Halfedge;
    };
  };

  // compute facet normal
  struct Facet_normal { 	// (functor)
    template <class	Facet>
    void operator()(Facet& f) {
      typename Facet::Normal sum = CGAL::NULL_VECTOR;
      typename Facet::Halfedge_around_facet_circulator h = f.facet_begin();
      do {
        typename Facet::Normal normal = CGAL::cross_product(
          h->next()->vertex()->point() - h->vertex()->point(),
          h->next()->next()->vertex()->point() - h->next()->vertex()->point());
        FT sqnorm = normal * normal;
        if (sqnorm != 0) {
          normal = normal / std::sqrt(sqnorm);
        }
        sum = sum + normal;
      } while (++h != f.facet_begin());
      FT	sqnorm = sum * sum;
      if (sqnorm != 0.0) {
        f.normal() = sum / std::sqrt(sqnorm);
      }
      else {
        f.normal() = CGAL::NULL_VECTOR;
        std::cerr << std::endl << red << "degenerated facet " << white;
      }
    }
  };

  // compute vertex	normal
  struct Vertex_normal {   //	(functor)
    template <class	Vertex>
    void operator()(Vertex&	v) {
      typename Vertex::Normal	normal = CGAL::NULL_VECTOR;
      typename Vertex::Halfedge_around_vertex_const_circulator	he = v.vertex_begin();
      typename Vertex::Halfedge_around_vertex_const_circulator	begin = he;
      CGAL_For_all(he, begin) {
        if (!he->is_border()) {
          normal = normal + he->facet()->normal();
        }
      }
      double	sqnorm = normal * normal;
      if (sqnorm != 0.0) {
        v.normal() = normal / (double)std::sqrt(sqnorm);
      }
      else {
        v.normal() = CGAL::NULL_VECTOR;
      }
    }
  };

  //***************************************************************************
  template <class Kernel, class Items>
  class Enriched_polyhedron : public CGAL::Polyhedron_3<Kernel, Items> {
  public:
    using typename CGAL::Polyhedron_3<Kernel, Items>::Halfedge_handle;
    using typename CGAL::Polyhedron_3<Kernel, Items>::Halfedge_const_handle;
    using typename CGAL::Polyhedron_3<Kernel, Items>::Halfedge_const_iterator;
    using typename CGAL::Polyhedron_3<Kernel, Items>::Halfedge_iterator;
    using typename CGAL::Polyhedron_3<Kernel, Items>::Vertex_const_handle;
    using typename CGAL::Polyhedron_3<Kernel, Items>::Vertex_const_iterator;
    using typename CGAL::Polyhedron_3<Kernel, Items>::Vertex_handle;
    using typename CGAL::Polyhedron_3<Kernel, Items>::Vertex_iterator;
    using typename CGAL::Polyhedron_3<Kernel, Items>::Facet_handle;
    using typename CGAL::Polyhedron_3<Kernel, Items>::Facet_const_handle;
    using typename CGAL::Polyhedron_3<Kernel, Items>::Facet_const_iterator;
    using typename CGAL::Polyhedron_3<Kernel, Items>::Edge_iterator;
    using typename CGAL::Polyhedron_3<Kernel, Items>::Edge_const_iterator;
    using typename CGAL::Polyhedron_3<Kernel, Items>::Facet_iterator;
    using typename CGAL::Polyhedron_3<Kernel, Items>::Halfedge_around_vertex_circulator;
    using typename CGAL::Polyhedron_3<Kernel, Items>::Halfedge_around_vertex_const_circulator;
    using typename CGAL::Polyhedron_3<Kernel, Items>::Halfedge_around_facet_circulator;
    using typename CGAL::Polyhedron_3<Kernel, Items>::Halfedge_around_facet_const_circulator;
    using typename CGAL::Polyhedron_3<Kernel, Items>::Point_const_iterator;
    using typename CGAL::Polyhedron_3<Kernel, Items>::Point_iterator;
    using CGAL::Polyhedron_3<Kernel, Items>::facets_begin;
    using CGAL::Polyhedron_3<Kernel, Items>::facets_end;
    using CGAL::Polyhedron_3<Kernel, Items>::vertices_begin;
    using CGAL::Polyhedron_3<Kernel, Items>::vertices_end;
    using CGAL::Polyhedron_3<Kernel, Items>::halfedges_begin;
    using CGAL::Polyhedron_3<Kernel, Items>::halfedges_end;
    using CGAL::Polyhedron_3<Kernel, Items>::points_begin;
    using CGAL::Polyhedron_3<Kernel, Items>::points_end;
    using CGAL::Polyhedron_3<Kernel, Items>::edges_begin;
    using CGAL::Polyhedron_3<Kernel, Items>::edges_end;
    using CGAL::Polyhedron_3<Kernel, Items>::size_of_vertices;
    using CGAL::Polyhedron_3<Kernel, Items>::size_of_halfedges;
    using CGAL::Polyhedron_3<Kernel, Items>::size_of_facets;
    using CGAL::Polyhedron_3<Kernel, Items>::join_vertex;
    using CGAL::Polyhedron_3<Kernel, Items>::join_facet;
    using CGAL::Polyhedron_3<Kernel, Items>::split_edge;
    using CGAL::Polyhedron_3<Kernel, Items>::split_facet;

    // 2D BVD
    typedef CGAL::Triangulation_vertex_base_2<Kernel> Vbb;
    typedef CGAL::Triangulation_hierarchy_vertex_base_2<Vbb> Vb;
    typedef CGAL::Triangulation_face_base_2<Kernel> Fb;
    typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
    // 3D BVD
    typedef CBvd<Kernel, Tds> Bvd;
    typedef typename Bvd::Vector_3 Vector;
    typedef typename Bvd::Vector_3 Normal;
    typedef typename Bvd::Line_3 Line;
    typedef typename Bvd::Segment_3 Segment;
    typedef typename Bvd::Triangle_3 Triangle;
    typedef typename Bvd::Plane_3 Plane;
    // Point list
    typedef typename Bvd::Point_list Point_list;
    typedef typename Bvd::Point_iter Point_iter;
    typedef typename Bvd::Point_const_iter Point_const_iter;
    // Color list
    typedef typename Bvd::Color_list Color_list;
    typedef typename Bvd::Color_iter Color_iter;
    typedef typename Bvd::Color_const_iter Color_const_iter;
    // Element list
    typedef typename std::list<Halfedge_handle> Halfedge_list;
    typedef typename std::list<Halfedge_handle>::iterator Halfedge_iter;
    typedef typename std::list<Halfedge_const_handle> Halfedge_const_list;
    typedef typename std::list<Halfedge_const_handle>::const_iterator
      Halfedge_const_iter;
    typedef std::list<Facet_handle> Facet_list;
    typedef typename std::list<Facet_handle>::iterator Facet_iter;
    typedef std::list<Facet_const_handle> Facet_const_list;
    typedef typename std::list<Facet_const_handle>::const_iterator
      Facet_const_iter;
    typedef std::list<Vertex_handle> Vertex_list;
    typedef typename std::list<Vertex_handle>::iterator Vertex_iter;
    typedef std::list<Vertex_const_handle> Vertex_const_list;
    typedef typename std::list<Vertex_const_handle>::const_iterator
      Vertex_const_iter;


    // compare functions
    struct Point_Comp {
      bool operator()(const Point &a, const Point &b) {
        Vector vec = a - b;
        if (vec * vec < SQUARED_MIN_VALUE) {
          return false;
        }
        else {
          return a < b;
        }
      }
    };

    // life cycle
    Enriched_polyhedron() {}
    virtual ~Enriched_polyhedron() {}

    // normals
    void calculate_normals(const std::string &name) {
      // normals: per facets, then per vertex
      std::cout << "Computing " << name << " normals...";
      std::for_each(facets_begin(), facets_end(), Facet_normal());
      std::for_each(vertices_begin(), vertices_end(), Vertex_normal());
      std::cout << "Done" << std::endl;
    }

    void calculate_local_normals(std::set<Facet_handle> *facets,
                                 std::set<Vertex_handle> *vertices) {
      for (auto it = facets->begin(); it != facets->end(); ++it) {
        calculate_facet_normal(*it);
      }
      for (auto it = vertices->begin(); it != vertices->end(); ++it) {
        calculate_vertex_normal(*it);
      }
    }

    // max errors
    void calculate_max_squared_errors() {
      reset_facet_max_squared_errors(-1.0);
      update_facet_in_max_squared_errors();
      update_facet_out_max_squared_errors();
      update_edge_in_max_squared_errors();
      update_edge_out_max_squared_errors();
      update_vertex_in_max_squared_errors();
      update_vertex_out_max_squared_errors();
    }

    void calculate_max_squared_errors(std::set<Facet_handle> *facets) const {
      for (auto fit = facets->begin(); fit != facets->end(); ++fit) {
        FT max_se = 0.0, se = 0.0;
        Facet_handle fh = *fit;
        // facet in links
        for (Link_iter_list_iter lit = fh->facet_in_links().begin();
          lit != fh->facet_in_links().end(); ++lit) {
          Link_list_iter it = *lit;
          se = CGAL::squared_distance(it->second.first, it->second.second);
          max_se = CGAL::max(max_se, se);
        }
        // edge in links
        for (Link_iter_list_iter lit = fh->edge_in_links().begin();
          lit != fh->edge_in_links().end(); ++lit) {
          Link_list_iter it = *lit;
          se = CGAL::squared_distance(it->second.first, it->second.second);
          max_se = CGAL::max(max_se, se);
        }
        // vertex in links
        for (Link_pointer_iter lit = fh->vertex_in_links().begin();
          lit != fh->vertex_in_links().end(); ++lit) {
          Link *link = *lit;
          se = CGAL::squared_distance(link->second.first, link->second.second);
          max_se = CGAL::max(max_se, se);
        }
        // facet out links
        for (Link_list_iter it = fh->facet_out_links().begin();
          it != fh->facet_out_links().end(); ++it) {
          se = CGAL::squared_distance(it->second.first, it->second.second);
          max_se = CGAL::max(max_se, se);
        }
        Halfedge_handle hh = fh->halfedge();
        for (int i = 0; i <= 2; ++i) {
          // edge out links
          Halfedge_handle hi = hh;
          if (hi->normal_dihedral() == -1.0) {
            hi = hi->opposite();
          }
          for (Link_list_iter it = hi->edge_out_links().begin();
            it != hi->edge_out_links().end(); ++it) {
            se = CGAL::squared_distance(it->second.first, it->second.second);
            max_se = CGAL::max(max_se, se);
          }
          // vertex out links
          Vertex_handle vh = hh->vertex();
          const Link &link = vh->vertex_out_link();
          se = CGAL::squared_distance(link.second.first, link.second.second);
          max_se = CGAL::max(max_se, se);
          hh = hh->next();
        }
        fh->max_squared_error() = max_se;
      }
    }

    Halfedge_handle get_maximal_error(FT *max_error) {
      Facet_iterator fi = facets_begin();
      Facet_handle max_error_facet = fi;
      FT max_se = max_error_facet->max_squared_error();
      ++fi;
      for (; fi != facets_end(); ++fi) {
        if (max_se < fi->max_squared_error()) {
          max_se = fi->max_squared_error();
          max_error_facet = fi;
        }
      }
      *max_error = CGAL::sqrt(max_se);
      return get_longest_halfedge(max_error_facet);
    }

    Halfedge_handle get_local_maximal_error(
        const std::set<Facet_handle> &facets, FT *max_error) {
      assert(!facets.empty());
      typename std::set<Facet_handle>::iterator it = facets.begin();
      Facet_handle fh = *it;
      FT max_se = fh->max_squared_error();
      Facet_handle max_error_facet = fh;
      ++it;
      for (; it != facets.end(); ++it) {
        fh = *it;
        if (max_se < fh->max_squared_error()) {
          max_se = fh->max_squared_error();
          max_error_facet = fh;
        }
      }
      *max_error = CGAL::sqrt(max_se);
      return get_longest_halfedge(max_error_facet);
    }

    // feature intensities
    void calculate_feature_intensities(const std::string &name,
        FT dihedral_theta, FT dihedral_delta, FT sum_theta, FT sum_delta,
        bool inherit_element_types, FT feature_control_delta) {
      std::cout << "Computing " << name << " feature intensities...";
      calculate_edge_feature_intensities(dihedral_theta, dihedral_delta);
      calculate_vertex_feature_intensities(sum_theta, sum_delta);
      if (inherit_element_types) {
        calculate_edge_classifications(feature_control_delta);
      }
      std::cout << "Done" << std::endl;
    }

    void update_local_feature_intensity(Vertex_handle vh,
        bool reset_normal_dihedral, FT dihedral_theta, FT dihedral_delta,
        FT sum_theta, FT sum_delta) {
      Vertex_list incident_vertices;
      incident_vertices.push_back(vh);
      // step 1: calcualte the edge feature intensities around it
      Halfedge_around_vertex_circulator vcirc = vh->vertex_begin();
      Halfedge_around_vertex_circulator vend = vcirc;
      CGAL_For_all(vcirc, vend) {
        calculate_edge_feature_intensity(vcirc, reset_normal_dihedral,
          dihedral_theta, dihedral_delta);
        incident_vertices.push_back(get_source_vertex(vcirc));
      }
      // step 2: calculate the vertex feature intensity
      for (Vertex_iter it = incident_vertices.begin();
        it != incident_vertices.end(); ++it) {
        calculate_vertex_feature_intensity(*it, sum_theta, sum_delta);
      }
    }

    // samples and links
    void clear_out_links() {
      // step 1: clear out links in facets
      for (Facet_iterator fi = facets_begin(); fi != facets_end(); ++fi) {
        fi->facet_out_links().clear();
      }
      // step 2: clear out links in halfedges
      for (Halfedge_iterator hi = halfedges_begin();
        hi != halfedges_end(); ++hi) {
        hi->edge_out_links().clear();
      }
      // step 3: clear out links in vertices (do not need to do anything)
    }

    void clear_in_link_iterators() {
      for (Facet_iterator fi = facets_begin(); fi != facets_end(); ++fi) {
        fi->facet_in_links().clear();
        fi->edge_in_links().clear();
        fi->vertex_in_links().clear();
      }
    }

    void clear_local_links(Halfedge_handle hh,
                           std::set<Facet_handle> *in_link_facets) const {
      // step 1: clear all in link_facets
      clear_local_in_links(in_link_facets);
      // step 2: clear the out links in facets that are incident to hh
      if (hh->normal_dihedral() == -1.0) {
        hh = hh->opposite();
      }
      hh->edge_out_links().clear();
      if (!hh->is_border()) {
        hh->facet()->facet_out_links().clear();
      }
      if (!hh->opposite()->is_border()) {
        hh->opposite()->facet()->facet_out_links().clear();
      }
    }

    void clear_local_links(Vertex_handle vh,
                           std::set<Facet_handle> *in_link_facets) const {
      // step 1: clear all in links in in_link_facets
      clear_local_in_links(in_link_facets);
      // step 2: clear the out links in one-ring facets of vh
      Halfedge_around_vertex_circulator vcirc = vh->vertex_begin();
      Halfedge_around_vertex_circulator vend = vcirc;
      CGAL_For_all(vcirc, vend) {
        Halfedge_handle hh = vcirc;
        if (!hh->is_border()) {                 // clear facet out links
          Facet_handle fh = hh->facet();
          fh->facet_out_links().clear();
        }
        if (hh->normal_dihedral() == -1.0) {    // clear edge out links
          hh = hh->opposite();
        }
        hh->edge_out_links().clear();
      }
    }

    void get_nb_samples_per_facet(int nb_samples_per_facets,
        int max_samples_per_area, int min_samples_per_triangle,
        SampleStrategy sample_strategy, const Facet_list &facets) {
      size_t nb_facets = facets.size();
      FT total_area = get_sum_area(facets);
      if (total_area < SQUARED_MIN_VALUE) {   // degenerated case
        for (auto it = facets.begin(); it != facets.end(); ++it) {
          Facet_handle fh = *it;
          fh->tag() = 0;
        }
        return;
      }
      if (sample_strategy == SampleStrategy::k_uniform) {
        // the number of samples per facet is proportional to its area
        size_t nb_samples = nb_samples_per_facets * nb_facets;
        for (auto it = facets.begin(); it != facets.end(); ++it) {
          Facet_handle fh = *it;
          FT facet_area = area(fh);
          int nb_facet_samples =
            static_cast<int>(nb_samples * facet_area / total_area);
          int nb_max_samples = max_samples_per_area * facet_area;
          nb_facet_samples = std::min(nb_facet_samples, nb_max_samples);
          nb_facet_samples = std::max(nb_facet_samples,
            min_samples_per_triangle);
          fh->tag() = nb_facet_samples;
        }
      }
      else {  // SampleStrategy::k_adaptive, samples per facet is the same
        for (auto it = facets.begin(); it != facets.end(); ++it) {
          Facet_handle fh = *it;
          std::set<Facet_handle> incident_facets;
          get_incident_facets(fh, &incident_facets);
          FT sum_area = 0.0;
          for (auto it2 = incident_facets.begin();
            it2 != incident_facets.end(); ++it2) {
            sum_area += area(*it2);
          }
          FT facet_area = area(fh);
          int nb_facet_samples = static_cast<int>(incident_facets.size() *
            nb_samples_per_facets * facet_area / sum_area);
          int nb_max_samples = max_samples_per_area * facet_area;
          nb_facet_samples = std::min(nb_facet_samples, nb_max_samples);
          nb_facet_samples = std::max(nb_facet_samples,
            min_samples_per_triangle);
          fh->tag() = nb_facet_samples;
        }
      }
    }

    void generate_random_samples(bool use_stratified_sampling, int nb_samples,
        int bvd_iteration_count, Facet_const_handle fh,
        Point_list *inner_samples, std::list<double> *feature_weights) const {
      // step 1: generate the samples
      inner_samples->clear();
      feature_weights->clear();
      Halfedge_const_handle hh = fh->halfedge();
      const Point &a = get_target_vertex(hh)->point();
      const Point &b = get_opposite_vertex(hh)->point();
      const Point &c = get_source_vertex(hh)->point();
      if (nb_samples <= 0) {
        return;
      }
      else if (nb_samples == 1) {
        inner_samples->push_back(CGAL::centroid(a, b, c));
      }
      else if (nb_samples == 2) {
        FT bc = CGAL::squared_distance(b, c);
        FT ca = CGAL::squared_distance(c, a);
        FT ab = CGAL::squared_distance(a, b);
        Point d;
        int longest_edge_index;
        if (bc > ca) {
          longest_edge_index = bc > ab ? 0 : 2;
        }
        else {
          longest_edge_index = ab < ca ? 1 : 2;
        }
        switch (longest_edge_index) {
        case 0:		// split in edge bc
          d = midpoint(fh->halfedge()->prev());
          inner_samples->push_back(CGAL::centroid(a, b, d));
          inner_samples->push_back(CGAL::centroid(a, d, c));
          break;
        case 1:		// split in edge ca
          d = midpoint(fh->halfedge());
          inner_samples->push_back(CGAL::centroid(b, c, d));
          inner_samples->push_back(CGAL::centroid(b, d, a));
          break;
        case 2:		// split in edge ab
          d = midpoint(fh->halfedge()->next());
          inner_samples->push_back(CGAL::centroid(c, a, d));
          inner_samples->push_back(CGAL::centroid(c, d, b));
          break;
        default:
          break;
        }
      }
      else {
        // step 1: generate the unique inner samples
        std::set<Point, Point_Comp> samples;
        Random random;
        srand(time(NULL));
        Vector ab = b - a, ac = c - a;  // edge vectors
        while (samples.size() < nb_samples) {
          FT u = 0.0, v = 0.0;
          u = random.random_value<double>(0.0, 1.0);
          u = 0.9 * u + 0.05;
          while (v == 0.0 || u + v == 1.0) {
            v = random.random_value<double>(0.0, 1.0);
            v = 0.9 * v + 0.05;
          }
          if (u + v > 1.0) {    // flip over diag if needed
            u = 1.0 - u, v = 1.0 - v;
          }
          samples.insert(a + u * ab + v * ac);
        }
        // step 2: convert to Point_list data structure
        Point_list all_samples(samples.begin(), samples.end());
        // step 3: relocate them using BVD iteration
        for (int i = 0; i < bvd_iteration_count; ++i) {
          Bvd bvd(triangle(fh));
          bvd.run(all_samples);
        }
        // step 3: pick the inner samples
        inner_samples->clear();
        for (auto it = all_samples.begin(); it != all_samples.end(); ++it) {
          inner_samples->push_back(*it);
        }

        //backup: we use disturbed_border_samples if !use_stratified sampling
        // step 1£º generate the unique disturbed border samples if necessary
        //FT disturb_ratio = 0.01;
        //std::set<Point, Point_Comp> disturbed_border_samples;
        //if (!use_stratified_sampling) {
        //  get_disturbed_border_samples(fh, disturb_ratio,
        //                               &disturbed_border_samples);
        //}
        //// step 2: generate the unique inner samples
        //std::set<Point, Point_Comp> samples(disturbed_border_samples);
        //Random random;
        //srand(time(NULL));
        //Vector ab = b - a, ac = c - a;  // edge vectors
        //while (samples.size() <
        //    nb_samples + disturbed_border_samples.size()) {
        //  FT u = 0.0, v = 0.0;
        //  u = random.random_value<double>(0.0, 1.0);
        //  u = 0.9 * u + 0.05;
        //  while (v == 0.0 || u + v == 1.0) {
        //    v = random.random_value<double>(0.0, 1.0);
        //    v = 0.9 * v + 0.05;
        //  }
        //  if (u + v > 1.0) {    // flip over diag if needed
        //    u = 1.0 - u, v = 1.0 - v;
        //  }
        //  samples.insert(a + u * ab + v * ac);
        //}
        //// step 3: convert to Point_list data structure
        //Point_list all_samples(samples.begin(), samples.end());
        //// step 4: relocate them using BVD iteration
        //for (int i = 0; i < bvd_iteration_count; ++i) {
        //  Bvd bvd(triangle(fh));
        //  bvd.run(all_samples, disturbed_border_samples);
        //}
        //// step 5: pick the inner samples
        //inner_samples->clear();
        //std::set<Point, Point_Comp>::iterator sit;
        //for (Point_iter it = all_samples.begin();
        //  it != all_samples.end(); ++it) {
        //  sit = disturbed_border_samples.find(*it);
        //  if (sit == disturbed_border_samples.end()) {
        //    inner_samples->push_back(*it);
        //  }
        //}
        //if (all_samples.size() !=
        //  inner_samples->size() + disturbed_border_samples.size()) {
        //  std::cout << std::endl << "all samples: " << all_samples.size()
        //            << " inner samples: " << inner_samples->size()
        //            << " disturbed_border_samples: "
        //            << disturbed_border_samples.size() << std::endl;
        //}
      }
      // step 2: calcualte the feature weights
      FT fi_a = get_target_vertex(hh)->feature_intensity();
      FT fi_b = get_opposite_vertex(hh)->feature_intensity();
      FT fi_c = get_source_vertex(hh)->feature_intensity();
      FT feature_weight = 0.0;
      // version 1: each sample has the same feature weight (efficient)
      /*feature_weight = (fi_a + fi_b + fi_c) / 3.0;
      for (Point_iter it = inner_samples->begin();
        it != inner_samples->end(); ++it) {
        feature_weights->push_back(feature_weight);
      }*/
      // version 2: each sample has different feature weights (accurate)
      FT facet_area = area(fh);
      FT weight_a, weight_b, weight_c;
      for (Point_iter it = inner_samples->begin();
        it != inner_samples->end(); ++it) {
        const Point &p = *it;
        weight_a = area(p, b, c) / facet_area;
        weight_b = area(p, c, a) / facet_area;
        weight_c = area(p, a, b) / facet_area;
        feature_weight = weight_a * fi_a + weight_b * fi_b + weight_c * fi_c;
        feature_weights->push_back(feature_weight);
      }
    }

    // split
    Halfedge_handle split_long_edge(const Point &new_point,
                                    Halfedge_handle hh) {
      Halfedge_handle hnew = split_edge(hh);
      hnew->vertex()->point() = new_point;
      if (!hnew->is_border()) {
        split_facet(hnew, hnew->next()->next());
      }
      if (!hnew->opposite()->is_border()) {
        split_facet(hnew->opposite()->next(),
          hnew->opposite()->next()->next()->next());
      }
      return hnew;
    }

    // collapse
    Vertex_handle collapse_short_edge(const Point &new_point,
                                      Halfedge_handle hh) {
      bool joined = join_facets_before_collapse(hh);
      if (!joined) {    // should never go here
        std::cout << red << "unable to join facets (premature ending)"
          << white << std::endl;
        return NULL;
      }
      Halfedge_handle hh_joined = join_vertex(hh);
      hh_joined->vertex()->point() = new_point;
      return hh_joined->vertex();
    }

    bool is_collapsible(Halfedge_const_handle hh) const {
      // precondition: hh is not a boundary halfedge (its opposite can be)
      if (hh->is_border()) {
        return false;
      }
      Halfedge_const_handle ho = hh->opposite();
      if (ho->is_border()) {              // border case
        Vertex_const_handle vs = get_target_vertex(ho->next());
        Vertex_const_handle vt = get_source_vertex(ho->prev());
        if (vs == vt) {                   // do not close a degree-3 hole
          return false;
        }
      }
      else {                              // inner case
        Vertex_const_handle vr = get_opposite_vertex(hh);
        Vertex_const_handle vs = get_opposite_vertex(ho);
        if (vr == vs) {
          return false;
        }
        // inner edge but two boundary vertices
        if (is_on_boundary(get_target_vertex(hh)) &&
          is_on_boundary(get_target_vertex(ho))) {
          return false;
        }
      }
      return check_link_test(hh);         // link condition
    }

    bool collapse_would_cause_wrinkle(const Halfedge_list &halfedges,
        const Point &new_point, Halfedge_const_handle hh) const {
      for (auto it = halfedges.begin(); it != halfedges.end(); ++it) {
        Halfedge_handle h = *it;
        const Point &start_point = get_source_vertex(h)->point();
        const Point &end_point = get_target_vertex(h)->point();
        const Point &old_point = get_opposite_vertex(h)->point();
        FT radian = get_radian(end_point, old_point, start_point);
        if (radian < MIN_VALUE || CGAL_PI - radian < MIN_VALUE) {
          return true;    // degenerate cases
        }
        Plane plane(end_point, old_point, start_point);
        Point projection = plane.projection(new_point);
        if (same_side(start_point, end_point, old_point, projection) < 0) {
          return true;
        }
      }
      return false;
    }

    bool predict_facets_after_collapse(Halfedge_handle hh,
                                       Halfedge_list *halfedges) {
      // whether facets compose a ring, halfedges are inserted in order
      // we use the halfedges to represent the facets
      Vertex_handle vp = get_source_vertex(hh);
      Vertex_handle vq = get_target_vertex(hh);
      if (hh->opposite()->is_border()) {
        if (!hh->next()->opposite()->is_border()) {
          Halfedge_handle h_add = hh->opposite()->prev()->opposite()->next();
          while (h_add != hh->prev()) {
            halfedges->push_back(h_add);
            h_add = h_add->next()->opposite()->next();
          }
        }
        if (!hh->prev()->opposite()->is_border()) {
          Halfedge_handle h_add = hh->prev()->opposite()->next();
          while (h_add != hh->opposite()->next()->opposite()->prev()) {
            halfedges->push_back(h_add);
            h_add = h_add->next()->opposite()->next();
          }
          halfedges->push_back(h_add);
        }
        return false;
      }
      else if (is_on_boundary(vp)) {
        Halfedge_around_vertex_circulator vcirc = vp->vertex_begin();
        Halfedge_around_vertex_circulator vend = vcirc;
        CGAL_For_all(vcirc, vend) {
          if (vcirc->is_border()) {
            break;
          }
        }
        Halfedge_handle h_add = vcirc->opposite()->next();
        while (h_add != hh->opposite()->prev()) {
          halfedges->push_back(h_add);
          h_add = h_add->next()->opposite()->next();
        }
        h_add = h_add->opposite()->next();
        while (h_add != hh->prev()) {
          halfedges->push_back(h_add);
          h_add = h_add->next()->opposite()->next();
        }
        if (!h_add->opposite()->is_border()) {
          h_add = h_add->opposite()->next();
          while (!h_add->next()->opposite()->is_border()) {
            halfedges->push_back(h_add);
            h_add = h_add->next()->opposite()->next();
          }
          halfedges->push_back(h_add);
        }
        return false;
      }
      else if (is_on_boundary(vq)) {
        Halfedge_around_vertex_circulator vcirc = vq->vertex_begin();
        Halfedge_around_vertex_circulator vend = vcirc;
        CGAL_For_all(vcirc, vend) {
          if (vcirc->is_border()) {
            break;
          }
        }
        Halfedge_handle h_add = vcirc->opposite()->next();
        while (h_add != hh->prev()) {
          halfedges->push_back(h_add);
          h_add = h_add->next()->opposite()->next();
        }
        h_add = h_add->opposite()->next();
        while (h_add != hh->opposite()->prev()) {
          halfedges->push_back(h_add);
          h_add = h_add->next()->opposite()->next();
        }
        if (!h_add->opposite()->is_border()) {
          h_add = h_add->opposite()->next();
          while (!h_add->next()->opposite()->is_border()) {
            halfedges->push_back(h_add);
            h_add = h_add->next()->opposite()->next();
          }
          halfedges->push_back(h_add);
        }
        return false;
      }
      else {
        Halfedge_handle h_add = hh->prev()->opposite()->next();
        while (h_add != hh->opposite()->prev()) {
          halfedges->push_back(h_add);
          h_add = h_add->next()->opposite()->next();
        }
        h_add = h_add->opposite()->next();
        while (h_add != hh->prev()) {
          halfedges->push_back(h_add);
          h_add = h_add->next()->opposite()->next();
        }
        return true;
      }
    }

    // flip
    bool is_flippable(bool inherit_element_types,
        EdgeFlipStrategy edge_flip_strategy, Halfedge_const_handle hh) const {
      // step 1: check whether it is border edge or the new_edge already exist
      if (hh->is_border_edge()) {
        return false;
      }
      Vertex_const_handle p = get_target_vertex(hh);
      Vertex_const_handle q = get_source_vertex(hh);
      Vertex_const_handle r = get_opposite_vertex(hh);		// safe
      Vertex_const_handle s = get_opposite_vertex(hh->opposite());
      if (are_neighbors(r, s)) {
        return false;
      }
      // step 2: check whether we are flippling a crease edge
      if (inherit_element_types && is_crease_edge(hh)) {
        return false;
      }
      // step 3: check whether flipping improves the valence or angle
      if (edge_flip_strategy == EdgeFlipStrategy::k_improve_valence) {
        return edge_flip_would_improve_valence(hh);
      }
      else {
        return edge_flip_would_improve_radian(hh);
      }
    }

    // relocate
    bool relocate_would_cause_wrinkle(const Point &new_point,
                                      Vertex_const_handle vh) const {
      Halfedge_around_vertex_const_circulator vcirc = vh->vertex_begin();
      Halfedge_around_vertex_const_circulator vend = vcirc;
      CGAL_For_all(vcirc, vend) {
        if (!vcirc->is_border()) {
          const Point &start_point = get_opposite_vertex(vcirc)->point();
          const Point &end_point = get_source_vertex(vcirc)->point();
          const Point &old_point = get_target_vertex(vcirc)->point();
          FT radian = get_radian(end_point, old_point, start_point);
          if (radian < MIN_VALUE || CGAL_PI - radian < MIN_VALUE) {
            return true;  // degenerate cases
          }
          Plane plane(end_point, old_point, start_point);
          Point projection = plane.projection(new_point);
          if (same_side(start_point, end_point, old_point, projection) < 0) {
            return true;
          }
        }
      }
      return false;
    }

    // vertex_handle access
    inline Vertex_handle get_target_vertex(Halfedge_handle h) {
      return h->vertex();
    }

    inline Vertex_const_handle get_target_vertex(
        Halfedge_const_handle h) const {
      return h->vertex();
    }

    inline Vertex_handle get_source_vertex(Halfedge_handle h) {
      return h->opposite()->vertex();
    }

    inline Vertex_const_handle get_source_vertex(
        Halfedge_const_handle h) const {
      return h->opposite()->vertex();
    }

    inline Vertex_handle get_opposite_vertex(Halfedge_handle h) {
      assert(!h->is_border());
      if (h->is_border()) {
        std::cout << "Error definition of oppsite vertex for border halfedge."
          << std::endl;
        return NULL;
      }
      else {
        return h->next()->vertex();
      }
    }

    inline Vertex_const_handle get_opposite_vertex(
        Halfedge_const_handle h) const {
      assert(!h->is_border());
      if (h->is_border()) {
        std::cout << "Error definition of opposite vertex for border halfedge."
          << std::endl;
        return NULL;
      }
      else {
        return h->next()->vertex();
      }
    }

    // halfedge_handle access
    Halfedge_handle get_shortest_halfedge(Facet_handle fh) {
      Halfedge_handle hh = fh->halfedge();
      Halfedge_handle shortest_hh = hh;
      FT shortest_squared_length = squared_length(shortest_hh);
      FT length = squared_length(hh->next());   // check the next halfedge
      if (length < shortest_squared_length) {
        shortest_squared_length = length;
        shortest_hh = hh->next();
      }
      length = squared_length(hh->prev());      // check the previous halfedge
      if (length < shortest_squared_length) {
        shortest_hh = hh->prev();
      }
      return shortest_hh;
    }

    Halfedge_const_handle get_shortest_halfedge(Facet_const_handle fh) const {
      Halfedge_const_handle hh = fh->halfedge();
      Halfedge_const_handle shortest_hh = hh;
      FT shortest_squared_length = squared_length(shortest_hh);
      FT length = squared_length(hh->next());   // check the next halfedge
      if (length < shortest_squared_length) {
        shortest_squared_length = length;
        shortest_hh = hh->next();
      }
      length = squared_length(hh->prev());      // check the previous halfedge
      if (length < shortest_squared_length) {
        shortest_hh = hh->prev();
      }
      return shortest_hh;
    }

    Halfedge_handle get_longest_halfedge(Facet_handle fh) {
      Halfedge_handle hh = fh->halfedge();
      Halfedge_handle longest_hh = hh;
      FT longest_squared_length = squared_length(longest_hh);
      FT length = squared_length(hh->next());   // check the next halfedge
      if (length > longest_squared_length) {
        longest_squared_length = length;
        longest_hh = hh->next();
      }
      length = squared_length(hh->prev());      // check the previous halfedge
      if (length > longest_squared_length) {
        longest_squared_length = length;
        longest_hh = hh->prev();
      }
      return longest_hh;
    }

    Halfedge_const_handle get_longest_halfedge(Facet_const_handle fh) const {
      Halfedge_const_handle hh = fh->halfedge();
      Halfedge_const_handle longest_hh = hh;
      FT longest_squared_length = squared_length(longest_hh);
      FT length = squared_length(hh->next());   // check the next halfedge
      if (length > longest_squared_length) {
        longest_squared_length = length;
        longest_hh = hh->next();
      }
      length = squared_length(hh->prev());      // check the previous halfedge
      if (length > longest_squared_length) {
        longest_squared_length = length;
        longest_hh = hh->prev();
      }
      return longest_hh;
    }

    Halfedge_handle get_minimal_radian(FT *minimal_radian) {
      *minimal_radian = CGAL_PI;
      Halfedge_handle minimal_radian_halfedge = NULL;
      for (Facet_iterator fi = facets_begin(); fi != facets_end(); ++fi) {
        FT radian = get_smallest_radian(fi);
        if (radian < *minimal_radian) {
          *minimal_radian = radian;
          minimal_radian_halfedge = get_shortest_halfedge(fi);
        }
      }
      return minimal_radian_halfedge;
    }

    Halfedge_const_handle get_minimal_radian(FT *minimal_radian) const {
      *minimal_radian = CGAL_PI;
      Halfedge_const_handle minimal_radian_halfedge = NULL;
      for (Facet_const_iterator fi = facets_begin();
        fi != facets_end(); ++fi) {
        FT radian = get_smallest_radian(fi);
        if (radian < *minimal_radian) {
          *minimal_radian = radian;
          minimal_radian_halfedge = get_shortest_halfedge(fi);
        }
      }
      return minimal_radian_halfedge;
    }

    Halfedge_handle longest_side_propagation(Halfedge_handle hh) {
      // precondition: hh is the longest edge in its incident facet
      if (hh->is_border() || hh->opposite()->is_border()) {
        return hh;
      }
      Facet_handle fh = hh->opposite()->facet();
      Halfedge_handle longest_in_neighbor = get_longest_halfedge(fh);
      if (longest_in_neighbor == hh->opposite()) {
        return hh;
      }
      else {
        return longest_side_propagation(longest_in_neighbor);
      }
    }

    // vertex property access
    bool is_on_boundary(Vertex_const_handle vh) const {
      Halfedge_around_vertex_const_circulator vcirc = vh->vertex_begin();
      Halfedge_around_vertex_const_circulator vend = vcirc;
      CGAL_For_all(vcirc, vend) {
        if (vcirc->is_border()) {   // only one halfedge will answer true
          return true;
        }
      }
      return false;
    }

    FT get_vertex_capacity(Vertex_handle vh) const {
      FT vertex_capacity = 0.0;
      const Point &p = vh->point();
      Halfedge_around_vertex_circulator vcirc = vh->vertex_begin();
      Halfedge_around_vertex_circulator vend = vcirc;
      CGAL_For_all(vcirc, vend) {
        if (!vcirc->is_border()) {
          Facet_handle fh = vcirc->facet();
          const Point p1 = midpoint(vcirc);
          const Point c = centroid(fh);       // option: use the circumcenter?
          const Point p2 = midpoint(vcirc->next());
          vertex_capacity += area(p1, p, c);
          vertex_capacity += area(c, p, p2);
        }
      }
      return vertex_capacity;
    }

    // edge property access
    inline FT get_opposite_radian(Halfedge_const_handle hh) const {
      if (hh->is_border()) {
        std::cout << "Error definition of opposite angle for a border."
          << std::endl;
        return -1;
      }
      else {
        return get_radian(get_target_vertex(hh)->point(),
                          get_opposite_vertex(hh)->point(),
                          get_source_vertex(hh)->point());
      }
    }

    inline FT get_opposite_angle(Halfedge_const_handle hh) const {
      FT radian = get_opposite_radian(hh);
      return to_angle(radian);
    }

    inline FT squared_length(Halfedge_const_handle hh) const {
      const Point &a = get_source_vertex(hh)->point();
      const Point &b = get_target_vertex(hh)->point();
      return CGAL::squared_distance(a, b);
    }

    inline FT length(Halfedge_const_handle hh) const {
      return CGAL::sqrt(squared_length(hh));
    }

    inline Point midpoint(Halfedge_const_handle hh) const {
      const Point &a = get_target_vertex(hh)->point();
      const Point &b = get_source_vertex(hh)->point();
      return CGAL::midpoint(a, b);
    }

    inline bool is_crease_edge(Halfedge_const_handle hh) const {
      if (hh->normal_dihedral() == -1.0) {
        hh = hh->opposite();
      }
      return hh->is_crease();
    }

    Point get_least_qem_point(Halfedge_handle hh) const {
      std::set<Facet_handle> facets;
      collect_one_ring_facets_incident_to_edge(hh, &facets);
      const Point &start_point = get_source_vertex(hh)->point();
      const Point mid_point = midpoint(hh);
      const Point &end_point = get_target_vertex(hh)->point();
      FT qem_start = get_sum_qem_value(facets, start_point);
      FT qem_mid = get_sum_qem_value(facets, mid_point);
      FT qem_end = get_sum_qem_value(facets, end_point);
      if (qem_mid <= qem_start) {
        return qem_mid <= qem_end ? mid_point : end_point;
      }
      else {
        return qem_start <= qem_end ? start_point : end_point;
      }
    }

    // facet property access
    FT get_smallest_radian(Facet_const_handle fh) const {
      Halfedge_const_handle shortest_hh = get_shortest_halfedge(fh);
      return get_opposite_radian(shortest_hh);
    }

    FT get_smallest_angle(Facet_const_handle fh) const {
      FT radian = get_smallest_radian(fh);
      return to_angle(radian);
    }

    FT get_largest_radian(Facet_const_handle fh) const {
      Halfedge_const_handle longest_hh = get_longest_halfedge(fh);
      return get_opposite_radian(longest_hh);
    }

    FT get_largest_angle(Facet_const_handle fh) const {
      FT radian = get_largest_radian(fh);
      return to_angle(radian);
    }

    inline FT area(Facet_const_handle fh) const {
      // precondition: the normal has been calcualted
      FT smallest_radian = get_smallest_radian(fh);
      if (smallest_radian < MIN_VALUE) {
        return 0.0;
      }
      else {
        return CGAL::sqrt(triangle(fh).squared_area());
      }
    }

    FT squared_distance(const Point &p, Facet_handle fh,
      Point *nearest_point) const {
      Halfedge_handle hh = fh->halfedge();
      const Point &a = hh->vertex()->point();
      const Point &b = hh->next()->vertex()->point();
      const Point &c = hh->prev()->vertex()->point();
      Plane plane(a, b, c);
      Point projection = plane.projection(p);
      if (point_in_triangle(a, b, c, projection)) {
        *nearest_point = projection;
        return CGAL::squared_distance(p, projection);
      }
      else {
        Segment ab(a, b), bc(b, c), ca(c, a);
        FT sd_ab = CGAL::squared_distance(p, ab);
        FT sd_bc = CGAL::squared_distance(p, bc);
        FT sd_ca = CGAL::squared_distance(p, ca);
        Segment closest_segment;
        if (sd_ab < sd_bc) {
          closest_segment = sd_ab < sd_ca ? ab : ca;
        }
        else {
          closest_segment = sd_bc < sd_ca ? bc : ca;
        }
        *nearest_point = get_nearest_point(closest_segment, p);
        return CGAL::min(sd_ab, CGAL::min(sd_bc, sd_ca));
      }
    }

    // partial property access
    FT get_local_minimal_radian(const std::set<Facet_handle> &facets) const {
      FT minimal_radian = CGAL_PI;
      for (auto it = facets.begin(); it != facets.end(); ++it) {
        Facet_const_handle fh = *it;
        FT radian = get_smallest_radian(fh);
        minimal_radian = CGAL::min(minimal_radian, radian);
      }
      return minimal_radian;
    }

    FT get_minimal_radian_around_vertex(Vertex_const_handle vh) const {
      FT minimal_radian = CGAL_PI;
      Halfedge_around_vertex_const_circulator vcirc = vh->vertex_begin();
      Halfedge_around_vertex_const_circulator vend = vcirc;
      CGAL_For_all(vcirc, vend) {
        if (!vcirc->is_border()) {
          Facet_const_handle fh = vcirc->facet();
          FT radian = get_smallest_radian(fh);
          minimal_radian = CGAL::min(minimal_radian, radian);
        }
      }
      return minimal_radian;
    }

    FT get_minimal_radian_incident_to_vertex(Vertex_const_handle vh) const {
      FT minimal_radian = CGAL_PI;
      Halfedge_around_vertex_const_circulator vcirc = vh->vertex_begin();
      Halfedge_around_vertex_const_circulator vend = vcirc;
      CGAL_For_all(vcirc, vend) {
        if (!vcirc->is_border()) {
          FT radian = get_radian(get_source_vertex(vcirc)->point(),
            get_target_vertex(vcirc)->point(),
            get_opposite_vertex(vcirc)->point());
          minimal_radian = CGAL::min(minimal_radian, radian);
        }
      }
      return minimal_radian;
    }

    FT ge_min_squared_distance_in_one_ring_facets(Vertex_handle vh) const {
      FT min_sd = std::numeric_limits<double>::max();
      Halfedge_around_vertex_circulator vcirc = vh->vertex_begin();
      Halfedge_around_vertex_circulator vend = vcirc;
      CGAL_For_all(vcirc, vend) {
        if (!vcirc->is_border()) {
          Line line(get_source_vertex(vcirc)->point(),
            get_opposite_vertex(vcirc)->point());
          FT sd = CGAL::squared_distance(vh->point(), line);
          if (sd < min_sd) {
            min_sd = sd;
          }
        }
      }
      return min_sd;
    }

    void reset_facet_tags(int value, const Facet_list &facets) {
      for (auto it = facets.begin(); it != facets.end(); ++it) {
        Facet_handle fh = *it;
        fh->tag() = value;
      }
    }

    Point barycenter(bool inherit_element_types, FT feature_control_delta,
                     Vertex_const_handle vh) const {
      Vector vec = CGAL::NULL_VECTOR;
      Point pivot = vh->point();
      FT denominator = 0.0;
      Halfedge_const_list effective_edges;
      VertexType vertex_type = get_vertex_type(inherit_element_types,
        feature_control_delta, vh, &effective_edges);
      if (vertex_type == VertexType::k_feature_vertex) {
        return pivot;
      }
      else if (vertex_type == VertexType::k_crease_vertex) {
        assert(effective_edges.size() == 2);
        Halfedge_const_iter it;
        // step 1: get the initial_point
        for (it = effective_edges.begin(); it != effective_edges.end(); ++it) {
          Halfedge_const_handle hh = *it;
          const Point &neighbor = get_source_vertex(hh)->point();
          vec = vec + (neighbor - pivot);
        }
        Point initial_point = pivot + vec / 2;
        // step 2: project the initial_point
        Point min_projection(DOUBLE_MAX, DOUBLE_MAX, DOUBLE_MAX);
        FT min_sd = DOUBLE_MAX;
        for (it = effective_edges.begin(); it != effective_edges.end(); ++it) {
          Segment segment(get_source_vertex(*it)->point(),
                          get_target_vertex(*it)->point());
          Point nearest_point = get_nearest_point(segment, initial_point);
          FT sd = CGAL::squared_distance(initial_point, nearest_point);
          if (sd < min_sd) {
            min_sd = sd;
            min_projection = nearest_point;
          }
        }
        return min_projection;
      }
      else {
        Halfedge_around_vertex_const_circulator vcirc = vh->vertex_begin();
        Halfedge_around_vertex_const_circulator vend = vcirc;
        CGAL_For_all(vcirc, vend) {
          const Point &neighbor = get_source_vertex(vcirc)->point();
          vec = vec + (neighbor - pivot);
        }
        return pivot + vec / vh->degree();
      }
    }

    Point average_center(Vertex_const_handle vh) const {
      // step 1: check whether it is a boundary vertex
      Halfedge_const_list boundary_edges;
      Halfedge_const_list inner_edges;
      Halfedge_around_vertex_const_circulator vcirc = vh->vertex_begin();
      Halfedge_around_vertex_const_circulator vend = vcirc;
      CGAL_For_all(vcirc, vend) {
        if (vcirc->is_border() || vcirc->opposite()->is_border()) {
          boundary_edges.push_back(vcirc);
        }
        inner_edges.push_back(vcirc);
      }
      // step 2: get the barycenter according to different cases
      const Point &pivot = vh->point();
      Vector vec = CGAL::NULL_VECTOR;
      if (!boundary_edges.empty()) {      // boundary case
        for (Halfedge_const_iter it = boundary_edges.begin();
          it != boundary_edges.end(); ++it) {
          const Point &neighbor = get_source_vertex(*it)->point();
          vec = vec + (neighbor - pivot);
        }
        return pivot + vec / boundary_edges.size();
      }
      else {                              // inner case
        for (Halfedge_const_iter it = inner_edges.begin();
          it != inner_edges.end(); ++it) {
          const Point &neighbor = get_source_vertex(*it)->point();
          vec = vec + (neighbor - pivot);
        }
        return pivot + vec / inner_edges.size();
      }
    }

    Point cvt_barycenter(bool inherit_element_types, FT feature_control_delta,
                         Vertex_const_handle vh) const {
      Vector vec = CGAL::NULL_VECTOR;
      Point pivot = vh->point();
      FT denominator = 0.0;
      Halfedge_const_list effective_edges;
      VertexType vertex_type = get_vertex_type(inherit_element_types,
        feature_control_delta, vh, &effective_edges);
      if (vertex_type == VertexType::k_feature_vertex) {
        return pivot;
      }
      else if (vertex_type == VertexType::k_crease_vertex) {
        assert(effective_edges.size() == 2);
        Halfedge_const_iter it;
        // step 1: get the initial_point
        for (it = effective_edges.begin(); it != effective_edges.end(); ++it) {
          vec = vec + (midpoint(*it) - pivot);
        }
        Point initial_point = pivot + vec / 2;
        // step 2: project the initial_point
        Point min_projection(DOUBLE_MAX, DOUBLE_MAX, DOUBLE_MAX);
        FT min_sd = DOUBLE_MAX;
        for (it = effective_edges.begin(); it != effective_edges.end(); ++it) {
          Segment segment(get_source_vertex(*it)->point(),
            get_target_vertex(*it)->point());
          Point nearest_point = get_nearest_point(segment, initial_point);
          FT sd = CGAL::squared_distance(initial_point, nearest_point);
          if (sd < min_sd) {
            min_sd = sd;
            min_projection = nearest_point;
          }
        }
        return min_projection;
      }
      else {  //vertex_type == VertexType::k_smooth_vertex
        Halfedge_around_vertex_const_circulator vcirc = vh->vertex_begin();
        Halfedge_around_vertex_const_circulator vend = vcirc;
        Vector vector;
        CGAL_For_all(vcirc, vend) {
          Halfedge_const_handle h1 = vcirc;
          Halfedge_const_handle h2 = h1->next();
          Halfedge_const_handle h3 = h1->prev();
          const Point &ph1 = get_target_vertex(h1)->point();
          const Point &ph2 = get_target_vertex(h2)->point();
          const Point &ph3 = get_target_vertex(h3)->point();
          const Point p = CGAL::centroid(ph1, ph2, ph3);
          const Point p1 = midpoint(h1);
          const Point p2 = midpoint(h2);
          const Point c1 = CGAL::centroid(pivot, p, p1);
          const Point c2 = CGAL::centroid(pivot, p2, p);
          const FT area1 = area(pivot, p, p1);
          const FT area2 = area(pivot, p2, p);
          vector = area1 * (c1 - CGAL::ORIGIN) + area2 * (c2 - CGAL::ORIGIN);
          vec = vec + vector;
          denominator += (area1 + area2);
        }
        if (denominator < SQUARED_MIN_VALUE) {
          return pivot;
        }
        else {
          return CGAL::ORIGIN + vec / denominator;
        }
      }
    }

    // vertex_handle collection
    void collect_vertices(const std::set<Facet_handle> &facets,
                          std::set<Vertex_handle> *vertices) const {
      for (auto it = facets.begin(); it != facets.end(); ++it) {
        Facet_handle fh = *it;
        Halfedge_handle hh = fh->halfedge();
        vertices->insert(hh->vertex());
        vertices->insert(hh->next()->vertex());
        vertices->insert(hh->prev()->vertex());
      }
    }

    void collect_indicent_vertices(Vertex_handle vh,
                                   std::set<Vertex_handle> *vertices) {
      Halfedge_around_vertex_circulator vcirc = vh->vertex_begin();
      Halfedge_around_vertex_circulator vend = vcirc;
      CGAL_For_all(vcirc, vend) {
        vertices->insert(get_source_vertex(vcirc));
      }
    }

    // halfedge_handle collection
    VertexType get_vertex_type(bool inherit_element_types,
        FT feature_control_delta, Vertex_const_handle vh,
        Halfedge_const_list *effective_edges) const {
      if (inherit_element_types) {
        get_crease_edges_in_one_ring(vh, effective_edges);
        if (effective_edges->size() <= 1) {
          return VertexType::k_smooth_vertex;
        }
        else if (effective_edges->size() == 2) {
          return VertexType::k_crease_vertex;
        }
        else {
          return VertexType::k_feature_vertex;
        }
      }
      else {
        get_effective_edges_in_one_ring(feature_control_delta, vh,
          effective_edges);
        if (effective_edges->size() == 0) {
          return VertexType::k_feature_vertex;
        }
        else if (effective_edges->size() == 2) {
          return VertexType::k_crease_vertex;
        }
        else {
          return VertexType::k_smooth_vertex;
        }
      }
    }

    // facet_handle collection
    void extend_facets(const std::set<Facet_handle> &one_ring_facets,
        int stencil_ring_size, std::set<Facet_handle> *extended_facets) const {
      extended_facets->clear();
      extended_facets->insert(one_ring_facets.begin(), one_ring_facets.end());
      for (int i = 0; i < stencil_ring_size; ++i) {
        extend_facets_by_one_ring(extended_facets);
      }
    }

    void collect_one_ring_facets_incident_to_vertex(Vertex_handle vh,
        std::set<Facet_handle> *facets) const {
      Halfedge_around_vertex_circulator vcirc = vh->vertex_begin();
      Halfedge_around_vertex_circulator vend = vcirc;
      CGAL_For_all(vcirc, vend) {
        if (!vcirc->is_border()) {
          facets->insert(vcirc->facet());
        }
      }
    }

    void collect_one_ring_facets_incident_to_edge(Halfedge_handle hh,
        std::set<Facet_handle> *facets) const {
      Vertex_handle vp = hh->opposite()->vertex();
      Vertex_handle vq = hh->vertex();
      collect_one_ring_facets_incident_to_vertex(vp, facets);
      collect_one_ring_facets_incident_to_vertex(vq, facets);
    }

    void collect_facets_incident_to_edge(Halfedge_handle hh,
        std::set<Facet_handle> *facets) const {
      if (!hh->is_border()) {
        facets->insert(hh->facet());
      }
      if (!hh->opposite()->is_border()) {
        facets->insert(hh->opposite()->facet());
      }
    }

    // polyhedron property access
    Bbox get_bounding_box() const {
      Bbox bbox = Bbox(DOUBLE_MAX, DOUBLE_MAX, DOUBLE_MAX,
        DOUBLE_MIN, DOUBLE_MIN, DOUBLE_MIN);
      Point_const_iterator pit = points_begin();
      bbox = (*pit).bbox();
      for (; pit != points_end(); ++pit) {
        bbox = bbox + pit->bbox();
      }
      return bbox;
    }

    FT get_diagonal_length() const {
      Bbox bbox = get_bounding_box();
      if (bbox.xmax() <= bbox.xmin()) {
        return -1.0;
      }
      else {        // invalid case
        FT x_length = bbox.xmax() - bbox.xmin();
        FT y_length = bbox.ymax() - bbox.ymin();
        FT z_length = bbox.zmax() - bbox.zmin();
        return std::sqrt(x_length * x_length + y_length * y_length +
          z_length * z_length);
      }
    }

    FT get_average_length() const {
      if (size_of_halfedges() == 0) {
        return 0;
      }
      else {
        FT sum_length = 0.0;
        for (Edge_const_iterator ei = edges_begin(); ei != edges_end(); ++ei) {
          sum_length += length(ei);
        }
        return sum_length * 2.0 / size_of_halfedges();
      }
    }


    void trace_properties(std::string name) {
      std::cout << std::endl;
      std::transform(name.begin(), name.end(), name.begin(), [](unsigned char c){ return std::toupper(c);});
      std::cout << yellow << name + " PROPERTIES" << white << std::endl;
      std::cout << size_of_vertices() << " vertices" << std::endl;
      std::cout << size_of_facets() << " facets" << std::endl;
      std::cout << size_of_halfedges() / 2 << " edges" << std::endl;
      std::cout << nb_boundaries() << " boundary(ies)" << std::endl;
      std::cout << nb_components() << " component(s)" << std::endl;
      trace_edge_length();
    }


    void trace_additional_properties(FT diagonal_length) const {
      if (size_of_facets() == 0) {
        return;
      }
      FT min_facet_area = DOUBLE_MAX;
      FT max_squared_error = 0.0;
      FT min_radian = CGAL_PI;
      FT max_radian = 0.0;
      FT avg_min_radian = 0.0;
      size_t min_out_link_count = MAX_VALUE;
      for (Facet_const_iterator fi = facets_begin();
        fi != facets_end(); ++fi) {
        FT facet_area = area(fi);
        min_facet_area = CGAL::min(min_facet_area, facet_area);
        min_out_link_count =
          CGAL::min(min_out_link_count, fi->facet_out_links().size());
        max_squared_error =
          CGAL::max(max_squared_error, fi->max_squared_error());
        FT smallest_radian = get_smallest_radian(fi);
        min_radian = CGAL::min(min_radian, smallest_radian);
        avg_min_radian += smallest_radian;
        FT largest_radian = get_largest_radian(fi);
        max_radian = CGAL::max(max_radian, largest_radian);
      }
      avg_min_radian /= size_of_facets();
      FT rms_error = get_rms_distance();
      FT max_error = CGAL::sqrt(max_squared_error);
      //FT diagonal_length = get_diagonal_length();
      std::cout << "Minimal quality: " << get_min_quality() << std::endl;
      std::cout << "Average quality: " << get_avg_quality() << std::endl;
      std::cout << "Minimal angle: " << to_angle(min_radian) << std::endl;
      std::cout << "Average minimal angle: "
        << to_angle(avg_min_radian) << std::endl;
      std::cout << "Maximal angle: " << to_angle(max_radian) << std::endl;
      std::cout << "Maximal error: " << max_error << std::endl;
      std::cout << "Maximal error ratio: "
        << max_error / diagonal_length * 100 << "%" << std::endl;
      std::cout << "RMS facet error " << rms_error << std::endl;
      std::cout << "RMS facet error ratio: "
        << rms_error / diagonal_length * 100 << "%" << std::endl;
      std::cout << "Angles smaller than 30â°: "
        << get_smaller_angle_ratio(30.0) * 100 << "%" << std::endl;
      std::cout << "Regular vertices ratio: "
        << get_regular_vertex_ratio() * 100 << "%" << std::endl;
      std::cout << "Minimal facet area: " << min_facet_area << std::endl;
      std::cout << "Minimal facet out link count: "
        << min_out_link_count << std::endl;
    }

    int get_facet_out_link_count() const {
      size_t facet_out_link_count = 0;
      for (Facet_const_iterator fi = facets_begin();
        fi != facets_end(); ++fi) {
        facet_out_link_count += fi->facet_out_links().size();
      }
      return static_cast<int>(facet_out_link_count);
    }

    int get_edge_out_link_count() const {
      size_t edge_out_link_count = 0;
      for (Edge_const_iterator ei = edges_begin(); ei != edges_end(); ++ei) {
        Halfedge_const_handle hh = ei;
        if (hh->normal_dihedral() == -1.0) {
          hh = hh->opposite();
        }
        edge_out_link_count += hh->edge_out_links().size();
      }
      return static_cast<int>(edge_out_link_count);
    }

    int get_vertex_out_link_count() const {
      return static_cast<int>(size_of_vertices());
    }

    void normalize(FT radius) {
      // step 1: calculate the bounding box
      Point_iterator pi = points_begin();
      Bbox bb = pi->bbox();
      for (; pi != points_end(); ++pi) {
        bb = bb + pi->bbox();
      }
      // step 2: get the center and radius
      FT x_center = (bb.xmin() + bb.xmax()) / 2.0;
      FT y_center = (bb.ymin() + bb.ymax()) / 2.0;
      FT z_center = (bb.zmin() + bb.zmax()) / 2.0;
      FT x_radius = (bb.xmax() - bb.xmin()) / 2.0;
      FT y_radius = (bb.ymax() - bb.ymin()) / 2.0;
      FT z_radius = (bb.zmax() - bb.zmin()) / 2.0;
      FT max_radius = max(max(x_radius, y_radius), z_radius);
      // step 3: transfer
      for (pi = points_begin(); pi != points_end(); ++pi) {
        (*pi) = Point(pi->x() - x_center,
          pi->y() - y_center,
          pi->z() - z_center);
      }
      // step 4: scale
      for (pi = points_begin(); pi != points_end(); ++pi) {
        (*pi) = Point(pi->x() * radius / max_radius,
          pi->y() * radius / max_radius,
          pi->z() * radius / max_radius);
      }
    }

    // IO
    void save_as(std::string file_name) const {
      size_t pos = file_name.find_last_of('.');
      if (pos == std::string::npos) {
        std::cout << "Invalid file name." << std::endl;
        return;
      }
      std::string extension = file_name.substr(pos);
      std::transform(extension.begin(), extension.end(),
                     extension.begin(), [](unsigned char c){ return std::tolower(c);});
      if (extension == ".off") {
        save_as_off(file_name);
      }
      //else if (extension == ".mesh") {
      //  save_as_mesh(file_name);    // save as .mesh for YANS test
      //}
      else {
        std::cout << "Invalid file name." << std::endl;
      }
    }

    void save_as_off(std::string file_name) const {
      std::ofstream ofs(file_name);
      CGAL::set_ascii_mode(ofs);
      ofs << "OFF\n" << size_of_vertices() << ' '
        << size_of_facets() << " 0\n";
      std::copy(points_begin(), points_end(),
        std::ostream_iterator<Point>(ofs, "\n"));
      for (Facet_const_iterator fi = facets_begin();
        fi != facets_end(); ++fi) {
        Halfedge_around_facet_const_circulator j = fi->facet_begin();
        CGAL_assertion(CGAL::circulator_size(j) >= 3);
        ofs << CGAL::circulator_size(j) << ' ';
        do {
          ofs << ' ' << std::distance(vertices_begin(), get_target_vertex(j));
        } while (++j != fi->facet_begin());
        ofs << "\n";
      }
      ofs.flush();
      ofs.close();
    }

    // utilities
    static inline FT to_radian(FT angle) {
      return angle * CGAL_PI / 180.0;
    }

    static inline FT to_angle(FT radian) {
      return radian * 180.0 / CGAL_PI;
    }

    static inline FT area(const Point &a, const Point &b, const Point &c) {
      if (CGAL::squared_distance(a, b) < SQUARED_MIN_VALUE ||
        CGAL::squared_distance(a, c) < SQUARED_MIN_VALUE ||
        CGAL::squared_distance(b, c) < SQUARED_MIN_VALUE) {
        return 0.0;
      }
      Triangle t(a, b, c);
      return CGAL::sqrt(t.squared_area());
    }

    static inline FT get_smallest_radian(const Point &a, const Point &b,
                                         const Point &c) {
      FT ab = CGAL::squared_distance(a, b);
      FT ac = CGAL::squared_distance(a, c);
      FT bc = CGAL::squared_distance(b, c);
      if (ab < ac) {
        return ab < bc ? get_radian(b, c, a) : get_radian(c, a, b);
      }
      else {
        return ac < bc ? get_radian(a, b, c) : get_radian(c, a, b);
      }
    }

    static FT get_radian(const Point &p1, const Point &p2, const Point &p3) {
      if (CGAL::squared_distance(p1, p2) < SQUARED_MIN_VALUE ||
        CGAL::squared_distance(p1, p3) < SQUARED_MIN_VALUE ||
        CGAL::squared_distance(p2, p3) < SQUARED_MIN_VALUE) {
        return 0.0;                         // degenerated case
      }
      Vector v1 = p1 - p2;				// v1
      FT v1_length = std::sqrt(v1 * v1);
      if (v1_length < MIN_VALUE) {
        return 0.0;
      }
      v1 = v1 / v1_length;
      Vector v2 = p3 - p2;				// v2
      FT v2_length = std::sqrt(v2 * v2);
      if (v2_length < MIN_VALUE) {
        return 0.0;
      }
      v2 = v2 / v2_length;
      FT cos_value = v1 * v2;
      if (cos_value > 1.0) {
        cos_value = 1.0;
      }
      if (cos_value < -1.0) {
        cos_value = -1.0;
      }
      return std::acos(cos_value);
    }

    static FT get_angle(const Point &p1, const Point &p2, const Point &p3) {
      FT radian = get_radian(p1, p2, p3);
      return to_angle(radian);
    }

    static FT get_radian(const Vector &v1, const Vector &v2) {
      if (v1 == CGAL::NULL_VECTOR || v2 == CGAL::NULL_VECTOR) {
        return 0.0;     // degenerated case
      }
      FT v1_length = std::sqrt(v1.squared_length());
      FT v2_length = std::sqrt(v2.squared_length());
      if (v1_length < MIN_VALUE || v2_length < MIN_VALUE) {
        return 0.0;
      }
      FT inner_product = (v1 * v2) / (v1_length * v2_length);
      if (inner_product > 1.0) {
        inner_product = 1.0;
      }
      if (inner_product < -1.0) {
        inner_product = -1.0;
      }
      return std::acos(inner_product);
    }

    static FT get_angle(const Vector &v1, const Vector &v2) {
      FT radian = get_radian(v1, v2);
      return to_angle(radian);
    }

    static int same_side(const Point &a, const Point &b, const Point &c,
                         const Point &p) {
      // determine whether c and p lays on the same side of ab
      // 1: the same side;
      // 0: c or p is on the line ab or projected on the line;
      // -1: the opposite side

      // precondition: a, b, c and p are on the same plane
      Vector ab = b - a;
      Vector ac = c - a;
      Vector ap = p - a;
      Vector v1 = CGAL::cross_product(ab, ac);
      Vector v2 = CGAL::cross_product(ab, ap);
      FT cp = v1 * v2;
      if (cp > 0) {
        return 1;
      }
      else if (cp == 0.0) {
        return 0;
      }
      else {
        return -1;
      }
    }

    static Point get_nearest_point(const Segment &s, const Point &p) {
      Line line = s.supporting_line();
      if (CGAL::squared_distance(s, p) == CGAL::squared_distance(line, p)) {
        return line.projection(p);
      }
      else {
        return CGAL::squared_distance(s.source(), p) <
          CGAL::squared_distance(s.target(), p) ? s.source() : s.target();
      }
    }

    static bool point_in_triangle(const Point &a, const Point &b, const Point &c,
      const Point &p) {
      return same_side(a, b, c, p) > 0 &&
        same_side(b, c, a, p) > 0 &&
        same_side(c, a, b, p) > 0;
    }

    // visual elements
    void compute_facets(DrawType draw_type, RenderType render_type,
        bool inherit_element_types, FT feature_control_delta,
        FT sum_theta_value, FT dihedral_theta_value,
        FT max_error_threshold_value, bool is_input,
        std::vector<float> *pos_faces, std::vector<float> *pos_face_normals,
        std::vector<float> *pos_face_colors,
        std::vector<float> *pos_boundaries,
        std::vector<float> *pos_samples) const {
      pos_faces->resize(0);
      pos_face_normals->resize(0);
      pos_face_colors->resize(0);
      pos_boundaries->resize(0);
      pos_samples->resize(0);
      switch (draw_type) {
      case DrawType::k_polyhedron:
        compute_polyhedron_facets(render_type, sum_theta_value,
          dihedral_theta_value, max_error_threshold_value, is_input,
          pos_faces, pos_face_normals, pos_face_colors);
        break;
      case DrawType::k_all_voronoi:
        compute_all_voronois(render_type, sum_theta_value,
          dihedral_theta_value, pos_faces, pos_face_normals,
          pos_face_colors, pos_boundaries, pos_samples);
        break;
      case DrawType::k_vertex_voronoi:
        compute_vertex_voronois(render_type, inherit_element_types,
          feature_control_delta, sum_theta_value, dihedral_theta_value,
          pos_faces, pos_face_normals, pos_face_colors,
          pos_boundaries, pos_samples);
        break;
      case DrawType::k_edge_voronoi:
        compute_edge_voronois(render_type, sum_theta_value,
          dihedral_theta_value, pos_faces, pos_face_normals,
          pos_face_colors, pos_boundaries, pos_samples);
        break;
      case DrawType::k_facet_voronoi:
        compute_facet_voronois(render_type, sum_theta_value,
          dihedral_theta_value, pos_faces, pos_face_normals,
          pos_face_colors, pos_boundaries, pos_samples);
        break;
      default:
        break;
      }
    }

    void compute_facet_samples(std::vector<float> *pos_samples) const {
      for (Facet_const_iterator fi = facets_begin();
        fi != facets_end(); ++fi) {
        for (Link_list_const_iter it = fi->facet_out_links().begin();
          it != fi->facet_out_links().end(); ++it) {
          const Point &p = it->second.first;
          compute_point(p, pos_samples);
        }
      }
    }

    void compute_edges(std::vector<float> *pos_edges) const {
      pos_edges->resize(0);
      for (Edge_const_iterator ei = edges_begin(); ei != edges_end(); ++ei) {
        compute_halfedge(ei, pos_edges);
      }
    }

    void compute_min_radian_edges(
        std::vector<float> *pos_min_radian_edges) const {
      pos_min_radian_edges->resize(0);
      FT min_radian = CGAL_PI;
      Halfedge_const_handle min_radian_hh = get_minimal_radian(&min_radian);
      if (min_radian != CGAL_PI) {
        compute_halfedge(min_radian_hh->next(), pos_min_radian_edges);
        compute_halfedge(min_radian_hh->prev(), pos_min_radian_edges);
      }
    }

    void compute_classified_edges(std::vector<float> *pos_normal_edges,
        std::vector<float> *pos_special_edges) const {
      pos_normal_edges->resize(0);
      pos_special_edges->resize(0);
      for (Edge_const_iterator ei = edges_begin(); ei != edges_end(); ++ei) {
        Halfedge_const_handle hh = ei;
        if (hh->normal_dihedral() == -1.0) {
          hh = hh->opposite();
        }
        if (hh->is_crease()) {
          compute_halfedge(hh, pos_special_edges);
        }
        else {
          compute_halfedge(hh, pos_normal_edges);
        }
      }
    }

    void compute_facet_start_points(
        std::vector<float> *pos_facet_start_point) const {
      pos_facet_start_point->clear();
      for (Facet_const_iterator fi = facets_begin();
        fi != facets_end(); ++fi) {
        for (Link_list_const_iter it = fi->facet_out_links().begin();
          it != fi->facet_out_links().end(); ++it) {
          compute_point(it->second.first, pos_facet_start_point);
        }
      }
    }

    void compute_facet_end_points(
        std::vector<float> *pos_facet_end_point) const {
      pos_facet_end_point->clear();
      for (Facet_const_iterator fi = facets_begin();
        fi != facets_end(); ++fi) {
        for (Link_list_const_iter it = fi->facet_out_links().begin();
          it != fi->facet_out_links().end(); ++it) {
          compute_point(it->second.second, pos_facet_end_point);
        }
      }
    }

    void compute_facet_links(
        std::vector<float> *pos_facet_links) const {
      pos_facet_links->clear();
      for (Facet_const_iterator fi = facets_begin();
        fi != facets_end(); ++fi) {
        for (Link_list_const_iter it = fi->facet_out_links().begin();
          it != fi->facet_out_links().end(); ++it) {
          const Point &p = it->second.first;
          const Point &q = it->second.second;
          compute_segment(p, q, pos_facet_links);
        }
      }
    }

    void compute_edge_start_points(
        std::vector<float> *pos_edge_start_points) const {
      pos_edge_start_points->clear();
      for (Edge_const_iterator ei = edges_begin(); ei != edges_end(); ++ei) {
        Halfedge_const_handle hh = ei;
        if (hh->normal_dihedral() == -1.0) {
          hh = hh->opposite();
        }
        for (auto it = hh->edge_out_links().begin();
          it != hh->edge_out_links().end(); ++it) {
          compute_point(it->second.first, pos_edge_start_points);
        }
      }
    }

    void compute_edge_end_points(
        std::vector<float> *pos_edge_end_points) const {
      pos_edge_end_points->clear();
      for (Edge_const_iterator ei = edges_begin(); ei != edges_end(); ++ei) {
        Halfedge_const_handle hh = ei;
        if (hh->normal_dihedral() == -1.0) {
          hh = hh->opposite();
        }
        for (auto it = hh->edge_out_links().begin();
          it != hh->edge_out_links().end(); ++it) {
          compute_point(it->second.second, pos_edge_end_points);
        }
      }
    }

    void compute_edge_links(
        std::vector<float> *pos_edge_links) const {
      pos_edge_links->clear();
      for (Edge_const_iterator ei = edges_begin(); ei != edges_end(); ++ei) {
        Halfedge_const_handle hh = ei;
        if (hh->normal_dihedral() == -1.0) {
          hh = hh->opposite();
        }
        for (auto it = hh->edge_out_links().begin();
          it != hh->edge_out_links().end(); ++it) {
          const Point &p = it->second.first;
          const Point &q = it->second.second;
          compute_segment(p, q, pos_edge_links);
        }
      }
    }

    void compute_vertex_start_points(
        std::vector<float> *pos_vertex_start_points) const {
      pos_vertex_start_points->clear();
      for (Vertex_const_iterator vi = vertices_begin();
        vi != vertices_end(); ++vi) {
        compute_point(vi->vertex_out_link().second.first,
          pos_vertex_start_points);
      }
    }

    void compute_vertex_end_points(
        std::vector<float> *pos_vertex_end_points) const {
      pos_vertex_end_points->clear();
      for (Vertex_const_iterator vi = vertices_begin();
        vi != vertices_end(); ++vi) {
        compute_point(vi->vertex_out_link().second.second,
          pos_vertex_end_points);
      }
    }

    void compute_vertex_links(
        std::vector<float> *pos_vertex_links) const {
      pos_vertex_links->clear();
      for (Vertex_const_iterator vi = vertices_begin();
        vi != vertices_end(); ++vi) {
        const Point &p = vi->vertex_out_link().second.first;
        const Point &q = vi->vertex_out_link().second.second;
        compute_segment(p, q, pos_vertex_links);
      }
    }

   private:
    // normals
    void calculate_facet_normal(Facet_handle fh) {
      Normal sum = CGAL::NULL_VECTOR;
      Halfedge_around_facet_circulator h = fh->facet_begin();
      do {
        Normal normal = CGAL::cross_product(
          get_opposite_vertex(h)->point() - get_target_vertex(h)->point(),
          get_source_vertex(h)->point() - get_opposite_vertex(h)->point());
        FT sqnorm = normal * normal;
        if (sqnorm != 0.0) {
          normal = normal / (double)std::sqrt(sqnorm);
        }
        sum = sum + normal;
      } while (++h != fh->facet_begin());
      FT sqnorm = sum * sum;
      if (sqnorm != 0.0) {
        fh->normal() = sum / std::sqrt(sqnorm);
      }
      else
      {
        std::cerr << red << "degenerated facet" << white << std::endl;
        fh->normal() = CGAL::NULL_VECTOR;
      }
    }

    void calculate_vertex_normal(Vertex_handle vh) {
      // precondition: the facet normals are computed
      Normal normal = CGAL::NULL_VECTOR;
      Halfedge_around_vertex_const_circulator he = vh->vertex_begin();
      Halfedge_around_vertex_const_circulator begin = he;
      CGAL_For_all(he, begin) {
        if (!he->is_border()) {
          normal = normal + he->facet()->normal();
        }
      }
      FT sqnorm = normal * normal;
      if (sqnorm != 0.0) {
        vh->normal() = normal / std::sqrt(sqnorm);
      }
      else {
        vh->normal() = CGAL::NULL_VECTOR;
      }
    }

    // links
    void clear_local_in_links(std::set<Facet_handle> *in_link_facets) const {
      for (typename std::set<Facet_handle>::iterator it = in_link_facets->begin();
        it != in_link_facets->end(); ++it) {
        Facet_handle fh = *it;
        fh->facet_in_links().clear();
        fh->edge_in_links().clear();
        fh->vertex_in_links().clear();
      }
    }

    // max errors
    void reset_facet_max_squared_errors(FT value) {
      for (Facet_iterator fi = facets_begin(); fi != facets_end(); ++fi) {
        fi->max_squared_error() = value;
      }
    }

    void update_facet_in_max_squared_errors() {
      FT se = 0.0;  // squared error
      for (Facet_iterator fi = facets_begin(); fi != facets_end(); ++fi) {
        for (Link_iter_list_iter it = fi->facet_in_links().begin();
          it != fi->facet_in_links().end(); ++it) {
          Link_list_iter llit = *it;
          const Link &link = *llit;
          se = CGAL::squared_distance(link.second.first, link.second.second);
          fi->max_squared_error() = CGAL::max(fi->max_squared_error(), se);
        }
      }
    }

    void update_facet_out_max_squared_errors() {
      FT se = 0.0;
      for (Facet_iterator fi = facets_begin(); fi != facets_end(); ++fi) {
        for (Link_list_iter it = fi->facet_out_links().begin();
          it != fi->facet_out_links().end(); ++it) {
          const Link &link = *it;
          se = CGAL::squared_distance(link.second.first, link.second.second);
          fi->max_squared_error() = CGAL::max(fi->max_squared_error(), se);
        }
      }
    }

    void update_edge_in_max_squared_errors() {
      FT se = 0.0;
      for (Facet_iterator fi = facets_begin(); fi != facets_end(); ++fi) {
        for (Link_iter_list_iter it = fi->edge_in_links().begin();
          it != fi->edge_in_links().end(); ++it) {
          Link_list_iter llit = *it;
          const Link &link = *llit;
          se = CGAL::squared_distance(link.second.first, link.second.second);
          fi->max_squared_error() = CGAL::max(fi->max_squared_error(), se);
        }
      }
    }

    void update_edge_out_max_squared_errors() {
      FT se = 0.0;
      for (Edge_iterator ei = edges_begin(); ei != edges_end(); ++ei) {
        Halfedge_handle hh = ei;
        if (hh->normal_dihedral() == -1.0) {
          hh = hh->opposite();
        }
        FT max_se = 0.0;
        for (Link_list_iter it = hh->edge_out_links().begin();
          it != hh->edge_out_links().end(); ++it) {
          const Link &link = *it;
          se = CGAL::squared_distance(link.second.first, link.second.second);
          max_se = CGAL::max(max_se, se);
        }
        Facet_handle fh = hh->facet();
        fh->max_squared_error() = CGAL::max(fh->max_squared_error(), max_se);
        if (!hh->opposite()->is_border()) {
          fh = hh->opposite()->facet();
          fh->max_squared_error() = CGAL::max(fh->max_squared_error(), max_se);
        }
      }
    }

    void update_vertex_in_max_squared_errors() {
      FT se = 0.0;
      for (Facet_iterator fi = facets_begin(); fi != facets_end(); ++fi) {
        for (auto it = fi->vertex_in_links().begin();
          it != fi->vertex_in_links().end(); ++it) {
          Link *link = *it;
          se = CGAL::squared_distance(link->second.first, link->second.second);
          fi->max_squared_error() = CGAL::max(fi->max_squared_error(), se);
        }
      }
    }

    void update_vertex_out_max_squared_errors() {
      for (Vertex_iterator vi = vertices_begin(); vi != vertices_end(); ++vi) {
        const Link &link = vi->vertex_out_link();
        FT se = CGAL::squared_distance(link.second.first, link.second.second);
        Halfedge_around_vertex_circulator vcirc = vi->vertex_begin();
        Halfedge_around_vertex_circulator vend = vcirc;
        CGAL_For_all(vcirc, vend) {
          if (!vcirc->is_border()) {
            Facet_handle fh = vcirc->facet();
            fh->max_squared_error() = CGAL::max(fh->max_squared_error(), se);
          }
        }
      }
    }

    // feature intensities
    void calculate_vertex_feature_intensities(FT sum_theta, FT sum_delta) {
      // precondition: the halfedge's normal_dihedrals have been updated
      for (Vertex_iterator vi = vertices_begin(); vi != vertices_end(); ++vi) {
        calculate_vertex_feature_intensity(vi, sum_theta, sum_delta);
      }
    }

    void calculate_vertex_feature_intensity(Vertex_handle vh,
                                            FT sum_theta, FT sum_delta) {
      bool on_boundary = false;
      FT sum_radian = 0.0;
      FT gaussian_curvature = 0.0;
      FT max_halfedge_dihedral = 0.0;
      Halfedge_around_vertex_circulator vcirc = vh->vertex_begin();
      Halfedge_around_vertex_circulator vend = vcirc;
      Halfedge_const_handle hh;
      CGAL_For_all(vcirc, vend) {
        if (vcirc->is_border()) {
          on_boundary = true;
          hh = vcirc->opposite();
          max_halfedge_dihedral = CGAL::max(max_halfedge_dihedral,
            hh->normal_dihedral());
        }
        else {
          hh = vcirc;
          sum_radian += get_radian(get_source_vertex(hh)->point(),
            get_target_vertex(hh)->point(),
            get_target_vertex(hh->next())->point());
          if (hh->normal_dihedral() == -1.0) {
            hh = hh->opposite();
          }
          max_halfedge_dihedral = CGAL::max(max_halfedge_dihedral,
            hh->normal_dihedral());
        }
      }
      if (on_boundary) {
        gaussian_curvature = CGAL::abs(CGAL_PI - sum_radian);
      }
      else {
        gaussian_curvature = CGAL::abs(2 * CGAL_PI - sum_radian);
      }
      gaussian_curvature /= sum_delta;
      gaussian_curvature = CGAL::min(gaussian_curvature, sum_theta * CGAL_PI);
      vh->gaussian_curvature() = gaussian_curvature;
      vh->max_halfedge_dihedral() = max_halfedge_dihedral;
    }

    void calculate_edge_feature_intensities(FT dihedral_theta,
                                            FT dihedral_delta) {
      for (Edge_iterator ei = edges_begin(); ei != edges_end(); ++ei) {
        calculate_edge_feature_intensity(ei, false,
          dihedral_theta, dihedral_delta);
      }
    }

    void calculate_edge_feature_intensity(Edge_iterator ei,
        bool reset_normal_dihedral, FT dihedral_theta, FT dihedral_delta) {
      // 1) normalized diheral = min(diheral_theta, dihedral / dihedral_delta);
      // 2) default value: dihedral_theta = 1 (PI), dihedral_delta = 0.5
      // 3) the range of the normalized dihedral: [0, dihedral_theta * PI]
      Halfedge_handle hh = ei;
      if (hh->is_border()) {
        hh = hh->opposite();
      }
      if (hh->opposite()->is_border()) {      // boundary case
        hh->normal_dihedral() = dihedral_theta * CGAL_PI;
        hh->opposite()->normal_dihedral() = -1.0;
        hh->opposite()->is_crease() = false;
      }
      else {                                  // inner case
        FT normal_dihedral = get_normal_dihedral(hh);
        normal_dihedral /= dihedral_delta;    // normalize the normal dihedral
        normal_dihedral = CGAL::min(normal_dihedral, dihedral_theta * CGAL_PI);
        if (reset_normal_dihedral) {
          if (hh->normal_dihedral() == -1.0) {
            hh->is_crease() = hh->opposite()->is_crease();  // reset is_crease
          }
          hh->normal_dihedral() = normal_dihedral;
        }
        else {
          if (hh->opposite()->normal_dihedral() != -1.0) {  // no need to reset
            hh = hh->opposite();
          }
          hh->normal_dihedral() = normal_dihedral;
        }
        hh->opposite()->normal_dihedral() = -1.0;
        hh->opposite()->is_crease() = false;
      }
    }

    void calculate_edge_classifications(FT feature_control_delta) {
      // precondition: edge and normal feature intensities has been computed
      for (Vertex_iterator vi = vertices_begin(); vi != vertices_end(); ++vi) {
        Halfedge_list effective_edges;
        get_effective_edges_in_one_ring(feature_control_delta, vi,
          &effective_edges);
        if (effective_edges.size() == 2) {
          for (Halfedge_iter hi = effective_edges.begin();
            hi != effective_edges.end(); ++hi) {
            Halfedge_handle hh = *hi;
            if (hh->normal_dihedral() == -1.0) {
              hh = hh->opposite();
            }
            hh->is_crease() = true;
          }
        }
      }
    }

    // elements properties
    bool is_regular(Vertex_const_handle vh) const {
      size_t degree = vh->vertex_degree();
      if (is_on_boundary(vh)) {
        return degree >= 3 && degree <= 5;
        //return degree == 4;
      }
      else {
        return degree >= 5 && degree <= 7;
        //return degree == 6;
      }
    }

    bool are_neighbors(Vertex_const_handle va, Vertex_const_handle vb) const {
      Halfedge_around_vertex_const_circulator vcirc = va->vertex_begin();
      Halfedge_around_vertex_const_circulator vend = vcirc;
      CGAL_For_all(vcirc, vend) {
        if (get_source_vertex(vcirc) == vb) {
          return true;
        }
      }
      return false;
    }

    FT get_normal_dihedral(Halfedge_handle hh) const {
      // get the dihedral of the normals between the two incident facets
      Facet_handle fh1 = hh->facet();
      Facet_handle fh2 = hh->opposite()->facet();
      FT cos_value = fh1->normal() * fh2->normal();
      if (cos_value > 1.0) {
        cos_value = 1.0;
      }
      else if (cos_value < -1.0) {
        cos_value = -1.0;
      }
      return std::acos(cos_value);       // expressed in radians [0, PI]
    }

    FT get_quality(Facet_const_handle fh) const {
      Halfedge_const_handle hh = fh->halfedge();
      FT s_t = area(fh);						              // area of the triangle
      FT h_t = length(get_longest_halfedge(fh));  // longest edge
      FT a = length(hh);							            // length of the first edge
      FT b = length(hh->next());					        // length of the second edge
      FT c = length(hh->prev());					        // length of the third edge
      FT p_t = (a + b + c) / 2;				            // halfedge perimeter
      if (h_t <= 0) {                             // invalid case
        return 0.0;
      }
      else {
        return 6.0 * s_t / (CGAL::sqrt(3.0) * p_t * h_t);
      }
    }

    Point centroid(Facet_const_handle fh) const {
      Halfedge_const_handle hh = fh->halfedge();
      const Point &a = hh->vertex()->point();
      const Point &b = hh->next()->vertex()->point();
      const Point &c = hh->prev()->vertex()->point();
      return CGAL::centroid(a, b, c);
    }

    Triangle triangle(Facet_const_handle fh) const {
      Halfedge_const_handle hh = fh->halfedge();
      const Point &a = hh->vertex()->point();
      const Point &b = hh->next()->vertex()->point();
      const Point &c = hh->prev()->vertex()->point();
      return Triangle(a, b, c);
    }

    FT get_sum_area(const Facet_list &facets) const {
      FT total_area = 0.0;
      for (auto it = facets.begin(); it != facets.end(); ++it) {
        total_area += area(*it);
      }
      return total_area;
    }

    FT get_sum_qem_value(const std::set<Facet_handle> &facets,
                         const Point &point) const {
      FT sum_qem_value = 0.0;
      for (auto it = facets.begin(); it != facets.end(); ++it) {
        Facet_handle fh = *it;
        if (area(fh) >= SQUARED_MIN_VALUE) {  // ignore too small facets
          Halfedge_handle hh = fh->halfedge();
          Plane plane(get_source_vertex(hh)->point(),
            get_target_vertex(hh)->point(),
            get_opposite_vertex(hh)->point());
          sum_qem_value += CGAL::squared_distance(point, plane);
        }
      }
      return sum_qem_value;
    }

    // polyhedron properties
    FT get_avg_quality() const {
      FT sum_quality = 0.0;
      for (Facet_const_iterator fi = facets_begin();
        fi != facets_end(); ++fi) {
        sum_quality += get_quality(fi);
      }
      if (size_of_facets() == 0) {    // invalid case
        return -1.0;
      }
      else {
        return sum_quality / size_of_facets();
      }
    }

    FT get_min_quality() const {
      FT min_quality = DOUBLE_MAX, quality;
      for (Facet_const_iterator fi = facets_begin();
        fi != facets_end(); ++fi) {
        quality = get_quality(fi);
        min_quality = CGAL::min(min_quality, quality);
      }
      return min_quality;
    }

    FT get_rms_distance() const {
      FT rms_distance = 0.0;
      size_t nb_samples = 0;
      for (Facet_const_iterator fi = facets_begin();
        fi != facets_end(); ++fi) {
        // facet in links
        for (Link_iter_list_const_iter it = fi->facet_in_links().begin();
          it != fi->facet_in_links().end(); ++it) {
          Link_list_const_iter lit = *it;
          const Point &start = lit->second.first;
          const Point &end = lit->second.second;
          rms_distance += CGAL::squared_distance(start, end);
        }
        nb_samples += fi->facet_in_links().size();
        // edge in links
        for (Link_iter_list_const_iter it = fi->edge_in_links().begin();
          it != fi->edge_in_links().end(); ++it) {
          Link_list_const_iter lit = *it;
          const Point &start = lit->second.first;
          const Point &end = lit->second.second;
          rms_distance += CGAL::squared_distance(start, end);
        }
        nb_samples += fi->edge_in_links().size();
        // vertex in links
        for (auto it = fi->vertex_in_links().begin();
          it != fi->vertex_in_links().end(); ++it) {
          Link *link = *it;
          const Point &start = link->second.first;
          const Point &end = link->second.second;
          rms_distance += CGAL::squared_distance(start, end);
        }
        nb_samples += fi->vertex_in_links().size();
      }
      // facet out links
      for (Facet_const_iterator fi = facets_begin();
        fi != facets_end(); ++fi) {
        for (Link_list_const_iter it = fi->facet_out_links().begin();
          it != fi->facet_out_links().end(); ++it) {
          const Point &start = it->second.first;
          const Point &end = it->second.second;
          rms_distance += CGAL::squared_distance(start, end);
        }
        nb_samples += fi->facet_out_links().size();
      }
      // edge out links
      for (Edge_const_iterator ei = edges_begin(); ei != edges_end(); ++ei) {
        Halfedge_const_handle hh = ei;
        if (hh->normal_dihedral() == -1.0) {
          hh = hh->opposite();
        }
        for (Link_list_const_iter it = hh->edge_out_links().begin();
          it != hh->edge_out_links().end(); ++it) {
          const Point &start = it->second.first;
          const Point &end = it->second.second;
          rms_distance += CGAL::squared_distance(start, end);
        }
        nb_samples += hh->edge_out_links().size();
      }
      // vertex out links
      for (Vertex_const_iterator vi = vertices_begin();
        vi != vertices_end(); ++vi) {
        const Point &start = vi->vertex_out_link().second.first;
        const Point &end = vi->vertex_out_link().second.second;
        rms_distance += CGAL::squared_distance(start, end);
      }
      nb_samples += size_of_vertices();
      if (nb_samples == 0) {      // invalid case
        return -1.0;
      }
      else {
        return CGAL::sqrt(rms_distance / nb_samples);
      }
    }

    void trace_edge_length() const {
      if (size_of_halfedges() == 0) {
        return;
      }
      FT sum_length = 0.0;
      FT min_length = MAX_VALUE;
      FT max_length = 0.0;
      int nb_edges = 0;
      for (Edge_const_iterator ei = edges_begin(); ei != edges_end(); ++ei) {
        const FT edge_length = length(ei);
        sum_length += edge_length;
        min_length = std::min(min_length, edge_length);
        max_length = std::max(max_length, edge_length);
        ++nb_edges;
      }
      std::cout << "Min edge length: " << min_length << std::endl;
      std::cout << "Max edge length: " << max_length << std::endl;
      std::cout << "Average edge length: "
        << sum_length / nb_edges << std::endl;
    }

    FT get_regular_vertex_ratio() const {
      int nb_regular_vertices = 0;
      for (Vertex_const_iterator vi = vertices_begin();
        vi != vertices_end(); ++vi) {
        nb_regular_vertices += is_regular(vi);
      }
      if (size_of_vertices() == 0) {    // invalid case
        return -1.0;
      }
      else {
        return static_cast<double>(nb_regular_vertices) / size_of_vertices();
      }
    }

    FT get_smaller_angle_ratio(FT angle) const {
      FT radian = to_radian(angle);
      int nb_facets = 0;
      for (Facet_const_iterator fi = facets_begin();
        fi != facets_end(); ++fi) {
        FT smallest_radian = get_smallest_radian(fi);
        nb_facets += smallest_radian < radian;
      }
      if (size_of_facets() == 0) {    // invalid case
        return -1.0;
      }
      else {
        return static_cast<double>(nb_facets) / size_of_facets();
      }
    }

    int nb_boundaries() {    // count #boundaries
      unsigned int nb = 0;
      tag_halfedges(0);
      Halfedge_iterator hi;
      for (hi = halfedges_begin(); hi != halfedges_end(); ++hi) {
        if (hi->is_border() && hi->tag() == 0) {
          ++nb;
          Halfedge_handle curr = hi;
          do {
            curr = curr->next();
            curr->tag() = 1;
          } while (curr != hi);
        }
      }
      return nb;
    }

    int nb_components() {    // count #components
      unsigned int nb = 0;
      tag_facets(0);
      for (Facet_iterator fi = facets_begin(); fi != facets_end(); ++fi) {
        if (fi->tag() == 0) {
          nb++;
          tag_component(fi, 0, 1);
        }
      }
      return nb;
    }

    void tag_component(Facet_handle pSeedFacet,
                       const int tag_free, const int tag_done) {
      pSeedFacet->tag() = tag_done;
      Facet_list facets;
      facets.push_front(pSeedFacet);
      while (!facets.empty()) {
        Facet_handle pFacet = facets.front();
        facets.pop_front();
        pFacet->tag() = tag_done;
        Halfedge_around_facet_circulator pHalfedge = pFacet->facet_begin();
        Halfedge_around_facet_circulator end = pHalfedge;
        CGAL_For_all(pHalfedge, end) {
          Facet_handle pNFacet = pHalfedge->opposite()->facet();
          if (pNFacet != NULL && pNFacet->tag() == tag_free) {
            facets.push_front(pNFacet);
            pNFacet->tag() = tag_done;
          }
        }
      }
    }

    void tag_facets(const int tag) {
      for (Facet_iterator fi = facets_begin(); fi != facets_end(); ++fi) {
        fi->tag() = tag;
      }
    }

    void tag_halfedges(const int tag) {
      for (Halfedge_iterator hi = halfedges_begin();
        hi != halfedges_end(); ++hi) {
        hi->tag() = tag;
      }
    }

    // element collections
    void get_crease_edges_in_one_ring(Vertex_const_handle vh,
        Halfedge_const_list *crease_edges) const {
      Halfedge_around_vertex_const_circulator vcirc = vh->vertex_begin();
      Halfedge_around_vertex_const_circulator vend = vcirc;
      CGAL_For_all(vcirc, vend) {
        Halfedge_const_handle hh = vcirc;
        if (hh->normal_dihedral() == -1.0) {
          hh = hh->opposite();
        }
        if (hh->is_crease()) {
          crease_edges->push_back(vcirc);
        }
      }
    }

    void get_incident_facets(Facet_handle fh,
                             std::set<Facet_handle> *facets) const {
      Halfedge_handle hh = fh->halfedge();
      // the first vertex
      Halfedge_around_vertex_circulator vcirc = hh->vertex()->vertex_begin();
      Halfedge_around_vertex_circulator vend = vcirc;
      CGAL_For_all(vcirc, vend) {
        if (!vcirc->is_border()) {
          facets->insert(vcirc->facet());
        }
      }
      // the second vertex
      vcirc = hh->next()->vertex()->vertex_begin();
      vend = vcirc;
      CGAL_For_all(vcirc, vend) {
        if (!vcirc->is_border()) {
          facets->insert(vcirc->facet());
        }
      }
      // the third vertex
      vcirc = hh->prev()->vertex()->vertex_begin();
      vend = vcirc;
      CGAL_For_all(vcirc, vend) {
        if (!vcirc->is_border()) {
          facets->insert(vcirc->facet());
        }
      }
    }

    void get_effective_edges_in_one_ring(FT feature_control_delta,
        Vertex_const_handle vh, Halfedge_const_list *effective_edges) const {
      FT f1 = vh->feature_intensity() * feature_control_delta;
      FT f2 = (vh->max_halfedge_dihedral() + 1) * feature_control_delta;
      Halfedge_around_vertex_const_circulator vcirc = vh->vertex_begin();
      Halfedge_around_vertex_const_circulator vend = vcirc;
      CGAL_For_all(vcirc, vend) {
        Vertex_const_handle v = get_source_vertex(vcirc);
        bool b1 = v->feature_intensity() >= f1;
        Halfedge_const_handle hh = vcirc;
        if (hh->normal_dihedral() == -1.0) {
          hh = hh->opposite();
        }
        bool b2 = (hh->normal_dihedral() + 1) >= f2;
        if (b1 && b2) {
          effective_edges->push_back(hh);
        }
      }
    }

    void get_effective_edges_in_one_ring(FT feature_control_delta,
        Vertex_handle vh, Halfedge_list *effective_edges) {
      FT f1 = vh->feature_intensity() * feature_control_delta;
      FT f2 = (vh->max_halfedge_dihedral() + 1) * feature_control_delta;
      Halfedge_around_vertex_circulator vcirc = vh->vertex_begin();
      Halfedge_around_vertex_circulator vend = vcirc;
      CGAL_For_all(vcirc, vend) {
        Vertex_handle v = get_source_vertex(vcirc);
        bool b1 = v->feature_intensity() >= f1;
        Halfedge_handle hh = vcirc;
        if (hh->normal_dihedral() == -1.0) {
          hh = hh->opposite();
        }
        bool b2 = (hh->normal_dihedral() + 1) >= f2;
        if (b1 && b2) {
          effective_edges->push_back(hh);
        }
      }
    }

    void extend_facets_by_one_ring(std::set<Facet_handle> *facets) const {
      std::set<Vertex_handle> vertices;
      collect_vertices(*facets, &vertices);
      for (auto it = vertices.begin(); it != vertices.end(); ++it) {
        collect_one_ring_facets_incident_to_vertex(*it, facets);
      }
    }

    // local operators
    bool check_link_test(Halfedge_const_handle h) const {
      // h = (p -> q)
      Vertex_const_handle vp = get_source_vertex(h);
      Vertex_const_handle vq = get_target_vertex(h);
      if (h->opposite()->is_border()) {
        // (p, q, r)
        Vertex_const_handle vr = get_opposite_vertex(h);
        Halfedge_around_vertex_const_circulator vcirc = vp->vertex_begin();
        Halfedge_around_vertex_const_circulator vend = vcirc;
        CGAL_For_all(vcirc, vend) {
          Vertex_const_handle v = get_source_vertex(vcirc);
          if (v == vq || v == vr) {
            continue;
          }
          if (are_neighbors(v, vq)) {
            return false;
          }
        }
        return true;
      }
      // (p, q, r)
      Vertex_const_handle vr = get_opposite_vertex(h);
      // (q, p, s)
      Vertex_const_handle vs = get_opposite_vertex(h->opposite());
      //
      Halfedge_around_vertex_const_circulator vcirc = vp->vertex_begin();
      Halfedge_around_vertex_const_circulator vend = vcirc;
      CGAL_For_all(vcirc, vend) {
        Vertex_const_handle v = get_source_vertex(vcirc);
        if (v == vq || v == vr || v == vs) {
          continue;
        }
        if (are_neighbors(v, vq)) {
          return false;
        }
      }
      return true;
    }

    bool join_facets_before_collapse(Halfedge_handle hh) {
      assert(!hh->is_border());
      if (hh->opposite()->is_border()) {
        return join_facets_before_collapse_boundary_case(hh);
      }
      else {
        return join_facets_before_collapse_inner_case(hh);
      }
    }

    bool join_facets_before_collapse_inner_case(Halfedge_handle hh) {
      Halfedge_handle qr = hh->next();
      Halfedge_handle rp = hh->prev();
      Halfedge_handle ps = hh->opposite()->next();
      Halfedge_handle sq = hh->opposite()->prev();
      Facet_handle frq = qr->opposite()->facet();
      Facet_handle fpr = rp->opposite()->facet();
      Facet_handle fsp = ps->opposite()->facet();
      Facet_handle fqs = sq->opposite()->facet();
      // check if faces are different
      if (valid_join_facet_pair(fsp, frq)) {
        if (valid_join_facet(ps) && valid_join_facet(qr)) {
          this->join_facet(ps);
          this->join_facet(qr);
          return true;
        }
      }
      if (valid_join_facet_pair(fsp, fpr)) {
        if (valid_join_facet(ps) && valid_join_facet(rp)) {
          this->join_facet(ps);
          this->join_facet(rp);
          return true;
        }
      }
      if (valid_join_facet_pair(fqs, frq)) {
        if (valid_join_facet(sq) && valid_join_facet(qr)) {
          this->join_facet(sq);
          this->join_facet(qr);
          return true;
        }
      }
      if (valid_join_facet_pair(fqs, fpr)) {
        if (valid_join_facet(sq) && valid_join_facet(rp)) {
          this->join_facet(sq);
          this->join_facet(rp);
          return true;
        }
      }
      std::cerr << red << "join facets failed" << white << std::endl;
      return false;
    }

    bool join_facets_before_collapse_boundary_case(Halfedge_handle hh) {
      CGAL_precondition(hh->opposite()->is_border());
      // given edge = pq
      Halfedge_handle qr = hh->next();
      Halfedge_handle rp = hh->next()->next();
      Facet_handle frq = qr->opposite()->facet();
      Facet_handle fpr = rp->opposite()->facet();
      if (frq != Facet_handle()) {
        if (valid_join_facet(qr)) {
          this->join_facet(qr);
          return true;
        }
      }
      else if (fpr != Facet_handle()) {
        if (valid_join_facet(rp)) {
          this->join_facet(rp);
          return true;
        }
      }
      // should never come here
      std::cerr << red << "Error in join_facets_before_collapse_boundary_case"
        << white << std::endl;
      return false;
    }

    bool valid_join_facet_pair(Facet_handle f1, Facet_handle f2) const {
      if (f1 == f2) {
        return false;
      }
      if (f1 == Facet_handle()) {
        return false;
      }
      if (f2 == Facet_handle()) {
        return false;
      }
      return true;
    }

    bool valid_join_facet(Halfedge_handle hh) const {
      if (circulator_size(hh->vertex_begin()) < 3) {
        return false;
      }
      if (circulator_size(hh->opposite()->vertex_begin()) < 3) {
        return false;
      }
      return true;
    }

    bool edge_flip_would_improve_valence(Halfedge_const_handle hh) const {
      // incident vertices
      Vertex_const_handle a = get_target_vertex(hh);
      Vertex_const_handle b = get_opposite_vertex(hh);
      Vertex_const_handle c = get_source_vertex(hh);
      Vertex_const_handle d = get_opposite_vertex(hh->opposite());
      // target valences
      int tar_var_a = is_on_boundary(a) ? 4 : 6;
      int tar_var_b = is_on_boundary(b) ? 4 : 6;
      int tar_var_c = is_on_boundary(c) ? 4 : 6;
      int tar_var_d = is_on_boundary(d) ? 4 : 6;
      // real valence
      int var_a = static_cast<int>(a->vertex_degree());
      int var_b = static_cast<int>(b->vertex_degree());
      int var_c = static_cast<int>(c->vertex_degree());
      int var_d = static_cast<int>(d->vertex_degree());
      int deviation_pre = std::abs(var_a - tar_var_a) +
        std::abs(var_b - tar_var_b) +
        std::abs(var_c - tar_var_c) +
        std::abs(var_d - tar_var_d);
      // after flipping, the valences of a, b, c and d change
      var_a--; var_c--; var_b++; var_d++;
      int deviation_post = std::abs(var_a - tar_var_a) +
        std::abs(var_b - tar_var_b) +
        std::abs(var_c - tar_var_c) +
        std::abs(var_d - tar_var_d);
      return deviation_pre > deviation_post;
    }

    bool edge_flip_would_improve_radian(Halfedge_const_handle hh) const {
      // precondition: hh or hh->opposite() cannot be boundary edges
      const Point &p = get_source_vertex(hh)->point();
      const Point &q = get_target_vertex(hh)->point();
      const Point &s = get_opposite_vertex(hh)->point();
      const Point &t = get_opposite_vertex(hh->opposite())->point();
      FT min_radian_before = get_smallest_radian(p, q, s);
      min_radian_before = std::min(min_radian_before,
        get_smallest_radian(q, p, t));
      FT min_radian_after = get_smallest_radian(t, s, p);
      min_radian_after = std::min(min_radian_after,
        get_smallest_radian(s, t, q));
      return min_radian_before < min_radian_after;
    }

    // visual element
    void compute_polyhedron_facets(RenderType render_type,
        FT sum_theta_value, FT dihedral_theta_value,
        FT max_error_threshold_value, bool is_input,
        std::vector<float> *pos_faces, std::vector<float> *pos_face_normals,
        std::vector<float> *pos_face_colors) const {
      FT min_value = 0;
      FT max_value = (sum_theta_value + 1) * (dihedral_theta_value + 1) - 1;
      Color input_color(150, 150, 200), remesh_color(200, 150, 150);
      for (Facet_const_iterator fi = facets_begin();
        fi != facets_end(); ++fi) {
        const Normal &normal = fi->normal();
        Halfedge_const_handle hh = fi->halfedge();
        Vertex_const_handle vh = get_target_vertex(hh);
        FT error = CGAL::sqrt(fi->max_squared_error());
        if (error > max_error_threshold_value) {
          error = max_error_threshold_value;
        }
        else {
          // if error == max_error_threshold_value, normalizeit to 0.5
          error = error / (2 * max_error_threshold_value);  // normalize
        }
        Color error_color = get_rgb_color(0, error, 1.0);
        for (int i = 0; i <= 2; ++i) {
          compute_vertex(vh, pos_faces);
          compute_normal(normal, pos_face_normals);
          switch (render_type) {
          case RenderType::k_plain_facets: {
            is_input ? compute_color(input_color, pos_face_colors) :
              compute_color(remesh_color, pos_face_colors);
            break;
          }
          case RenderType::k_ifi_facets: {
            Color color = get_vertex_sample_normalized_color(vh,
              RenderType::k_ifi_facets, min_value, max_value, 240);
            compute_color(color, pos_face_colors);
            break;
          }
          case RenderType::k_mr_facets: {
            compute_color(error_color, pos_face_colors);
            break;
          }
          default: {
            break;
          }
          }
          hh = hh->next();
          vh = get_target_vertex(hh);
        }
      }
    }

    void compute_all_voronois(RenderType render_type,
        FT sum_theta_value, FT dihedral_theta_value,
        std::vector<float> *pos_faces, std::vector<float> *pos_face_normals,
        std::vector<float> *pos_face_colors,
        std::vector<float> *pos_boundaries,
        std::vector<float> *pos_samples) const {
      // step 1: compute all the minimal value and maximal value
      FT min_value = 0.0, max_value = 0.0;
      switch (render_type) {
      case RenderType::k_feature_intensity:
        max_value = (sum_theta_value + 1) * (dihedral_theta_value + 1) - 1;
        break;
      case RenderType::k_capacity:
        max_value = 2.0;
        break;
      case RenderType::k_weight:
        max_value = get_max_sample_weight();
        break;
      default:
        break;
      }
      // step 2: compute all sample cells and cell boundaries
      for (Facet_const_iterator fi = facets_begin();
        fi != facets_end(); ++fi) {
        const Normal &normal = fi->normal();
        Point_list samples;
        Color_list colors;
        get_all_sample_normalized_colors(fi, render_type, min_value, max_value,
          240, &samples, &colors);
        Bvd bvd(triangle(fi));
        bvd.compute_voronoi_cells_and_boundaries(samples, normal, colors,
          pos_faces, pos_face_normals, pos_face_colors, pos_boundaries);
      }
      // step 3: compute all samples
      compute_all_samples(pos_samples);
    }

    void compute_vertex_voronois(RenderType render_type,
        bool inherit_element_types, FT feature_control_delta,
        FT sum_theta_value, FT dihedral_theta_value,
        std::vector<float> *pos_faces, std::vector<float> *pos_face_normals,
        std::vector<float> *pos_face_colors,
        std::vector<float> *pos_boundaries,
        std::vector<float> *pos_samples) const {
      // step 1: comput the minimal and maximal value
      FT min_value = 0.0, max_value = 0.0;
      switch (render_type) {
      case RenderType::k_classifications:
        break;
      case RenderType::k_gaussian_curvature:
        max_value = sum_theta_value;
        break;
      case RenderType::k_maximal_halfedge_dihedral:
        max_value = dihedral_theta_value;
        break;
      case RenderType::k_feature_intensity:
        max_value = (sum_theta_value + 1) * (dihedral_theta_value + 1) - 1;
        break;
      case RenderType::k_capacity:
        //max_value = 2 * get_max_sample_capacity();
        max_value = 2 * get_max_vertex_sample_capacity();
        break;
      case RenderType::k_weight:
        //max_value = get_max_sample_weight();
        max_value = get_max_vertex_sample_weight();
        break;
      default:
        break;
      }
      // step 2: compute the cell and cell boundaries
      Color color;
      for (Vertex_const_iterator vi = vertices_begin();
        vi != vertices_end(); ++vi) {
        if (render_type == RenderType::k_classifications) {
          color = get_vertex_classification_color(inherit_element_types,
            feature_control_delta, vi);
        }
        else {
          color = get_vertex_sample_normalized_color(
            vi, render_type, min_value, max_value, 240);
        }
        const Point &p = vi->point();
        Halfedge_around_vertex_const_circulator vcirc = vi->vertex_begin();
        Halfedge_around_vertex_const_circulator vend = vcirc;
        CGAL_For_all(vcirc, vend) {
          if (!vcirc->is_border()) {
            // step 1: compute the facets
            Facet_const_handle fh = vcirc->facet();
            const Normal &normal = fh->normal();
            const Point p1 = midpoint(vcirc);
            const Point c = centroid(fh);
            const Point p2 = midpoint(vcirc->next());
            compute_triangle(p1, p, c, normal, color,
              pos_faces, pos_face_normals, pos_face_colors);
            compute_triangle(c, p, p2, normal, color,
              pos_faces, pos_face_normals, pos_face_colors);
            // step 2: compute the edges
            compute_segment(p1, c, pos_boundaries);
            compute_segment(c, p2, pos_boundaries);
          }
        }
      }
      // step 3: compute vertices
      compute_vertices(pos_samples);
    }

    void compute_edge_voronois(RenderType render_type,
        FT sum_theta_value, FT dihedral_theta_value,
        std::vector<float> *pos_faces, std::vector<float> *pos_face_normals,
        std::vector<float> *pos_face_colors,
        std::vector<float> *pos_boundaries,
        std::vector<float> *pos_samples) const {
      if (render_type == RenderType::k_normal_dihedral) {
        // only render two triangles for each edge. Do not render samples
        compute_edge_normal_dihedrals(dihedral_theta_value,
          pos_faces, pos_face_normals, pos_face_colors, pos_boundaries);
      }
      else {
        // we need to render triangles for each sample. render the samples
        compute_edge_sample_properties(render_type, sum_theta_value,
          dihedral_theta_value, pos_faces, pos_face_normals,
          pos_face_colors, pos_boundaries, pos_samples);
      }
    }

    void compute_facet_voronois(RenderType render_type,
        FT sum_theta_value, FT dihedral_theta_value,
        std::vector<float> *pos_faces, std::vector<float> *pos_face_normals,
        std::vector<float> *pos_face_colors,
        std::vector<float> *pos_boundaries,
        std::vector<float> *pos_samples) const {
      // step 1: compute the minimal and maximal values
      FT min_value = 0.0, max_value = 0.0;
      switch (render_type) {
      case RenderType::k_feature_intensity:
        max_value = (sum_theta_value + 1) * (dihedral_theta_value + 1) - 1;
        break;
      case RenderType::k_capacity:
        //max_value = 2 * get_max_sample_capacity();
        max_value = 2 * get_max_facet_sample_capacity();
        break;
      case RenderType::k_weight:
        //max_value = get_max_sample_weight();
        max_value = get_max_facet_sample_weight();
        break;
      default:
        break;
      }
      // step 2: compute facet voronoi cell and boundaries
      for (Facet_const_iterator fi = facets_begin();
        fi != facets_end(); ++fi) {
        const Normal &normal = fi->normal();
        Point_list samples;
        Color_list colors;
        get_facet_sample_normalized_colors(fi, render_type, min_value,
          max_value, 240, &samples, &colors);
        Bvd bvd(triangle(fi));
        bvd.compute_voronoi_cells_and_boundaries(samples, normal, colors,
          pos_faces, pos_face_normals, pos_face_colors, pos_boundaries);
      }
      // step 3: compute samples
      compute_facet_samples(pos_samples);
    }

    void compute_edge_normal_dihedrals(FT dihedral_theta_value,
        std::vector<float> *pos_faces, std::vector<float> *pos_face_normals,
        std::vector<float> *pos_face_colors,
        std::vector<float> *pos_boundaries) const {
      FT min_value = 0.0, max_value = dihedral_theta_value;
      Color color;
      // step 1: compute the faces
      for (Edge_const_iterator ei = edges_begin(); ei != edges_end(); ++ei) {
        Halfedge_const_handle hh = ei;
        if (hh->normal_dihedral() == -1.0) {
          hh = hh->opposite();
        }
        FT value = hh->normal_dihedral() / max_value; // normalized value
        color = get_rgb_color(240, value, 1.0);
        Facet_const_handle fh = hh->facet();            // the first triangle
        const Normal &normal = fh->normal();
        const Point &a = get_target_vertex(hh)->point();
        const Point b = centroid(fh);
        const Point &c = get_source_vertex(hh)->point();
        compute_triangle(a, b, c, normal, color,
          pos_faces, pos_face_normals, pos_face_colors);
        if (!hh->opposite()->is_border()) {             // the second triangle
          hh = hh->opposite();
          fh = hh->facet();
          const Normal &normal = fh->normal();
          const Point &a = get_target_vertex(hh)->point();
          const Point b = centroid(fh);
          const Point &c = get_source_vertex(hh)->point();
          compute_triangle(a, b, c, normal, color,
            pos_faces, pos_face_normals, pos_face_colors);
        }
      }
      // step 2: compute the edges
      for (Facet_const_iterator fi = facets_begin();
        fi != facets_end(); ++fi) {
        Halfedge_const_handle hh = fi->halfedge();
        const Point &a = get_source_vertex(hh)->point();
        const Point &b = get_target_vertex(hh)->point();
        const Point &c = get_opposite_vertex(hh)->point();
        Point d = centroid(fi);
        compute_segment(d, a, pos_boundaries);
        compute_segment(d, b, pos_boundaries);
        compute_segment(d, c, pos_boundaries);
      }
    }

    void compute_edge_sample_properties(RenderType render_type,
        FT sum_theta_value, FT dihedral_theta_value,
        std::vector<float> *pos_faces, std::vector<float> *pos_face_normals,
        std::vector<float> *pos_face_colors,
        std::vector<float> *pos_boundaries,
        std::vector<float> *pos_samples) const {
      // step 1: get the minimal and maximal values
      FT min_value = 0.0, max_value = 0.0;
      switch (render_type) {
      case RenderType::k_feature_intensity:
        max_value = (sum_theta_value + 1) * (dihedral_theta_value + 1) - 1;
        break;
      case RenderType::k_capacity:
        //max_value = 2 * get_max_sample_capacity();
        max_value = 2 * get_max_edge_sample_capacity();
        break;
      case RenderType::k_weight:
        //max_value = get_max_sample_weight();
        max_value = get_max_edge_sample_weight();
        break;
      default:
        break;
      }
      // step 2: compute the faces and samples
      for (Edge_const_iterator ei = edges_begin(); ei != edges_end(); ++ei) {
        Halfedge_const_handle hh = ei;
        if (hh->normal_dihedral() == -1.0) {
          hh = hh->opposite();
        }
        // step 2.1: compute the samples
        for (Link_list_const_iter it = hh->edge_out_links().begin();
          it != hh->edge_out_links().end(); ++it) {
          compute_point(it->second.first, pos_samples);
        }
        // step 2.2: compute the faces and edges
        Link_list_const_iter first = hh->edge_out_links().begin();
        Link_list_const_iter second = first;
        ++second;
        Point_list points;
        points.push_back(get_source_vertex(hh)->point());
        for (; second != hh->edge_out_links().end(); ++first, ++second) {
          points.push_back(CGAL::midpoint(first->second.first,
            second->second.second));
        }
        points.push_back(get_target_vertex(hh)->point());
        Color_list colors;
        get_edge_sample_normalized_colors(hh, render_type, min_value,
          max_value, 240, &colors);
        Facet_const_handle fh = hh->facet(); // samples of the first triangle
        Point c = centroid(fh);
        const Normal &normal = fh->normal();
        Point_iter it1 = points.begin();
        Point_iter it2 = it1;
        ++it2;
        Color_const_iter cit = colors.begin();
        for (; it2 != points.end(); ++it1, ++it2, ++cit) {
          if (it1 != points.begin()) {
            compute_segment(c, *it1, pos_boundaries);
          }
          compute_triangle(*it1, *it2, c, normal, *cit,
            pos_faces, pos_face_normals, pos_face_colors);
        }
        if (!hh->opposite()->is_border()) { // samples of the second triangle
          fh = hh->opposite()->facet();
          c = centroid(fh);
          const Normal &normal = fh->normal();
          it1 = points.begin();
          it2 = it1;
          ++it2;
          cit = colors.begin();
          for (; it2 != points.end(); ++it1, ++it2, ++cit) {
            if (it1 != points.begin()) {
              compute_segment(c, *it1, pos_boundaries);
            }
            compute_triangle(*it1, *it2, c, normal, *cit,
              pos_faces, pos_face_normals, pos_face_colors);
          }
        }
      }
      // step 3: compute the rest of the edges
      for (Facet_const_iterator fi = facets_begin();
        fi != facets_end(); ++fi) {
        Point c = centroid(fi);
        Halfedge_const_handle hh = fi->halfedge();
        for (int i = 0; i <= 2; ++i) {
          compute_segment(c, get_target_vertex(hh)->point(), pos_boundaries);
          hh = hh->next();
        }
      }
    }

    void compute_all_samples(std::vector<float> *pos_samples) const {
      // step 1: facet out samples
      compute_facet_samples(pos_samples);
      // step 2: edge out samples
      for (Edge_const_iterator ei = edges_begin(); ei != edges_end(); ++ei) {
        Halfedge_const_iterator hh = ei;
        if (hh->normal_dihedral() == -1.0) {
          hh = hh->opposite();
        }
        for (Link_list_const_iter it = hh->edge_out_links().begin();
          it != hh->edge_out_links().end(); ++it) {
          const Point &p = it->second.first;
          compute_point(p, pos_samples);
        }
      }
      // step 3: vertex out samples
      for (Vertex_const_iterator vi = vertices_begin();
        vi != vertices_end(); ++vi) {
        const Point &p = vi->vertex_out_link().second.first;
        compute_point(p, pos_samples);
      }
    }

    void compute_vertices(std::vector<float> *pos_samples) const {
      for (Vertex_const_iterator vi = vertices_begin();
        vi != vertices_end(); ++vi) {
        const Point &p = vi->point();
        compute_point(p, pos_samples);
      }
    }

    void inline compute_triangle(const Point &p1, const Point &p2,
        const Point &p3, const Normal &normal, const Color &color,
        std::vector<float> *pos_faces, std::vector<float> *pos_face_normals,
        std::vector<float> *pos_face_colors) const {
      compute_triangle_point(p1, normal, color,
        pos_faces, pos_face_normals, pos_face_colors);
      compute_triangle_point(p2, normal, color,
        pos_faces, pos_face_normals, pos_face_colors);
      compute_triangle_point(p3, normal, color,
        pos_faces, pos_face_normals, pos_face_colors);
    }

    void inline compute_triangle_point(const Point &p, const Normal &normal,
        const Color &color, std::vector<float> *pos_faces,
        std::vector<float> *pos_face_normals,
        std::vector<float> *pos_face_colors) const {
      compute_point(p, pos_faces);
      compute_normal(normal, pos_face_normals);
      compute_color(color, pos_face_colors);
    }

    void inline compute_halfedge(Halfedge_const_handle hh,
                                 std::vector<float> *pos) const {
      const Point &p = get_source_vertex(hh)->point();
      const Point &q = get_target_vertex(hh)->point();
      compute_segment(p, q, pos);
    }

    void inline compute_segment(const Point &p, const Point &q,
                                std::vector<float> *pos) const {
      compute_point(p, pos);
      compute_point(q, pos);
    }

    void inline compute_vertex(Vertex_const_handle vh,
                               std::vector<float> *pos) const {
      compute_point(vh->point(), pos);
    }

    void inline compute_point(const Point &point,
                              std::vector<float> *pos) const {
      pos->push_back(point.x());
      pos->push_back(point.y());
      pos->push_back(point.z());
    }

    void inline compute_normal(const Normal &normal,
                               std::vector<float> *pos) const {
      pos->push_back(normal.x());
      pos->push_back(normal.y());
      pos->push_back(normal.z());
    }

    void inline compute_color(const Color &color,
                              std::vector<float> *pos) const {
      pos->push_back(color.red() / 255.0f);
      pos->push_back(color.green() / 255.0f);
      pos->push_back(color.blue() / 255.0f);
    }

    FT get_max_sample_capacity() const {
      FT max_capacity = 0.0;
      max_capacity = CGAL::max(max_capacity, get_max_facet_sample_capacity());
      max_capacity = CGAL::max(max_capacity, get_max_edge_sample_capacity());
      max_capacity = CGAL::max(max_capacity, get_max_vertex_sample_capacity());
      return max_capacity;
    }

    FT get_max_facet_sample_capacity() const {
      FT max_capacity = 0.0;
      for (Facet_const_iterator fi = facets_begin();
        fi != facets_end(); ++fi) {
        FT capacity = area(fi) / fi->facet_out_links().size();
        max_capacity = CGAL::max(max_capacity, capacity);
      }
      return max_capacity;
    }

    FT get_max_edge_sample_capacity() const {
      FT max_capacity = 0.0;
      for (Edge_const_iterator ei = edges_begin(); ei != edges_end(); ++ei) {
        Halfedge_const_handle hh = ei;
        if (hh->normal_dihedral() == -1.0) {
          hh = hh->opposite();
        }
        FT facet_area = area(hh->facet());
        if (!hh->opposite()->is_border()) {
          facet_area += area(hh->opposite()->facet());
        }
        FT capacity = facet_area / (3 * hh->edge_out_links().size());
        max_capacity = CGAL::max(max_capacity, capacity);
      }
      return max_capacity;
    }

    FT get_max_vertex_sample_capacity() const {
      FT max_capacity = 0.0;
      FT weight = 0.0, feature_intensity = 0.0;
      for (Vertex_const_iterator vi = vertices_begin();
        vi != vertices_end(); ++vi) {
        weight = vi->vertex_out_link().first;
        feature_intensity = vi->feature_intensity();
        max_capacity = CGAL::max(max_capacity, weight / feature_intensity);
      }
      return max_capacity;
    }

    FT get_max_sample_weight() const {
      FT max_weight = 0.0;
      max_weight = CGAL::max(max_weight, get_max_facet_sample_weight());
      max_weight = CGAL::max(max_weight, get_max_edge_sample_weight());
      max_weight = CGAL::max(max_weight, get_max_vertex_sample_weight());
      return max_weight;
    }

    FT get_max_facet_sample_weight() const {
      FT max_weight = 0.0;
      for (Facet_const_iterator fi = facets_begin();
        fi != facets_end(); ++fi) {
        for (Link_list_const_iter it = fi->facet_out_links().begin();
          it != fi->facet_out_links().end(); ++it) {
          max_weight = CGAL::max(max_weight, it->first);
        }
      }
      return max_weight;
    }

    FT get_max_edge_sample_weight() const {
      FT max_weight = 0.0;
      for (Halfedge_const_iterator hi = edges_begin();
        hi != edges_end(); ++hi) {
        if (hi->normal_dihedral() != -1.0) {
          for (Link_list_const_iter it = hi->edge_out_links().begin();
            it != hi->edge_out_links().end(); ++it) {
            max_weight = CGAL::max(max_weight, it->first);
          }
        }
      }
      return max_weight;
    }

    FT get_max_vertex_sample_weight() const {
      FT max_weight = 0.0;
      for (Vertex_const_iterator vi = vertices_begin();
        vi != vertices_end(); ++vi) {
        max_weight = CGAL::max(max_weight, vi->vertex_out_link().first);
      }
      return max_weight;
    }

    void get_all_sample_normalized_colors(Facet_const_handle fh,
        RenderType render_type, FT min_value, FT max_value,
        FT h, Point_list *samples, Color_list *colors) const {
      // calculate the samples and their weights
      std::map<Point, double, Point_Comp> disturbed_border_samples;
      get_disturbed_border_samples_with_weights(fh, 0.01, render_type,
        &disturbed_border_samples);
      for (Link_list_const_iter it = fh->facet_out_links().begin();
        it != fh->facet_out_links().end(); ++it) {
        if (disturbed_border_samples.find(it->second.first) ==
          disturbed_border_samples.end()) {
          switch (render_type) {
          case RenderType::k_feature_intensity:
            disturbed_border_samples[it->second.first] = it->first - 1.0;
            break;
          case RenderType::k_weight:
            disturbed_border_samples[it->second.first] = it->first;
            break;
          case RenderType::k_capacity:
            disturbed_border_samples[it->second.first] = 1.0;
            break;
          default:
            break;
          }
        }
      }
      // normalize the weights, and convert to colors
      for (auto it = disturbed_border_samples.begin();
        it != disturbed_border_samples.end(); ++it) {
        if (it->second > max_value) {
          it->second = max_value;
        }
        if (it->second < min_value) {
          it->second = min_value;
        }
        it->second = (it->second - min_value) / (max_value - min_value);
        samples->push_back(it->first);
        colors->push_back(get_rgb_color(h, it->second, 1.0));
      }
    }

    void get_facet_sample_normalized_colors(Facet_const_handle fh,
        RenderType render_type, FT min_value, FT max_value,
        FT h, Point_list *samples, Color_list *colors) const {
        FT capacity = area(fh) / fh->facet_out_links().size();
      // step 1: calculate the original values
      std::list<double> values;
      for (Link_list_const_iter it = fh->facet_out_links().begin();
        it != fh->facet_out_links().end(); ++it) {
        samples->push_back(it->second.first);
        switch (render_type) {
        case RenderType::k_feature_intensity:
          values.push_back(it->first / capacity - 1.0);
          break;
        case RenderType::k_capacity:
          values.push_back(capacity);
          break;
        case RenderType::k_weight:
          values.push_back(it->first);
          break;
        default:
          break;
        }
      }
      // step 2: normalize the values, and convert to colors
      for (auto it = values.begin(); it != values.end(); ++it) {
        if (*it > max_value) {
          *it = max_value;
        }
        if (*it < min_value) {
          *it = min_value;
        }
        *it = (*it - min_value) / (max_value - min_value);
        colors->push_back(get_rgb_color(h, *it, 1.0));
      }
    }

    void get_edge_sample_normalized_colors(Halfedge_const_handle hh,
        RenderType render_type, FT min_value, FT max_value,
        FT h, Color_list *colors) const {
        FT facet_area = area(hh->facet());      // calculate the facet area
      if (!hh->opposite()->is_border()) {
        facet_area += area(hh->opposite()->facet());
      }
      FT capacity = facet_area / (3 * hh->edge_out_links().size());
      std::list<double> values;
      // step 1: calculate the original values
      for (Link_list_const_iter it = hh->edge_out_links().begin();
        it != hh->edge_out_links().end(); ++it) {
        switch (render_type) {
        case RenderType::k_feature_intensity:
          values.push_back(it->first / capacity - 1.0);
          break;
        case RenderType::k_capacity:
          values.push_back(capacity);
          break;
        case RenderType::k_weight:
          values.push_back(it->first);
          break;
        default:
          break;
        }
      }
      // step 2: normalize the values, and convert to colors
      for (auto it = values.begin(); it != values.end(); ++it) {
        if (*it > max_value) {
          *it = max_value;
        }
        if (*it < min_value) {
          *it = min_value;
        }
        *it = (*it - min_value) / (max_value - min_value);
        colors->push_back(get_rgb_color(h, *it, 1.0));
      }
    }

    Color get_vertex_sample_normalized_color(Vertex_const_handle vh,
        RenderType render_type, FT min_value, FT max_value, FT h) const {
      FT value = 0.0;
      switch (render_type) {
      case RenderType::k_gaussian_curvature:
        value = vh->gaussian_curvature();
        break;
      case RenderType::k_maximal_halfedge_dihedral:
        value = vh->max_halfedge_dihedral();
        break;
      case RenderType::k_ifi_facets:
      case RenderType::k_feature_intensity:
        value = vh->feature_intensity() - 1.0;
        break;
      case RenderType::k_capacity: {
        FT weight = vh->vertex_out_link().first;
        FT feature_intensity = vh->feature_intensity();
        value = weight / feature_intensity;
        break;
      }
      case RenderType::k_weight:
        value = vh->vertex_out_link().first;
        break;
      default:
        break;
      }
      if (value > max_value) {
        value = max_value;
      }
      if (value < min_value) {
        value = min_value;
      }
      value = (value - min_value) / (max_value - min_value);  // normalize
      return get_rgb_color(h, value, 1.0);
    }

    Color get_vertex_classification_color(bool inherit_element_types,
        FT feature_control_delta, Vertex_const_handle vh) const {
      Halfedge_const_list effective_edges;
      VertexType vertex_type = get_vertex_type(
        inherit_element_types, feature_control_delta, vh, &effective_edges);
      if (vertex_type == VertexType::k_feature_vertex) {
        return get_rgb_color(300, 0.8, 1.0);      // feature vertex: purple;
      }
      else if (vertex_type == VertexType::k_crease_vertex) {
        return get_rgb_color(180, 0.8, 1.0);      // crease vertex: cyan;
      }
      else {
        return get_rgb_color(60, 0.8, 1.0);       // smooth vertex: yellow.
      }
    }

    Color get_rgb_color(FT h, FT s, FT v) const {
      // h(1~360), s(0~1), v(0~1)
      int hi = h / 60;
      FT f = h / 60 - hi, p = v * (1.0 - s);
      FT q = v * (1.0 - f * s), t = v * (1.0 - (1 - f) * s);
      v *= 255, t *= 255, p *= 255, q *= 255;
      switch (hi) {
      case 0:
        return Color(v, t, p);
      case 1:
        return Color(q, v, p);
      case 2:
        return Color(p, v, t);
      case 3:
        return Color(p, q, v);
      case 4:
        return Color(t, p, v);
      case 5:
        return Color(v, p, q);
      default:
        return Color(0, 0, 0);
      }
    }

    void get_disturbed_border_samples_with_weights(Facet_const_handle fh,
        FT disturb_ratio, RenderType render_type,
        std::map<Point, double, Point_Comp> *disturbed_border_samples) const {
      // this version is used for rendering
      Halfedge_const_handle hh = fh->halfedge();
      FT value = 0.0;
      for (int i = 0; i <= 2; ++i) {
        // step 1: add the vertex sample
        Vertex_const_handle vh = get_opposite_vertex(hh);
        Vector vec = midpoint(hh) - vh->point();
        vec = vec * disturb_ratio;
        if (render_type == RenderType::k_capacity) {
          value = 1.0;
        }
        else if (render_type == RenderType::k_feature_intensity) {
          value = vh->vertex_out_link().first - 1.0;
        }
        else {
          value = vh->vertex_out_link().first;
        }
        disturbed_border_samples->insert(std::make_pair(
          vh->vertex_out_link().second.first + vec, value));
        // step 2: add the edge samples
        Halfedge_const_handle h = hh;
        if (hh->normal_dihedral() == -1.0) {
          h = h->opposite();
        }
        for (Link_list_const_iter it = h->edge_out_links().begin();
          it != h->edge_out_links().end(); ++it) {
          const Point &p = it->second.first;
          vec = vh->point() - p;
          vec = vec * disturb_ratio;
          if (render_type == RenderType::k_capacity) {
            value = 1.0;
          }
          else if (render_type == RenderType::k_feature_intensity) {
            value = it->first - 1.0;
          }
          else {
            value = it->first;
          }
          disturbed_border_samples->insert(std::make_pair(p + vec, value));
        }
        hh = hh->next();
      }
    }
  };
}

#endif  // ENRICHEDPOLYHEDRON_H_

/*void get_disturbed_border_samples(
  Facet_const_handle fh, FT disturb_ratio,
  std::set<Point, Point_Comp> *disturbed_border_samples) const {
  // this version is used for sampling

  Halfedge_const_handle hh = fh->halfedge();
  for (int i = 0; i <= 2; ++i) {
    Vertex_const_handle vh = get_opposite_vertex(hh);
    Vector vec = midpoint(hh) - vh->point();
    vec = vec * disturb_ratio;
    disturbed_border_samples->insert(vh->vertex_out_link().second.first + vec);
    Halfedge_const_handle h = hh;
    if (hh->normal_dihedral() == -1.0) {
      h = h->opposite();
    }
    for (Link_list_const_iter it = h->edge_out_links().begin();
      it != h->edge_out_links().end(); ++it) {
      const Point &p = it->second.first;
      vec = vh->point() - p;
      vec = vec * disturb_ratio;
      disturbed_border_samples->insert(p + vec);
    }
    hh = hh->next();
  }
}

void render_edge_voronoi_cell_boundaries(RenderType render_type) const {
  ::glBegin(GL_LINES);
  for (Facet_const_iterator fi = facets_begin();
    fi != facets_end(); ++fi) {
    const Point c = centroid(fi);
    Halfedge_const_handle hh = fi->halfedge();
    for (int i = 0; i <= 2; ++i) {
      // step 1: draw the edge feature intensity boundaries
      const Point &a = get_source_vertex(hh)->point();
      const Point &b = get_target_vertex(hh)->point();
      ::glVertex3d(a[0], a[1], a[2]);
      ::glVertex3d(c[0], c[1], c[2]);
      ::glVertex3d(b[0], b[1], b[2]);
      ::glVertex3d(c[0], c[1], c[2]);
      // step 2: draw the edge sample boundaries if necessary
      if (render_type != RenderType::k_normal_dihedral) {
        Halfedge_const_handle h = hh;
        if (h->normal_dihedral() == -1.0) {
          h = h->opposite();
        }
        Link_list_const_iter it1 = h->edge_out_links().begin();
        Link_list_const_iter it2 = it1;
        ++it2;
        for (; it2 != h->edge_out_links().end(); ++it1, ++it2) {
          const Point p = CGAL::midpoint(it1->second.first,
            it2->second.first);
          ::glVertex3d(p[0], p[1], p[2]);
          ::glVertex3d(c[0], c[1], c[2]);
        }
      }
      hh = hh->next();
    }
  }
  ::glEnd();
}

void render_facet_voronoi_cell_boundaries() const {
  for (Facet_const_iterator fi = facets_begin();
    fi != facets_end(); ++fi) {
    Point_list samples;
    for (Link_list_const_iter it = fi->facet_out_links().begin();
      it != fi->facet_out_links().end(); ++it) {
      samples.push_back(it->second.first);
    }
    Bvd bvd(triangle(fi));
    bvd.draw_voronoi_boundaries(samples);
  }
}

void render_all_samples() const {
  render_facet_samples();
  render_edge_samples();
  render_vertices();
}

void render_edges(RenderType render_type) const {
  ::glBegin(GL_LINES);
  for (Edge_const_iterator ei = edges_begin(); ei != edges_end(); ++ei) {
    Halfedge_const_handle hh = ei;
    if (hh->normal_dihedral() == -1.0) {
      hh = hh->opposite();
    }
    const Point &p = get_source_vertex(hh)->point();
    const Point &q = get_target_vertex(hh)->point();
    if (render_type == RenderType::k_classifications) {
      if (hh->is_crease()) {
        ::glColor3ub(200, 0, 0);
      }
      else {
        ::glColor3ub(0, 0, 200);
      }
    }
    else {
      ::glColor3ub(0, 0, 0);
    }
    ::glVertex3d(p[0], p[1], p[2]);
    ::glVertex3d(q[0], q[1], q[2]);
  }
  ::glEnd();
}

void render_all_voronoi(RenderType render_type,
  bool render_polyhedron_edges,
  FT sum_theta_value,
  FT dihedral_theta_value) const {
  // step 1: render all the sample cells
  ::glEnable(GL_LIGHTING);
  ::glEnable(GL_LIGHT0);
  ::glEnable(GL_POLYGON_OFFSET_FILL);
  ::glPolygonOffset(1.0f, 1.0f);
  render_all_voronoi_cells(render_type,
    sum_theta_value, dihedral_theta_value);
  ::glDisable(GL_POLYGON_OFFSET_FILL);
  ::glDisable(GL_LIGHTING);
  // step 2: render all the sample cell boundaries
  ::glColor3ub(0, 200, 0);
  ::glLineWidth(1.0f);
  render_all_voronoi_cell_boundaries();
  // step 3: render the edges
  if (render_polyhedron_edges) {
    ::glLineWidth(2.0f);
    render_edges(render_type);
  }
  // step 4: render the sample points
  ::glColor3ub(0, 200, 0);
  ::glPointSize(3.0f);
  render_all_samples();
}

void render_edge_voronoi_cells(RenderType render_type,
  FT sum_theta_value,
  FT dihedral_theta_value) const {
  // step 1: get the minimal and maximal values for color
  FT min_value = DOUBLE_MAX, max_value = DOUBLE_MIN;
  if (render_type == RenderType::k_normal_dihedral) {
    min_value = 0.0;
    max_value = dihedral_theta_value;
  }
  else if (render_type == RenderType::k_feature_intensity) {
    min_value = 0.0;
    max_value = (sum_theta_value + 1) * (dihedral_theta_value + 1) - 1;
  }
  else {  // RenderType::k_capacity or RenderType::k_weight
    min_value = 0.0;
    for (Edge_const_iterator ei = edges_begin(); ei != edges_end(); ++ei) {
      Halfedge_const_handle hh = ei;
      if (hh->normal_dihedral() == -1.0) {
        hh = hh->opposite();
      }
      FT facet_area = area(hh->facet());
      if (!hh->opposite()->is_border()) {
        facet_area += area(hh->opposite()->facet());
      }
      FT capacity = facet_area / (3 * hh->edge_out_links().size());
      if (render_type == RenderType::k_capacity) {
        max_value = CGAL::max(max_value, 2 * capacity);
      }
      else {             // RenderType::k_weight
        for (Link_list_const_iter it = hh->edge_out_links().begin();
          it != hh->edge_out_links().end(); ++it) {
          max_value = CGAL::max(max_value, it->first);
        }
      }
    }
  }
  // step 2: render the voronoi cells
  Color color;
  ::glBegin(GL_TRIANGLES);
  for (Edge_const_iterator ei = edges_begin(); ei != edges_end(); ++ei) {
    Halfedge_const_handle hh = ei;
    if (hh->normal_dihedral() == -1.0) {
      hh = hh->opposite();
    }
    if (render_type == RenderType::k_normal_dihedral) {
      // for edge feanture intensity, we only need to render two triangles
      FT value =
        (hh->normal_dihedral() - min_value) / (max_value - min_value);
      color = get_rgb_color(240, value, 1.0);
      Facet_const_handle fh = hh->facet();  // the first triangle
      const Normal &n = fh->normal();
      const Point &a = get_target_vertex(hh)->point();
      const Point b = centroid(fh);
      const Point &c = get_source_vertex(hh)->point();
      ::glColor3ub(color.r(), color.g(), color.b());
      ::glNormal3d(n[0], n[1], n[2]);
      ::glVertex3d(a[0], a[1], a[2]);
      ::glVertex3d(b[0], b[1], b[2]);
      ::glVertex3d(c[0], c[1], c[2]);
      if (!hh->opposite()->is_border()) {   // the second triangle if exists
        hh = hh->opposite();
        Facet_const_handle fh = hh->facet();
        const Normal &n = fh->normal();
        const Point &a = get_target_vertex(hh)->point();
        const Point b = centroid(fh);
        const Point &c = get_source_vertex(hh)->point();
        ::glColor3ub(color.r(), color.g(), color.b());
        ::glNormal3d(n[0], n[1], n[2]);
        ::glVertex3d(a[0], a[1], a[2]);
        ::glVertex3d(b[0], b[1], b[2]);
        ::glVertex3d(c[0], c[1], c[2]);
      }
    }
    else {  // for others, we need to render for each edge sample
      // step 2.1: get the midpoints of neighboring edge samples
      Link_list_const_iter first = hh->edge_out_links().begin();
      Link_list_const_iter second = first;
      ++second;
      Point_list points;
      points.push_back(get_source_vertex(hh)->point());
      for (; second != hh->edge_out_links().end(); ++first, ++second) {
        points.push_back(CGAL::midpoint(first->second.first,
          second->second.second));
      }
      points.push_back(get_target_vertex(hh)->point());
      // step 2.2: get the normalized values
      Color_list colors;
      get_edge_normalized_render_colors(hh, render_type, min_value,
        max_value, 240, &colors);
      Facet_const_handle fh = hh->facet();
      Point c = centroid(fh);
      const Normal &n = fh->normal();
      Point_const_iter it1 = points.begin();
      Point_const_iter it2 = it1;
      ++it2;
      Color_const_iter cit = colors.begin();
      for (; it2 != points.end(); ++it1, ++it2, ++cit) {
        ::glColor3ub(cit->r(), cit->g(), cit->b());
        ::glNormal3d(n[0], n[1], n[2]);
        ::glVertex3d((*it1)[0], (*it1)[1], (*it1)[2]);
        ::glVertex3d((*it2)[0], (*it2)[1], (*it2)[2]);
        ::glVertex3d(c[0], c[1], c[2]);
      }
      if (!hh->opposite()->is_border()) {
        fh = hh->opposite()->facet();
        c = centroid(fh);
        const Normal &n = fh->normal();
        it1 = points.begin();
        it2 = it1;
        ++it2;
        cit = colors.begin();
        for (; it2 != points.end(); ++it1, ++it2, ++cit) {
          ::glColor3ub(cit->r(), cit->g(), cit->b());
          ::glNormal3d(n[0], n[1], n[2]);
          ::glVertex3d((*it1)[0], (*it1)[1], (*it1)[2]);
          ::glVertex3d((*it2)[0], (*it2)[1], (*it2)[2]);
          ::glVertex3d(c[0], c[1], c[2]);
        }
      }
    }
  }
  ::glEnd();
}

void render_edge_voronoi(RenderType render_type,
  bool render_polyhedron_edges,
  FT sum_theta_value,
  FT dihedral_theta_value) const {
  // step 1: render the edge (samples) voronoi cells
  ::glEnable(GL_LIGHTING);
  ::glEnable(GL_LIGHT0);
  ::glEnable(GL_POLYGON_OFFSET_FILL);
  ::glPolygonOffset(1.0f, 1.0f);
  render_edge_voronoi_cells(render_type,
    sum_theta_value, dihedral_theta_value);
  ::glDisable(GL_POLYGON_OFFSET_FILL);
  ::glDisable(GL_LIGHTING);
  // step 2: render the edge voronoi boundaries
  ::glColor3ub(0, 200, 0);
  ::glLineWidth(1.0f);
  render_edge_voronoi_cell_boundaries(render_type);
  // step 4: render the edges
  if (render_polyhedron_edges) {
    ::glColor3ub(0, 0, 0);
    render_edges(render_type);
  }
  // step 3: render the samples
  ::glColor3ub(0, 200, 0);
  ::glPointSize(3.0);
  render_edge_samples();
}

void render_vertex_voronoi_cells(RenderType render_type,
  bool inherit_element_types,
  FT feature_control_delta,
  FT sum_theta_value,
  FT dihedral_theta_value) {
  // step 1: get the minimal and maximal values for color
  FT min_value = 0.0, max_value = DOUBLE_MIN;
  if (render_type == RenderType::k_feature_intensity) {
    max_value = (sum_theta_value + 1) * (dihedral_theta_value + 1) - 1;
  }
  else if (render_type == RenderType::k_gaussian_curvature) {
    max_value = sum_theta_value;
  }
  else if (render_type == RenderType::k_maximal_halfedge_dihedral) {
    max_value = dihedral_theta_value;
  }
  else {
    FT weight = 0.0, feature_intensity = 0.0;
    for (Vertex_const_iterator vi = vertices_begin();
      vi != vertices_end(); ++vi) {
      switch (render_type) {
      case RenderType::k_capacity:
        weight = vi->vertex_out_link().first;
        feature_intensity = vi->feature_intensity();
        max_value = CGAL::max(max_value, 2 * weight / feature_intensity);
        break;
      case RenderType::k_weight:
        max_value = CGAL::max(max_value, vi->vertex_out_link().first);
        break;
      default:
        break;
      }
    }
  }
  // step 2: render the voronoi cells
  Color color;
  ::glBegin(GL_TRIANGLES);
  for (Vertex_const_iterator vi = vertices_begin();
      vi != vertices_end(); ++vi) {
    const Point &p = vi->point();
    Halfedge_around_vertex_const_circulator vcirc = vi->vertex_begin();
    Halfedge_around_vertex_const_circulator vend = vcirc;
    if (render_type == RenderType::k_classifications) {
      color = get_vertex_classification_color(inherit_element_types,
        feature_control_delta, vi);
    }
    else {
      color = get_vertex_normalized_render_color(vi, render_type, min_value,
        max_value, 240);
    }
    ::glColor3ub(color.r(), color.g(), color.b());
    CGAL_For_all(vcirc, vend) {
      if (!vcirc->is_border()) {
        Facet_const_handle fh = vcirc->facet();
        const Point p1 = midpoint(vcirc);
        const Point c = centroid(fh);
        const Point p2 = midpoint(vcirc->next());
        const Normal &n = fh->normal();
        ::glNormal3d(n[0], n[1], n[2]);
        ::glVertex3d(p1[0], p1[1], p1[2]);	// the first triangle
        ::glVertex3d(p[0], p[1], p[2]);
        ::glVertex3d(c[0], c[1], c[2]);
        ::glVertex3d(c[0], c[1], c[2]);		  // the second triangle
        ::glVertex3d(p[0], p[1], p[2]);
        ::glVertex3d(p2[0], p2[1], p2[2]);
      }
    }
  }
  ::glEnd();
}

void render_vertex_voronoi_cell_boundaries() const {
  for (Vertex_const_iterator vi = vertices_begin();
    vi != vertices_end(); ++vi) {
    const Point &p = vi->point();
    Halfedge_around_vertex_const_circulator vcirc = vi->vertex_begin();
    Halfedge_around_vertex_const_circulator vend = vcirc;
    ::glBegin(GL_LINE_STRIP);
    const Point p1 = midpoint(vcirc);
    ::glVertex3d(p1[0], p1[1], p1[2]);
    CGAL_For_all(vcirc, vend) {
      if (!vcirc->is_border()) {
        Facet_const_handle fh = vcirc->facet();
        const Point p2 = centroid(fh);
        const Point p3 = midpoint(vcirc->next());
        ::glVertex3d(p2[0], p2[1], p2[2]);
        ::glVertex3d(p3[0], p3[1], p3[2]);
      }
      else {
        const Point p4 = midpoint(vcirc->next());
        ::glVertex3d(p[0], p[1], p[2]);
        ::glVertex3d(p4[0], p4[1], p4[2]);
      }
    }
    ::glEnd();
  }
}

void render_vertex_voronoi(RenderType render_type,
  bool inherit_element_types,
  bool render_polyhedron_edges,
  FT feature_control_delta,
  FT sum_theta_value,
  FT dihedral_theta_value) {
  // step 1: render the voronoi (sample) cells
  ::glEnable(GL_POLYGON_OFFSET_FILL);
  ::glPolygonOffset(1.0f, 1.0f);
  ::glEnable(GL_LIGHTING);
  ::glEnable(GL_LIGHT0);
  render_vertex_voronoi_cells(render_type, inherit_element_types,
    feature_control_delta, sum_theta_value, dihedral_theta_value);
  ::glDisable(GL_LIGHTING);
  ::glDisable(GL_POLYGON_OFFSET_FILL);
  ::glLineWidth(1.0f);
  if (render_type != RenderType::k_normal_dihedral &&
    render_type != RenderType::k_classifications) {
    // step 2: render the voronoi cell boundaries
    ::glColor3ub(0, 200, 0);
    render_vertex_voronoi_cell_boundaries();
    // step 3: render the voronoi cell points
    ::glPointSize(3.0f);
    render_vertices();
  }
  // step 3: render the edges
  if (render_polyhedron_edges) {
    render_edges(render_type);
  }
}

void render_facet_voronoi(RenderType render_type,
  bool render_polyhedron_edges,
  FT sum_theta_value,
  FT dihedral_theta_value) const {
  // step 1: render the facet sample cells
  ::glEnable(GL_LIGHTING);
  ::glEnable(GL_LIGHT0);
  ::glEnable(GL_POLYGON_OFFSET_FILL);
  ::glPolygonOffset(1.0f, 1.0f);
  render_facet_voronoi_cells(render_type,
    sum_theta_value, dihedral_theta_value);
  ::glDisable(GL_POLYGON_OFFSET_FILL);
  ::glDisable(GL_LIGHTING);
  // step 2: render the facet sample cell boundaries
  ::glColor3ub(0, 200, 0);
  ::glLineWidth(1.0f);
  render_facet_voronoi_cell_boundaries();
  // step 3: render the edges
  if (render_polyhedron_edges) {
    ::glLineWidth(2.0f);
    render_edges(render_type);
  }
  // step 4: render the facet samples
  ::glColor3ub(0, 200, 0);
  ::glPointSize(3.0f);
  render_facet_samples();
}

void render_facet_voronoi_cells(RenderType render_type,
  FT sum_theta_value,
  FT dihedral_theta_value) const {
  // step 1: get the minimal and maximal values for color
  FT min_value = 0.0, max_value = DOUBLE_MIN;
  if (render_type == RenderType::k_feature_intensity) {
    max_value = (sum_theta_value + 1) * (dihedral_theta_value + 1) - 1;
  }
  else {  // RenderType::k_capacity or RenderType::k_weight
    for (Facet_const_iterator fi = facets_begin();
      fi != facets_end(); ++fi) {
      FT capacity = area(fi) / fi->facet_out_links().size();
      if (render_type == RenderType::k_capacity) {
        max_value = CGAL::max(max_value, 2 * capacity);
      }
      else {
        for (Link_list_const_iter it = fi->facet_out_links().begin();
          it != fi->facet_out_links().end(); ++it) {
          max_value = CGAL::max(max_value, it->first);
        }
      }
    }
  }
  // step 2: render the facet voronoi cells
  Color color;
  for (Facet_const_iterator fi = facets_begin();
    fi != facets_end(); ++fi) {
    const Normal n = fi->normal();
    Point_list samples;
    Color_list colors;
    get_facet_normalized_render_colors(fi, render_type, min_value,
      max_value, 240, &samples, &colors);
    ::glNormal3d(n[0], n[1], n[2]);
    Bvd bvd(triangle(fi));
    bvd.draw_voronoi_cells(samples, colors);
  }
}

void render_vertices() const {
  ::glBegin(GL_POINTS);
  for (Vertex_const_iterator vi = vertices_begin();
    vi != vertices_end(); ++vi) {
    const Point &p = vi->point();
    ::glVertex3d(p[0], p[1], p[2]);
  }
  ::glEnd();
}

void render_edge_samples() const {
  ::glBegin(GL_POINTS);
  for (Edge_const_iterator ei = edges_begin(); ei != edges_end(); ++ei) {
    Halfedge_const_handle hh = ei;
    if (hh->normal_dihedral() == -1.0) {
      hh = hh->opposite();
    }
    for (Link_list_const_iter it = hh->edge_out_links().begin();
      it != hh->edge_out_links().end(); ++it) {
      const Point &p = it->second.first;
      ::glVertex3d(p[0], p[1], p[2]);
    }
  }
  ::glEnd();
}

void render_facet_samples() const {
  ::glBegin(GL_POINTS);
  for (Facet_const_iterator fi = facets_begin();
    fi != facets_end(); ++fi) {
    for (Link_list_const_iter it = fi->facet_out_links().begin();
      it != fi->facet_out_links().end(); ++it) {
      const Point &p = it->second.first;
      ::glVertex3d(p[0], p[1], p[2]);
    }
  }
  ::glEnd();
}

void render_all_voronoi_cells(RenderType render_type,
  FT sum_theta_value,
  FT dihedral_theta_value) const {
  // step 1: get the minimal and maximal values for color
  FT min_value = 0.0, max_value = DOUBLE_MIN;
  if (render_type == RenderType::k_feature_intensity) {
    max_value = (sum_theta_value + 1) * (dihedral_theta_value + 1) - 1;
  }
  else if (render_type == RenderType::k_capacity) {
    max_value = 2.0;
  }
  else {                  // RenderType::k_weight
    // the maximial value will certainly be in vertex sample
    for (Vertex_const_iterator vi = vertices_begin();
      vi != vertices_end(); ++vi) {
      max_value = CGAL::max(max_value, vi->vertex_out_link().first);
    }
  }
  // step 2: render all the sample cells
  Color color;
  for (Facet_const_iterator fi = facets_begin();
    fi != facets_end(); ++fi) {
    const Normal n = fi->normal();
    Point_list samples;
    Color_list colors;
    get_all_normalized_render_colors(fi, render_type, min_value, max_value,
      240, &samples, &colors);
    ::glNormal3d(n[0], n[1], n[2]);
    Bvd bvd(triangle(fi));
    bvd.draw_voronoi_cells(samples, colors);
  }
}

void render_all_voronoi_cell_boundaries() const {
  for (Facet_const_iterator fi = facets_begin();
    fi != facets_end(); ++fi) {
    // facet samples
    Point_list samples;
    for (Link_list_const_iter it = fi->facet_out_links().begin();
      it != fi->facet_out_links().end(); ++it) {
      samples.push_back(it->second.first);
    }
    // edge samples and vertex samples
    std::set<Point, Point_Comp> disturbed_border_samples;
    get_disturbed_border_samples(fi, 0.01, &disturbed_border_samples);
    for (auto it = disturbed_border_samples.begin();
      it != disturbed_border_samples.end(); ++it) {
      samples.push_back(*it);
    }
    Bvd bvd(triangle(fi));
    bvd.draw_voronoi_boundaries(samples);
  }
}

VertexType calculate_vertex_type(Vertex_const_handle vh,
  FT feature_control_delta,
  Halfedge_const_list *effective_edges) const {
  // effective_edges will contain all the one-ring effective edges
  // if effective_edges is empty, vh cannot be relocated
  //    (vh might be a feature_vertex or an end of a crease)
  // if effective_edges has 2 elements, vh can be relocated along crease
  //    (vh should be a crease_vertex)
  // if effective_edges has more elements, vh can be relocated freely
  //    (vh should be a smooth_vertex, or might be a boundary vertex)
  effective_edges->clear();
  std::multimap<FT, Halfedge_const_handle> edge_map;
  Halfedge_around_vertex_const_circulator vcirc = vh->vertex_begin();
  Halfedge_around_vertex_const_circulator vend = vcirc;
  int effective_vertex_count = 0;
  int degree = static_cast<int>(vh->degree());
  FT f1 = vh->feature_intensity() * feature_control_delta;
  FT f2 = (vh->max_halfedge_dihedral() + 1) * feature_control_delta;
  CGAL_For_all(vcirc, vend) {
    Vertex_const_handle v = get_source_vertex(vcirc);
    effective_edges->push_back(vcirc);
    bool b1 = v->feature_intensity() >= f1;
    FT normal_dihedral = vcirc->normal_dihedral() == -1.0 ?
      vcirc->opposite()->normal_dihedral() : vcirc->normal_dihedral();
    bool b2 = (normal_dihedral + 1) >= f2;
    if (b1 && b2) {
      ++effective_vertex_count;
      edge_map.insert(std::make_pair(v->feature_intensity(), vcirc));
    }
  }
  if (effective_vertex_count == 0) {
    effective_edges->clear();
    return VertexType::k_feature_vertex;
  }
  else if (effective_vertex_count >= 3) {
    return VertexType::k_smooth_vertex;
  }
  else {
    effective_edges->clear();
    std::multimap<FT, Halfedge_const_handle>::reverse_iterator it;
    it = edge_map.rbegin();
    for (int i = 0; i < effective_vertex_count; ++i) {
      effective_edges->push_back(it->second);
      ++it;
    }
    return VertexType::k_crease_vertex;
  }
}*/