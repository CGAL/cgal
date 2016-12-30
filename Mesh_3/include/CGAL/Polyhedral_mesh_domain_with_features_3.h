// Copyright (c) 2009-2010 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Stéphane Tayeb
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#ifndef CGAL_POLYHEDRAL_MESH_DOMAIN_WITH_FEATURES_3_H
#define CGAL_POLYHEDRAL_MESH_DOMAIN_WITH_FEATURES_3_H

#include <CGAL/Mesh_3/config.h>

#include <CGAL/Random.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/Mesh_domain_with_polyline_features_3.h>
#include <CGAL/Mesh_polyhedron_3.h>

#include <CGAL/Mesh_3/Detect_polylines_in_polyhedra.h>
#include <CGAL/Mesh_3/Polyline_with_context.h>
#include <CGAL/Mesh_3/Detect_features_in_polyhedra.h>

#include <CGAL/enum.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/helpers.h>
#include <boost/graph/filtered_graph.hpp>

#include <boost/iterator/transform_iterator.hpp>
#include <boost/foreach.hpp>

#include <string>
#include <vector>
#include <fstream>


namespace CGAL {

namespace internal {
namespace Mesh_3 {

template <typename Kernel>
struct Angle_tester
{
  template <typename vertex_descriptor, typename Graph>
  bool operator()(vertex_descriptor& v, const Graph& g) const
  {
    typedef typename boost::graph_traits<Graph>::out_edge_iterator out_edge_iterator;
    if (out_degree(v, g) != 2)
      return true;
    else
    {
      out_edge_iterator out_edge_it, out_edges_end;
      boost::tie(out_edge_it, out_edges_end) = out_edges(v, g);

      vertex_descriptor v1 = target(*out_edge_it++, g);
      vertex_descriptor v2 = target(*out_edge_it++, g);
      CGAL_assertion(out_edge_it == out_edges_end);

      const typename Kernel::Point_3& p = v->point();
      const typename Kernel::Point_3& p1 = v1->point();
      const typename Kernel::Point_3& p2 = v2->point();

      return (CGAL::angle(p1, p, p2) == CGAL::ACUTE);
    }
  }
};

template <typename Polyhedron>
struct Is_featured_edge {
  const Polyhedron* polyhedron;
  Is_featured_edge() : polyhedron(0) {} // required by boost::filtered_graph
  Is_featured_edge(const Polyhedron& polyhedron) : polyhedron(&polyhedron) {}

  bool operator()(typename boost::graph_traits<Polyhedron>::edge_descriptor e) const {
    return halfedge(e, *polyhedron)->is_feature_edge();
  }
}; // end Is_featured_edge<Polyhedron>

template <typename Polyhedron>
struct Is_border_edge {
  const Polyhedron* polyhedron;
  Is_border_edge() : polyhedron(0) {} // required by boost::filtered_graph
  Is_border_edge(const Polyhedron& polyhedron) : polyhedron(&polyhedron) {}

  bool operator()(typename boost::graph_traits<Polyhedron>::edge_descriptor e) const {
    return is_border(halfedge(e, *polyhedron), *polyhedron) ||
      is_border(opposite(halfedge(e, *polyhedron), *polyhedron), *polyhedron);
  }
}; // end Is_featured_edge<Polyhedron>

template<typename Polyhedral_mesh_domain, typename Polyline_with_context>
struct Extract_polyline_with_context_visitor
{
  typedef typename Polyhedral_mesh_domain::Polyhedron Polyhedron;
  std::vector<Polyline_with_context>& polylines;
  const Polyhedron& polyhedron;

  Extract_polyline_with_context_visitor
  (const Polyhedron& polyhedron,
   typename std::vector<Polyline_with_context>& polylines)
    : polylines(polylines), polyhedron(polyhedron)
  {}

  void start_new_polyline()
  {
    polylines.push_back(Polyline_with_context());
  }

  void add_node(typename Polyhedron::Vertex_handle vd)
  {
    if(polylines.back().polyline_content.empty()) {
      polylines.back().polyline_content.push_back(vd->point());
    }
  }

  void add_edge(typename boost::graph_traits<Polyhedron>::edge_descriptor ed)
  {
    Polyline_with_context& polyline = polylines.back();
    typename Polyhedron::Halfedge_handle he = halfedge(ed, polyhedron);
    CGAL_assertion(!polyline.polyline_content.empty());
    if(polyline.polyline_content.back() != he->vertex()->point()) {
      polyline.polyline_content.push_back(he->vertex()->point());
    } else if(polyline.polyline_content.back() !=
         he->opposite()->vertex()->point())
    { // if the edge is zero-length, it is ignored
        polyline.polyline_content.push_back(he->opposite()->vertex()->point());
    }
    typename Polyhedral_mesh_domain::Surface_patch_index_generator generator;
    if(!is_border(he, polyhedron)) {
      polyline.context.adjacent_patches_ids.insert(generator(he->face()));
    }
    he = he->opposite();
    if(!is_border(he, polyhedron)) {
      polyline.context.adjacent_patches_ids.insert(generator(he->face()));
    }
  }

  void end_polyline()
  {
    // ignore degenerated polylines
    if(polylines.back().polyline_content.size() < 2)
      polylines.resize(polylines.size() - 1);
    // else {
    //   std::cerr << "Polyline with " << polylines.back().polyline_content.size()
    //             << " vertices, incident to "
    //             << polylines.back().context.adjacent_patches_ids.size()
    //             << " patches:\n ";
    //   for(auto p: polylines.back().polyline_content)
    //     std::cerr << " " << p;
    //   std::cerr << "\n";
    // }
  }
};


} // end CGAL::internal::Mesh_3
} // end CGAL::internal

/**
 * @class Polyhedral_mesh_domain_with_features_3
 *
 *
 */
template < class IGT_,
           class Polyhedron_ = typename Mesh_polyhedron_3<IGT_>::type,
           class TriangleAccessor=Triangle_accessor_3<Polyhedron_,IGT_>,
           class Use_patch_id_tag = Tag_true,
           class Use_exact_intersection_construction_tag = Tag_true >
class Polyhedral_mesh_domain_with_features_3
  : public Mesh_domain_with_polyline_features_3<
      Polyhedral_mesh_domain_3< Polyhedron_,
                                IGT_,
                                TriangleAccessor,
                                Use_patch_id_tag,
                                Use_exact_intersection_construction_tag > >
{
  typedef Mesh_domain_with_polyline_features_3<
    Polyhedral_mesh_domain_3<
      Polyhedron_, IGT_, TriangleAccessor,
      Use_patch_id_tag, Use_exact_intersection_construction_tag > > Base;

public:
  typedef Polyhedron_ Polyhedron;

  // Index types
  typedef typename Base::Index                Index;
  typedef typename Base::Corner_index         Corner_index;
  typedef typename Base::Curve_segment_index  Curve_segment_index;
  typedef typename Base::Surface_patch_index  Surface_patch_index;
  typedef typename Base::Subdomain_index      Subdomain_index;
  
  // Backward compatibility
#ifndef CGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX
  typedef Surface_patch_index                 Surface_index;
#endif // CGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX

  typedef typename Base::R         R;
  typedef typename Base::Point_3   Point_3;
  typedef typename Base::FT        FT;
  
  typedef CGAL::Tag_true           Has_features;

  typedef std::vector<Point_3> Bare_polyline;
  typedef Mesh_3::Polyline_with_context<Surface_patch_index, Curve_segment_index,
                                        Bare_polyline > Polyline_with_context;
  /// Constructors
  Polyhedral_mesh_domain_with_features_3(const Polyhedron& p,
    CGAL::Random* p_rng = NULL);

  CGAL_DEPRECATED Polyhedral_mesh_domain_with_features_3(const std::string& filename,
    CGAL::Random* p_rng = NULL);

  // The following is needed, because otherwise, when a "const char*" is
  // passed, the constructors templates are a better match, than the
  // constructor with `std::string`.
  CGAL_DEPRECATED Polyhedral_mesh_domain_with_features_3(const char* filename,
    CGAL::Random* p_rng = NULL);

  // Inherited constructors
  template <typename T2>
  Polyhedral_mesh_domain_with_features_3(const Polyhedron& p,
                                         const T2& b,
                                         CGAL::Random* p_rng = NULL)
    : Base(p, b)
    , borders_detected_(false)
  {
    poly_.resize(1);
    poly_[0] = p;
    this->set_random_generator(p_rng);
  }

  template <typename InputPolyhedraPtrIterator>
  Polyhedral_mesh_domain_with_features_3(InputPolyhedraPtrIterator begin,
                                         InputPolyhedraPtrIterator end,
                                         CGAL::Random* p_rng = NULL)
    : Base(begin, end)
    , borders_detected_(false)
  {
    if (begin != end) {
      for (; begin != end; ++begin)
        poly_.push_back(**begin);
    }
    else
      poly_.push_back(**begin);
    this->set_random_generator(p_rng);
  }

  template <typename InputPolyhedraPtrIterator>
  Polyhedral_mesh_domain_with_features_3(InputPolyhedraPtrIterator begin,
                                         InputPolyhedraPtrIterator end,
                                         const Polyhedron& bounding_polyhedron,
                                         CGAL::Random* p_rng = NULL)
    : Base(begin, end, bounding_polyhedron)
    , borders_detected_(false)
  {
    if (begin != end) {
      for (; begin != end; ++begin)
        poly_.push_back(**begin);
    }
    else
      poly_.push_back(**begin);
    this->set_random_generator(p_rng);
  }

  /// Destructor
  ~Polyhedral_mesh_domain_with_features_3() {}

  /// Detect features
  void initialize_ts(Polyhedron& p);

  void detect_features(FT angle_in_degree, std::vector<Polyhedron>& p);
  void detect_features(FT angle_in_degree = FT(60)) { detect_features(angle_in_degree, poly_); }

  void detect_borders(const std::vector<Polyhedron>& p);
  void detect_borders() { detect_borders(poly_); };

private:
  template <typename Edge_predicate>
  void add_features_from_split_graph_into_polylines(const Polyhedron& poly,
                                                    const Edge_predicate& pred);

  std::vector<Polyhedron> poly_;
  bool borders_detected_;

private:
  // Disabled copy constructor & assignment operator
  typedef Polyhedral_mesh_domain_with_features_3 Self;
  Polyhedral_mesh_domain_with_features_3(const Self& src);
  Self& operator=(const Self& src);

};  // end class Polyhedral_mesh_domain_with_features_3


template < typename GT_, typename P_, typename TA_,
           typename Tag_, typename E_tag_>
Polyhedral_mesh_domain_with_features_3<GT_,P_,TA_,Tag_,E_tag_>::
Polyhedral_mesh_domain_with_features_3(const Polyhedron& p,
                                       CGAL::Random* p_rng)
  : Base()
  , borders_detected_(false)
{
  poly_.resize(1);
  poly_[0] = p;
  this->add_primitives(poly_[0]);
  this->set_random_generator(p_rng);
}


template < typename GT_, typename P_, typename TA_,
           typename Tag_, typename E_tag_>
CGAL_DEPRECATED Polyhedral_mesh_domain_with_features_3<GT_,P_,TA_,Tag_,E_tag_>::
Polyhedral_mesh_domain_with_features_3(const char* filename,
                                       CGAL::Random* p_rng)
  : Base()
  , borders_detected_(false)
{
  // Create input polyhedron
  std::ifstream input(filename);
  poly_.resize(1);
  input >> poly_[0];
  this->add_primitives(poly_[0]);
  this->set_random_generator(p_rng);
}


template < typename GT_, typename P_, typename TA_,
           typename Tag_, typename E_tag_>
CGAL_DEPRECATED Polyhedral_mesh_domain_with_features_3<GT_,P_,TA_,Tag_,E_tag_>::
Polyhedral_mesh_domain_with_features_3(const std::string& filename,
                                       CGAL::Random* p_rng)
  : Base()
  , borders_detected_(false)
{
  // Create input polyhedron
  std::ifstream input(filename.c_str());
  poly_.resize(1);
  input >> poly_[0];
  this->add_primitives(poly_[0]);
  this->set_random_generator(p_rng);
}


template < typename GT_, typename P_, typename TA_,
           typename Tag_, typename E_tag_>
void
Polyhedral_mesh_domain_with_features_3<GT_,P_,TA_,Tag_,E_tag_>::
initialize_ts(Polyhedron& p)
{
  std::size_t ts = 0;
  for(typename Polyhedron::Vertex_iterator v = p.vertices_begin(),
      end = p.vertices_end() ; v != end ; ++v)
  {
    v->set_time_stamp(ts++);
  }
  for(typename Polyhedron::Facet_iterator fit = p.facets_begin(),
       end = p.facets_end() ; fit != end ; ++fit )
  {
    fit->set_time_stamp(ts++);
  }
  for(typename Polyhedron::Halfedge_iterator hit = p.halfedges_begin(),
       end = p.halfedges_end() ; hit != end ; ++hit )
  {
    hit->set_time_stamp(ts++);
  }
}




template <typename Graph>
void dump_graph_edges(std::ostream& out, const Graph& g)
{
  typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<Graph>::edge_descriptor edge_descriptor;

  out.precision(17);
  BOOST_FOREACH(edge_descriptor e, edges(g))
  {
    vertex_descriptor s = source(e, g);
    vertex_descriptor t = target(e, g);
    out << "2 " << s->point() << " " << t->point() << "\n";
  }
}

template <typename Graph>
void dump_graph_edges(const char* filename, const Graph& g)
{
  std::ofstream out(filename);
  dump_graph_edges(out, g);
}

template < typename GT_, typename P_, typename TA_,
           typename Tag_, typename E_tag_>
void
Polyhedral_mesh_domain_with_features_3<GT_,P_,TA_,Tag_,E_tag_>::
detect_features(FT angle_in_degree, std::vector<Polyhedron>& poly)
{
  CGAL_assertion(!borders_detected_);
  if (borders_detected_)
    return;//prevent from not-terminating

  BOOST_FOREACH(Polyhedron& p, poly)
  {
    initialize_ts(p);

    // Get sharp features
    Mesh_3::detect_features(p,angle_in_degree);

    internal::Mesh_3::Is_featured_edge<Polyhedron> is_featured_edge(p);

    add_features_from_split_graph_into_polylines(p, is_featured_edge);
  }
  borders_detected_ = true;/*done by Mesh_3::detect_features*/
}

template < typename GT_, typename P_, typename TA_,
           typename Tag_, typename E_tag_>
template <typename Edge_predicate>
void
Polyhedral_mesh_domain_with_features_3<GT_,P_,TA_,Tag_,E_tag_>::
add_features_from_split_graph_into_polylines(const Polyhedron& p,
                                             const Edge_predicate& pred)
{
  typedef boost::filtered_graph<Polyhedron,
                                Edge_predicate > Featured_edges_graph;
  Featured_edges_graph graph(p, pred);

  std::vector<Polyline_with_context> polylines;

  internal::Mesh_3::Extract_polyline_with_context_visitor<
    Polyhedral_mesh_domain_with_features_3,
    Polyline_with_context
    > visitor(p, polylines);
  internal::Mesh_3::Angle_tester<GT_> angle_tester;
  split_graph_into_polylines(graph, visitor, angle_tester);

  this->add_features_with_context(polylines.begin(),
                                  polylines.end());

}


template < typename GT_, typename P_, typename TA_,
           typename Tag_, typename E_tag_>
void
Polyhedral_mesh_domain_with_features_3<GT_,P_,TA_,Tag_,E_tag_>::
detect_borders(const std::vector<Polyhedron>& poly)
{
  if (borders_detected_)
    return;//border detection has already been done

  BOOST_FOREACH(const Polyhedron& p, poly)
  {
    internal::Mesh_3::Is_border_edge<Polyhedron> is_border_edge(p);
    add_features_from_split_graph_into_polylines(p, is_border_edge);
  }

  borders_detected_ = true;
}


} //namespace CGAL


#endif // CGAL_POLYHEDRAL_MESH_DOMAIN_WITH_FEATURES_3_H
