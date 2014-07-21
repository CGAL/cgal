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

#include <boost/iterator/transform_iterator.hpp>

#include <string>
#include <vector>
#include <fstream>


namespace CGAL {

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
  
  typedef Polyhedron_ Polyhedron;
  
public:
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

  /// Constructors
  Polyhedral_mesh_domain_with_features_3(const Polyhedron& p,
    CGAL::Random* p_rng = NULL);
  Polyhedral_mesh_domain_with_features_3(const std::string& filename,
    CGAL::Random* p_rng = NULL);

  // The following is needed, because otherwise, when a "const char*" is
  // passed, the constructors templates are a better match, than the
  // constructor with `std::string`.
  Polyhedral_mesh_domain_with_features_3(const char* filename,
    CGAL::Random* p_rng = NULL);

  template <typename T1, typename T2>
  Polyhedral_mesh_domain_with_features_3(const T1& a, const T2& b,
                                         CGAL::Random* p_rng = NULL)
    : Base(a, b)
  { this->set_random_generator(p_rng); }


  template <typename T1, typename T2, typename T3>
  Polyhedral_mesh_domain_with_features_3(const T1& a, const T2& b, const T3& c,
                                         CGAL::Random* p_rng = NULL)
    : Base(a, b, c)
  { this->set_random_generator(p_rng); }

  /// Destructor
  ~Polyhedral_mesh_domain_with_features_3() {}

  /// Detect features
  void initialize_ts(Polyhedron& p);

  void detect_features(FT angle_in_degree, Polyhedron& p);
  void detect_features(FT angle_in_degree = FT(60)) { detect_features(angle_in_degree, polyhedron_); }

private:
  Polyhedron polyhedron_;

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
  , polyhedron_(p)
{
  this->add_primitives(polyhedron_);
  this->set_random_generator(p_rng);
}

template < typename GT_, typename P_, typename TA_,
           typename Tag_, typename E_tag_>
Polyhedral_mesh_domain_with_features_3<GT_,P_,TA_,Tag_,E_tag_>::
Polyhedral_mesh_domain_with_features_3(const char* filename,
                                       CGAL::Random* p_rng)
  : Base()
  , polyhedron_()
{
  // Create input polyhedron
  std::ifstream input(filename);
  input >> polyhedron_;
  this->add_primitives(polyhedron_);
  this->set_random_generator(p_rng);
}

template < typename GT_, typename P_, typename TA_,
           typename Tag_, typename E_tag_>
Polyhedral_mesh_domain_with_features_3<GT_,P_,TA_,Tag_,E_tag_>::
Polyhedral_mesh_domain_with_features_3(const std::string& filename,
                                       CGAL::Random* p_rng)
  : Base()
  , polyhedron_()
{
  // Create input polyhedron
  std::ifstream input(filename.c_str());
  input >> polyhedron_;
  this->add_primitives(polyhedron_);
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


template < typename GT_, typename P_, typename TA_,
           typename Tag_, typename E_tag_>
void
Polyhedral_mesh_domain_with_features_3<GT_,P_,TA_,Tag_,E_tag_>::
detect_features(FT angle_in_degree, Polyhedron& p)
{
  initialize_ts(p);
  // Get sharp features
  Mesh_3::detect_features(p,angle_in_degree);
  
  // Get polylines
  typedef std::vector<Point_3> Bare_polyline;
  typedef Mesh_3::Polyline_with_context<Surface_patch_index, Curve_segment_index,
    Bare_polyline > Polyline;
  
  std::vector<Polyline> polylines;
  typedef std::back_insert_iterator<std::vector<Polyline> > Output_iterator;

  Mesh_3::detect_polylines<Polyhedron,Polyline,Output_iterator>(
    &p, std::back_inserter(polylines));
    
  // Insert polylines in domain
  Mesh_3::Extract_bare_polyline<Polyline> extractor;
  
  this->add_features(
    boost::make_transform_iterator(polylines.begin(),extractor),
    boost::make_transform_iterator(polylines.end(),extractor));
}

} //namespace CGAL


#endif // CGAL_POLYHEDRAL_MESH_DOMAIN_WITH_FEATURES_3_H
