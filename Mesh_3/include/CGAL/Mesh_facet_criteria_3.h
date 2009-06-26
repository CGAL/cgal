// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Mesh_default_facet_criteria_3 class. See class description.
//******************************************************************************

#ifndef MESH_DEFAULT_FACET_CRITERIA_3_H_
#define MESH_DEFAULT_FACET_CRITERIA_3_H_


#include <CGAL/Surface_mesher/Standard_criteria.h>
//#include <CGAL/Mesh_3/Facet_on_surface_criterion.h>

//#include <boost/optional.hpp>
#include <CGAL/Mesh_3/mesh_standard_facet_criteria.h>


namespace CGAL {

/**
 * @class Mesh_default_facet_criteria_3
 *
 *  Default facet criteria to drive Mesh_3 process
 */
template<class Tr>
class Mesh_default_facet_criteria_3
{
public:
  /// Types
  typedef Tr Triangulation;
  typedef typename Tr::Facet Facet;
  typedef typename Tr::Cell_handle Cell_handle;
  typedef typename Tr::Geom_traits::FT FT;

  typedef Surface_mesher::Refine_criterion<Tr> Criterion;
  typedef Surface_mesher::Standard_criteria<Criterion> Criteria;
  typedef typename Criteria::Quality Facet_quality;
  typedef boost::optional<Facet_quality> Facet_badness;

  /// Constructor
  Mesh_default_facet_criteria_3(const FT angle_bound,
                                const FT radius_bound,
                                const FT distance_bound)
    : curvature_size_criterion_(distance_bound)
    , uniform_size_criterion_(radius_bound)
    , aspect_ratio_criterion_(angle_bound)
//    , facet_on_surface_criterion_()
  {
    criterion_vector_.push_back (&aspect_ratio_criterion_);
    criterion_vector_.push_back (&uniform_size_criterion_);
    criterion_vector_.push_back (&curvature_size_criterion_);
//    criterion_vector_.push_back (&facet_on_surface_criterion_);

    criteria_.set_criteria(criterion_vector_);
  }

  /// Destructor
  ~Mesh_default_facet_criteria_3() { };

//  /**
//   * Tests if f is good or bad with respect to criteria
//   * @param[in] f the facet to test
//   * @param [out] q the quality of the facet if \c is bad
//   * @return \c true if \c f is a bad facet
//   */
//  bool is_bad (const Facet& f, Quality& q) const
//  {
//    return criteria_.is_bad(f, q);
//  }

  Facet_quality operator()() const
  {
    return criteria_.default_quality();
  }

  Facet_badness operator()(const Facet& f) const
  {
    Facet_quality q;
    if ( criteria_.is_bad(f,q) )
      return Facet_badness(q);
    else
      return Facet_badness();
  }

  Facet_badness operator()(const Cell_handle& c, const int i) const
  {
    return this->operator()(Facet(c,i));
  }

  /**
   * Adds a new criterion to the criterion used to test if a facet is bad
   * @param criterion the criterion to add
   *
   * TODO: fix const
   */
  void add_criterion(Criterion& criterion)
  {
    criterion_vector_.push_back(&criterion);
    criteria_.set_criteria(criterion_vector_);
  }

private:
  /// bound on Hausdorff distance does not play any role if bigger than
  /// the square of the Uniform_size_criterion
  Surface_mesher::Curvature_size_criterion<Tr> curvature_size_criterion_;

  /// bound on radii of surface Delaunay balls
  Surface_mesher::Uniform_size_criterion<Tr> uniform_size_criterion_;

  /// lower bound on minimum angle in degrees
  Surface_mesher::Aspect_ratio_criterion<Tr> aspect_ratio_criterion_;

  /// Ensure that each vertex of surface facet have dimension 2
//  Mesh_3::Facet_on_surface_criterion<Tr> facet_on_surface_criterion_;

  /// Criterion vector
  std::vector<Criterion*> criterion_vector_;
  /// Criteria
  Criteria criteria_;

private:
  // Disabled copy constructor & assignment operator
  typedef Mesh_default_facet_criteria_3<Tr> Self;
  //Mesh_default_facet_criteria_3(const Self& src);
  //Self& operator=(const Self& src);

};  // end class Mesh_default_facet_criteria_3



template<typename Tr, typename Visitor_ = Mesh_3::Facet_criterion_visitor<Tr> >
class Mesh_facet_criteria_3
{
  typedef Visitor_ Visitor;
  typedef Mesh_3::Criteria<Tr,Visitor> Criteria;

  typedef typename Tr::Facet Facet;
  typedef typename Tr::Geom_traits::FT FT;

  typedef Mesh_facet_criteria_3<Tr> Self;

public:
  typedef typename Visitor::Facet_quality Facet_quality;
  typedef typename Visitor::Facet_badness Facet_badness;

  /**
   * @brief Constructor
   */
  Mesh_facet_criteria_3(const FT& angle_bound,
                        const FT& radius_bound,
                        const FT& distance_bound)
  {
    typedef Mesh_3::Aspect_ratio_criterion<Tr,Visitor> Aspect_criterion;
    typedef Mesh_3::Uniform_size_criterion<Tr,Visitor> Uniform_size_criterion;
    typedef Mesh_3::Curvature_size_criterion<Tr,Visitor> Curvature_criterion;
    typedef Mesh_3::Facet_on_surface_criterion<Tr,Visitor> On_surface_criterion;

    if ( 0 != angle_bound )
      criteria_.add(new Aspect_criterion(angle_bound));

    if ( 0 != radius_bound )
      criteria_.add(new Uniform_size_criterion(radius_bound));

    if ( 0 != distance_bound )
      criteria_.add(new Curvature_criterion(distance_bound));

    criteria_.add(new On_surface_criterion());
  }

  /// Destructor
  ~Mesh_facet_criteria_3() { }

   /**
   * @brief returns the badness of facet \c facet
   * @param facet the facet
   * @return the badness of \c facet
   */
  Facet_badness operator()(const Facet& facet) const
  {
    return criteria_(facet);
  }

private:
  Criteria criteria_;
};  // end class Mesh_facet_criteria_3

}  // end namespace CGAL

#endif // MESH_DEFAULT_FACET_CRITERIA_3_H_
