// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
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
//
//******************************************************************************
// File Description :
// class Labeled_mesh_domain_3. See class description.
//******************************************************************************

#ifndef CGAL_POLYHEDRON_DEMO_LABELED_MESH_DOMAIN_3_H
#define CGAL_POLYHEDRON_DEMO_LABELED_MESH_DOMAIN_3_H

#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Random.h>
#include "Image_type.h"

#include <boost/type_traits.hpp>

namespace CGAL {

/**
 * \class Polyhedron_demo_labeled_mesh_domain_3
 * LabeledDomain must be a Labeled_mesh_domain_3
 */
template<class LabeledDomain, class Image = CGAL::Image_3>
class Polyhedron_demo_labeled_mesh_domain_3
  : public LabeledDomain
{
public:
  typedef LabeledDomain Base;
  typedef Polyhedron_demo_labeled_mesh_domain_3 Domain;

  //-------------------------------------------------------
  // Index Types
  //-------------------------------------------------------
  /// Type of indexes for cells of the input complex
  typedef typename Base::Subdomain_index       Subdomain_index;
  typedef boost::optional<Subdomain_index>     Subdomain;

  /// Type of indexes for surface patch of the input complex
  typedef int                                  Surface_patch_index;
  typedef boost::optional<Surface_patch_index> Surface_patch;

  typedef typename Base::Point_3   Point_3;
  typedef typename Base::Segment_3 Segment_3;

  /// Type of indexes to characterize the lowest dimensional face of the input
  /// complex on which a vertex lie
  typedef int Index;
  typedef CGAL::cpp11::tuple<Point_3, Index, int> Intersection;

  //constructors
  Polyhedron_demo_labeled_mesh_domain_3(
    const typename Base::Fct& f,
    const typename Base::Bbox_3& bbox,
    const typename Base::FT& error_bound = Base::FT(1e-3),
    CGAL::Random* p_rng = NULL)
    : Base(f, bbox, error_bound, p_rng)
  {}

  Polyhedron_demo_labeled_mesh_domain_3(
    const Image& img,
    const typename Base::FT& error_bound = Base::FT(1e-3),
    CGAL::Random* p_rng = NULL)
    : Base(img, error_bound, p_rng)
  {}

  /**
   * Returns the index to be stored in a vertex lying on the surface identified
   * by \c index.
   */
  Index index_from_surface_patch_index(
    const Surface_patch_index& index) const
  {
    return Index(index);
  }

  Index index_from_surface_patch_index(
    const typename Base::Surface_patch_index& index_pair) const
  {
    return Index(index_pair.first * 1000 + index_pair.second);
  }

  Index index_from_surface_patch_index(
    const typename Base::Index& index) const
  {
    return index_from_surface_patch_index(surface_patch_index(index));
  }

  Index index_from_index(
    const typename Base::Index& index) const
  {
    if (const typename Base::Surface_patch_index* index_pair =
      boost::get<typename Base::Surface_patch_index>(&index))
        return Index(index_pair->first * 1000 + index_pair->second);

    else if (const typename Base::Subdomain_index* sub_index =
      boost::get<typename Base::Subdomain_index>(&index))
        return Index(*sub_index);

    return Index(-1);//~error
  }

  /**
   * Returns the index to be stored in a vertex lying in the subdomain
   * identified by \c index.
   */
  Index index_from_subdomain_index(const Subdomain_index& index) const
  { return Index(index); }

  /**
   * Returns the \c Surface_patch_index of the surface patch
   * where lies a vertex with dimension 2 and index \c index.
   */
  Surface_patch_index
    surface_patch_index(const Index& index) const
  {
    return index;
  }

  Surface_patch_index
    surface_patch_index(const typename Base::Index& index) const
  {
    typename Base::Surface_patch_index index_pair =
      boost::get<typename Base::Surface_patch_index>(index);
    return Surface_patch_index(index_pair.first * 1000 + index_pair.second);
  }

  typename Base::Surface_patch_index
    surface_patch_index_base(const Index& index) const
  {
    return typename Base::Surface_patch_index(index / 1000,
                                              index % 1000);
  }

  /**
   * Returns the index of the subdomain containing a vertex
   *  with dimension 3 and index \c index.
   */
  Subdomain_index subdomain_index(const Index& index) const
  { return index; }
  
  // -----------------------------------
  // Backward Compatibility
  // -----------------------------------
#ifndef CGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX
  typedef Surface_patch_index   Surface_index;
  
  Index index_from_surface_index(const Surface_index& index) const
  { return index_from_surface_patch_index(index); }
  
  Surface_index surface_index(const Index& index) const
  { return surface_patch_index(index); }
#endif // CGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX
  // -----------------------------------
  // End backward Compatibility
  // -----------------------------------


  struct Construct_initial_points
  {
    Construct_initial_points(const Domain& domain)
      : r_domain_(domain) {}

    template<class OutputIterator>
    OutputIterator operator()(OutputIterator pts, const int n = 12) const
    {
      typename Base::Construct_initial_points cip =
        r_domain_.Base::construct_initial_points_object();

      std::vector<std::pair<Point_3, typename Base::Index> > initial_points;
      cip(std::back_inserter(initial_points), n);

      for (std::size_t i = 0; i < initial_points.size(); ++i)
      {
        std::pair<Point_3, typename Base::Index> p = initial_points[i];
        *pts++ = std::make_pair(p.first,
                                r_domain_.surface_patch_index(p.second));
      }
      return pts;
    }

  private:
    const Domain& r_domain_;
  };

  /// Returns Construct_initial_points object
  Construct_initial_points construct_initial_points_object() const
  {
    return Construct_initial_points(*this);
  }


  struct Construct_intersection
  {
    Construct_intersection(const Domain& domain)
      : r_domain_(domain) {}

    template<typename Query>
    Intersection operator()(const Query& q) const
    {
      typename Base::Construct_intersection
        ci = r_domain_.Base::construct_intersection_object();
      typename Base::Intersection bi = ci(q);

      return CGAL::cpp11::make_tuple(
        CGAL::cpp11::get<0>(bi),
        r_domain_.index_from_index(CGAL::cpp11::get<1>(bi)),
        CGAL::cpp11::get<2>(bi));
    }

  private:
    const Domain& r_domain_;
  };

  /// Returns Construct_intersection object
  Construct_intersection construct_intersection_object() const
  {
    return Construct_intersection(*this);
  }


private:
  ///// Returns Surface_patch_index from \c i and \c j
  //typename Base::Surface_patch_index
  //  make_surface_index(const Subdomain_index i,
  //                     const Subdomain_index j) const
  //{
  //  if ( i < j )
  //    return typename Base::Surface_patch_index(i, j);
  //  else
  //    return typename Base::Surface_patch_index(j, i);
  //}

  Surface_patch_index
    make_surface_index(const Subdomain_index i,
                       const Subdomain_index j) const
  {
    if (i < j)
      return (i * 1000 + j);
    else
      return (j * 1000 + i);
  }

};  // end class Labeled_mesh_domain_3

}  // end namespace CGAL

#endif // CGAL_POLYHEDRON_DEMO_LABELED_MESH_DOMAIN_3_H

