// Copyright (c) 2002,2011 Utrecht University (The Netherlands).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Hans Tangelder (<hanst@cs.uu.nl>)


#ifndef CGAL_MANHATTAN_DISTANCE_ISO_BOX_POINT_H
#define CGAL_MANHATTAN_DISTANCE_ISO_BOX_POINT_H

#include <CGAL/license/Spatial_searching.h>


#include <CGAL/result_of.h>
#include <CGAL/Kd_tree_rectangle.h>
#include <CGAL/internal/Get_dimension_tag.h>
#include <vector>

namespace CGAL {

  template <class SearchTraits>
  class Manhattan_distance_iso_box_point {
    SearchTraits traits;
  public:

    typedef typename SearchTraits::Point_d Point_d;
    typedef typename SearchTraits::Iso_box_d Iso_box_d;
    typedef typename SearchTraits::FT    FT;
    typedef Iso_box_d                  Query_item;
    typedef typename internal::Get_dimension_tag<SearchTraits>::Dimension Dimension;
    typedef typename CGAL::cpp11::result_of<typename SearchTraits::Construct_max_vertex_d(Query_item)>::type Max_vertex;

    Manhattan_distance_iso_box_point(const SearchTraits& traits_=SearchTraits()):traits(traits_) {}


    // obsolete as we no longer store dimension Manhattan_distance_iso_box_point(const int d) : the_dimension(d) {}

    inline FT transformed_distance(const Query_item& q, const Point_d& p) const {
                FT distance = FT(0);
                typename SearchTraits::Construct_cartesian_const_iterator_d construct_it=
                  traits.construct_cartesian_const_iterator_d_object();
                typename SearchTraits::Construct_min_vertex_d construct_min_vertex;
                typename SearchTraits::Construct_max_vertex_d construct_max_vertex;
                Max_vertex maxv = construct_max_vertex(q);
                typename SearchTraits::Cartesian_const_iterator_d qmaxit = construct_it(maxv),
                  qe = construct_it(maxv,1), qminit = construct_it(construct_min_vertex(q)),
                  pit = construct_it(p);
                for (; qmaxit != qe; ++pit,++qmaxit,++qminit) {
                        if ((*pit)>(*qmaxit)) distance +=
                        (*pit)-(*qmaxit);
                        else if ((*pit)<(*qminit)) distance +=
                        (*qminit)-(*pit);
                }
                return distance;
    }


    inline FT min_distance_to_rectangle(const Query_item& q,
                                              const Kd_tree_rectangle<FT,Dimension>& r) const {
                FT distance = FT(0);
                typename SearchTraits::Construct_cartesian_const_iterator_d construct_it=
                  traits.construct_cartesian_const_iterator_d_object();
                typename SearchTraits::Construct_min_vertex_d construct_min_vertex;
                typename SearchTraits::Construct_max_vertex_d construct_max_vertex;
                Max_vertex maxv = construct_max_vertex(q);
                typename SearchTraits::Cartesian_const_iterator_d qmaxit = construct_it(maxv),
                  qe = construct_it(maxv,1), qminit = construct_it(construct_min_vertex(q));
                for (unsigned int i = 0; qmaxit != qe; ++ qmaxit, ++qminit, ++i)  {
                        if (r.min_coord(i)>(*qmaxit))
                          distance +=(r.min_coord(i)-(*qmaxit));
                        if (r.max_coord(i)<(*qminit))
                          distance += ((*qminit)-r.max_coord(i));
                }
                return distance;
        }

     inline FT min_distance_to_rectangle(const Query_item& q,
                                              const Kd_tree_rectangle<FT,Dimension>& r,std::vector<FT>& dists) const {
                FT distance = FT(0);
                typename SearchTraits::Construct_cartesian_const_iterator_d construct_it=
                  traits.construct_cartesian_const_iterator_d_object();
                typename SearchTraits::Construct_min_vertex_d construct_min_vertex;
                typename SearchTraits::Construct_max_vertex_d construct_max_vertex;
                Max_vertex maxv = construct_max_vertex(q);
                typename SearchTraits::Cartesian_const_iterator_d qmaxit = construct_it(maxv),
                  qe = construct_it(maxv,1), qminit = construct_it(construct_min_vertex(q));
                for (unsigned int i = 0; qmaxit != qe; ++ qmaxit, ++qminit, ++i)  {
                        if (r.min_coord(i)>(*qmaxit)) {
                          dists[i]=(r.min_coord(i)-(*qmaxit));
                          distance +=dists[i];
                        }
                        if (r.max_coord(i)<(*qminit)) {
                          dists[i]=((*qminit)-r.max_coord(i));
                          distance += dists[i];
                        }
                }
                return distance;
        }

    inline
    FT
    max_distance_to_rectangle(const Query_item& q,
                              const Kd_tree_rectangle<FT,Dimension>& r) const {
      FT distance=FT(0);
      typename SearchTraits::Construct_cartesian_const_iterator_d construct_it=
        traits.construct_cartesian_const_iterator_d_object();
      typename SearchTraits::Construct_min_vertex_d construct_min_vertex;
      typename SearchTraits::Construct_max_vertex_d construct_max_vertex;
      Max_vertex maxv = construct_max_vertex(q);
      typename SearchTraits::Cartesian_const_iterator_d qmaxit = construct_it(maxv),
        qe = construct_it(maxv,1), qminit = construct_it(construct_min_vertex(q));
      for (unsigned int i = 0; qmaxit != qe; ++ qmaxit, ++qminit, ++i)  {
        if ( r.max_coord(i)-(*qminit) >(*qmaxit)-r.min_coord(i) )
          distance += (r.max_coord(i)-(*qminit));
        else
          distance += ((*qmaxit)-r.min_coord(i));
      }
      return distance;
    }

     inline
    FT
    max_distance_to_rectangle(const Query_item& q,
                              const Kd_tree_rectangle<FT,Dimension>& r,std::vector<FT>& dists) const {
      FT distance=FT(0);
      typename SearchTraits::Construct_cartesian_const_iterator_d construct_it=
        traits.construct_cartesian_const_iterator_d_object();
      typename SearchTraits::Construct_min_vertex_d construct_min_vertex;
      typename SearchTraits::Construct_max_vertex_d construct_max_vertex;
      Max_vertex maxv = construct_max_vertex(q);
      typename SearchTraits::Cartesian_const_iterator_d qmaxit = construct_it(maxv),
        qe = construct_it(maxv,1), qminit = construct_it(construct_min_vertex(q));
      for (unsigned int i = 0; qmaxit != qe; ++ qmaxit, ++qminit, ++i)  {
        if ( r.max_coord(i)-(*qminit) >(*qmaxit)-r.min_coord(i) )  {
          dists[i]=(r.max_coord(i)-(*qminit));
          distance += dists[i];
        }
        else {
          dists[i]=((*qmaxit)-r.min_coord(i));
          distance += dists[i];
        }
      }
      return distance;
    }

  inline
  FT
  transformed_distance(FT d) const
  {
    return d;
  }

  inline
  FT
  inverse_of_transformed_distance(FT d) const
  {
    return d;
  }

}; // class Manhattan_distance_iso_box_point

} // namespace CGAL
#endif // MANHATTAN_DISTANCE_ISO_BOX_POINT_H
