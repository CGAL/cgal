// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// $Source: /CVSROOT/CGAL/Packages/Envelope_3/include/CGAL/Envelope_caching_traits_3.h,v $
// $Revision$ $Date$
// $Name:  $
//
// Author(s)     : Michal Meyerovitch     <gorgymic@post.tau.ac.il>

#ifndef CGAL_ENV_CACHING_TRAITS_3_H
#define CGAL_ENV_CACHING_TRAITS_3_H

/*! \file
 * The caching traits class for the envelope package.
 */

#include <map>

CGAL_BEGIN_NAMESPACE

/*!
 * \class A decorator for a traits class of the envelope divide-and-conquer
 * algorithm, adding a caching ability.
 */
template <class EnvTraits_>
class Env_caching_traits_3 : public EnvTraits_
{
public:
  typedef EnvTraits_                                     Base_traits_3;
  typedef Env_caching_traits_3<Base_traits_3>            Self;
  
  typedef typename Base_traits_3::Point_2                Point_2;
  typedef typename Base_traits_3::Curve_2                Curve_2;
  typedef typename Base_traits_3::X_monotone_curve_2     X_monotone_curve_2;
  typedef typename Base_traits_3::Surface_3              Surface_3;
  typedef typename Base_traits_3::Xy_monotone_surface_3  Xy_monotone_surface_3;
 
protected:
  typedef std::pair<Curve_2, Intersection_type>          Intersection_curve;

  // caching for intersections
  typedef std::list<Object>                              Intersections_list;
  typedef std::pair<Xy_monotone_surface_3,
                    Xy_monotone_surface_3>               Surface_pair;

  struct Less_surface_pair
  {
    bool operator() (const Surface_pair& sp1,
                     const Surface_pair& sp2) const
    {
      // Compare the pairs of IDs lexicographically.
      return (sp1.first < sp2.first ||
              (sp1.first == sp2.first && sp1.second < sp2.second));
    }
  };

  // caching for construct_projected_intersections traits method's result
  typedef std::map<Surface_pair,
                   Intersections_list,
                   Less_surface_pair>                  Intersections_cache;

  // caching for compare_distance_to_envelope traits method's result
  typedef std::map<Surface_pair,
                   Comparison_result,
                   Less_surface_pair>                  Compare_cache;

public:

  class Construct_projected_intersections_2
  {
  protected:
    Self& parent;
    Intersections_cache& inter_cache;
    unsigned int& intersections_number;

  public:
    Construct_projected_intersections_2(Self* p,
                                        Intersections_cache& inter,
                                        unsigned int& inter_num)
      : parent(*p), inter_cache(inter),
        intersections_number(inter_num)
    {}
    
    // insert into OutputIterator all the (2d) projections on the xy plane of
    // the intersection objects between the 2 surfaces
    // the data type of OutputIterator is Object
    template <class OutputIterator>
    OutputIterator
    operator()(const Xy_monotone_surface_3& s1,
               const Xy_monotone_surface_3& s2,
               OutputIterator o)
    {
      Intersections_list                      inter_list;

      typename Intersections_cache::iterator  cache_iter;
      Surface_pair                            spair(s1, s2);

      // search for this pair in the cache
      cache_iter = inter_cache.find(spair);
      if (cache_iter == inter_cache.end())
      {
        parent.Base_traits_3::construct_projected_intersections_2_object()
                                  (s1, s2, std::back_inserter(inter_list));
        inter_cache[spair] = inter_list;

        // update statistics
        if (inter_list.size() > 0)
          ++intersections_number;
      }
      else
      {
        inter_list = (*cache_iter).second;
      }

      // report intersections
      typename Intersections_list::const_iterator  iter;
      for (iter = inter_list.begin(); iter != inter_list.end(); ++iter)
        *o++ = *iter;
      return o;
    }  
  };  

  /*! Get a Construct_projected_intersections_2 functor object. */
  Construct_projected_intersections_2
  construct_projected_intersections_2_object()
  {                                                                
    return Construct_projected_intersections_2(this, inter_cache,
                                               intersections_number);
  }
   
  class Compare_z_at_xy_3
  {
  protected:
    Self& parent;
    Compare_cache& compare_cache;
  public:
    Compare_z_at_xy_3(Self* p, Compare_cache& comp_cache)
      : parent(*p), compare_cache(comp_cache)
    {}

    // check which of the surfaces is closer to the envelope at the xy
    // coordinates of point (i.e. lower if computing the lower envelope, or
    // upper if computing the upper envelope)
    // precondition: the surfaces are defined in point
    Comparison_result operator()(const Point_2& p,
                                 const Xy_monotone_surface_3& s1,
                                 const Xy_monotone_surface_3& s2)
    {
      return cmp_z_at_xy(p,s1,s2);
    }

    // check which of the surfaces is closer to the envelope at the xy
    // coordinates of cv (i.e. lower if computing the lower envelope, or upper
    // if computing the upper envelope)
    // precondition: the surfaces are defined in all points of cv, and the 
    //               answer is the same for each of these points
    Comparison_result operator()(const X_monotone_curve_2& cv,
                                 const Xy_monotone_surface_3& s1,
                                 const Xy_monotone_surface_3& s2)
    {
      return cmp_z_at_xy(cv,s1,s2);
    }

  protected:
    template <class Geometry>
    Comparison_result
    cmp_z_at_xy (Geometry& g,
                 const Xy_monotone_surface_3& s1,
                 const Xy_monotone_surface_3& s2)
    {
      // check that s1 and s2 do not intersect
      Intersections_list                      inter_list;
      parent.construct_projected_intersections_2_object()
        (s1, s2,
         std::back_inserter(inter_list));
      
      // if they do not intersect, we can use the cache
      if (inter_list.empty())
      {
        // then we should check the cache of compare_distance, and only if
        // it is empty we use the traits
        Surface_pair                        spair(s1, s2);

        typename Compare_cache::iterator    comp_cache_iter;
        comp_cache_iter = compare_cache.find(spair);
        if (comp_cache_iter == compare_cache.end())
        {
          Comparison_result res =
              parent.Base_traits_3::compare_z_at_xy_3_object()
                (g, s1, s2);
          compare_cache[spair] = res;
          return res;
        }
        else
          return (*comp_cache_iter).second;
      }
      else
        // they intersect, so must use the traits
        return (parent.Base_traits_3::compare_z_at_xy_3_object()
                (g, s1, s2));
    }
  };
   
  /*! Get a Compare_z_at_xy_3 functor object. */
  Compare_z_at_xy_3
  compare_z_at_xy_3_object()
  {
    return Compare_z_at_xy_3(this, compare_cache);
  }

  /*! Default constructor. */
  Env_caching_traits_3 () :
    Base_traits_3(), 
    intersections_number(0)
  {}

  /*! Constructor with envelope type. */
  Env_caching_traits_3 (Envelope_type t) : 
    Base_traits_3(t),
    intersections_number(0)
  {}

  virtual ~Env_caching_traits_3()
  {
    reset();
  }

  // methods for benchmarks
  void reset()
  {
    intersections_number = 0;
    inter_cache.clear();
    compare_cache.clear();
  }

protected:

  mutable Intersections_cache inter_cache;
  mutable unsigned int intersections_number;

  mutable Compare_cache       compare_cache;
};

CGAL_END_NAMESPACE

#endif
