// Copyright (c) 2002,2011 Utrecht University (The Netherlands).
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
//
// Author(s)     : Hans Tangelder (<hanst@cs.uu.nl>)

// Note: Use p=0 to denote the weighted Linf-distance 
// For 0<p<1 Lp is not a metric

#ifndef CGAL_WEIGHTED_MINKOWSKI_DISTANCE_H
#define CGAL_WEIGHTED_MINKOWSKI_DISTANCE_H

#include <CGAL/license/Spatial_searching.h>


#include <cmath>
#include <vector>

#include <CGAL/array.h>
#include <CGAL/number_utils.h>
#include <CGAL/Kd_tree_rectangle.h>
#include <CGAL/internal/Get_dimension_tag.h>

namespace CGAL {
  namespace internal {
    template<class T, class Dim>
      struct Array_or_vector_selector {
	typedef std::vector<T> type;
	static void resize(type&v, std::size_t d) { v.resize(d); }
      };
    template<class T, int D>
      struct Array_or_vector_selector<T, Dimension_tag<D> > {
	typedef cpp11::array<T,D> type;
	static void resize(type&, std::size_t CGAL_assertion_code(d)) { CGAL_assertion(d==D); }
      };
  }

  template <class SearchTraits>
  class Weighted_Minkowski_distance {
    SearchTraits traits;
    public:

    typedef typename SearchTraits::Point_d Point_d;
    typedef Point_d                        Query_item;
    typedef typename SearchTraits::FT      FT;
    typedef typename internal::Get_dimension_tag<SearchTraits>::Dimension Dimension;
    typedef internal::Array_or_vector_selector<FT,Dimension> Weight_vector_traits;
    typedef typename Weight_vector_traits::type Weight_vector;

    private:

    typedef typename SearchTraits::Cartesian_const_iterator_d Coord_iterator;
    FT power; 

    Weight_vector the_weights;

    public:


    // default constructor
    Weighted_Minkowski_distance(const SearchTraits& traits_=SearchTraits())
      : traits(traits_),power(2) 
    {}

    Weighted_Minkowski_distance(const int d,const SearchTraits& traits_=SearchTraits()) 
      : traits(traits_),power(2), the_weights(d)
    {
      for (int i = 0; i < d; ++i) the_weights[i]=FT(1);
    }

    //default copy constructor and destructor
    

    template <class InputIterator>
    Weighted_Minkowski_distance (FT pow, int dim,
				 InputIterator begin,
				 InputIterator CGAL_assertion_code(end),
                                 const SearchTraits& traits_=SearchTraits()) 
      : traits(traits_),power(pow)
    {
      CGAL_assertion(power >= FT(0));
      Weight_vector_traits::resize(the_weights, dim);
      for (int i = 0; i < dim; ++i){
	the_weights[i] = *begin;
	++begin;
	CGAL_assertion(the_weights[i]>=FT(0));
      }
      CGAL_assertion(begin == end);
    }


    inline FT transformed_distance(const Query_item& q, const Point_d& p) const {
        return transformed_distance(q,p, Dimension());
    }

    //Dynamic version for runtime dimension
    inline 
    FT 
    transformed_distance(const Query_item& q, const Point_d& p, Dynamic_dimension_tag) const
    {
      FT distance = FT(0);
      typename SearchTraits::Construct_cartesian_const_iterator_d construct_it=
        traits.construct_cartesian_const_iterator_d_object();
      Coord_iterator qit = construct_it(q),
	             qe = construct_it(q,1), 
	             pit = construct_it(p);
      if (power == FT(0)) {
	for (unsigned int i = 0; qit != qe; ++qit, ++i)
	  if (the_weights[i] * CGAL::abs((*qit) - (*pit)) > distance)
	    distance = the_weights[i] * CGAL::abs((*qit)-(*pit));
      }
      else
	for (unsigned int i = 0; qit != qe; ++qit, ++i)
	  distance += 
	    the_weights[i] * std::pow(CGAL::abs((*qit)-(*pit)),power);
      return distance;
    }

    //Generic version for DIM > 3
    template <int DIM>
    inline FT 
    transformed_distance(const Query_item& q, const Point_d& p, Dimension_tag<DIM>) const
    {
      FT distance = FT(0);
      typename SearchTraits::Construct_cartesian_const_iterator_d construct_it=
        traits.construct_cartesian_const_iterator_d_object();
      Coord_iterator qit = construct_it(q),
	             qe = construct_it(q,1), 
	             pit = construct_it(p);
      if (power == FT(0)) {
	for (unsigned int i = 0; qit != qe; ++qit, ++i)
	  if (the_weights[i] * CGAL::abs((*qit) - (*pit)) > distance)
	    distance = the_weights[i] * CGAL::abs((*qit)-(*pit));
      }
      else
	for (unsigned int i = 0; qit != qe; ++qit, ++i)
	  distance += 
	    the_weights[i] * std::pow(CGAL::abs((*qit)-(*pit)),power);
      return distance;
    }

    //DIM = 2 loop unrolled
    inline FT 
    transformed_distance(const Query_item& q, const Point_d& p, Dimension_tag<2>) const
    {
      FT distance = FT(0);
      typename SearchTraits::Construct_cartesian_const_iterator_d construct_it=
        traits.construct_cartesian_const_iterator_d_object();
      Coord_iterator qit = construct_it(q),
	             pit = construct_it(p);
      if (power == FT(0)) {
	  if (the_weights[0] * CGAL::abs((*qit) - (*pit)) > distance)
	    distance = the_weights[0] * CGAL::abs((*qit)-(*pit));
          qit++;pit++;
          if (the_weights[1] * CGAL::abs((*qit) - (*pit)) > distance)
	    distance = the_weights[1] * CGAL::abs((*qit)-(*pit));
      }
      else{
	  distance += 
	    the_weights[0] * std::pow(CGAL::abs((*qit)-(*pit)),power);
          qit++;pit++;
          distance += 
	    the_weights[1] * std::pow(CGAL::abs((*qit)-(*pit)),power);
      }
      return distance;
    }

    //DIM = 3 loop unrolled
    inline FT 
    transformed_distance(const Query_item& q, const Point_d& p, Dimension_tag<3>) const
    {
      FT distance = FT(0);
      typename SearchTraits::Construct_cartesian_const_iterator_d construct_it=
        traits.construct_cartesian_const_iterator_d_object();
      Coord_iterator qit = construct_it(q),
	             pit = construct_it(p);
      if (power == FT(0)) {
	  if (the_weights[0] * CGAL::abs((*qit) - (*pit)) > distance)
	    distance = the_weights[0] * CGAL::abs((*qit)-(*pit));
          qit++;pit++;
          if (the_weights[1] * CGAL::abs((*qit) - (*pit)) > distance)
	    distance = the_weights[1] * CGAL::abs((*qit)-(*pit));
          qit++;pit++;
          if (the_weights[2] * CGAL::abs((*qit) - (*pit)) > distance)
	    distance = the_weights[2] * CGAL::abs((*qit)-(*pit));
      }
      else{
	  distance += 
	    the_weights[0] * std::pow(CGAL::abs((*qit)-(*pit)),power);
          qit++;pit++;
          distance += 
	    the_weights[1] * std::pow(CGAL::abs((*qit)-(*pit)),power);
          qit++;pit++;
          distance += 
	    the_weights[2] * std::pow(CGAL::abs((*qit)-(*pit)),power);
      }
      return distance;
    }

    inline 
    FT 
    min_distance_to_rectangle(const Query_item& q,
			      const Kd_tree_rectangle<FT,Dimension>& r) const 
    {
      FT distance = FT(0);
      typename SearchTraits::Construct_cartesian_const_iterator_d construct_it=
        traits.construct_cartesian_const_iterator_d_object();
      Coord_iterator qit = construct_it(q), qe = construct_it(q,1);
      if (power == FT(0))
	{
	  for (unsigned int i = 0; qit != qe; ++qit, ++i) {
	    if (the_weights[i]*(r.min_coord(i) - 
				(*qit)) > distance)
	      distance = the_weights[i] * (r.min_coord(i)-
					   (*qit));
	    if (the_weights[i] * ((*qit) - r.max_coord(i)) > 
		distance)
	      distance = the_weights[i] * 
		((*qit)-r.max_coord(i));
	  }
	}
      else
	{
	  for (unsigned int i = 0; qit != qe; ++qit, ++i) {
	    if ((*qit) < r.min_coord(i))
	      distance += the_weights[i] * 
		std::pow(r.min_coord(i)-(*qit),power);
	    if ((*qit) > r.max_coord(i))
	      distance += the_weights[i] * 
		std::pow((*qit)-r.max_coord(i),power);
	  }
	};
      return distance;
    }

    inline 
    FT 
    min_distance_to_rectangle(const Query_item& q,
			      const Kd_tree_rectangle<FT,Dimension>& r,std::vector<FT>& dists) const {
      FT distance = FT(0);
      typename SearchTraits::Construct_cartesian_const_iterator_d construct_it=
        traits.construct_cartesian_const_iterator_d_object();
      Coord_iterator qit = construct_it(q), qe = construct_it(q,1);
      if (power == FT(0))
	{
	  for (unsigned int i = 0; qit != qe; ++qit, ++i) {
	    if (the_weights[i]*(r.min_coord(i) - 
				(*qit)) > distance){
              dists[i] = (r.min_coord(i)-
		(*qit));
	      distance = the_weights[i] * dists[i];
            }
	    if (the_weights[i] * ((*qit) - r.max_coord(i)) > 
		distance){
                  dists[i] = 
		((*qit)-r.max_coord(i));
	      distance = the_weights[i] * dists[i];
            }
	  }
	}
      else
	{
	  for (unsigned int i = 0; qit != qe; ++qit, ++i) {
	    if ((*qit) < r.min_coord(i)){
              dists[i] = r.min_coord(i)-(*qit);
	      distance += the_weights[i] * 
		std::pow(dists[i],power);
            }
	    if ((*qit) > r.max_coord(i)){
              dists[i] = (*qit)-r.max_coord(i);
	      distance += the_weights[i] * 
		std::pow(dists[i],power);
            }
	  }
	};
      return distance;
    }

    inline 
    FT
    max_distance_to_rectangle(const Query_item& q,
			      const Kd_tree_rectangle<FT,Dimension>& r) const {
      FT distance=FT(0);
      typename SearchTraits::Construct_cartesian_const_iterator_d construct_it=
        traits.construct_cartesian_const_iterator_d_object();
      Coord_iterator qit = construct_it(q), qe = construct_it(q,1);
      if (power == FT(0))
	{
	  for (unsigned int i = 0; qit != qe; ++qit, ++i) {
	    if ((*qit) >= (r.min_coord(i) + 
			 r.max_coord(i))/FT(2.0)) {
	      if (the_weights[i] * ((*qit) - 
				    r.min_coord(i)) > distance)
		distance = the_weights[i] * 
		  ((*qit)-r.min_coord(i));
	      else
		if (the_weights[i] * 
		    (r.max_coord(i) - (*qit)) > distance)
		  distance = the_weights[i] * 
		    ( r.max_coord(i)-(*qit));
            }
	  }
	}
      else
	{
	  for (unsigned int i = 0; qit != qe; ++qit, ++i) {
	    if ((*qit) <= (r.min_coord(i)+r.max_coord(i))/FT(2.0))
	      distance += the_weights[i] * std::pow(r.max_coord(i)-(*qit),power);
	    else
	      distance += the_weights[i] * std::pow((*qit)-r.min_coord(i),power);
	  }
	};
      return distance;
    }

     inline 
    FT
    max_distance_to_rectangle(const Query_item& q,
			      const Kd_tree_rectangle<FT,Dimension>& r,std::vector<FT>& dists) const {
      FT distance=FT(0);
      typename SearchTraits::Construct_cartesian_const_iterator_d construct_it=
        traits.construct_cartesian_const_iterator_d_object();
      Coord_iterator qit = construct_it(q), qe = construct_it(q,1);
      if (power == FT(0))
	{
	  for (unsigned int i = 0; qit != qe; ++qit, ++i) {
	    if ((*qit) >= (r.min_coord(i) + 
			 r.max_coord(i))/FT(2.0)) {
	      if (the_weights[i] * ((*qit) - 
				    r.min_coord(i)) > distance){
                dists[i] = (*qit)-r.min_coord(i);
		distance = the_weights[i] * 
		  (dists[i]);
              }
	      else
		if (the_weights[i] * 
		    (r.max_coord(i) - (*qit)) > distance){
                      dists[i] =  r.max_coord(i)-(*qit);
		  distance = the_weights[i] * 
		    (dists[i]);
                }
            }
	  }
	}
      else
	{
	  for (unsigned int i = 0; qit != qe; ++qit, ++i) {
	    if ((*qit) <= (r.min_coord(i)+r.max_coord(i))/FT(2.0)){
              dists[i] = r.max_coord(i)-(*qit);
	      distance += the_weights[i] * std::pow(dists[i],power);
            }
	    else{
              dists[i] = (*qit)-r.min_coord(i);
	      distance += the_weights[i] * std::pow(dists[i],power);
            }
	  }
	};
      return distance;
    }
    
    inline 
    FT 
    new_distance(FT dist, FT old_off, FT new_off,
		 int cutting_dimension)  const 
    {
      FT new_dist;
      if (power == FT(0))
	{
	  if (the_weights[cutting_dimension]*CGAL::abs(new_off) 
	      > dist) 
	    new_dist= 
	      the_weights[cutting_dimension]*CGAL::abs(new_off);
	  else new_dist=dist;
	}
      else
	{
	  new_dist = dist + the_weights[cutting_dimension] * 
	    (std::pow(CGAL::abs(new_off),power)-std::pow(CGAL::abs(old_off),power));
	}
      return new_dist;
    }
    
    inline 
    FT 
    transformed_distance(FT d) const 
    {
      if (power <= FT(0)) return d;
      else return std::pow(d,power);
      
    }
    
    inline 
    FT 
    inverse_of_transformed_distance(FT d) const 
    {
      if (power <= FT(0)) return d;
      else return std::pow(d,1/power);
      
    }

  }; // class Weighted_Minkowski_distance

} // namespace CGAL

#endif // CGAL_WEIGHTED_MINKOWSKI_DISTANCE_H
