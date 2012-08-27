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
// 
//
// Author(s)     : Hans Tangelder (<hanst@cs.uu.nl>)


#ifndef CGAL_FUZZY_ISO_BOX_H
#define CGAL_FUZZY_ISO_BOX_H
#include <CGAL/Kd_tree_rectangle.h>
#include <boost/mpl/if.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>
#include <CGAL/Search_traits_adapter.h>

namespace CGAL {

  namespace internal{
    template <class SearchTraits,class Point>
    struct Is_from_point_from_adapter_traits{
      typedef boost::false_type type;
    };
    
    
    template <class K,class PM,class Base,class Point>
    struct Is_from_point_from_adapter_traits<Search_traits_adapter<K,PM,Base>,Point>{
      typedef typename boost::is_same<Point,typename Base::Point_d> type;
    };
  } //namespace internal
  
  template <class SearchTraits>
  class Fuzzy_iso_box{
    SearchTraits traits;
    
    public:

    typedef typename SearchTraits::Point_d Point_d;
    typedef typename SearchTraits::Iso_box_d Iso_box_d;
    typedef typename SearchTraits::FT FT;
    typedef typename SearchTraits::Construct_min_vertex_d Construct_min_vertex_d;
    typedef typename SearchTraits::Construct_max_vertex_d Construct_max_vertex_d;
    typedef typename SearchTraits::Cartesian_const_iterator_d Cartesian_const_iterator_d;
    typedef typename SearchTraits::Construct_cartesian_const_iterator_d Construct_cartesian_const_iterator_d;

    private:

    typename boost::remove_cv< 
      typename boost::remove_reference< typename Construct_min_vertex_d::result_type >::type 
      >::type min, max;
    Cartesian_const_iterator_d min_begin, max_begin;
    FT eps;
    unsigned int dim;

    //constructor implementation
    template <class Point,class Construct_iso_box_d>
    void construct(const Point& p, const Point& q)
    {
      Construct_cartesian_const_iterator_d construct_it=traits.construct_cartesian_const_iterator_d_object();
      Cartesian_const_iterator_d begin = construct_it(p),
        end = construct_it(p,1);
      dim = static_cast<unsigned int>(end - begin);

      Iso_box_d box = Construct_iso_box_d()(p,q);
      Construct_min_vertex_d construct_min_vertex_d;
      Construct_max_vertex_d construct_max_vertex_d;
      min = construct_min_vertex_d(box);
      max = construct_max_vertex_d(box);
      min_begin = construct_it(min);
      max_begin = construct_it(max);
    }
   
    public:

    	// default constructor
    	Fuzzy_iso_box(const SearchTraits& traits_=SearchTraits()):traits(traits_) {}

        // constructor
	Fuzzy_iso_box(const Point_d& p, const Point_d& q, FT epsilon=FT(0),const SearchTraits& traits_=SearchTraits()) 
	  : traits(traits_), eps(epsilon)
        {
          construct<Point_d,typename SearchTraits::Construct_iso_box_d>(p,q);
        }
        
        //additional constructor if SearchTraits = Search_traits_adapter
        template <class Point>
	Fuzzy_iso_box(const Point& p,const Point&q,FT epsilon=FT(0),const SearchTraits& traits_=SearchTraits(),
          typename boost::enable_if<typename internal::Is_from_point_from_adapter_traits<SearchTraits,Point>::type>::type* = 0) 
	  : traits(traits_), eps(epsilon)
        {
          construct<Point,typename SearchTraits::Base::Construct_iso_box_d>(p,q);
        }

        bool contains(const Point_d& p) const {	
	  Construct_cartesian_const_iterator_d construct_it=traits.construct_cartesian_const_iterator_d_object();
	  Cartesian_const_iterator_d pit = construct_it(p);
	  Cartesian_const_iterator_d minit= min_begin, maxit = max_begin; 
		for (unsigned int i = 0; i < dim; ++i, ++pit, ++minit, ++maxit) {
			if ( ((*pit) < (*minit)) || ((*pit) >= (*maxit)) ) return false;
		}
		return true; 
        }

	bool inner_range_intersects(const Kd_tree_rectangle<FT>& rectangle) const { 
	  Cartesian_const_iterator_d minit= min_begin, maxit = max_begin;   
 		for (unsigned int i = 0; i < dim; ++i, ++minit, ++maxit) {
        		if ( ((*maxit)-eps < rectangle.min_coord(i)) 
			|| ((*minit)+eps >= rectangle.max_coord(i)) ) return false;
    		}
    		return true;                                     
	}


	bool outer_range_contains(const Kd_tree_rectangle<FT>& rectangle) const {  
	  Cartesian_const_iterator_d minit= min_begin, maxit = max_begin;   
    		for (unsigned int i = 0; i < dim; ++i, ++minit, ++maxit) {
        		if (  ((*maxit)+eps < rectangle.max_coord(i) ) 
			|| ((*minit)-eps >= rectangle.min_coord(i)) ) return false;
    		}
    		return true;
  	} 

	

	

  }; // class Fuzzy_iso_box

} // namespace CGAL
#endif // FUZZY_ISO_BOX_H
