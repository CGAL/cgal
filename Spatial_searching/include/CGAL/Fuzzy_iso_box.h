// Copyright (c) 2002,2011 Utrecht University (The Netherlands).
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
// Author(s)     : Hans Tangelder (<hanst@cs.uu.nl>)


#ifndef CGAL_FUZZY_ISO_BOX_H
#define CGAL_FUZZY_ISO_BOX_H
#include <CGAL/Kd_tree_rectangle.h>

namespace CGAL {

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

    typename Construct_min_vertex_d::result_type min, max;
    Cartesian_const_iterator_d min_begin, max_begin;
    FT eps;
    unsigned int dim;

    public:

    	// default constructor
    	Fuzzy_iso_box(const SearchTraits& traits_=SearchTraits()):traits(traits_) {}

        // constructor
	Fuzzy_iso_box(const Point_d& p, const Point_d& q, FT epsilon=FT(0),const SearchTraits& traits_=SearchTraits()) 
	  : traits(traits_), eps(epsilon)
        {
	  Construct_cartesian_const_iterator_d construct_it=traits.construct_cartesian_const_iterator_d_object();
	  Cartesian_const_iterator_d begin = construct_it(p),
	    end = construct_it(p,1);
	  dim = static_cast<unsigned int>(end - begin);

	  Iso_box_d box = typename SearchTraits::Construct_iso_box_d()(p,q);
	  Construct_min_vertex_d construct_min_vertex_d;
	  Construct_max_vertex_d construct_max_vertex_d;
	  min = construct_min_vertex_d(box);
	  max = construct_max_vertex_d(box);
	  min_begin = construct_it(min);
	  max_begin = construct_it(max);
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
