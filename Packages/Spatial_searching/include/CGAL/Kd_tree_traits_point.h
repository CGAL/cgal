// Copyright (c) 2002  Utrecht University (The Netherlands).
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Authors       : Hans Tangelder (<hanst@cs.uu.nl>)


#ifndef CGAL_KD_TREE_TRAITS_POINT_H
#define CGAL_KD_TREE_TRAITS_POINT_H
#include <CGAL/Splitters.h>

namespace CGAL {

  template <class Point_, 
            class Splitter=Sliding_midpoint<Point_> >
  class Kd_tree_traits_point {

  public:
    typedef Point_ Point;
    typedef Point** Point_iterator;
    typedef typename Kernel_traits<Point>::Kernel::FT NT;
    typedef typename Splitter::Container Container;
    typedef typename Splitter::Separator Separator;
    
  private:

    unsigned int the_bucket_size;
    NT the_aspect_ratio;
    bool use_extended_nodes_option;

  public:

        //default constructor

	Kd_tree_traits_point() {
		the_bucket_size = 5;		
		the_aspect_ratio = NT(3);
		use_extended_nodes_option = true;
        }
               
	Kd_tree_traits_point(unsigned int bucket_size, 
			     NT aspect_ratio=NT(3), 
			     bool use_extended_nodes=true) {
		the_bucket_size = bucket_size;
		the_aspect_ratio = aspect_ratio;
		use_extended_nodes_option = use_extended_nodes;
	}

    	NT aspect_ratio() const {return the_aspect_ratio;}
	
        unsigned int bucket_size() const {return the_bucket_size;}
	bool use_extended_nodes() const {return use_extended_nodes_option;}

	// split c0 in c0 and c1
    	void split(Separator& sep,
    	           Container& c0, Container& c1) 
	{
		Splitter S;
		S(sep,c0,c1,the_aspect_ratio);
    	}

  };

  
} // namespace CGAL
#endif // KD_TREE_TRAITS_POINT_H
