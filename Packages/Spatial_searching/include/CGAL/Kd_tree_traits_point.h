// ======================================================================
//
// Copyright (c) 2002 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Kd_tree_traits_point.h
// package       : ASPAS
// revision      : 3.0
// revision_date : 2003/07/10 
// authors       : Hans Tangelder (<hanst@cs.uu.nl>)
// maintainer    : Hans Tangelder (<hanst@cs.uu.nl>)
// coordinator   : Utrecht University
//
// ======================================================================


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
		the_bucket_size = 100;		
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
