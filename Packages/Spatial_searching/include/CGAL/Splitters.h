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
// release       : $CGAL_Revision: CGAL-2.5-I-99 $
// release_date  : $CGAL_Date: 2003/05/23 $
//
// file          : include/CGAL/Splitters.h
// package       : ASPAS (3.12)
// maintainer    : Hans Tangelder <hanst@cs.uu.nl>
// revision      : 3.0
// revision_date : 2003/07/10 
// authors       : Hans Tangelder (<hanst@cs.uu.nl>)
// coordinator   : Utrecht University
//
// ======================================================================

// Defines rules used for constructing a split node. That is, it implements, 
// in several ways, the concept
// Boxtree_splitter<NT>.

#ifndef CGAL_SPLITTERS_H
#define CGAL_SPLITTERS_H
#include <CGAL/Point_container.h>
#include <CGAL/Plane_separator.h>
namespace CGAL {

template <class Point, class Container_=Point_container<Point>, 
          class Separator_=Plane_separator<typename Kernel_traits<Point>::Kernel::FT>  >
class Median_of_max_spread {
public:
  typedef typename Kernel_traits<Point>::Kernel::FT NT;
  typedef Container_ Container;
  typedef Separator_ Separator;
 
  void operator() (Separator& sep, Container& c0,
  			     Container& c1, 
  			     NT Aspect_ratio=NT(3)) {        
        sep=Separator(c0.max_tight_span_coord(),NT(0));
        sep.set_cutting_value(c0.median(sep.cutting_dimension()));
        c0.split_container(c1,sep,true);
  }
};

template <class Point, class Container_=Point_container<Point>,
	  class Separator_=Plane_separator<typename Kernel_traits<Point>::Kernel::FT>  > 
class Fair {
public:
  typedef typename Kernel_traits<Point>::Kernel::FT NT;
  typedef Container_ Container;
  typedef Separator_ Separator;

  void operator() (Separator& sep, Container& c0,
  			     Container& c1, 
  			     NT Aspect_ratio=NT(3)) {
	// find legal cut with max spread
        sep=Separator(c0.max_tight_span_coord_balanced(Aspect_ratio),
				NT(0));
        sep.set_cutting_value(c0.balanced_fair(sep.cutting_dimension(),
				Aspect_ratio));
        c0.split_container(c1,sep);
  }
};

template <class Point,  class Container_=Point_container<Point>,
	  class Separator_=Plane_separator<typename Kernel_traits<Point>::Kernel::FT>  >
class Sliding_fair {
public:
  typedef typename Kernel_traits<Point>::Kernel::FT NT;
  typedef Container_ Container;
  typedef Separator_ Separator;

  void operator() (Separator& sep, Container& c0,
  			     Container& c1, 
  			     NT Aspect_ratio=NT(3))  {
    // find legal cut with max spread
   
    sep=Separator(c0.max_tight_span_coord_balanced(Aspect_ratio),
			    NT(0));
    
    sep.set_cutting_value(c0.balanced_sliding_fair(sep.cutting_dimension(),
			 Aspect_ratio));
    c0.split_container(c1,sep,true);
  }
};


template <class Point,  class Container_=Point_container<Point>,
	  class Separator_=Plane_separator<typename Kernel_traits<Point>::Kernel::FT>  >
class Sliding_midpoint {
public:
  typedef typename Kernel_traits<Point>::Kernel::FT NT;
  typedef Container_ Container;
  typedef Separator_ Separator;

  void operator() (Separator& sep, Container& c0,
  			     Container& c1, 
  			     NT Aspect_ratio=NT(3))
  {
        sep=Separator(c0.max_span_coord(),
              (c0.max_span_upper() + c0.max_span_lower())/NT(2));
	NT max_span_lower = 
	c0.tight_bounding_box().min_coord(c0.max_span_coord());
	NT max_span_upper = 
	c0.tight_bounding_box().max_coord(c0.max_span_coord());
	if (max_span_upper <= sep.cutting_value()) {
		sep.set_cutting_value(max_span_upper); 
	};
	if (max_span_lower >= sep.cutting_value()) {
		sep.set_cutting_value(max_span_lower); 
	};
	c0.split_container(c1,sep,true);
  }
};

template <class Point,  class Container_=Point_container<Point>,
	  class Separator_=Plane_separator<typename Kernel_traits<Point>::Kernel::FT>  >
class Median_of_rectangle {
public:
  typedef typename Kernel_traits<Point>::Kernel::FT NT;
  typedef Container_ Container;
  typedef Separator_ Separator;

  void operator() (Separator& sep, Container& c0,
  			     Container& c1, 
  			     NT Aspect_ratio=NT(3))
  {
    sep=Separator(c0.max_span_coord(),NT(0));
    sep.set_cutting_value(c0.median(sep.cutting_dimension()));
    c0.split_container(c1,sep,true);
  }
};

template <class Point,  class Container_=Point_container<Point>,
	  class Separator_=Plane_separator<typename Kernel_traits<Point>::Kernel::FT>  >
class Midpoint_of_max_spread {
public:
  typedef typename Kernel_traits<Point>::Kernel::FT NT;
  typedef Container_ Container;
  typedef Separator_ Separator;

  void operator() (Separator& sep, Container& c0,
  			     Container& c1, 
  			     NT Aspect_ratio=NT(3))
  {
    sep= Separator(c0.max_tight_span_coord(),
    (c0.max_tight_span_upper() + c0.max_tight_span_lower())/NT(2));
    c0.split_container(c1,sep);
  }
};

template <class Point, class Container_=Point_container<Point>,
	  class Separator_=Plane_separator<typename Kernel_traits<Point>::Kernel::FT>  >
class Midpoint_of_rectangle {
public:
  typedef typename Kernel_traits<Point>::Kernel::FT NT;
  typedef Container_ Container;
  typedef Separator_ Separator;
  void operator() (Separator& sep, Container& c0,
  			     Container& c1, 
  			     NT Aspect_ratio=NT(3))
  {
    sep = Separator(c0.max_span_coord(),
              (c0.max_span_upper() + c0.max_span_lower())/NT(2));
    c0.split_container(c1,sep);          
  }
 
};

} // namespace CGAL
#endif // CGAL_SPLITTERS



