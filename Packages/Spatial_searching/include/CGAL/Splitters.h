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
// Boxtree_splitter<FT>.

#ifndef CGAL_SPLITTERS_H
#define CGAL_SPLITTERS_H
#include <CGAL/Point_container.h>
#include <CGAL/Plane_separator.h>
namespace CGAL {

  template <class FT>
  class Splitter_base {
  private:
    unsigned int the_bucket_size;
    FT the_aspect_ratio;

public:
    Splitter_base(unsigned int bucket_size = 3, 
		  FT aspect_ratio = FT(3))
      : the_bucket_size(bucket_size),
	the_aspect_ratio(aspect_ratio)
    {}

    FT 
    aspect_ratio() const {
      return the_aspect_ratio;
    }
    
    unsigned int 
    bucket_size() const {
      return the_bucket_size;
    }

  };


template <class SearchTraits, class Container_=Point_container<SearchTraits>, 
          class Separator_=Plane_separator<typename SearchTraits::FT>  >
class Median_of_max_spread 
  : public Splitter_base<typename SearchTraits::FT> 
{

  typedef Splitter_base<typename SearchTraits::FT> Base;
public:
  typedef typename SearchTraits::FT FT;
  typedef Container_ Container;
  typedef Separator_ Separator;
 
  Median_of_max_spread()
    : Base()
  {}

  Median_of_max_spread(unsigned int bucket_size)
    : Base(bucket_size)
  {}

  void operator() (Separator& sep, 
		   Container& c0,
		   Container& c1) {        
        sep=Separator(c0.max_tight_span_coord(),FT(0));
        sep.set_cutting_value(c0.median(sep.cutting_dimension()));
        c0.split(c1,sep,true);
  }
};

template <class SearchTraits, class Container_=Point_container<SearchTraits>,
	  class Separator_=Plane_separator<typename SearchTraits::FT>  > 
class Fair
  : public Splitter_base<typename SearchTraits::FT> {

  typedef Splitter_base<typename SearchTraits::FT> Base;
public:
  typedef typename SearchTraits::FT FT;
  typedef Container_ Container;
  typedef Separator_ Separator;

  Fair()
    : Base()
  {}

  Fair(unsigned int bucket_size, 
       FT aspect_ratio=FT(3))
    : Base(bucket_size, aspect_ratio)
  {}

  void operator() (Separator& sep, Container& c0,
  			     Container& c1) {
	// find legal cut with max spread
        sep=Separator(c0.max_tight_span_coord_balanced(aspect_ratio()),
				FT(0));
        sep.set_cutting_value(c0.balanced_fair(sep.cutting_dimension(),
				aspect_ratio()));
        c0.split(c1,sep);
  }
};

template <class SearchTraits,  class Container_=Point_container<SearchTraits>,
	  class Separator_=Plane_separator<typename SearchTraits::FT>  >
class Sliding_fair
  : public Splitter_base<typename SearchTraits::FT> {

  typedef Splitter_base<typename SearchTraits::FT> Base;

public:
  typedef typename SearchTraits::FT FT;
  typedef Container_ Container;
  typedef Separator_ Separator;

  Sliding_fair()
    : Base()
  {}

  Sliding_fair(unsigned int bucket_size, 
	       FT aspect_ratio=FT(3))
    : Base(bucket_size, aspect_ratio)
  {}

  void operator() (Separator& sep, Container& c0,
  			     Container& c1)  {
    // find legal cut with max spread
   
    sep=Separator(c0.max_tight_span_coord_balanced(aspect_ratio()),
			    FT(0));
    
    sep.set_cutting_value(c0.balanced_sliding_fair(sep.cutting_dimension(),
			 aspect_ratio()));
    c0.split(c1,sep,true);
  }
};


template <class SearchTraits,  class Container_=Point_container<SearchTraits>,
	  class Separator_=Plane_separator<typename SearchTraits::FT>  >
class Sliding_midpoint
  : public Splitter_base<typename SearchTraits::FT> {

  typedef Splitter_base<typename SearchTraits::FT> Base;

public:
  typedef typename SearchTraits::FT FT;
  typedef Container_ Container;
  typedef Separator_ Separator;

  Sliding_midpoint()
    : Base()
  {}

  Sliding_midpoint(unsigned int bucket_size)
    : Base(bucket_size)
  {}

  void operator() (Separator& sep, Container& c0,
  			     Container& c1)
  {
        sep=Separator(c0.max_span_coord(),
              (c0.max_span_upper() + c0.max_span_lower())/FT(2));
	FT max_span_lower = 
	c0.tight_bounding_box().min_coord(c0.max_span_coord());
	FT max_span_upper = 
	c0.tight_bounding_box().max_coord(c0.max_span_coord());
	if (max_span_upper <= sep.cutting_value()) {
		sep.set_cutting_value(max_span_upper); 
	};
	if (max_span_lower >= sep.cutting_value()) {
		sep.set_cutting_value(max_span_lower); 
	};
	c0.split(c1,sep,true);
  }
};

template <class SearchTraits,  class Container_=Point_container<SearchTraits>,
	  class Separator_=Plane_separator<typename SearchTraits::FT>  >
class Median_of_rectangle
  : public Splitter_base<typename SearchTraits::FT> {

  typedef Splitter_base<typename SearchTraits::FT> Base;

public:
  typedef typename SearchTraits::FT FT;
  typedef Container_ Container;
  typedef Separator_ Separator;


  Median_of_rectangle()
    : Base()
  {}

  Median_of_rectangle(unsigned int bucket_size)
    : Base(bucket_size)
  {}

  void operator() (Separator& sep, Container& c0,
  			     Container& c1)
  {
    sep=Separator(c0.max_span_coord(),FT(0));
    sep.set_cutting_value(c0.median(sep.cutting_dimension()));
    c0.split(c1,sep,true);
  }
};

template <class SearchTraits,  class Container_=Point_container<SearchTraits>,
	  class Separator_=Plane_separator<typename SearchTraits::FT>  >
class Midpoint_of_max_spread
  : public Splitter_base<typename SearchTraits::FT> {

  typedef Splitter_base<typename SearchTraits::FT> Base;

public:
  typedef typename SearchTraits::FT FT;
  typedef Container_ Container;
  typedef Separator_ Separator;


  Midpoint_of_max_spread()
    : Base()
  {}

  Midpoint_of_max_spread(unsigned int bucket_size)
    : Base(bucket_size)
  {}

  void operator() (Separator& sep, Container& c0,
  			     Container& c1)
  {
    sep= Separator(c0.max_tight_span_coord(),
    (c0.max_tight_span_upper() + c0.max_tight_span_lower())/FT(2));
    c0.split(c1,sep);
  }
};

template <class SearchTraits, class Container_=Point_container<SearchTraits>,
	  class Separator_=Plane_separator<typename SearchTraits::FT>  >
class Midpoint_of_rectangle
  : public Splitter_base<typename SearchTraits::FT> {

  typedef Splitter_base<typename SearchTraits::FT> Base;
public:
  typedef typename SearchTraits::FT FT;
  typedef Container_ Container;
  typedef Separator_ Separator;


  Midpoint_of_rectangle()
    : Base()
  {}

  Midpoint_of_rectangle(unsigned int bucket_size)
    : Base(bucket_size)
  {}

  void operator() (Separator& sep, Container& c0,
  			     Container& c1)
  {
    sep = Separator(c0.max_span_coord(),
		    (c0.max_span_upper() + c0.max_span_lower())/FT(2));
    c0.split(c1,sep);          
  }
 
};

} // namespace CGAL
#endif // CGAL_SPLITTERS



