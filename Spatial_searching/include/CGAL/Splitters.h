// Copyright (c) 2002-2011 Utrecht University (The Netherlands).
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
    aspect_ratio() const 
    {
      return the_aspect_ratio;
    }
    
    unsigned int 
    bucket_size() const 
    {
      return the_bucket_size;
    }

  };


  template <class SearchTraits,
            class Separator_= Plane_separator<typename SearchTraits::FT>  >
  class Median_of_max_spread 
    : public Splitter_base<typename SearchTraits::FT> 
  {

    typedef Splitter_base<typename SearchTraits::FT> Base;
  public:
    typedef typename SearchTraits::FT FT;
    typedef Point_container<SearchTraits> Container;
    typedef Separator_ Separator;
    
    Median_of_max_spread()
      : Base()
    {}

    Median_of_max_spread(unsigned int bucket_size)
      : Base(bucket_size)
    {}
    
    void 
    operator() (Separator& sep, 
		Container& c0,
		Container& c1) const
    {        
      sep=Separator(c0.max_tight_span_coord(),FT(0));
      sep.set_cutting_value(c0.median(sep.cutting_dimension()));
      c0.split(c1,sep,true);
    }
  };

  template <class SearchTraits,
	    class Separator_= Plane_separator<typename SearchTraits::FT>  > 
  class Fair
    : public Splitter_base<typename SearchTraits::FT> 
  {

    typedef Splitter_base<typename SearchTraits::FT> Base;
  public:
    typedef typename SearchTraits::FT FT;
    typedef Point_container<SearchTraits> Container;
    typedef Separator_ Separator;
    
    Fair()
      : Base()
    {}

    Fair(unsigned int bucket_size, 
	 FT aspect_ratio=FT(3))
      : Base(bucket_size, aspect_ratio)
    {}

    void 
    operator()(Separator& sep, Container& c0, Container& c1) const
    {
      // find legal cut with max spread
      sep=Separator(c0.max_tight_span_coord_balanced(this->aspect_ratio()),
		    FT(0));
      sep.set_cutting_value(c0.balanced_fair(sep.cutting_dimension(),
					     this->aspect_ratio()));
      c0.split(c1,sep);
    }
  };

  template <class SearchTraits,
           class Separator_= Plane_separator<typename SearchTraits::FT>  >
  class Sliding_fair
    : public Splitter_base<typename SearchTraits::FT> 
  {
    
    typedef Splitter_base<typename SearchTraits::FT> Base;
    
  public:
    typedef typename SearchTraits::FT FT;
    typedef Point_container<SearchTraits> Container;
    typedef Separator_ Separator;
    
    Sliding_fair()
      : Base()
    {}
    
    Sliding_fair(unsigned int bucket_size, 
		 FT aspect_ratio=FT(3))
      : Base(bucket_size, aspect_ratio)
    {}
    
    void 
    operator() (Separator& sep, Container& c0, Container& c1)  const
    {
      // find legal cut with max spread
      
      sep = Separator(c0.max_tight_span_coord_balanced(this->aspect_ratio()),
		      FT(0));
      
      sep.set_cutting_value(c0.balanced_sliding_fair(sep.cutting_dimension(),
						     this->aspect_ratio()));
      c0.split(c1,sep,true);
    }
  };


  template <class SearchTraits,
	    class Separator_= Plane_separator<typename SearchTraits::FT>  >
  class Sliding_midpoint
    : public Splitter_base<typename SearchTraits::FT> 
  {
    
    typedef Splitter_base<typename SearchTraits::FT> Base;
    
  public:
    typedef typename SearchTraits::FT FT;
    typedef Point_container<SearchTraits> Container;
    typedef Separator_ Separator;
    
    Sliding_midpoint()
      : Base()
    {}
    
    Sliding_midpoint(unsigned int bucket_size)
      : Base(bucket_size)
    {}

    void 
    operator()(Separator& sep, Container& c0, Container& c1) const
    {
      CGAL_assertion(c0.is_valid());
      CGAL_assertion(c1.is_valid());
      sep = Separator(c0.max_span_coord(),
		      (c0.max_span_upper() + c0.max_span_lower())/FT(2));

      FT max_span_lower = 
	c0.tight_bounding_box().min_coord(c0.max_span_coord());

      CGAL_assertion(max_span_lower >= c0.max_span_lower());
      FT max_span_upper = 
	c0.tight_bounding_box().max_coord(c0.max_span_coord());
      
      CGAL_assertion(max_span_upper <= c0.max_span_upper());
      if (max_span_upper <= sep.cutting_value()) {
	sep.set_cutting_value(max_span_upper); 
      };      
      if (max_span_lower >= sep.cutting_value()) {    
	sep.set_cutting_value(max_span_lower); 
      };   
      c0.split(c1,sep,true); 
    }
  };
  
  template <class SearchTraits,
	    class Separator_= Plane_separator<typename SearchTraits::FT>  >
  class Median_of_rectangle
    : public Splitter_base<typename SearchTraits::FT> 
  {
    
    typedef Splitter_base<typename SearchTraits::FT> Base;
    
  public:
    typedef typename SearchTraits::FT FT;
    typedef Point_container<SearchTraits> Container;
    typedef Separator_ Separator;
    
    
    Median_of_rectangle()
      : Base()
    {}
    
    Median_of_rectangle(unsigned int bucket_size)
      : Base(bucket_size)
    {}

    void 
    operator() (Separator& sep, Container& c0, Container& c1) const
    {
      sep = Separator(c0.max_span_coord(),FT(0));
      sep.set_cutting_value(c0.median(sep.cutting_dimension()));
      c0.split(c1,sep,true);
    }
  };

  template <class SearchTraits,
	   class Separator_= Plane_separator<typename SearchTraits::FT>  >
  class Midpoint_of_max_spread
    : public Splitter_base<typename SearchTraits::FT> 
  {
    typedef Splitter_base<typename SearchTraits::FT> Base;
    
  public:
    typedef typename SearchTraits::FT FT;
    typedef Point_container<SearchTraits> Container;
    typedef Separator_ Separator;
    

    Midpoint_of_max_spread()
      : Base()
    {}
    
    Midpoint_of_max_spread(unsigned int bucket_size)
      : Base(bucket_size)
    {}
    
    void 
    operator()(Separator& sep, Container& c0, Container& c1) const
    {
      sep = Separator(c0.max_tight_span_coord(),
		      (c0.max_tight_span_upper() + c0.max_tight_span_lower())/FT(2));
      c0.split(c1,sep);
    }
  };
  
  template <class SearchTraits,
	   class Separator_= Plane_separator<typename SearchTraits::FT>  >
  class Midpoint_of_rectangle
  : public Splitter_base<typename SearchTraits::FT> 
  {

    typedef Splitter_base<typename SearchTraits::FT> Base;
  public:
    typedef typename SearchTraits::FT FT;
    typedef Point_container<SearchTraits> Container;
    typedef Separator_ Separator;
    
    
    Midpoint_of_rectangle()
      : Base()
    {}
    
    Midpoint_of_rectangle(unsigned int bucket_size)
      : Base(bucket_size)
    {}
    
    void 
    operator()(Separator& sep, Container& c0, Container& c1) const
    {
      sep = Separator(c0.max_span_coord(),
		      (c0.max_span_upper() + c0.max_span_lower())/FT(2));
      c0.split(c1,sep);          
    }
    
  };

} // namespace CGAL
#endif // CGAL_SPLITTERS
