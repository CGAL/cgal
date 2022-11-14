// Copyright (c) 2002-2011 Utrecht University (The Netherlands).
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

// Defines rules used for constructing a split node. That is, it implements,
// in several ways, the concept
// Boxtree_splitter<FT>.

#ifndef CGAL_SPLITTERS_H
#define CGAL_SPLITTERS_H

#include <CGAL/license/Spatial_searching.h>

#include <CGAL/Point_container.h>
#include <CGAL/Plane_separator.h>

namespace CGAL {

  template <class FT>
  class Splitter_base {
  private:
    unsigned int the_bucket_size;
    FT the_aspect_ratio;

  public:
          //default bucket_size should be 10
    Splitter_base(unsigned int bucket_size = 10,
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
          int cutdim = c0.max_span_coord();

          //Bugfix: avoid linear tree in degenerated cases
          if(c0.tight_bounding_box().min_coord(cutdim) != c0.tight_bounding_box().max_coord(cutdim)){
                  sep = Separator(cutdim,
                                  (c0.max_span_upper() + c0.max_span_lower())/FT(2));
          }
          else{
                  cutdim = c0.max_tight_span_coord();
                  sep = Separator(cutdim,
                                  (c0.max_tight_span_upper() + c0.max_tight_span_lower())/FT(2));
          }

      FT max_span_lower =
        c0.tight_bounding_box().min_coord(cutdim);

      CGAL_assertion(max_span_lower >= c0.bounding_box().min_coord(cutdim));
      FT max_span_upper =
        c0.tight_bounding_box().max_coord(cutdim);

      CGAL_assertion(max_span_upper <= c0.bounding_box().max_coord(cutdim));
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
