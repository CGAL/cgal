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

  template<
  class SearchTraits,
  class Separator_ = Plane_separator<typename SearchTraits::FT> >
  class Balanced_splitter : public Splitter_base<typename SearchTraits::FT> {
    using Base = Splitter_base<typename SearchTraits::FT>;

  public:
    using FT = typename SearchTraits::FT;
    using Container = Point_container<SearchTraits>;
    using Separator = Separator_;
    using Key_compare = typename Container::Key_compare;
    using FTP = typename Container::FTP;
    using Point_d = typename SearchTraits::Point_d;

    Balanced_splitter() : Base() { }
    Balanced_splitter(const unsigned int bucket_size) : Base(bucket_size) { }

    void print_reference(
      const std::size_t dim,
      const std::vector<FTP>& reference) const {

      std::cout << std::endl;
      for (std::size_t i = 0; i < reference.size(); ++i) {
        std::cout << *(reference[i]) << std::endl;
        for (std::size_t k = 0; k < dim; ++k) {
          std::cout << reference[i][k] << " ";
        }
        std::cout << std::endl;
      }
    }

    void print_references(
      const std::size_t dim,
      const std::vector< std::vector<FTP> >& references,
      const std::string name) const {

      std::cout << std::endl << "--- " << name << " : " << std::endl;
      for (const auto& reference : references) {
        print_reference(dim, reference);
      }
    }

    void balanced_split(Container& c0, Container& c1, Separator& sep) const {

      CGAL_assertion(!c0.is_last_call());
      CGAL_assertion(!c1.is_last_call());

      CGAL_assertion(c0.dimension() == c1.dimension());
      CGAL_assertion(c0.is_valid());

      const std::size_t dim = c0.get_dim();
      const std::size_t depth = c0.get_depth();

      const long start = c0.get_start();
      const long end = c0.get_end();

      const auto& ref_ptr = c0.get_ref_ptr();
      auto& references = *(ref_ptr);

      // std::cout << "dim: "   << dim   << std::endl;
      // std::cout << "depth: " << depth << std::endl;
      // std::cout << "start: " << start << std::endl;
      // std::cout << "end: "   << end   << std::endl;

      CGAL_assertion(dim   != static_cast<std::size_t>(-1));
      CGAL_assertion(depth != static_cast<std::size_t>(-1));

      CGAL_assertion(start != -1);
      CGAL_assertion(end   != -1);

      const long median = start + ((end - start) / 2);
      const std::size_t axis = depth % dim;
      CGAL_assertion(median >= 0);

      // std::cout << "median: " << median << std::endl;

      const Key_compare compare_keys;
      std::vector<FTP> tmp(references[0].size());
      CGAL_assertion_msg(end > start + 2, "ERROR: LAST CALL MUST BE AVOIDED!");

      // Copy the first ref to save it.
      for (long i = start; i <= end; ++i) {
        tmp[i] = references[0][i];
      }

      // Apply permutations for dim - 1 axis.
      long lower, upper;
      for (std::size_t k = 1; k < dim; ++k) {
        lower = start - 1;
        upper = median;

        for (long i = start; i <= end; ++i) { // always skip median key
          const FT compare = compare_keys(
            references[k][i], references[0][median], axis, dim);
          if (compare < FT(0)) {
            references[k-1][++lower] = references[k][i];
          } else if (compare > FT(0)) {
            references[k-1][++upper] = references[k][i];
          }
        }
      }
      CGAL_assertion(lower >= 0 && upper >= 0);

      // Copy back.
      for (long i = start; i <= end; ++i) {
        references[dim-1][i] = tmp[i];
      }

      // print_references(dim, references, "DEPTH " + std::to_string(depth));

      SearchTraits traits;
      auto construct_it = traits.construct_cartesian_const_iterator_d_object();
      const auto it = std::partition(c0.begin(), c0.end(),
        [&](const Point_d* pt) {
          const auto p = construct_it(*pt);
          const auto q = references[dim-1][median];
          const FT compare = compare_keys(p, q, axis, dim);
          return compare < FT(0);
        }
      );

      c1.set_range(c0.begin(), it);
      c0.set_range(it, c0.end());

      CGAL_assertion(c0.size() > 0);
      CGAL_assertion(c1.size() > 0);

      // std::cout << "c0 size: " << c0.size() << std::endl;
      // std::cout << "c1 size: " << c1.size() << std::endl;

      // std::cout << std::endl << "c0 container, depth: " << depth << std::endl;
      // for (const auto& item : c0) {
      //   std::cout << *item << std::endl;
      // }
      // std::cout << std::endl << "c1 container, depth: " << depth << std::endl;
      // for (const auto& item : c1) {
      //   std::cout << *item << std::endl;
      // }

      c1.set_ref_ptr(ref_ptr);
      c1.set_data(dim, depth + 1, start, lower);
      c0.set_data(dim, depth + 1, median + 1, upper);

      // Adjust boxes and other data.
      // TODO: The cutting value and split_coord obtained from the new algorithm
      // do not satisfy the is_valid() criteria. Maybe it is correct and due to
      // a different type of partition but maybe we do something wrong. For the moment,
      // we recompute these values for tight bboxes and use midpoint as ref. The question is:
      // does it affect the search anyhow?

      c0.bounding_box(). template update_from_point_pointers<
      typename SearchTraits::Construct_cartesian_const_iterator_d>(c0.begin(), c0.end(), construct_it);
      c0.tight_bounding_box(). template update_from_point_pointers<
      typename SearchTraits::Construct_cartesian_const_iterator_d>(c0.begin(), c0.end(), construct_it);

      c1.bounding_box(). template update_from_point_pointers<
      typename SearchTraits::Construct_cartesian_const_iterator_d>(c1.begin(), c1.end(), construct_it);
      c1.tight_bounding_box(). template update_from_point_pointers<
      typename SearchTraits::Construct_cartesian_const_iterator_d>(c1.begin(), c1.end(), construct_it);

      c0.set_built_coordinate(c0.max_span_coord());
      c1.set_built_coordinate(c0.max_span_coord());

      sep.set_cutting_dimension(c0.max_span_coord());
      const FT cutting_value = (c0.max_span_upper() + c0.max_span_lower()) / FT(2); // midpoint
      sep.set_cutting_value(cutting_value);

      // std::cout << "c0.bbox: " << std::endl;
      // c0.bounding_box().print(std::cout);
      // std::cout << "c0.tbox: " << std::endl;
      // c0.tight_bounding_box().print(std::cout);
      // std::cout << std::endl;

      // std::cout << "c1.bbox: " << std::endl;
      // c1.bounding_box().print(std::cout);
      // std::cout << "c1.tbox: " << std::endl;
      // c1.tight_bounding_box().print(std::cout);
      // std::cout << std::endl;

      CGAL_assertion(c0.is_valid());
      CGAL_assertion(c1.is_valid());
    }

    void operator()(
      Separator& sep, Container& c0, Container& c1) const {

      CGAL_assertion(c0.is_valid());
      CGAL_assertion(c1.is_valid());
      this->balanced_split(c0, c1, sep);
    }
  };

} // namespace CGAL
#endif // CGAL_SPLITTERS
