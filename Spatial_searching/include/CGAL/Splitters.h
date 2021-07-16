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
    using Ref_pair = typename Container::Ref_pair;
    using Point_d = typename SearchTraits::Point_d;

    Balanced_splitter() : Base() { }
    Balanced_splitter(const unsigned int bucket_size) : Base(bucket_size) { }

    struct Balanced_cmp {

      const long m_median;
      const std::size_t m_dim;
      const std::size_t m_axis;
      const std::vector<Ref_pair>& m_reference;
      const Key_compare compare_keys;

      Balanced_cmp(
        const std::vector<Ref_pair>& reference,
        const std::size_t axis, const long median, const std::size_t dim) :
      m_median(median), m_dim(dim), m_axis(axis), m_reference(reference)
      { }

      // TODO: Can we avoid this linear search?
      bool operator()(const Point_d* pt) const {

        FTP p; bool duplicates_found = true;
        for (const auto& ref_pair : m_reference) {
          if (ref_pair.second == pt) {
            p = ref_pair.first;
            duplicates_found = false;
            break;
          }
        }
        if (duplicates_found) {
          CGAL_assertion_msg(false, "ERROR: DUPLICATES ARE FOUND!");
          return false;
        }

        const FT compare = compare_keys(
          p, m_reference[m_median].first, m_axis, m_dim);
        return compare < FT(0);
      }
    };

    template<typename References>
    void print_references(const References& references) const {

      for (const auto& reference : (*references)) {
        std::cout << std::endl;
        for (const auto& ref_pair : reference) {
          std::cout << *(ref_pair.first) << std::endl;
        }
      }
    }

    void balanced_split(Container& c0, Container& c1, Separator& sep) const {

      CGAL_assertion(c0.dimension() == c1.dimension());
      CGAL_assertion(c0.is_valid());

      const std::size_t dim = c0.get_dim();
      const std::size_t depth = c0.get_depth();

      const long start = c0.get_start();
      const long end = c0.get_end();

      const auto& references = c0.get_references();

      // NEW CODE!
      CGAL_assertion(dim   != static_cast<std::size_t>(-1));
      CGAL_assertion(depth != static_cast<std::size_t>(-1));

      CGAL_assertion(start != -1);
      CGAL_assertion(end   != -1);

      const long median = start + ((end - start) / 2);
      const std::size_t axis = depth % dim;
      CGAL_assertion(median >= 0);

      // std::cout << "data: "   << c0.size() << std::endl;
      // std::cout << "start: "  << start     << std::endl;
      // std::cout << "end: "    << end       << std::endl;
      // std::cout << "median: " << median    << std::endl;

      const Key_compare compare_keys;
      std::vector<FTP> tmp((*references)[0].size());
      CGAL_assertion(end > start + 2);
      for (long i = start; i <= end; ++i) { // copy the first ref to save it
        tmp[i] = (*references)[0][i].first;
      }

      long lower, upper;
      for (std::size_t i = 1; i < dim; ++i) { // apply permutations for dim - 1 axis
        lower = start - 1;
        upper = median;

        for (long j = start; j <= end; ++j) { // always skip median key
          const FT compare = compare_keys(
            (*references)[i][j].first, (*references)[0][median].first, axis, dim);
          if (compare < FT(0)) {
            (*references)[i-1][++lower].first = (*references)[i][j].first;
          } else if (compare > FT(0)) {
            (*references)[i-1][++upper].first = (*references)[i][j].first;
          }
        }
      }
      CGAL_assertion(lower >= 0 && upper >= 0);

      for (long i = start; i <= end; ++i) { // copy back
        (*references)[dim-1][i].first = tmp[i];
      }

      // std::cout << "DEPTH: " << depth << std::endl;
      // print_references(references);

      ///////////////////////

      // OLD CODE!

      // c1.bounding_box() = c0.boundig_box();
      // const int split_coord = static_cast<int>(axis); // TODO: Is it correct?
      // const FT cutting_value = ((*references)[dim-1][median].first)[axis]; // TODO: Is it correct?

      // std::cout << "split coord: "   << split_coord   << std::endl;
      // std::cout << "cutting value: " << cutting_value << std::endl;

      // const int split_coord = sep.cutting_dimension();
      // const FT cutting_value = sep.cutting_value();

      // c0.set_built_coordinate(split_coord);
      // c1.set_built_coordinate(split_coord);

      // auto construct_it = traits.construct_cartesian_const_iterator_d_object();
      // Cmp<Traits> cmp(split_coord, cutting_value, construct_it);
      // const iterator it = std::partition(c0.begin(), c0.end(), cmp);

      ///////////////////////

      // NEW VERSION!
      const Balanced_cmp cmp((*references)[dim-1], axis, median, dim);
      const auto it = std::partition(c0.begin(), c0.end(), cmp);

      c1.set_range(c0.begin(), it);
      c0.set_range(it, c0.end());

      ///////////////////////

      // std::cout << "pre size c0: " << c0.size() << std::endl;
      // std::cout << "pre size c1: " << c1.size() << std::endl;

      c1.set_references(references);
      c1.set_data(dim, depth + 1, start, lower);
      c0.set_data(dim, depth + 1, median + 1, upper);

      // Adjust boxes and other data.
      // TODO: The cutting value and split_coord obtained from the new algorithm
      // do not satisfy the is_valid() criteria. Maybe it is correct and due to
      // a different type of partition but maybe we do something wrong. For the moment,
      // we recompute these values for tight bboxes and use midpoint as ref. The question is:
      // does it affect the search anyhow?

      SearchTraits traits;
      auto construct_it = traits.construct_cartesian_const_iterator_d_object();
      // c0.bounding_box().set_lower_bound(split_coord, cutting_value); // old version
      c0.bounding_box(). template update_from_point_pointers<
      typename SearchTraits::Construct_cartesian_const_iterator_d>(c0.begin(), c0.end(), construct_it);
      c0.tight_bounding_box(). template update_from_point_pointers<
      typename SearchTraits::Construct_cartesian_const_iterator_d>(c0.begin(), c0.end(), construct_it);

      // c1.bounding_box().set_upper_bound(split_coord, cutting_value); // old version
      c1.bounding_box(). template update_from_point_pointers<
      typename SearchTraits::Construct_cartesian_const_iterator_d>(c1.begin(), c1.end(), construct_it);
      c1.tight_bounding_box(). template update_from_point_pointers<
      typename SearchTraits::Construct_cartesian_const_iterator_d>(c1.begin(), c1.end(), construct_it);

      c0.set_built_coordinate(c0.max_span_coord());
      c1.set_built_coordinate(c0.max_span_coord());

      sep.set_cutting_dimension(c0.max_span_coord());
      const FT cutting_value = (c0.max_span_upper() + c0.max_span_lower()) / FT(2); // midpoint
      // const FT cutting_value = c0.median(sep.cutting_dimension()); // median
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

      // CGAL_assertion_msg(false, "TODO: FINISH BALANCED SPLIT!");
    }

    void operator()(
      Separator& sep, Container& c0, Container& c1) const {

      CGAL_assertion(c0.is_valid());
      CGAL_assertion(c1.is_valid());

      this->balanced_split(c0, c1, sep);
      // CGAL_assertion_msg(false, "TODO: FINISH BALANCED SPLITTER!");
    }
  };

} // namespace CGAL
#endif // CGAL_SPLITTERS
