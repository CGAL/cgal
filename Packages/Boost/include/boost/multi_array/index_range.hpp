// Copyright 2002 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software 
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Boost.MultiArray Library
//  Authors: Ronald Garcia
//           Jeremy Siek
//           Andrew Lumsdaine
//  See http://www.boost.org/libs/multi_array for documentation.

#ifndef BOOST_INDEX_RANGE_RG071801_HPP
#define BOOST_INDEX_RANGE_RG071801_HPP

#include <boost/config.hpp>
#include <utility>
#include <boost/limits.hpp>

// For representing intervals, also with stride.
// A degenerate range is a range with one element.

// Thanks to Doug Gregor for the really cool idea of using the
// comparison operators to express various interval types!

// Internally, we represent the interval as half-open.

namespace boost {
namespace detail {
namespace multi_array {

  template <typename Index,typename SizeType>
  class index_range {
  public:
    typedef Index index;
    typedef SizeType size_type;

    index_range()
    {
      start_ = from_start();
      finish_ = to_end();
      stride_ = 1;
      degenerate_ = false;
    }

    explicit index_range(index pos)
    {
      start_ = pos;
      finish_ = pos;
      stride_ = 1;
      degenerate_ = true;
    }

    explicit index_range(index start, index finish, index stride=1)
      : start_(start), finish_(finish), stride_(stride),
        degenerate_(false)
    { }


    // These are for chaining assignments to an index_range
    index_range& start(index s) {
      start_ = s;
      degenerate_ = (start_ == finish_);
      return *this;
    }

    index_range& finish(index f) {
      finish_ = f;
      degenerate_ = (start_ == finish_);
      return *this;
    }

    index_range& stride(index s) { stride_ = s; return *this; }

    index start() const
    { 
      return start_; 
    }

    index get_start(index low_index_range = 0) const
    { 
      if (start_ == from_start())
        return low_index_range;
      return start_; 
    }

    index finish() const
    {
      return finish_;
    }

    index get_finish(index high_index_range = 0) const
    {
      if (finish_ == to_end())
        return high_index_range;
      return finish_;
    }

    size_type size(index recommended_length = 0) const
    {
      if ((start_ == from_start()) || (finish_ == to_end()))
        return recommended_length;
      else 
        return (finish_ - start_) / stride_;
    }

    index stride() const { return stride_; }

    bool is_ascending_contiguous() const
    {
      return (start_ < finish_) && is_unit_stride();
    }

    void set_index_range(index start, index finish, index stride=1)
    {
      start_ = start;
      finish_ = finish;
      stride_ = stride;
    }

    static index_range all() 
    { return index_range(from_start(), to_end(), 1); }

    bool is_unit_stride() const
    { return stride_ == 1; }

    bool is_degenerate() const { return degenerate_; }

    index_range operator-(index shift) const
    { 
      return index_range(start_ - shift, finish_ - shift, stride_); 
    }

    index_range operator+(index shift) const
    { 
      return index_range(start_ + shift, finish_ + shift, stride_); 
    }

    index operator[](unsigned i) const
    {
      return start_ + i * stride_;
    }

    index operator()(unsigned i) const
    {
      return start_ + i * stride_;
    }

    // add conversion to std::slice?

  private:
    static index from_start()
      { return (std::numeric_limits<index>::min)(); }

    static index to_end()
      { return (std::numeric_limits<index>::max)(); }
  public:
    index start_, finish_, stride_;
    bool degenerate_;
  };

  // Express open and closed interval end-points using the comparison
  // operators.

  // left closed
  template <typename Index, typename SizeType>
  inline index_range<Index,SizeType>
  operator<=(Index s, const index_range<Index,SizeType>& r)
  {
    return index_range<Index,SizeType>(s, r.finish(), r.stride());
  }

  // left open
  template <typename Index, typename SizeType>
  inline index_range<Index,SizeType>
  operator<(Index s, const index_range<Index,SizeType>& r)
  {
    return index_range<Index,SizeType>(s + 1, r.finish(), r.stride());
  }

  // right open
  template <typename Index, typename SizeType>
  inline index_range<Index,SizeType>
  operator<(const index_range<Index,SizeType>& r, Index f)
  {
    return index_range<Index,SizeType>(r.start(), f, r.stride());
  }

  // right closed
  template <typename Index, typename SizeType>
  inline index_range<Index,SizeType>
  operator<=(const index_range<Index,SizeType>& r, Index f)
  {
    return index_range<Index,SizeType>(r.start(), f + 1, r.stride());
  }

} // namespace multi_array
} // namespace detail  
} // namespace boost

#endif // BOOST_INDEX_RANGE_RG071801_HPP
