// Copyright (c) 2012 GeometryFactory (France). All rights reserved.
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
// $URL: 
// $Id: 
// 
// Author(s)     : Philipp Moeller

#ifndef CGAL_SKIPLIST_H
#define CGAL_SKIPLIST_H

#include <utility>
#include <cassert>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/sequenced_index.hpp>
#include <boost/multi_index/random_access_index.hpp>

namespace CGAL {

/// The skiplist maintains two iterator ranges on the set of elements:
/// the all view and the skip view. The all view at all times contains
/// all elements in the Skiplist. The skip view can be modified and
/// items can be removed from it.
///
/// insert will add elements to the all_view.
/// 
/// @tparam T the value_type to store in the Skiplist
template<typename T>
class Skiplist
{
public:
  typedef T value_type;
  
private:
  typedef boost::multi_index_container<
    T,
    boost::multi_index::indexed_by<
      boost::multi_index::sequenced<> // all elements
      , boost::multi_index::random_access<> // skipping elements
      >
    > container_type;

  typedef typename container_type::template nth_index<0>::type all_index;
  typedef typename container_type::template nth_index<1>::type skip_index;
public:
  typedef typename all_index::iterator all_iterator;
  typedef std::pair<all_iterator, all_iterator> all_range;
  typedef typename skip_index::iterator skip_iterator;
  typedef std::pair<skip_iterator, skip_iterator> skip_range;

  Skiplist() : skip_end_(boost::get<1>(c).end()) {}

  /// Construct a Skiplist from the range [begin,end)
  /// @postcond Both views are equal to the range [begin,end)
  template<typename InputIterator>
  Skiplist(InputIterator begin, InputIterator end)
  {
    boost::get<0>(c).insert(all_begin(), begin, end);
    skip_end_ = boost::get<1>(c).end();
  }

  all_iterator all_begin() const { return boost::get<0>(c).begin(); }
  all_iterator all_end() const { return boost::get<0>(c).end(); }

  skip_iterator skip_begin() const { return boost::get<1>(c).begin(); }
  skip_iterator skip_end() const { return skip_end_; }

  all_range
  all_elements() const { return std::make_pair(all_begin(), all_end()); }

  skip_range
  skip_elements() const { return std::make_pair(skip_begin(), skip_end()); }

  /// The elements pointed to by it are no longer in the range
  /// [skip_begin(), skip_end()).
  ///
  /// @precond it and to_skip(it) are valid iterators
  ///
  void skip(all_iterator it)
  {
    skip_iterator sit = c.template project<1>(it);
    boost::get<1>(c).relocate(skip_end_, sit);
    skip_end_ = sit;
  }
  
  void unskip(skip_iterator pos, all_iterator it)
  {
    skip_iterator m = to_skip(it);
    if(m == skip_end_)
      ++skip_end_;
    boost::get<1>(c).relocate(pos, m);
  }

  /// The elements pointed to by [begin, end) are no longer in the
  /// range [skip_begin(), skip_end()). If an element in [begin, end)
  /// already is removed from the skip view, it will stay removed.
  ///
  /// @precond [begin,end) is a slice of the range [all_begin(), all_end())
  ///
  void skip(all_iterator begin, all_iterator end)
  {
    if(end == all_end()) {
      skip_end_ = c.template project<1>(begin);
    } else {
      skip_iterator sbegin = c.template project<1>(begin);
      skip_iterator send = c.template project<1>(end);
      boost::get<1>(c).relocate(skip_end_, sbegin, send);
      skip_end_ = sbegin;
    }
  }

  
  /// Convert the argument from an skip_iterator to the corresponding
  /// all_iterator.
  ///
  /// @precond is_skipped(it) == false
  all_iterator to_all(skip_iterator it) const
  {
    // special treatment for the end
    if(it == skip_end_) {
      return all_end();
    }
    return c.template project<0>(it);
  }

  /// Convert the argument from an all_iterator to the corresponding
  /// skip_iterator.
  ///
  /// @precond is_skipped(it) == false and it is a valid dereferencable iterator in the range [all_begin(), all_end())
  skip_iterator to_skip(all_iterator it) const
  {
    return c.template project<1>(it);
  }

  /// Check if an all_iterator has been skipped.
  /// 
  /// @param it a valid all_iterator
  ///
  /// @return true if the element pointed to by it is not in the range [skip_begin(), skip_end())
  bool is_skipped(all_iterator it) const
  {
    return skip_end_ <= c.template project<1>(it);
  }

  /// Adds an element to the end of both views in Skiplist.
  void push_back(const value_type& t)
  {
    insert(all_end(), t);
    unskip(skip_end(), boost::prior(all_end()));
  }

  /// Adds an element to the front of both views in Skiplist.
  void push_front(const value_type& t)
  {
    insert(all_begin(), t);
    unskip(skip_begin(), all_begin());
  }

  /// Insert \c t before \c pos in the all_view. \t will not be inserted into the skip view.
  /// @returns an iterator to the inserted element.
  all_iterator insert(all_iterator pos, const value_type& t)
  {
    std::pair<all_iterator, bool> p = boost::get<0>(c).insert(pos, t);
    assert(p.second);
    // we are in the non-skipped state and need to readjust the end
    if(skip_end_ == boost::get<1>(c).end())
      --skip_end_;
    
    return p.first;
  }


  /// Insert \c t before \c pos in the all_view. \t will be inserted into the skip view.
  /// @returns an iterator to the inserted element.
  skip_iterator insert(skip_iterator pos, const value_type& t)
  {
    std::pair<skip_iterator, bool> p = boost::get<1>(c).insert(pos, t);
    assert(p.second);
    boost::get<0>(c).relocate(to_all(pos), boost::prior(all_end()));
    return p.first;
  }

  /// Insert the range [begin,end) into the all view. If the container
  /// is empty() the range will also be visible in the skip view.
  template<typename InputIterator>
  void insert(all_iterator pos, InputIterator begin, InputIterator end)
  {
    // we are in the non-skipped state and need to readjust the end
    if(skip_end_ == boost::get<1>(c).end()) {
      if(skip_size() != 0) {
        --skip_end_;
        boost::get<0>(c).insert(pos, begin, end);
      } else {
        boost::get<0>(c).insert(pos, begin, end);
        skip_end_ = boost::get<1>(c).end();
      }
    } else {
      boost::get<0>(c).insert(pos, begin, end);
    }
  }

  typename all_index::size_type 
  all_size() const { return  boost::get<0>(c).size(); }

  typename skip_index::size_type 
  skip_size() const { return std::distance(skip_begin(), skip_end()); }
  
  bool empty() const { return c.empty(); }

  friend void swap(Skiplist& a, Skiplist& b) {
    std::swap(a.c, b.c);
    std::swap(a.skip_end_, b.skip_end_);
  }

  /// Reset the container.
  /// @postcond *this.empty() == true
  void clear()
  {
    c.clear();
    skip_end_ = boost::get<1>(c).end();
  }

  template<typename Modifier>
  bool modify(all_iterator it, Modifier m)
  {
    return boost::get<0>(c).modify(it, m);
  }

  template<typename Modifier>
  bool modify(skip_iterator it, Modifier m)
  {
    return boost::get<1>(c).modify(it, m);
  }

private:
  container_type c;
  skip_iterator skip_end_;
};
} // CGAL


#endif /* CGAL_SKIPLIST_H */
