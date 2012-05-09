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

#include <boost/intrusive/list.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/iterator/transform_iterator.hpp>

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
  typedef T        value_type;
  typedef T&       reference;
  typedef const T& const_reference;
  typedef T*       pointer;
  typedef const T* const_pointer;
  
private:
  struct Node {
    explicit Node(const T& t) : t_(t) {}
    value_type t_;
    boost::intrusive::list_member_hook<> skip_hook_;
    boost::intrusive::list_member_hook<> all_hook_;
  };

  struct Node_disposer {
    void operator()(Node* p) const { delete p; }
  };

  static value_type& extract(Node& n) { return n.t_; }
  static const value_type& extract_const(const Node& n) { return n.t_; }

  typedef boost::intrusive::member_hook<
    Node,
    boost::intrusive::list_member_hook<>,
    &Node::skip_hook_>
  SkipOption;
  typedef boost::intrusive::member_hook<
    Node,
    boost::intrusive::list_member_hook<>,
    &Node::all_hook_>
  AllOption;

  typedef boost::intrusive::list<Node, SkipOption> skip_list;
  typedef boost::intrusive::list<Node, AllOption> all_list;
public:
  // assume both lists have the same size_type/ptr_diff
  typedef typename all_list::size_type size_type;
  typedef typename all_list::difference_type difference_type;

  typedef boost::transform_iterator<
    value_type& (*)(Node&)
    , typename all_list::iterator > all_iterator;
  typedef boost::transform_iterator<
    value_type& (*)(Node&)
    , typename skip_list::iterator > skip_iterator;
  
  typedef boost::iterator_range<all_iterator>  all_range;
  typedef boost::iterator_range<skip_iterator> skip_range;

  Skiplist() {}
  ~Skiplist() 
  {
    clear();
  }

  /// Construct a Skiplist from the range [begin,end)
  /// @postcond Both views are equal to the range [begin,end)
  template<typename InputIterator>
  Skiplist(InputIterator begin, InputIterator end)
  {
    while(begin != end) {
      push_back(*begin++);
    }
  }

  all_iterator all_begin() 
  { 
    return boost::make_transform_iterator(all_.begin(), &extract);
  }
  all_iterator all_end() 
  { 
    return  boost::make_transform_iterator(all_.end(), &extract);
  }

  skip_iterator skip_begin()
  {
    return boost::make_transform_iterator(skip_.begin(), &extract);
  }
  skip_iterator skip_end()
  {
    return  boost::make_transform_iterator(skip_.end(), &extract);
  }

  all_range
  all_elements()
  {
    return boost::make_iterator_range(all_begin(), all_end());
  }

  skip_range
  skip_elements()
  {
    return boost::make_iterator_range(skip_begin(), skip_end());
  }

  /// The elements pointed to by it are no longer in the range
  /// [skip_begin(), skip_end()).
  ///
  /// @precond it and to_skip(it) are valid iterators
  ///
  void skip(all_iterator it)
  {
    skip_.erase(skip_.iterator_to((*it.base())));
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
      skip_.erase(skip_.iterator_to(*begin.base()), skip_.end());
    } else {
      skip_.erase(skip_.iterator_to(*begin.base()), skip_.iterator_to(*end.base()));
    }
  }

  void unskip(skip_iterator pos, all_iterator it)
  {
    skip_.insert(pos.base(), *(it.base()));
  }

  /// Convert the argument from an skip_iterator to the corresponding
  /// all_iterator.
  ///
  /// @precond is_skipped(it) == false
  all_iterator to_all(skip_iterator it)
  {
    if(it == skip_end()) {
      return all_end();
    } else {
      return boost::make_transform_iterator(all_.iterator_to(*(it.base())), &extract);
    }
  }

  /// Convert the argument from an all_iterator to the corresponding
  /// skip_iterator.
  ///
  /// @precond is_skipped(it) == false and it is a valid dereferencable iterator in the range [all_begin(), all_end())
  skip_iterator to_skip(all_iterator it) const
  {
    if(it == all_end()) {
      return all_end();
    } else {
      return boost::make_transform_iterator(skip_.iterator_to(*(it.base())), &extract);
    }
  }

  /// Check if an all_iterator has been skipped.
  ///
  /// @param it a valid all_iterator
  ///
  /// @return true if the element pointed to by it is not in the range [skip_begin(), skip_end())
  bool is_skipped(all_iterator it) const
  {
    return !(it.base()->skip_hook_.is_linked());
  }

  /// Adds an element to the end of both views in Skiplist.
  void push_back(const value_type& t)
  {
    all_.push_back(*new Node(t));
    skip_.push_back(all_.back());
  }

  /// Adds an element to the front of both views in Skiplist.
  void push_front(const value_type& t)
  {
    all_.push_front(*new Node(t));
    skip_.push_front(all_.front());
  }

  /// Insert \c t before \c pos in the all_view. \t will not be inserted into the skip view.
  /// @returns an iterator to the inserted element.
  all_iterator insert(all_iterator pos, const value_type& t)
  {
    return boost::make_transform_iterator(all_.insert(pos.base(), *new Node(t)), &extract);
  }

  // /// Insert \c t before \c pos in the all_view. \t will be inserted into the skip view.
  // /// @returns an iterator to the inserted element.
  // skip_iterator insert(skip_iterator pos, const value_type& t)
  // {
  // }

  /// Insert the range [begin,end) into the all view. If the container
  /// is empty() the range will also be visible in the skip view.
  template<typename InputIterator>
  void insert(all_iterator pos, InputIterator begin, InputIterator end)
  {
    if(all_.empty()) {
      while(begin != end) {
        push_back(*begin++);
      }
    } else {
      while(begin != end) {
        pos = insert(pos, *begin++);
      }
    }
  }

  size_type
  all_size() const { return all_.size(); }

  size_type
  skip_size() const { return skip_.size(); }

  bool empty() const { return all_.empty(); }

  friend void swap(Skiplist& a, Skiplist& b) {
  }

  /// Reset the container.
  /// @postcond *this.empty() == true
  void clear()
  {
    skip_.clear();
    all_.clear_and_dispose(Node_disposer());
  }

private:
  all_list all_;
  skip_list skip_;
};
} // CGAL


#endif /* CGAL_SKIPLIST_H */
