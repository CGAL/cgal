// Copyright (c) 2012 GeometryFactory (France). All rights reserved.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Philipp Moeller

#ifndef CGAL_SKIPLIST_H
#define CGAL_SKIPLIST_H

#include <boost/intrusive/list.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/iterator/iterator_adaptor.hpp>

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
    const T& get() const { return t_; }
    T& get() { return t_; }
  private:
    value_type t_;
  public:
    boost::intrusive::list_member_hook<> skip_hook;
    boost::intrusive::list_member_hook<> all_hook;
  };

  struct Node_disposer {
    void operator()(Node* p) const { delete p; }
  };

  typedef boost::intrusive::member_hook<
    Node,
    boost::intrusive::list_member_hook<>,
    &Node::skip_hook>
  SkipOption;
  typedef boost::intrusive::member_hook<
    Node,
    boost::intrusive::list_member_hook<>,
    &Node::all_hook>
  AllOption;

  typedef boost::intrusive::list<Node, SkipOption> skip_list;
  typedef boost::intrusive::list<Node, AllOption> all_list;
public:
  // assume both lists have the same size_type/ptr_diff
  typedef typename all_list::size_type size_type;
  typedef typename all_list::difference_type difference_type;

  struct all_iterator :
    public boost::iterator_adaptor<
      all_iterator
    , typename all_list::iterator
    , T
    >
  {
  public:
    all_iterator() {}
    all_iterator(typename all_list::iterator it)
      : all_iterator::iterator_adaptor_(it) {}
  private:
    friend class boost::iterator_core_access;
    T& dereference() const { return this->base()->get(); }
  };

  struct skip_iterator :
    public boost::iterator_adaptor<
      skip_iterator
    , typename skip_list::iterator
    , T
    >
  {
  public:
    skip_iterator() {}
    skip_iterator(typename skip_list::iterator it)
      : skip_iterator::iterator_adaptor_(it) {}
    operator all_iterator() const { return all_list::s_iterator_to(*this->base()); }
  private:
    friend class boost::iterator_core_access;
    T& dereference() const { return this->base()->get(); }
  };

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

  /// The semantics of front and back try to be consistent with
  /// push_back and push_front. Questionable.

  /// Returns front of the skiplist
  const_reference
  front() const { return skip_.front().get(); }

  /// Returns front of the skiplist
  reference
  front() { return skip_.front().get(); }

  /// Returns back of the skiplist
  const_reference
  back() const { return skip_.back().get(); }

  /// Returns back of the skiplist
  reference
  back() { return skip_.back().get(); }

  all_iterator all_begin()
  {
    return all_.begin();
  }
  all_iterator all_end()
  {
    return all_.end();
  }

  skip_iterator skip_begin()
  {
    return skip_.begin();
  }
  skip_iterator skip_end()
  {
    return skip_.end();
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
  void skip(skip_iterator it)
  {
    skip_.erase(it.base());
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

  /// Check if an all_iterator has been skipped.
  ///
  /// @param it a valid all_iterator
  ///
  /// @return true if the element pointed to by it is not in the range [skip_begin(), skip_end())
  bool is_skipped(all_iterator it) const
  {
    return !(it.base()->skip_hook.is_linked());
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

  void pop_back()
  {
    all_.pop_back();
    skip_.pop_back();
  }

  /// Insert \c t before \c pos in the all_view. \t will not be inserted into the skip view.
  /// @returns an skip_iterator to the inserted element.
  all_iterator insert(all_iterator pos, const value_type& t)
  {
    return all_.insert(pos.base(), *new Node(t));
  }

  /// Insert \c t before \c pos in the all_view. \t will be inserted into the skip view.
  /// @returns an iterator to the inserted element.
  skip_iterator insert(skip_iterator pos, const value_type& t)
  {
    all_iterator it = insert(static_cast<all_iterator>(pos), t);
    return skip_.insert(pos.base(), *it.base());
  }

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

  /// Drop the contents of iterator \c it from both views.
  all_iterator erase(all_iterator it)
  {
    if(!is_skipped(it)) {
      skip_.erase(skip_.iterator_to(*it.base()));
    }

    return all_.erase_and_dispose(it.base(), Node_disposer());
  }

  void splice(skip_iterator pos, Skiplist& other,
              skip_iterator first, skip_iterator last)
  {
    all_iterator alllast = last == other.skip_end() ?
      other.all_.end() : other.all_.iterator_to(*last.base());
    all_iterator allpos = pos == skip_end() ?
      all_end() : all_.iterator_to(*pos.base());

    all_.splice(allpos.base(), other.all_,
                other.all_.iterator_to(*first.base()), alllast.base());
    skip_.splice(pos.base(), other.skip_, first.base(), last.base());
  }

  size_type
  all_size() const { return all_.size(); }

  size_type
  skip_size() const { return skip_.size(); }

  bool empty() const { return all_.empty(); }

  void swap(Skiplist& other)
  {
    using std::swap;
    this->all_.swap(other.all_);
    this->skip_.swap(other.skip_);
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

template<typename T>
void swap(Skiplist<T>& a, Skiplist<T>& b)
{
  a.swap(b);
}

} // CGAL


#endif /* CGAL_SKIPLIST_H */
