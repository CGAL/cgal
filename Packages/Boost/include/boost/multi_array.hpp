// Copyright 2002 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software 
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Boost.MultiArray Library
//  Authors: Ronald Garcia
//           Jeremy Siek
//           Andrew Lumsdaine
//  See http://www.boost.org/libs/multi_array for documentation.

#ifndef BOOST_MULTI_ARRAY_RG071801_HPP
#define BOOST_MULTI_ARRAY_RG071801_HPP

//
// multi_array.hpp - contains the multi_array class template
// declaration and definition
//

#include "boost/multi_array/base.hpp"
#include "boost/multi_array/collection_concept.hpp"
#include "boost/multi_array/copy_array.hpp"
#include "boost/multi_array/iterator.hpp"
#include "boost/multi_array/subarray.hpp"
#include "boost/multi_array/multi_array_ref.hpp"
#include "boost/multi_array/algorithm.hpp"
#include "boost/array.hpp"
#include "boost/type_traits.hpp"
#include <algorithm>
#include <cstddef>
#include <functional>
#include <numeric>
#include <vector>



namespace boost {
  namespace detail {
    namespace multi_array {
      struct populate_index_ranges {
        multi_array_types::index_range
        operator()(multi_array_types::index base,
                   multi_array_types::size_type extent) {
          return multi_array_types::index_range(base,base+extent);
        }
      };
    } //namespace multi_array
  } // namespace detail

template<typename T, std::size_t NumDims,
  typename Allocator>
class multi_array :
  public multi_array_ref<T,NumDims>
{
  typedef multi_array_ref<T,NumDims> super_type;
public:
  typedef typename super_type::value_type value_type;
  typedef typename super_type::reference reference;
  typedef typename super_type::const_reference const_reference;
  typedef typename super_type::iterator iterator;
  typedef typename super_type::const_iterator const_iterator;
  typedef typename super_type::reverse_iterator reverse_iterator;
  typedef typename super_type::const_reverse_iterator const_reverse_iterator;
  typedef typename super_type::element element;
  typedef typename super_type::size_type size_type;
  typedef typename super_type::difference_type difference_type;
  typedef typename super_type::index index;
  typedef typename super_type::extent_range extent_range;


  template <std::size_t NDims>
  struct const_array_view {
    typedef boost::detail::multi_array::const_multi_array_view<T,NDims> type;
  };

  template <std::size_t NDims>
  struct array_view {
    typedef boost::detail::multi_array::multi_array_view<T,NDims> type;
  };

  explicit multi_array() :
    super_type((T*)initial_base_) {
    allocate_space();
  }
    
  template <class ExtentList>
  explicit multi_array(
      ExtentList const& extents
#ifdef BOOST_NO_FUNCTION_TEMPLATE_ORDERING
    , typename detail::multi_array::disable_non_sub_array<ExtentList>::type* = 0
#endif 
  ) :
    super_type((T*)initial_base_,extents) {
    boost::function_requires<
      detail::multi_array::CollectionConcept<ExtentList> >();
    allocate_space();
  }
    
  template <class ExtentList>
  explicit multi_array(ExtentList const& extents,
                       const general_storage_order<NumDims>& so) :
    super_type((T*)initial_base_,extents,so) {
    boost::function_requires<
      detail::multi_array::CollectionConcept<ExtentList> >();
    allocate_space();
  }

  template <class ExtentList>
  explicit multi_array(ExtentList const& extents,
                       const general_storage_order<NumDims>& so,
                       Allocator const& alloc) :
    super_type((T*)initial_base_,extents,so), allocator_(alloc) {
    boost::function_requires<
      detail::multi_array::CollectionConcept<ExtentList> >();
    allocate_space();
  }


  explicit multi_array(const detail::multi_array
                       ::extent_gen<NumDims>& ranges) :
    super_type((T*)initial_base_,ranges) {

    allocate_space();
  }


  explicit multi_array(const detail::multi_array
                       ::extent_gen<NumDims>& ranges,
                       const general_storage_order<NumDims>& so) :
    super_type((T*)initial_base_,ranges,so) {

    allocate_space();
  }


  explicit multi_array(const detail::multi_array
                       ::extent_gen<NumDims>& ranges,
                       const general_storage_order<NumDims>& so,
                       Allocator const& alloc) :
    super_type((T*)initial_base_,ranges,so), allocator_(alloc) {

    allocate_space();
  }

  multi_array(const multi_array& rhs) :
  super_type(rhs), allocator_(rhs.allocator_) {
    allocate_space();
    boost::copy_n(rhs.base_,rhs.num_elements(),base_);
  }

  template <typename OPtr>
  multi_array(const detail::multi_array::
              const_sub_array<T,NumDims,OPtr>& rhs) :
    super_type(rhs) {
    allocate_space();
    std::copy(rhs.begin(),rhs.end(),this->begin());
  }

  // For some reason, gcc 2.95.2 doesn't pick the above template
  // member function when passed a subarray, so i was forced to
  // duplicate the functionality here...
  multi_array(const detail::multi_array::
              sub_array<T,NumDims>& rhs) :
    super_type(rhs) {
    allocate_space();
    std::copy(rhs.begin(),rhs.end(),this->begin());
  }
    
  // Since assignment is a deep copy, multi_array_ref
  // contains all the necessary code.
  template <typename ConstMultiArray>
  multi_array& operator=(const ConstMultiArray& other) {
    super_type::operator=(other);
    return *this;
  }

  multi_array& operator=(const multi_array& other) {
    if (&other != this) {
      super_type::operator=(other);
    }
    return *this;
  }


  multi_array& resize(const detail::multi_array
                      ::extent_gen<NumDims>& ranges) {


    // build a multi_array with the specs given
    multi_array new_array(ranges);


    // build a view of tmp with the minimum extents

    // Get the minimum extents of the arrays.
    boost::array<size_type,NumDims> min_extents;

    const size_type& (*min)(const size_type&, const size_type&) =
      std::min;
    std::transform(new_array.extent_list_.begin(),new_array.extent_list_.end(),
                   this->extent_list_.begin(),
                   min_extents.begin(),
                   min);


    // typedef boost::array<index,NumDims> index_list;
    // Build index_gen objects to create views with the same shape

    // these need to be separate to handle non-zero index bases
    typedef detail::multi_array::index_gen<NumDims,NumDims> index_gen;
    index_gen old_idxes;
    index_gen new_idxes;

    std::transform(new_array.index_base_list_.begin(),
                   new_array.index_base_list_.end(),
                   min_extents.begin(),old_idxes.ranges_.begin(),
                   detail::multi_array::populate_index_ranges());

    std::transform(this->index_base_list_.begin(),
                   this->index_base_list_.end(),
                   min_extents.begin(),new_idxes.ranges_.begin(),
                   detail::multi_array::populate_index_ranges());

    // Build same-shape views of the two arrays
    typename
      multi_array::BOOST_NESTED_TEMPLATE array_view<NumDims>::type view_old = (*this)[old_idxes];
    typename
      multi_array::BOOST_NESTED_TEMPLATE array_view<NumDims>::type view_new = new_array[new_idxes];

    // Set the right portion of the new array
    view_new = view_old;

    using std::swap;
    // Swap the internals of these arrays.
    swap(this->super_type::base_,new_array.super_type::base_);
    swap(this->storage_,new_array.storage_);
    swap(this->extent_list_,new_array.extent_list_);
    swap(this->stride_list_,new_array.stride_list_);
    swap(this->index_base_list_,new_array.index_base_list_);
    swap(this->origin_offset_,new_array.origin_offset_);
    swap(this->directional_offset_,new_array.directional_offset_);
    swap(this->num_elements_,new_array.num_elements_);
    swap(this->allocator_,new_array.allocator_);
    swap(this->base_,new_array.base_);
    swap(this->allocated_elements_,new_array.allocated_elements_);

    return *this;
  }


  ~multi_array() {
    deallocate_space();
  }

private:
  void allocate_space() {
    typename Allocator::const_pointer no_hint=0;
    base_ = allocator_.allocate(this->num_elements(),no_hint);
    this->set_base_ptr(base_);
    allocated_elements_ = this->num_elements();
    std::uninitialized_fill_n(base_,allocated_elements_,T());
  }

  void deallocate_space() {
    if(base_) {
      for(T* i = base_; i != base_+allocated_elements_; ++i)
        allocator_.destroy(i);
      allocator_.deallocate(base_,allocated_elements_);
    }
  }

  typedef boost::array<size_type,NumDims> size_list;
  typedef boost::array<index,NumDims> index_list;

  Allocator allocator_;
  T* base_;
  size_type allocated_elements_;
  enum {initial_base_ = 0};
};

} // namespace boost

#endif // BOOST_MULTI_ARRAY_RG071801_HPP
