// Copyright (c) 2014  INRIA Sophia-Antipolis (France)
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
// $URL: $
// $Id: $
//
//
// Author(s)     : Clement Jamin

#ifndef CGAL_TC_UTILITIES_H
#define CGAL_TC_UTILITIES_H

#include <CGAL/basic.h>
#include <CGAL/Dimension.h>

#include <set>
#include <vector>
#include <atomic> // CJTODO: this is C++11 => use boost.Atomic (but it's too recent) 
                  // or tbb::atomic (works for doubles, but not officially)

namespace CGAL {
namespace Tangential_complex_ {

  // Provides copy constructors to std::atomic so that
  // it can be used in a vector
  template <typename T>
  struct Atomic_wrapper 
    : public std::atomic<T>
  {
    typedef std::atomic<T> Base;

    Atomic_wrapper() {}
    Atomic_wrapper(const T &t) : Base(t) {}
    Atomic_wrapper(const std::atomic<T> &a) : Base(a.load()) {}
    Atomic_wrapper(const Atomic_wrapper &other) : Base(other.load())
    {}
    
    Atomic_wrapper &operator=(const T &other)
    {
      Base::store(other);
      return *this;
    }
    Atomic_wrapper &operator=(const std::atomic<T> &other)
    {
      Base::store(other.load());
      return *this;
    }
    Atomic_wrapper &operator=(const Atomic_wrapper &other)
    {
      Base::store(other.load());
      return *this;
    }
  };

  /*template <typename T>
  struct Atomic_wrapper
  {
    std::atomic<T> _a;

    Atomic_wrapper()
      :_a()
    {}

    Atomic_wrapper(const std::atomic<T> &other)
      :_a(other.load())
    {}

    Atomic_wrapper(const Atomic_wrapper &other)
      :_a(other._a.load())
    {}

    Atomic_wrapper(const T &other)
      :_a(other)
    {}

    Atomic_wrapper &operator=(const std::atomic<T> &other)
    {
      _a.store(other._a.load());
      return *this;
    }

    Atomic_wrapper &operator=(const Atomic_wrapper &other)
    {
      _a.store(other._a.load());
      return *this;
    }

    Atomic_wrapper &operator=(const T &other)
    {
      _a.store(other);
      return *this;
    }

    operator T() const
    {
      return _a.load();
    }

    operator std::atomic<T>() const
    {
      return _a;
    }
  };*/


  template <typename K>
  std::vector<typename K::Vector_d>
  compute_gram_schmidt_basis(
    std::vector<typename K::Vector_d> const& input_basis,
    K const& kernel)
  {
    typedef typename K::FT            FT;
    typedef typename K::Vector_d      Vector;
    typedef std::vector<Vector>       Basis;

    // Kernel functors
    typename K::Squared_length_d        sqlen      = kernel.squared_length_d_object();
    typename K::Scaled_vector_d         scaled_vec = kernel.scaled_vector_d_object();
    typename K::Scalar_product_d        inner_pdct = kernel.scalar_product_d_object();
    typename K::Difference_of_vectors_d diff_vec   = kernel.difference_of_vectors_d_object();

    Basis output_basis;

    typename Basis::const_iterator inb_it = input_basis.begin();
    typename Basis::const_iterator inb_it_end = input_basis.end();
    for (int i = 0 ; inb_it != inb_it_end ; ++inb_it, ++i)
    {
      Vector u = *inb_it;

      typename Basis::iterator outb_it = output_basis.begin();
      for (int j = 0 ; j < i ; ++j)
      {
        Vector const& ej = *outb_it;
        Vector u_proj = scaled_vec(ej, inner_pdct(u, ej));
        u = diff_vec(u, u_proj);
      }

      output_basis.push_back(
        scaled_vec(u, FT(1)/CGAL::sqrt(sqlen(u))));
    }

    return output_basis;
  }
  
  // Compute all the k-combinations of elements
  // Output_iterator::value_type must be std::set<std::size_t> >
  template <typename Elements_container, typename Output_iterator>
  void combinations(const Elements_container elements, int k,
                        Output_iterator combinations)
  {
    std::size_t n = elements.size();
    std::vector<bool> booleans(n, false);
    std::fill(booleans.begin() + n - k, booleans.end(), true);
    do
    {
      std::set<std::size_t> combination;
      Elements_container::const_iterator it_elt = elements.begin();
      for (std::size_t i = 0 ; i < n ; ++i, ++it_elt)
      {
        if (booleans[i])
          combination.insert(*it_elt);
      }
      *combinations++ = combination;

    } while (std::next_permutation(booleans.begin(), booleans.end()));
  }

} // namespace Tangential_complex_
} //namespace CGAL

#endif // CGAL_TC_UTILITIES_H
