// Copyright (c) 2015  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Olivier Devillers
//               : Pedro Machado Manhaes de Castro

#ifndef CGAL_HILBERT_SORT_ON_SPHERE_H
#define CGAL_HILBERT_SORT_ON_SPHERE_H

#include <CGAL/hilbert_sort.h>
#include <CGAL/Hilbert_sort_on_sphere_3.h>

namespace CGAL {

namespace internal {

template <class RandomAccessIterator, class Kernel, class Policy>
void hilbert_sort_on_sphere (RandomAccessIterator begin,
                   RandomAccessIterator end,
       const Kernel &k,
                   Policy,
                   typename Kernel::Point_3 *,
                   double sq_r,
                   const typename Kernel::Point_3 &p)
{
  boost::rand48 random;
  boost::random_number_generator<boost::rand48> rng(random);
  std::random_shuffle(begin,end, rng);
  (Hilbert_sort_on_sphere_3<Kernel, Policy> (k,sq_r,p))(begin, end);
}

} //end of namespace internal

template <class RandomAccessIterator>
void hilbert_sort_on_sphere (RandomAccessIterator begin, RandomAccessIterator end,
  double sq_r = 1.0,
  const typename CGAL::Kernel_traits<typename std::iterator_traits<RandomAccessIterator>::value_type>::Kernel::Point_3 &p =
	    typename CGAL::Kernel_traits<typename std::iterator_traits<RandomAccessIterator>::value_type>::Kernel::Point_3(0,0,0))
{

    typedef std::iterator_traits<RandomAccessIterator> ITraits;
    typedef typename ITraits::value_type               value_type;
    typedef CGAL::Kernel_traits<value_type>            KTraits;
    typedef typename KTraits::Kernel                   Kernel;

    internal::hilbert_sort_on_sphere(begin, end, Kernel(), Hilbert_sort_median_policy(),
				  static_cast<value_type *> (0), sq_r, p);

}

template <class RandomAccessIterator>
void hilbert_sort_on_sphere (RandomAccessIterator begin, RandomAccessIterator end, Hilbert_sort_median_policy policy,
  double sq_r = 1.0,
  const typename CGAL::Kernel_traits<typename std::iterator_traits<RandomAccessIterator>::value_type>::Kernel::Point_3 &p =
    typename CGAL::Kernel_traits<typename std::iterator_traits<RandomAccessIterator>::value_type>::Kernel::Point_3(0,0,0))
{
    typedef std::iterator_traits<RandomAccessIterator> ITraits;
    typedef typename ITraits::value_type               value_type;
    typedef CGAL::Kernel_traits<value_type>            KTraits;
    typedef typename KTraits::Kernel                   Kernel;

    internal::hilbert_sort_on_sphere(begin, end, Kernel(), policy,
				  static_cast<value_type *> (0), sq_r, p);
}


template <class RandomAccessIterator>
void hilbert_sort_on_sphere (RandomAccessIterator begin, RandomAccessIterator end, Hilbert_sort_middle_policy policy,
  double sq_r = 1.0,
  const typename CGAL::Kernel_traits<typename std::iterator_traits<RandomAccessIterator>::value_type>::Kernel::Point_3 &p =
    typename CGAL::Kernel_traits<typename std::iterator_traits<RandomAccessIterator>::value_type>::Kernel::Point_3(0,0,0))
{
    typedef std::iterator_traits<RandomAccessIterator> ITraits;
    typedef typename ITraits::value_type               value_type;
    typedef CGAL::Kernel_traits<value_type>            KTraits;
    typedef typename KTraits::Kernel                   Kernel;

    internal::hilbert_sort_on_sphere(begin, end, Kernel(), policy,
				  static_cast<value_type *> (0), sq_r, p);
}


template <class RandomAccessIterator, class Kernel, class Policy>
void hilbert_sort_on_sphere (RandomAccessIterator begin, RandomAccessIterator end, const Kernel &k, Policy policy,
	double sq_r = 1.0, const typename Kernel::Point_3 &p = typename Kernel::Point_3(0,0,0))
{
    typedef std::iterator_traits<RandomAccessIterator> ITraits;
    typedef typename ITraits::value_type               value_type;

    internal::hilbert_sort_on_sphere(begin, end,
			   k, policy, static_cast<value_type *> (0), sq_r, p);
}

} // end of namespace CGAL

#endif // CGAL_HILBERT_SORT_ON_SPHERE_H
