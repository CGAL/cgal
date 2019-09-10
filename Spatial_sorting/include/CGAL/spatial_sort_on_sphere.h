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
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Olivier Devillers
//               : Pedro Machado Manhaes de Castro


#ifndef CGAL_SPATIAL_SORT_ON_SPHERE_H
#define CGAL_SPATIAL_SORT_ON_SPHERE_H

#include <CGAL/spatial_sort.h>
#include <CGAL/hilbert_sort_on_sphere.h>

namespace CGAL {

namespace internal {

template <class RandomAccessIterator, class Policy, class Kernel>
void spatial_sort_on_sphere (
                   RandomAccessIterator begin, RandomAccessIterator end,
                   const Kernel &k,
                   Policy /*policy*/,
                   typename Kernel::Point_3 *,
                   double sq_r,
                   const typename Kernel::Point_3 &p,
                   std::ptrdiff_t threshold_hilbert,
                   std::ptrdiff_t threshold_multiscale,
                   double ratio)
{
  typedef Hilbert_sort_on_sphere_3<Kernel, Policy> Sort;
    typedef std::iterator_traits<RandomAccessIterator> ITraits;
    typedef typename ITraits::difference_type Diff_t;

    boost::rand48 random;
    boost::random_number_generator<boost::rand48, Diff_t> rng(random);
    CGAL::cpp98::random_shuffle(begin,end, rng);

            if (threshold_hilbert==0) threshold_hilbert=4;
            if (threshold_multiscale==0) threshold_multiscale=16;
            if (ratio==0.0) ratio=0.25;

    (Multiscale_sort<Sort> (Sort (k, sq_r, p, threshold_hilbert),
                            threshold_multiscale, ratio)) (begin, end);
}

} // end of namespace internal

template <class RandomAccessIterator, class Policy, class Kernel>
void spatial_sort_on_sphere (RandomAccessIterator begin, RandomAccessIterator end,
                   const Kernel &k,
		   Policy policy,
		   double sq_r=1.0,
		   const typename Kernel::Point_3 &p = typename Kernel::Point_3(0,0,0),
		   std::ptrdiff_t threshold_hilbert=0,
		   std::ptrdiff_t threshold_multiscale=0,
		   double ratio=0.0)
{
    typedef std::iterator_traits<RandomAccessIterator> ITraits;
    typedef typename ITraits::value_type               value_type;

    internal::spatial_sort_on_sphere(begin, end, k, policy, static_cast<value_type *> (0),
			   sq_r,p, threshold_hilbert,threshold_multiscale,ratio);
}

template <class RandomAccessIterator>
void spatial_sort_on_sphere (RandomAccessIterator begin, RandomAccessIterator end,
		   Hilbert_sort_median_policy policy,
		   double sq_r=1.0,
		   const typename CGAL::Kernel_traits<typename std::iterator_traits<RandomAccessIterator>::value_type>::Kernel::Point_3 &p =
		     typename CGAL::Kernel_traits<typename std::iterator_traits<RandomAccessIterator>::value_type>::Kernel::Point_3(0,0,0),
		   std::ptrdiff_t threshold_hilbert=0,
		   std::ptrdiff_t threshold_multiscale=0,
		   double ratio=0.0)
{
    typedef std::iterator_traits<RandomAccessIterator> ITraits;
    typedef typename ITraits::value_type               value_type;
    typedef CGAL::Kernel_traits<value_type>            KTraits;
    typedef typename KTraits::Kernel                   Kernel;

    spatial_sort_on_sphere (begin, end, Kernel(), policy, sq_r, p,
		  threshold_hilbert,threshold_multiscale,ratio);
}
template <class RandomAccessIterator>
void spatial_sort_on_sphere (RandomAccessIterator begin, RandomAccessIterator end,
		   Hilbert_sort_middle_policy policy,
		   double sq_r=1.0,
		   const typename CGAL::Kernel_traits<typename std::iterator_traits<RandomAccessIterator>::value_type>::Kernel::Point_3 &p =
		     typename CGAL::Kernel_traits<typename std::iterator_traits<RandomAccessIterator>::value_type>::Kernel::Point_3(0,0,0),
		   std::ptrdiff_t threshold_hilbert=0,
		   std::ptrdiff_t threshold_multiscale=0,
		   double ratio=0.0)
{
    typedef std::iterator_traits<RandomAccessIterator> ITraits;
    typedef typename ITraits::value_type               value_type;
    typedef CGAL::Kernel_traits<value_type>            KTraits;
    typedef typename KTraits::Kernel                   Kernel;

    spatial_sort_on_sphere (begin, end, Kernel(), policy, sq_r,p,
		  threshold_hilbert,threshold_multiscale,ratio);
}


template <class RandomAccessIterator, class Kernel>
void spatial_sort_on_sphere (RandomAccessIterator begin, RandomAccessIterator end,
		   const Kernel &k,
		   double sq_r = 1.0,
           const typename Kernel::Point_3 &p = typename Kernel::Point_3(0,0,0),
           std::ptrdiff_t threshold_hilbert=0,
		   std::ptrdiff_t threshold_multiscale=0,
		   double ratio=0.0)
{
    spatial_sort_on_sphere (begin, end, k,
		  Hilbert_sort_median_policy(), sq_r, p,
		  threshold_hilbert,threshold_multiscale,ratio);
}

template <class RandomAccessIterator>
void spatial_sort_on_sphere (RandomAccessIterator begin, RandomAccessIterator end,
           double sq_r = 1.0,
		   const typename CGAL::Kernel_traits<typename std::iterator_traits<RandomAccessIterator>::value_type>::Kernel::Point_3 &p =
	       typename CGAL::Kernel_traits<typename std::iterator_traits<RandomAccessIterator>::value_type>::Kernel::Point_3(0,0,0),
		   std::ptrdiff_t threshold_hilbert=0,
		   std::ptrdiff_t threshold_multiscale=0,
		   double ratio=0.0)
{
    spatial_sort_on_sphere (begin, end,
		  Hilbert_sort_median_policy(), sq_r,p,
		  threshold_hilbert,threshold_multiscale,ratio);
}

} // end of namespace CGAL

#endif // CGAL_SPATIAL_SORT_ON_SPHERE_H
