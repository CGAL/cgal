// Copyright (c) 2011  INRIA Sophia-Antipolis (France).
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
// Author(s)     :  Olivier Devillers

#ifndef CGAL_HILBERT_SORT_MEDIAN_d_H
#define CGAL_HILBERT_SORT_MEDIAN_d_H

#include <CGAL/config.h>
#include <functional>
#include <cstddef>
#include <iterator>
#include <vector>
#include <CGAL/Hilbert_sort_base.h>

namespace CGAL {

namespace internal {

    template <class K>
    struct Hilbert_cmp_d
        : public CGAL::binary_function<typename K::Point_d,
                                      typename K::Point_d, bool>
    {
        typedef typename K::Point_d Point;
        K k;
	int axe;
	bool orient;
        Hilbert_cmp_d (int a, bool o, const K &_k = K()) 
	  : k(_k), axe(a),  orient(o) {}
        bool operator() (const Point &p, const Point &q) const
        {
	  return (orient  ? (k.less_coordinate_d_object() (q,p,axe) )
	        	  : (k.less_coordinate_d_object() (p,q,axe) ));
        }
    };

}

template <class K>
class Hilbert_sort_median_d
{
public:
    typedef K Kernel;
    typedef typename Kernel::Point_d Point;
    typedef std::vector< bool > Starting_position;

private:
    Kernel _k;
    std::ptrdiff_t _limit;
    mutable int _dimension;
    mutable int two_to_dim;

    struct Cmp : public internal::Hilbert_cmp_d<Kernel>
    { Cmp (int a, bool dir, const Kernel &k) 
	  : internal::Hilbert_cmp_d<Kernel> (a,dir,k) {} };

public:
    Hilbert_sort_median_d(const Kernel &k = Kernel(), std::ptrdiff_t limit = 1)
        : _k(k), _limit (limit)
    {}

    template <class RandomAccessIterator>
    void sort (RandomAccessIterator begin, RandomAccessIterator end,
	       Starting_position start, int direction) const
   {
     if (end - begin <= _limit) return;

     int nb_directions = _dimension;
     int nb_splits     = two_to_dim;

     if ( (end-begin) < (two_to_dim/2) ) { // not many points
       nb_splits = 1;
       nb_directions = 0;
       while ( (end-begin) > nb_splits) {
	 ++nb_directions;
	 nb_splits *= 2;        // compute 2^nb_directions
       }
     }

     std::vector<RandomAccessIterator> places(nb_splits +1);
     std::vector<int>                  dir   (nb_splits +1);
     places[0]=begin;
     places[nb_splits]=end;

     int last_dir = (direction + nb_directions) % _dimension;
     int current_dir = direction;
     int current_level_step =nb_splits;
     do{
       int half_step = current_level_step/2;
       int left=0;
       int middle = half_step;
       int right=current_level_step;
       bool orient = start[current_dir];
       do{
	 dir[middle]    = current_dir; 
	 places[middle] = internal::hilbert_split 
	   (places[left], places[right], Cmp (current_dir,orient,_k));
	 left =right;
	 right+=current_level_step;
	 middle+=current_level_step;
	 orient = ! orient;
       }while( left< nb_splits);
       current_level_step = half_step;
       current_dir = (current_dir +1) % _dimension;
     }while (current_dir != last_dir);

     if ( end-begin < two_to_dim) return; // less than 2^dim points

     /////////////start recursive calls
     last_dir = (direction + _dimension -1) % _dimension;
     // first step is special
     sort( places[0], places[1], start, last_dir);

     for(int i=1; i<two_to_dim-1; i +=2){
       sort( places[i  ], places[i+1], start, dir[i+1]);
       sort( places[i+1], places[i+2], start, dir[i+1]);
       start[dir[i+1]] = !  start[dir[i+1]];
       start[last_dir] = !  start[last_dir];
     }

     //last step is special
     sort( places[two_to_dim-1], places[two_to_dim], start, last_dir);
    }


    template <class RandomAccessIterator>
    void operator() (RandomAccessIterator begin, RandomAccessIterator end) const
    {
      _dimension = _k.point_dimension_d_object()(*begin);
      two_to_dim = 1;
      Starting_position start(_dimension);

      typename std::iterator_traits<RandomAccessIterator>::difference_type N=end-begin;
      N*=2;
      for (int i=0; i<_dimension; ++i) 	start[i]=false; 	// we start below in all coordinates
      for (int i=0; i<_dimension; ++i) {
	two_to_dim *= 2;        // compute 2^_dimension
	N/=2;
	if (N==0) break;  // not many points, this number of dimension is enough
      }


      // we start with  direction 0;
      sort (begin, end, start, 0);
    }
};




} // namespace CGAL

#endif//CGAL_HILBERT_SORT_MEDIAN_d_H
