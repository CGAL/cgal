// Copyright (c) 2010  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
// Author(s)     :  Olivier Devillers

#ifndef CGAL_HILBERT_SORT_MEDIAN_d_H
#define CGAL_HILBERT_SORT_MEDIAN_d_H

#include <CGAL/basic.h>
#include <functional>
#include <cstddef>
#include <CGAL/Hilbert_sort_base.h>

namespace CGAL {

namespace internal {

    template <class K>
    struct Hilbert_cmp_d
        : public std::binary_function<typename K::Point_d,
                                      typename K::Point_d, bool>
    {
        typedef typename K::Point_d Point;
        K k;
	int axe;
	bool orient;
        Hilbert_cmp_d (int a, bool o, const K &_k = K()) 
	  : axe(a), orient(o), k(_k) {}
        bool operator() (const Point &p, const Point &q) const
        {
	  return (orient  ? (k.compute_coordinate_d_object() (p,axe) 
			     > k.compute_coordinate_d_object() (q,axe))
		          : (k.compute_coordinate_d_object() (p,axe) 
			     < k.compute_coordinate_d_object() (q,axe)) ) ;
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
      _dimension = begin->dimension();
      two_to_dim = 1;
      Starting_position start(_dimension);


      for (int i=0; i<_dimension; ++i) {
	start[i]=false; 	// we start below in all coordinates
	two_to_dim *= 2;        // compute 2^_dimension
	if (two_to_dim*2 <= 0) {
	  CGAL_assertion(end-begin < two_to_dim);//too many points in such dim
	  break;
	}
      }

      // we start with  direction 0;
      sort (begin, end, start, 0);
    }
};




} // namespace CGAL

#endif//CGAL_HILBERT_SORT_MEDIAN_d_H
