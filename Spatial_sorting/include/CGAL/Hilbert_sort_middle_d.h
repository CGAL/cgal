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
//
// Author(s)     :  Olivier Devillers

#ifndef CGAL_HILBERT_SORT_MIDDLE_d_H
#define CGAL_HILBERT_SORT_MIDDLE_d_H

#include <CGAL/basic.h>
#include <functional>
#include <cstddef>
#include <CGAL/Hilbert_sort_middle_base.h>

namespace CGAL {

namespace internal {

    template <class K>
    struct Fixed_hilbert_cmp_d
        : public std::binary_function<typename K::Point_d,
                                      typename K::Point_d, bool>
    {
        typedef typename K::Point_d Point;
        K k;
	int axe;
	bool orient;
	double value;
        Fixed_hilbert_cmp_d (int a, bool o, double v, const K &_k = K()) 
	  :  k(_k), axe(a), orient(o), value(v) {}
        bool operator() (const Point &p) const
        {
	  return (orient  
		  ? ( to_double( k.compute_coordinate_d_object() (p,axe) ) >  value)
		  : ( to_double( k.compute_coordinate_d_object() (p,axe) ) <= value));
        }
    };

}

template <class K>
class Hilbert_sort_middle_d
{
public:
    typedef K Kernel;
    typedef typename Kernel::Point_d Point;
    typedef std::vector< bool > Starting_position;
    typedef std::vector< double > Corner;

private:
    Kernel _k;
    std::ptrdiff_t _limit;
    mutable int _dimension;
    mutable int two_to_dim;

    struct Cmp : public internal::Fixed_hilbert_cmp_d<Kernel>
    { Cmp (int a, bool dir, double v, const Kernel &k) 
	: internal::Fixed_hilbert_cmp_d<Kernel> (a,dir,v,k) {} };

public:
    Hilbert_sort_middle_d (const Kernel &k = Kernel(), std::ptrdiff_t limit = 1)
        : _k(k), _limit (limit)
    {}

    template <class RandomAccessIterator>
    void sort (RandomAccessIterator begin, RandomAccessIterator end,
	       Starting_position start, int direction, 
	       Corner mini, Corner maxi) const
   {
     if (end - begin <= _limit) return;

     Corner med(_dimension);
     for( int i=0; i<_dimension; ++i) med[i]=(mini[i]+maxi[i])/2;
     Corner cmin=mini,cmax=med;

     std::vector<RandomAccessIterator> places(two_to_dim +1);
     std::vector<int>                  dir   (two_to_dim +1);
     places[0]=begin;
     places[two_to_dim]=end;

     int last_dir = (direction + _dimension) % _dimension;
     int current_dir = direction;
     int current_level_step =two_to_dim;
     do{
       int half_step = current_level_step/2;
       int left=0;
       int middle = half_step;
       int right=current_level_step;
       bool orient = start[current_dir];
       do{
	 dir[middle]    = current_dir; 
	 places[middle] = internal::fixed_hilbert_split 
                             (places[left], places[right], 
			      Cmp (current_dir,orient,med[current_dir],_k));
	 left =right;
	 right+=current_level_step;
	 middle+=current_level_step;
	 orient = ! orient;
       }while( left< two_to_dim);
       current_level_step = half_step;
       current_dir = (current_dir +1) % _dimension;
     }while (current_dir != last_dir);

     /////////////start recursive calls
     last_dir = (direction + _dimension -1) % _dimension;
     // first step is special
     sort( places[0], places[1], start, last_dir,cmin,cmax);
     cmin[last_dir] = med[last_dir];
     cmax[last_dir] = maxi[last_dir];
     

     for(int i=1; i<two_to_dim-1; i +=2){
       //std::cout<<i<<";"<<start[0]<<start[1]<<start[2]<<start[3]<<"/"<<dir[i+1]<<std::endl;
       sort( places[i  ], places[i+1], start, dir[i+1],cmin,cmax);
       cmax[ dir[i+1] ] =  (cmin[ dir[i+1]]==mini[ dir[i+1]])
	                    ? maxi[ dir[i+1] ] : mini[ dir[i+1] ];
       cmin[ dir[i+1] ] =  med[ dir[i+1] ];

       sort( places[i+1], places[i+2], start, dir[i+1],cmin,cmax);
       cmin[ dir[i+1] ] =  cmax[ dir[i+1] ];
       cmax[ dir[i+1] ] =  med[ dir[i+1] ];
       cmax[ last_dir ] = (cmax[last_dir]==maxi[last_dir])
	                    ? mini[ last_dir ] : maxi[ last_dir ];
       start[dir[i+1]] = !  start[dir[i+1]];
       start[last_dir] = !  start[last_dir];
     }

     //last step is special
     sort( places[two_to_dim-1], places[two_to_dim], start, last_dir,cmin,cmax);
    }


    template <class RandomAccessIterator>
    void operator() (RandomAccessIterator begin, RandomAccessIterator end) const
    {
      _dimension = _k.point_dimension_d_object()(*begin);
      two_to_dim = 1;
      Starting_position start(_dimension);
      Corner mini(_dimension),maxi(_dimension);

      for (int i=0; i<_dimension; ++i) 
	mini[i]=maxi[i]=to_double( _k.compute_coordinate_d_object() (*begin,i) );
      for(RandomAccessIterator it=begin+1; it<end; ++it){
	for (int i=0; i<_dimension; ++i){
	  double d=  to_double( _k.compute_coordinate_d_object() (*it,i) );
	  if (d < mini[i]) mini[i] = d;
	  if (d > maxi[i]) maxi[i] = d;
	}
      }


      for (int i=0; i<_dimension; ++i) {
	start[i]=false; 	// we start below in all coordinates
	two_to_dim *= 2;        // compute 2^_dimension
	if (two_to_dim*2 <= 0) {
	  CGAL_assertion(end-begin < two_to_dim);//too many points in such dim
	  break;
	}
      }



      // we start with  direction 0;
      sort (begin, end, start, 0, mini, maxi);
    }
};




} // namespace CGAL

#endif//CGAL_HILBERT_SORT_MIDDLE_d_H
