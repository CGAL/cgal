// Copyright (c) 2013 INRIA Sophia-Antipolis (France).
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
// $URL$
// $Id$
//
//
// Author(s)     : Yannick Verdié, Clément Jamin
//
//******************************************************************************
// File Description :
//
//
//******************************************************************************

#ifndef CGAL_EFFICIENT_RANSAC_H
#define CGAL_EFFICIENT_RANSAC_H

#include "Types.h"
#include "Octree.h"
#include "Primitive.h"
#include "Plane.h"
#include "Sphere.h"
#include "Cylinder.h"

//for octree ------------------------------
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_iso_box.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/basic.h>
#include <CGAL/Search_traits_adapter.h>
#include <boost/iterator/zip_iterator.hpp>
#include <CGAL/bounding_box.h> //----------

#include <utility> // defines std::pair
#include <set>	    // defines std::set
#include <list>
#include <vector>
#include <fstream>
#include <sstream>
#include <random>
#define  _USE_MATH_DEFINES
#include <math.h>


//boost --------------
#include <boost/iterator/counting_iterator.hpp>
#include <boost/bind/make_adaptable.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <functional>
#include <boost/foreach.hpp>
#include <boost/tuple/tuple.hpp>
//---------------------

#define D18


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

namespace CGAL {

  namespace Efficient_ransac {
#ifdef DEBUG
#define printd(...) printf(__VA_ARGS__)
#else
#define printd(...) printf(__VA_ARGS__)
#endif

#ifdef D1
#define printd1(...) printf(__VA_ARGS__)
#else
#define printd1(...)
#endif
    
    //m for members
    //l for local
    //c for constants/readonlys
    //p for pointer (and pp for pointer to pointer)
    //v for volatile
    //s for static
    //i for indexes and iterators
    //e for eventstemplate <typename Kernel, class inputDataType>

    template <typename Kernel, class inputDataType>
    class valid_points
    {
      //--------------------------------------------typedef
      typedef typename Kernel::FT FT;

      // basic geometric types
      typedef typename Kernel::Point_3 Point;
      typedef typename Kernel::Point_2 Point_2d;
      typedef typename Kernel::Vector_3 Vector;
      typedef typename Kernel::Plane_3 Plane_3;
      typedef typename boost::tuple<Point, int> Point_and_int;

      typedef typename Kernel::Iso_cuboid_3 Iso_cuboid_3;
      typedef typename Kernel::Iso_rectangle_2 bbox_2d;

      //more advance container for geometric types
      typedef typename std::vector<Point>::const_iterator Point_const_iterator;
      typedef typename std::pair<Point, Vector> Pwn;
      typedef typename boost::tuple<Pwn, int> Pwn_and_int;
      typedef typename std::vector<Pwn> Pwn_vector;
      typedef typename Pwn_vector::iterator Pwn_iterator;
      typedef typename std::vector<Pwn_and_int>::iterator Pwn_and_int_iterator;

      //for kdTree
      typedef typename CGAL::Search_traits_3<Kernel> Traits_base;
      typedef typename CGAL::Search_traits_adapter<Point_and_int, My_point_property_map<Kernel>, Traits_base>    Traits;
      typedef typename CGAL::Kd_tree<Traits> Tree;
      typedef typename CGAL::Fuzzy_sphere<Traits> Fuzzy_sphere;

      typedef Primitive_ab<Kernel, inputDataType> Primitive;

      typedef boost::function<Point(Pwn&)> getPointFunc;
      typedef boost::transform_iterator<getPointFunc, Pwn_iterator> point_iterator;

      typedef boost::function<Vector(Pwn&)> getNormalFunc;
      typedef boost::transform_iterator<getNormalFunc, Pwn_iterator> normal_iterator;

    private:
      point_iterator m_point;
      normal_iterator m_normal;
      const FT m_epsilon;
      const FT m_ndev;
      const Plane_3 m_plane;

    public:	 
      //template <typename T> T operator ( ) (T& _elem ){return _elem;} 

      // Constructor initializes the value to multiply by

      valid_points  ( const Plane_3 _p, const FT _gam, const FT _al, const point_iterator _point, const normal_iterator _Val) : m_point(_point), m_normal ( _Val ), m_epsilon(_gam), m_plane(_p), m_ndev(_al){}

      int operator ( ) (int& _elem ) 
      {
        if (CGAL::squared_distance ( m_plane, *(m_point+_elem)) > m_epsilon*m_epsilon) return -1;

        normal_iterator it  = 	(m_normal +_elem);
        //if ( abs( *(it) * m_plane.orthogonal_vector () ) <  m_ndev*sqrtf(it->squared_length()*m_plane.orthogonal_vector ().squared_length()))  return -1;
        if ( (*(it) * m_plane.orthogonal_vector ())*( *(it) * m_plane.orthogonal_vector () ) <  m_ndev*m_ndev*(it->squared_length()*m_plane.orthogonal_vector ().squared_length()))  return -1;


        return _elem;
      };	  	
    };

    template <typename Kernel, class inputDataType>
    class Efficient_ransac
    {
    public:
      typedef typename std::vector<inputDataType>::iterator inputIterator;
      //--------------------------------------------typedef
      typedef typename Kernel::FT FT;

      // basic geometric types
      typedef typename Kernel::Point_3 Point;
      typedef typename Kernel::Vector_3 Vector;
      typedef typename std::pair<Point, Vector> Pwn;
      typedef typename boost::tuple<Point, int> Point_and_int;

      typedef typename Kernel::Iso_cuboid_3 Iso_cuboid_3;

      //more advance container for geometric types
      typedef typename std::vector<Point>::const_iterator Point_const_iterator;
      typedef typename boost::tuple<inputDataType, int> Pwn_and_int;

      typedef typename std::vector<Pwn_and_int>::iterator Pwn_and_int_iterator;

      //for kdTree
      typedef typename CGAL::Search_traits_3<Kernel> Traits_base;
      typedef typename CGAL::Search_traits_adapter<Point_and_int, My_point_property_map<Kernel>, Traits_base>    Traits;
      typedef typename CGAL::Kd_tree<Traits> Tree;
      typedef typename CGAL::Fuzzy_sphere<Traits> Fuzzy_sphere;

      typedef Primitive_ab<Kernel, inputDataType> Primitive;
      typedef Plane<Kernel, inputDataType> DebugPlane;
      typedef Octree<Kernel, DirectPointAccessor<inputDataType>, inputDataType> DirectOctree;
      typedef Octree<Kernel, IndexedPointAccessor<inputDataType>, inputDataType> IndexedOctree;

      //boost
      typedef boost::function<Point(inputDataType&)> getPointFunc;
      typedef boost::transform_iterator<getPointFunc, inputIterator> point_iterator;
      point_iterator m_it_Point, m_last_it_Point;

      typedef boost::function<Vector(inputDataType&)> getNormalFunc;
      typedef boost::transform_iterator<getNormalFunc, inputIterator> normal_iterator;
      normal_iterator m_it_Normal;
      //--------------------------------------------typedef




      //--------------------------------------------Functions
    public:
      Efficient_ransac() {};
      ~Efficient_ransac() {/*if (m_global_tree) delete m_global_tree;*/}
      Efficient_ransac(inputIterator first, inputIterator beyond);

      void /*std::pair<std::vector<Primitive *>, std::vector<int> >*/ run();		 //the primitive and the index of the point sorted
      void addPrimitives(Primitive* p) {m_list_sought_primitives.insert(p->type());delete p;}
      void setOptions(Option_RANSAC _m) {m_options = _m;}

    private:
      //Kernel:: Iso_cuboid_3 getCell(int l_index_point_p1);
      int getLevelOctree();
      Primitive* getBestCandidate(std::vector<Primitive* >& l_list_candidates, const int _SizeP, const point_iterator _it_Point, const normal_iterator _it_Normal);
      inline FT getRadiusSphere_from_level(int l_level_octree) {return m_max_radiusSphere_Octree/powf(2, l_level_octree);}
      bool improveBound(Primitive *candidate, const int _SizeP, const point_iterator _it_Points, const normal_iterator _it_Normal, unsigned int max_subset, unsigned int min_points);
      void rescoreCandidate(Primitive *candidate);

      inline FT StopProbability(FT _sizeC, FT _np, FT _dC, FT _l) const
      {
        //printf("stop without thr %f\n", 	 std::pow(1.f - _sizeC/ (_np * _l * 4), _dC));
        return std::min(std::pow(1.f - _sizeC / (_np * _l * 4), _dC), 1.);		//4 is (1 << (m_reqSamples - 1))) with m_reqSamples=3 (min number of points to create a candidate)
      }

      //--------------------------------------------Functions


      //--------------------------------------------Variables
    public:
      static const unsigned int scm_max_depth_octree = 4;

    private:
      Option_RANSAC m_options;

      //Tree* m_global_tree;
      //std::vector<Tree*> m_subsets_tree;
      DirectOctree **m_direct_octrees;
      IndexedOctree *m_global_octree;
      std::vector<int> m_availableOctreeSizes;
      unsigned int m_num_subsets;
      std::vector<int> m_shapeIndex; // maps index into points to assigned extracted primitive
      unsigned int m_numAvailablePoints;
      //std::vector<CGAL::Efficient_ransac::Octree<Kernel, DirectPointAccessor<Pri>> m_subsets_octree;
      //std::pair<std::vector<Tree*>*, std::vector<int> > m_subsets_tree_int;
      std::vector<int> m_index_subsets;  //give the index of the subset of point i

      //std::vector<std::vector<std::vector<voxel> > > voxelize_space;
      //std::vector<voxel> voxelize_space;
      //FT m_size_voxel; int m_width; int m_height;  int m_depth;
      //Iso_cuboid_3 m_bbox;
      std::set<int> m_list_sought_primitives;
      //std::vector<int> m_index_available_points;
      inputIterator m_it_Point_Normal, m_last_it_Point_Normal; //copy iterator of input points

      FT m_max_radiusSphere_Octree;
      std::vector<FT> m_level_weighting;  	//sum must be 1
      //--------------------------------------------Variables
    };

    //int current(0);
    //int increment () { return current++; }


    //requirement: object inputDataType has x(), y(), and z() implemented

    template <typename Kernel, class inputDataType>
    Efficient_ransac<Kernel, inputDataType>::Efficient_ransac(inputIterator first, inputIterator beyond)
    {
      srand ( time(NULL) );

      m_global_octree = new IndexedOctree(first, beyond, 0);
      m_numAvailablePoints = beyond - first;
      
      m_it_Point_Normal = first;
      m_last_it_Point_Normal = beyond;

      //Index of available points
      //sequence 0, 1 , ..., nbPoints - 1
/*
      int count = 0;
      for(inputIterator it = first; it != beyond;it++)
        m_index_available_points.push_back(count++);*/

      //std::iota(m_index_available_points.begin(), m_index_available_points.end(), 0);


      //init octree -------------------------
      getPointFunc f = boost::bind( &Pwn::first, _1);
      m_it_Point = boost::make_transform_iterator(first, f);
      m_last_it_Point  = boost::make_transform_iterator(beyond, f);

      getNormalFunc f2 = boost::bind( &Pwn::second, _1);
      m_it_Normal = boost::make_transform_iterator(first, f2);

      //create subsets ------------------------------------------------
      //how many subsets ?
      m_num_subsets = std::max( (int) std::floor(std::logf(m_numAvailablePoints)/std::logf(2) )-9, 2);
      printd("number of subtree: %d\n", m_num_subsets);

      // SUBSET GENERATION ->
      // approach with increasing subset sizes -> replace with octree later on
      std::random_device rd;
      std::mt19937 gen(rd());
      std::uniform_int_distribution<> dis(0, 1);
      
      inputIterator end = beyond - 1;
      int remainingPoints = m_numAvailablePoints;

      m_availableOctreeSizes.resize(m_num_subsets);
      m_direct_octrees = new DirectOctree *[m_num_subsets];
      std::cout << "subSetSizes: ";
      for (int s = m_num_subsets - 1;s>=0;--s) {
        int subsetSize = remainingPoints;
        std::vector<unsigned int> indices(subsetSize);
        if (s) {
          subsetSize >>= 1;
          for (unsigned int i = 0;i<subsetSize;i++) {
            int index = dis(gen);
            index = index + (i<<1);
            index = (index >= remainingPoints) ? remainingPoints - 1 : index;
            indices[i] = index;
          }

          // move points to the end of the point vector
          for (int i = subsetSize - 1;i>=0;i--) {
            inputDataType tmp = (*end);
            *end = first[indices[i]];
            first[indices[i]] = tmp;
            end--;
          }
          m_direct_octrees[s] = new DirectOctree(end + 1, end + subsetSize + 1, remainingPoints - subsetSize);
        }
        else
          m_direct_octrees[0] = new DirectOctree(first, first + (subsetSize), 0);

        std::cout << subsetSize << " ";

        m_availableOctreeSizes[s] = subsetSize;
        m_direct_octrees[s]->createTree();

        remainingPoints -= subsetSize;
      }

      m_global_octree = new IndexedOctree(first, beyond);
      m_global_octree->createTree();
      m_global_octree->verify();

      // <- approach with increasing subset sizes

      // YANNICKs version ->
      //create a random vector of number in [0, l_nb_subsets-1] thus as this indicates to which subset the point belongs to.
      /*srand(time(NULL));
      std::random_device rd2;
      std::mt19937 gen2(rd2());
      std::uniform_int_distribution<> dis2(0, m_num_subsets-1);
      std::generate_n(back_inserter(m_index_subsets), m_index_available_points.size(), [&]{ return  dis2(gen2); });

      //create the tree for each subset
      std::vector<std::vector<Point_and_int> > l_vector_subsets(m_num_subsets);
      count = 0;
      for(point_iterator it =  m_it_Point; it !=  m_last_it_Point;it++, count++)
        l_vector_subsets[m_index_subsets[count]].push_back(boost::make_tuple(*it, count));
      
      // SUBSET GENERATION <-

      std::vector<int> l_size_subsets;
      for (int i=0;i<m_num_subsets;i++)
      {
        m_subsets_tree.push_back(new Tree(l_vector_subsets[i].begin(), l_vector_subsets[i].end()) );
        l_size_subsets.push_back(m_subsets_tree[i]->size());
        std::cout << " " << m_subsets_tree[i]->size();
      }
      std::cout << std::endl;

      //group subset with number of available points in the subset (need for the scoring)
      m_subsets_tree_int	= std::make_pair(&m_subsets_tree, l_size_subsets);
      //---------------------------------------------------------------


      //bounding box and parameters for octree  
      m_bbox = bounding_box (m_it_Point, m_last_it_Point);

      // widen by very small delta
      m_bbox = Iso_cuboid_3(m_bbox.xmin() - 0.000001, m_bbox.ymin() - 0.000001, m_bbox.zmin() - 0.000001, m_bbox.xmax() + 0.000001, m_bbox.ymax() + 0.000001, m_bbox.zmax() + 0.000001);

      int l_max_sizeCell = ceil(max (max ( m_bbox.xmax() - m_bbox.xmin() , m_bbox.ymax() - m_bbox.ymin()), m_bbox.zmax() - m_bbox.zmin() ));			
      m_max_radiusSphere_Octree = l_max_sizeCell * sqrtf(3)/2 ;
      m_level_weighting = std::vector<FT>(scm_max_depth_octree, 1.f/scm_max_depth_octree); 
      //-------------------------------------

      //init voxelize for attaching point to primitive at the end.
      m_size_voxel = sqrtf(2)*m_options.m_epsilon;		//round down

      m_width = int(ceil( (m_bbox.xmax () - m_bbox.xmin ())/m_size_voxel)) ; 
      m_height= int(ceil( (m_bbox.ymax () - m_bbox.ymin ())/m_size_voxel)) ;
      m_depth = int(ceil( (m_bbox.zmax () - m_bbox.zmin ())/m_size_voxel)) ;

      //voxelize_space =  std::vector< voxel >(m_width*m_height*m_depth);
      voxelize_space = std::vector< std::vector< std::vector< voxel > > >(m_width , 
        std::vector< std::vector< voxel > >  (m_height , 
        std::vector< voxel >				   (m_depth )));

      //voxelize_space.reserve(m_index_available_points.size()); 
      count = 0;
      //now put indices of input data in the voxel
      for(point_iterator it = m_it_Point; it != m_last_it_Point;it++)
      { 
        //voxelize_space[int((it->x() - m_bbox.xmin ())/m_size_voxel) + int((it->y() - m_bbox.ymin ())/m_size_voxel)*m_width + int((it->z() - m_bbox.zmin ())/m_size_voxel)*m_width*m_height ].push_back(count++); 
        voxelize_space[(it->x() - m_bbox.xmin ())/m_size_voxel][(it->y() - m_bbox.ymin ())/m_size_voxel][(it->z() - m_bbox.zmin ())/m_size_voxel].push_back(count++); 
      }*/

      printd("init Ransac done\n");
    };	 

    template <typename Kernel, class inputDataType>
    int Efficient_ransac<Kernel, inputDataType>::getLevelOctree()
    {
/*
      int l_level_octree = 0;
      FT selected = FT(rand())/RAND_MAX;	 //between [0, 1]

      FT l_sum_CDF = m_level_weighting[0];
      while(selected > l_sum_CDF) {
        l_sum_CDF += m_level_weighting[l_level_octree++];
        if (l_level_octree == m_level_weighting.size() -1)
          break;
      };*/

      return abs(getRandomInt()) % (m_global_octree->maxLevel() + 1);
    };	

    template <typename Kernel, class inputDataType>
    Primitive_ab<Kernel, inputDataType>* Efficient_ransac<Kernel, inputDataType>::getBestCandidate(std::vector<Primitive* >& l_list_candidates, const int _SizeP, const point_iterator _it_Point, const normal_iterator _it_Normal)
    {
      //printf("in getBestCandidate\n");
      if (l_list_candidates.size() == 1)	return l_list_candidates.back();

      int index_worse_candidate = 0;
      bool improved = true;
      while (index_worse_candidate < l_list_candidates.size()-1 && improved)		//quit if find best candidate or no more improvement possible
      {
        improved = false;
        std::sort(l_list_candidates.begin()+index_worse_candidate, l_list_candidates.end(), [](Primitive* a, Primitive* b)
        {
          return a->maxBound() < b->maxBound();
          //return a->ExpectedValue() < b->ExpectedValue();
        });

        /*
        for (int i=0;i<  l_list_candidates.size();i++)
        l_list_candidates.at(i)->info();

        system("pause");
        */

        //printf("best min/max = %f %f | min/max second %f %f\n", 
        //l_list_candidates.at(l_list_candidates.size()-1)->minBound(), 	l_list_candidates.at(l_list_candidates.size()-1)->maxBound(), 
        //l_list_candidates.at(l_list_candidates.size()-2)->minBound(), 	l_list_candidates.at(l_list_candidates.size()-2)->maxBound() );
        //l_list_candidates.back()->info();  l_list_candidates.at(l_list_candidates.size()-2)->info();

        //refine the best one 
        improveBound(l_list_candidates.back(), _SizeP, _it_Point, _it_Normal, m_num_subsets, m_options.m_minNbPoints);

        int position_stop;
        //printf("-- sorted here %d, exp %f %f\n", l_list_candidates.size(), l_list_candidates.back()->minBound(), l_list_candidates.at( l_list_candidates.size()-1)->ExpectedValue());


        //printf("best min/max = %f %f | min/max second %f %f\n", 
        //l_list_candidates.at(l_list_candidates.size()-1)->minBound(), 	l_list_candidates.at(l_list_candidates.size()-1)->maxBound(), 
        //l_list_candidates.at(l_list_candidates.size()-2)->minBound(), 	l_list_candidates.at(l_list_candidates.size()-2)->maxBound() );
        //l_list_candidates.back()->info();  l_list_candidates.at(l_list_candidates.size()-2)->info();
        //system("pause");
        //Take all those intersecting the best one
        for (position_stop= l_list_candidates.size()-2;position_stop>=index_worse_candidate;position_stop--)
        {
          if  (l_list_candidates.back()->minBound() > l_list_candidates.at(position_stop)->maxBound() ) break;//the intervals do not overlaps anymore

          if  (l_list_candidates.at(position_stop)->maxBound() <= m_options.m_minNbPoints) break;  //the following candidate doesnt have enough points !

          //printf("%d intersect (over %d)", position_stop, l_list_candidates.size());
          //if we reach this point, there is an overlaps between best one and position_stop
          //so request refining bound on position_stop
          improved |= improveBound(l_list_candidates.at(position_stop), _SizeP, _it_Point, _it_Normal, m_num_subsets, m_options.m_minNbPoints);//this include the next subset for computing bounds, -> the score is then returned by ExpectedValue()	

          //printf("%d improved %d\n", 	position_stop, 	improved);
          //if (!improved) break; //all the subsets have been used, but no improvement, quit !

          //test again after refined
          if  (l_list_candidates.back()->minBound() > l_list_candidates.at(position_stop)->maxBound() ) break;//the intervals do not overlaps anymore
        }


        //remove the condidates not processed so far
        //l_list_candidates.erase(l_list_candidates.begin(), l_list_candidates.begin()+position_stop+1);
        index_worse_candidate = position_stop;
      }

      if (index_worse_candidate == l_list_candidates.size()-1 && improved) printf("delete everything, should never happen (%d) !!!!!!\n", index_worse_candidate);
      //printf("best with value %f\n", l_list_candidates.at(0)->ExpectedValue());
      //l_list_candidates.back()->info();
      // system("pause");
      //printf("out getBestCandidate\n");



      return l_list_candidates.back();
    };
    
    template <typename Kernel, class inputDataType>
    bool Efficient_ransac<Kernel, inputDataType>::improveBound(Primitive *candidate, const int _SizeP, const point_iterator _it_Points, const normal_iterator _it_Normal, unsigned int max_subset, unsigned int min_points) {
      candidate->m_nb_subset_used = (candidate->m_nb_subset_used >= m_num_subsets) ? m_num_subsets - 1 : candidate->m_nb_subset_used;
      
      //what it does is add another subset and recompute lower and upper bound
      //the next subset to include is provided by m_nb_subset_used

      //1. sum nb points of previous subset, 	  and score
      int l_sum_score = candidate->m_score;
      int l_nb_total_points_subsets = 0;
      for (int i=0;i<candidate->m_nb_subset_used;i++)
        l_nb_total_points_subsets += m_availableOctreeSizes[i];	   //no some points are forbidden now

      //2. need score of new subset as well as sum of the score of the previous considered subset
      int l_new_score = 0;

      do
      {
        //l_new_score = candidate->score(_t.first, candidate->m_query_ball, _it_Points, _it_Normal, is_available);
        l_new_score = m_direct_octrees[candidate->m_nb_subset_used]->score(candidate, m_it_Point_Normal, m_shapeIndex, m_options.m_epsilon, m_options.m_normalThresh);
        //candidate->m_Indices.clear();
        //l_new_score3 = m_direct_octrees[candidate->m_nb_subset_used]->fullScore(candidate, m_shapeIndex, m_options.m_epsilon, m_options.m_normalThresh);
        candidate->m_score += l_new_score;

        //int tmp = score(octree)
        //printf("subset %d exp value %f\n", m_nb_subset_used, l_new_score);
        //candidate->m_score_by_subset.push_back(l_new_score);	//add to vector of score

        l_nb_total_points_subsets += m_availableOctreeSizes[candidate->m_nb_subset_used]; //add point new subset
        l_sum_score += l_new_score;

        candidate->m_nb_subset_used++;
      }
      while (l_sum_score < min_points && candidate->m_nb_subset_used < m_num_subsets);
      // while (l_new_score == 0 && m_nb_subset_used < _t.first->size());YANNICK	   

      if (l_new_score == 0) {
        //std::cout << "c" << std::endl;
        //printd(info().c_str());
        return false;
      };

      //printf("->>>>>  %d %d\n" , l_nb_total_points_subsets, _SizeP);
      candidate->computeBound(l_nb_total_points_subsets, _SizeP);//estimate the bound

      //m_nb_subset_used++;
      //std::cout << info() << " d" << std::endl;
      return true;
    }
/*

    template <typename Kernel, class inputDataType>
    bool Efficient_ransac<Kernel, inputDataType>::rescoreCandidate(Primitive *candidate) {
      candidate->m_nb_subset_used = (candidate->m_nb_subset_used >= m_num_subsets) ? m_num_subsets - 1 : candidate->m_nb_subset_used;

      //what it does is add another subset and recompute lower and upper bound
      //the next subset to include is provided by m_nb_subset_used

      //1. sum nb points of previous subset, 	  and score
      int l_sum_score = 0;
      int l_nb_total_points_subsets = 0;
      candidate->m_score_by_subset.clear();
      for (int i=0;i<candidate->m_nb_subset_used;i++)
      {
        l_nb_total_points_subsets += m_direct_octrees[i]->size();	   //no some points are forbidden now
        l_sum_score += candidate->m_score_by_subset[i] ;
      }

      //2. need score of new subset as well as sum of the score of the previous considered subset
      int l_new_score = 0;

      do
      {
        //l_new_score = candidate->score(_t.first, candidate->m_query_ball, _it_Points, _it_Normal, is_available);
        l_new_score = m_direct_octrees[candidate->m_nb_subset_used]->score(candidate, m_shapeIndex, m_options.m_epsilon, m_options.m_normalThresh);
        //candidate->m_Indices.clear();
        //l_new_score3 = m_direct_octrees[candidate->m_nb_subset_used]->fullScore(candidate, m_shapeIndex, m_options.m_epsilon, m_options.m_normalThresh);
        candidate->m_score_by_subset.push_back(l_new_score);

        //int tmp = score(octree)
        //printf("subset %d exp value %f\n", m_nb_subset_used, l_new_score);
        //candidate->m_score_by_subset.push_back(l_new_score);	//add to vector of score

        l_nb_total_points_subsets += m_direct_octrees[candidate->m_nb_subset_used]->size(); //add point new subset
        l_sum_score += l_new_score;

        candidate->m_nb_subset_used++;
      }
      while (l_sum_score < min_points && candidate->m_nb_subset_used < m_num_subsets);
      // while (l_new_score == 0 && m_nb_subset_used < _t.first->size());YANNICK	   

      if (l_new_score == 0) {
        //std::cout << "c" << std::endl;
        //printd(info().c_str());
        return false;
      };

      //printf("->>>>>  %d %d\n" , l_nb_total_points_subsets, _SizeP);
      candidate->computeBound(l_nb_total_points_subsets, _SizeP, l_sum_score);//estimate the bound

      //m_nb_subset_used++;
      //std::cout << info() << " d" << std::endl;
      return true;
    }
*/

    template <typename Kernel, class inputDataType>
    void Efficient_ransac<Kernel, inputDataType>::run()
    {
      //no primitives added, exit
      if (m_list_sought_primitives.size() == 0) return;

      // initializing the shape index
      m_shapeIndex.resize(m_numAvailablePoints, -1);

      std::vector<Primitive* > l_result;  // Y <- 0	   (no result)
      std::vector<Primitive* > l_list_candidates;  // C <- 0	   (no candidate)

      int requiredSamples = 4;

/*
      for (unsigned int i = 0;i<m_list_sought_primitives.size();i++) {
        requiredSamples = std::max<int>(requiredSamples, m_list_sought_primitives[i])
      }*/

      int l_index_point_p1;
      int l_level_octree ;		//level 0 is full size, i.e 
      FT l_radius_Sphere;

      FT bestExp = 0;

      int numInvalid = 0; // number of points that have been assigned to a shape

      //somehow, we also need a vector<bool> to tell if the points is available
      //std::vector<bool> l_points_available(m_index_available_points.size(), true);
      //m_shapeIndex.resize(m_numAvailablePoints, -1);

      int nbNewCandidates = 0;
      int nbFailedCandidates = 0;  bool forceExit = false;

      printf("Starting...\n");
      do	  //main loop
      {
        bestExp = 0;

        //nbNewCandidates = 0;
        //TODO: do this in batch of n candidate, and thus can be parallelized in n threads
        do	  //generate candidate
        {
          //1. pick a point p1 randomly among available points
          std::set<int> indices;
          bool done = false;

          do {
            l_index_point_p1 = getRandomInt() % m_numAvailablePoints;
            //done = m_global_octree->getPointsInCellContainingPoint((m_it_Point_Normal + l_index_point_p1)->first, getLevelOctree(), indices, m_shapeIndex);
            done = m_global_octree->drawSamplesFromCellContainingPoint((m_it_Point_Normal + l_index_point_p1)->first, getLevelOctree(), indices, m_shapeIndex, requiredSamples);
          } while (m_shapeIndex[l_index_point_p1] != -1 || !done);
          //l_index_point_p1 = m_index_available_points[numInvalid + getRandomInt() % (m_index_available_points.size()-numInvalid)];

          // indices are already filtered by shapeIndex
/*
          std::set<int> l_list_index_selected;
          l_list_index_selected.insert(l_index_point_p1);
          do {
            l_list_index_selected.insert(getRandomInt() % indices.size());
          } while (l_list_index_selected.size() < requiredSamples);*/

          // select number of necessary points from indices

          //printf("îck %d\n", 	l_index_point_p1);

          //2. pick a cell size randomly among cell that include the point p1
          //in our implementation, we use sphere centered on p1, thus no need to check that the cell includes p1

/*          l_level_octree = getLevelOctree();		//level 0 is full size, i.e 

          l_radius_Sphere = getRadiusSphere_from_level(l_level_octree);
          Fuzzy_sphere s_query(*(m_it_Point + l_index_point_p1), l_radius_Sphere, 0.1);
          std::vector<Point_and_int> l_points;

          //for(int i = 0;i< m_num_subsets;i++)
            //m_subsets_tree[i]->search(std::back_inserter(l_points), s_query);

          //3. Pick 2 more points in the l_points and give it to each primitives to create candidate
          int l_nb_try = 0;
          l_list_index_selected.insert(l_index_point_p1 );

          while (l_list_index_selected.size() < 3 &&  l_nb_try++ < l_points.size())
          {
            int id_p2 = boost::get<1>(l_points[getRandomInt() % l_points.size()]);	//pick a point in the sphere randomly
            if (m_shapeIndex[id_p2] == -1)												//is it free ?
              l_list_index_selected.insert(id_p2);									//put in a set (avoid duplicate with p1)

            //printf("id picked : %d %d\n", id_p2, int(l_points_available[id_p2]));
          }
          if (l_list_index_selected.size() < 3 &&  l_nb_try < l_points.size()) {printd("error !, impossible to pick 3 random points\n"); return;};*/
          
          //add candidate for each type of primitives
          for(std::set<int>::iterator it =  m_list_sought_primitives.begin(); it != m_list_sought_primitives.end(); it++)		  //1 type
          {
            Primitive *p = Primitive::create(*it, m_options.m_epsilon, m_options.m_normalThresh);
            p->compute(	indices, m_it_Point_Normal);	//compute the primitive and says if the candidate is valid 

            /*
            Vector n = ((DebugPlane *)p)->m_plane.orthogonal_vector();
            n = n * (1.0 / sqrt(n.squared_length()));
            auto it2 = l_list_index_selected.begin();

            if (abs(n.x()) > 0.9999 || abs(n.y()) > 0.9999 || abs(n.z()) > 0.9999) {
            printd("%s, %i %i %i", p->info().c_str(), *it2++, *it2++, *it2++);
            }*/

            if (p->isValid()) {
              //printf("in valid\n");  
              //NOTE: octree NEED to discard points already assigned !
              //improveBound(p, m_subsets_tree_int, (m_index_available_points.size() - numInvalid), m_it_Point, m_it_Normal, l_points_available, 1, 500);//this include the next subset for computing bounds, -> the score is then returned by ExpectedValue()
              improveBound(p, m_numAvailablePoints - numInvalid, m_it_Point, m_it_Normal, 1, 500);//this include the next subset for computing bounds, -> the score is then returned by ExpectedValue()

              // debug out ->
              /*
              Vector n = ((DebugPlane *)p)->m_plane.orthogonal_vector();
              if (abs(n.z()) > 0.9) {
              n = n * (1.0 / sqrt(n.squared_length()));
              FT d = (CGAL::ORIGIN - ((DebugPlane *)p)->m_point_on_primitive) * n;
              if (abs(d) > 0.5);
              }*/
              // <- debug out
              //printd(p->info().c_str());

              //printf("(%f %f    ----- %f)\n", p->minBound(), p->maxBound(), p->getLastScore());

              //update level sampling
              //paper does that way
              //(*sampleLevelScores)[node->Level()].first += cand.ExpectedValue();
              //++(*sampleLevelScores)[node->Level()].second;

              //evaluate the candidate
              if(p->maxBound() >= m_options.m_minNbPoints) {
                //printf("add new candidate %f %f\n", p->minBound(), p->maxBound());
                if (bestExp < p->ExpectedValue()) bestExp = p->ExpectedValue();
                l_list_candidates.push_back(p);
                nbNewCandidates++;
              }
              else {
                nbFailedCandidates++;
                delete p;
              }        
            }
            else {
              nbFailedCandidates++;
              delete p;
            }
          }

          if (nbFailedCandidates >= 100000)
            forceExit = true;
          //printf("%d %d\n", nbNewCandidates, nbFailedCandidates);

          // printf("%d %f ->>> %f %f\n", nbNewCandidates, bestExp, 
          //	    StopProbability(bestExp, (m_index_available_points.size()-numInvalid), nbNewCandidates, scm_max_depth_octree), 
          //		StopProbability(m_options.m_minNbPoints, (m_index_available_points.size()-numInvalid), nbNewCandidates, scm_max_depth_octree));

        } while( !forceExit
          && StopProbability(bestExp, m_numAvailablePoints - numInvalid, nbNewCandidates /*l_list_candidates.size()*/, scm_max_depth_octree)                 > m_options.m_probability 
          && StopProbability(m_options.m_minNbPoints, m_numAvailablePoints - numInvalid, nbNewCandidates /*l_list_candidates.size()*/, scm_max_depth_octree) > m_options.m_probability);
        // end of generate candidate

        printd1("#candidates: %i\n", l_list_candidates.size());

        if (forceExit) break;
        // printf("criteria stop:%f\n", 
        //StopProbability(bestExp, (m_index_available_points.size()-numInvalid), nbNewCandidates /*l_list_candidates.size()*/, scm_max_depth_octree)                 , 
        //StopProbability(m_options.m_minNbPoints, m_index_available_points.size() - numInvalid, nbNewCandidates /*l_list_candidates.size()*/, scm_max_depth_octree) );

        // printf("size list candidates %d (best%f)\n", l_list_candidates.size(), bestExp);
        // for (int i=0;i<  l_list_candidates.size();i++)
        //	 l_list_candidates[i]->info();

        //now get the best candidate in the current set of all candidates
        //Note that the function sort the candidates: the best candidate is always the last element of the vector
        Primitive *best_Candidate = getBestCandidate(l_list_candidates, m_numAvailablePoints - numInvalid, m_it_Point, m_it_Normal);

        //if the bestCandidate is good enough (proba of overlook lower than m_options.m_probability)
        if (StopProbability(best_Candidate->ExpectedValue(), (m_numAvailablePoints - numInvalid), 
          nbNewCandidates/*l_list_candidates.size()*/, scm_max_depth_octree) <= m_options.m_probability)
        {
          //First, put best Candidate in first in the vector
          //printf("find best candidate\n");
          //1 find the points attached to the primitives (return immediately if points already attached)

          // Yannick is doing the connected component via the voxel grid
          best_Candidate->m_indices.clear();
          m_global_octree->score(best_Candidate, m_it_Point_Normal, m_shapeIndex, 3 * m_options.m_epsilon, m_options.m_normalThresh);
          best_Candidate->connectedComponent(m_it_Point_Normal, m_options.m_bitmapEpsilon, m_global_octree->m_center, m_global_octree->m_width);
          //best_Candidate->attachPoints(m_it_Point, m_it_Normal, voxelize_space, m_bbox, m_size_voxel, l_points_available);

          printd1("best candidate: ");
          printd1(best_Candidate->info().c_str());

          //we keep it
          if (best_Candidate->getPointsIndices()->size() >=  m_options.m_minNbPoints) 
          {
            l_list_candidates.back() = NULL;	//put null like that when we delete the vector, the object is not lost (pointer saved in bestCandidate)

            //1. add best candidate to final result.
            l_result.push_back(best_Candidate);

            //2. remove the points
            //2.1 update boolean
            std::vector<int> *indices_points_best_candidate = best_Candidate->getPointsIndices();
            //#pragma omp parallel for
            for (int i=0;i<indices_points_best_candidate->size();i++)
            {
              if (m_shapeIndex[indices_points_best_candidate->at(i)] != -1) {
                std::cout << "shapeIndex doppelt zugewiesen!" << std::endl;
              }
              m_shapeIndex[indices_points_best_candidate->at(i)] = l_result.size() - 1;
              numInvalid++;

              bool exactlyOnce = true;

              for (unsigned int j = 0;j<m_num_subsets;j++) {
                if (m_direct_octrees[j] && m_direct_octrees[j]->m_root) {
                  int offset = m_direct_octrees[j]->offset();
                  if (offset <= indices_points_best_candidate->at(i) && (indices_points_best_candidate->at(i) - offset) < m_direct_octrees[j]->size()) {
                    if (!exactlyOnce) {
                      std::cout << "adjusting available octree sizes failed! (twice)" << std::endl;
                    }
                    exactlyOnce = false;
                    m_availableOctreeSizes[j]--;
                  }
                }
              }

              if (exactlyOnce) {
                std::cout << "adjusting available octree sizes failed! (never)" << std::endl;
              }
              //l_points_available[	indices_points_best_candidate->at(i) ] = false;
              //m_subsets_tree_int.second[m_index_subsets[indices_points_best_candidate->at(i)]]--;	//the subset indicated by m_index_subsets has one point less
            }

            //2.2 swap index from m_index_available_points
            //std::vector<int >::iterator bound = std::partition (m_index_available_points.begin() + numInvalid, m_index_available_points.end(), [&l_points_available](int const &p){return !l_points_available[p];});
            //printf("nbInvalid %d %d\n", numInvalid , numInvalid +	std::distance(m_index_available_points.begin() + numInvalid, bound));
            //numInvalid += std::distance(m_index_available_points.begin() + numInvalid, bound);

            //2.3 block also the points for the subtrees        

            nbNewCandidates--;
            nbFailedCandidates = 0;
            bestExp = 0;
            //nbNewCandidates = l_result.size();

            //3. refit	(not implemented yet)
            //best_Candidate->LSfit();

            //4. save primitive
            std::cout << "extracted primitive: " << best_Candidate->info().c_str() << std::endl;
            std::stringstream ss;
            ss << best_Candidate->type_str() << "_" << l_result.size() << ".off";
            best_Candidate->save(ss.str().c_str(), m_it_Point_Normal);
            //system("pause");
          }

          std::vector<int> subsetSizes(m_num_subsets);
          subsetSizes[0] = m_availableOctreeSizes[0];
          for (unsigned int i = 1;i<m_num_subsets;i++) {
            subsetSizes[i] = subsetSizes[i-1] + m_availableOctreeSizes[i];
          }


          //5. delete ALL the candidate and empty the list
          //INSTEAD DELETE COMMON POINTS
#pragma omp parallel for
          for (int i=0;i< l_list_candidates.size()-1;i++)   //-1 because we keep BestCandidate
          {
            if (l_list_candidates[i]) {
              l_list_candidates[i]->updatePoints(m_shapeIndex);
              if (l_list_candidates[i]->score() < m_options.m_minNbPoints) {
                delete l_list_candidates[i];
                l_list_candidates[i] = NULL;
              }
              else l_list_candidates[i]->computeBound(subsetSizes[l_list_candidates[i]->m_nb_subset_used - 1], m_numAvailablePoints - numInvalid);
            }
          }

          int start = 0, end = l_list_candidates.size() - 1;
          while (start < end) {
            while (l_list_candidates[start] && start < end) start++;
            while (!l_list_candidates[end] && start < end) end--;
            if (!l_list_candidates[start] && l_list_candidates[end] && start < end) {
              l_list_candidates[start] = l_list_candidates[end];
              l_list_candidates[end] = NULL;
              start++;
              end--;
            }
          }

          l_list_candidates.resize(end);

          //l_list_candidates = std::vector<Primitive*>(0);

          /*
          // << << << << << << << here what I understood for the paper BUT it is slow, so for now try nomething else	>>>>>>>>>>>>>

          best_Candidate->save("plane.off", m_it_Point_Normal);
          system("pause");

          //compute for the other candidate or delete it
          printf("nb candidate %d\n", 	 l_list_candidates.size());
          //#pragma omp parallel for
          for (int i=0;i< l_list_candidates.size()-1;i++)
          {
          l_list_candidates[i]->attachPoints(m_it_Point, m_it_Normal, voxelize_space, m_bbox, m_size_voxel);
          }
          //-------------------------------------------------------------------------------------------


          //2. remove the points
          //instead of removing, we can move the point in front, count the nb of point attached to a candidate (numInvalid), and disregard those points during future iteration
          //2.1 First we update l_points_available and will use it for 2.2
          std::vector<int> *indices_points_best_candidate = best_Candidate->getPointsIndices();
          #pragma omp parallel for
          for (int i=0;i<indices_points_best_candidate->size();i++)
          {
          l_points_available[	indices_points_best_candidate->at(i) ] = false;
          }

          //2.2 even better, we are going to swap index from m_index_available_points
          //we process only the part that has not been partitioned before
          //here, all indices (from m_index_available_points) that are false in l_points_available are on the left partition
          std::vector<int >::iterator bound = std::partition (m_index_available_points.begin() + numInvalid, m_index_available_points.end(), [&l_points_available](int const &p){return !l_points_available[p];});
          numInvalid += std::distance(m_index_available_points.begin() + numInvalid, bound);

          //3. remove the candidates that have at least one point in common with bestCandidate
          #pragma omp parallel for
          for (int i=0;i< l_list_candidates.size()-1;i++)
          {
          bool deleteMe = l_list_candidates[i]->hasCommonPoints(indices_points_best_candidate);
          if (deleteMe)
          {
          delete l_list_candidates[i];
          l_list_candidates[i] = NULL;
          }
          }

          //and resize the vector of candidate afterward
          std::vector<Primitive* >::iterator bound2 = std::partition (l_list_candidates.begin(), l_list_candidates.end(), [](Primitive const *c){return !c==NULL;});
          l_list_candidates.resize(std::distance(l_list_candidates.begin(), bound2));


          //4. refit	(not implemented yet)
          //best_Candidate->LSfit();


          best_Candidate->save("plane.off", m_it_Point_Normal);
          system("pause");
          */

        }

      }
      while( 	 
        StopProbability(m_options.m_minNbPoints, m_numAvailablePoints - numInvalid, nbNewCandidates /*l_list_candidates.size()*/, scm_max_depth_octree) > m_options.m_probability
        && FT(m_numAvailablePoints - numInvalid) >= m_options.m_minNbPoints
        );


      //printf("criteria stop:%f %f\n", 
      //	StopProbability(m_options.m_minNbPoints, m_index_available_points.size() - numInvalid, nbNewCandidates /*l_list_candidates.size()*/, scm_max_depth_octree), 
      //	  FT(m_index_available_points.size() - numInvalid) );


      std::cout << "run Ransac done" << std::endl;

      //return std::pair<std::vector<Primitive *>, std::vector<int> > (std::vector<Primitive *>(), std::vector<int>());
    };
  }
}


#endif
