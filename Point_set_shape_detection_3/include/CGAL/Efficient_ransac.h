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

#include "Cylinder.h"
#include "Plane.h"
#include "Sphere.h"

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

extern int octTime;
extern int runTime;
extern int candGenCount;
extern int candGenTime;
extern int findBestCount;
extern int findBestTime;
extern int improveCount;
extern int improveTime;


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

namespace CGAL {

  namespace Efficient_ransac {
#ifdef DEBUG
#define printd(...) printf(__VA_ARGS__)
#else
#define printd(...) printf(__VA_ARGS__)
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
    class Efficient_ransac
    {
    public:
      typedef typename std::vector<inputDataType>::iterator inputIterator;
      //--------------------------------------------typedef
      typedef typename Kernel::FT FT;

      // basic geometric types
      typedef typename Kernel::Point_3 Point;
      typedef typename Kernel::Vector_3 Vector;
      typedef typename Kernel::Line_3 Line;
      typedef typename std::pair<Point, Vector> Pwn;

      //more advance container for geometric types
      typedef typename std::vector<Point>::const_iterator Point_const_iterator;

      typedef Primitive_ab<Kernel, inputDataType> Primitive;
      typedef Octree<Kernel, DirectPointAccessor<inputDataType>, inputDataType> DirectOctree;
      typedef Octree<Kernel, IndexedPointAccessor<inputDataType>, inputDataType> IndexedOctree;

      //--------------------------------------------typedef

    public:
      Efficient_ransac() {};
      ~Efficient_ransac() {/*if (m_global_tree) delete m_global_tree;*/}
      Efficient_ransac(inputIterator first, inputIterator beyond);

      void /*std::pair<std::vector<Primitive *>, std::vector<int> >*/ run(bool save = false);		 //the primitive and the index of the point sorted
      void addPrimitives(Primitive* p) {m_list_sought_primitives.insert(p->type());delete p;}
      const std::vector<Primitive *> &getExtractedPrimitives() {return m_extractedPrimitives;}
      int unassignedPoints() {return m_numAvailablePoints;}
      void setOptions(Option_RANSAC _m) {m_options = _m;}

    private:
      int getLevelOctree();
      Primitive* getBestCandidate(std::vector<Primitive* >& l_list_candidates, const int _SizeP);
      inline FT getRadiusSphere_from_level(int l_level_octree) {return m_max_radiusSphere_Octree/powf(2, l_level_octree);}
      bool improveBound(Primitive *candidate, const int _SizeP, unsigned int max_subset, unsigned int min_points);
      void rescoreCandidates(std::vector<Primitive *> &candidates);

      inline FT StopProbability(FT _sizeC, FT _np, FT _dC, FT _l) const
      {
        //printf("stop without thr %f\n", 	 std::pow(1.f - _sizeC/ (_np * _l * 4), _dC));
        return std::min<float>(std::pow(1.f - _sizeC / (_np * _l * 4), _dC), 1.);		//4 is (1 << (m_reqSamples - 1))) with m_reqSamples=3 (min number of points to create a candidate)
      }
      static bool candComp(const Primitive* a, const Primitive* b) {
        return a->ExpectedValue() < b->ExpectedValue();
      }

      //--------------------------------------------Functions


      //--------------------------------------------Variables
    public:
      static const unsigned int scm_max_depth_octree = 4;

    private:
      Option_RANSAC m_options;

      DirectOctree **m_direct_octrees;
      IndexedOctree *m_global_octree;
      std::vector<int> m_availableOctreeSizes;
      unsigned int m_num_subsets;
      std::vector<int> m_shapeIndex; // maps index into points to assigned extracted primitive
      unsigned int m_numAvailablePoints;
      std::vector<int> m_index_subsets;  //give the index of the subset of point i

      std::vector<Primitive *> m_extractedPrimitives;

      std::set<int> m_list_sought_primitives;
      inputIterator m_it_Point_Normal, m_last_it_Point_Normal; //copy iterator of input points

      FT m_max_radiusSphere_Octree;
      std::vector<FT> m_level_weighting;  	//sum must be 1
      //--------------------------------------------Variables
    };
    
    template <typename Kernel, class inputDataType>
    Efficient_ransac<Kernel, inputDataType>::Efficient_ransac(inputIterator first, inputIterator beyond)
    {
      srand ( time(NULL) );

      clock_t s = clock();

      m_global_octree = new IndexedOctree(first, beyond, 0);
      m_numAvailablePoints = beyond - first;
      
      m_it_Point_Normal = first;
      m_last_it_Point_Normal = beyond;

      //create subsets ------------------------------------------------
      //how many subsets ?
      m_num_subsets = std::max( (int) std::floor(std::logf(m_numAvailablePoints)/std::logf(2) )-9, 2);
      printd("number of subtrees: %d\n", m_num_subsets);

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
       octTime = clock() - s;
      //m_global_octree->verify();

      printd("init Ransac done\n");
    };	 

    template <typename Kernel, class inputDataType>
    int Efficient_ransac<Kernel, inputDataType>::getLevelOctree()
    {
      return abs(getRandomInt()) % (m_global_octree->maxLevel() + 1);
    };	

    template <typename Kernel, class inputDataType>
    Primitive_ab<Kernel, inputDataType>* Efficient_ransac<Kernel, inputDataType>::getBestCandidate(std::vector<Primitive* >& l_list_candidates, const int _SizeP) {
      if (l_list_candidates.size() == 1)	return l_list_candidates.back();
      clock_t s = clock();

      int index_worse_candidate = 0;
      bool improved = true;
      while (index_worse_candidate < l_list_candidates.size()-1 && improved)		//quit if find best candidate or no more improvement possible
      {
        improved = false;
        std::sort(l_list_candidates.begin()+index_worse_candidate, l_list_candidates.end(), [](Primitive* a, Primitive* b)
        {
          return a->maxBound() < b->maxBound();
        });


        //refine the best one 
        improveBound(l_list_candidates.back(), _SizeP, m_num_subsets, m_options.m_minNbPoints);

        int position_stop;

        //Take all those intersecting the best one
        for (position_stop= l_list_candidates.size()-2;position_stop>=index_worse_candidate;position_stop--)
        {
          if  (l_list_candidates.back()->minBound() > l_list_candidates.at(position_stop)->maxBound() ) break;//the intervals do not overlaps anymore

          if  (l_list_candidates.at(position_stop)->maxBound() <= m_options.m_minNbPoints) break;  //the following candidate doesnt have enough points !

          //if we reach this point, there is an overlaps between best one and position_stop
          //so request refining bound on position_stop
          improved |= improveBound(l_list_candidates.at(position_stop), _SizeP, m_num_subsets, m_options.m_minNbPoints);//this include the next subset for computing bounds, -> the score is then returned by ExpectedValue()	
          
          //test again after refined
          if  (l_list_candidates.back()->minBound() > l_list_candidates.at(position_stop)->maxBound() ) break;//the intervals do not overlaps anymore
        }

        index_worse_candidate = position_stop;
      }

      if (index_worse_candidate == l_list_candidates.size()-1 && improved) printf("delete everything, should never happen (%d) !!!!!!\n", index_worse_candidate);

      findBestTime += clock() - s;
      findBestCount++;

      return l_list_candidates.back();
    };
    
    template <typename Kernel, class inputDataType>
    bool Efficient_ransac<Kernel, inputDataType>::improveBound(Primitive *candidate, const int _SizeP, unsigned int max_subset, unsigned int min_points) {
      if (candidate->m_nb_subset_used >= max_subset)
        return false;
      if (candidate->m_nb_subset_used >= m_num_subsets)
        return false;

      improveCount++;
      clock_t s = clock();

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
      int l_new_sampled_points = 0;

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
        l_new_sampled_points += m_availableOctreeSizes[candidate->m_nb_subset_used];

        l_sum_score += l_new_score;

        candidate->m_nb_subset_used++;
      }
      while (l_new_sampled_points < min_points && candidate->m_nb_subset_used < m_num_subsets);

      candidate->m_score = candidate->m_indices.size();

//       if (l_new_score == 0) {
//         return false;
//       };

      candidate->computeBound(l_nb_total_points_subsets, _SizeP);//estimate the bound

      improveTime += clock() - s;

      return true;
    }

    template <typename Kernel, class inputDataType>
    void Efficient_ransac<Kernel, inputDataType>::run(bool save) {
      //no primitives added, exit
      if (m_list_sought_primitives.size() == 0) return;

      clock_t s = clock();

      // initializing the shape index
      m_shapeIndex.resize(m_numAvailablePoints, -1);

      std::vector<Primitive* > l_result;  // Y <- 0	   (no result)
      std::vector<Primitive* > l_list_candidates;  // C <- 0	   (no candidate)

      int requiredSamples = 4;

      int l_index_point_p1;
      int l_level_octree ;		//level 0 is full size, i.e 
      FT l_radius_Sphere;

      FT bestExp = 0;

      int numInvalid = 0; // number of points that have been assigned to a shape

      int nbNewCandidates = 0;
      int nbFailedCandidates = 0;
      bool forceExit = false;

      printf("Starting...\n");
      do	  //main loop
      {
        bestExp = 0;

        clock_t sti = clock();
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
          
          nbNewCandidates++;
          
          //add candidate for each type of primitives
          for(std::set<int>::iterator it =  m_list_sought_primitives.begin(); it != m_list_sought_primitives.end(); it++)	{
            Primitive *p = Primitive::create(*it, m_options.m_epsilon, m_options.m_normalThresh);
            p->compute(	indices, m_it_Point_Normal);	//compute the primitive and says if the candidate is valid 


            if (p->isValid()) {
              improveBound(p, m_numAvailablePoints - numInvalid, 1, 500);//this include the next subset for computing bounds, -> the score is then returned by ExpectedValue()

              //evaluate the candidate
              if(p->maxBound() >= m_options.m_minNbPoints) {
                if (bestExp < p->ExpectedValue()) bestExp = p->ExpectedValue();
                l_list_candidates.push_back(p);
              }
              else {
                nbFailedCandidates++;
                nbNewCandidates--;
                delete p;
              }        
            }
            else {
              nbFailedCandidates++;
              delete p;
            }
          }

          candGenCount++;

          if (nbFailedCandidates >= 10000)
            forceExit = true;

        } while( !forceExit
          && StopProbability(bestExp, m_numAvailablePoints - numInvalid, nbNewCandidates /*l_list_candidates.size()*/, scm_max_depth_octree)                 > m_options.m_probability 
          && StopProbability(m_options.m_minNbPoints, m_numAvailablePoints - numInvalid, nbNewCandidates /*l_list_candidates.size()*/, scm_max_depth_octree) > m_options.m_probability);
        // end of generate candidate
        candGenTime += clock() - sti;

        if (forceExit) break;

        if (l_list_candidates.empty())
          continue;

        //now get the best candidate in the current set of all candidates
        //Note that the function sort the candidates: the best candidate is always the last element of the vector
        Primitive *best_Candidate = getBestCandidate(l_list_candidates, m_numAvailablePoints - numInvalid);
        if (!best_Candidate)
          continue;

        best_Candidate->m_indices.clear();
        int s = m_global_octree->score(best_Candidate, m_it_Point_Normal, m_shapeIndex, 3 * m_options.m_epsilon, m_options.m_normalThresh);
        int cc = best_Candidate->connectedComponent(m_it_Point_Normal, m_options.m_bitmapEpsilon, m_global_octree->m_center, m_global_octree->m_width);

        //if the bestCandidate is good enough (proba of overlook lower than m_options.m_probability)
        if (StopProbability(best_Candidate->ExpectedValue(), (m_numAvailablePoints - numInvalid), 
          nbNewCandidates/*l_list_candidates.size()*/, scm_max_depth_octree) <= m_options.m_probability)
        {
          //we keep it
          if (best_Candidate->getPointsIndices()->size() >=  m_options.m_minNbPoints) 
          {
            l_list_candidates.back() = NULL;	//put null like that when we delete the vector, the object is not lost (pointer saved in bestCandidate)

            //1. add best candidate to final result.
            l_result.push_back(best_Candidate);

            //2. remove the points
            //2.1 update boolean
            std::vector<int> *indices_points_best_candidate = best_Candidate->getPointsIndices();

            for (int i=0;i<indices_points_best_candidate->size();i++)
            {
              if (m_shapeIndex[indices_points_best_candidate->at(i)] != -1) {
                std::cout << "shapeIndex already assigned!" << std::endl;
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
            }


            //2.3 block also the points for the subtrees        

            nbNewCandidates--;
            nbFailedCandidates = 0;
            bestExp = 0;
            //nbNewCandidates = l_result.size();

            //3. refit	(not implemented yet)
            //best_Candidate->LSfit();

            //4. save primitive
            if (save) {
              std::cout << "extracted primitive: " << best_Candidate->info().c_str() << std::endl;
              std::stringstream ss;
              ss << "ours_" << best_Candidate->type_str() << "_" << l_result.size() << ".ply";
              best_Candidate->save(ss.str().c_str(), m_it_Point_Normal);
            }
            m_extractedPrimitives.push_back(best_Candidate);
          }

          std::vector<int> subsetSizes(m_num_subsets);
          subsetSizes[0] = m_availableOctreeSizes[0];
          for (unsigned int i = 1;i<m_num_subsets;i++) {
            subsetSizes[i] = subsetSizes[i-1] + m_availableOctreeSizes[i];
          }


          //5. Remove points from candidates common with extracted primitive
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
        }

      }
      while( 	 
        StopProbability(m_options.m_minNbPoints, m_numAvailablePoints - numInvalid, nbNewCandidates /*l_list_candidates.size()*/, scm_max_depth_octree) > m_options.m_probability
        && FT(m_numAvailablePoints - numInvalid) >= m_options.m_minNbPoints
        );
      std::cout << "run Ransac done" << std::endl;

      m_numAvailablePoints -= numInvalid;

      runTime = clock() - s;
    };
  }
}


#endif
