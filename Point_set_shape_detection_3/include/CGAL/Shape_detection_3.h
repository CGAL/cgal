// Copyright (c) 2013-2015 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Sven Oesau, Yannick Verdié, Clément Jamin, Pierre Alliez
//
//******************************************************************************
// File Description :
//
//
//******************************************************************************

#ifndef CGAL_SHAPE_DETECTION_3_H
#define CGAL_SHAPE_DETECTION_3_H

#include "Octree.h"
#include "Shape_base.h"

//for octree ------------------------------
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_iso_box.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/basic.h>
#include <CGAL/Search_traits_adapter.h>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/iterator/filter_iterator.hpp>
#include <CGAL/bounding_box.h> //----------

#include <utility> // defines std::pair
#include <set>	    // defines std::set
#include <list>
#include <vector>
#include <fstream>
#include <sstream>
#include <random>
#define  _USE_MATH_DEFINES
#include <cmath>


//boost --------------
#include <boost/iterator/counting_iterator.hpp>
#include <boost/bind/make_adaptable.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <functional>
#include <boost/foreach.hpp>
#include <boost/tuple/tuple.hpp>
//---------------------


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

/*! 
  \file Shape_detection_3.h
*/

namespace CGAL {
    /*!
     \brief Traits class for definition of types.
     */

  template <class Gt,
            class Input_it,
            class P_pmap,
            class N_pmap>
  struct Shape_detection_traits_3 {
    typedef Gt Geom_traits;     ///< Geometric types for definition of point, vector types, etc.
    typedef Input_it Input_iterator;  ///< Random access iterator type used for providing input data to the method.
    typedef P_pmap Point_pmap;    ///< Property map to access point location from input data.
    typedef N_pmap Normal_pmap;   ///< Property map to access normal vector from input data.
  };

  /*!
\brief Implementation of a RANSAC method for shape detection.

Given a point set in 3D space with unoriented normals, sampled on surfaces,
the method detects sets of connected points on the surface of primitive shapes.
Each input point is assigned to either none or at most one detected primitive shape.
This implementation follows the algorithm published by Schnabel
et al. in 2007 \cgalCite{Schnabel07}.

\tparam Shape detection traits class. 

*/
  template <class Sd_traits>
  class Shape_detection_3 {
  public:

    /// \cond SKIP_IN_MANUAL
    struct Filter_unassigned_points {
      Filter_unassigned_points() : m_shapeIndex() {}
      Filter_unassigned_points(const std::vector<int> &shapeIndex)
        : m_shapeIndex(shapeIndex) {}

      bool operator()(int x) {
        if (x < m_shapeIndex.size())
          return m_shapeIndex[x] == -1;
        else return true; // to prevent infinite incrementing
      }
      std::vector<int> m_shapeIndex;
    };
    /// \endcond

    /// \name Types 
    /// @{
    typedef typename Sd_traits::Input_iterator Input_iterator; ///< random access iterator for input data.
    typedef typename Sd_traits::Geom_traits::FT FT; ///< number type.
    typedef typename Sd_traits::Geom_traits::Point_3 Point; ///< point type.
    typedef typename Sd_traits::Geom_traits::Vector_3 Vector; ///< vector type.
    typedef typename Sd_traits::Point_pmap Point_pmap; ///< property map to access the location of an input point.
    typedef typename Sd_traits::Normal_pmap Normal_pmap; ///< property map to access the unoriented normal of an input point.

    typedef Shape_base<Sd_traits> Shape; ///< shape type.
    typedef typename std::vector<Shape *>::const_iterator Shape_iterator; ///< iterator for extracted shapes.
    typedef boost::filter_iterator<Filter_unassigned_points, boost::counting_iterator<std::size_t> > Point_index_iterator; ///< iterator for indices of points.
    /// @}
    
    /// \name Parameters 
    /// @{
      /*!
       \brief parameters for Shape_detection_3.
       */
    struct Parameters {
      Parameters() : probability(0.01), min_points(SIZE_MAX), epsilon(-1), normal_threshold(0.9), cluster_epsilon(-1) {}
      FT probability;         ///< Probability to control search endurance. Default value 0.05.
      size_t min_points;      ///< Minimum number of points of a shape. Default value 1% of total number of input points.
      FT epsilon;             ///< Maximum tolerance Euclidian distance from a point and a shape. Default value 1% of bounding box diagonal.
      FT normal_threshold;	  ///< Maximum tolerance normal deviation from a point's normal to the normal on shape at projected point. Default value 0.9 (around 25 degrees).
      FT cluster_epsilon;	    ///< Maximum distance between points to be considered connected. Default value 1% of bounding box diagonal.
    };
    /// @}

  private:
    typedef internal::Octree<internal::DirectPointAccessor<Sd_traits> >
      Direct_octree;
    typedef internal::Octree<internal::IndexedPointAccessor<Sd_traits> >
      Indexed_octree;
    //--------------------------------------------typedef

  public:

  /// \name Initialization
  /// @{


  /*! 
    Constructor to provide random access iterators over the input data and
    property maps to access point locations and unoriented normals.
    Internal data structures depending on the input data are constructed.
  */ 
    Shape_detection_3(Input_iterator first, ///< iterator over the first input point.
      Input_iterator beyond, ///< past-the-end iterator over the input points.
      Point_pmap point_pmap, ///< property map to access the position of an input point.
      Normal_pmap normal_pmap ///< property map to access the unoriented normal of an input point.
      ) : m_rng(std::random_device()()) {

      m_point_pmap = point_pmap;
      m_normal_pmap = normal_pmap;

      m_inputIterator_first = first;
      m_inputIterator_beyond = beyond;
      
      m_global_octree = new Indexed_octree(first, beyond, 0);
      m_numAvailablePoints = beyond - first;

      //create subsets ------------------------------------------------
      //how many subsets ?
      m_num_subsets = (std::max<size_t>)((size_t)
        std::floor(std::log(double(m_numAvailablePoints))/std::log(2.))-9, 2);

      // SUBSET GENERATION ->
      // approach with increasing subset sizes -> replace with octree later on
      Input_iterator last = beyond - 1;
      size_t remainingPoints = m_numAvailablePoints;

      m_availableOctreeSizes.resize(m_num_subsets);
      m_direct_octrees = new Direct_octree *[m_num_subsets];
      std::cout << "subSetSizes: ";
      for (int s = m_num_subsets - 1;s >= 0;--s) {
        size_t subsetSize = remainingPoints;
        std::vector<size_t> indices(subsetSize);
        if (s) {
          subsetSize >>= 1;
          for (size_t i = 0;i<subsetSize;i++) {
            size_t index = m_rng() % 2;
            index = index + (i<<1);
            index = (index >= remainingPoints) ? remainingPoints - 1 : index;
            indices[i] = index;
          }

          // move points to the end of the point vector
          for (int i = subsetSize - 1;i >= 0;i--) {
            typename std::iterator_traits<Input_iterator>::value_type
              tmp = (*last);
            *last = first[indices[size_t(i)]];
            first[indices[size_t(i)]] = tmp;
            last--;
          }
          m_direct_octrees[s] = new Direct_octree(last + 1,
                                                  last + subsetSize + 1,
                                                  remainingPoints - subsetSize);
        }
        else
          m_direct_octrees[0] = new Direct_octree(first,
                                                  first + (subsetSize), 
                                                  0);

        std::cout << subsetSize << " ";

        m_availableOctreeSizes[s] = subsetSize;
        m_direct_octrees[s]->createTree();

        remainingPoints -= subsetSize;
      }

      m_global_octree = new Indexed_octree(first, beyond);
      m_global_octree->createTree();
    }	 

    /// \cond SKIP_IN_MANUAL
    ~Shape_detection_3() {
      if (m_global_octree)
        delete m_global_octree;

      if (m_direct_octrees) {
        for (size_t i = 0;i<m_num_subsets;i++) {
          delete m_direct_octrees[i];
        }

        delete [] m_direct_octrees;
      }

      if (m_extractedShapes.size()) {
        for (size_t i = 0;i<m_extractedShapes.size();i++)
          delete m_extractedShapes[i];
      }

      for (size_t i = 0;i<m_shapeFactories.size();i++)
        delete m_shapeFactories[i];
    }
    /// \endcond

  /*!
    This function registers a shape type for detection.
  */ 
    template <class S>
    void add_shape_factory(Shape_factory<S> *factory///< Factory for shape type defined by 'Shape_factory<"shape_type">'
      ) {
      m_shapeFactories.push_back(factory);
    }

    /// @}

    /// \name Detection 
    /// @{
    /*! 
      This function initiates the shape detection. Shape types to be detected
      must be registered before with 'add_shape_factory'.
    */ 
    void detect(
      const Parameters &options = Parameters()///< Parameters for shape detection.
                ) {
      // no shape types for detection, exit
      if (m_shapeFactories.size() == 0) return;

      // use bounding box diagonal as reference for default values
      Bbox_3 bbox = m_global_octree->boundingBox();
      FT bboxDiagonal = sqrt((bbox.xmax() - bbox.xmin()) * (bbox.xmax() - bbox.xmin()) + (bbox.ymax() - bbox.ymin()) * (bbox.ymax() - bbox.ymin()) + (bbox.zmax() - bbox.zmin()) * (bbox.zmax() - bbox.zmin()));

      m_options = options;

      // epsilon or cluster_epsilon have been set by the user? if not, derive from bounding box diagonal
      m_options.epsilon = (m_options.epsilon < 0) ? bboxDiagonal * 0.01 : m_options.epsilon;
      m_options.cluster_epsilon = (m_options.cluster_epsilon < 0) ? bboxDiagonal * 0.01 : m_options.cluster_epsilon;

      // minimum number of points has been set?
      m_options.min_points = (m_options.min_points >= m_numAvailablePoints) ? (size_t)((FT)0.001 * m_numAvailablePoints) : m_options.min_points;
      
      // initializing the shape index
      m_shapeIndex.resize(m_numAvailablePoints, -1);

      // list of all randomly drawn candidates
      // with the minimum number of points
      std::vector<Shape *> candidates;

      int requiredSamples = 4;

      int firstSample; // first sample for RANSAC

      FT bestExp = 0;

      // number of points that have been assigned to a shape
      int numInvalid = 0;

      int nbNewCandidates = 0;
      int nbFailedCandidates = 0;
      bool forceExit = false;

      printf("Starting...\n");
      do { // main loop
        bestExp = 0;

        do {  //generate candidate
          //1. pick a point p1 randomly among available points
          std::set<size_t> indices;
          bool done = false;
          do {
            do 
              firstSample = m_rng() % m_numAvailablePoints;
              while (m_shapeIndex[firstSample] != -1);
              
            done = m_global_octree->drawSamplesFromCellContainingPoint(
              get(m_point_pmap, 
                  *(m_inputIterator_first + firstSample)),
              selectRandomOctreeLevel(),
              indices,
              m_shapeIndex,
              requiredSamples);

          } while (m_shapeIndex[firstSample] != -1 || !done);

          nbNewCandidates++;

          //add candidate for each type of primitives
          for(std::vector<internal::Shape_factory_base *>::iterator it =
            m_shapeFactories.begin(); it != m_shapeFactories.end(); it++)	{
            Shape *p = (Shape *) (*it)->create();
            //compute the primitive and says if the candidate is valid
            p->compute(indices,
                       m_inputIterator_first,
                       m_point_pmap,
                       m_normal_pmap,
                       m_options.epsilon, 
                       m_options.normal_threshold);

            if (p->is_valid()) {
              improveBound(p, m_numAvailablePoints - numInvalid, 1, 500);

              //evaluate the candidate
              if(p->max_bound() >= m_options.min_points) {
                if (bestExp < p->expected_value())
                  bestExp = p->expected_value();

                candidates.push_back(p);
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
          
          if (nbFailedCandidates >= 10000)
            forceExit = true;

        } while( !forceExit
          && StopProbability(bestExp,
                             m_numAvailablePoints - numInvalid, 
                             nbNewCandidates,
                             m_global_octree->maxLevel()) 
                > m_options.probability
          && StopProbability(m_options.min_points,
                             m_numAvailablePoints - numInvalid, 
                             nbNewCandidates,
                             m_global_octree->maxLevel())
                > m_options.probability);
        // end of generate candidate

        if (forceExit) {
          std::cout << "force exit" << std::endl;
          break;
        }

        if (candidates.empty())
          continue;

        // Now get the best candidate in the current set of all candidates
        // Note that the function sorts the candidates:
        //  the best candidate is always the last element of the vector

        Shape *best_Candidate = 
          getBestCandidate(candidates, m_numAvailablePoints - numInvalid);

        if (!best_Candidate)
          continue;

        best_Candidate->m_indices.clear();

        best_Candidate->m_score = 
          m_global_octree->score(best_Candidate,
                                 m_shapeIndex,
                                 3 * m_options.epsilon,
                                 m_options.normal_threshold);

        best_Candidate->connected_component(m_options.cluster_epsilon/*, 
                                            m_global_octree->m_center,
                                            m_global_octree->m_width*/);

        if (StopProbability(best_Candidate->expected_value(),
                            (m_numAvailablePoints - numInvalid),
                            nbNewCandidates,
                            m_global_octree->maxLevel())
                <= m_options.probability) {
          //we keep it
          if (best_Candidate->assigned_points().size() >=
                       m_options.min_points) {
            candidates.back() = NULL;

            //1. add best candidate to final result.
            m_extractedShapes.push_back(best_Candidate);

            //2. remove the points
            //2.1 update boolean
            const std::vector<size_t> &indices_points_best_candidate =
              best_Candidate->assigned_points();

            for (size_t i = 0;i<indices_points_best_candidate.size();i++) {
              if (m_shapeIndex[indices_points_best_candidate.at(i)] != -1) {
                std::cout << "shapeIndex already assigned!" << std::endl;
              }
              m_shapeIndex[indices_points_best_candidate.at(i)] =
                m_extractedShapes.size() - 1;

              numInvalid++;

              bool exactlyOnce = true;

              for (size_t j = 0;j<m_num_subsets;j++) {
                if (m_direct_octrees[j] && m_direct_octrees[j]->m_root) {
                  size_t offset = m_direct_octrees[j]->offset();

                  if (offset <= indices_points_best_candidate.at(i) &&
                      (indices_points_best_candidate.at(i) - offset) 
                      < m_direct_octrees[j]->size()) {
                    if (!exactlyOnce) {
                      std::cout << "adjusting available octree\
                                   sizes failed! (twice)" << std::endl;
                    }
                    exactlyOnce = false;
                    m_availableOctreeSizes[j]--;
                  }
                }
              }

              if (exactlyOnce) {
                std::cout << std::endl << "adjusting available octree\
                                          sizes failed! (never)" << std::endl;
              }
            }

            //2.3 block also the points for the subtrees        

            nbNewCandidates--;
            nbFailedCandidates = 0;
            bestExp = 0;
          }

          std::vector<size_t> subsetSizes(m_num_subsets);
          subsetSizes[0] = m_availableOctreeSizes[0];
          for (size_t i = 1;i<m_num_subsets;i++) {
            subsetSizes[i] = subsetSizes[i-1] + m_availableOctreeSizes[i];
          }


          //3. Remove points from candidates common with extracted primitive
          //#pragma omp parallel for
          for (size_t i=0;i< candidates.size()-1;i++)
          {
            if (candidates[i]) {
              candidates[i]->update_points(m_shapeIndex);
              if (candidates[i]->score() < m_options.min_points) {
                delete candidates[i];
                candidates[i] = NULL;
              }
              else 
                candidates[i]->compute_bound(
                   subsetSizes[candidates[i]->m_nb_subset_used - 1],
                   m_numAvailablePoints - numInvalid);
            }
          }

          size_t start = 0, end = candidates.size() - 1;
          while (start < end) {
            while (candidates[start] && start < end) start++;
            while (!candidates[end] && start < end) end--;
            if (!candidates[start] && candidates[end] && start < end) {
              candidates[start] = candidates[end];
              candidates[end] = NULL;
              start++;
              end--;
            }
          }

          candidates.resize(end);
        }

      }
      while(StopProbability(m_options.min_points,
                            m_numAvailablePoints - numInvalid,
                            nbNewCandidates, 
                            m_global_octree->maxLevel())
               > m_options.probability
        && FT(m_numAvailablePoints - numInvalid) >= m_options.min_points);

      m_numAvailablePoints -= numInvalid;
    }

    /// @}

    /// \name Access
    /// @{
      /*!
       Number of detected shapes.
       */
    std::size_t number_of_shapes() const {
      return m_extractedShapes.size();
    }
      
    /*!
      Iterator to the first detected shape.
      The order of the shapes is the order of the detection.
      Depending on the chosen probability for the detection
      the shapes are ordered with decreasing size.
    */

    Shape_iterator shapes_begin() const {
      return m_extractedShapes.begin();
    }
      
    /*!
      Past-the-end shape iterator.
    */
    Shape_iterator shapes_end() const {
      return m_extractedShapes.end();
    }
      
    /*! 
      Number of points not assigned to a shape.
    */ 
    std::size_t number_of_unassigned_points() {
      return m_numAvailablePoints;
    }
    
    /*! 
      Provides iterator to the index of the first point in
      the input data that has not been assigned to a shape.
    */ 
    Point_index_iterator unassigned_points_begin() {
      Filter_unassigned_points fup(m_shapeIndex);

      return boost::make_filter_iterator<Filter_unassigned_points>(
        fup,
        boost::counting_iterator<int>(0),
        boost::counting_iterator<int>(m_shapeIndex.size()));
    }
       
    /*! 
      Provides the past-the-end iterator
      for the indices of the unassigned points.
    */ 
    Point_index_iterator unassigned_points_end() {
      Filter_unassigned_points fup(m_shapeIndex);

      return boost::make_filter_iterator<Filter_unassigned_points>(
        fup,
        boost::counting_iterator<int>(0),
        boost::counting_iterator<int>(m_shapeIndex.size())).end();
    }
    /// @}

  private:
    int selectRandomOctreeLevel() {
      return m_rng() % (m_global_octree->maxLevel() + 1);
    }

    Shape* getBestCandidate(std::vector<Shape* >& candidates,
                            const int _SizeP) {
      if (candidates.size() == 1)
        return candidates.back();

      int index_worse_candidate = 0;
      bool improved = true;

      while (index_worse_candidate < (int)candidates.size() - 1 && improved) {
        improved = false;

        typename Shape::Compare_by_max_bound comp;

        std::sort(candidates.begin() + index_worse_candidate,
                  candidates.end(),
                  comp);

        //refine the best one 
        improveBound(candidates.back(),
                     _SizeP, m_num_subsets,
                     m_options.min_points);

        int position_stop;

        //Take all those intersecting the best one, check for equal ones
        for (position_stop = candidates.size() - 2;
             position_stop >= index_worse_candidate;
             position_stop--) {
          if (candidates.back()->min_bound() >
            candidates.at(position_stop)->max_bound())
            break;//the intervals do not overlaps anymore

          if (candidates.at(position_stop)->max_bound()
              <= m_options.min_points)
            break;  //the following candidate doesnt have enough points!

          //if we reach this point, there is an overlap
          //  between best one and position_stop
          //so request refining bound on position_stop
          improved |= improveBound(candidates.at(position_stop),
                                   _SizeP,
                                   m_num_subsets,
                                   m_options.min_points);

          //test again after refined
          if (candidates.back()->min_bound() >
            candidates.at(position_stop)->max_bound()) 
            break;//the intervals do not overlaps anymore
        }

        index_worse_candidate = position_stop;
      }

      return candidates.back();
    }

    inline FT getRadiusSphere_from_level(int l_level_octree) {
      return m_max_radiusSphere_Octree/powf(2, l_level_octree);
    }

    bool improveBound(Shape *candidate,
                      const int _SizeP,
                      size_t max_subset,
                      size_t min_points) {
      if (candidate->m_nb_subset_used >= max_subset)
        return false;

      if (candidate->m_nb_subset_used >= m_num_subsets)
        return false;

      candidate->m_nb_subset_used =
        (candidate->m_nb_subset_used >= m_num_subsets) ? 
        m_num_subsets - 1 : candidate->m_nb_subset_used;

      //what it does is add another subset and recompute lower and upper bound
      //the next subset to include is provided by m_nb_subset_used

      size_t numPointsEvaluated = 0;
      for (size_t i=0;i<candidate->m_nb_subset_used;i++)
        numPointsEvaluated += m_availableOctreeSizes[i];

      // need score of new subset as well as sum of
      // the score of the previous considered subset
      size_t newScore = 0;
      size_t newSampledPoints = 0;

      do {
        newScore = m_direct_octrees[candidate->m_nb_subset_used]->score(
          candidate, 
          m_shapeIndex, 
          m_options.epsilon,
          m_options.normal_threshold);

        candidate->m_score += newScore;
        
        numPointsEvaluated += 
          m_availableOctreeSizes[candidate->m_nb_subset_used];

        newSampledPoints +=
          m_availableOctreeSizes[candidate->m_nb_subset_used];

        candidate->m_nb_subset_used++;
      } while (newSampledPoints < min_points &&
        candidate->m_nb_subset_used < m_num_subsets);

      candidate->m_score = candidate->m_indices.size();

      candidate->compute_bound(numPointsEvaluated, _SizeP);

      return true;
    }

    inline FT StopProbability(FT _sizeC, FT _np, FT _dC, FT _l) const {
      return (std::min<FT>)(std::pow(1.f - _sizeC / (_np * _l * 3), _dC), 1.);
    }

    static bool candComp(const Shape* a, const Shape* b) {
      return a->expectedValue() < b->expected_value();
    }
    
  private:
    Parameters m_options;

    std::mt19937 m_rng;

    Direct_octree **m_direct_octrees;
    Indexed_octree *m_global_octree;
    std::vector<int> m_availableOctreeSizes;
    size_t m_num_subsets;

    // maps index into points to assigned extracted primitive
    std::vector<int> m_shapeIndex;
    size_t m_numAvailablePoints; 

    //give the index of the subset of point i
    std::vector<int> m_index_subsets;

    std::vector<Shape *> m_extractedShapes;

    std::vector<internal::Shape_factory_base *> m_shapeFactories;

    // iterators of input data
    Input_iterator m_inputIterator_first, m_inputIterator_beyond; 
    Point_pmap m_point_pmap;
    Normal_pmap m_normal_pmap;

    FT m_max_radiusSphere_Octree;
    std::vector<FT> m_level_weighting;  // sum must be 1
  };
}


#endif
