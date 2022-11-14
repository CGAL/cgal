// Copyright (c) 2015 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sven Oesau, Yannick Verdie, Clément Jamin, Pierre Alliez
//

#ifndef CGAL_SHAPE_DETECTION_EFFICIENT_RANSAC_H
#define CGAL_SHAPE_DETECTION_EFFICIENT_RANSAC_H

#include <CGAL/license/Shape_detection.h>

#include <CGAL/Random.h>

#include <CGAL/Shape_detection/Efficient_RANSAC/Octree.h>
#include <CGAL/Shape_detection/Efficient_RANSAC/Shape_base.h>
#include <CGAL/Shape_detection/Efficient_RANSAC/Plane.h>

// for octree ------------------------------
#include <boost/iterator/filter_iterator.hpp>
#include <CGAL/bounding_box.h>
#include <CGAL/Iterator_range.h>
//----------

#include <vector>
#include <cmath>
#include <limits>
#include <fstream>
#include <sstream>
#include <functional>

// boost --------------
#include <CGAL/boost/iterator/counting_iterator.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
//---------------------

namespace CGAL {
namespace Shape_detection {

/*!
  \ingroup PkgShapeDetectionRANSAC

  \brief Shape detection algorithm based on the RANSAC method.

  Given a point set in 3D space with unoriented normals, sampled on surfaces,
  this class enables to detect subsets of connected points lying on the surface of primitive shapes.
  Each input point is assigned to either none or at most one detected primitive
  shape. The implementation follows \cgalCite{schnabel2007efficient}.

  \tparam Traits must be a model of `EfficientRANSACTraits`.
*/
template<class Traits>
class Efficient_RANSAC {
public:

  /// \cond SKIP_IN_MANUAL
  struct Filter_unassigned_points {
    Filter_unassigned_points() : m_shape_index(dummy) {}

    Filter_unassigned_points(const std::vector<int> &shapeIndex)
            : m_shape_index(shapeIndex) {}

    bool operator()(std::size_t x) {
      if (x < m_shape_index.size())
        return m_shape_index[x] == -1;
      else return true; // to prevent infinite incrementing
    }

    const std::vector<int> &m_shape_index;
    std::vector<int> dummy;
  };

  typedef boost::filter_iterator<Filter_unassigned_points,
          boost::counting_iterator<std::size_t, boost::use_default, std::ptrdiff_t> > Point_index_iterator;
  ///< iterator for indices of points.
  /// \endcond

  /// \name Types
  /// @{
  /// \cond SKIP_IN_MANUAL
  typedef typename Traits::Input_range::iterator Input_iterator;
  typedef typename Traits::FT FT; ///< number type.
  typedef typename Traits::Point_3 Point; ///< point type.
  typedef typename Traits::Vector_3 Vector; ///< vector type.
  /// \endcond

  typedef typename Traits::Input_range Input_range;
  ///< Model of the concept `Range` with random access iterators, providing input points and normals
  /// through the following two property maps.

  typedef typename Traits::Point_map Point_map;
  ///< Property map to access the location of an input point.
  typedef typename Traits::Normal_map Normal_map;
  ///< Property map to access the unoriented normal of an input point.
  typedef Shape_base<Traits> Shape; ///< Shape type.
  typedef Plane<Traits> Plane_shape; ///< %Plane shape type.

#ifdef DOXYGEN_RUNNING
  typedef unspecified_type Shape_range;
  ///< `Iterator_range` with a bidirectional constant iterator type with value type `boost::shared_ptr<Shape>`.
  typedef unspecified_type Plane_range;
  ///< `Iterator_range` with a bidirectional constant iterator type with value type `boost::shared_ptr<Plane_shape>`.
#else

  struct Shape_range : public Iterator_range<
          typename std::vector<boost::shared_ptr<Shape> >::const_iterator> {
    typedef Iterator_range<
            typename std::vector<boost::shared_ptr<Shape> >::const_iterator> Base;

    Shape_range(boost::shared_ptr<std::vector<boost::shared_ptr<Shape> > >
                extracted_shapes) : Base(make_range(extracted_shapes->begin(),
                                                    extracted_shapes->end())), m_extracted_shapes(extracted_shapes) {}

  private:
    boost::shared_ptr<std::vector<boost::shared_ptr<Shape> > >
            m_extracted_shapes; // keeps a reference to the shape vector
  };

  struct Plane_range : public Iterator_range<
          typename std::vector<boost::shared_ptr<Plane_shape> >::const_iterator> {
    typedef Iterator_range<
            typename std::vector<boost::shared_ptr<Plane_shape> >::const_iterator> Base;

    Plane_range(boost::shared_ptr<std::vector<boost::shared_ptr<Plane_shape> > >
                extracted_shapes) : Base(make_range(extracted_shapes->begin(),
                                                    extracted_shapes->end())), m_extracted_shapes(extracted_shapes) {}

  private:
    boost::shared_ptr<std::vector<boost::shared_ptr<Plane_shape> > >
            m_extracted_shapes; // keeps a reference to the shape vector
  };

#endif

#ifdef DOXYGEN_RUNNING
  typedef unspecified_type Point_index_range;
  ///< `Iterator_range` with a bidirectional iterator with value type `std::size_t`
  ///  as indices into the input data that has not been assigned to a shape.
  ///  As this range class has no `size()` method, the method
  ///  `Efficient_RANSAC::number_of_unassigned_points()` is provided.
#else
  typedef Iterator_range<Point_index_iterator>
          Point_index_range;
#endif

  /// @}

  /// \name Parameters
  /// @{
  /*!
   Parameters for the shape detection algorithm. They are explained in detail
   in Section \ref Shape_detection_RANSACParameters  of the User Manual.
   */
  struct Parameters {
    Parameters()
            : probability((FT) 0.01), min_points((std::numeric_limits<std::size_t>::max)()), epsilon(-1),
              normal_threshold((FT) 0.9), cluster_epsilon(-1) {}

    /*!
      Probability to control search endurance.
      %Default value is 0.05.

      A lower probability provides a higher reliability and determinism at the cost
      of longer running time due to a higher search endurance.

      It must belong to the interval [0, 1].
    */
    FT probability;

    /*!
      Minimum number of points in a shape.
      %Default value is 1% of total number of input points.

      It must belong to the interval [0, +inf).
    */
    std::size_t min_points;

    /*!
      Maximum acceptable Euclidean distance between a point and a shape.
      %Default value is 1% of the bounding box diagonal.

      It must belong to the interval [0, +inf).
    */
    FT epsilon;

    /*!
      Maximum threshold on the dot product between the estimated
      shape's normal and the point's normal, that is the cosine of the angle (cos(25°) = 0.9).
      %Default value is 0.9 (around 25 degrees).

      It must belong to the interval [0, 1].
    */
    FT normal_threshold;

    /*!
      Maximum acceptable Euclidean distance between points, which are assumed to be neighbors.
      %Default value is 1% of the bounding box diagonal.

      It must belong to the interval [0, +inf).
    */
    FT cluster_epsilon;
  };
  /// @}

private:

  typedef internal::RANSAC_octree<Traits> Direct_octree;
  typedef internal::RANSAC_octree<Traits> Indexed_octree;

  //--------------------------------------------typedef

  // Creates a function pointer for instancing shape instances.
  template<class ShapeT>
  static Shape *factory() {
    return new ShapeT;
  }

public:

  /// \name Initialization
  /// @{

  /*!
    Constructs an empty shape detection object.
  */
  Efficient_RANSAC(Traits t = Traits())
          : m_traits(t), m_direct_octrees(nullptr), m_global_octree(nullptr), m_num_subsets(0),
            m_num_available_points(0), m_num_total_points(0), m_valid_iterators(false) {
  }

  /*!
    Releases all memory allocated by this instance including shapes.
  */
  ~Efficient_RANSAC() {
    clear();
  }

  /*!
    Retrieves the traits class.
   */
  const Traits &
  traits() const {
    return m_traits;
  }

  /*!
    Retrieves the point property map.
  */
  const Point_map &point_map() const { return m_point_pmap; }

  /*!
    Retrieves the normal property map.
  */
  const Normal_map &normal() const { return m_normal_pmap; }

  Input_iterator input_iterator_first() const {
    return m_input_iterator_first;
  }

  Input_iterator input_iterator_beyond() const {
    return m_input_iterator_beyond;
  }

  /*!
    Sets the input data. The range must stay valid
    until the detection has been performed and the access to the
    results is no longer required. The data in the input is reordered by the methods
    `detect()` and `preprocess()`. This function first calls `clear()`.
  */
  void set_input(
          Input_range &input_range,
          ///< Range of input data.
          Point_map point_map = Point_map(),
          ///< Property map to access the position of an input point.
          Normal_map normal_map = Normal_map()
          ///< Property map to access the normal of an input point.
  ) {
    m_point_pmap = point_map;
    m_normal_pmap = normal_map;

    m_input_iterator_first = input_range.begin();
    m_input_iterator_beyond = input_range.end();

    clear();

    m_extracted_shapes =
            boost::make_shared<std::vector<boost::shared_ptr<Shape> > >();

    m_num_available_points = m_num_total_points = std::distance(
            m_input_iterator_first, m_input_iterator_beyond);

    m_valid_iterators = true;

  }

  /*!
    Registers the shape type `ShapeType` in the detection engine that must inherit from `Shape_base`.
    For example, for registering a plane as detectable shape, you should call
    `ransac.add_shape_factory< Shape_detection::Plane<Traits> >();`. Note
    that if your call is within a template, you should add the `template`
    keyword just before `add_shape_factory`:
    `ransac.template add_shape_factory< Shape_detection::Plane<Traits> >();`.
  */
  template<class Shape_type>
  void add_shape_factory() {
    m_shape_factories.push_back(factory<Shape_type>);
  }

  /*!
    Constructs internal data structures required for the shape detection.
    These structures only depend on the input data, i.e. the points and
    normal vectors. This method is called by `detect()`, if it was not called
    before by the user.
  */
  bool preprocess() {

    if (m_num_total_points == 0)
      return false;

    // Generation of subsets
    m_num_subsets = (std::size_t) (std::max<std::ptrdiff_t>)((std::ptrdiff_t)
                                                                     std::floor(std::log(double(m_num_total_points)) /
                                                                                std::log(2.)) - 9, 2);

    // SUBSET GENERATION ->
    // approach with increasing subset sizes -> replace with octree later on
    Input_iterator last = m_input_iterator_beyond - 1;
    std::size_t remainingPoints = m_num_total_points;

    m_available_octree_sizes.resize(m_num_subsets);
    m_direct_octrees = new Direct_octree *[m_num_subsets];
    for (int s = int(m_num_subsets) - 1; s >= 0; --s) {
      std::size_t subsetSize = remainingPoints;
      std::vector<std::size_t> indices(subsetSize);
      if (s) {
        subsetSize >>= 1;
        for (std::size_t i = 0; i < subsetSize; i++) {
          std::size_t index = get_default_random()(2);
          index = index + (i << 1);
          index = (index >= remainingPoints) ? remainingPoints - 1 : index;
          indices[i] = index;
        }

        // move points to the end of the point vector
        std::size_t j = subsetSize;
        do {
          j--;
          typename std::iterator_traits<Input_iterator>::value_type
                  tmp = (*last);
          *last = m_input_iterator_first[indices[std::size_t(j)]];
          m_input_iterator_first[indices[std::size_t(j)]] = tmp;
          last--;
        } while (j > 0);
        m_direct_octrees[s] = new Direct_octree(
                m_traits, last + 1,
                last + subsetSize + 1,
                m_point_pmap,
                remainingPoints - subsetSize);
      } else
        m_direct_octrees[0] = new Direct_octree(
                m_traits, m_input_iterator_first,
                m_input_iterator_first + (subsetSize),
                m_point_pmap,
                0);

      m_available_octree_sizes[s] = subsetSize;
      m_direct_octrees[s]->refine(m_options.cluster_epsilon);

      remainingPoints -= subsetSize;
    }

    m_global_octree = new Indexed_octree(
            m_traits, m_input_iterator_first, m_input_iterator_beyond,
            m_point_pmap
    );
    m_global_octree->refine(m_options.cluster_epsilon);

    return true;
  }

  /// @}

  /// \name Memory Management
  /// @{
  /*!
    Removes all shape types registered for detection.
   */
  void clear_shape_factories() {
    m_shape_factories.clear();
  }

  /*!
    Frees memory allocated for the internal search structures but keeps the detected shapes.
    It invalidates the range retrieved using `unassigned_points()`.
   */
  void clear_octrees() {
    // If there is no data yet, there are no data structures.
    if (!m_valid_iterators)
      return;

    if (m_global_octree) {
      delete m_global_octree;
      m_global_octree = nullptr;
    }

    if (m_direct_octrees) {
      for (std::size_t i = 0; i < m_num_subsets; i++)
        delete m_direct_octrees[i];
      delete[] m_direct_octrees;

      m_direct_octrees = nullptr;
    }

    m_num_subsets = 0;
  }

  /*!
    Calls `clear_octrees()` and removes all detected shapes.
    All internal structures are cleaned, including formerly detected shapes.
    Thus iterators and ranges retrieved through `shapes()`, `planes()` and `indices_of_unassigned_points()`
    are invalidated.
  */
  void clear() {
    // If there is no data yet, there are no data structures.
    if (!m_valid_iterators)
      return;

    std::vector<int>().swap(m_shape_index);

    m_extracted_shapes =
            boost::make_shared<std::vector<boost::shared_ptr<Shape> > >();

    m_num_available_points = m_num_total_points;

    clear_octrees();
    clear_shape_factories();
  }
  /// @}

  /// \name Detection
  /// @{

  /*!
    Performs the shape detection. Shape types considered during the detection
    are those registered using `add_shape_factory()`.

    \param options parameters for shape detection

    \param callback can be omitted if the algorithm should be run
    without any callback. It is called regularly when the algorithm
    is running: the current advancement (between 0.0 and 1.0) is
    passed as parameter. If it returns `true`, then the algorithm
    continues its execution normally; if it returns `false`, the
    algorithm is stopped. Note that this interruption may leave the
    class in an invalid state.

    \return `true` if shape types have been registered and
            input data has been set. Otherwise, `false` is returned.
  */
  bool detect(const Parameters &options = Parameters(),
              const std::function<bool(double)> &callback
              = std::function<bool(double)>()) {

    m_options = options;

    // No shape types for detection or no points provided, exit
    if (m_shape_factories.size() == 0 ||
        (m_input_iterator_beyond - m_input_iterator_first) == 0)
      return false;

    if (m_num_subsets == 0 || m_global_octree == 0) {
      if (!preprocess())
        return false;
    }

    if (callback && !callback(0.)) {
      clear_octrees();
      clear_shape_factories();
      return false;
    }

    // Reset data structures possibly used by former search
    m_extracted_shapes =
            boost::make_shared<std::vector<boost::shared_ptr<Shape> > >();
    m_num_available_points = m_num_total_points;

    for (std::size_t i = 0; i < m_num_subsets; i++) {
      m_available_octree_sizes[i] = m_direct_octrees[i]->size();
    }

    // Use bounding box diagonal as reference for default values
    Bbox_3 bbox = m_global_octree->boundingBox();
    FT bbox_diagonal = (FT) CGAL::sqrt(
            (bbox.xmax() - bbox.xmin()) * (bbox.xmax() - bbox.xmin())
            + (bbox.ymax() - bbox.ymin()) * (bbox.ymax() - bbox.ymin())
            + (bbox.zmax() - bbox.zmin()) * (bbox.zmax() - bbox.zmin()));

    // Epsilon or cluster_epsilon have been set by the user?
    // If not, derive from bounding box diagonal
    m_options.epsilon = (m_options.epsilon < 0)
                        ? bbox_diagonal * (FT) 0.01 : m_options.epsilon;

    m_options.cluster_epsilon = (m_options.cluster_epsilon < 0)
                                ? bbox_diagonal * (FT) 0.01 : m_options.cluster_epsilon;

    // Minimum number of points has been set?
    m_options.min_points =
      (m_options.min_points == (std::numeric_limits<std::size_t>::max)()) ?
      (std::size_t)((FT)0.01 * m_num_available_points) :
      m_options.min_points;
    m_options.min_points = (m_options.min_points < 10) ? 10 : m_options.min_points;

    // Initializing the shape index
    m_shape_index.assign(m_num_available_points, -1);

    if (m_options.min_points > m_num_available_points)
      return true;

    // List of all randomly drawn candidates
    // with the minimum number of points
    std::vector<Shape *> candidates;

    // Identifying minimum number of samples
    m_required_samples = 0;
    for (std::size_t i = 0; i < m_shape_factories.size(); i++) {
      Shape *tmp = (Shape *) m_shape_factories[i]();
      m_required_samples = (std::max<std::size_t>)(m_required_samples, tmp->minimum_sample_size());
      delete tmp;
    }

    std::size_t first_sample; // first sample for RANSAC

    FT best_expected = 0;

    // number of points that have been assigned to a shape
    std::size_t num_invalid = 0;

    std::size_t generated_candidates = 0;
    std::size_t failed_candidates = 0;
    std::size_t limit_failed_candidates = (std::max)(std::size_t(10000),
                                                     std::size_t(m_input_iterator_beyond
                                                                 - m_input_iterator_first)
                                                     / std::size_t(100));

    bool force_exit = false;
    bool keep_searching = true;

    do { // main loop
      best_expected = 0;

      if (keep_searching)
        do {
            // Search (remaining_points / min_points) shapes (max 200 per iteration, min 1)
            std::size_t search_number
              = (std::min)(std::size_t(200),
                           (std::max)(std::size_t((m_num_available_points - num_invalid) / double(m_options.min_points)),
                                      std::size_t(1)));
            for (std::size_t nb = 0; nb < search_number; ++ nb)
            {
              // Generate candidates
              //1. pick a point p1 randomly among available points
              std::set<std::size_t> indices;
              bool done = false;
              do {
                do
                  first_sample = get_default_random()(
                    static_cast<unsigned int>(m_num_available_points));
                while (m_shape_index[first_sample] != -1);

                done = drawSamplesFromCellContainingPoint
                  (m_global_octree,
                   get(m_point_pmap,
                       *(m_input_iterator_first + first_sample)),
                   select_random_octree_level(),
                   indices,
                   m_shape_index,
                   m_required_samples);

                if (callback && !callback(num_invalid / double(m_num_total_points))) {
                  clear(num_invalid, candidates);
                  return false;
                }

              } while (m_shape_index[first_sample] != -1 || !done);

              generated_candidates++;

              //add candidate for each type of primitives
              bool candidate_success = false;
              for(typename std::vector<Shape *(*)()>::iterator it =
                    m_shape_factories.begin(); it != m_shape_factories.end(); it++)        {
                if (callback && !callback(num_invalid / double(m_num_total_points))) {
                  clear(num_invalid, candidates);
                  return false;
                }
                Shape *p = (Shape *) (*it)();
                //compute the primitive and says if the candidate is valid
                p->compute(indices,
                           m_input_iterator_first,
                           m_traits,
                           m_point_pmap,
                           m_normal_pmap,
                           m_options.epsilon,
                           m_options.normal_threshold);

                if (p->is_valid()) {
                  improve_bound(p, m_num_available_points - num_invalid, 1, 500);

                  //evaluate the candidate
                  if(p->max_bound() >= m_options.min_points && p->score() > 0) {
                    if (best_expected < p->expected_value())
                      best_expected = p->expected_value();

                    candidates.push_back(p);
                    candidate_success = true;
                  }
                  else {
                    delete p;
                  }
                }
                else {
                  delete p;
                }
              }
              if (!candidate_success)
                ++ failed_candidates;

            }

            if (failed_candidates >= limit_failed_candidates)
            {
              force_exit = true;
            }

          keep_searching = (stop_probability(m_options.min_points,
                                             m_num_available_points - num_invalid,
                                             generated_candidates, m_global_octree->maxLevel())
                            > m_options.probability);
        } while (!force_exit
                 && stop_probability((std::size_t) best_expected,
                                     m_num_available_points - num_invalid,
                                     generated_candidates,
                                     m_global_octree->maxLevel())
                    > m_options.probability
                 && keep_searching);
        // end of generate candidate

      if (force_exit) {
        break;
      }

      if (candidates.empty())
        continue;

      // Now get the best candidate in the current set of all candidates
      // Note that the function sorts the candidates:
      //  the best candidate is always the last element of the vector

      Shape *best_candidate =
              get_best_candidate(candidates, m_num_available_points - num_invalid);

      if (callback && !callback(num_invalid / double(m_num_total_points))) {
        clear(num_invalid, candidates);
        return false;
      }

      // If search is done and the best candidate is too small, we are done.
      if (!keep_searching && best_candidate->m_score < m_options.min_points)
        break;

      if (!best_candidate)
        continue;

      best_candidate->m_indices.clear();

      best_candidate->m_score =
              score(m_global_octree,
                    best_candidate,
                    m_shape_index,
                    FT(3) * m_options.epsilon,
                    m_options.normal_threshold);

      best_expected = static_cast<FT>(best_candidate->m_score);

      best_candidate->connected_component(best_candidate->m_indices,
                                          m_options.cluster_epsilon);

      if (callback && !callback(num_invalid / double(m_num_total_points))) {
        clear(num_invalid, candidates);
        return false;
      }
      // check score against min_points and clear out candidates if too low
      if (best_candidate->indices_of_assigned_points().size() <
          m_options.min_points) {
        if (!(best_candidate->indices_of_assigned_points().empty()))
          for (std::size_t i = 0; i < candidates.size() - 1; i++) {
            if (best_candidate->is_same(candidates[i])) {
              delete candidates[i];
              candidates[i] = nullptr;
            }
          }

        candidates.back() = nullptr;
        delete best_candidate;
        best_candidate = nullptr;

        if (callback && !callback(num_invalid / double(m_num_total_points))) {
          clear(num_invalid, candidates);
          return false;
        }

        // Trimming candidates list
        std::size_t empty = 0, occupied = 0;
        while (empty < candidates.size()) {
          while (empty < candidates.size() && candidates[empty]) empty++;

          if (empty >= candidates.size())
            break;

          if (occupied < empty)
            occupied = empty + 1;

          while (occupied < candidates.size() && !candidates[occupied])
            occupied++;

          if (occupied >= candidates.size())
            break;
          candidates[empty] = candidates[occupied];
          candidates[occupied] = nullptr;
          empty++;
          occupied++;
        }

        candidates.resize(empty);

        if (callback && !callback(num_invalid / double(m_num_total_points))) {
          clear(num_invalid, candidates);
          return false;
        }
      } else if (stop_probability((std::size_t) best_candidate->expected_value(),
                                  (m_num_available_points - num_invalid),
                                  generated_candidates,
                                  m_global_octree->maxLevel())
                 <= m_options.probability) {

        // Remove candidate from list
        candidates.back() = nullptr;

        //1. add best candidate to final result.
        m_extracted_shapes->push_back(
                boost::shared_ptr<Shape>(best_candidate));

        if (callback && !callback(num_invalid / double(m_num_total_points))) {
          clear(num_invalid, candidates);
          return false;
        }

        //2. remove the points
        const std::vector<std::size_t> &indices_points_best_candidate =
                best_candidate->indices_of_assigned_points();

        // update generated candidates to reflect removal of points
        generated_candidates = std::size_t(std::pow(1.f - (indices_points_best_candidate.size() /
                                                           float(m_num_available_points - num_invalid)), 3.f)
                                           * generated_candidates);

        //2.3 Remove the points from the subtrees
        for (std::size_t i = 0; i < indices_points_best_candidate.size(); i++) {
          m_shape_index[indices_points_best_candidate.at(i)] =
                  int(m_extracted_shapes->size()) - 1;

          num_invalid++;

          for (std::size_t j = 0; j < m_num_subsets; j++) {
            if (m_direct_octrees[j]) {
              std::size_t offset = m_direct_octrees[j]->offset();

              if (offset <= indices_points_best_candidate.at(i) &&
                  (indices_points_best_candidate.at(i) - offset)
                  < m_direct_octrees[j]->size()) {
                m_available_octree_sizes[j]--;
              }
            }
          }
        }

        failed_candidates = 0;
        best_expected = 0;

        if (callback && !callback(num_invalid / double(m_num_total_points))) {
          clear(num_invalid, candidates);
          return false;
        }

        std::vector<std::size_t> subset_sizes(m_num_subsets);
        subset_sizes[0] = m_available_octree_sizes[0];
        for (std::size_t i = 1; i < m_num_subsets; i++) {
          subset_sizes[i] = subset_sizes[i - 1] + m_available_octree_sizes[i];
        }


        //3. Remove points from candidates common with extracted primitive
        //#pragma omp parallel for
        best_expected = 0;
        for (std::size_t i = 0; i < candidates.size() - 1; i++) {
          if (candidates[i]) {
            candidates[i]->update_points(m_shape_index);
            candidates[i]->compute_bound(
                    subset_sizes[candidates[i]->m_nb_subset_used - 1],
                    m_num_available_points - num_invalid);

            if (candidates[i]->max_bound() < m_options.min_points) {
              delete candidates[i];
              candidates[i] = nullptr;
            } else {
              best_expected = (candidates[i]->expected_value() > best_expected) ?
                              candidates[i]->expected_value() : best_expected;
            }
          }
        }

        if (callback && !callback(num_invalid / double(m_num_total_points))) {
          clear(num_invalid, candidates);
          return false;
        }

        std::size_t start = 0, end = candidates.size() - 1;
        while (start < end) {
          while (candidates[start] && start < end) start++;
          while (!candidates[end] && start < end) end--;
          if (!candidates[start] && candidates[end] && start < end) {
            candidates[start] = candidates[end];
            candidates[end] = nullptr;
            start++;
            end--;
          }
        }

        if (candidates[end]) end++;

        candidates.resize(end);
      } else if (!keep_searching)
        ++generated_candidates;

      if (callback && !callback(num_invalid / double(m_num_total_points))) {
        clear(num_invalid, candidates);
        return false;
      }

      keep_searching = (stop_probability(m_options.min_points,
                                         m_num_available_points - num_invalid,
                                         generated_candidates,
                                         m_global_octree->maxLevel())
                        > m_options.probability);
    } while ((keep_searching
              && FT(m_num_available_points - num_invalid) >= m_options.min_points)
             || best_expected >= m_options.min_points);

    // Clean up remaining candidates.
    clear_candidates(num_invalid, candidates);
    return true;
  }

  /// @}

  /// \name Access
  /// @{
  /*!
    Returns an `Iterator_range` with a bidirectional iterator with value type
    `boost::shared_ptr<Shape>` over the detected shapes in the order of detection.
    Depending on the chosen probability
    for the detection, the shapes are ordered with decreasing size.
  */
  Shape_range shapes() const {
    return Shape_range(m_extracted_shapes);
  }

  /*!
    Returns an `Iterator_range` with a bidirectional iterator with
    value type `boost::shared_ptr<Plane_shape>` over only the
    detected planes in the order of detection.  Depending on the
    chosen probability for the detection, the planes are ordered
    with decreasing size.
  */
  Plane_range planes() const {
    boost::shared_ptr<std::vector<boost::shared_ptr<Plane_shape> > > planes
            = boost::make_shared<std::vector<boost::shared_ptr<Plane_shape> > >();

    for (std::size_t i = 0; i < m_extracted_shapes->size(); ++i) {
      boost::shared_ptr<Plane_shape> pshape
              = boost::dynamic_pointer_cast<Plane_shape>((*m_extracted_shapes)[i]);

      // Ignore all shapes other than plane
      if (pshape != boost::shared_ptr<Plane_shape>())
        planes->push_back(pshape);
    }
    return Plane_range(planes);
  }

  /*!
    Number of points not assigned to a shape.
  */
  std::size_t number_of_unassigned_points() const {
    return m_num_available_points;
  }

  /*!
    Returns an `Iterator_range` with a bidirectional iterator with value type `std::size_t`
    as indices into the input data that has not been assigned to a shape.
  */
  Point_index_range indices_of_unassigned_points() {
    Filter_unassigned_points fup(m_shape_index);

    Point_index_iterator p1 =
            boost::make_filter_iterator<Filter_unassigned_points>(
                    fup,
                    boost::counting_iterator<std::size_t, boost::use_default, std::ptrdiff_t>(0),
                    boost::counting_iterator<std::size_t, boost::use_default, std::ptrdiff_t>(m_shape_index.size()));

    return make_range(p1, Point_index_iterator(p1.end()));
  }
  /// @}

private:
  void clear(
    const std::size_t num_invalid, std::vector<Shape *>& candidates) {

    clear_octrees();
    clear_shape_factories();
    clear_candidates(num_invalid, candidates);
  }

  void clear_candidates(
    const std::size_t num_invalid, std::vector<Shape *>& candidates) {

    for (std::size_t i = 0; i < candidates.size(); i++) {
      delete candidates[i];
    }
    candidates.resize(0);
    m_num_available_points -= num_invalid;
  }

  int select_random_octree_level() {
    auto upper_bound = static_cast<unsigned int>(m_global_octree->maxLevel() + 1);
    return (int) get_default_random()(upper_bound);
  }

  Shape *get_best_candidate(std::vector<Shape *> &candidates,
                            const std::size_t num_available_points) {

    if (candidates.size() == 1)
      return candidates.back();

    int index_worse_candidate = 0;
    bool improved = true;

    while (index_worse_candidate < (int) candidates.size() - 1 && improved) {
      improved = false;

      typename Shape::Compare_by_max_bound comp;

      std::sort(candidates.begin() + index_worse_candidate,
                candidates.end(),
                comp);

      //refine the best one
      improve_bound(candidates.back(),
                    num_available_points, m_num_subsets,
                    m_options.min_points);

      int position_stop;

      //Take all those intersecting the best one, check for equal ones
      for (position_stop = int(candidates.size()) - 1;
           position_stop > index_worse_candidate;
           position_stop--) {
        if (candidates.back()->min_bound() >
            candidates.at(position_stop)->max_bound())
          break;//the intervals do not overlaps anymore

        if (candidates.at(position_stop)->max_bound()
            <= m_options.min_points)
          break;  //the following candidate doesn't have enough points!

        //if we reach this point, there is an overlap
        //  between best one and position_stop
        //so request refining bound on position_stop
        improved |= improve_bound(candidates.at(position_stop),
                                  num_available_points,
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

  bool improve_bound(Shape *candidate,
                     std::size_t num_available_points,
                     std::size_t max_subset,
                     std::size_t min_points) {

    if (candidate->m_nb_subset_used >= max_subset)
      return false;

    if (candidate->m_nb_subset_used >= m_num_subsets)
      return false;

    candidate->m_nb_subset_used =
            (candidate->m_nb_subset_used >= m_num_subsets) ?
            m_num_subsets - 1 : candidate->m_nb_subset_used;

    //what it does is add another subset and recompute lower and upper bound
    //the next subset to include is provided by m_nb_subset_used

    std::size_t num_points_evaluated = 0;
    for (std::size_t i = 0; i < candidate->m_nb_subset_used; i++)
      num_points_evaluated += m_available_octree_sizes[i];

    // need score of new subset as well as sum of
    // the score of the previous considered subset
    std::size_t new_score = 0;
    std::size_t new_sampled_points = 0;

    do {
      new_score =
              score(m_direct_octrees[candidate->m_nb_subset_used],
                    candidate,
                    m_shape_index,
                    m_options.epsilon,
                    m_options.normal_threshold);

      candidate->m_score += new_score;

      num_points_evaluated +=
              m_available_octree_sizes[candidate->m_nb_subset_used];

      new_sampled_points +=
              m_available_octree_sizes[candidate->m_nb_subset_used];

      candidate->m_nb_subset_used++;
    } while (new_sampled_points < min_points &&
             candidate->m_nb_subset_used < m_num_subsets);

    candidate->m_score = candidate->m_indices.size();

    candidate->compute_bound(num_points_evaluated, num_available_points);

    return true;
  }

  inline FT stop_probability(std::size_t largest_candidate, std::size_t num_pts, std::size_t num_candidates, std::size_t octree_depth) const {
    return (std::min<FT>)((FT)std::pow(FT(1) - FT(largest_candidate)
                                       / (FT(num_pts) * FT(octree_depth+1)
                                          * FT(1 << (m_required_samples - 1))),
                                       int(num_candidates)), FT(1));
  }

  template<class Octree>
  std::size_t score(const Octree *octree,
                    Shape *candidate,
                    std::vector<int> &shapeIndex,
                    FT epsilon,
                    FT normal_threshold) {

    typedef typename Octree::Node Cell;

    std::stack<Cell> stack;
    stack.push(octree->root());

    while (!stack.empty()) {
      Cell cell = stack.top();
      stack.pop();

      FT width = octree->width() / (1 << (cell.depth()));

      FT diag = CGAL::sqrt(FT(3) * width * width) + epsilon;

      FT dist = candidate->squared_distance(octree->barycenter(cell));

      if (dist > (diag * diag))
        continue;

      // differ between full or partial overlap?
      // if full overlap further traversal of this branch is not necessary
      if (cell.is_leaf()) {
        std::vector<std::size_t> indices;
        indices.reserve(cell.size());
        for (std::size_t i = 0; i < cell.size(); i++) {
          if (shapeIndex[octree->index(cell, i)] == -1) {
            indices.push_back(octree->index(cell, i));
          }
        }

        candidate->cost_function(epsilon,
                                 normal_threshold,
                                 indices);
      } else {

        if (!cell.is_leaf()) {
          for (std::size_t i = 0; i < 8; i++) {
            if (!cell[i].empty())
              stack.push(cell[i]);
          }
        }
      }

    }

    return candidate->m_indices.size();
  }


  template<class Octree>
  const typename Octree::Node node_containing_point(const Octree *octree, const Point &p, std::size_t level) {

    // Find the node containing the point
    typename Octree::Node cur = octree->root();
    while (!cur.is_null() && cur.depth() < level) {

      // If cur is a leaf node, its child is null
      if (cur.is_leaf())
        return typename Octree::Node();

      // If that child is empty, return null
      if (cur.empty())
        return typename Octree::Node();

      // Determine the coordinate of the child
      Point center = octree->barycenter(cur);
      std::bitset<3> coordinate;
      coordinate[0] = center.x() <= p.x();
      coordinate[1] = center.y() <= p.y();
      coordinate[2] = center.z() <= p.z();

      // Otherwise, return the correct child of cur
      cur = cur[coordinate.to_ulong()];

    }

    return cur;
  }

  template<class Octree>
  bool drawSamplesFromCellContainingPoint(const Octree *octree,
                                          const Point &p,
                                          std::size_t level,
                                          std::set<std::size_t> &indices,
                                          const std::vector<int> &shapeIndex,
                                          std::size_t requiredSamples) {

    typedef typename Octree::Node Cell;

    const Cell cur = node_containing_point(octree, p, level);

    // Stop if the node we need doesn't exist
    if (cur.is_null())
      return false;

    // Count point indices that map to -1 in the shape index
    std::size_t enough = 0;
    for (auto j : cur) {
      if (shapeIndex[j] == -1)
        enough++;
      if (enough >= requiredSamples)
        break;
    }

    // Make sure we found enough samples
    if (enough < requiredSamples)
      return false;

    do {
      std::size_t p = CGAL::get_default_random().
              uniform_int<std::size_t>(0, cur.size() - 1);
      std::size_t j = octree->index(cur, p);

      if (shapeIndex[j] == -1)
        indices.insert(j);
    } while (indices.size() < requiredSamples);

    return true;
  }

private:
  Parameters m_options;

  // Traits class.
  Traits m_traits;

  // Octrees build on input data for quick shape evaluation and
  // sample selection within an octree cell.
  Direct_octree **m_direct_octrees;
  Indexed_octree *m_global_octree;
  std::vector<std::size_t> m_available_octree_sizes;
  std::size_t m_num_subsets;

  // maps index into points to assigned extracted primitive
  std::vector<int> m_shape_index;
  std::size_t m_num_available_points;
  std::size_t m_num_total_points;
  std::size_t m_required_samples;

  //give the index of the subset of point i
  std::vector<int> m_index_subsets;

  boost::shared_ptr<std::vector<boost::shared_ptr<Shape> > > m_extracted_shapes;

  std::vector<Shape *(*)()> m_shape_factories;

  // iterators of input data
  bool m_valid_iterators;
  Input_iterator m_input_iterator_first, m_input_iterator_beyond;
  Point_map m_point_pmap;
  Normal_map m_normal_pmap;
};


}

}

#endif // CGAL_SHAPE_DETECTION_EFFICIENT_RANSAC_H
