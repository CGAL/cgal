// Copyright (c) 2005,2008,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Shlomo Golubev   <golubevs@post.tau.ac.il>

#ifndef CGAL_ARR_LANDMARKS_SPECIFIED_POINTS_GENERATOR_H
#define CGAL_ARR_LANDMARKS_SPECIFIED_POINTS_GENERATOR_H

#include <CGAL/license/Arrangement_on_surface_2.h>


/*! \file
 * Definition of the Arr_lm_specified_points_generator<Arrangement> template.
 */

#include <CGAL/Arr_point_location/Arr_lm_generator_base.h>

namespace CGAL {

/*! \class Arr_landmarks_specified_points_generator
 * A generator for the landmarks point-locatoion class, which uses
 * specified set of points as its set of landmarks.
*/
template <typename Arrangement_,
          typename Nearest_neighbor_ =
            Arr_landmarks_nearest_neighbor<Arrangement_> >
class Arr_landmarks_specified_points_generator :
    public Arr_landmarks_generator_base<Arrangement_, Nearest_neighbor_>
{
public:
  typedef Arrangement_                                  Arrangement_2;
  typedef Nearest_neighbor_                             Nearest_neighbor;

private:
  typedef Arr_landmarks_generator_base<Arrangement_2, Nearest_neighbor>
                                                        Base;
  typedef Arr_landmarks_vertices_generator<Arrangement_2, Nearest_neighbor>
                                                        Self;

public:
  typedef typename Arrangement_2::Geometry_traits_2     Geometry_traits_2;
  typedef typename Arrangement_2::Vertex_const_handle   Vertex_const_handle;
  typedef typename Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Arrangement_2::Face_const_handle     Face_const_handle;
  typedef typename Arrangement_2::Vertex_handle         Vertex_handle;
  typedef typename Arrangement_2::Halfedge_handle       Halfedge_handle;
  typedef typename Arrangement_2::Face_handle           Face_handle;
  typedef typename Arrangement_2::Vertex_const_iterator Vertex_const_iterator;

  typedef typename Arrangement_2::Point_2               Point_2;
  typedef typename Arrangement_2::X_monotone_curve_2    X_monotone_curve_2;

  typedef typename Base::Points_set                     Points_set;

  typedef typename Nearest_neighbor::NN_Point_2         NN_Point_2;
  typedef std::list<NN_Point_2>                         NN_Point_list;

protected:
  typedef Arr_traits_basic_adaptor_2<Geometry_traits_2> Traits_adaptor_2;
  typedef typename Base::PL_result_type                 PL_result_type;
  typedef std::pair<Point_2, PL_result_type>            PL_pair;
  typedef std::vector<PL_pair>                          Pairs_set;

  // Data members:
  const Traits_adaptor_2* m_traits;  // Its associated traits object.
  Points_set              m_points;  // container of the specified points
  Pairs_set               lm_pairs;
  int                     num_landmarks;

private:
  /*! Copy constructor - not supported. */
  Arr_landmarks_specified_points_generator(const Self&);

  /*! Assignment operator - not supported. */
  Self& operator=(const Self&);

public:
  /*! Constructor.
   * Create landmarks in the points that are given to it.
   */
  Arr_landmarks_specified_points_generator(const Arrangement_2& arr,
                                           const Points_set points) :
    Base(arr),
    m_traits(static_cast<const Traits_adaptor_2*>(arr.geometry_traits())),
    m_points(points),
    num_landmarks(points.size())
  { build_landmark_set(); }

  /*! Constructor. from an arrangement.
   * \param arr (in) The arrangement.
   */
  Arr_landmarks_specified_points_generator(const Arrangement_2& arr) :
    Base(arr),
    m_traits(static_cast<const Traits_adaptor_2*> (arr.geometry_traits())),
    num_landmarks(1)
  {
    //this constructor creates a single landmark in the origin
    m_points.push_back(Point_2(0,0));
    build_landmark_set();
  }

  virtual void _create_points_set(Points_set & /* points */)
  {
    std::cerr << "should not reach here!"<< std::endl;
    CGAL_error();
  }

  /*! Create the landmark set, using all arrangement vertices.
   */
  void build_landmark_set()
  {
    lm_pairs.clear();
    locate(*(this->arrangement()), m_points.begin(), m_points.end(),
           std::back_inserter(lm_pairs));

    // Go over the container of the specified points and insert them as
    // landmarks.
    NN_Point_list                   nnp_list;
    typename Points_set::iterator   pt_it;
    typename Pairs_set::iterator    pairs_it;
    for (pt_it = m_points.begin(); pt_it != m_points.end(); ++pt_it) {
      for (pairs_it = lm_pairs.begin();
           pairs_it != lm_pairs.end() && (*pairs_it).first != (*pt_it);
           ++pairs_it) {};
      if ((*pairs_it).first == (*pt_it)) {
        nnp_list.push_back (NN_Point_2 ((*pt_it),(*pairs_it).second));
      }
    }

    // Update the search structure.
    this->nn.clear();
    this->nn.init (nnp_list.begin(), nnp_list.end());

    this->num_small_not_updated_changes = 0;
    this->updated = true;
  }

  /*!
   * Clear the landmark set.
   */
  void clear_landmark_set()
  {
    this->nn.clear();

    num_landmarks = 0;
    this->num_small_not_updated_changes = 0;
    this->updated = false;
  }

  /*!
   * Get the nearest neighbor (landmark) to the given point.
   * \param q The query point.
   * \param obj Output: The location of the nearest landmark point in the
   *                    arrangement (a vertex, halfedge, or face handle).
   * \return The nearest landmark point.
   */
  virtual Point_2 closest_landmark(const Point_2& q, PL_result_type& obj)
  {
    CGAL_assertion(this->updated);
    return (this->nn.find_nearest_neighbor(q, obj));
  }
};

} //namespace CGAL

#endif
