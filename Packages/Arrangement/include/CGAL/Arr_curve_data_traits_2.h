// Copyright (c) 1997  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$
// $Name$
//
// Author(s)     : Ron Wein          <wein@post.tau.ac.il>
//                 Efi Fogel         <efif@post.tau.ac.il>
#ifndef CGAL_CURVE_DATA_TRAITS_2_H
#define CGAL_CURVE_DATA_TRAITS_2_H

#include <list>

CGAL_BEGIN_NAMESPACE

/*! A generic traits class for maintaining an arrangement of curves that have
 * an extra data field. This traits class is templated with a Data class an
 * an ordinary traits class which is also used as a based traits class to
 * inherit from. It extracts the original Curve_2 and X_monotone_curve_2 types
 * from the ordinary traits class, and redefines them to have Data as an extra
 * field.
 *
 * The Data field is updated when the curves are converted from Curve_2 to
 * X_monotone_curve_2, and when the X_monotone_curve_2 curves are split. All
 * other methods are inherited from the base ordinary traits class.
 */

template <class Traits_, class Data_>
class Arr_curve_data_traits_2 : public Traits_ 
{
public:
  typedef Traits_                               Traits;
  typedef Data_                                 Data;
  typedef typename Traits::Curve_2              Org_curve_2;
  typedef typename Traits::X_monotone_curve_2   Org_x_monotone_curve_2;
  typedef typename Traits::Point_2              Point_2;

  /*!
   * Representation of an input curve with an addtional data field.
   */
  class Curve_2 : public Org_curve_2 
  {
  private:
    Data m_data;

  public:

    /*!
     * Default constructor.
     */
    Curve_2 ()
    {}
    
    /*!
     * Construct a curve from an original curve and a data object.
     * \param cv The original curve.
     * \param data The data object.
     */ 
    Curve_2 (const Org_curve_2 & cv, const Data & data) :
      Org_curve_2(cv),
      m_data(data)
    {}

    /*!
     * Get the data.
     * \return The data object associated with the curve.
     */
    const Data & get_data () const
    {
      return m_data;
    }

    /*!
     * Set the curve data.
     * \param data The data object to be associated with the curve.
     */
    void set_data (const Data & data)
    {
      m_data = data;
      return;
    }
  };
  
  /*!
   * Representation of an x-monotone cuvre. As this curve may represent
   * an overlapping section of several input curves, we store a list of data
   * objects with it.
   */
  class X_monotone_curve_2 : public Org_x_monotone_curve_2 
  {
  private:

    std::list<Data>  m_data_list;

  public:

    /*!
     * Default constructor.
     */
    X_monotone_curve_2()
    {}
    
    /*!
     * Construct a curve from an original x-monotone curve and a data object.
     * \param cv The original x-monotone curve.
     * \param data The data object.
     */ 
    X_monotone_curve_2 (const Org_x_monotone_curve_2 & cv, const Data & data) :
      Org_x_monotone_curve_2(cv),
      m_data_list()
    {
      m_data_list.push_back(data);
    }

    /*!
     * Construct a curve from an original x-monotone curve and a range of 
     * data objects.
     * \param cv The original x-monotone curve.
     * \param begin A begin iterator for the data range.
     * \param end A past-the-end iterator for the data range.
     */
    template <class InputIterator>
    X_monotone_curve_2 (const Org_x_monotone_curve_2 & cv, 
			const InputIterator& begin, const InputIterator& end) :
      Org_x_monotone_curve_2(cv),
      m_data_list()
    {
      InputIterator    it;

      for (it = begin; it != end; it++)
	m_data_list.push_back(*it);
    }

    /*!
     * Get the number of data objects associated with the x-monotne curve.
     */
    int number_of_data_objects () const
    {
      return (m_data_list.size());
    }

    /*!
     * Get the first data object associated with the curve.
     * \pre number_of_data_objects() is not 0.
     */
    const Data& get_data () const
    {
      CGAL_precondition (m_data_list.size() > 0);

      return (m_data_list.front());
    }

    /*!
     * Get the data iterators (const version).
     */
    typename std::list<Data>::const_iterator begin_data () const
    {
      return (m_data_list.begin());
    }

    typename std::list<Data>::const_iterator end_data () const
    {
      return (m_data_list.end());
    }

    /*!
     * Get the data iterators (non-const version).
     */
    typename std::list<Data>::iterator begin_data ()
    {
      return (m_data_list.begin());
    }

    typename std::list<Data>::iterator end_data ()
    {
      return (m_data_list.end());
    }

    /*!
     * Add a data object to the curve.
     * \param data The additional data object.
     */
    void add_data (const Data & data)
    {
      m_data_list.push_back (data);
      return;
    }

    /*!
     * Set a data object to the curve.
     * \param data The data object to set.
     */
    void set_data (const Data & data)
    {
      clear_data();
      add_data(data);
      return;
    }

    /*!
     * Add a range of data objects to the curve.
     * \param begin A begin iterator for the data range.
     * \param end A past-the-end iterator for the data range.
     */
    template <class InputIterator>
    void add_data (const InputIterator & begin, const InputIterator & end)
    {
      InputIterator    it;

      for (it = begin; it != end; it++)
	m_data_list.push_back(*it);

      return;
    }

    /*!
     * Clear the data objects.
     */
    void clear_data ()
    {
      m_data_list.clear();
      return;
    }
  };

  // For backward compatibility:
  typedef Point_2                               Point;
  typedef X_monotone_curve_2                    X_curve;
  typedef Curve_2                               Curve;
  
public:

  /*!
   * Default constructor.
   */
  Arr_curve_data_traits_2 () : 
    Traits() {}
  
  /*!
   * Cut the given curve into x-monotone subcurves and insert them to the
   * given output iterator. While segments are x_monotone, still need to pass
   * them out.
   * \param cv The curve.
   * \param o The output iterator.
   * \return The past-the-end iterator.
   */
  template<class OutputIterator>
  OutputIterator curve_make_x_monotone(const Curve_2 & cv,
                                       OutputIterator o) const
  {
    // Make the original curve x-monotone.
    std::list<Org_x_monotone_curve_2>  org_x_curves;
    
    Traits::curve_make_x_monotone (cv, std::back_inserter(org_x_curves));

    // Attach the data to each of the resulting x-monotone curves.
    typename std::list<Org_x_monotone_curve_2>::const_iterator it;

    for (it = org_x_curves.begin(); it != org_x_curves.end(); it++) 
      *o++ = X_monotone_curve_2 (*it, cv.get_data());
    
    return o;
  } 

  /*! 
   * Split a given curve at a given split point into two sub-curves.
   * \param cv The curve to split
   * \param c1 The output first part of the split curve. 
   *           Its source is the source of the original curve.
   * \param c2 The output second part of the split curve.
   *           Its target is the target of the original curve.
   * \param p The split point.
   * \pre p lies on cv but is not one of its end-points.
   */
  void curve_split(const X_monotone_curve_2 & cv, 
		   X_monotone_curve_2 & c1, X_monotone_curve_2 & c2, 
                   const Point_2 & p) const
  {
    // Split the original curve.
    Org_x_monotone_curve_2 org_c1, org_c2;
    Traits::curve_split (cv, org_c1, org_c2, p);

    // Attach data.
    c1 = X_monotone_curve_2 (org_c1, cv.begin_data(), cv.end_data());
    c2 = X_monotone_curve_2 (org_c2, cv.begin_data(), cv.end_data());

    return;
  }

  /*!
   * Find the nearest intersection of the two given curves to the right of 
   * a given reference point.
   * Nearest is defined as the lexicographically nearest point, not including 
   * the point reference point itself.
   * If the intersection of the two curves is an X_monotone_curve_2, that is,
   * there is an overlapping subcurve, that contains the reference point in
   * its x-range, the function should return an X_monotone_curve_2 whose 
   * interior is strictly to the right of the reference point (that is, whose
   * left endpoint is the projection of the reference point onto the 
   * overlapping subcurve).
   * \param cv1 The first curve.
   * \param cv2 The second curve.
   * \param p The refernece point.
   * \return An empty object if there is no intersection to the right of p.
   *         An object wrapping a Point_2 in case of a simple intersection.
   *         An object wrapping an X_monotone_curve_2 in case of an overlap.
   */
  Object nearest_intersection_to_right (const X_monotone_curve_2 & cv1,
				        const X_monotone_curve_2 & cv2,
					const Point_2 & p) const
  {
    // Use the base method and resturn the result as is in case there is
    // no intersection or if the intersection is a single point.
    CGAL::Object res = Traits::nearest_intersection_to_right (cv1, cv2, p);
    Point_2      ip;

    if (res.is_empty() || CGAL::assign (ip, res))
      return (res);

    // In case we have an overlapping curve, attach data from both cv1 and cv2
    // to it.
    Org_x_monotone_curve_2  iocv;
    bool                    assign_success;
    
    assign_success = CGAL::assign (iocv, res);
    CGAL_assertion (assign_success);
    
    X_monotone_curve_2      icv (iocv, cv1.begin_data(), cv1.end_data());

    icv.add_data (cv2.begin_data(), cv2.end_data());
    return (CGAL::make_object (icv));
  }

  /*!
   * Find the nearest intersection of the two given curves to the left of 
   * a given reference point.
   * Nearest is defined as the lexicographically nearest point, not including 
   * the point reference point itself.
   * If the intersection of the two curves is an X_monotone_curve_2, that is,
   * there is an overlapping subcurve, that contains the reference point in
   * its x-range, the function should return an X_monotone_curve_2 whose 
   * interior is strictly to the left of the reference point (that is, whose
   * right endpoint is the projection of the reference point onto the 
   * overlapping subcurve).
   * \param cv1 The first curve.
   * \param cv2 The second curve.
   * \param p The refernece point.
   * \return An empty object if there is no intersection to the left of p.
   *         An object wrapping a Point_2 in case of a simple intersection.
   *         An object wrapping an X_monotone_curve_2 in case of an overlap.
   */
  Object nearest_intersection_to_left (const X_monotone_curve_2 & cv1,
				       const X_monotone_curve_2 & cv2,
				       const Point_2 & p) const
  {
    // Use the base method and resturn the result as is in case there is
    // no intersection or if the intersection is a single point.
    CGAL::Object res = Traits::nearest_intersection_to_left (cv1, cv2, p);
    Point_2      ip;

    if (res.is_empty() || CGAL::assign (ip, res))
      return (res);

    // In case we have an overlapping curve, attach data from both cv1 and cv2
    // to it.
    Org_x_monotone_curve_2  iocv;
    bool                    assign_success;
    
    assign_success = CGAL::assign (iocv, res);
    CGAL_assertion (assign_success);
    
    X_monotone_curve_2      icv (iocv, cv1.begin_data(), cv1.end_data());

    icv.add_data (cv2.begin_data(), cv2.end_data());
    return (CGAL::make_object (icv));
  }
};

CGAL_END_NAMESPACE

#endif
