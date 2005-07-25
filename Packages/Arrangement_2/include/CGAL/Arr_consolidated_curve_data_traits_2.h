// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Ron Wein          <wein@post.tau.ac.il>
//                 Efi Fogel         <efif@post.tau.ac.il>
#ifndef CGAL_ARR_CONSOLIDATED_CURVE_DATA_TRAITS_2_H
#define CGAL_ARR_CONSOLIDATED_CURVE_DATA_TRAITS_2_H

/*! \file
 * Definition of the Arr_consolidated_curve_data_traits_2<Traits,Data> class.
 */

#include <list>

CGAL_BEGIN_NAMESPACE

/*! \class
 * A generic traits class for maintaining an arrangement of curves that have
 * an extra data field. This traits class is templated with a Data class an
 * an ordinary traits class which is also used as a based traits class to
 * inherit from. It extracts the original Curve_2 and X_monotone_curve_2 types
 * from the ordinary traits class, and redefines them to have Data as an extra
 * field in the Curve_2 type, and a container of Data objects for the extended
 * X_monotone_curve_2 type.
 * The Data field is updated when the curves are converted from Curve_2 to
 * X_monotone_curve_2, and when the X_monotone_curve_2 curves are split.
 * When two x-monotone curves overlap, their data containers are consolidated
 * and attached to the resulting subcurve.
 * All other functors are inherited from the base ordinary traits class.
 */
template <class Traits_, class Data_>
class Arr_consolidated_curve_data_traits_2 : public Traits_ 
{
public:

  typedef Traits_                                   Base_traits;
  typedef Data_                                     Data;
  typedef typename Base_traits::Curve_2             Base_curve_2;
  typedef typename Base_traits::X_monotone_curve_2  Base_x_monotone_curve_2;
  typedef typename Base_traits::Point_2             Point_2;

  typedef typename Base_traits::Has_left_category   Has_left_category;
  typedef typename Base_traits::Has_merge_category  Base_has_merge_category;
  typedef Tag_true                                  Has_merge_category;

  /*!
   * Representation of an input curve with an addtional data field.
   */
  class Curve_2 : public Base_curve_2 
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
    Curve_2 (const Base_curve_2& cv, const Data& data) :
      Base_curve_2(cv),
      m_data(data)
    {}

    /*!
     * Get the data.
     * \return The data object associated with the curve.
     */
    const Data& get_data () const
    {
      return m_data;
    }

    /*!
     * Set the curve data.
     * \param data The data object to be associated with the curve.
     */
    void set_data (const Data& data)
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
  class X_monotone_curve_2 : public Base_x_monotone_curve_2 
  {
  private:
    typedef std::list<Data>                          Data_container;

  public:
    typedef typename Data_container::iterator        Data_iterator;
    typedef typename Data_container::const_iterator  Data_const_iterator;

  private:

    Data_container  m_data_list;

  public:

    /*!
     * Default constructor.
     */
    X_monotone_curve_2()
    {}

    /*!
     * Construct a curve from an original x-monotone with no data object. 
     * \param cv The original x-monotone curve.
     */
    X_monotone_curve_2 (const Base_x_monotone_curve_2& cv) :
      Base_x_monotone_curve_2(cv),
      m_data_list()
    {}

    /*!
     * Construct a curve from an original x-monotone curve and a data object.
     * \param cv The original x-monotone curve.
     * \param data The data object.
     */ 
    X_monotone_curve_2 (const Base_x_monotone_curve_2& cv,
                        const Data& data) :
      Base_x_monotone_curve_2(cv),
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
    X_monotone_curve_2 (const Base_x_monotone_curve_2& cv, 
                        const InputIterator& begin,
                        const InputIterator& end) :
      Base_x_monotone_curve_2(cv),
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
    Data_const_iterator data_begin () const
    {
      return (m_data_list.begin());
    }

    Data_const_iterator data_end () const
    {
      return (m_data_list.end());
    }

    /*!
     * Get the data iterators (non-const version).
     */
    Data_iterator data_begin ()
    {
      return (m_data_list.begin());
    }

    Data_iterator data_end ()
    {
      return (m_data_list.end());
    }

    /*!
     * Check whether the given curve has the same data.
     */
    bool has_same_data (const X_monotone_curve_2& cv) const
    {
      return (m_data_list == cv.m_data_list);
    }

    /*!
     * Add a data object to the curve.
     * \param data The additional data object.
     */
    void add_data (const Data& data)
    {
      m_data_list.push_back (data);
      return;
    }

    /*!
     * Remove a data object from the curve.
     * \param data The data object to be removed.
     */
    void remove_data (const Data& data)
    {
      m_data_list.remove (data);
      return;
    }

    /*!
     * Set a data object to the curve.
     * \param data The data object to set.
     */
    void set_data (const Data& data)
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
    void add_data (const InputIterator& begin, const InputIterator& end)
    {
      InputIterator    it;

      for (it = begin; it != end; ++it)
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



public:
  
  /// \name Construction.
  //@{

  /*! Default constructor. */
  Arr_consolidated_curve_data_traits_2 ()
  {}
  
  /*! Constructor from a base-traits class. */
  Arr_consolidated_curve_data_traits_2 (const Base_traits& traits) :
    Base_traits (traits)
  {}
  //@}

  /// \name Overriden functors.
  //@{

  class Make_x_monotone_2
  {
  private:
    Base_traits    *base;

  public:

    /*! Constructor. */
    Make_x_monotone_2 (Base_traits *_base) :
      base (_base)
    {}
    
    /*!
     * Cut the given curve into x-monotone subcurves and insert them to the
     * given output iterator. As segments are always x_monotone, only one
     * x-monotone curve will be contained in the iterator.
     * \param cv The curve.
     * \param oi The output iterator, whose value-type is Object. The output
     *           may be either X_monotone_curve_2 objects or Point_2 objects
     *           (in case the input curve contains isolated points).
     * \return The past-the-end iterator.
     */
    template<class OutputIterator>
    OutputIterator operator() (const Curve_2& cv, OutputIterator oi) const
    {
      // Make the original curve x-monotone.
      std::list<Object>   base_objects;
    
      base->make_x_monotone_2_object() (cv,
                                        std::back_inserter (base_objects));

      // Go over the returned objects and attach the data to each of the
      // resulting x-monotone curves.
      typename std::list<Object>::const_iterator  it;
      const Base_x_monotone_curve_2              *base_x_curve;

      for (it = base_objects.begin(); it != base_objects.end(); ++it)
      { 
        base_x_curve = object_cast<Base_x_monotone_curve_2> (&(*it));
        if (base_x_curve != NULL)
        {
          // The current object is a base x-monotone curve: Attach data to it.
          *oi = make_object (X_monotone_curve_2 (*base_x_curve,
                                                 cv.get_data()));
        }
        else
        {
          CGAL_assertion (object_cast<Point_2> (&(*it)) != NULL);

          // The current object is an isolated point: Leave it as is.
          *oi = *it;
        }
        ++oi;
      }

      return (oi);
    }
  };

  /*! Get a Make_x_monotone_2 functor object. */
  Make_x_monotone_2 make_x_monotone_2_object ()
  {
    return Make_x_monotone_2 (this);
  }

  class Split_2
  {
  private:
    Base_traits    *base;

  public:

    /*! Constructor. */
    Split_2 (Base_traits *_base) :
      base (_base)
    {}

    /*!
     * Split a given x-monotone curve at a given point into two sub-curves.
     * \param cv The curve to split
     * \param p The split point.
     * \param c1 Output: The left resulting subcurve (p is its right endpoint).
     * \param c2 Output: The right resulting subcurve (p is its left endpoint).
     * \pre p lies on cv but is not one of its end-points.
     */
    void operator() (const X_monotone_curve_2& cv, const Point_2 & p,
                     X_monotone_curve_2& c1, X_monotone_curve_2& c2) const
    {
      // Split the original curve.
      base->split_2_object() (cv, p,
                              c1, c2);

      // Attach data to the split curves.
      c1.clear_data();
      c1.add_data (cv.data_begin(), cv.data_end());

      c2.clear_data();
      c2.add_data (cv.data_begin(), cv.data_end());

      return;
    }
  };

  /*! Get a Split_2 functor object. */
  Split_2 split_2_object ()
  {
    return Split_2 (this);

  }

  class Intersect_2
  {
  private:
    Base_traits    *base;

  public:

    /*! Constructor. */
    Intersect_2 (Base_traits *_base) :
      base (_base)
    {}

    /*!
     * Find the intersections of the two given curves and insert them to the
     * given output iterator. As two segments may itersect only once, only a
     * single will be contained in the iterator.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param oi The output iterator.
     * \return The past-the-end iterator.
     */
    template<class OutputIterator>
    OutputIterator operator() (const X_monotone_curve_2& cv1,
                               const X_monotone_curve_2& cv2,
                               OutputIterator oi) const
    {
      // Use the base functor to obtain all intersection objects.
      std::list<Object>                      base_list;

      base->intersect_2_object() (cv1, cv2,
                                  std::back_inserter (base_list));

      // Stop if the list is empty:
      if (base_list.empty())
        return (oi);

      // Go over all intersection objects and prepare the output.
      typename std::list<Object>::iterator   curr;
      const Base_x_monotone_curve_2         *base_cv;

      for (curr = base_list.begin(); curr != base_list.end(); ++curr)
      {
        if ((base_cv = object_cast<Base_x_monotone_curve_2>(&(*curr))) != NULL)
        {
          // The current intersection object is an overlapping x-monotone
          // curve: First attach data from the first curve.
          X_monotone_curve_2  cv (*base_cv,
                                  cv1.data_begin(), cv1.data_end());

          // Add data from the second curve.
          cv.add_data (cv2.data_begin(), cv2.data_end());

          // Output the curve:
          *oi = make_object (cv);
        }
        else
        {
          // The current intersection object is an intersection point:
          // Copy it as is.
          *oi = *curr;
        }

        ++oi;
      }

      return (oi);
    }
  };

  /*! Get an Intersect_2 functor object. */
  Intersect_2 intersect_2_object () 
  {
    return Intersect_2 (this);
  }

  class Are_mergeable_2
  {
  private:
    const Base_traits    *base;

  public:

    /*! Constructor. */
    Are_mergeable_2 (const Base_traits *_base) :
      base (_base)
    {}

    /*!
     * Check whether it is possible to merge two given x-monotone curves.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \return (true) if the two curves are mergeable; (false) otherwise.
     */
    bool operator() (const X_monotone_curve_2& cv1,
                     const X_monotone_curve_2& cv2) const
    {
      // In case the two base curves are not mergeable, the extended curves
      // are not mergeable as well.
      if (!_are_mergeable_base_imp (cv1, cv2,
                                    Base_has_merge_category()))
      {
        return (false);
      }

      // Make sure that the data attached to both curves is the same.
      return (cv1.has_same_data (cv2));
    }

  private:

    /*!
     * Implementation of the base predicate in case the HasMerge tag is true.
     */
    bool _are_mergeable_base_imp (const X_monotone_curve_2& cv1,
                                  const X_monotone_curve_2& cv2,
                                  Tag_true) const
    {
      return (base->are_mergeable_2_object() (cv1, cv2));      
    }

    /*!
     * Implementation of the base predicate in case the HasMerge tag is false.
     */
    bool _are_mergeable_base_imp (const X_monotone_curve_2& ,
                                  const X_monotone_curve_2& ,
                                  Tag_false) const
    {
      // Curve merging is not supported:
      return (false);
    }
  };
  
  /*! Get an Are_mergeable_2 functor object. */
  Are_mergeable_2 are_mergeable_2_object () const
  {
    return Are_mergeable_2 (this);
  }

  class Merge_2
  {
  private:
    Base_traits    *base;

  public:

    /*! Constructor. */
    Merge_2 (Base_traits *_base) :
      base (_base)
    {}

    /*!
     * Merge two given x-monotone curves into a single curve (segment).
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param c Output: The merged curve.
     * \pre The two curves are mergeable.
     */
    void operator() (const X_monotone_curve_2& cv1,
                     const X_monotone_curve_2& cv2,
                     X_monotone_curve_2& c) const
    {
      // The function is implemented based on the base Has_merge category.
      _merge_imp (cv1, cv2, c, Base_has_merge_category());
    }

  private:

    /*!
     * Implementation of the operator() in case the HasMerge tag is true.
     */
    void _merge_imp (const X_monotone_curve_2& cv1,
                     const X_monotone_curve_2& cv2,
                     X_monotone_curve_2& c,
                     Tag_true) const
    {      
      // Merge the two base curve.
      Base_x_monotone_curve_2  base_cv;

      base->merge_2_object() (cv1, cv2,
                              base_cv);

      // Attach data from one of the curves.
      CGAL_precondition (cv1.has_same_data (cv2));

      c = X_monotone_curve_2 (base_cv,
                              cv1.data_begin(), cv1.data_end());
      return;
    }

    /*!
     * Implementation of the operator() in case the HasMerge tag is false.
     */
    void _merge_imp (const X_monotone_curve_2& ,
                     const X_monotone_curve_2& ,
                     X_monotone_curve_2& ,
                     Tag_false) const
    {
      // This function should never be called!
      CGAL_assertion_msg (false,
                          "Merging curves is not supported.");
    }
  };

  /*! Get a Merge_2 functor object. */
  Merge_2 merge_2_object ()
  {
    return Merge_2 (this);
  }

  class Construct_x_monotone_curve_2
  {
  private:
    const Base_traits    *base;

  public:

    /*! Constructor. */
    Construct_x_monotone_curve_2 (const Base_traits *_base) :
      base (_base)
    {}

    /*!
     * Return an x-monotone curve connecting the two given endpoints.
     * \param p The first point.
     * \param q The second point.
     * \pre p and q must not be the same.
     * \return An x-monotone curve connecting p and q.
     */
    X_monotone_curve_2 operator() (const Point_2& p,
                                   const Point_2& q) const
    {
      Base_x_monotone_curve_2  base_cv =
        base->construct_x_monotone_curve_2_object() (p, q);

      return (X_monotone_curve_2 (base_cv));
    }
  };

  /*! Get a Construct_x_monotone_curve_2 functor object. */
  Construct_x_monotone_curve_2 construct_x_monotone_curve_2_object () const
  {
    return Construct_x_monotone_curve_2 (this);
  }
  //@}

};

CGAL_END_NAMESPACE

#endif

