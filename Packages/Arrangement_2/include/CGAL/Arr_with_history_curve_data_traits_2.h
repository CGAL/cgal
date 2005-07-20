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
// Author(s)     : Baruch Zukerman   <baruchzu@post.tau.ac.il>
//                 Efi Fogel         <efif@post.tau.ac.il>
#ifndef CGAL_ARR_WITH_HISTORY_CURVE_DATA_TRAITS_2_H
#define CGAL_ARR_WITH_HISTORY_CURVE_DATA_TRAITS_2_H

#include <CGAL/Arr_consolidated_curve_data_traits_2.h>

#include <list>

template <class T_Traits>
class Arr_with_history_curve_data_traits_2 :
  public Arr_consolidated_curve_data_traits_2<T_Traits,
                                              typename T_Traits::Curve_2 *>
{
private:
  typedef T_Traits                              Traits;
  typedef typename Traits::Curve_2              Base_curve_2;
  typedef Arr_consolidated_curve_data_traits_2<T_Traits,Base_curve_2 *>
    Base;

  std::list<Base_curve_2>                       m_base_curves;

public:
  class Make_x_monotone_2 {
  private:
    Base_traits * base;
    
  public:

    /*! Constructor. */
    Make_x_monotone_2(Base_traits *_base) : base (_base) {}
    
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
    OutputIterator operator()(const Curve_2 & cv, OutputIterator oi)
    {
      m_base_curves.push_back(static_cast<Base_curve_2>(cv));
      cv.set_data(&(m_base_curves.back()));
      oi = base->make_x_monotone_2_object()(cv, oi);
      return oi;
    }
  };

  /*! Get a Make_x_monotone_2 functor object. */
  Make_x_monotone_2 make_x_monotone_2_object ()
  {
    return Make_x_monotone_2 (this);
  }  
};

#endif
