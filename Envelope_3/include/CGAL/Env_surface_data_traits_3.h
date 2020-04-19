// Copyright (c) 2006  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Ron Wein          <wein@post.tau.ac.il>
#ifndef CGAL_ENV_SURFACE_DATA_TRAITS_3_H
#define CGAL_ENV_SURFACE_DATA_TRAITS_3_H

#include <CGAL/license/Envelope_3.h>


/*! \file
 * Definition of the env_surface_data_traits_3<> class template.
 */

#include <CGAL/Arr_geometry_traits/Curve_data_aux.h>
#include <list>

namespace CGAL {

/*! \class
 * A generic traits class for maintaining the envelope of surfaces that have
 * an extra data field. This traits class is templated with an ordinary traits
 * class, which is also used as a based traits class to inherit from.
 * It can attach data objects to Surface_3 and to Xy_monotone_surface_3 objects
 * (possibly of two different types).
 */
template <class Traits_, class XyMonotoneSurfaceData_,
          class SurfaceData_ = XyMonotoneSurfaceData_,
          class Convert_ = _Default_convert_func<SurfaceData_,
                                                 XyMonotoneSurfaceData_> >
class Env_surface_data_traits_3 : public Traits_
{
public:

  typedef Traits_                                  Base_traits_3;
  typedef XyMonotoneSurfaceData_                   Xy_monotone_surface_data;
  typedef SurfaceData_                             Surface_data;
  typedef Convert_                                 Convert;

  typedef typename Base_traits_3::Surface_3        Base_surface_3;
  typedef typename Base_traits_3::Xy_monotone_surface_3
                                                   Base_xy_monotone_surface_3;

  // Representation of a surface with an addtional data field:
  typedef _Curve_data_ex<Base_surface_3, Surface_data>  Surface_3;

  // Representation of an xy-monotone surface with an addtional data field:
  typedef _Curve_data_ex<Base_xy_monotone_surface_3,
                         Xy_monotone_surface_data>      Xy_monotone_surface_3;

public:

  /// \name Construction.
  //@{

  /*! Default constructor. */
  Env_surface_data_traits_3 ()
  {}

  /*! Constructor from a base-traits class. */
  Env_surface_data_traits_3 (const Base_traits_3 & traits) :
    Base_traits_3 (traits)
  {}
  //@}

  /// \name Overriden functors.
  //@{

  class Make_xy_monotone_3
  {
  private:
    const Base_traits_3 * base;

  public:

    /*! Constructor. */
    Make_xy_monotone_3 (const Base_traits_3 * _base) : base (_base)
    {}

    /*!
     * Subdivide the given surface into xy-monotone surfaces and insert them
     * to the given output iterator.
     * \param S The surface.
     * \param oi The output iterator,
     *           whose value-type is Xy_monotone_surface_2.
     * \return The past-the-end iterator.
     */
    template<class OutputIterator>
    OutputIterator operator() (const Surface_3& S, bool is_lower,
                               OutputIterator oi) const
    {
      // Make the original surface xy-monotone.
      std::list<Base_xy_monotone_surface_3>                     xy_surfaces;
      typename std::list<Base_xy_monotone_surface_3>::iterator  xys_it;

      base->make_xy_monotone_3_object()
        (S, is_lower, std::back_inserter (xy_surfaces));

      // Attach the data to each of the resulting xy-monotone surfaces.
      for (xys_it = xy_surfaces.begin(); xys_it != xy_surfaces.end(); ++xys_it)
      {
        *oi = Xy_monotone_surface_3 (*xys_it,
                                     Convert() (S.data()));
        ++oi;
      }

      return (oi);
    }
  };

  /*! Get a Make_xy_monotone_3 functor object. */
  Make_xy_monotone_3 make_xy_monotone_3_object() const
  {
    return Make_xy_monotone_3 (this);
  }

  //@}

};

} //namespace CGAL

#endif
