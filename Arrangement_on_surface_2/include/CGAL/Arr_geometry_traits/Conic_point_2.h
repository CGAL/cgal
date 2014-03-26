// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
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
// Author(s)     : Ron Wein <wein@post.tau.ac.il>

#ifndef CGAL_CONIC_POINT_2_H
#define CGAL_CONIC_POINT_2_H

/*! \file
 * Header file for the _Conic_point_2<Alg_kernel> class.
 */

#include <list>
#include <CGAL/assertions.h>

namespace CGAL {

/*!
 * \class A class that stores additional information with the point's 
 * coordinates, namely the conic IDs of the generating curves.
 */
template <class Alg_kernel_>
class _Conic_point_2 : public Alg_kernel_::Point_2
{
public:

  typedef Alg_kernel_                       Alg_kernel;
  typedef typename Alg_kernel::Point_2      Base;
  typedef _Conic_point_2<Alg_kernel>        Self;
    
  typedef typename Alg_kernel::FT           Algebraic;

  /*! \class
   * Representation of an ID of a conic arc.
   */
  class Conic_id
  {
  private:

    unsigned int   index;           // The index of the conic arc.

  public:

    /*! Default constructor. */
    Conic_id () :
      index (0)
    {}

    /*! Constructor. */
    Conic_id (unsigned int ind) :
      index (ind)
    {
      CGAL_precondition (ind != 0);
    }

    /*! Check if the ID is valid. */
    bool is_valid () const
    {
      return (index != 0);
    }
    
    /*! Equality operator. */
    bool operator== (const Conic_id& id) const
    {
      return (index == id.index);
    }

    /*! Inequality operator. */
    bool operator!= (const Conic_id& id) const
    {
      return (index != id.index);
    }

    /*! Less-than operator. */
    bool operator< (const Conic_id& id) const
    {
      return (index < id.index);
    }

    /*! Greater-than operator. */
    bool operator> (const Conic_id& id) const
    {
      return (index > id.index);
    }
  };
        
private:

  typedef std::list<Conic_id>                          Ids_container;
  typedef typename std::list<Conic_id>::const_iterator Ids_iterator;

  Ids_container   conic_ids;       // The IDs of the generating conics.

 public:

  /// \name Constructors.
  //@{

  /*! Default constructors. */
  _Conic_point_2 () :
    Base()
  {}

  /*! Constrcutor from the base class. */
  _Conic_point_2 (const Base& p) :
    Base (p)
  {}

  /*! Constructor with homegeneous coordinates. */
  _Conic_point_2 (const Algebraic& hx, 
		  const Algebraic& hy,
		  const Algebraic& hz) :
    Base (hx, hy, hz)
  {}

  /*! Constructor with Cartesian coordinates. */
  _Conic_point_2 (const Algebraic& x, const Algebraic& y) :
    Base (x, y)
  {}
  //@}

  /// \name Maintaining the generating conic IDs.
  //@{

  /*! Add a generating conic ID. */
  void set_generating_conic (const Conic_id& id)
  {
    if (id.is_valid())
      conic_ids.push_back (id);

    return;
  }

  /*! Check if the given conic generates the point. */
  bool is_generating_conic (const Conic_id& id) const
  {
    if (! id.is_valid())
      return (false);

    Ids_iterator       it;

    for (it = conic_ids.begin(); it != conic_ids.end(); ++it)
    {
      if (*it == id)
        return (true);
    }

    return (false);
  }
  //@}

};

} //namespace CGAL

#endif
