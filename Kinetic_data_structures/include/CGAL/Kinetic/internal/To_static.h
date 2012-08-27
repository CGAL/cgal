// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_KINETIC_TO_STATIC_H
#define CGAL_KINETIC_TO_STATIC_H
#include <CGAL/Kinetic/basic.h>
#include <CGAL/Simple_cartesian.h>

namespace CGAL { namespace Kinetic {

//! a functor that turns a moving object into a static object. It needs to pick the static type.
/*!  This is the default implementation which looks for a class named
  Static_traits in the moving object type (the first template
  argument) and uses the Static_type field of that traits when
  instantiated with the static kernel to figure out the return
  type. It then calls the to_static member method of the class.
*/
template <class Arg,
	  class SK >
class To_static
{
  typedef typename Arg::template Static_traits<SK> Traits;
public:
  //! The way time is represented.
  typedef typename SK::FT Time;

  //! Construct it with a static kernel
  To_static(const SK &sk=SK()): sk_(sk), t_(-666666666) {
#ifndef NDEBUG
    initialized_=false;
#endif
  }

  typedef typename Traits::Static_type result_type;
  typedef Arg argument_type;
  //! Convert an appropriate moving object to a static object
  result_type operator()(const argument_type &arg) const
  {
#ifndef NDEBUG
    if (!initialized_) {
      std::cerr << "You must set time before using a snapshot.\n";
      CGAL_precondition(initialized_);
    }
#endif
    return Traits::to_static(arg, time(), sk_);
  }
  //! What this believes the time to be.
  const Time& time() const
  {
    /*#ifndef NDEBUG
      if (!initialized_) {
      std::cerr << "You must set time before using a snapshot.\n";
      CGAL_precondition(initialized_);
      }
      #endif    */
    return t_;
  }
  //! Set the time.
  void set_time(const Time &t) {
    t_=t;
#ifndef NDEBUG
    initialized_=true;
#endif
  }
protected:
  SK sk_;
  Time t_;
#ifndef NDEBUG
  bool initialized_;
#endif
};
} } //namespace CGAL::Kinetic
#endif
