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

#ifndef CGAL_KINETIC_INSERT_EVENT_H_
#define CGAL_KINETIC_INSERT_EVENT_H_
#include <CGAL/Kinetic/basic.h>
#include <CGAL/Kinetic/internal/To_static.h>
#include <CGAL/Kinetic/Event_base.h>

namespace CGAL { namespace Kinetic {

//! An event to insert a single object into a MovingObjectTable

template <class MOT>
class Insert_event: public Event_base<int*>
{
  typedef typename MOT::Handle Pointer;
  typedef typename MOT::Data Object;
public:
  //! Construct it with the time, the object and where to put it
  Insert_event(const Object &obj,
	       Pointer mot):mot_(mot),
			    obj_(obj){}
  void process() {
    CGAL_LOG(Log::SOME, "Inserting object.\n");
    mot_->insert(obj_);
  }
  void* kds() const {return NULL;}

  std::ostream& write(std::ostream &out) const
  {
    out << " I" << obj_ ;
    return out;
  }
protected:
  Pointer mot_;
  Object obj_;
};

template <class MH>
std::ostream &operator<<(std::ostream &out, const Insert_event< MH> &moi)
{
  moi.write(out);
  return out;
}


} } //namespace CGAL::Kinetic
#endif
