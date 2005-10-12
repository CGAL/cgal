// Copyright (c) 2005  Stanford University (USA).
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
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_KDS_INSERT_EVENT_H_
#define CGAL_KDS_INSERT_EVENT_H_
#include <CGAL/KDS/basic.h>
#include <CGAL/KDS/internal/To_static.h>

CGAL_KDS_BEGIN_NAMESPACE;


//! An event to insert a single object into a MovingObjectTable

template <class MOT>
class Insert_event {
  typedef typename MOT::Pointer Pointer;
  typedef typename MOT::Data Object;
public:
  //! Construct it with the time, the object and where to put it
  Insert_event(const Object &obj, 
	       Pointer mot):mot_(mot),
			    obj_(obj){}
  template <class T>
  void process(const T&) {
    CGAL_KDS_LOG(LOG_SOME, "Inserting object.\n");
    mot_->insert(obj_);
  }

   void write(std::ostream &out) const {
     out << " I" << obj_ ;
  }
protected:
  Pointer mot_;
  Object obj_;
};

template <class MH>
std::ostream &operator<<(std::ostream &out, const Insert_event< MH> &moi){
  moi.write(out);
  return out;
}



CGAL_KDS_END_NAMESPACE;
#endif
