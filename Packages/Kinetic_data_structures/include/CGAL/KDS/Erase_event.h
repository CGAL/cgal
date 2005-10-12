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

#ifndef CGAL_KDS_MOVING_OBJECT_ERASER_H_
#define CGAL_KDS_MOVING_OBJECT_ERASER_H_
#include <CGAL/KDS/basic.h>

CGAL_KDS_BEGIN_NAMESPACE;


//! Delete a single moving object from the MOT at a particular time. 
/*!
  Note that this class has not been used. 
*/
template <class MOT>
class Erase_event{
  typedef typename MOT::Pointer Pointer;
  typedef typename MOT::Key Key;
public:
  Erase_event(Key k, 
			Pointer mot):mot_(mot),
				     k_(k){}
  template <class T>
  void process(const T&) {
    CGAL_KDS_LOG(LOG_SOME,"Deleting object.\n");
    mot_->erase(k_);
  }
  void write(std::ostream &out) const {
    out << "E" << k_;
  }
protected:
  Pointer mot_;
  Key k_;
};

template <class MH>
std::ostream &operator<<(std::ostream &out, const Erase_event<MH> &moi){
  moi.write(out);
  return out;
}



CGAL_KDS_END_NAMESPACE;
#endif
