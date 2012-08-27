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

#ifndef CGAL_KDS_EVENT_BASE_H
#define CGAL_KDS_EVENT_BASE_H
#include <CGAL/Kinetic/basic.h>
namespace CGAL { namespace Kinetic {

template <class KDS_ptr>
class Event_base {
public:
  typedef KDS_ptr KDS_handle;
  Event_base(KDS_ptr pt): kds_(pt){}
  Event_base(): kds_(NULL){}
	
  KDS_handle kds() const {
    return kds_;
  }
  KDS_handle kds() {
    return kds_;
  }
  template <class Key>
  CGAL::Comparison_result compare_concurrent(Key a,
					     Key b) const {
    if (a < b) return CGAL::SMALLER;
    else if (b < a) return CGAL::LARGER;
    else return CGAL::EQUAL;
    //return CGAL::compare(a,b);
  }	
  template <class Key>
  bool merge_concurrent(Key ,
			Key ) {
    return false;
  }
  template <class K>
  void audit(K ) const {}
  std::ostream &write(std::ostream &out) const {
    return out << "Event does not implement write";
  }
protected:
  KDS_handle kds_;
};

} } //namespace CGAL::Kinetic
#endif
