// Copyright (c) 2014 GeometryFactory Sarl (France)
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
//
//
// Author(s)     : Jane Tournois

#ifndef CGAL_TIME_STAMPER_H
#define CGAL_TIME_STAMPER_H

namespace CGAL {

template <typename T>
struct Time_stamper
{
  Time_stamper()
   : time_stamp_(0) {}

  Time_stamper(const Time_stamper& ts)
   : time_stamp_(ts.time_stamp_) {}

  void set_time_stamp(T* pt) {
    pt->set_time_stamp(time_stamp_++);
  }

  static bool less(T* p_t1, T* p_t2) {
    if(p_t1 == NULL)      return (p_t2 != NULL);
    else if(p_t2 == NULL) return false;
    else                  return p_t1->time_stamp() < p_t2->time_stamp();
  }

  void reset() {
    time_stamp_ = 0;
  }
private:
  std::size_t time_stamp_;
}; // end class template Time_stamper<T>

template <typename T>
struct No_time_stamp
{
public:
  void set_time_stamp(T* pt)  {}
  static bool less(T* p_t1, T* p_t2) {
    return p_t1 < p_t2;
  }
  void reset()                {}
}; // end class template No_time_stamp<T>

} //end of CGAL namespace

#endif // CGAL_TIME_STAMPER_H
