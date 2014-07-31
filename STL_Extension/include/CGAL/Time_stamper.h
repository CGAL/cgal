// Copyright (c) 2014 GeometryFactory Sarl (France)
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
//
//
// Author(s)     : Jane Tournois

#include <CGAL/Has_timestamp.h>

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

  static bool less(const T* p_t1, const T* p_t2) {
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
  void set_time_stamp(T*)  {}
  static bool less(const T* p_t1,const T* p_t2) {
    return p_t1 < p_t2;
  }
  void reset()                {}
}; // end class template No_time_stamp<T>

// That class template is an auxiliary class.  It has a
// specialization for the case where `T::Has_timestamp` does not exists.
// The non-specialized template, when `T::Has_timestamp` exists, derives
// from `Time_stamper<T>` or `No_time_stamp<T>` depending on the
// value of the Boolean constant `T::Has_timestamp`.
// The declaration of that class template requires `T` to be a complete type.
template <class T, bool has_timestamp = internal::Has_timestamp<T>::value>
struct Get_time_stamper{
  typedef Time_stamper<T> type;
};

// Specialization when `T::Has_timestamp` does not exist, derives from
// `TimeStamper_`, or from `No_time_stamp<T>`.
template <class T>
struct Get_time_stamper<T,false>{
  typedef No_time_stamp<T> type;
};

// Implementation of the timestamp policy. It is very important that the
// declaration of that class template does not require `T` to be a complete
// type.  That way, the declaration of a pointer of type `Time_stamper_impl<T, Ts>
// in `Compact_container` for example is possible with an incomplete type.
template <class T>
struct Time_stamper_impl : public Get_time_stamper<T>::type {};

} //end of namespace CGAL

#endif // CGAL_TIME_STAMPER_H
