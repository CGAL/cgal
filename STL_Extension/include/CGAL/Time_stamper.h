// Copyright (c) 2014 GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Jane Tournois

#ifndef CGAL_TIME_STAMPER_H
#define CGAL_TIME_STAMPER_H

#include <CGAL/Has_timestamp.h>
#include <CGAL/atomic.h>

namespace CGAL {

namespace internal {

constexpr size_t rounded_down_log2(size_t n)
{
  return ( (n<2) ? 0 : 1+rounded_down_log2(n/2));
}
} // namespace internal

template <typename T>
struct Time_stamper
{
  static void initialize_time_stamp(T* pt) {
    pt->set_time_stamp(std::size_t(-1));
  }

  template <typename time_stamp_t>
  static void set_time_stamp(T* pt, time_stamp_t& time_stamp_) {
    if(pt->time_stamp() == std::size_t(-1)) {
      const std::size_t new_ts = time_stamp_++;
      pt->set_time_stamp(new_ts);
    }
    else {
      // else: the time stamp is re-used

      // Enforces that the time stamp is greater than the current value.
      // That is used when a TDS_3 is copied: in that case, the
      // time stamps are copied from the old element to the new one,
      // but the time stamper does not know.
#ifdef CGAL_NO_ATOMIC
      time_stamp_ = (std::max)(time_stamp_, pt->time_stamp() + 1);
#else
      // How to atomically update a maximum value?
      //   https://stackoverflow.com/a/16190791/1728537
      const std::size_t new_value = pt->time_stamp() + 1;
      std::size_t prev_value = time_stamp_;
      while(prev_value < new_value &&
            !time_stamp_.compare_exchange_weak(prev_value, new_value))
        ;
#endif // atomic
    }
  }

  static std::size_t time_stamp(const T* pt)
  {
    if(pt == nullptr){
      return std::size_t(-1);
    }
    return pt->time_stamp();
  }

  static std::size_t hash_value(const T* p) {
    if(nullptr == p)
      return std::size_t(-1);
    else
      return p->time_stamp();
  }

  static bool less(const T* p_t1, const T* p_t2) {
    if(p_t1 == nullptr)      return (p_t2 != nullptr);
    else if(p_t2 == nullptr) return false;
    else {
      CGAL_assertion((p_t1 == p_t2) == (time_stamp(p_t1) == time_stamp(p_t2)));
      return time_stamp(p_t1) < time_stamp(p_t2);
    }
  }
}; // end class template Time_stamper<T>

template <typename T>
struct No_time_stamp
{
public:
  template <typename time_stamp_t>
  static void set_time_stamp(T*, time_stamp_t&)  {}
  static bool less(const T* p_t1,const T* p_t2) {
    return p_t1 < p_t2;
  }

  static void initialize_time_stamp(T*) {
  }

  static std::size_t time_stamp(const T*)
  {
    return 0;
  }

  static std::size_t hash_value(const T* p) {

    constexpr std::size_t shift = internal::rounded_down_log2(sizeof(T));
    return reinterpret_cast<std::size_t>(p) >> shift;
  }
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

struct Hash_handles_with_or_without_timestamps
{
  template <typename Handle>
  std::size_t operator()(const Handle h) const
  {
    typedef typename std::iterator_traits<Handle>::value_type Type;

    return Get_time_stamper<Type>::type::hash_value(&*h);
  }
};

} //end of namespace CGAL

#endif // CGAL_TIME_STAMPER_H
