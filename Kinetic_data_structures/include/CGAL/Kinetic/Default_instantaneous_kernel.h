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

#ifndef CGAL_CARTESIAN_DEFAULT_INSTANTANEOUS_KERNEL_H
#define CGAL_CARTESIAN_DEFAULT_INSTANTANEOUS_KERNEL_H

#include <CGAL/Kinetic/basic.h>
#include <CGAL/Kinetic/internal/Instantaneous_adaptor.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>
#include <map>
#include <iostream>
#include <CGAL/Kinetic/Ref_counted.h>
#include <CGAL/Kinetic/internal/To_static.h>

#define CGAL_MSA(Pred, pred, Arg, d) typedef Instantaneous_adaptor<typename Static_kernel::Pred##_##d, typename Kinetic_kernel::Pred##_##d, Rep, Arg> Pred##_##d; \
  Pred##_##d pred##_##d##_object() const				\
  {									\
    typename Static_kernel::Pred##_##d sp= rep_->static_kernel().pred##_##d##_object();	\
    typename Kinetic_kernel::Pred##_##d kp= rep_->kinetic_kernel().pred##_##d##_object();	\
    return Pred##_##d(rep_, sp, kp);		\
  }

#define CGAL_TSO(name) typedef typename Static_kernel::name name

namespace CGAL { namespace Kinetic {

template <class CIK>
class Default_instantaneous_kernel_rep: public Ref_counted<Default_instantaneous_kernel_rep<CIK> >
{
public:
  typedef typename CIK::Traits Traits;
  typedef typename Traits::Static_kernel Static_kernel;
  typedef typename Traits::Kinetic_kernel Kinetic_kernel;

  typedef typename Kinetic_kernel::Point_1::template Static_traits<Static_kernel> Static_traits_point_1;
  typedef typename Kinetic_kernel::Point_2::template Static_traits<Static_kernel> Static_traits_point_2;
  typedef typename Kinetic_kernel::Point_3::template Static_traits<Static_kernel> Static_traits_point_3;
  typedef typename Kinetic_kernel::Weighted_point_3::template Static_traits<Static_kernel> Static_traits_weighted_point_3;

  typedef typename Static_kernel::FT NT;
  typedef typename CIK::Traits::Simulator::Time Time;

  Default_instantaneous_kernel_rep(Traits tr): tr_(tr) {
    initialized_=false;
    time_is_nt_=false;
  }
  
  template <class T>
  void set_time(const T &t, bool limit) {
    if (!initialized_) {
      time_is_nt_=false;
      time_=t;
    } else {
      if ((!time_is_nt_ && time_ != t) || time_is_nt_) {
        time_is_nt_=false;
        time_=t;
        cache_1_.clear();
        cache_2_.clear();
        cache_3_.clear();
        cache_w3_.clear();
      }
    }
    initialized_=true;
    after_=limit;
  }

  void set_time(const NT &t, bool limit)
  {
    if (initialized_ && ((time_is_nt_ && time_nt_ != t) || !time_is_nt_)) {
      cache_1_.clear();
      cache_2_.clear();
      cache_3_.clear();
      cache_w3_.clear();
    }
    time_is_nt_=true;
    time_nt_=t;
    time_= Time(t);

    initialized_=true;
    after_=limit;
  }

  bool time_after() const {
    return after_;
  }

  bool initialized() const {
    return initialized_;
  }
  bool time_is_nt() const {
    return time_is_nt_;
  }

  const NT & time_as_nt() const
  {
    CGAL_precondition(initialized_);
    CGAL_precondition(time_is_nt());
    return time_nt_;
  }
  
  const Time &time() const
  {
    CGAL_precondition(initialized_);
    return time_;
  }
  
  void check_static_object() const {
#ifndef NDEBUG
    if (!initialized_) {
      std::cerr << "The InstantaneousKernel (or one of its predicates) was\n";
      std::cerr << "used without the time being set. This probably is the sign\n";
      std::cerr << "of misusing it--specifically make sure you store a copy\n";
      std::cerr << "from the SimulatorTraits and get predicates from it.\n";
    }
#endif
    if (!time_is_nt()) {
      std::cerr << "You can only compute static objects when the IK current\n";
      std::cerr << "time is an FT, rather than a root.\n";
      CGAL_error();
    }
  }

  typedef typename CIK::Traits::Active_points_1_table::Data::template Static_traits<Static_kernel> Static_traits_1;
  typedef typename CIK::Traits::Active_points_2_table::Data::template Static_traits<Static_kernel> Static_traits_2;
  typedef typename CIK::Traits::Active_points_3_table::Data::template Static_traits<Static_kernel> Static_traits_3;
  typedef typename CIK::Traits::Active_weighted_points_3_table::Data::template Static_traits<Static_kernel> Weighted_static_traits_3;

  const typename Static_traits_1::Static_type&
  static_object(typename CIK::Point_1 k) const {
    check_static_object();
    if (cache_1_.find(k) == cache_1_.end()) {
      cache_1_[k]= Static_traits_1::to_static(tr_.active_points_1_table_handle()->at(k),
                                              time_nt_, static_kernel());
    }
    return cache_1_[k];
  }

  const typename Static_traits_2::Static_type&
  static_object(typename CIK::Point_2 k) const {
    check_static_object();
    if (cache_2_.find(k) == cache_2_.end()) {
      cache_2_[k]= Static_traits_2::to_static(tr_.active_points_2_table_handle()->at(k),
                                              time_nt_, static_kernel());
    }
    return cache_2_[k];
  }

  const typename Static_traits_3::Static_type&
  static_object(typename CIK::Point_3 k) const {
    check_static_object();
    if (cache_3_.find(k) == cache_3_.end()) {
      cache_3_[k]= Static_traits_3::to_static(tr_.active_points_3_table_handle()->at(k),
                                              time_nt_, static_kernel());
    }
    return cache_3_[k];
  }

  const typename Weighted_static_traits_3::Static_type&
  static_object(typename CIK::Weighted_point_3 k) const {
    check_static_object();
    if (cache_w3_.find(k) == cache_w3_.end()) {
      cache_w3_[k]= Weighted_static_traits_3::to_static(tr_.active_weighted_points_3_table_handle()->at(k),
                time_nt_, static_kernel());
    }
    return cache_w3_[k];
  }

  const typename CIK::Traits::Active_points_1_table::Data&
  kinetic_object(typename CIK::Point_1 k) const {
    return tr_.active_points_1_table_handle()->at(k);
  }

  const typename CIK::Traits::Active_points_2_table::Data&
  kinetic_object(typename CIK::Point_2 k) const {
    return tr_.active_points_2_table_handle()->at(k);
  }

  const typename CIK::Traits::Active_points_3_table::Data&
  kinetic_object(typename CIK::Point_3 k) const {
    return tr_.active_points_3_table_handle()->at(k);
  }

  const typename CIK::Traits::Active_weighted_points_3_table::Data&
  kinetic_object(typename CIK::Weighted_point_3 k) const {
    return tr_.active_weighted_points_3_table_handle()->at(k);
  }
  
  const Static_kernel& static_kernel() const
  {
    return tr_.static_kernel_object();
  }
  
  const Kinetic_kernel& kinetic_kernel() const
  {
    return tr_.kinetic_kernel_object();
  }

protected:
  mutable bool initialized_;
  bool time_is_nt_;
  typename CIK::Traits tr_;
  mutable std::map<typename CIK::Point_1,
                   typename Static_traits_1::Static_type> cache_1_;
  mutable std::map<typename CIK::Point_2,
                   typename Static_traits_2::Static_type> cache_2_;
  mutable std::map<typename CIK::Point_3,
                   typename Static_traits_3::Static_type> cache_3_;
  mutable std::map<typename CIK::Weighted_point_3,
                   typename Weighted_static_traits_3::Static_type> cache_w3_;
  NT time_nt_;
  Time time_;
  bool after_;
};

template <class Traitst >
class Default_instantaneous_kernel
{
  typedef Default_instantaneous_kernel<Traitst> This;
public:
  typedef Traitst Traits;
  typedef Default_instantaneous_kernel_rep< This>  Rep;
  typedef typename Traits::Static_kernel Static_kernel;
  typedef typename Traits::Kinetic_kernel Kinetic_kernel;
  typedef typename Static_kernel::FT NT;
  typedef typename Traits::Simulator::Time Time;

  CGAL_static_assertion((boost::is_convertible<NT, Time>::value));
  CGAL_static_assertion((boost::is_convertible<Time, typename Kinetic_kernel::Certificate::Time>::value));
  CGAL_static_assertion((boost::is_convertible<typename Kinetic_kernel::Certificate::Time,
  Time>::value));

  Default_instantaneous_kernel(const Traits &tr):
    rep_(new Rep(tr)) {
  }
  template <class N>
  void set_time(const N &cur_time) const
  {
    rep_->set_time(cur_time, false);
  }

  void set_time(const Time &cur_time) const
  {
    rep_->set_time(cur_time, false);
  }

  template <class N>
  void set_time_to_after(const N &cur_time) const
  {
    rep_->set_time(cur_time, true);
  }

  void set_time_to_after(const Time &cur_time) const
  {
    rep_->set_time(cur_time, true);
  }
 
  bool time_is_nt() const
  {
    return rep_->time_is_nt();
  }

  const NT & time_as_nt() const
  {
    return rep_->time_as_nt();
  }
  
  const Time & time() const
  {
    return rep_->time();
  }

  bool has_time() const
  {
    return rep_->initialized();
  }

  typedef typename Static_kernel::RT RT;
  typedef typename Static_kernel::FT FT;

  typedef typename Traits::Active_points_1_table::Key Point_1;
  typedef typename Traits::Active_points_2_table::Key Point_2;
  typedef typename Traits::Active_points_3_table::Key Point_3;
  typedef typename Traits::Active_weighted_points_3_table::Key Weighted_point_3;

  struct Current_coordinates
  {
    Current_coordinates(typename Rep::Handle rep): rep_(rep){}
   
    const FT  & operator()(Point_1 k) const {
      return rep_->static_object(k);
    }
    const typename Static_kernel::Point_2  & operator()(Point_2 k) const {
      return rep_->static_object(k);
    }
    const typename Static_kernel::Point_3  & operator()(Point_3 k) const {
      return rep_->static_object(k);
    }
    const typename Static_kernel::Weighted_point_3  & operator()(Weighted_point_3 k) const {
      return rep_->static_object(k);
    }
    typename Rep::Handle rep_;
  };

  Current_coordinates current_coordinates_object() const {
    return Current_coordinates(rep_);
  }

  template <class T>
  class Compare_static {
  public:
    typedef CGAL::Comparison_result result_type;
    typedef T first_argument_type;
    typedef T second_argument_type;
    result_type operator()(const T &a, const T&b) const {
      return CGAL::compare(a,b);
    }
  };

  typedef Instantaneous_adaptor<Compare_static<RT>, typename Kinetic_kernel::Compare_x_1, Rep, Point_1> Compare_x_1;
  Compare_x_1 compare_x_1_object() const
  {
    Compare_static<NT> sp;
    return Compare_x_1(rep_, sp, rep_->kinetic_kernel().compare_x_1_object());
  }

  struct Construct_point_2
  {
    template<typename>
    struct result {
      typedef Point_2 type;
    };

    template<typename F>
    struct result<F(Point_2)> {
      typedef const Point_2& type;
    };

    const Point_2& operator()(const Point_2& p) const { return p; }
  };

  struct Construct_point_3
  {
    template<typename>
    struct result {
      typedef Point_3 type;
    };

    template<typename F>
    struct result<F(Point_3)> {
      typedef const Point_3& type;
    };

    const Point_3& operator()(const Point_3& p) const { return p; }
    Point_3 operator()(const Weighted_point_3& wp) const { return Point_3(wp.index()); }
  };

  Construct_point_2 construct_point_2_object() const { return Construct_point_2(); }
  Construct_point_3 construct_point_3_object() const { return Construct_point_3(); }

  CGAL_MSA(Side_of_oriented_circle,side_of_oriented_circle, Point_2, 2);
  CGAL_MSA(Orientation,orientation, Point_2, 2);
  CGAL_MSA(Compare_x, compare_x, Point_2, 2);
  CGAL_MSA(Compare_y, compare_y, Point_2, 2);
  CGAL_MSA(Compare_distance, compare_distance, Point_2, 2);
  CGAL_MSA(Compare_distance, compare_distance, Point_3, 3);
  CGAL_TSO(Segment_2);
  CGAL_TSO(Triangle_2);

  CGAL_MSA(Side_of_oriented_sphere,side_of_oriented_sphere, Point_3, 3);
  CGAL_MSA(Orientation,orientation, Point_3, 3);
  CGAL_MSA(Compare_x,compare_x, Point_3, 3);
  CGAL_MSA(Compare_y,compare_y, Point_3, 3);
  CGAL_MSA(Compare_z,compare_z, Point_3, 3);
  CGAL_MSA(Compare_xyz,compare_xyz, Point_3, 3);
  /*CGAL_MSA(Less_x, less_x, Point_3, 3);
  CGAL_MSA(Less_y, less_y, Point_3, 3);
  CGAL_MSA(Less_z, less_z, Point_3, 3);*/
  CGAL_MSA(Coplanar_orientation, coplanar_orientation, Point_3, 3);
  CGAL_MSA(Coplanar_side_of_bounded_circle, coplanar_side_of_bounded_circle, Point_3, 3);
  CGAL_MSA(Equal, equal, Point_3, 3);
  CGAL_MSA(Power_side_of_oriented_power_sphere,power_side_of_oriented_power_sphere, Weighted_point_3, 3);

  CGAL_TSO(Segment_3);
  CGAL_TSO(Triangle_3);
  CGAL_TSO(Tetrahedron_3);
  CGAL_TSO(Line_3);
  CGAL_TSO(Ray_3);
  CGAL_TSO(Object_3);
  CGAL_TSO(Plane_3);

protected:
  typename Rep::Handle rep_;
};
#undef CGAL_MSA
#undef CGAL_TSO
} } //namespace CGAL::Kinetic
#endif
