// Copyright (c) 1997  Tel-Aviv University (Israel).
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
// Author(s)     : Oren Nechushtan <theoren@math.tau.ac.il>

#ifndef CGAL_PM_ADVANCED_BOUNDING_BOX
template <class T_>
Bounding_box_base* init_default_bounding_box(T_*) const
{
  return new Pm_unbounding_box<Self>;
}

#ifdef CGAL_PM_STRAIGHT_EXACT_TRAITS_H
template <class R_>
Bounding_box_base* init_default_bounding_box(Pm_straight_traits_2<R_>*)
const
{
  return new Pm_dynamic_open_bounding_box<Self>;
}
template <class R_>
Bounding_box_base* 
init_default_bounding_box(const Pm_straight_traits_2<R_>*)
const
{
  return new Pm_dynamic_open_bounding_box<Self>;
}
template <class R_>
Bounding_box_base* init_default_bounding_box(Pm_straight_traits_2<R_>*)
{
  return new Pm_dynamic_open_bounding_box<Self>;
}
template <class R_>
Bounding_box_base* 
init_default_bounding_box(const Pm_straight_traits_2<R_>*)
{
  return new Pm_dynamic_open_bounding_box<Self>;
}
#endif //CGAL_PM_STRAIGHT_EXACT_TRAITS_H
#else // CGAL_PM_ADVANCED_BOUNDING_BOX
template <class T_>
Bounding_box_base* init_default_bounding_box(T_*) const
{
  return new Pm_dynamic_open_bounding_box<Self>;
}

// special initializers for backward compatability.
#ifdef CGAL_PM_SEGMENT_EXACT_TRAITS_H
template <class R_> 
Bounding_box_base* init_default_bounding_box(
  Pm_segment_traits_2<R_>* t_) const
{
  return new Pm_unbounding_box<Self>;
}
template <class R_> 
Bounding_box_base* init_default_bounding_box(
  const Pm_segment_traits_2<R_>* t_) const
{
  return new Pm_unbounding_box<Self>;
}
template <class R_> 
Bounding_box_base* init_default_bounding_box(
  Pm_segment_traits_2<R_>* t_) 
{
  return new Pm_unbounding_box<Self>;
}
template <class R_> 
Bounding_box_base* init_default_bounding_box(
  const Pm_segment_traits_2<R_>* t_) 
{
  return new Pm_unbounding_box<Self>;
}
#endif
#ifdef CGAL_PM_SEGMENT_EPSILON_TRAITS_H
template <class R_> 
Bounding_box_base* init_default_bounding_box(
  const Pm_segment_epsilon_traits<R_>*  t_) const
{
  return new Pm_unbounding_box<Self>;
}
template <class R_> 
Bounding_box_base* init_default_bounding_box(
  Pm_segment_epsilon_traits<R_>*  t_) const
{
  return new Pm_unbounding_box<Self>;
}
template <class R_> 
Bounding_box_base* init_default_bounding_box(
  const Pm_segment_epsilon_traits<R_>*  t_)
{
  return new Pm_unbounding_box<Self>;
}
template <class R_> 
Bounding_box_base* init_default_bounding_box(
  Pm_segment_epsilon_traits<R_>*  t_)
{
  return new Pm_unbounding_box<Self>;
}
#endif
#ifdef CGAL_PM_LEDA_SEGMENT_EXACT_TRAITS_H
template <class R_> 
Bounding_box_base* init_default_bounding_box(
  const Pm_leda_segment_traits_2<R_>* t_) const
{
  return new Pm_unbounding_box<Self>;
}
template <class R_> 
Bounding_box_base* init_default_bounding_box(
  Pm_leda_segment_traits_2<R_>* t_) const
{
  return new Pm_unbounding_box<Self>;
}
template <class R_> 
Bounding_box_base* init_default_bounding_box(
  const Pm_leda_segment_traits_2<R_>* t_) 
{
  return new Pm_unbounding_box<Self>;
}
template <class R_> 
Bounding_box_base* init_default_bounding_box(
  Pm_leda_segment_traits_2<R_>* t_) 
{
  return new Pm_unbounding_box<Self>;
}
#endif
#ifdef CGAL_PM_TRAITS_CHECKER_H
template <class T1_,class T2_,class C_> 
Bounding_box_base* init_default_bounding_box(
  const Pm_traits_checker<T1_,T2_,C_>*  t_) const
{
  return init_default_bounding_box((const T1_* ) t_);
}
template <class T1_,class T2_,class C_> 
Bounding_box_base* init_default_bounding_box(
  Pm_traits_checker<T1_,T2_,C_>*  t_) const
{
  return init_default_bounding_box((T1_* ) t_);
}
template <class T1_,class T2_,class C_> 
Bounding_box_base* init_default_bounding_box(
  const Pm_traits_checker<T1_,T2_,C_>*  t_) 
{
  return init_default_bounding_box((const T1_* ) t_);
}
template <class T1_,class T2_,class C_> 
Bounding_box_base* init_default_bounding_box(
  Pm_traits_checker<T1_,T2_,C_>*  t_) 
{
  return init_default_bounding_box((T1_* ) t_);
}
#endif
#endif // CGAL_PM_ADVANCED_BOUNDING_BOX



