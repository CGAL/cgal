// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 1999, October 13
//
// file          : include/CGAL/Planar_map_2/Bounding_box_special_initializer.h
// package       : pm (4.08)
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Oren Nechushtan <theoren@math.tau.ac.il>
//
//
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter       : 
// ======================================================================

#ifndef CGAL_PM_ADVANCED_BOUNDING_BOX
template <class T_>
Bounding_box_base* init_default_bounding_box(T_*) const
{
  return new Pm_unbounding_box<Self>;
}

#ifdef CGAL_PM_STRAIGHT_EXACT_TRAITS_H
template <class R_>
Bounding_box_base* init_default_bounding_box(Pm_straight_exact_traits<R_>*)
const
{
  return new Pm_dynamic_open_bounding_box<Self>;
}
template <class R_>
Bounding_box_base* 
init_default_bounding_box(const Pm_straight_exact_traits<R_>*)
const
{
  return new Pm_dynamic_open_bounding_box<Self>;
}
template <class R_>
Bounding_box_base* init_default_bounding_box(Pm_straight_exact_traits<R_>*)
{
  return new Pm_dynamic_open_bounding_box<Self>;
}
template <class R_>
Bounding_box_base* 
init_default_bounding_box(const Pm_straight_exact_traits<R_>*)
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
  Pm_segment_exact_traits<R_>* t_) const
{
  return new Pm_unbounding_box<Self>;
}
template <class R_> 
Bounding_box_base* init_default_bounding_box(
  const Pm_segment_exact_traits<R_>* t_) const
{
  return new Pm_unbounding_box<Self>;
}
template <class R_> 
Bounding_box_base* init_default_bounding_box(
  Pm_segment_exact_traits<R_>* t_) 
{
  return new Pm_unbounding_box<Self>;
}
template <class R_> 
Bounding_box_base* init_default_bounding_box(
  const Pm_segment_exact_traits<R_>* t_) 
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
  const Pm_leda_segment_exact_traits<R_>* t_) const
{
  return new Pm_unbounding_box<Self>;
}
template <class R_> 
Bounding_box_base* init_default_bounding_box(
  Pm_leda_segment_exact_traits<R_>* t_) const
{
  return new Pm_unbounding_box<Self>;
}
template <class R_> 
Bounding_box_base* init_default_bounding_box(
  const Pm_leda_segment_exact_traits<R_>* t_) 
{
  return new Pm_unbounding_box<Self>;
}
template <class R_> 
Bounding_box_base* init_default_bounding_box(
  Pm_leda_segment_exact_traits<R_>* t_) 
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



