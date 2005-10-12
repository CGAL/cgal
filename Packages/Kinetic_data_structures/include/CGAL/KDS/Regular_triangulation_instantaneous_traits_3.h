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

#ifndef CGAL_REGULAR_CARTESIAN_INSTANTANEOUS_KERNEL_H
#define  CGAL_REGULAR_CARTESIAN_INSTANTANEOUS_KERNEL_H
#include <CGAL/KDS/Cartesian_instantaneous_kernel.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/KDS/internal/Kernel/Cartesian_moving_weighted_point_3.h>

#define CGAL_MSAW(Pred, pred, d) typedef Instantaneous_adaptor<typename Parent::Static_kernel::Pred##_##d, This, Weighted_point_3> Pred##_##d;\
Pred##_##d pred##_##d##_object() const {\
typename Parent::Static_kernel::Pred##_##d sp= Parent::rep_->static_kernel().pred##_##d##_object();\
return Pred##_##d(*this, sp);\
}



CGAL_KDS_BEGIN_NAMESPACE

//! A kernel that allows CGAL datastructures to act on a snapshot of a kinetic data structure. 
/*!  This version acts as a three dimensional regular triangulation
  traits. Use this if you want to compute a regular triangulation.

*/
template <class MPT, class SK = CGAL::Regular_triangulation_euclidean_traits_3<CGAL::Simple_cartesian<typename MPT::Object::Coordinate::NT> > > 
class Regular_triangulation_instantaneous_traits_3:
  public Cartesian_instantaneous_kernel<MPT, SK> {
  typedef Regular_triangulation_instantaneous_traits_3<MPT, SK> This;
  typedef Cartesian_instantaneous_kernel<MPT, SK> Parent;
public:
  
  Regular_triangulation_instantaneous_traits_3(typename MPT::Const_pointer mot,
					       const typename Parent::Static_kernel &sk= typename Parent::Static_kernel()): Parent(mot, sk){
  }
  //! Use a key from the moving point table as the point primitive.
  /*!
    The keys will get transformed into actual geometry when the geometry is needed.
  */
  typedef typename MPT::Key Weighted_point_3;
  //! A bare point is represented by the geometry.
  /*!  It shouldn't actually be used by the algorithm or data
    structure wrapping this kernel for anything other than debugging
    or temporaries.
  */
  typedef typename Parent::Static_kernel::Bare_point Bare_point;
  CGAL_MSAW(Power_test,power_test, 3);
  
  void write_point(const Weighted_point_3 p, std::ostream &out) const {
    out << p << " = " << Parent::rep_->static_object(p) << std::endl;;
  }
};
#undef CGAL_MSAW

CGAL_KDS_END_NAMESPACE
#endif
