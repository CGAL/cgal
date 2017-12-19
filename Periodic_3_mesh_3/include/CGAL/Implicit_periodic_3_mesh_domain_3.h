// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
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
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Mikhail Bogdanov
//
#ifndef CGAL_PERIODIC_3_MESH_3_IMPLICIT_PERIODIC_3_MESH_DOMAIN_3_H
#define CGAL_PERIODIC_3_MESH_3_IMPLICIT_PERIODIC_3_MESH_DOMAIN_3_H

#include <CGAL/license/Periodic_3_mesh_3.h>

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4180) // qualifier applied to function type has no meaning; ignored
#endif

#include <CGAL/Periodic_3_mesh_3/config.h>

#include <CGAL/Labeled_periodic_3_mesh_domain_3.h>
#include <CGAL/Implicit_to_labeling_function_wrapper.h>

namespace CGAL {
/**
 * @class Implicit_periodic_3_mesh_domain_3
 *
 * Implements mesh_traits for a domain defined as the negative values of
 * an implicit function.
 */
template<class Function,
         class BGT,
         class Wrapper = Implicit_to_labeling_function_wrapper<Function,BGT> >
class Implicit_periodic_3_mesh_domain_3
 : public Labeled_periodic_3_mesh_domain_3<Wrapper, BGT >
{
public:
  /// Base type
  typedef Labeled_periodic_3_mesh_domain_3<Wrapper, BGT> Base;

  /// Public types
  typedef typename Base::Bbox_3 Bbox_3;
  typedef typename Base::FT FT;
  typedef BGT Geom_traits;
  typedef typename BGT::Iso_cuboid_3 Iso_cuboid_3;

  /**
   * Constructor
   * @param f the function which negative values defines the domain
   * @param cuboid the fundamental domain
   * @param error_bound the error bound relative to the sphere radius
   */
  Implicit_periodic_3_mesh_domain_3(const Function& f,
                                    const Iso_cuboid_3& cuboid,
                                    FT error_bound = FT(1e-6))
    : Base(Wrapper(f), cuboid, error_bound)
  { }

  /// Destructor
  virtual ~Implicit_periodic_3_mesh_domain_3() { }

private:
  // Disabled copy constructor & assignment operator
  typedef Implicit_periodic_3_mesh_domain_3<Function,BGT> Self;
  Implicit_periodic_3_mesh_domain_3(const Self& src);
  Self& operator=(const Self& src);

};

} // end namespace CGAL

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

#endif // CGAL_PERIODIC_3_MESH_3_IMPLICIT_PERIODIC_3_MESH_DOMAIN_3_H
