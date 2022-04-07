// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stéphane Tayeb
//
//******************************************************************************
// File Description :
// class Implicit_mesh_domain_3. See class description.
//******************************************************************************

#ifndef CGAL_IMPLICIT_MESH_DOMAIN_3_H
#define CGAL_IMPLICIT_MESH_DOMAIN_3_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Implicit_to_labeling_function_wrapper.h>
#include <CGAL/Random.h>

namespace CGAL {


/**
 * @class Implicit_mesh_domain_3
 *
 * Implements mesh_traits for a domain defined as the negative values of
 * an implicit function.
 */
template<class Function_,
  class BGT,
  class Wrapper = Implicit_to_labeling_function_wrapper<Function_,BGT> >
class
CGAL_DEPRECATED_MSG
( "The class template `CGAL::Implicit_mesh_domain_3` is now deprecated. "
  "Use the static member function template "
  "`Labeled_mesh_domain_3<K>::create_implicit_image_mesh_domain` instead.")
Implicit_mesh_domain_3
 : public Labeled_mesh_domain_3<BGT>
{
public:
  /// Base type
  typedef Labeled_mesh_domain_3<BGT> Base;

  /// Public types
  typedef typename Base::Sphere_3 Sphere_3;
  typedef typename Base::FT FT;
  typedef BGT Geom_traits;

  /**
   * Constructor
   * @param f the function which negative values defines the domain
   * @param bounding_sphere a bounding sphere of the domain
   * @param error_bound the error bound relative to the sphere radius
   */
  Implicit_mesh_domain_3(Function_ f,
                         const Sphere_3& bounding_sphere,
                         const FT& error_bound = FT(1e-6),
                         CGAL::Random* p_rng = nullptr)
    : Base(Wrapper(f), bounding_sphere, error_bound,
           Null_subdomain_index(), p_rng)  {}

  /// Destructor
  virtual ~Implicit_mesh_domain_3() {}

  using Base::bbox;
private:
  // Disabled copy constructor & assignment operator
  typedef Implicit_mesh_domain_3<Function_,BGT> Self;
  Implicit_mesh_domain_3(const Self& src);
  Self& operator=(const Self& src);

};  // end class Implicit_mesh_domain_3


}  // end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_IMPLICIT_MESH_DOMAIN_3_H
