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

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4180) // qualifier applied to function type has no meaning; ignored
#endif

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
template<class Function,
  class BGT,
  class Wrapper = Implicit_to_labeling_function_wrapper<Function,BGT> >
class Implicit_mesh_domain_3
 : public Labeled_mesh_domain_3<Wrapper, BGT >
{
public:
  /// Base type
  typedef Labeled_mesh_domain_3<Wrapper, BGT> Base;

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
  Implicit_mesh_domain_3(const Function& f,
                         const Sphere_3& bounding_sphere,
                         const FT& error_bound = FT(1e-6),
                         CGAL::Random* p_rng = NULL)
    : Base(Wrapper(f), bounding_sphere, error_bound, p_rng)  {}

  /// Destructor
  virtual ~Implicit_mesh_domain_3() {}


private:
  // Disabled copy constructor & assignment operator
  typedef Implicit_mesh_domain_3<Function,BGT> Self;
  Implicit_mesh_domain_3(const Self& src);
  Self& operator=(const Self& src);

};  // end class Implicit_mesh_domain_3


}  // end namespace CGAL

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

#endif // CGAL_IMPLICIT_MESH_DOMAIN_3_H
