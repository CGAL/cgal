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
#include <CGAL/Mesh_3/Null_subdomain_index.h>
#include <CGAL/Random.h>

namespace CGAL {


/*!
\ingroup PkgMesh3Domains

\deprecated The class template `Implicit_mesh_domain_3` is deprecated
since CGAL-4.13, in favor of the class template `Labeled_mesh_domain_3` and
its static function
`Labeled_mesh_domain_3::create_implicit_mesh_domain()`.

The class `Implicit_mesh_domain_3` implements a domain whose bounding surface is
described
implicitly as the zero level set of a function.
The domain to be discretized is assumed to be the domain where
the function has negative values.
This class is a model of the concept `MeshDomain_3`.


\tparam Function_ provides the definition of the function.
This parameter stands for a model of the concept
`ImplicitFunction` described in the
surface mesh generation package.
The number types `Function::FT`
and `BGT::FT` are required to match.

\tparam BGT is a geometric traits which provides the basic operations to implement
intersection tests and computations
through a bisection method. This parameter must be instantiated
with a model of the concept `BisectionGeometricTraits_3`.

The constructor of `Implicit_mesh_domain_3`
takes as argument a bounding sphere which is required
to circumscribe the surface and to have its center inside the
domain.
This domain constructs intersection points between
the surface and segments/rays/lines by bisection. It needs an
`error_bound` such that the bisection process is stopped
when the query segment is smaller than the error bound.
The `error_bound` passed as argument to the domain constructor
is a relative error bound expressed as a ratio to the bounding sphere radius.

\cgalModels `MeshDomain_3`

\sa `BisectionGeometricTraits_3`
\sa `CGAL::make_mesh_3()`.

*/
template<class Function_
         ,class BGT
#ifndef DOXYGEN_RUNNING
         ,class Wrapper = Implicit_to_labeling_function_wrapper<Function_,BGT>
#endif
         >
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

  /// \name Creation
  /// @{

  /*!
    @param f is the object of type `Function_` that represents the implicit
    surface.

    @param bounding_sphere is a bounding sphere of the implicit surface. The
    value of `f` at the sphere center `c` must be
    negative: \f$ f(c)<0\f$.

    @param error_bound is the relative error bound
    used to compute intersection points between the implicit surface
    and query segments. The
    bisection is stopped when the length of the intersected
    segment is less than the product of `error_bound` by the
    radius of `bounding_sphere`.
  */
  Implicit_mesh_domain_3(Function_ f
                         ,const Sphere_3& bounding_sphere
                         ,const FT& error_bound = FT(1e-6)
#ifndef DOXYGEN_RUNNING
                         ,CGAL::Random* p_rng = nullptr
#endif
                         )
    : Base(parameters::function = Wrapper(f), parameters::bounding_object = bounding_sphere, parameters::relative_error_bound = error_bound,
           parameters::null_subdomain_index = Null_subdomain_index(), parameters::p_rng = p_rng)  {}

  /// @}

  // Destructor
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
