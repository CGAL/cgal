// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
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
// $URL: svn+ssh://mbogdanov@scm.gforge.inria.fr/svn/cgal/trunk/Mesh_3/include/CGAL/Implicit_mesh_domain_3.h $
// $Id: Implicit_mesh_domain_3.h 52705 2009-10-23 10:27:15Z stayeb $
//
//
// Author(s)     : Mikhail Bogdanov
//
//******************************************************************************
// File Description :
// class Periodic_implicit_mesh_domain_3. See class description.
//******************************************************************************

#ifndef CGAL_PERIODIC_IMPLICIT_MESH_DOMAIN_3_H
#define CGAL_PERIODIC_IMPLICIT_MESH_DOMAIN_3_H


#include <CGAL/Periodic_labeled_mesh_domain_3.h>
#include <CGAL/Implicit_to_labeling_function_wrapper.h>

#include <CGAL/Periodic_mesh_3/config.h>

namespace CGAL {


/**
 * @class Periodic_implicit_mesh_domain_3
 *
 * Implements mesh_traits for a domain defined as the negative values of
 * an implicit function.
 */
template<class Function,
  class BGT,
  class Wrapper = Implicit_to_labeling_function_wrapper<Function,BGT> >
class Periodic_implicit_mesh_domain_3
 : public Periodic_labeled_mesh_domain_3<Wrapper, BGT >
{
public:
  /// Base type
  typedef Periodic_labeled_mesh_domain_3<Wrapper, BGT> Base;

  /// Public types
  typedef typename Base::Bbox_3 Bbox_3;
  typedef typename Base::FT FT;
  typedef BGT Geom_traits;
  typedef typename BGT::Iso_cuboid_3 Iso_cuboid_3;

  /**
   * Constructor
   * @param f the function which negative values defines the domain
   * @param cuboid a bounding box of the domain, periodic domain
   * @param error_bound the error bound relative to the sphere radius
   */
  Periodic_implicit_mesh_domain_3(const Function& f,
                         const Iso_cuboid_3& cuboid,
                         FT error_bound = FT(1e-3))
    : Base(Wrapper(f), cuboid, error_bound) {}

  /// Destructor
  virtual ~Periodic_implicit_mesh_domain_3() {}

private:
  // Disabled copy constructor & assignment operator
  typedef Periodic_implicit_mesh_domain_3<Function,BGT> Self;
  Periodic_implicit_mesh_domain_3(const Self& src);
  Self& operator=(const Self& src);

};  // end class Periodic_implicit_mesh_domain_3


}  // end namespace CGAL

#endif // CGAL_PERIODIC_IMPLICIT_MESH_DOMAIN_3_H
