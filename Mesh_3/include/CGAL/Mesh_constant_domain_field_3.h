// Copyright (c) 2010 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description :
//******************************************************************************

#ifndef CGAL_MESH_3_MESH_CONSTANT_DOMAIN_FIELD_3_H
#define CGAL_MESH_3_MESH_CONSTANT_DOMAIN_FIELD_3_H

#include <CGAL/license/Mesh_3.h>


#include <map>
#include <utility>
#include <CGAL/STL_Extension/internal/Has_nested_type_Bare_point.h>
#include <boost/mpl/eval_if.hpp>
#include <boost/mpl/identity.hpp>

namespace CGAL {

/*!
* @ingroup PkgMesh3DomainFields
* The class `Mesh_constant_domain_field_3` is a model of concept `MeshDomainField_3`. It provides
* a constant field accessible using queries on 3D-points.
*
* The class `Mesh_constant_domain_field_3` can also be customized through `set_size()` operations to become
* a piecewise constant field, i.e. a sizing field with a constant size on each subpart
* of the domain.
*
* @tparam GT is the geometric traits class. It must match the type `Triangulation::Geom_traits`,
* where `Triangulation` is the nested type of the model of `MeshComplex_3InTriangulation_3` used
* in the meshing process.
*
* @tparam Index_ is the type of index of the vertices of the triangulation.
* It must match the type `%Index` of the model of `MeshDomain_3` used in the meshing process.
*
* @cgalModels{MeshDomainField_3}
*/
template <typename GT, typename Index_>
class Mesh_constant_domain_field_3
{
public:
  /// \name Types
  /// @{

  /*!
  * Numerical type.
  */
  typedef typename GT::FT         FT;

  /*!

  * Point type.
  */
#ifdef DOXYGEN_RUNNING
  typedef GT::Point_3 Point_3;
#else
  typedef typename boost::mpl::eval_if_c<
      internal::Has_nested_type_Bare_point<GT>::value,
      typename internal::Bare_point_type<GT>,
      boost::mpl::identity<typename GT::Point_3>
    >::type                       Point_3;
  #endif

  /*!
  * Type of index of the vertices of the triangulation.
  */
  typedef Index_                  Index;

  /// @}

private:
  // Map to store field values
  typedef std::map<std::pair<int,Index>,FT> Values;

public:
  /// \name Creation
  /// @{

  /*!
  * Builds a constant domain field with size `size`.
  */
  Mesh_constant_domain_field_3(const FT& size) : d_(size) {}

  /// @}

  /// \name Operations
  /// @{

  /*!
  * Returns the size of query points of dimension `dim` and index `index`.
  */
  FT operator()(const Point_3&, const int dim, const Index& index) const
  {
    typename Values::const_iterator it = values_.find(std::make_pair(dim,index));
    if ( it != values_.end() ) { return it->second; }

    return d_;
  }


  /*!
  * Sets the size such that `operator()` for any query point
  * of dimension `dim` and index `index` returns `size`.
  */
  void set_size(const FT& size, const int dim, const Index& index)
  {
    values_.insert(std::make_pair(std::make_pair(dim,index),size));
  }
  /// @}

private:
  // default value
  FT d_;
  // Other values
  Values values_;
};

} // end namespace CGAL

#endif // CGAL_MESH_3_MESH_CONSTANT_DOMAIN_FIELD_3_H
