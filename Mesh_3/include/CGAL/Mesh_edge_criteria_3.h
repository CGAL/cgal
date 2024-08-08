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

#ifndef CGAL_MESH_EDGE_CRITERIA_3_H
#define CGAL_MESH_EDGE_CRITERIA_3_H

#include <CGAL/license/Mesh_3.h>


#include <CGAL/Mesh_constant_domain_field_3.h>
#include <CGAL/Mesh_3/Is_mesh_domain_field_3.h>
#include <type_traits>
#include <cfloat>

namespace CGAL {
namespace Mesh_3 {
namespace internal {

  // Those two classes are designed to handle dynamic initialization of
  // Sizing_field type (using named parameters of make_mesh_3 for example)
  template < typename FT_, typename Point_, typename Index_ >
  class Sizing_field_interface
  {
  public:
    typedef FT_     FT;
    typedef Point_  Point_3;
    typedef Index_  Index;

    virtual ~Sizing_field_interface() {}

    virtual FT operator()(const Point_3& p,
                          const int dim,
                          const Index& index) const = 0;

    virtual Sizing_field_interface* clone() const = 0;
  };

  template < typename Sizing_field,
             typename FT,
             typename Point_3,
             typename Index>
  struct Sizing_field_container
    : public Sizing_field_interface < FT,
                                      Point_3,
                                      Index >
  {
    typedef Sizing_field_interface < FT,
                                     Point_3,
                                     Index > Base;

    typedef Sizing_field_container<Sizing_field, FT, Point_3, Index> Self;

  public:
    Sizing_field_container(const Sizing_field& s) : s_(s) {}
    virtual ~Sizing_field_container() {}

    virtual FT operator()(const Point_3& p,
                          const int dim,
                          const Index& index) const
    {
      return s_(p,dim,index);
    }

    virtual Base* clone() const
    {
      return new Self(*this);
    }

  private:
    Sizing_field s_;
  };

} // end namespace internal
} // end namespace Mesh_3


/*!
\ingroup PkgMesh3MeshClasses

The function object class `Mesh_edge_criteria_3` is a model of `MeshEdgeCriteria_3`. It
provides bounds for the size and approximation criteria.

\cgalModels{MeshEdgeCriteria_3}

\sa `MeshCriteriaWithFeatures_3`
\sa `CGAL::Mesh_criteria_3<Tr>`
\sa `MeshDomainField_3`
*/
template < typename Tr >
class Mesh_edge_criteria_3
{
private:
  typedef Mesh_edge_criteria_3 Self;
  typedef typename Tr::Geom_traits GT;

public:

  /// \name Types
  /// @{
  /*!
  Numerical type.
  */
  typedef typename Tr::Geom_traits::FT  FT;
  typedef typename Tr::Vertex::Index    Index;
  typedef typename Tr::Bare_point       Point_3;

  /// @}


  /// \name Creation
  /// @{

  /*!
  * returns an object to serve as criteria for edges.
  *
  * \param length_bound is an upper bound
  * for the length of the edges which are used to discretize the curves.
  * It can be a functional or a constant.
  * \param min_length_bound is a desired lower bound
  * for the length of the edges which are used to discretize the curves.
  * Only edges that are longer than this bound will be refined. Using
  * this lower bound can be handy on some domains, but using it may
  * break all the surface topology guarantees of the meshing algorithm.
  * It is not guaranteed to be exactly respected in the output mesh.
  * \param distance_bound is an upper bound for the distance from the
  * edge to the corresponding 1D feature.
  * It can be a functional or a constant.
  *
  * \note If one parameter is set to 0, then its corresponding criterion is ignored.
  *
  * @tparam SizingField scalar or model of `MeshDomainField_3`
  * @tparam DistanceField scalar or model of `MeshDomainField_3`
  */
  template < typename SizingField, typename DistanceField = FT >
  Mesh_edge_criteria_3(const SizingField& length_bound,
                       const FT& min_length_bound = 0,
                       const DistanceField& distance_bound = FT(0))
      : min_length_bound_(min_length_bound)
  {
    init_p_size(length_bound,
                  Mesh_3::Is_mesh_domain_field_3<Tr, SizingField>());
    init_distance_bound(distance_bound,
                  Mesh_3::Is_mesh_domain_field_3<Tr, DistanceField>());
  }

  /// @}

#ifndef DOXYGEN_RUNNING
  Mesh_edge_criteria_3(const Self& rhs)
    : p_size_(rhs.p_size_->clone())
    , min_length_bound_(rhs.min_length_bound_)
    , distance_bound_(rhs.distance_bound_ == nullptr
                      ? nullptr
                      : rhs.distance_bound_->clone())
  {}

  /// Destructor
  ~Mesh_edge_criteria_3()
  {
    delete p_size_;
    if(distance_bound_ != nullptr)
      delete distance_bound_;
  }

  /// Returns size of tuple (p,dim,index)
  FT sizing_field(const Point_3& p, const int dim, const Index& index) const
  {
    const FT s = (*p_size_)(p, dim, index);
    if (min_length_bound_ == FT(0))
      return s;
    else
      return (std::max)(s, min_length_bound_);
  }

  FT distance_field(const Point_3& p, const int dim, const Index& index) const
  {
    if (distance_bound_ == nullptr)
      return FT(0);
    return (*distance_bound_)(p,dim,index);
  }

public:
  const FT& min_length_bound() const
  {
    return min_length_bound_;
  }
  bool has_distance_field() const
  {
    return distance_bound_ != nullptr;
  }

#endif

private:
  typedef Mesh_3::internal::Sizing_field_interface<FT,Point_3,Index>
    Sizing_field_interface;

  void init_p_size(const FT& length_bound, Tag_false)
  {
    p_size_ = new Mesh_3::internal::Sizing_field_container<
        Mesh_constant_domain_field_3<GT,Index> ,
        FT,
        Point_3,
        Index>(length_bound);
  }

  template <typename SizingField>
  void init_p_size(const SizingField& length_bound, Tag_true)
  {
    p_size_ = new Mesh_3::internal::Sizing_field_container<
        SizingField,
        FT,
        Point_3,
        Index>(length_bound);
  }

  void init_distance_bound(const FT& distance_bound, Tag_false)
  {
    if (distance_bound == 0.)
      distance_bound_ = nullptr;
    else
      distance_bound_ = new Mesh_3::internal::Sizing_field_container<
        Mesh_constant_domain_field_3<GT,Index> ,
        FT,
        Point_3,
        Index>(distance_bound);
  }

  template <typename DistanceField>
  void init_distance_bound(const DistanceField& distance_bound, Tag_true)
  {
    distance_bound_ = new Mesh_3::internal::Sizing_field_container<
        DistanceField,
        FT,
        Point_3,
        Index>(distance_bound);
  }

  // A pointer to Sizing_field_interface to handle dynamic wrapping of
  // real SizingField type
  Sizing_field_interface* p_size_;
  const FT min_length_bound_;
  Sizing_field_interface* distance_bound_;
};


} // end namespace CGAL

#endif // CGAL_MESH_EDGE_CRITERIA_3_H
