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
// Author(s)     : Stéphane Tayeb, Aymeric PELLE
//
//******************************************************************************
// File Description :
// class Labeled_mesh_domain_3. See class description.
//******************************************************************************

#ifndef CGAL_LABELED_MESH_DOMAIN_3_H
#define CGAL_LABELED_MESH_DOMAIN_3_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Mesh_3/config.h>

#include <CGAL/Bbox_3.h>
#include <CGAL/point_generators_3.h>
#include <memory>
#include <CGAL/tuple.h>
#include <CGAL/Origin.h>

#include <functional>
#include <type_traits>

#include <CGAL/SMDS_3/internal/Handle_IO_for_pair_of_int.h>
#include <CGAL/SMDS_3/internal/indices_management.h>

// support for `CGAL::Image_3`
#include <CGAL/Image_3.h>
#include <CGAL/Mesh_3/Image_to_labeled_function_wrapper.h>
#include <CGAL/Mesh_3/Image_plus_weights_to_labeled_function_wrapper.h>

// support for implicit functions
#include <CGAL/Implicit_to_labeling_function_wrapper.h>
#include <CGAL/Named_function_parameters.h>
#ifdef CGAL_MESH_3_VERBOSE
#  include <boost/format.hpp>
#endif
#include <optional>

#include <CGAL/Mesh_3/Null_subdomain_index.h>
#include <CGAL/Mesh_domain_with_polyline_features_3.h>

#include <CGAL/Mesh_3/polylines_to_protect.h>

namespace CGAL {
namespace Mesh_3 {
namespace internal {

  struct Identity {
    template <typename T>
    const T& operator()(const T& x) { return x; }
  };

  template<typename T>
  struct Greater_than {
    typedef T argument_type;
    typedef bool result_type;
    Greater_than(const T& second) : second(second) {}
    bool operator()(const T& first) const {
      return std::greater<T>()(first, second);
    }
    T second;
  };

  struct Do_not_delete {
    template <typename T> void operator()(T*) const { }
  };

  // Returns a box enclosing image `im`
  inline Bbox_3 compute_bounding_box(const Image_3& im)
  {
    return Bbox_3(-1+im.tx(),-1+im.ty(),-1+im.tz(),
                  double(im.xdim())*im.vx()+im.tx()+1,
                  double(im.ydim())*im.vy()+im.ty()+1,
                  double(im.zdim())*im.vz()+im.tz()+1);
  }

  template <typename Image_values_to_subdom_indices>
  struct Create_gray_image_values_to_subdomain_indices {
    typedef Image_values_to_subdom_indices type;
    template <typename FT>
    type operator()(Image_values_to_subdom_indices functor, const FT&) const {
      return functor;
    }
  };

  // specialization for `Null_functor`: create the default functor
  template <>
  struct Create_gray_image_values_to_subdomain_indices<Null_functor> {
    typedef Mesh_3::internal::Greater_than<double> type;
    template <typename FT>
    type operator()(Null_functor, const FT& iso_value) const {
      return type(iso_value);
    }
  };

  template <typename Image_values_to_subdom_indices>
  struct Create_labeled_image_values_to_subdomain_indices {
    typedef Image_values_to_subdom_indices type;
    type operator()(Image_values_to_subdom_indices functor) const {
      return functor;
    }
  };

  // specialization for `Null_functor`: create the default functor
  template <>
  struct Create_labeled_image_values_to_subdomain_indices<Null_functor> {
    typedef Identity type;
    type operator()(Null_functor) const {
      return type();
    }
  };

  // Detect_features_in_domain
  template<typename Point, typename DetectFunctor>
  struct Detect_features_in_domain {

    std::vector<std::vector<Point>>
    operator()(const CGAL::Image_3& image, CGAL::Image_3& weights, DetectFunctor functor) const {
#if defined(BOOST_MSVC) && (BOOST_MSVC < 1910) //before msvc2017
      if(weights.is_valid())
        return functor.operator()<Point>(image, weights);
      else
        return functor.operator()<Point>(image);
#else
      if(weights.is_valid())
        return functor.template operator()<Point>(image, weights);
      else
        return functor.template operator()<Point>(image);
#endif
    }
  };
  // specialization for `Null_functor`: create the default functor
  template<typename Point>
  struct Detect_features_in_domain<Point, Null_functor> {
    std::vector<std::vector<Point>>
    operator()(const CGAL::Image_3&, CGAL::Image_3&, Null_functor) const {
      return std::vector<std::vector<Point>>();
    }
  };

  template<typename Point, typename DetectFunctor>
  std::vector<std::vector<Point>>
    detect_features(const CGAL::Image_3& image, CGAL::Image_3& weights, DetectFunctor functor)
  {
    Detect_features_in_domain<Point, DetectFunctor> detector;
    return detector(image, weights, functor);
  }

  template<bool WithFeatures>
  struct Add_features_in_domain {
    template<typename MeshDomain, typename InputFeatureRange, typename DetectFunctor>
    void operator()(const CGAL::Image_3&, CGAL::Image_3&, MeshDomain&, const InputFeatureRange&, DetectFunctor)
    {}
  };

  template<>
  struct Add_features_in_domain<true>
  {
    template<typename MeshDomain, typename InputFeatureRange, typename DetectFunctor>
    void operator()(const CGAL::Image_3& image,
                    CGAL::Image_3& weights,
                    MeshDomain& domain,
                    const InputFeatureRange& input_features,
                    DetectFunctor functor)
    {
      using P = typename MeshDomain::Point_3;
      auto detected_feature_range
        = CGAL::Mesh_3::internal::detect_features<P>(image, weights, functor);

      CGAL::merge_and_snap_polylines(image, detected_feature_range, input_features);

      if (!input_features.empty())
        domain.add_features(input_features.begin(), input_features.end());
      domain.add_features(detected_feature_range.begin(), detected_feature_range.end());
    }
  };

} // end namespace CGAL::Mesh_3::internal
} // end namespace CGAL::Mesh_3

#ifndef DOXYGEN_RUNNING
template <typename Subdomain_index>
struct Construct_pair_from_subdomain_indices {
  typedef std::pair<Subdomain_index, Subdomain_index> result_type;

  result_type operator()(Subdomain_index a, Subdomain_index b) const {
    return result_type(a, b);
  }
}; // end class template Construct_pair_from_subdomain_indices
#endif

namespace details
{

template <typename Geom_traits,
          typename Subdomain_index,
          typename Surface_patch_index_>
class Labeled_mesh_domain_3_impl
{
protected:
  typedef Surface_patch_index_ Surface_patch_index;
  typedef typename Geom_traits::Point_3 Point_3;
  typedef typename Geom_traits::Sphere_3 Sphere_3;
  typedef typename Geom_traits::Iso_cuboid_3 Iso_cuboid_3;
  typedef typename Geom_traits::FT FT;
  typedef std::shared_ptr<CGAL::Random> CGAL_Random_share_ptr_t;
  // Returns squared error bound from `bbox` and `error`
  FT squared_error_bound(const Iso_cuboid_3& bbox, const FT& error) const
  {
    typename Geom_traits::Compute_squared_distance_3 squared_distance =
                             Geom_traits().compute_squared_distance_3_object();
    return squared_distance((bbox.min)(), (bbox.max)())*error*error/4;
  }

  static Iso_cuboid_3 iso_cuboid(const Bbox_3& bbox)
  {
    const Point_3 p_min(bbox.xmin(), bbox.ymin(), bbox.zmin());
    const Point_3 p_max(bbox.xmax(), bbox.ymax(), bbox.zmax());

    return Iso_cuboid_3(p_min,p_max);
  }

  static Iso_cuboid_3 iso_cuboid(const typename Geom_traits::Sphere_3& sphere)
  {
    return iso_cuboid(sphere.bbox());
  }

  static Iso_cuboid_3 iso_cuboid(const typename Geom_traits::Iso_cuboid_3& c)
  {
    return c;
  }

  static Construct_pair_from_subdomain_indices<Subdomain_index>
  construct_pair_functor() {
    return Construct_pair_from_subdomain_indices<Subdomain_index>();
  }

  template <typename Function,
            typename Bounding_object,
            typename Null,
            typename Construct_surface_patch_index>
  Labeled_mesh_domain_3_impl(const Function& f,
                                     const Bounding_object& bounding,
                                     const FT& error_bound,
                                     Construct_surface_patch_index cstr_s_p_i,
                                     Null null,
                                     CGAL::Random* p_rng)
    : function_(f)
    , bbox_(iso_cuboid(bounding))
    , cstr_s_p_index(cstr_s_p_i)
    , null(null)
    , p_rng_(p_rng == 0 ?
             CGAL_Random_share_ptr_t(new CGAL::Random(0)) :
             CGAL_Random_share_ptr_t(p_rng, Mesh_3::internal::Do_not_delete()))
    , squared_error_bound_(squared_error_bound(bbox_,error_bound))
  {}

  // The function which answers subdomain queries
  typedef std::function<Subdomain_index(const Point_3&)> Function;
  Function function_;
  // The bounding box
  const Iso_cuboid_3 bbox_;

  typedef std::function<
    Surface_patch_index(Subdomain_index,
                        Subdomain_index)> Construct_surface_patch_index;
  Construct_surface_patch_index cstr_s_p_index;
  // The functor that decides which sub-domain indices correspond to the
  // outside of the domain.
  typedef std::function<bool(Subdomain_index)> Null;
  Null null;
  // The random number generator used by Construct_initial_points
  CGAL_Random_share_ptr_t p_rng_;
  // Error bound relative to sphere radius
  FT squared_error_bound_;
}; // Labeled_mesh_domain_3_impl

} // namespace details

/*!
\ingroup PkgMesh3Domains

\brief The class `Labeled_mesh_domain_3` implements indexed domains.

This class is a model of concept `MeshDomain_3`.

Any boundary facet is labeled <a,b>, with a<b, where a and b are the
tags of its incident subdomains.
Thus, a boundary facet of the domain is labeled <0,b>, where b!=0.

This class includes a <em>labeling function</em> that provides the index of the subdomain in which any
query point lies. An intersection between a segment and bounding
surfaces is detected when both segment endpoints are associated with different
values of subdomain indices. The intersection is then constructed by bisection.
The bisection stops when the query segment is shorter than an error bound
`e` given by the product of the
length of the diagonal of the bounding box (in world coordinates), or the radius of the bounding sphere, and
a relative error bound passed as argument to the constructor of `Labeled_mesh_domain_3`.

This class has a constructor taking a labeling function. It has also three
static template member functions that act as named constructors:
<ul>
<li>`create_gray_image_mesh_domain()`, to create a domain from a 3D gray image,
<li>`create_labeled_image_mesh_domain()`, to create a domain from a 3D labeled image, and
<li>`create_implicit_mesh_domain()`, to create a domain from an implicit function.
</ul>

\tparam BGT is a geometric traits class that provides
the basic operations to implement
intersection tests and intersection computations
through a bisection method. This parameter must be instantiated
with a model of the concept `BisectionGeometricTraits_3`.

\cgalHeading{Labeling function}

A labeling function `f` must return `0` if the point isn't located in any subdomain. The return type of labeling functions is an integer.

Let `p` be a Point.
<ul>
<li>`f(p)=0` means that `p` is outside domain.</li>
<li>`f(p)=a`, `a!=0` means that `p` is inside subdomain `a`.</li>
</ul>
`CGAL::Implicit_multi_domain_to_labeling_function_wrapper` is a good candidate for this template parameter
if there are several components to mesh.

The function type can be any model of the concept `Callable` compatible with the signature
`Subdomain_index(const %Point_3&)`: it can be a function, a function object, a lambda expression...
that takes a `%Point_3` as argument, and returns a type convertible to `Subdomain_index`.

\cgalModels{MeshDomain_3}

\sa `CGAL::Implicit_multi_domain_to_labeling_function_wrapper`
\sa `CGAL::make_mesh_3()`

*/
template<class BGT,
         class Subdomain_index_ = int,
         class Surface_patch_index_ = std::pair<Subdomain_index_,
                                                Subdomain_index_> >
class Labeled_mesh_domain_3
#ifndef DOXYGEN_RUNNING
  : protected details::Labeled_mesh_domain_3_impl<BGT,
                                                  Subdomain_index_,
                                                  Surface_patch_index_>
#endif
{
public:
  //-------------------------------------------------------
  // Index Types
  //-------------------------------------------------------
  // Type of indexes for cells of the input complex
/// \name Types
///@{
  /// The subdomain index of this model of `MeshDomain_3`
  typedef Subdomain_index_                  Subdomain_index;
  //
#ifdef DOXYGEN_RUNNING
  /// The type of object that stores the function using type-erasure.
  typedef std::function< Subdomain_index(const Point_3 &)> Labeling_function;
///@}

/// \name Types imported from the geometric traits class
///@{
  /// The point type of the geometric traits class
  typedef typename Geom_traits::Point_3      Point_3;
  /// The sphere type of the geometric traits class
  typedef typename Geom_traits::Sphere_3     Sphere_3;
  /// The iso-cuboid type of the geometric traits class
  typedef typename Geom_traits::Iso_cuboid_3 Iso_cuboid_3;
  /// The number type (a field type) of the geometric traits class
  typedef typename Geom_traits::FT           FT;
///@}
#else // DOXYGEN_RUNNING
  typedef std::optional<Subdomain_index>  Subdomain;

  // Type of indexes for cells of the input complex
  typedef Surface_patch_index_                  Surface_patch_index;
  typedef std::optional<Surface_patch_index>  Surface_patch;

  // Type of indexes to characterize the lowest dimensional face of the input
  // complex on which a vertex lie
  typedef typename CGAL::Mesh_3::internal::
    Index_generator<Subdomain_index, Surface_patch_index>::Index Index;

private:
  typedef details::Labeled_mesh_domain_3_impl<BGT,
                                              Subdomain_index,
                                              Surface_patch_index
                                             > Impl_details;
  typedef typename Impl_details::Null     Null;
  typedef typename Impl_details::Construct_surface_patch_index
                                          Construct_surface_patch_index;
  typedef typename Impl_details::Function Function;

public:
  // Geometric object types
  typedef typename BGT::Point_3    Point_3;
  typedef typename BGT::Segment_3  Segment_3;
  typedef typename BGT::Ray_3      Ray_3;
  typedef typename BGT::Line_3     Line_3;
  typedef typename BGT::Vector_3   Vector_3;
  typedef typename BGT::Sphere_3   Sphere_3;
  typedef CGAL::Bbox_3             Bbox_3;

protected:
  typedef typename BGT::Iso_cuboid_3 Iso_cuboid_3;

public:
  // Kernel_traits compatibility
  typedef BGT R;
  // access Function type from inherited class
  typedef Function Fct;

  typedef std::tuple<Point_3,Index,int> Intersection;


  typedef typename BGT::FT FT;
  typedef BGT Geom_traits;
  using Impl_details::construct_pair_functor;
#endif // DOXYGEN_RUNNING

/// \name Creation
/// @{
  /*!  \brief Construction from a function, a bounding object and a relative error bound.
   *
   * \tparam Function a type compatible with `Labeling_function`
   * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
   * \tparam Bounding_object either a bounding sphere (of type `Sphere_3`), a bounding box (type `Bbox_3`),
   *                         or a bounding `Iso_cuboid_3`
   *
   * \param function the labeling function
   * \param bounding_object the bounding object bounding the meshable space.
   * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below:
   *
   * \cgalNamedParamsBegin
   *   \cgalParamNBegin{relative_error_bound}
   *      \cgalParamDescription{the relative error bound used to compute intersection points between the implicit surface and query segments.
   *                            The bisection is stopped when the length of the intersected segment is less than the product
   *                            of `relative_error_bound` by the diameter of the bounding object.}
   *      \cgalParamDefault{FT(1e-3)}
   *   \cgalParamNEnd
   * \cgalNamedParamsEnd
   *
   * \cgalHeading{Example}
   * From the example (\ref Mesh_3/mesh_implicit_domains_2.cpp):
   * \snippet Mesh_3/mesh_implicit_domains_2.cpp Domain creation
   */
  template<typename Function, typename Bounding_object, typename CGAL_NP_TEMPLATE_PARAMETERS>
  Labeled_mesh_domain_3(const Function& function,
                        const Bounding_object& bounding_object,
                        const CGAL_NP_CLASS& np = parameters::default_values()
#ifndef DOXYGEN_RUNNING
                        , typename std::enable_if<!is_named_function_parameter<Function>>::type* = nullptr
#endif // DOXYGEN_RUNNING
                        )
  :Impl_details(function,
                bounding_object,
                parameters::choose_parameter(parameters::get_parameter(np, internal_np::error_bound), FT(1e-3)),
                parameters::choose_parameter(parameters::get_parameter(np, internal_np::surface_patch_index), construct_pair_functor()),
                parameters::choose_parameter(parameters::get_parameter(np, internal_np::null_subdomain_index_param), Null_subdomain_index()),
                parameters::choose_parameter(parameters::get_parameter(np, internal_np::rng), nullptr))
  {}
///@}

#ifndef DOXYGEN_RUNNING
  template<typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT>
  Labeled_mesh_domain_3(const CGAL_NP_CLASS& np)
  :Impl_details(parameters::get_parameter(np, internal_np::function_param),
                parameters::get_parameter(np, internal_np::bounding_object_param),
                parameters::choose_parameter(parameters::get_parameter(np, internal_np::error_bound), FT(1e-3)),
                parameters::choose_parameter(parameters::get_parameter(np, internal_np::surface_patch_index), construct_pair_functor()),
                parameters::choose_parameter(parameters::get_parameter(np, internal_np::null_subdomain_index_param), Null_subdomain_index()),
                parameters::choose_parameter(parameters::get_parameter(np, internal_np::rng), nullptr))
  {}

  // Overload handling parameters passed with operator=
  template<typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_1,
           typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_2,
           typename ... NP>
  Labeled_mesh_domain_3(const CGAL_NP_CLASS_1&  np1,
                        const CGAL_NP_CLASS_2&  np2,
                        const NP& ... nps)
    : Labeled_mesh_domain_3(internal_np::combine_named_parameters(np1, np2, nps...))
  {}


#ifndef CGAL_NO_DEPRECATED_CODE
  template<typename Function, typename Bounding_object>
#if !defined(BOOST_MSVC)
  CGAL_DEPRECATED
#endif // BOOST_MSVC
  Labeled_mesh_domain_3(const Function& function,
                        const Bounding_object& bounding_object,
                        double error_bound,
                        typename std::enable_if<!is_named_function_parameter<Function>>::type* = nullptr)
  : Labeled_mesh_domain_3(function,
                          bounding_object,
                          parameters::relative_error_bound(error_bound))
  {}
#endif // CGAL_NO_DEPRECATED_CODE
#endif // DOXYGEN_RUNNING

/// \name Creation of domains from 3D images
/// @{
  /*!
   * \brief Construction from a 3D gray image
   *
   * This static method is a <em>named constructor</em>. It constructs a domain
   * described by a 3D gray image. A 3D gray image is a grid of voxels,
   * where each voxel is associated with a gray level value.  Unless otherwise specified by the parameter `image_values_to_subdom_indices`, the domain to
   * be discretized is the union of voxels that lie inside a surface
   * described by an isolevel value, called \a isovalue. The voxels lying
   * inside the domain have gray level values that are larger than the
   * isovalue.
   *
   * The value of voxels is interpolated to a gray level value at any query point.
   *
   * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
   *
   * \param image_ the input 3D image.
   * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below:
   *
   * \cgalNamedParamsBegin
   *   \cgalParamNBegin{iso_value}
   *     \cgalParamDescription{the isovalue, inside
   *                           `image`, of the surface describing the boundary of the object to be
   *                            meshed.}
   *     \cgalParamDefault{0}
   *   \cgalParamNEnd
   *   \cgalParamNBegin{image_values_to_subdomain_indices}
   *     \cgalParamDescription{a function or a function object, compatible with the signature
   *                           `Subdomain_index(double)`. This function returns the subdomain index
   *                           corresponding to a pixel value. If this parameter is used, then the
   *                           parameter `iso_value` is ignored.}
   *     \cgalParamDefault{Null_functor()}
   *   \cgalParamNEnd
   *   \cgalParamNBegin{value_outside}
   *      \cgalParamDescription{the value attached to voxels
   *                            outside of the domain to be meshed. It should be lower than
   *                           `iso_value`.}
   *      \cgalParamDefault{0}
   *   \cgalParamNEnd
   *   \cgalParamNBegin{relative_error_bound}
   *      \cgalParamDescription{ is the relative error
   *                             bound, relative to the diameter of the box of the image.}
   *      \cgalParamDefault{FT(1e-3)}
   *   \cgalParamNEnd
   * \cgalNamedParamsEnd
   *
   * \cgalHeading{Examples}
   *
   * From the example (\ref Mesh_3/mesh_3D_gray_image.cpp):
   *
   * \snippet Mesh_3/mesh_3D_gray_image.cpp Domain creation
   *
   * From the example (\ref Mesh_3/mesh_3D_gray_vtk_image.cpp):
   *
   * \snippet Mesh_3/mesh_3D_gray_vtk_image.cpp Domain creation
   *
   */
  template<typename CGAL_NP_TEMPLATE_PARAMETERS>
  static Labeled_mesh_domain_3 create_gray_image_mesh_domain(const CGAL::Image_3& image_, const CGAL_NP_CLASS& np = parameters::default_values())
  {
    using parameters::get_parameter;
    using parameters::choose_parameter;
    auto iso_value_ = choose_parameter(get_parameter(np, internal_np::iso_value_param), 0);
    auto value_outside_ = choose_parameter(get_parameter(np, internal_np::voxel_value), 0);
    FT relative_error_bound_ = choose_parameter(get_parameter(np, internal_np::error_bound), FT(1e-3));
    auto image_values_to_subdomain_indices_ = choose_parameter(get_parameter(np, internal_np::image_subdomain_index), Null_functor());
    CGAL::Random* p_rng_ = choose_parameter(get_parameter(np, internal_np::rng), nullptr);
    auto null_subdomain_index_ = choose_parameter(get_parameter(np, internal_np::null_subdomain_index_param), Null_functor());
    auto construct_surface_patch_index_ = choose_parameter(get_parameter(np, internal_np::surface_patch_index), Null_functor());
    namespace p = CGAL::parameters;
    return Labeled_mesh_domain_3
              (p::function = create_gray_image_wrapper
                       (image_,
                        iso_value_,
                        image_values_to_subdomain_indices_,
                        value_outside_),
               p::bounding_object = Mesh_3::internal::compute_bounding_box(image_),
               p::relative_error_bound = relative_error_bound_,
               p::p_rng = p_rng_,
               p::null_subdomain_index =
                       create_null_subdomain_index(null_subdomain_index_),
               p::construct_surface_patch_index =
                       create_construct_surface_patch_index(construct_surface_patch_index_));

  }

  /*!
   * \brief Construction from a 3D labeled image
   *
   * This static method is a <em>named constructor</em>. It constructs a
   * domain described by a 3D labeled image. A 3D labeled image is a grid
   * of voxels, where each voxel is associated with an index (a subdomain
   * index) characterizing the subdomain in which the voxel lies. The
   * domain to be discretized is the union of voxels that have non-zero
   * values.
   *
   * \returns either a `Labeled_mesh_domain_3`,
   *   or a `Mesh_domain_with_polyline_features_3<Labeled_mesh_domain_3>`
   *   depending on whether one or more of the named parameters
   *   `features_detector` and `input_features` are provided.
   *
   * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
   * \param image_ the input 3D image.
   * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below:
   *
   * \cgalNamedParamsBegin
   *   \cgalParamNBegin{weights}
   *     \cgalParamDescription{a reference to an input 3D image that provides
   *                           weights associated to each voxel (the word type is `unsigned char`,
   *                           and the voxels values are integers between 0 and 255).
   *                           The weights image can be generated with `CGAL::Mesh_3::generate_label_weights()`.
   *                           Its dimensions must be the same as the dimensions of `parameters::image`.}
   *     \cgalParamDefault{CGAL::Image_3()}
   *     \cgalParamType{CGAL::Image_3&}
   *     \cgalParamExtra{if `features_detector` is provided, `weights` should be modified accordingly.
   *                     The available functors described in See \ref PkgMesh3FeatureDetection
   *                     implement the necessary modifications.}
   *     \cgalParamExtra{if `input_features` is provided, `weights` should be modified accordingly
   *                     to keep consistency of the output `MeshDomainWithFeatures_3`}
   *   \cgalParamNEnd
   *   \cgalParamNBegin{value_outside}
   *     \cgalParamDescription{the value attached to voxels
   *                           outside of the domain to be meshed. It should be lower than
   *                           `iso_value`.}
   *     \cgalParamDefault{0}
   *   \cgalParamNEnd
   *
   *   \cgalParamNBegin{relative_error_bound}
   *     \cgalParamDescription{ is the relative error
   *                            bound, relative to the diameter of the box of the image.}
   *     \cgalParamDefault{FT(1e-3)}
   *   \cgalParamNEnd
   *
   *   \cgalParamNBegin{features_detector}
   *    \cgalParamDescription{ a functor that implements
   *      `std::vector<std::vector<Point>> operator()(const Image_3& img) const`,
   *      and `std::vector<std::vector<Point>> operator()(const Image_3& img, Image_3& weights) const`,
   *      where `%Point` matches the mesh domain point type,
   *      that both return a range of detected polyline features for feature protection.
   *      Only one implementation is used, depending on whether the named parameter `weights`
   *      is provided or not.
   *      Polyline features are added to the domain for further feature protection.
   *      See \ref PkgMesh3FeatureDetection for available functors.}
   *    \cgalParamDefault{CGAL::Null_functor()}
   *    \cgalParamExtra{The return type of the function depends on whether this parameter
   *                    or `input_features` are provided or not.}
   *    \cgalParamExtra{If `weights` is provided, it must either be adapted to the detected features,
   *                    or post-processed during feature detection to keep consistency
   *                    of the output `MeshDomainWithFeatures_3`.
   *                    Available functors implement the necessary modifications.}
   *   \cgalParamNEnd
   *
   *   \cgalParamNBegin{input_features}
   *    \cgalParamDescription{ a `Range` of polyline features, represented as `Range`s of `Point_3`.
   *         Polyline features are added to the domain for further feature protection.
   *         Input polyline features must be different from the detected features
   *         and can intersect only at vertices, if they do. Otherwise,
   *         the meshing process may not terminate.}
   *    \cgalParamDefault{`std::vector<std::vector<Point_3>>()`}
   *    \cgalParamExtra{The return type of the function depends on whether this parameter
                        or `input_features` are provided or not.}
   *    \cgalParamExtra{It is recommended to pass a const-reference for this parameter,
   *                    possibly using `std::cref(polylines_range)` to avoid useless copies.}
   *    \cgalParamExtra{If `weights` is provided, it must be adapted to the input features,
   *                    to keep consistency of the output `MeshDomainWithFeatures_3`}
   *   \cgalParamNEnd
   *
   * \cgalNamedParamsEnd
   *
   * \cgalHeading{Example}
   *
   * From the example (\ref Mesh_3/mesh_3D_image.cpp):
   *
   * \snippet Mesh_3/mesh_3D_image.cpp Domain creation
   *
   * From the example (\ref Mesh_3/mesh_3D_weighted_image.cpp),
   * where the labeled image is used with a precomputed 3D image of weights :
   *
   * \snippet Mesh_3/mesh_3D_weighted_image.cpp Domain creation
   *
   * From the example (\ref Mesh_3/mesh_3D_image_with_detection_of_features.cpp)
   * where the features are detected in `image`:
   *
   * \snippet Mesh_3/mesh_3D_image_with_detection_of_features.cpp Domain creation
   *
   * From the example (\ref Mesh_3/mesh_3D_image_with_input_features.cpp)
   * where the features are provided by the user:
   *
   * \snippet Mesh_3/mesh_3D_image_with_input_features.cpp Domain creation
   *
   */
  template<typename CGAL_NP_TEMPLATE_PARAMETERS>
  static auto
  create_labeled_image_mesh_domain(const CGAL::Image_3& image_, const CGAL_NP_CLASS& np = parameters::default_values())
  {
    using parameters::get_parameter;
    using parameters::get_parameter_reference;
    using parameters::choose_parameter;

    auto iso_value_ = choose_parameter(get_parameter(np, internal_np::iso_value_param), 0);
    auto value_outside_ = choose_parameter(get_parameter(np, internal_np::voxel_value), 0);
    FT relative_error_bound_ = choose_parameter(get_parameter(np, internal_np::error_bound), FT(1e-3));
    auto image_values_to_subdomain_indices_ = choose_parameter(get_parameter(np, internal_np::image_subdomain_index), Null_functor());
    CGAL::Random* p_rng_ = choose_parameter(get_parameter(np, internal_np::rng), nullptr);
    auto null_subdomain_index_ = choose_parameter(get_parameter(np, internal_np::null_subdomain_index_param), Null_functor());
    auto construct_surface_patch_index_ = choose_parameter(get_parameter(np, internal_np::surface_patch_index), Null_functor());

    using Image_ref_type = typename internal_np::Lookup_named_param_def<internal_np::weights_param_t,
                                                                        CGAL_NP_CLASS,
                                                                        CGAL::Image_3>::reference;
    CGAL::Image_3 no_weights_;
    Image_ref_type weights_ = choose_parameter(get_parameter_reference(np, internal_np::weights_param), no_weights_);
    auto features_detector_ = choose_parameter(get_parameter(np, internal_np::features_detector_param), Null_functor());

    using Default_input_features = std::vector<std::vector<typename Labeled_mesh_domain_3::Point_3>>;
    using Input_features_ref_type = typename internal_np::Lookup_named_param_def<internal_np::input_features_param_t,
                                                                                 CGAL_NP_CLASS,
                                                                                 Default_input_features>::reference;
    Default_input_features empty_vec;
    Input_features_ref_type input_features_
          = choose_parameter(get_parameter_reference(np, internal_np::input_features_param), empty_vec);

    CGAL_USE(iso_value_);
    namespace p = CGAL::parameters;

    const bool use_weights = weights_.is_valid();
    auto image_wrapper = use_weights
      ? create_weighted_labeled_image_wrapper(image_,
                                              weights_,
                                              image_values_to_subdomain_indices_,
                                              value_outside_)
      : create_labeled_image_wrapper(image_,
                                     image_values_to_subdomain_indices_,
                                     value_outside_);

    // warning : keep Return_type consistent with actual return type
    const bool no_features
      =  CGAL::parameters::is_default_parameter<CGAL_NP_CLASS, internal_np::features_detector_param_t>::value
      && CGAL::parameters::is_default_parameter<CGAL_NP_CLASS, internal_np::input_features_param_t>::value;
    using Return_type = std::conditional_t <
      no_features,
      Labeled_mesh_domain_3,
      Mesh_domain_with_polyline_features_3<Labeled_mesh_domain_3>
    >;

    Return_type domain
      (p::function = image_wrapper,
       p::bounding_object = Mesh_3::internal::compute_bounding_box(image_),
       p::relative_error_bound = relative_error_bound_,
       p::p_rng = p_rng_,
       p::null_subdomain_index =
               create_null_subdomain_index(null_subdomain_index_),
       p::construct_surface_patch_index =
               create_construct_surface_patch_index(construct_surface_patch_index_));

    // features
    Mesh_3::internal::Add_features_in_domain<!no_features>()
      (image_, weights_, domain, input_features_, features_detector_);

    return domain;
  }
/// @}

#ifndef DOXYGEN_RUNNING
  template<typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT>
  static Labeled_mesh_domain_3 create_gray_image_mesh_domain(const CGAL_NP_CLASS& np)
  {
    static_assert(!parameters::is_default_parameter<CGAL_NP_CLASS, internal_np::image_3_param_t>::value, "Value for required parameter not found");
    using parameters::get_parameter;
    using parameters::get_parameter_reference;
    using parameters::choose_parameter;
    const CGAL::Image_3& image_ = get_parameter_reference(np,internal_np::image_3_param);
    auto iso_value_ = choose_parameter(get_parameter(np, internal_np::iso_value_param), 0);
    auto value_outside_ = choose_parameter(get_parameter(np, internal_np::voxel_value), 0);
    FT relative_error_bound_ = choose_parameter(get_parameter(np, internal_np::error_bound), FT(1e-3));
    auto image_values_to_subdomain_indices_ = choose_parameter(get_parameter(np, internal_np::image_subdomain_index), Null_functor());
    CGAL::Random* p_rng_ = choose_parameter(get_parameter(np, internal_np::rng), nullptr);
    auto null_subdomain_index_ = choose_parameter(get_parameter(np, internal_np::null_subdomain_index_param), Null_functor());
    auto construct_surface_patch_index_ = choose_parameter(get_parameter(np, internal_np::surface_patch_index), Null_functor());
    namespace p = CGAL::parameters;
    return Labeled_mesh_domain_3
            (p::function = create_gray_image_wrapper
                     (image_,
                      iso_value_,
                      image_values_to_subdomain_indices_,
                      value_outside_),
             p::bounding_object = Mesh_3::internal::compute_bounding_box(image_),
             p::relative_error_bound = relative_error_bound_,
             p::p_rng = p_rng_,
             p::null_subdomain_index =
                     create_null_subdomain_index(null_subdomain_index_),
             p::construct_surface_patch_index =
                     create_construct_surface_patch_index(construct_surface_patch_index_));

  }

  // Overload handling parameters passed with operator=
  template<typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_1,
           typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_2,
           typename ... NP>
  static Labeled_mesh_domain_3 create_gray_image_mesh_domain(const CGAL::Image_3& image_,
                                                             const CGAL_NP_CLASS_1&  np1,
                                                             const CGAL_NP_CLASS_2&  np2,
                                                             const NP& ... nps)
  {
    return create_gray_image_mesh_domain(image_, internal_np::combine_named_parameters(np1, np2, nps...));
  }
  // Overload handling parameters passed with operator=
  template<typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_1,
           typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_2,
           typename ... NP>
  static Labeled_mesh_domain_3 create_gray_image_mesh_domain(const CGAL_NP_CLASS_1&  np1,
                                                             const CGAL_NP_CLASS_2&  np2,
                                                             const NP& ... nps)
  {
    return create_gray_image_mesh_domain(internal_np::combine_named_parameters(np1, np2, nps...));
  }

  template<typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT>
  static auto create_labeled_image_mesh_domain(const CGAL_NP_CLASS& np)
  {
    static_assert(!parameters::is_default_parameter<CGAL_NP_CLASS, internal_np::image_3_param_t>::value, "Value for required parameter not found");
    using parameters::get_parameter_reference;
    const CGAL::Image_3& image_ = get_parameter_reference(np,internal_np::image_3_param);
    return create_labeled_image_mesh_domain(image_, np);
  }

  // Overload handling parameters passed with operator=
  template<typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_1,
           typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_2,
           typename ... NP>
  static auto create_labeled_image_mesh_domain(const CGAL::Image_3& image_,
                                               const CGAL_NP_CLASS_1&  np1,
                                               const CGAL_NP_CLASS_2&  np2,
                                               const NP& ... nps)
  {
    return create_labeled_image_mesh_domain(image_, internal_np::combine_named_parameters(np1, np2, nps...));
  }
  // Overload handling parameters passed with operator=
  template<typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_1,
           typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_2,
           typename ... NP>
  static auto create_labeled_image_mesh_domain(const CGAL_NP_CLASS_1& np1,
                                               const CGAL_NP_CLASS_2& np2,
                                               const NP& ... nps)
  {
    return create_labeled_image_mesh_domain(internal_np::combine_named_parameters(np1, np2, nps...));
  }
#endif

/// \name Creation of domains from implicit functions
/// @{

  /*!
   * \brief Construction from an implicit function
   *
   * This static method is a <em>named constructor</em>. It constructs a domain
   * whose bounding surface is described implicitly as the zero level set of a
   * function.  The domain to be discretized is assumed to be the domain where
   * the function has negative values.
   *
   * The method takes as argument a bounding sphere which is required to
   * circumscribe the surface and to have its center inside the domain.
   *
   * \tparam Function a type compatible with the signature `FT(Point_3)`: it takes a point as argument,
   *                  and returns a scalar value. That object must be model of `CopyConstructible`
   * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
   * \tparam Bounding_object either a bounding sphere (of type `Sphere_3`), a bounding box (type `Bbox_3`),
   *                         or a bounding `Iso_cuboid_3`
   *
   * \param function the implicit function
   * \param bounding_object object boundint the meshable domain and its center is inside the domain.
   * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below:
   *
   * \cgalNamedParamsBegin
   *   \cgalParamNBegin{relative_error_bound}
   *     \cgalParamDescription{ is the relative error
   *                            bound, relative to the diameter of the box of the image.}
   *     \cgalParamDefault{FT(1e-3)}
   *   \cgalParamNEnd
   * \cgalNamedParamsEnd
   *
   * \cgalHeading{Examples}
   *
   * From the example (\ref Mesh_3/mesh_implicit_sphere.cpp):
   *
   * \snippet Mesh_3/mesh_implicit_sphere.cpp Domain creation
   *
   * From the example (\ref Mesh_3/mesh_implicit_sphere_variable_size.cpp):
   *
   * \snippet Mesh_3/mesh_implicit_sphere_variable_size.cpp Domain creation
   *
   */
  template<typename Function, typename Bounding_object, typename CGAL_NP_TEMPLATE_PARAMETERS>
  static Labeled_mesh_domain_3 create_implicit_mesh_domain(const Function& function,
                                                           const Bounding_object& bounding_object,
                                                           const CGAL_NP_CLASS& np = parameters::default_values()
#ifndef DOXYGEN_RUNNING
                                                           , typename std::enable_if<!is_named_function_parameter<Function>>::type* = nullptr
#endif
)
  {
    using parameters::get_parameter;
    using parameters::choose_parameter;
    FT relative_error_bound_ = choose_parameter(get_parameter(np, internal_np::error_bound), FT(1e-3));
    CGAL::Random* p_rng_ = choose_parameter(get_parameter(np, internal_np::rng), nullptr);
    auto null_subdomain_index_ = choose_parameter(get_parameter(np, internal_np::null_subdomain_index_param), Null_functor());
    auto construct_surface_patch_index_ = choose_parameter(get_parameter(np, internal_np::surface_patch_index), Null_functor());
    namespace p = CGAL::parameters;
    return Labeled_mesh_domain_3
            (p::function = make_implicit_to_labeling_function_wrapper<BGT>(function),
             p::bounding_object = bounding_object,
             p::relative_error_bound = relative_error_bound_,
             p::p_rng = p_rng_,
             p::null_subdomain_index =
                     create_null_subdomain_index(null_subdomain_index_),
             p::construct_surface_patch_index =
                     create_construct_surface_patch_index(construct_surface_patch_index_));
  }
/// @}
#ifndef DOXYGEN_RUNNING
  template<typename CGAL_NP_TEMPLATE_PARAMETERS>
  static Labeled_mesh_domain_3 create_implicit_mesh_domain(const CGAL_NP_CLASS& np = parameters::default_values())
  {
    static_assert(!parameters::is_default_parameter<CGAL_NP_CLASS, internal_np::function_param_t>::value, "Value for required parameter not found");
    static_assert(!parameters::is_default_parameter<CGAL_NP_CLASS, internal_np::bounding_object_param_t>::value, "Value for required parameter not found");

    using parameters::get_parameter;
    return create_implicit_mesh_domain(parameters::get_parameter(np, internal_np::function_param),
                                       parameters::get_parameter(np, internal_np::bounding_object_param),
                                       np);
  }

  // Overload handling parameters passed with operator=
  template<typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_1,
           typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_2,
           typename ... NP>
  static Labeled_mesh_domain_3 create_implicit_mesh_domain(const CGAL_NP_CLASS_1&  np1,
                                                           const CGAL_NP_CLASS_2&  np2,
                                                           const NP& ... nps)
  {
    return create_implicit_mesh_domain(internal_np::combine_named_parameters(np1, np2, nps...));
  }


#ifndef CGAL_NO_DEPRECATED_CODE
  template<typename SubdomainIndex = Null_functor,
           typename NullSubdomainIndex = Null_functor,
           typename ConstructSurfacePatchIndex = Null_functor>
  CGAL_DEPRECATED
  static Labeled_mesh_domain_3
  create_gray_image_mesh_domain(const CGAL::Image_3& image_,
                                double iso_value,
                                double value_outside=0,
                                double relative_error_bound = 1e-3,
                                CGAL::Random* rng = nullptr,
                                SubdomainIndex image_values_to_subdom_indices = SubdomainIndex(),
                                NullSubdomainIndex null_subdomain_index_ = NullSubdomainIndex(),
                                ConstructSurfacePatchIndex construct_surface_patch_index_ = ConstructSurfacePatchIndex())
  {
    return create_gray_image_mesh_domain(image_, parameters::iso_value(iso_value)
                                                            .image_values_to_subdomain_indices(image_values_to_subdom_indices)
                                                            .value_outside(value_outside)
                                                            .relative_error_bound(relative_error_bound)
                                                            .p_rng(rng).null_subdomain_index(null_subdomain_index_)
                                                            .construct_surface_patch_index(construct_surface_patch_index_));
  }

  template<typename SubdomainIndex = Null_functor,
           typename NullSubdomainIndex = Null_functor,
           typename ConstructSurfacePatchIndex = Null_functor>
  CGAL_DEPRECATED
  static Labeled_mesh_domain_3
  create_labeled_image_mesh_domain(const CGAL::Image_3& image_,
                                   double relative_error_bound,
                                   const CGAL::Image_3& weights_ = CGAL::Image_3(),
                                   int value_outside=0,
                                   CGAL::Random* rng = nullptr,
                                   SubdomainIndex image_values_to_subdom_indices = SubdomainIndex(),
                                   NullSubdomainIndex null_subdomain_index_ = NullSubdomainIndex(),
                                   ConstructSurfacePatchIndex construct_surface_patch_index_ = ConstructSurfacePatchIndex())
  {
    return create_labeled_image_mesh_domain(image_, parameters::weights(weights_)
                                                               .image_values_to_subdomain_indices(image_values_to_subdom_indices)
                                                               .value_outside(value_outside)
                                                               .relative_error_bound(relative_error_bound)
                                                               .p_rng(rng)
                                                               .null_subdomain_index(null_subdomain_index_)
                                                               .construct_surface_patch_index(construct_surface_patch_index_));
  }
#endif


  /*
   * Constructs  a set of `n` points on the surface, and output them to
   *  the output iterator `pts` whose value type is required to be
   *  `std::pair<Points_3, Index>`.
   */
  struct Construct_initial_points
  {
    Construct_initial_points(const Labeled_mesh_domain_3& domain)
      : r_domain_(domain) {}

    template<class OutputIterator>
    OutputIterator operator()(OutputIterator pts, const int n = 12) const;

  private:
    const Labeled_mesh_domain_3& r_domain_;
  };

  // Returns Construct_initial_points object
  Construct_initial_points construct_initial_points_object() const
  {
    return Construct_initial_points(*this);
  }

  /*
   * Returns a bounding box of the domain
   */
  Bbox_3 bbox() const {
    return this->bbox_.bbox();
  }

  /*
   * Returns true if point `p` is in the domain. If `p` is in the
   *  domain, the parameter index is set to the index of the subdomain
   *  including $p$. It is set to the default value otherwise.
   */
  struct Is_in_domain
  {
    Is_in_domain(const Labeled_mesh_domain_3& domain) : r_domain_(domain) {}

    Subdomain operator()(const Point_3& p) const
    {
      // null(f(p)) means p is outside the domain
      Subdomain_index index = (r_domain_.function_)(p);
      if ( r_domain_.null(index) )
        return Subdomain{};
      else
        return Subdomain{ index };
    }
  private:
    const Labeled_mesh_domain_3& r_domain_;
  };

  // Returns Is_in_domain object
  Is_in_domain is_in_domain_object() const { return Is_in_domain(*this); }

  /*
   * Returns `true` if the element `type` intersect properly any of the
   * surface patches describing the either the domain boundary or some
   * subdomain boundary.
   * `Type` is either `Segment_3`, `Ray_3` or `Line_3`.
   * Parameter index is set to the index of the intersected surface patch
   * if `true` is returned and to the default `Surface_patch_index`
   * value otherwise.
   */
  struct Do_intersect_surface
  {
    Do_intersect_surface(const Labeled_mesh_domain_3& domain)
      : r_domain_(domain) {}

    Surface_patch operator()(const Segment_3& s) const
    {
      return this->operator()(s.source(), s.target());
    }

    Surface_patch operator()(const Ray_3& r) const
    {
      return clip_to_segment(r);
    }

    Surface_patch operator()(const Line_3& l) const
    {
      return clip_to_segment(l);
    }

  private:
    // Returns true if points `a` and `b` do not belong to the same subdomain
    // `index` is set to the surface index of subdomains f(a), f(b)
    Surface_patch operator()(const Point_3& a, const Point_3& b) const
    {
      // If f(a) != f(b), then [a,b] intersects some surface. Here we consider
      // [a,b] intersects surface_patch labeled <f(a),f(b)> (or <f(b),f(a)>).
      // It may be false, further rafinement will improve precision
      const Subdomain_index value_a = r_domain_.function_(a);
      const Subdomain_index value_b = r_domain_.function_(b);

      if ( value_a != value_b ) {
        if( r_domain_.null(value_a) && r_domain_.null(value_b) )
          return Surface_patch();
        else
          return Surface_patch(r_domain_.make_surface_index(value_a, value_b));
      }
      else
        return Surface_patch();
    }

    /*
     * Clips  `query` to a segment `s`, and call `operator()(s)`
     */
    template<typename Query>
    Surface_patch clip_to_segment(const Query& query) const
    {
      const auto clipped = CGAL::intersection(query, r_domain_.bbox_);
      if(clipped)
        if(const Segment_3* s = std::get_if<Segment_3>(&*clipped))
          return this->operator()(*s);

      return Surface_patch();
    }

  private:
    const Labeled_mesh_domain_3& r_domain_;
  };

  // Returns Do_intersect_surface object
  Do_intersect_surface do_intersect_surface_object() const
  {
    return Do_intersect_surface(*this);
  }

  /*
   * Returns a point in the intersection of the primitive `type`
   * with some boundary surface.
   * `Type1` is either `Segment_3`, `Ray_3` or `Line_3`.
   * The integer `dimension` is set to the dimension of the lowest
   * dimensional face in the input complex containing the returned point, and
   * `index` is set to the index to be stored at a mesh vertex lying
   * on this face.
   */
  struct Construct_intersection
  {
    Construct_intersection(const Labeled_mesh_domain_3& domain)
      : r_domain_(domain) {}

    Intersection operator()(const Segment_3& s) const
    {
#ifndef CGAL_MESH_3_NO_LONGER_CALLS_DO_INTERSECT_3
      CGAL_precondition(r_domain_.do_intersect_surface_object()(s) != std::nullopt);
#endif // NOT CGAL_MESH_3_NO_LONGER_CALLS_DO_INTERSECT_3
      return this->operator()(s.source(),s.target());
    }

    Intersection operator()(const Ray_3& r) const
    {
      return clip_to_segment(r);
    }

    Intersection operator()(const Line_3& l) const
    {
      return clip_to_segment(l);
    }

  private:
    /*
     * Returns a point in the intersection of [a,b] with the surface
     *  `a` must be the source point, and `b` the out point. It's important
     * because it drives bisection cuts.
     * Indeed, the returned point is the first intersection from  `[a,b]`
     * with a subdomain surface.
     */
    Intersection operator()(const Point_3& a, const Point_3& b) const
    {
      // Functors
      typename BGT::Compute_squared_distance_3 squared_distance =
                                      BGT().compute_squared_distance_3_object();
      typename BGT::Construct_midpoint_3 midpoint =
                                      BGT().construct_midpoint_3_object();

      // Non const points
      Point_3 p1 = a;
      Point_3 p2 = b;
      Point_3 mid = midpoint(p1, p2);

      // Cannot be const: those values are modified below.
      Subdomain_index value_at_p1 = r_domain_.function_(p1);
      Subdomain_index value_at_p2 = r_domain_.function_(p2);
      Subdomain_index value_at_mid = r_domain_.function_(mid);

      // If both extremities are in the same subdomain,
      // there is no intersection.
      // This should not happen...
      if( value_at_p1 == value_at_p2 )
      {
        return Intersection();
      }
      if( r_domain_.null(value_at_p1) && r_domain_.null(value_at_p2) ) {
        return Intersection();
      }

      // Else lets find a point (by bisection)
      // Bisection ends when the point is near than error bound from surface
      while(true)
      {
        // If the two points are enough close, then we return midpoint
        if ( squared_distance(p1, p2) < r_domain_.squared_error_bound_ )
        {
          CGAL_assertion(value_at_p1 != value_at_p2 &&
             ! ( r_domain_.null(value_at_p1) && r_domain_.null(value_at_p2) ) );
          const Surface_patch_index sp_index =
            r_domain_.make_surface_index(value_at_p1, value_at_p2);
          const Index index = r_domain_.index_from_surface_patch_index(sp_index);
          return Intersection(mid, index, 2);
        }

        // Else we must go on
        // Here we consider that p1(a) is the source point. Thus, we keep p1 and
        // change p2 if f(p1)!=f(p2).
        // That allows us to find the first intersection from a of [a,b] with
        // a surface.
        if ( value_at_p1 != value_at_mid &&
             ! ( r_domain_.null(value_at_p1) && r_domain_.null(value_at_mid) ) )
        {
          p2 = mid;
          value_at_p2 = value_at_mid;
        }
        else
        {
          p1 = mid;
          value_at_p1 = value_at_mid;
        }

        mid = midpoint(p1, p2);
        value_at_mid = r_domain_.function_(mid);
      }
    }

    // Clips  `query` to a segment `s`, and call `operator()(s)`
    template<typename Query>
    Intersection clip_to_segment(const Query& query) const
    {
      const auto clipped = CGAL::intersection(query, r_domain_.bbox_);
      if(clipped)
        if(const Segment_3* s = std::get_if<Segment_3>(&*clipped))
          return this->operator()(*s);

      return Intersection();
    }

  private:
    const Labeled_mesh_domain_3& r_domain_;
  };

  // Returns Construct_intersection object
  Construct_intersection construct_intersection_object() const
  {
    return Construct_intersection(*this);
  }

  /*
   * Returns the index to be stored in a vertex lying on the surface identified
   * by `index`.
   */
  Index index_from_surface_patch_index(const Surface_patch_index& index) const
  { return Index(index); }

  /*
   * Returns the index to be stored in a vertex lying in the subdomain
   * identified by `index`.
   */
  Index index_from_subdomain_index(const Subdomain_index& index) const
  { return Index(index); }

  /*
   * Returns the `Surface_patch_index` of the surface patch
   * where lies a vertex with dimension 2 and index `index`.
   */
  Surface_patch_index surface_patch_index(const Index& index) const
  { return Mesh_3::internal::get_index<Surface_patch_index>(index); }

  /*
   * Returns the index of the subdomain containing a vertex
   *  with dimension 3 and index `index`.
   */
  Subdomain_index subdomain_index(const Index& index) const
  { return Mesh_3::internal::get_index<Subdomain_index>(index); }

  // -----------------------------------
  // Backward Compatibility
  // -----------------------------------
#ifndef CGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX
  typedef Surface_patch_index   Surface_index;

  Index index_from_surface_index(const Surface_index& index) const
  { return index_from_surface_patch_index(index); }

  Surface_index surface_index(const Index& index) const
  { return surface_patch_index(index); }
#endif // CGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX
  // -----------------------------------
  // End backward Compatibility
  // -----------------------------------

protected:
  // Returns Surface_patch_index from `i` and `j`
  Surface_patch_index make_surface_index(const Subdomain_index i,
                                         const Subdomain_index j) const
  {
    if(i < j)
      return this->cstr_s_p_index(i, j);
    else
      return this->cstr_s_p_index(j, i);
  }

  // Returns the bounding sphere of an Iso_cuboid_3
  Sphere_3 bounding_sphere(const Iso_cuboid_3& bbox) const
  {
    typename BGT::Construct_sphere_3 sphere = BGT().construct_sphere_3_object();
    return sphere((bbox.min)(), (bbox.max)());
  }

  template <typename Image_word_type,
            typename FT, typename FT2, typename Functor>
  static
  Function
  create_gray_image_wrapper_with_known_word_type
  (const CGAL::Image_3& image,
   const FT& iso_value,
   const Functor& image_values_to_subdomain_indices,
   const FT2& value_outside)
  {
    using Mesh_3::internal::Create_gray_image_values_to_subdomain_indices;
    typedef Create_gray_image_values_to_subdomain_indices<Functor> C_i_v_t_s_i;
    typedef typename C_i_v_t_s_i::type Image_values_to_subdomain_indices;
    Image_values_to_subdomain_indices transform_fct =
      C_i_v_t_s_i()(image_values_to_subdomain_indices, iso_value);

    typedef Mesh_3::Image_to_labeled_function_wrapper<Image_word_type,
                                                      double,
                                                      Subdomain_index,
                                                      false>           Wrapper;
    return Wrapper(image,
                   transform_fct,
                   value_outside) ;
  }

  template <typename FT, typename FT2, typename Functor>
  static
  Function
  create_gray_image_wrapper(const CGAL::Image_3& image,
                            const FT& iso_value,
                            const Functor& image_values_to_subdomain_indices,
                            const FT2& value_outside)
  {
    CGAL_IMAGE_IO_CASE(image.image(),
       return create_gray_image_wrapper_with_known_word_type<Word>
                       (image,
                        iso_value,
                        image_values_to_subdomain_indices,
                        value_outside);
                       );
    CGAL_error_msg("This place should never be reached, because it would mean "
                   "the image word type is a type that is not handled by "
                   "CGAL_ImageIO.");
    return Function();
  }

  template <typename Image_word_type,
            typename FT, typename Functor>
  static
  Function
  create_labeled_image_wrapper_with_known_word_type
  (const CGAL::Image_3& image,
   const Functor& image_values_to_subdomain_indices,
   const FT& value_outside)
  {
    using Mesh_3::internal::Create_labeled_image_values_to_subdomain_indices;
    typedef Create_labeled_image_values_to_subdomain_indices<Functor> C_i_v_t_s_i;
    typedef typename C_i_v_t_s_i::type Image_values_to_subdomain_indices;
    Image_values_to_subdomain_indices transform_fct =
      C_i_v_t_s_i()(image_values_to_subdomain_indices);

    typedef Mesh_3::Image_to_labeled_function_wrapper<Image_word_type,
                                                      int,
                                                      Subdomain_index> Wrapper;
    return Wrapper(image,
                   transform_fct,
                   transform_fct(value_outside));
  }

  template <typename Image_word_type,
            typename FT, typename Functor>
  static
  Function
  create_weighted_labeled_image_wrapper_with_know_word_type
  (const CGAL::Image_3& image,
   const CGAL::Image_3& weights,
   const Functor& image_values_to_subdomain_indices,
   const FT& value_outside)
  {
    using Mesh_3::internal::Create_labeled_image_values_to_subdomain_indices;
    typedef Create_labeled_image_values_to_subdomain_indices<Functor> C_i_v_t_s_i;
    typedef typename C_i_v_t_s_i::type Image_values_to_subdomain_indices;
    Image_values_to_subdomain_indices transform_fct =
      C_i_v_t_s_i()(image_values_to_subdomain_indices);

    typedef Mesh_3::Image_plus_weights_to_labeled_function_wrapper<
      Image_word_type,
      int, //interpolation type
      unsigned char, // Weights_type,
      Subdomain_index> Wrapper;
    return Wrapper(image,
                   weights,
                   transform_fct,
                   transform_fct(value_outside));
  }

  template <typename FT, typename Functor>
  static
  Function
  create_labeled_image_wrapper(const CGAL::Image_3& image,
                               const Functor& image_values_to_subdomain_indices,
                               const FT& value_outside)
  {
    CGAL_IMAGE_IO_CASE(image.image(),
       return create_labeled_image_wrapper_with_known_word_type<Word>
                       (image,
                        image_values_to_subdomain_indices,
                        value_outside);
                       );
    CGAL_error_msg("This place should never be reached, because it would mean "
                   "the image word type is a type that is not handled by "
                   "CGAL_ImageIO.");
    return Function();
  }

  template <typename FT, typename Functor>
  static
  Function
  create_weighted_labeled_image_wrapper(const CGAL::Image_3& image,
                                        const CGAL::Image_3& weights,
                                        const Functor& image_values_to_subdomain_indices,
                                        const FT& value_outside)
  {
    CGAL_IMAGE_IO_CASE(image.image(),
      return create_weighted_labeled_image_wrapper_with_know_word_type<Word>
                        (image,
                         weights,
                         image_values_to_subdomain_indices,
                         value_outside);
                        );
    CGAL_error_msg("This place should never be reached, because it would mean "
      "the image word type is a type that is not handled by "
      "CGAL_ImageIO.");
    return Function();
  }

  static
  Construct_surface_patch_index
  create_construct_surface_patch_index(const Null_functor&) {
    return Impl_details::construct_pair_functor();
  }

  template <typename Functor>
  static
  Construct_surface_patch_index
  create_construct_surface_patch_index(const Functor& functor) {
    return functor;
  }

  static Null create_null_subdomain_index(const Null_functor&) {
    return Null_subdomain_index();
  }

  template <typename Functor>
  static Null create_null_subdomain_index(const Functor& functor) {
    return functor;
  }

public:
  // Returns bounding box
  const Iso_cuboid_3& bounding_box() const { return this->bbox_; }
#endif //DOXYGEN_RUNNING

};  // end class Labeled_mesh_domain_3

#ifndef DOXYGEN_RUNNING
//-------------------------------------------------------
// Method implementation
//-------------------------------------------------------
template<class BGT, class Subdomain_index, class Surface_patch_index>
template<class OutputIterator>
OutputIterator
Labeled_mesh_domain_3<BGT, Subdomain_index, Surface_patch_index>::
Construct_initial_points::operator()(OutputIterator pts,
                                     const int nb_points) const
{
  // Create point_iterator on and in bounding_sphere
  typedef Random_points_on_sphere_3<Point_3> Random_points_on_sphere_3;
  typedef Random_points_in_sphere_3<Point_3> Random_points_in_sphere_3;

  const FT squared_radius = BGT().compute_squared_radius_3_object()(
      r_domain_.bounding_sphere(r_domain_.bbox_));

  const double radius = std::sqrt(CGAL::to_double(squared_radius));

  CGAL::Random& rng = *(r_domain_.p_rng_);
  Random_points_on_sphere_3 random_point_on_sphere(radius, rng);
  Random_points_in_sphere_3 random_point_in_sphere(radius, rng);

  // Get some functors
  typename BGT::Construct_segment_3 segment_3 =
                              BGT().construct_segment_3_object();
  typename BGT::Construct_vector_3 vector_3 =
                              BGT().construct_vector_3_object();
  typename BGT::Construct_translated_point_3 translate =
                              BGT().construct_translated_point_3_object();
  typename BGT::Construct_center_3 center = BGT().construct_center_3_object();

  // Get translation from origin to sphere center
  Point_3 center_pt = center(r_domain_.bounding_sphere(r_domain_.bbox_));
  const Vector_3 sphere_translation = vector_3(CGAL::ORIGIN, center_pt);

  // Create nb_point points
  int n = nb_points;
#ifdef CGAL_MESH_3_VERBOSE
  std::cerr << "construct initial points (nb_points: " << nb_points << ")\n";
#endif
  while ( 0 != n )
  {
    // Get a random segment
    const Point_3 random_point = translate(*random_point_on_sphere,
                                           sphere_translation);
    const Segment_3 random_seg = segment_3(center_pt, random_point);

    // Add the intersection to the output if it exists
    Surface_patch surface = r_domain_.do_intersect_surface_object()(random_seg);
    if ( surface )
    {
      const Point_3 intersect_pt = std::get<0>(
          r_domain_.construct_intersection_object()(random_seg));
      *pts++ = std::make_pair(intersect_pt,
                              r_domain_.index_from_surface_patch_index(*surface));
      --n;

#ifdef CGAL_MESH_3_VERBOSE
      std::cerr << boost::format("\r             \r"
                                 "%1%/%2% initial point(s) found...")
                   % (nb_points - n)
                   % nb_points;
#endif
    }
    else
    {
      // Get a new random point into sphere as center of object
      // It may be necessary if the center of the domain is empty, e.g. torus
      // In general case, it is good for input point dispersion
      ++random_point_in_sphere;
      center_pt = translate(*random_point_in_sphere, sphere_translation);
    }
    ++random_point_on_sphere;
  }

#ifdef CGAL_MESH_3_VERBOSE
  std::cerr << "\n";
#endif
  return pts;
}
#endif

}  // end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_LABELED_MESH_DOMAIN_3_H
