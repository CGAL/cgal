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

#include <CGAL/Mesh_3/internal/Handle_IO_for_pair_of_int.h>
#include <CGAL/Mesh_3/internal/indices_management.h>

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
#include <boost/optional.hpp>

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

  /// Returns a box enclosing image \c im
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

} // end namespace CGAL::Mesh_3::internal
} // end namespace CGAL::Mesh_3

struct Null_subdomain_index {
  template <typename T>
  bool operator()(const T& x) const { return 0 == x; }
};

template <typename Subdomain_index>
struct Construct_pair_from_subdomain_indices {
  typedef std::pair<Subdomain_index, Subdomain_index> result_type;

  result_type operator()(Subdomain_index a, Subdomain_index b) const {
    return result_type(a, b);
  }
}; // end class template Construct_pair_from_subdomain_indices

template <typename Geom_traits,
          typename Subdomain_index,
          typename Surface_patch_index_>
class Labeled_mesh_domain_3_impl_details
{
protected:
  typedef Surface_patch_index_ Surface_patch_index;
  typedef typename Geom_traits::Point_3 Point_3;
  typedef typename Geom_traits::Sphere_3 Sphere_3;
  typedef typename Geom_traits::Iso_cuboid_3 Iso_cuboid_3;
  typedef typename Geom_traits::FT FT;
  typedef std::shared_ptr<CGAL::Random> CGAL_Random_share_ptr_t;
  /// Returns squared error bound from \c bbox and \c error
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
  Labeled_mesh_domain_3_impl_details(const Function& f,
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

  /// The function which answers subdomain queries
  typedef std::function<Subdomain_index(const Point_3&)> Function;
  Function function_;
  /// The bounding box
  const Iso_cuboid_3 bbox_;

  typedef std::function<
    Surface_patch_index(Subdomain_index,
                        Subdomain_index)> Construct_surface_patch_index;
  Construct_surface_patch_index cstr_s_p_index;
  /// The functor that decides which sub-domain indices correspond to the
  /// outside of the domain.
  typedef std::function<bool(Subdomain_index)> Null;
  Null null;
  /// The random number generator used by Construct_initial_points
  CGAL_Random_share_ptr_t p_rng_;
  /// Error bound relative to sphere radius
  FT squared_error_bound_;
}; // Labeled_mesh_domain_3_impl_details

/**
 * \class Labeled_mesh_domain_3
 *
 * Function f must take his values into N.
 * Let p be a Point.
 *  - f(p)=0 means that p is outside domain.
 *  - f(p)=a, a!=0 means that p is inside subdomain a.
 *
 *  Any boundary facet is labelled <a,b>, with a<b, where a and b are the
 *  tags of it's incident subdomain.
 *  Thus, a boundary facet of the domain is labelled <0,b>, where b!=0.
 */
template<class BGT,
         class Subdomain_index_ = int,
         class Surface_patch_index_ = std::pair<Subdomain_index_,
                                                Subdomain_index_> >
class Labeled_mesh_domain_3 :
    protected Labeled_mesh_domain_3_impl_details<BGT,
                                                 Subdomain_index_,
                                                 Surface_patch_index_>
{
public:
  //-------------------------------------------------------
  // Index Types
  //-------------------------------------------------------
  /// Type of indexes for cells of the input complex
  typedef Subdomain_index_                  Subdomain_index;
  typedef boost::optional<Subdomain_index>  Subdomain;

  /// Type of indexes for cells of the input complex
  typedef Surface_patch_index_                  Surface_patch_index;
  typedef boost::optional<Surface_patch_index>  Surface_patch;

  /// Type of indexes to characterize the lowest dimensional face of the input
  /// complex on which a vertex lie
  typedef typename CGAL::Mesh_3::internal::
    Index_generator<Subdomain_index, Surface_patch_index>::Index Index;

private:
  typedef Labeled_mesh_domain_3_impl_details<BGT,
                                             Subdomain_index,
                                             Surface_patch_index
                                             > Impl_details;
  typedef typename Impl_details::Null     Null;
  typedef typename Impl_details::Construct_surface_patch_index
                                          Construct_surface_patch_index;
  typedef typename Impl_details::Function Function;

public:
  /// Geometric object types
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

  template<typename Func, typename BoundingBox, typename CGAL_NP_TEMPLATE_PARAMETERS>
  Labeled_mesh_domain_3(const Func& f, const BoundingBox& bbox, const CGAL_NP_CLASS& np = parameters::default_values())
  :Impl_details(f,bbox,
                parameters::choose_parameter(parameters::get_parameter(np, internal_np::error_bound), FT(1e-3)),
                parameters::choose_parameter(parameters::get_parameter(np, internal_np::surface_patch_index), construct_pair_functor()),
                parameters::choose_parameter(parameters::get_parameter(np, internal_np::null_subdomain_index_param), Null_subdomain_index()),
                parameters::choose_parameter(parameters::get_parameter(np, internal_np::rng), nullptr))
  {
  }
#ifndef CGAL_NO_DEPRECATED_CODE
        template<typename Func, typename BoundingBox, typename ... NP_PACK>
        Labeled_mesh_domain_3(const Func& f, const BoundingBox& bbox, const NP_PACK& ...nps):Labeled_mesh_domain_3(f, bbox, internal_np::combine_named_parameters(nps...))
        {
        }
#endif //CGAL_NO_DEPRECATED_CODE

  /**
   * Backward-compatibility constructors, with `null_subdomain_index` as
   * fourth parameter.
   * @{
   */
  Labeled_mesh_domain_3(const Function& f,
                        const Sphere_3& bounding_sphere,
                        const FT& error_bound = FT(1e-3),
                        Null null = Null_subdomain_index(),
                        CGAL::Random* p_rng = nullptr)
    : Impl_details(f, bounding_sphere,
                   error_bound,
                   construct_pair_functor(),
                   null, p_rng) {}

  Labeled_mesh_domain_3(const Function& f,
                        const Bbox_3& bbox,
                        const FT& error_bound = FT(1e-3),
                        Null null = Null_subdomain_index(),
                        CGAL::Random* p_rng = nullptr)
    : Impl_details(f, bbox,
                   error_bound,
                   construct_pair_functor(),
                   null, p_rng) {}

  Labeled_mesh_domain_3(const Function& f,
                        const Iso_cuboid_3& bbox,
                        const FT& error_bound = FT(1e-3),
                        Null null = Null_subdomain_index(),
                        CGAL::Random* p_rng = nullptr)
    : Impl_details(f, bbox, error_bound,
                   construct_pair_functor(),
                   null, p_rng)
  {}
  /**
   * @}
   */

  /// Named constructors
  /// @{
        /*!
      \brief Construction from a 3D gray image

      This static method is a <em>named constructor</em>. It constructs a domain
      described by a 3D gray image. A 3D gray image is a grid of voxels,
      where each voxel is associated with a gray level value.  Unless otherwise specified by the parameter `image_values_to_subdom_indices`, the domain to
      be discretized is the union of voxels that lie inside a surface
      described by an isolevel value, called \a isovalue. The voxels lying
      inside the domain have gray level values that are larger than the
      isovalue.

      The value of voxels is interpolated to a gray level value at any query point.

      This constructor uses named parameters . They can be specified in any order.
         \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

         \param image the input 3D image. Must be a `CGAL::Image_3` object.

         \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below:

         \cgalNamedParamsBegin
           \cgalParamNBegin{iso_value}
             \cgalParamDescription{the isovalue, inside
                                  `image`, of the surface describing the boundary of the object to be
                                   meshed.}
             \cgalParamDefault{0}

           \cgalParamNBegin{image_values_to_subdomain_indices}
             \cgalParamDescription{a function or
                                   a function object, compatible with the signature
                                   `Subdomain_index(double)`. This function returns the subdomain index
                                   corresponding to a pixel value. If this parameter is used, then the
                                   parameter `iso_value` is ignored..}
             \cgalParamDefault{Null_functor()}

           \cgalParamNBegin{value_outside}
             \cgalParamDescription{the value attached to voxels
                                   outside of the domain to be meshed. It should be lower than
                                  `iso_value`.}
             \cgalParamDefault{0}

           \cgalParamNBegin{relative_error_bound}
             \cgalParamDescription{ is the relative error
                                    bound, relative to the diameter of the box of the image.}
             \cgalParamDefault{FT(1e-3)}


         \cgalNamedParamsEnd

      \cgalHeading{Examples}

      From the example (\ref Mesh_3/mesh_3D_gray_image.cpp), where the name
      of the parameters is not specified, as they are given is the same
      order as the parameters definition:

      \snippet Mesh_3/mesh_3D_gray_image.cpp Domain creation

      From the example (\ref Mesh_3/mesh_3D_gray_vtk_image.cpp):

      \snippet Mesh_3/mesh_3D_gray_vtk_image.cpp Domain creation

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
    CGAL::Random* p_rng_ = choose_parameter(get_parameter(np, internal_np::rng), (CGAL::Random*)(0));
    auto null_subdomain_index_ = choose_parameter(get_parameter(np, internal_np::null_subdomain_index_param), Null_functor());
    auto construct_surface_patch_index_ = choose_parameter(get_parameter(np, internal_np::surface_patch_index), Null_functor());
    namespace p = CGAL::parameters;
    return Labeled_mesh_domain_3
              (create_gray_image_wrapper
                       (image_,
                        iso_value_,
                        image_values_to_subdomain_indices_,
                        value_outside_),
               Mesh_3::internal::compute_bounding_box(image_),
               p::relative_error_bound = relative_error_bound_,
               p::p_rng = p_rng_,
               p::null_subdomain_index =
                       create_null_subdomain_index(null_subdomain_index_),
               p::construct_surface_patch_index =
                       create_construct_surface_patch_index(construct_surface_patch_index_));

  }
#ifndef CGAL_NO_DEPRECATED_CODE
        template<typename... NP_Pack>
        static Labeled_mesh_domain_3 create_gray_image_mesh_domain(const CGAL::Image_3& image_, const NP_Pack& ...nps)
        {
            return create_gray_image_mesh_domain(image_, internal_np::combine_named_parameters(nps...));
        }
#endif //CGAL_NO_DEPRECATED_CODE


        /*!
         * \brief Construction from a 3D labeled image

        This static method is a <em>named constructor</em>. It constructs a
        domain described by a 3D labeled image. A 3D labeled image is a grid
        of voxels, where each voxel is associated with an index (a subdomain
        index) characterizing the subdomain in which the voxel lies. The
        domain to be discretized is the union of voxels that have non-zero
        values.

        This constructor uses named parameters . They can be specified in any order.
        \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

        \param image the input 3D image. Must be a `CGAL::Image_3` object.

        \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below:

          \cgalNamedParamsBegin
           \cgalParamNBegin{weights}
             \cgalParamDescription{an input 3D image that provides
                                   weights associated to each voxel (the word type is `unsigned char`,
                                   and the voxels values are integers between 0 and 255).
                                   The weights image can be generated with `CGAL::Mesh_3::generate_label_weights()`.
                                   Its dimensions must be the same as the dimensions of `parameters::image`.}
             \cgalParamDefault{CGAL::Image_3()}

           \cgalParamNBegin{value_outside}
             \cgalParamDescription{the value attached to voxels
                                   outside of the domain to be meshed. It should be lower than
                                  `iso_value`.}
             \cgalParamDefault{0}

           \cgalParamNBegin{relative_error_bound}
             \cgalParamDescription{ is the relative error
                                    bound, relative to the diameter of the box of the image.}
             \cgalParamDefault{FT(1e-3)}
         \cgalNamedParamsEnd

\cgalHeading{Example}

From the example (\ref Mesh_3/mesh_3D_image.cpp):

\snippet Mesh_3/mesh_3D_image.cpp Domain creation

From the example (\ref Mesh_3/mesh_3D_weighted_image.cpp),
where the labeled image is used with a precomputed 3D image of weights :

\snippet Mesh_3/mesh_3D_weighted_image.cpp Domain creation

 */
        template<typename CGAL_NP_TEMPLATE_PARAMETERS>
        static Labeled_mesh_domain_3 create_labeled_image_mesh_domain(const CGAL::Image_3& image_, const CGAL_NP_CLASS& np = parameters::default_values())
        {
            using parameters::get_parameter;
            using parameters::choose_parameter;
            auto iso_value_ = choose_parameter(get_parameter(np, internal_np::iso_value_param), FT(0));
            auto value_outside_ = choose_parameter(get_parameter(np, internal_np::voxel_value), FT(0));
            FT relative_error_bound_ = choose_parameter(get_parameter(np, internal_np::error_bound), FT(1e-3));
            auto image_values_to_subdomain_indices_ = choose_parameter(get_parameter(np, internal_np::image_subdomain_index), Null_functor());
            CGAL::Random* p_rng_ = choose_parameter(get_parameter(np, internal_np::rng), (CGAL::Random*)(0));
            auto null_subdomain_index_ = choose_parameter(get_parameter(np, internal_np::null_subdomain_index_param), Null_functor());
            auto construct_surface_patch_index_ = choose_parameter(get_parameter(np, internal_np::surface_patch_index), Null_functor());
            CGAL::Image_3 weights_ = choose_parameter(get_parameter(np, internal_np::weights_param), CGAL::Image_3());

            namespace p = CGAL::parameters;
            if (weights_.is_valid())
            {
                return Labeled_mesh_domain_3
                        (create_weighted_labeled_image_wrapper
                                 (image_,
                                  weights_,
                                  image_values_to_subdomain_indices_,
                                  value_outside_),
                         Mesh_3::internal::compute_bounding_box(image_),
                         p::relative_error_bound = relative_error_bound_,
                         p::p_rng = p_rng_,
                         p::null_subdomain_index =
                                 create_null_subdomain_index(null_subdomain_index_),
                         p::construct_surface_patch_index =
                                 create_construct_surface_patch_index(construct_surface_patch_index_));
            }
            else
            {
                return Labeled_mesh_domain_3
                        (create_labeled_image_wrapper
                                 (image_,
                                  image_values_to_subdomain_indices_,
                                  value_outside_),
                         Mesh_3::internal::compute_bounding_box(image_),
                         p::relative_error_bound = relative_error_bound_,
                         p::p_rng = p_rng_,
                         p::null_subdomain_index =
                                 create_null_subdomain_index(null_subdomain_index_),
                         p::construct_surface_patch_index =
                                 create_construct_surface_patch_index(construct_surface_patch_index_));
            }
        }
#ifndef CGAL_NO_DEPRECATED_CODE
        template<typename... NP_Pack>
        static Labeled_mesh_domain_3 create_labeled_image_mesh_domain(const CGAL::Image_3& image_, const NP_Pack& ...nps)
        {
            return create_labeled_image_mesh_domain(image_, internal_np::combine_named_parameters(nps...));
        }
#endif //CGAL_NO_DEPRECATED_CODE



/// \name Creation of domains from implicit functions

/*!
\brief Construction from an implicit function

This static method is a <em>named constructor</em>. It constructs a domain
whose bounding surface is described implicitly as the zero level set of a
function.  The domain to be discretized is assumed to be the domain where
the function has negative values.

The method takes as argument a bounding sphere which is required to
circumscribe the surface and to have its center inside the domain.

This constructor uses named parameters (from the <em>Boost Parameter
Library</em>). They can be specified in any order.

        \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

         \param function  the implicit function, compatible with the signature `FT(Point_3)`: it takes a point as argument,
                          and returns a scalar value. That object must be model of `CopyConstructible`.
         \param bounding_object the bounding object is either a bounding sphere (of type `Sphere_3`), a bounding box (type
                                `Bbox_3`), or a bounding `Iso_cuboid_3`. It must bounds the surface, and
                                its center must be inside the domain.

         \param np an optional sequence of \ref bgl_namedparameters "Named Parameters".

\cgalHeading{Examples}

From the example (\ref Mesh_3/mesh_implicit_sphere.cpp), where the name of
the parameters is not specified, as they are given is the same order as the
parameters definition:

\snippet Mesh_3/mesh_implicit_sphere.cpp Domain creation

From the example (\ref Mesh_3/mesh_implicit_sphere_variable_size.cpp):

\snippet Mesh_3/mesh_implicit_sphere_variable_size.cpp Domain creation

 */

  template<typename Function, typename BoundingObject, typename CGAL_NP_TEMPLATE_PARAMETERS>
  static Labeled_mesh_domain_3 create_implicit_mesh_domain(Function function_, BoundingObject bounding_object_, const CGAL_NP_CLASS& np = parameters::default_values())
  {
      using parameters::get_parameter;
      using parameters::choose_parameter;
      FT relative_error_bound_ = choose_parameter(get_parameter(np, internal_np::error_bound), FT(1e-3));
      CGAL::Random* p_rng_ = choose_parameter(get_parameter(np, internal_np::rng), (CGAL::Random*)(0));
      auto null_subdomain_index_ = choose_parameter(get_parameter(np, internal_np::null_subdomain_index_param), Null_functor());
      auto construct_surface_patch_index_ = choose_parameter(get_parameter(np, internal_np::surface_patch_index), Null_functor());
      namespace p = CGAL::parameters;
      return Labeled_mesh_domain_3
              (make_implicit_to_labeling_function_wrapper<BGT>(function_),
               bounding_object_,
               p::relative_error_bound = relative_error_bound_,
               p::p_rng = p_rng_,
               p::null_subdomain_index =
                       create_null_subdomain_index(null_subdomain_index_),
               p::construct_surface_patch_index =
                       create_construct_surface_patch_index(construct_surface_patch_index_));
  }

#ifndef CGAL_NO_DEPRECATED_CODE
  template<typename Function, typename BoundingObject, typename... NP_Pack>
  static Labeled_mesh_domain_3 create_implicit_mesh_domain(Function function_, BoundingObject bounding_object_, const NP_Pack& ...nps)
  {
      return create_implicit_mesh_domain(function_, bounding_object_, internal_np::combine_named_parameters(nps...));
  }
#endif //CGAL_NO_DEPRECATED_CODE

  /// @}

  /**
   * Constructs  a set of \ccc{n} points on the surface, and output them to
   *  the output iterator \ccc{pts} whose value type is required to be
   *  \ccc{std::pair<Points_3, Index>}.
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

  /// Returns Construct_initial_points object
  Construct_initial_points construct_initial_points_object() const
  {
    return Construct_initial_points(*this);
  }

  /**
   * Returns a bounding box of the domain
   */
  Bbox_3 bbox() const {
    return this->bbox_.bbox();
  }

  /**
   * Returns true if point~\ccc{p} is in the domain. If \ccc{p} is in the
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
        return Subdomain();
      else
        return Subdomain(index);
    }
  private:
    const Labeled_mesh_domain_3& r_domain_;
  };

  /// Returns Is_in_domain object
  Is_in_domain is_in_domain_object() const { return Is_in_domain(*this); }

  /**
   * Returns true is the element \ccc{type} intersect properly any of the
   * surface patches describing the either the domain boundary or some
   * subdomain boundary.
   * \ccc{Type} is either \ccc{Segment_3}, \ccc{Ray_3} or \ccc{Line_3}.
   * Parameter index is set to the index of the intersected surface patch
   * if \ccc{true} is returned and to the default \ccc{Surface_patch_index}
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
    /// Returns true if points \c a & \c b do not belong to the same subdomain
    /// \c index is set to the surface index of subdomains f(a), f(b)
    Surface_patch operator()(const Point_3& a, const Point_3& b) const
    {
      // If f(a) != f(b), then [a,b] intersects some surface. Here we consider
      // [a,b] intersects surface_patch labelled <f(a),f(b)> (or <f(b),f(a)>).
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

    /**
     * Clips \c query to a segment \c s, and call operator()(s)
     */
    template<typename Query>
    Surface_patch clip_to_segment(const Query& query) const
    {
      const auto clipped = CGAL::intersection(query, r_domain_.bbox_);
      if(clipped)
        if(const Segment_3* s = boost::get<Segment_3>(&*clipped))
          return this->operator()(*s);

      return Surface_patch();
    }

  private:
    const Labeled_mesh_domain_3& r_domain_;
  };

  /// Returns Do_intersect_surface object
  Do_intersect_surface do_intersect_surface_object() const
  {
    return Do_intersect_surface(*this);
  }

  /**
   * Returns a point in the intersection of the primitive \ccc{type}
   * with some boundary surface.
   * \ccc{Type1} is either \ccc{Segment_3}, \ccc{Ray_3} or \ccc{Line_3}.
   * The integer \ccc{dimension} is set to the dimension of the lowest
   * dimensional face in the input complex containing the returned point, and
   * \ccc{index} is set to the index to be stored at a mesh vertex lying
   * on this face.
   */
  struct Construct_intersection
  {
    Construct_intersection(const Labeled_mesh_domain_3& domain)
      : r_domain_(domain) {}

    Intersection operator()(const Segment_3& s) const
    {
#ifndef CGAL_MESH_3_NO_LONGER_CALLS_DO_INTERSECT_3
      CGAL_precondition(r_domain_.do_intersect_surface_object()(s));
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
    /**
     * Returns a point in the intersection of [a,b] with the surface
     * \c a must be the source point, and \c b the out point. It's important
     * because it drives bisection cuts.
     * Indeed, the returned point is the first intersection from \c [a,b]
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

    /// Clips \c query to a segment \c s, and call operator()(s)
    template<typename Query>
    Intersection clip_to_segment(const Query& query) const
    {
      const auto clipped = CGAL::intersection(query, r_domain_.bbox_);
      if(clipped)
        if(const Segment_3* s = boost::get<Segment_3>(&*clipped))
          return this->operator()(*s);

      return Intersection();
    }

  private:
    const Labeled_mesh_domain_3& r_domain_;
  };

  /// Returns Construct_intersection object
  Construct_intersection construct_intersection_object() const
  {
    return Construct_intersection(*this);
  }

  /**
   * Returns the index to be stored in a vertex lying on the surface identified
   * by \c index.
   */
  Index index_from_surface_patch_index(const Surface_patch_index& index) const
  { return Index(index); }

  /**
   * Returns the index to be stored in a vertex lying in the subdomain
   * identified by \c index.
   */
  Index index_from_subdomain_index(const Subdomain_index& index) const
  { return Index(index); }

  /**
   * Returns the \c Surface_patch_index of the surface patch
   * where lies a vertex with dimension 2 and index \c index.
   */
  Surface_patch_index surface_patch_index(const Index& index) const
  { return boost::get<Surface_patch_index>(index); }

  /**
   * Returns the index of the subdomain containing a vertex
   *  with dimension 3 and index \c index.
   */
  Subdomain_index subdomain_index(const Index& index) const
  { return boost::get<Subdomain_index>(index); }

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
  /// Returns Surface_patch_index from \c i and \c j
  Surface_patch_index make_surface_index(const Subdomain_index i,
                                         const Subdomain_index j) const
  {
    if(i < j)
      return this->cstr_s_p_index(i, j);
    else
      return this->cstr_s_p_index(j, i);
  }

  /// Returns the bounding sphere of an Iso_cuboid_3
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
  /// Returns bounding box
  const Iso_cuboid_3& bounding_box() const { return this->bbox_; }

};  // end class Labeled_mesh_domain_3

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


}  // end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_LABELED_MESH_DOMAIN_3_H
