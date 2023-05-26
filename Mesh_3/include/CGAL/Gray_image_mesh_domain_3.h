// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// Copyright (c) 2012  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stephane Tayeb, Laurent Rineau
//

#ifndef CGAL_GRAY_IMAGE_MESH_DOMAIN_3_H
#define CGAL_GRAY_IMAGE_MESH_DOMAIN_3_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Random.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Mesh_3/Image_to_labeled_function_wrapper.h>
#include <CGAL/Bbox_3.h>

namespace CGAL {


/*!
\ingroup PkgMesh3Domains

\deprecated The class template `Gray_image_mesh_domain_3` is deprecated
since CGAL-4.13, in favor of the class template `Labeled_mesh_domain_3` and
its static function
`Labeled_mesh_domain_3::create_gray_image_mesh_domain()`.

The class `Gray_image_mesh_domain_3` implements a domain described by a 3D
gray image. A 3D gray image is a grid of voxels,
where each voxel is associated with a gray level value.
This class is a model of the concept `MeshDomain_3`.
The domain to be discretized is the union of voxels that lie inside a surface
described by an isolevel value, called \a isovalue. The voxels lying inside the
domain have gray level values that are larger than the isovalue.

This class includes a member function that provides, by interpolation,
a gray level value at any query point.
An intersection between a segment and bounding
surfaces is detected when both segment endpoints are associated with gray level
values which are on both sides of the isovalue.
The intersection is then constructed by bisection.
The bisection stops when the query segment is shorter than a given error bound
`e`. This error bound is given by `e=d`\f$ \times\f$`bound` where `d` is the
length of the diagonal of the bounding box (in world coordinates) and
`bound` is the argument passed to the constructor of `Labeled_image_mesh_domain_3`.


\tparam Image is the type of the input image.
This parameter must be `CGAL::Image_3`.

\tparam BGT is a geometric traits class which provides
the basic operations to implement
intersection tests and intersection computations
through a bisection method. This parameter must be instantiated
with a model of the concept `BisectionGeometricTraits_3`.

\tparam Image_word_type is the data type encoded in the `Image`
input file

\cgalModels `MeshDomain_3`

\sa `BisectionGeometricTraits_3`
\sa `CGAL::make_mesh_3()`.

*/

template<class Image,
         class BGT,
         typename Image_word_type_ = float,
         typename Transform = Mesh_3::internal::Greater_than<double>,
         typename Subdomain_index = int>
class
CGAL_DEPRECATED_MSG
( "The class template `CGAL::Gray_image_mesh_domain_3` is now deprecated. "
  "Use the static member function template "
  "`Labeled_mesh_domain_3<K>::create_gray_image_mesh_domain` instead.")
Gray_image_mesh_domain_3
  : public Labeled_mesh_domain_3<BGT, Subdomain_index>
{
public:
  typedef Image_word_type_ Image_word_type;
  typedef Mesh_3::Image_to_labeled_function_wrapper<Image_word_type,
                                                    double,
                                                    Subdomain_index,
                                                    false>           Wrapper;

  typedef Labeled_mesh_domain_3<BGT, Subdomain_index>                Base;

  typedef typename Base::Sphere_3 Sphere_3;
  typedef typename Base::FT FT;
  typedef BGT Geom_traits;
  typedef CGAL::Bbox_3 Bbox_3;

  /// \name Creation
  /// @{

  /*!
    Construction from an image.
    The object to be meshed is described by the voxels that have a gray-level
    value higher than the input isovalue.
    @param image the input image
    @param iso_value the isovalue, inside `image`,
    of the surface describing the boundary of the object to be meshed.
    @param value_outside the value attached to voxels outside of the domain
    to be meshed. It should be lower than `iso_value`
    @param error_bound is relative to the size of the image.
  */
  Gray_image_mesh_domain_3(const Image& image
                           ,const Image_word_type iso_value
                           ,const Image_word_type value_outside = 0.
                           ,const FT& error_bound = FT(1e-3)
#ifndef DOXYGEN_RUNNING
                           ,CGAL::Random* p_rng = nullptr
#endif
                           )
    : Base(parameters::function = Wrapper(image,
                   Transform(iso_value),
                   Transform(iso_value)(value_outside)),
           parameters::bounding_object = Mesh_3::internal::compute_bounding_box(image),
           parameters::relative_error_bound = error_bound,
           parameters::p_rng = p_rng)
  {
    CGAL_assertion(Transform(iso_value)(value_outside) == 0);
  }

  /// @}

  Gray_image_mesh_domain_3(const Image& image,
                           const Transform& transform,
                           const Image_word_type value_outside = 0.,
                           const FT& error_bound = FT(1e-3),
                           CGAL::Random* p_rng = nullptr)
    : Base(parameters::function = Wrapper(image, transform, transform(value_outside)),
           parameters::bounding_object = Mesh_3::internal::compute_bounding_box(image),
           parameters::relative_error_bound = error_bound,
           parameters::p_rng = p_rng)
  {
    CGAL_assertion(transform(value_outside) == 0);
  }

  // Destructor
  virtual ~Gray_image_mesh_domain_3() {}
};

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_GRAY_IMAGE_MESH_DOMAIN_3_H
