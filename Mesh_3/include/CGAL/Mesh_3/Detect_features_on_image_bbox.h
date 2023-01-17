// Copyright (c) 2022 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau, Jane Tournois
//
//******************************************************************************
//
//******************************************************************************

#ifndef CGAL_MESH_3_DETECT_FEATURES_ON_IMAGE_BBOX_H
#define CGAL_MESH_3_DETECT_FEATURES_ON_IMAGE_BBOX_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/Mesh_3/polylines_to_protect.h>

#include <vector>


namespace CGAL
{
namespace Mesh_3
{
namespace internal
{
template<typename Word_type, typename Mesh_domain>
std::vector<std::vector<typename Mesh_domain::Point_3>>
detect_features_on_bbox_with_know_word_type(const CGAL::Image_3& image,
                                            Mesh_domain& domain)
{
  using Point_3 = typename Mesh_domain::Point_3;
  using Polyline_type = std::vector<Point_3>;
  using Polylines = std::vector<Polyline_type>;

  Polylines polylines_on_bbox;
  CGAL::polylines_to_protect<Point_3, Word_type>(image, polylines_on_bbox);

  return polylines_on_bbox;
}

template<typename Mesh_domain>
std::vector<std::vector<typename Mesh_domain::Point_3>>
detect_features_on_bbox(const CGAL::Image_3& image, Mesh_domain& domain)
{
  CGAL_IMAGE_IO_CASE(image.image(),
    return detect_features_on_bbox_with_know_word_type<Word>(image, domain)
  );
  CGAL_error_msg("This place should never be reached, because it would mean "
                 "the image word type is a type that is not handled by "
                 "CGAL_ImageIO.");
  return std::vector<std::vector<typename Mesh_domain::Point_3>>();
}

}// namespace internal

/*!
* \ingroup PkgMesh3FeatureDetection
*
* Functor for feature detection in labeled images.
*/
struct Detect_features_on_image_bbox
{
public:

  /*!
  * detects and constructs the polylines that lie at the
  * intersection of two subdomains and the bounding box of the input labeled image.
  * The constructed polylines are added to `domain` for further feature protection.
  *
  * \tparam Mesh_domain class model of `MeshDomainWithFeatures_3`
  *
  * \param image the input image
  * \param domain the mesh domain to be enriched with polyline features
  */
  template<typename Mesh_domain>
  std::vector<std::vector<typename Mesh_domain::Point_3>>
  operator()(const CGAL::Image_3& image, Mesh_domain& domain)
  {
    return internal::detect_features_on_bbox(image, domain);
  }
};


}//end namespace Mesh_3
}//end namespace CGAL


#endif //CGAL_MESH_3_DETECT_FEATURES_ON_IMAGE_BBOX_H
