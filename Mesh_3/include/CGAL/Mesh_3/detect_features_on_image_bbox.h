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
void detect_features_on_bbox_with_know_word_type(const CGAL::Image_3& image,
                                                 Mesh_domain& domain)
{
  using Gt = typename Mesh_domain::R;
  using Point_3 = typename Gt::Point_3;
  using Polyline_type = std::vector<Point_3>;
  using Polylines = std::vector<Polyline_type>;

  Polylines polylines_on_bbox;
  CGAL::polylines_to_protect<Point_3, Word_type>(image, polylines_on_bbox);

  domain.add_features(polylines_on_bbox.begin(), polylines_on_bbox.end());
}

template<typename Mesh_domain>
void detect_features_on_bbox(const CGAL::Image_3& image, Mesh_domain& domain)
{
  CGAL_IMAGE_IO_CASE(image.image(),
    return detect_features_on_bbox_with_know_word_type<Word>(image, domain)
  );
  CGAL_error_msg("This place should never be reached, because it would mean "
                 "the image word type is a type that is not handled by "
                 "CGAL_ImageIO.");
}

}// namespace internal

struct Detect_features_on_image_bbox
{
  template<typename Mesh_domain>
  void operator()(const CGAL::Image_3& image, Mesh_domain& domain)
  {
    internal::detect_features_on_bbox(image, domain);
  }
};


}//end namespace Mesh_3
}//end namespace CGAL


#endif //CGAL_MESH_3_DETECT_FEATURES_ON_IMAGE_BBOX_H
