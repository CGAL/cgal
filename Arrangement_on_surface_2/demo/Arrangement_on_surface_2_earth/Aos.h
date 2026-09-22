// Copyright (c) 2023, 2024 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Engin Deniz Diktas <denizdiktas@gmail.com>

#ifndef AOS_H
#define AOS_H

#include <map>
#include <memory>
#include <vector>
#include <qvector3d.h>

#include "Kml_reader.h"

class Aos {
public:
  using Approx_arc = std::vector<QVector3D>;
  using Approx_arcs = std::vector<Approx_arc>;
  //using Arr_handle = void*;
  using Arr_handle = std::shared_ptr<void>;

  static Approx_arc get_approx_identification_curve(double error);

  // this constructs some sample arcs manually (used to check the visual output)
  static Approx_arcs get_approx_arcs(double error);

  // this is used to construct approximate arcs from KML-Placemark data
  static Approx_arcs get_approx_arcs(const Kml::Placemark& placemark,
                                     double error);

  static void get_curves();

  static void check(const Kml::Placemarks& placemarks);

  // Extended check: here we use the extended version of the DCEL to understand
  // what is going on with the additional vertices and faces.
  // INPUT: GIS data
  // OUTPUT: vertices created during arrangement construction
  static std::vector<QVector3D> ext_check(const Kml::Placemarks& placemarks);

  // Same as above function but works by identifying the duplicate nodes
  static std::vector<QVector3D> ext_check_id_based(Kml::Placemarks& placemarks);

  // This function is used to find the newly-created faces during arrangement
  // construction. It uses the extended DCEL structure (see: above 2 functions).
  // NOTE: The need for this function arises due to the following observation:
  // when constructing the arrangement from our data-set I observed a newly
  // created face which is not defined explicitly in the data-set. My intiution
  // tells me that this is the Caspian sea, because it is the only land-locked
  // polygon whose boundaries are not defined in the data-set and its boundaries
  // are defined indirecly by its surrounding countries.
  static Approx_arcs find_new_faces(Kml::Placemarks& placemarks);

  // save the arrangement created with EPEC
  static void save_arr(Kml::Placemarks& placemarks,
                       const std::string& file_name);

  // loads the arrangement created by the save_arr function
  // NOTE: This one loads the arrangement as "Country_arr" type
  static Arr_handle load_arr(const std::string& file_name);

  static Arr_handle construct(Kml::Placemarks& placemarks);
  //static std::vector<QVector3D> get_triangles(Arr_handle arrh);

  //using Country_triangles_map = std::map<std::string, std::vector<QVector3D>>;
  //static Country_triangles_map get_triangles_by_country(Arr_handle arrh);

  using Country_color_map = std::map<std::string, std::size_t>;
  static Country_color_map get_color_mapping(Arr_handle arrh);

  static std::string locate_country(Arr_handle arrh, const QVector3D& point);

  // this will get the approximate arcs of face-edges from the arrangement
  // NOTE: this is similar to "get_approx_arcs(KML::Placemarks&, double)" above!
  static Approx_arcs get_approx_arcs_from_faces_edges(Arr_handle arrh,
                                                      double error);
};

#endif
