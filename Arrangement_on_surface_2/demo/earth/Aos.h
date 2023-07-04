
#ifndef AOS_H
#define AOS_H

#include <vector>
#include <qvector3d.h>

#include "Kml_reader.h"


class Aos
{
public:
  using Approx_arcs = std::vector<std::vector<QVector3D>>;

  // this constructs some sample arcs manually (used to check the visual output)
  static Approx_arcs get_approx_arcs(double error);

  // this is used to construct approximate arcs from KML-Placemark data
  static Approx_arcs get_approx_arcs(const Kml::Placemark& placemark, double error);


  static void get_curves();

  static void check(const Kml::Placemarks& placemarks);

  // Extended check: here we use the extended version of the DCEL to understand
  // what is going on with the additional vertices and faces.
  // INPUT: GIS data
  // OUTPUT: vertices created during arrangement construction
  static std::vector<QVector3D> ext_check(const Kml::Placemarks& placemarks);
};


#endif
