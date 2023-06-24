
#ifndef GEODESIC_ARCS_H
#define GEODESIC_ARCS_H

#include <vector>
#include <qvector3d.h>

#include "Kml_reader.h"


class Geodesic_arcs
{
public:
  using Approx_arcs = std::vector<std::vector<QVector3D>>;

  
  Approx_arcs get_approx_arcs(double error);
  
  // generate approximate arcs from KML data
  Approx_arcs get_approx_arcs(const Kml::Placemarks& placemarks, double error);
};


#endif
