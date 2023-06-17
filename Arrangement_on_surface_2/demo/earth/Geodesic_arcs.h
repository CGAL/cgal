
#ifndef GEODESIC_ARCS_H
#define GEODESIC_ARCS_H

#include <vector>
#include <qvector3d.h>


class Geodesic_arcs
{
public:

  std::vector<std::vector<QVector3D>> get_approximate_arcs(double error);
  
  //Line_strip_approx  get_approximate_arcs(double error);
};


#endif
