
#ifndef AOS_H
#define AOS_H

#include <vector>
#include <qvector3d.h>

#include "Kml_reader.h"


class Aos
{
public:

  static void check(const Kml::Placemarks& placemarks);

  // Extended check: here we use the extended version of the DCEL to understand
  // what is going on with the additional vertices and faces.
  // INPUT: GIS data
  // OUTPUT: vertices created during arrangement construction
  static void ext_check(const Kml::Placemarks& placemarks);
};


#endif
