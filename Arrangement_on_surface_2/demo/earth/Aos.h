
#ifndef AOS_H
#define AOS_H

#include <vector>
#include <qvector3d.h>

#include "Kml_reader.h"


class Aos
{
public:

  static void check(const Kml::Placemarks& placemarks);
};


#endif
