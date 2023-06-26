
#ifndef KML_READER_H
#define KML_READER_H

#include <string>
#include <vector>

#include <qvector3d.h>

class Kml
{
public:
  struct Node
  {
    double lon, lat;

    bool operator == (const Node& r) const  
    { 
      return  (lon == r.lon) && (lat == r.lat);
    }

    QVector3D get_coords_3d(const double r=1.0) const;
  };

  struct LinearRing
  {
    std::vector<Node> nodes;
  };

 
  struct Placemark
  {
    std::vector<LinearRing> polygons;
    std::string name;
  };

  using Nodes = std::vector<Node>;
  using Placemarks = std::vector<Placemark>;


  static Placemarks read(const std::string& file_name);

  static Nodes get_duplicates(const Placemarks& placemarks);
};


#endif
