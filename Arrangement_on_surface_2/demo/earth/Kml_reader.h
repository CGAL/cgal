
#ifndef KML_READER_H
#define KML_READER_H

#include <string>
#include <vector>


class Kml
{
public:
  struct Node
  {
    double lon, lat;
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

  using Placemarks = std::vector<Placemark>;


  static Placemarks read(const std::string& file_name);
};


#endif
