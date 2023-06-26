
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

    bool operator == (const Node& r) 
    { 
      return  (lon == r.lon) && (lat == r.lat);
    }
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

  static void check_duplicates(const Placemarks& placemarks);
};


#endif
