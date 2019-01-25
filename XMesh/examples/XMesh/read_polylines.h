#ifndef READ_POLYLINES_H
#define READ_POLYLINES_H

#include <cstddef>
#include <vector>
#include <fstream>

template <typename Point_3>
bool read_polylines(const char* fname,
                    std::vector<std::vector<Point_3> >& polylines)
{
  std::ifstream ifs(fname);
  if(ifs.bad()) return false;
  std::size_t n;
  while(ifs >> n) {
    polylines.resize(polylines.size()+1);
    std::vector<Point_3>& polyline = polylines.back();
    while(n-- != 0) {
      Point_3 p;
      ifs >> p;
      if(ifs.fail()) return false;
      polyline.push_back(p);
    }
  }
  if(ifs.bad()) return false;
  else return ifs.eof();
}

#endif // READ_POLYLINES_H
