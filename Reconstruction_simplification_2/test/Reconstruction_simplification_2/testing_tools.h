#include <fstream>
#include <vector>
#include <string>

template <class PointMassList, class Point>
void load_xy_file(const std::string& filename, PointMassList& points)
{
  std::ifstream ifs(filename.c_str());
  Point point;
  while (ifs >> point)
    points.push_back(std::make_pair(point, 1));

  ifs.close();
}

template <class Point>
void load_xy_file_points(const std::string& fileName, std::vector<Point>& points)
{
  std::ifstream ifs(fileName.c_str());
  Point point;
  while (ifs >> point)
  {
    points.push_back(point);
  }
  ifs.close();
}
