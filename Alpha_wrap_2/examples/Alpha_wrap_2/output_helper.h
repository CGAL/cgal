#ifndef CGAL_ALPHA_WRAP_2_EXAMPLES_OUTPUT_HELPER_H
#define CGAL_ALPHA_WRAP_2_EXAMPLES_OUTPUT_HELPER_H

#include <CGAL/IO/WKT.h>

#include <fstream>
#include <string>

// @todo wkt cannot be opened with the demo.....
// @todo multipolygon would make most sense? see also "purge inner CCs" + arrange_offset_polygons from SLS2

std::string generate_output_name(std::string input_name,
                                 const double alpha,
                                 const double offset)
{
  input_name = input_name.substr(input_name.find_last_of("/") + 1, input_name.length() - 1);
  input_name = input_name.substr(0, input_name.find_last_of("."));

#ifdef CGAL_AW2_OUTPUT_TO_WKT
  std::string output_name = input_name
                            + "_" + std::to_string(static_cast<int>(alpha))
                            + "_" + std::to_string(static_cast<int>(offset)) + "-wrap.wkt";
#else
  std::string output_name = input_name
                            + "_" + std::to_string(static_cast<int>(alpha))
                            + "_" + std::to_string(static_cast<int>(offset)) + "-wrap.txt";
#endif

  return output_name;
}

// write a range of Polygon_2 items as polylines
template <typename Polygons>
void write_polygons(const std::string& filename,
                    const Polygons& polygons)
{
  std::ofstream out(filename.c_str());
  out.precision(17);

#ifdef CGAL_AW2_OUTPUT_TO_WKT
  for(const auto& p : polygons)
    write_polygon_WKT(out, p);
#else
  for(const auto& p : polygons)
  {
    out << p.size();
    for(const auto& pt : p)
      out << " " << pt << " 0";
    out << std::endl;
  }
#endif
}

#endif // CGAL_ALPHA_WRAP_2_EXAMPLES_OUTPUT_HELPER_H
