#ifndef POLYGON_IO_H
#define POLYGON_IO_H

#include <typedefs.h>
#include <CGAL/IO/Window_stream.h>
#include <LEDA/string.h>
#include <LEDA/polygon.h>
#include <LEDA/point.h>

void convert_leda_polygon_to_CGAL(const leda_polygon &ledaP,
                                  Polygon_2 &p);

void convert_point_list_to_leda_poly(const Polygon_2& p,leda_polygon &ledaP);

void input_polygon(Polygon_2& p, CGAL::Window_stream& W);

template <class OutputStream>
void label_vertices(OutputStream& OS, const Polygon_2& poly);

template <class OutputStream>
void draw_a_polygon(OutputStream& OS, const Polygon_2& poly, CGAL::Color color);

template <class OutputStream>
void draw_poly_list(OutputStream& OS, const std::list<Polygon_2>& poly_list, 
                    int line_width, CGAL::Color color,  int& count,
                    bool show_coords);

bool read_poly_from_file(Polygon_2& polygon, leda_string file_name);

bool save_partition_p_list_to_file(const std::list<Polygon_2>& partition_p_list,
                                   leda_string& current_file_name, 
                                   leda_string file_suffix);

bool save_poly_to_file(const Polygon_2& polygon, leda_string file_name);

#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#include <polygon_io.C>
#endif // CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION

#endif // POLYGON_IO_H
