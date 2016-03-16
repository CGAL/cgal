// Copyright (c) 2016  INRIA Sophia-Antipolis (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: $
// $Id: $
//
//
// Author(s)     : Clement Jamin

#ifndef CGAL_TC_WRITE_RIB_FILE_H
#define CGAL_TC_WRITE_RIB_FILE_H

#include "utilities.h"

#include <CGAL/basic.h>

#include <tuple>
#include <string>

namespace CGAL {
namespace Tangential_complex_ {
  
template <typename PointRandomAccessRange, typename SimplexRange>
class RIB_exporter
{
  typedef typename PointRandomAccessRange::value_type  Point;
public:

  typedef std::tuple<double, double, double, double>  Color; // RGBA

  // Constructor
  RIB_exporter(
    PointRandomAccessRange const& points,
    SimplexRange const& simplices,
    std::ofstream &out,
    std::string const& rendered_image_filename = "export.tif",
    bool is_preview = false, // low-quality
    int image_width = 1920,
    int image_height = 1080,
    Color const& triangle_color = std::make_tuple(0., 1., 1., 1.),
    bool ambient_light = true,
    double ambient_intensity = 0.20,
    bool shadow = true,
    double shadow_intensity = 0.85)
  : m_points(points),
    m_simplices(simplices),
    m_out(out),
    m_rendered_image_filename(rendered_image_filename),
    m_is_preview(is_preview),
    m_image_width(image_width),
    m_image_height(image_height),
    m_current_color(0., 0., 0., 0.),
    m_current_alpha(1),
    m_triangle_color(triangle_color),
    m_ambient_light(ambient_light),
    m_ambient_intensity(ambient_intensity),
    m_shadow(shadow),
    m_shadow_intensity(shadow_intensity)
  {
    m_out.precision(8);
  }

  void write_file()
  {
    write_header();
    write_lights();
    write_triangles();

    m_out << "WorldEnd\n";
  }

private:

  void write_header()
  {
    m_out << "Option \"searchpath\" \"shader\" "
      "\".:./shaders:%PIXIE_SHADERS%:%PIXIEHOME%/shaders\"\n";

    if (m_is_preview)
    {
      m_out << "Attribute \"visibility\" \"specular\" 1\n"
        << "Attribute \"visibility\" \"transmission\" 1\n\n";
    }

    m_out << "Display \"" << m_rendered_image_filename << "\" \"file\" \"rgb\"\n";

    if (!m_is_preview)
    {
      m_out << "Format " << m_image_width << " " << m_image_height << " 1\n";
    }
    else
    {
      double ratio = double(m_image_height) / double(m_image_width);

      int width = (ratio < 1.) ? 300 : int(300. / ratio);
      int height = (ratio < 1.) ? int(ratio * 300.) : 300;

      m_out << "Format " << width << " " << height << " 1\n";
    }


    if (m_image_width > m_image_height)
    {
      double ratio = double(m_image_height) / double(m_image_width);
      m_out << "ScreenWindow -1 1 " << -ratio << " " << ratio << "\n";
    }
    else if (m_image_height > m_image_width)
    {
      double ratio = double(m_image_width) / double(m_image_height);
      m_out << "ScreenWindow " << -ratio << " " << ratio << " -1 1\n";
    }

    m_out << "Projection \"perspective\" \"fov\" 45\n"
      << "Translate 0 0 15\n"
      << "PixelSamples 4 4\n"
      << "PixelFilter \"catmull-rom\" 3 3\n"
      << "ShadingInterpolation \"smooth\"\n"
      << "Rotate 180 0 0 1\n"
      << "WorldBegin\n";
  }


  void write_lights()
  {
    if (!m_is_preview)
    {
      // ShadowLight
      m_out << "LightSource \"shadowdistant\" 1 \"from\" [0 0 0] \"to\" [0 0 1]"
        << " \"shadowname\" \"raytrace\" \"intensity\" "
        << m_shadow_intensity << "\n";

      // Ambient light
      m_out << "LightSource \"ambientlight\" 2 \"intensity\" "
        << m_ambient_intensity << "\n";
    }
    else
    {
      m_out << "LightSource \"distantLight\" 1 \"from\" [0 0 0] \"to\" [0 0 1]"
        << " \"intensity\" 0.85\n";
    }

    // Background light
    m_out << "LightSource \"ambientlight\" 99 \"intensity\" 1\n";

    // Turn background light OFF
    turn_background_light(false);
  }

  void turn_background_light(bool turn_on)
  {
    if (!turn_on)
    {
      m_out << "Illuminate 1 1" << std::endl;
      if (!m_is_preview) 
        m_out << "Illuminate 2 1" << std::endl;
      m_out << "Illuminate 99 0" << std::endl;
    }
    else
    {
      m_out << "Illuminate 1 0" << std::endl;
      if (!m_is_preview) 
        m_out << "Illuminate 2 0" << std::endl;
      m_out << "Illuminate 99 1" << std::endl;
    }
  }

  // CJTODO
  void write_background(const Color& color)
  {
    write_turn_background_light(false);

    /*m_out << "Surface \"constant\"" << std::endl;
    write_color(color, false);

    double corner = zmax_ * 2.;
    double depth_pos = zmax_ * 1.6;

    m_out << "Polygon \"P\" [";
    m_out << " " << -corner << " " << -corner << " " << depth_pos << " ";
    m_out << " " << corner << " " << -corner << " " << depth_pos << " ";
    m_out << " " << corner << " " << corner << " " << depth_pos << " ";
    m_out << " " << -corner << " " << corner << " " << depth_pos << " ";
    m_out << "]" << std::endl;*/
  }


  void write_color(Color const& color, bool use_transparency)
  {
    if (m_current_color == color)
      return;

    m_current_color = color;

    // Write opacity data
    if (use_transparency)
      write_opacity(std::get<3>(color));

    // Write color data
    m_out << "Color [ " << std::get<0>(color) << " " << std::get<1>(color) 
      << " " << std::get<2>(color) << " ]\n";
  }

  void write_opacity(const double alpha)
  {
    if (m_current_alpha == alpha)
      return;

    m_current_alpha = alpha;

    // Write opacity data
    m_out << "Opacity " << alpha << " " << alpha << " " << alpha << std::endl;
  }

  void write_point(Point const& p)
  {
    m_out << " " << p[0] << " " << p[1] << " " << p[2] << " ";
  }

  void write_triangles()
  {
    m_out << "Surface \"plastic\" \"Ka\" 0.65 \"Kd\" 0.85 \"Ks\" 0.25 \"roughness\" 0.1" << std::endl;

    for (auto simplex : m_simplices)
    {
      std::vector<std::set<std::size_t> > triangles;
      // Get the triangles composing the simplex
      combinations(simplex, 3, std::back_inserter(triangles));
      for (auto const& t : triangles)
        write_triangle(t);
    }
  }

  template <typename PointIndexRange>
  void write_triangle(PointIndexRange const& t)
  {
    // Color
    write_color(m_triangle_color, true);

    // Triangle
    m_out << "Polygon \"P\" [";
    for (auto idx : t)
      write_point(m_points[idx]);
    m_out << "]" << std::endl;

    // Edges (will be drawn later on)
    /*add_edge(p, q, edge_color);
    add_edge(p, r, edge_color);
    add_edge(q, r, edge_color);

    // Vertices (will be drawn later on)
    add_vertex(p, edge_color);
    add_vertex(q, edge_color);
    add_vertex(r, edge_color);*/
  }

  //===========================================================================
  
  PointRandomAccessRange const& m_points;
  SimplexRange const& m_simplices;
  std::ofstream &m_out;
  std::string m_rendered_image_filename;
  bool m_is_preview;
  int m_image_width;
  int m_image_height;
  Color m_current_color;
  Color m_triangle_color;
  double m_current_alpha;
  bool m_ambient_light;
  double m_ambient_intensity;
  bool m_shadow;
  double m_shadow_intensity;
};

} // namespace Tangential_complex_
} //namespace CGAL

#endif // CGAL_TC_WRITE_RIB_FILE_H
