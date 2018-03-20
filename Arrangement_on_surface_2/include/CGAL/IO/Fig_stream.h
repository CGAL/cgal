// Copyright (c) 2006,2007,2008,2009,2010,2011 Tel-Aviv University (Israel).
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
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
// 
// Author(s)     : Ron Wein           <wein@post.tau.ac.il>

#ifndef CGAL_FIG_STREAM_H
#define CGAL_FIG_STREAM_H

#include <CGAL/license/Arrangement_on_surface_2.h>


#include <CGAL/basic.h>
#include <CGAL/Polygon_2.h>

#include <vector>
#include <fstream>
#include <cstdio>

namespace CGAL {

/*!
 * FIG colors.
 */
enum Fig_color
{
  // Predefined colors:
  FIG_BLACK = 0,
  FIG_BLUE = 1,
  FIG_GREEN = 2,
  FIG_CYAN = 3,
  FIG_RED = 4,
  FIG_MAGENTA = 5,
  FIG_YELLOW = 6,
  FIG_WHITE = 7,
  FIG_BLUE_1 = 8, FIG_BLUE_2 = 9, FIG_BLUE_3 = 10, FIG_BLUE_4 = 11,
  FIG_GREEN_1 = 12, FIG_GREEN_2 = 13, FIG_GREEN_3 = 14,
  FIG_CYAN_1 = 15, FIG_CYAN_2 = 16, FIG_CYAN_3 = 17,
  FIG_RED_1 = 18, FIG_RED_2 = 19, FIG_RED_3 = 20,
  FIG_MAGENTA_1 = 21, FIG_MAGENTA_2 = 22, FIG_MAGENTA_3 = 23,
  FIG_BROWN_1 = 24, FIG_BROWN_2 = 25, FIG_BROWN_3 = 26,
  FIG_PINK_1 = 27, FIG_PINK_2 = 28, FIG_PINK_3 = 29, FIG_PINK_4 = 30,
  FIG_GOLD = 31,

  // User-defined colors:
  FIG_FIRST_USER_DEFINED_COLOR = 32,
  FIG_LAST_USER_DEFINED_COLOR = 543
};

/*!
 * FIG line styles.
 */
enum Fig_line_style
{
  FIG_SOLID = 0,
  FIG_DASHED = 1,
  FIG_DOTTED = 2,
  FIG_DASH_DOTTED = 3,
  FIG_DASH_DOUBLE_DOTTED = 4,
  FIG_DASH_TRIPLE_DOTTED = 5
};

#define FIG_DEFAULT_STYLE_VALUE   4.0
    
/*!
 * FIG fill styles.
 */
enum Fig_fill_style
{
  FIG_NOT_FILLED = -1,
  FIG_FILL_BLACK = 0,
  /// Values from 1 to 19 are shades of the color from darker to lighter.
  FIG_FILLED = 20,
  /// Values from 21 to 39 are tints of the color from the color to white.
  FIG_FILL_WHITE = 40,
  FIG_LEFT_DIAG_30DEG = 41,
  FIG_RIGHT_DIAG_30DEG = 42,
  FIG_CROSS_DIAG_30DEG = 43,
  FIG_LEFT_DIAG_45DEG = 44,
  FIG_RIGHT_DIAG_45DEG = 45,
  FIG_CROSS_DIAG_45DEG = 46,
  FIG_HORIZONTAL_BRICKS = 47,
  FIG_VERTICAL_BRICKS = 48,
  FIG_HORIZONTAL_LINES = 49,
  FIG_VERTICAL_LINES = 50,
  FIG_CROSS_LINES = 51,
  FIG_HORIZONTAL_RIGHT_SHINGLES = 52,
  FIG_HORIZONTAL_LEFT_SHINGLES = 53,
  FIG_VERTICAL_RIGHT_SHINGLES = 54,
  FIG_VERTICAL_LEFT_SHINGLES = 55,
  FIG_FISH_SCALES = 56,
  FIG_SMALL_FISH_SCALES = 57,
  FIG_CIRCLES = 58,
  FIG_HEXAGONS = 59,
  FIG_OCTAGONS = 60,
  FIG_HORIZONTAL_TIRE_TREADS = 61,
  FIG_VERTICAL_TIRE_TREADS = 62
};

/*!
 * FIG arrow types.
 */
enum Fig_arrow_type
{
  FIG_STICK = 0,
  FIG_TRIANGLE = 1,
  FIG_INDENTED_BUTT = 2,
  FIG_POINTED_BUTT = 3
};

/*!
 * Arrow modes (not based on the FIG format).
 */
enum Fig_arrow_mode
{
  FIG_NO_ARROW,
  FIG_FORWARD_ARROW,
  FIG_BACKWARD_ARROW,
  FIG_BOTH_ARROWS
};

/*!
 * Point styles (not based on the FIG format).
 */
enum Fig_point_style
{
  FIG_CROSS,
  FIG_PLUS,
  FIG_CIRCLE,
  FIG_DISC,
  FIG_SQUARE,
  FIG_BOX,
  FIG_RHOMBUS,
  FIG_DIAMOND
};

/*!
 * FIG fonts.
 */
enum Fig_font
{
  FIG_ROMAN = 1,
  FIG_BOLD = 2,
  FIG_ITALIC = 3,
  FIG_SANS_SERIF = 4,
  FIG_TYPEWRITER = 5
};

/*!
 * Depth constants.
 */
enum Fig_depth
{
  FIG_MIN_DEPTH = 0,
  FIG_DEFAULT_DEPTH = 50,
  FIG_MAX_DEPTH = 99
};

/*!
 * \class A class for writing geometric objects in a FIG format (version 3.2).
 * For more details, see: http://www.xfig.org/userman/fig-format.html
 */
template <class Kernel_>
class Fig_stream
{
public:

  typedef Kernel_                                Kernel;
  
  // Define the kernel objects.
  typedef typename Kernel::FT                    NT;
  typedef typename Kernel::Point_2               Point_2;
  typedef typename Kernel::Segment_2             Segment_2;
  typedef typename Kernel::Ray_2                 Ray_2;
  typedef typename Kernel::Line_2                Line_2;
  typedef typename Kernel::Triangle_2            Triangle_2;
  typedef typename Kernel::Iso_rectangle_2       Iso_rectangle_2;
  typedef Polygon_2<Kernel>                      Polygon_2;
  typedef typename Kernel::Circle_2              Circle_2;

protected:

  // Data members:
  std::ofstream           _ofile;       // The output file.

  Iso_rectangle_2         _bound_rect;  // A rectangle bounding the workspace.
  int                     _width;       // Figure width (in pixels).
  int                     _height;      // Figure height (in pixels).
  double                  _scale;       // Scaling factor.

  int                     _depth;

  Fig_color               _color;
  int                     _line_width;
  Fig_line_style          _line_style;
  double                  _style_value;

  Fig_color               _fill_color;
  Fig_fill_style          _fill_style;

  Fig_point_style         _point_style;
  NT                      _point_size;

  Fig_arrow_mode          _arrow_mode;
  Fig_arrow_type          _arrow_type;
  NT                      _arrow_width;
  NT                      _arrow_height;

  Fig_font                _font;
  int                     _font_size;

  bool                    colors[FIG_LAST_USER_DEFINED_COLOR + 1];

  // Kernel functors.
  typename Kernel::Intersect_2       intersect_func;

private:

  // Copy constructor and assignment operator - not supported.
  Fig_stream (const Fig_stream<Kernel>& );
  const Fig_stream<Kernel>& operator= (const Fig_stream<Kernel>& );

public:

  /// \name Constructors and destructor.
  //@{

  /*!
   * Default constructor.
   */
  Fig_stream () :
    _width (0),
    _height (0),
    _scale (0),
    _depth (FIG_DEFAULT_DEPTH),
    _color (FIG_BLACK),
    _line_width (1),
    _line_style (FIG_SOLID),
    _style_value (FIG_DEFAULT_STYLE_VALUE),
    _fill_color (FIG_WHITE),
    _fill_style (FIG_NOT_FILLED),
    _point_style (FIG_DISC),
    _arrow_mode (FIG_NO_ARROW),
    _arrow_type (FIG_STICK),
    _font (FIG_ROMAN),
    _font_size (12)
  {
    // Reset all colors.
    _reset_colors ();

    // Construct the necessary kernel functors.
    Kernel        ker;

    intersect_func = ker.intersect_2_object();
  }

  /*!
   * Constructor.
   * \param filename The name of the output FIG file.
   * \param rect A rectangle bounding the logical drawing area.
   * \param width The physical width of the figure (in FIG units).
   * \param height The physical height of the figure (in FIG units).
   * \pre The bounding rectangle is valid and the physical dimensions are 
   *      both positive.
   */
  Fig_stream (const char *filename,
              const Iso_rectangle_2& rect,
              const int& width = 12000,
              const int& height = 12000) :
    _width (0),
    _height (0),
    _scale (0),
    _depth (FIG_DEFAULT_DEPTH),
    _color (FIG_BLACK),
    _line_width (1),
    _line_style (FIG_SOLID),
    _style_value (FIG_DEFAULT_STYLE_VALUE),
    _fill_color (FIG_WHITE),
    _fill_style (FIG_NOT_FILLED),
    _point_style (FIG_DISC),
    _arrow_mode (FIG_NO_ARROW),
    _arrow_type (FIG_STICK),
    _font (FIG_ROMAN),
    _font_size (12)
  {
    // Reset all colors.
    _reset_colors ();

    // Construct the necessary kernel functors.
    Kernel        ker;

    intersect_func = ker.intersect_2_object();

    // Open the file.
    open (filename, rect, width, height);
  }

  /*!
   * Destructor.
   */
  virtual ~Fig_stream ()
  {
    _ofile.close();
  }
  //@}

  /// \name Openning and closing the file.
  //@{

  /*!
   * Check whether the file is open.
   */
  bool is_open ()
  {
    return (_ofile.is_open());
  }

  /*!
   * Open a FIG file.
   * \param filename The name of the output FIG file.
   * \param rect A rectangle bounding the logical drawing area.
   * \param width The physical width of the figure (in FIG units).
   * \param height The physical height of the figure (in FIG units).
   * \pre The bounding rectangle is valid and the physical dimensions are 
   *      both positive.
   * \return Whether the file was successfully opened.
   */
  bool open (const char *filename,
             const Iso_rectangle_2& rect,
             const int& width = 12000,
             const int& height = 12000)
  {
    CGAL_precondition (width > 0);
    CGAL_precondition (height > 0);
    CGAL_precondition (rect.xmax() > rect.xmin());
    CGAL_precondition (rect.ymax() > rect.ymin());

    // Reset all colors.
    _reset_colors ();

    // Close the current output file, if necessary.
    if (_ofile.is_open())
      _ofile.close();

    // Open the output file.
    _ofile.open (filename);

    // Set the logical and physical dimensions.
    _bound_rect = rect;
    _width = width;
    _height = height;

    // Compute the scale.
    const double x_scale = width / CGAL::to_double(rect.xmax() - rect.xmin());
    const double y_scale = height / CGAL::to_double(rect.ymax() - rect.ymin());
    _scale = (x_scale < y_scale) ? x_scale : y_scale;

    // Set the default point size and arrow dimensions.
    _point_size = (rect.xmax() - rect.xmin()) / NT(500);
    _arrow_width = _point_size;
    _arrow_height = 2*_point_size;

    // End here if the file is not opened.
    if (! _ofile.is_open())
      return (false);

    // Write the FIG header.
    _ofile << "#FIG 3.2" << std::endl;
    _ofile << "Landscape" << std::endl;
    _ofile << "Center" << std::endl;
    _ofile << "Inches" << std::endl;
    _ofile << "Letter" << std::endl;  
    _ofile << "100.00" << std::endl;
    _ofile << "Single" << std::endl;
    _ofile << "-2" << std::endl;
    _ofile << "1200 2" << std::endl;

    return (true);
  }

  /*!
   * Close the FIG file.
   */
  void close ()
  {
    if (_ofile.is_open())
      _ofile.close();

    // Reset all colors.
    _reset_colors ();
  }
  //@}

  /// \name Accessing drawing properties.
  //@{
  
  /*!
   * Get the workspace bounding rectangle.
   */
  Iso_rectangle_2 bounding_rect() const
  {
      return (_bound_rect);
  }

/*   /\*! */
/*    * Set the workspace bounding rectangle. */
/*    *\/ */
/*   void set_bounding_rect(const Iso_rectangle_2& rect) */
/*   { */
/*       _bound_rect = rect; */
/*   } */


  /*!
   * Get the physical width of the fig
   */  
  int width () const
  {
    return _width;
  }

  /*!
   * Get the physical height of the fig
   */  
  int height () const
  {
    return _height;
  }

  /*!
   * Get the depth.
   */
  int depth () const
  {
    return (depth);
  }

  /*!
   * Get the color.
   */
  Fig_color color () const
  {
    return (_color);
  }

  /*!
   * Get the line width.
   */
  int line_width () const
  {
    return (line_width);
  }

  /*!
   * Get the line style.
   */
  Fig_line_style line_style () const
  {
    return (_line_style);
  }

  /*!
   * Get the style value.
   */
  double style_value () const
  {
    return (_style_value);
  }

  /*!
   * Get the fill color.
   */
  Fig_color fill_color () const
  {
    return (_fill_color);
  }

  /*!
   * Get the fill style.
   */
  Fig_fill_style fill_style () const
  {
    return (_fill_style);
  }

  /*!
   * Get the point style.
   */
  Fig_point_style point_style () const
  {
    return (_point_style);
  }

  /*!
   * Get the point size.
   */
  const NT& point_size () const
  {
    return (_point_size);
  }

  /*!
   * Get the arrow drawing mode (this mode is relevent when drawing segments,
   * polylines, circular arcs or splines).
   */
  Fig_arrow_mode arrow_mode () const
  {
    return (_arrow_mode);
  }

  /*!
   * Get the arrow type.
   */
  Fig_arrow_type arrow_type () const
  {         
    return (_arrow_type);
  }

  /*!
   * Get the arrow width.
   */
  const NT& arrow_width () const
  {
    return (_arrow_width);
  }

  /*!
   * Get the arrow height.
   */
  const NT& arrow_height () const
  {
    return (_arrow_height);
  }

  /*!
   * Get the font.
   */
  Fig_font font () const
  {
    return (_font);
  }

  /*!
   * Get the font size.
   */
  int font_size () const
  {
    return (_font_size);
  }
  //@}

  /// \name Set the drawing properties.
  //@{

  /*!
   * Set the depth.
   */
  void set_depth (const int& depth)
  {
    if (depth < static_cast<int>(FIG_MIN_DEPTH))
      _depth = static_cast<int>(FIG_MIN_DEPTH);
    else if (depth > static_cast<int>(FIG_MAX_DEPTH))
      _depth = static_cast<int>(FIG_MAX_DEPTH);
    else
      _depth = depth;

    return;
  }
  
  /*!
   * Set the color.
   * \pre The color must be defined.
   */
  void set_color (const Fig_color& color)
  {
    CGAL_precondition (color_defined (color));

    if (color_defined (color))
        _color = color;
    
    return;
  }

  /*!
   * Set the line width.
   */
  void set_line_width (const unsigned int& width)
  {
    _line_width = static_cast<int>(width);
    return;
  }

  /*!
   * Set the line style.
   */
  void set_line_style (const Fig_line_style& style)
  {
    _line_style = style;
    return;
  }

  /*!
   * Set the style value.
   */
  void set_style_value (const double& val)
  {
    CGAL_precondition (val > 0);

    _style_value = val;
    return;
  }

  /*!
   * Set the fill color.
   * \pre The color must be defined.
   */
  void set_fill_color (const Fig_color& color)
  {
    CGAL_precondition (color_defined (color));

    if (color_defined (color))
        _fill_color = color;

    return;
  }

  /*!
   * Set the fill style.
   */
  void set_fill_style (const Fig_fill_style& style)
  {
    _fill_style = style;
    return;
  }

  /*!
   * Set the point style.
   */
  void set_point_style (const Fig_point_style& style)
  {
    _point_style = style;
    return;
  }

  /*!
   * Set the point size.
   */
  void set_point_size (const NT& size)
  {
    _point_size = CGAL::abs(size);
    return;
  }

  /*!
   * Set the arrow drawing mode. This mode will be applied when drawing
   * segments, polylines, circular arcs or splines.
   */
  void set_arrow_mode (const Fig_arrow_mode& mode)
  {
    _arrow_mode = mode;
    return;
  }

  /*!
   * Set the arrow type.
   */
  void set_arrow_type (const Fig_arrow_type& type)
  {         
    _arrow_type = type;
    return;
  }

  /*!
   * Set the arrow width.
   */
  void set_arrow_width (const NT& width)
  {
    _arrow_width = CGAL::abs(width);
    return;
  }

  /*!
   * Get the arrow height.
   */
  void set_arrow_height (const NT& height)
  {
    _arrow_height = CGAL::abs(height);
    return;
  }

  /*!
   * Set the font.
   */
  void set_font (const Fig_font& font)
  {
    _font = font;
    return;
  }

  /*!
   * Set the font size.
   */
  void set_font_size (const unsigned int& size)
  {
    _font_size = static_cast<int>(size);
    return;
  }
  //@}

  /// \name Defining colors.
  //@{

  /*!
   * Check if a color is defined.
   */
  bool color_defined (const Fig_color& color) const
  {
    int     col = static_cast<int>(color);

    if (col < 0 || col > FIG_LAST_USER_DEFINED_COLOR)
        return (false);

    return (colors[col]);
  }

  /*!
   * Add a user-defined color.
   * Use this function after openning the FIG stream and before writing any
   * other object (i.e. before calling the write_<object> () functions).
   * \param color The color.
   * \param r The red component (0 - 255).
   * \param g The green component (0 - 255).
   * \param b The blue component (0 - 255).
   * \pre The color must be undefined.
   */
  void define_color (const Fig_color& color,
                     const unsigned char& r, 
                     const unsigned char& g, 
                     const unsigned char& b)
  {
    CGAL_precondition (color_defined (color));
    CGAL_precondition (_ofile.is_open());

    if (color_defined (color))
      return;

    // Prepare a string desribing the color.
    char    color_desc [10];

    sprintf ("#%02x%02x%02x", r, g, b);

    // Write the color to the FIG file.
    _ofile << "0 "                        // Desginates a color pseudo-object.
           << static_cast<int>(color) << ' '
           << color_desc << std::endl;

    // Mark that the color is now defined.
    colors[static_cast<int>(color)] = true;

    return;
  }

  //@}

  /// \name Writing objects.
  //@{

  /*!
   * Write a point.
   */
  void write_point (const Point_2& p)
  {
    CGAL_precondition (_ofile.is_open());

    //is the point outside the iso-rectangle?
    if(_bound_rect.has_on_unbounded_side(p))
      return;

    switch (_point_style)
    {
    case (FIG_CROSS):
    case (FIG_PLUS):
    {
      // Draw two segments intersecting at p.
      Point_2  s1, t1;
      Point_2  s2, t2;

      if (_point_style == FIG_PLUS)
      {
        // Draw a '+'.
        s1 = Point_2 (p.x() - _point_size, p.y());
        t1 = Point_2 (p.x() + _point_size, p.y());
        s2 = Point_2 (p.x(), p.y() - _point_size);
        t2 = Point_2 (p.x(), p.y() + _point_size);
      }
      else
      {
        // Draw an 'x'.
        s1 = Point_2 (p.x() - _point_size, p.y() - _point_size);
        t1 = Point_2 (p.x() + _point_size, p.y() + _point_size);
        s2 = Point_2 (p.x() - _point_size, p.y() + _point_size);
        t2 = Point_2 (p.x() + _point_size, p.y() - _point_size);
      }

      // Draw solid lines with width 1.
      _write_segment (Segment_2(s1, t1),
                      _color, 1, FIG_SOLID, _style_value,
                      false);
      _write_segment (Segment_2(s2, t2),
                      _color, 1, FIG_SOLID, _style_value,
                      false);

      break;
    }

    case (FIG_CIRCLE):
    {
      // Draw an empty circle (use a solid line with width 1).
      _write_ellipse (p,
                      CGAL::square(_point_size), 
                      CGAL::square(_point_size),
                      _color, 1, FIG_SOLID, _style_value,
                      FIG_WHITE, FIG_NOT_FILLED);

      break;
    }

    case (FIG_DISC):
    {
      // Draw a filled disc.
      _write_ellipse (p,
                      CGAL::square(_point_size), 
                      CGAL::square(_point_size),
                      _color, 1, FIG_SOLID, _style_value,
                      _color, FIG_FILLED);

      break;
    }

    case (FIG_SQUARE):
    case (FIG_BOX):
    case (FIG_RHOMBUS):
    case (FIG_DIAMOND):
    {
      // Prepare the rectangle vertices.
      std::vector<Point_2>   vertices (4);

      if (_point_style == FIG_SQUARE || _point_style == FIG_BOX)
      {
        vertices[0] = Point_2 (p.x() - _point_size, p.y() - _point_size);
        vertices[1] = Point_2 (p.x() - _point_size, p.y() + _point_size);
        vertices[2] = Point_2 (p.x() + _point_size, p.y() + _point_size);
        vertices[3] = Point_2 (p.x() + _point_size, p.y() - _point_size);
      }
      else
      {
        vertices[0] = Point_2 (p.x(), p.y() - _point_size);
        vertices[1] = Point_2 (p.x() + _point_size, p.y());
        vertices[2] = Point_2 (p.x(), p.y() + _point_size);
        vertices[3] = Point_2 (p.x() - _point_size, p.y());
      }

      if (_point_style == FIG_SQUARE || _point_style == FIG_RHOMBUS)
      {
        // Draw an empty rectangular shape (use a solid line with width 1).
        _write_polygon (4, vertices.begin(), vertices.end(),
                        _color, 1, FIG_SOLID, _style_value,
                        FIG_WHITE, FIG_NOT_FILLED);
      }
      else
      {
        // Draw a filled rectangular shape.
        _write_polygon (4, vertices.begin(), vertices.end(),
                        _color, 1, FIG_SOLID, _style_value,
                        _color, FIG_FILLED);
      }

      break;
      }
    }

    return;
  }

  /*!
   * Write a segment.
   */
  void write_segment (const Segment_2& seg)
  {
    CGAL_precondition (_ofile.is_open());

    // Clip the ray using the bounding rectangle.
    CGAL::Object    obj = intersect_func (_bound_rect, seg);
    Segment_2       clipped_seg;

    // Draw only the clipped segment (draw nothing if the segment does not
    // intersect the bounding rectangle).
    if (CGAL::assign (clipped_seg, obj))
    {
      _write_segment (clipped_seg,
                      _color, _line_width, _line_style, _style_value,
                      false);
    }

    return;
  }

  /*!
   * Write a ray.
   */
  void write_ray (const Ray_2& ray)
  {
    CGAL_precondition (_ofile.is_open());

    // Clip the ray using the bounding rectangle.
    CGAL::Object    obj = intersect_func (_bound_rect, ray);
    Segment_2       seg;

    // Draw only the clipped segment (draw nothing if the ray does not
    // intersect the bounding rectangle).
    if (CGAL::assign (seg, obj))
    {
      _write_segment (seg,
                      _color, _line_width, _line_style, _style_value,
                      false);
    }

    return;
  }

  /*!
   * Write a line.
   */
  void write_line (const Line_2& line)
  {
    CGAL_precondition (_ofile.is_open());

    // Clip the ray using the bounding rectangle.
    CGAL::Object    obj = intersect_func (_bound_rect, line);
    Segment_2       seg;

    // Draw only the clipped segment (draw nothing if the ray does not
    // intersect the bounding rectangle).
    if (CGAL::assign (seg, obj))
    {
      _write_segment (seg,
                      _color, _line_width, _line_style, _style_value,
                      false);
    }

    return;
  }

  /*!
   * Write a triangle.
   */
  void write_triangle (const Triangle_2& tri)
  {
    CGAL_precondition (_ofile.is_open());

    std::vector<Point_2>   vertices(3);

    vertices[0] = tri.vertex(0);
    vertices[1] = tri.vertex(1);
    vertices[2] = tri.vertex(2);

    _write_polygon (3, vertices.begin(), vertices.end(),
                    _color, _line_width, _line_style, _style_value,
                    _fill_color, _fill_style);
    return;
  }

  /*!
   * Write an iso-rectangle.
   */
  void write_rectangle (const Iso_rectangle_2& rect)
  {
    CGAL_precondition (_ofile.is_open());

    std::vector<Point_2>   vertices(4);

    vertices[0] = rect.vertex(0);
    vertices[1] = rect.vertex(1);
    vertices[2] = rect.vertex(2);
    vertices[3] = rect.vertex(3);

    _write_polygon (4, vertices.begin(), vertices.end(),
                    _color, _line_width, _line_style, _style_value,
                    _fill_color, _fill_style);
    return;
  }

  /*!
   * Write a polyline.
   * \param begin An iterator of the control points (of type Point_2).
   * \param end A past-the-end iterator for the control points.
   */
  template <class Input_iterator>
  void write_polyline (const Input_iterator& begin, const Input_iterator& end)
  {
    CGAL_precondition (_ofile.is_open());

    _write_polyline (begin, end,
                     _color, _line_width, _line_style, _style_value, false);
    return;
  }

  /*!
   * Write a polygon.
   */
  void write_polygon (const Polygon_2& pgn)
  {
    CGAL_precondition (_ofile.is_open());

    _write_polygon (pgn.size(),
                    pgn.vertices_begin(), pgn.vertices_end(),
                    _color, _line_width, _line_style, _style_value,
                    _fill_color, _fill_style);
    return;
  }

  /*!
   * Write a circle.
   */
  void write_circle (const Circle_2& circ)
  {
    CGAL_precondition (_ofile.is_open());

    _write_ellipse (circ.center(),
                    circ.squared_radius(), circ.squared_radius(),
                    _color, _line_width, _line_style, _style_value,
                    _fill_color, _fill_style);
    return;
  }

  /*!
   * Write a canonical ellipse.
   * \param center The center of the ellipse (the intersection of its axes).
   * \param r1_squared The squared length of the axis parallel to the x-axis.
   * \param r2_squared The squared length of the axis parallel to the y-axis.
   */
  void write_ellipse (const Point_2& center,
                      const NT& r1_squared, const NT& r2_squared)
  {
    CGAL_precondition (_ofile.is_open());

    _write_ellipse (center,
                    r1_squared, r2_squared,
                    _color, _line_width, _line_style, _style_value,
                    _fill_color, _fill_style);
    return;
  }

  /*!
   * Write a circular arc.
   * \param p1 The source point of the arc.
   * \param p2 A midpoint on the arc.
   * \param p3 The target point of the arc.
   * \pre The three points are not collinear.
   */
  void write_circular_arc (const Point_2& p1,
                           const Point_2& p2,
                           const Point_2& p3)
  {
    CGAL_precondition (_ofile.is_open());
    
    _write_arc (p1, p2, p3,
                _color, _line_width, _line_style, _style_value);
    return;
  }

  /*!
   * Write a spline.
   * \param begin An iterator of the control points (of type Point_2).
   * \param end A past-the-end iterator for the control points.
   * \param factor A shape factor for the spline: A value in the range [-1,1],
   *               where negative values are used for interpolated splines
   *               and positive values for approximated splines.
   *               The default value if 1.
   */
  template <class Input_iterator>
  void write_spline (const Input_iterator& begin, const Input_iterator& end,
                     const float& factor = 1)
  {
    CGAL_precondition (_ofile.is_open());

    if (begin == end)
      return;

    // Normalize the shape factor.
    float shape_factor;

    if (factor > 1)
      shape_factor = 1;
    else if (factor < -1)
      shape_factor = -1;
    else
      shape_factor = factor;

    _write_spline (begin, end, shape_factor,
                   _color, _line_width, _line_style, _style_value);
  }

  /*!
   * Write a text box.
   * \param pos The lower-left corner of the text box.
   * \param text The text to write.
   * \param angle The angle (in radians) that the text forms with the x-axis
   *              (0 by default).
   */
  void write_text (const Point_2& pos,
		   const char *text,
		   const double& angle = 0)
  {
    CGAL_precondition (_ofile.is_open());

    if (text == NULL || strlen(text) == 0)
      return;

    _write_text (pos, 
		 reinterpret_cast<const unsigned char*>(text), strlen(text),
		 angle,
		 _color, _font, _font_size);
    return;
  }
  //@}

  /// \name Setting the draw properties via the << operator.
  //@{

  /*!
   * Set the depth.
   */
  Fig_stream& operator<< (const Fig_depth& depth)
  {
    set_depth (static_cast<int>(depth));
    return (*this);
  }
  
  /*!
   * Set the color.
   */
  Fig_stream& operator<< (const Fig_color& color)
  {
    set_color (color);
    return (*this);
  }

  /*!
   * Set the line style.
   */
  Fig_stream& operator<< (const Fig_line_style& style)
  {
    set_line_style (style);
    return (*this);
  }

  /*!
   * Set the fill style.
   */
  Fig_stream& operator<< (const Fig_fill_style& style)
  {
    set_fill_style (style);
    return (*this);
  }

  /*!
   * Set the point style.
   */
  Fig_stream& operator<< (const Fig_point_style& style)
  {
    set_point_style (style);
    return (*this);
  }

  /*!
   * Set the arrow drawing mode. This mode will be applied when drawing
   * segments, polylines, circular arcs or splines.
   */
  Fig_stream& operator<< (const Fig_arrow_mode& mode)
  {
    set_arrow_mode (mode);
    return (*this);
  }

  /*!
   * Set the arrow type.
   */
  Fig_stream& operator<< (const Fig_arrow_type& type)
  {         
    set_arrow_type (type);
    return (*this);
  }

  /*!
   * Set the font.
   */
  Fig_stream& operator<< (const Fig_font& font)
  {
    set_font (font);
    return (*this);
  }
  //@}

  /// \name Drawing objects via the << operator.
  //@{

  /*!
   * Write a point.
   */
  Fig_stream& operator<< (const Point_2& p)
  {
    write_point (p);
    return (*this);
  }

  /*!
   * Write a line segment.
   */
  Fig_stream& operator<< (const Segment_2& seg)
  {
    write_segment (seg);
    return (*this);
  }

  /*!
   * Write a ray.
   */
  Fig_stream& operator<< (const Ray_2& ray)
  {
    write_ray (ray);
    return (*this);
  }

  /*!
   * Write a line.
   */
  Fig_stream& operator<< (const Line_2& line)
  {
    write_line (line);
    return (*this);
  }
  
  /*!
   * Write a triangle.
   */
  Fig_stream& operator<< (const Triangle_2& tri)
  {
    write_triangle (tri);
    return (*this);
  }

  /*!
   * Write a rectangle.
   */
  Fig_stream& operator<< (const Iso_rectangle_2& rect)
  {
    write_rectangle (rect);
    return (*this);
  }

  /*!
   * Write a polygon.
   */
  Fig_stream& operator<< (const Polygon_2& pgn)
  {
    write_polygon (pgn);
    return (*this);
  }

  /*!
   * Write a circle.
   */
  Fig_stream& operator<< (const Circle_2& circ)
  {
    write_circle (circ);
    return (*this);
  }
  //@}

protected:

  /*!
   * Convert a point to FIG units.
   */
  void _convert_point (const Point_2& p,
                       int& ix, int& iy) const
  {
    ix = static_cast<int> (_scale * 
                           CGAL::to_double(p.x() - _bound_rect.xmin()));
    iy = static_cast<int> (_scale * 
                           CGAL::to_double( _bound_rect.ymax() - p.y()));
    return;
  }

  /*!
   * Write a segment.
   */
  void _write_segment (const Segment_2& seg,
                       const Fig_color&      line_color,
                       const int&            line_width,
                       const Fig_line_style& line_style,
                       const double&         style_value,
                       const bool&      draw_arrows)
  {
    // Convert the segment to a polyline with two points and write it.
    std::vector<Point_2>   points (2);

    points[0] = seg.source();
    points[1] = seg.target();

    _write_polyline (points.begin(), points.end(),
                     line_color, line_width, line_style, style_value,
                     draw_arrows);
    return;
  }

  /*!
   * Write a polyline.
   */
  template <class Input_iterator>
  void _write_polyline (Input_iterator begin, Input_iterator end,
                        const Fig_color&      line_color,
                        const int&            line_width,
                        const Fig_line_style& line_style,
                        const double&         style_value,
                        const bool&      draw_arrows)
  {
    // Check if we should draw arrows.
    bool    forward_arrow = false;
    bool    backward_arrow = false;

    if (draw_arrows)
    {
      forward_arrow = (_arrow_mode == FIG_FORWARD_ARROW) ||
                      (_arrow_mode == FIG_BOTH_ARROWS);
      backward_arrow = (_arrow_mode == FIG_BACKWARD_ARROW) ||
                       (_arrow_mode == FIG_BOTH_ARROWS);
    }

    // Count the number of points in the spline.
    int     n_points = std::distance (begin, end);

    // Write the segment properties.
    _ofile << "2 1 "                      // Desginate a polyline.
           << line_style << ' ' 
           << line_width << ' ' 
           << line_color << ' '
           << FIG_WHITE << ' '            // Fill color (dummy).
           << _depth << ' ' 
           << "0 "                        // Pen style (not in use, always 0).
           << FIG_NOT_FILLED << ' '
           << style_value << ' '
           << "0 "                        // Join style (always 0).
           << "0 "                        // Cap style (always 0).
           << "-1 "                       // Radius (not in use for lines).
           << forward_arrow << ' '
           << backward_arrow << ' '
           << n_points << std::endl;

    // Write the points defining the polyline.
    bool             is_first = true;
    int              ix, iy;

    while (begin != end)
    {
      if (is_first)
      {
        _ofile << '\t';
        is_first = false;
      }
      else
      {
        _ofile << ' ';
      }

      _convert_point (*begin, ix, iy);
      _ofile << ix << ' ' << iy;

      begin++;
    }
    _ofile << std::endl;

    // Write the arrows, if necessary.
    if (forward_arrow)
        _write_arrow_line ();

    if (backward_arrow)
        _write_arrow_line ();

    return;
  }

  /*!
   * Write a polygon, reprsented as a range of points.
   */
  template <class Input_iterator>
  void _write_polygon (const int n_points,
                       Input_iterator begin, Input_iterator end,
                       const Fig_color&      line_color,
                       const int             line_width,
                       const Fig_line_style& line_style,
                       const double&         style_value,
                       const Fig_color&      fill_color,
                       const Fig_fill_style& fill_style)
  {
    // Write the polyline properties.
    _ofile << "2 3 "                      // Desginate a polygon.
           << line_style << ' ' 
           << line_width << ' ' 
           << line_color << ' '
           << fill_color << ' '
           << _depth << ' ' 
           << "0 "                        // Pen style (not in use, always 0).
           << fill_style << ' '
           << style_value << ' '
           << "0 "                        // Join style (always 0).
           << "0 "                        // Cap style (always 0).
           << "-1 "                       // Radius (not in use for lines).
           << "0 "                        // No forward arrow.
           << "0 "                        // No backward arrow.
           << n_points + 1 << std::endl;

    // Write the points.
    bool             is_first = true;
    Point_2          first = *begin;
    int              ix, iy;

    while (begin != end)
    {
      if (is_first)
      {
        _ofile << '\t';
        is_first = false;
      }
      else
      {
        _ofile << ' ';
      }

      _convert_point (*begin, ix, iy);
      _ofile << ix << ' ' << iy;

      begin++;
    }

    // Write the first point again.
    _convert_point (first, ix, iy);
    _ofile << ' ' << ix << ' ' << iy << std::endl;

    return;
  }

  /*!
   * Write an ellipse.
   */
  void _write_ellipse (const Point_2& center,
                       const NT& squared_radius_x,
                       const NT& squared_radius_y,
                       const Fig_color&      line_color,
                       const int&            line_width,
                       const Fig_line_style& line_style,
                       const double&         style_value,
                       const Fig_color&      fill_color,
                       const Fig_fill_style& fill_style) 
  {
    

    // Write the ellipse properties.
    _ofile << "1 1 "                      // Desginate an ellipse.
           << line_style << ' ' 
           << line_width << ' ' 
           << line_color << ' '
           << fill_color << ' '
           << _depth << ' ' 
           << "0 "                        // Pen style (not in use, always 0).
           << fill_style << ' '
           << style_value << ' '
           << "1 "                        // Direction (always 1).
           << "0.000 ";                   // Angle (in radians).

    // Write the center point.
    int     ix, iy;

    _convert_point (center, ix, iy);
    _ofile << ' ' << ix << ' ' << iy;

    // Write the radii.
    int  rx = static_cast<int> (_scale * 
                                std::sqrt(CGAL::to_double(squared_radius_x)));
    int  ry = static_cast<int> (_scale * 
                                std::sqrt(CGAL::to_double(squared_radius_y)));

    _ofile << ' ' << rx << ' ' << ry;

    // Write the start point (the center) and the end point (one corner of
    // the rectangles bounding the ellipse).
    _ofile << ' ' << ix << ' ' << iy
           << ' ' << (ix + rx) << ' ' << (iy +ry) << std::endl;

    return;
  }

  /*!
   * Write an arc.
   */
  void _write_arc (const Point_2& p1, const Point_2& p2, const Point_2& p3,
                   const Fig_color&      line_color,
                   const int&            line_width,
                   const Fig_line_style& line_style,
                   const double&         style_value)
  {
    // Check if we should draw arrows.
    bool    forward_arrow;
    bool    backward_arrow;

    forward_arrow = (_arrow_mode == FIG_FORWARD_ARROW) ||
                  (_arrow_mode == FIG_BOTH_ARROWS);
    backward_arrow = (_arrow_mode == FIG_BACKWARD_ARROW) ||
                   (_arrow_mode == FIG_BOTH_ARROWS);

    // Construct the supporting circle of the arc and use its center and
    // orientation.
    Circle_2   circ (p1, p2, p3);
    int        orient = (circ.orientation() == CGAL::CLOCKWISE) ? 0 : 1;

    // Write the arc properties.
    _ofile << "5 1 "                      // Desginate an open arc.
           << line_style << ' ' 
           << line_width << ' ' 
           << line_color << ' '
           << FIG_WHITE << ' '            // Fill color (dummy).
           << _depth << ' ' 
           << "0 "                        // Pen style (not in use, always 0).
           << FIG_NOT_FILLED << ' '
           << style_value << ' '
           << "0 "                        // Cap style (always 0).
           << orient << ' '
           << forward_arrow << ' '
           << backward_arrow << ' ';

    // Write the center of the circle.
    int     ix, iy;

    _convert_point (circ.center(), ix, iy);
    _ofile << ix << ' ' << iy;

    // Write the three points defining the arc.
    _convert_point (p1, ix, iy);
    _ofile << ' ' << ix << ' ' << iy;
    
    _convert_point (p2, ix, iy);
    _ofile << ' ' << ix << ' ' << iy;
    
    _convert_point (p3, ix, iy);
    _ofile << ' ' << ix << ' ' << iy << std::endl;

    // Write the arrows, if necessary.
    if (forward_arrow)
        _write_arrow_line ();

    if (backward_arrow)
        _write_arrow_line ();

    return;
  }

  /*!
   * Write a spline
   */
  template <class Input_iterator>
  void _write_spline (Input_iterator begin, Input_iterator end,
                      const float& factor,
                      const Fig_color&      line_color,
                      const int&            line_width,
                      const Fig_line_style& line_style,
                      const double&         style_value)
  {
    // Check if we should draw arrows.
    bool    forward_arrow;
    bool    backward_arrow;

    forward_arrow = (_arrow_mode == FIG_FORWARD_ARROW) ||
                  (_arrow_mode == FIG_BOTH_ARROWS);
    backward_arrow = (_arrow_mode == FIG_BACKWARD_ARROW) ||
                   (_arrow_mode == FIG_BOTH_ARROWS);

    // Count the number of points in the spline.
    int     n_points = std::distance (begin, end);

    // Write the spline properties.
    _ofile << "3 0 "                      // Desginate an open spline.
           << line_style << ' ' 
           << line_width << ' ' 
           << line_color << ' '
           << FIG_WHITE << ' '            // Fill color (dummy).
           << _depth << ' ' 
           << "0 "                        // Pen style (not in use, always 0).
           << FIG_NOT_FILLED << ' '
           << style_value << ' '
           << "0 "                        // Cap style (always 0).
           << forward_arrow << ' '
           << backward_arrow << ' '
           << n_points << std::endl;

    // Write the points defining the spline.
    bool             is_first = true;
    int              ix, iy;

    while (begin != end)
    {
      if (is_first)
      {
        _ofile << '\t';
        is_first = false;
      }
      else
      {
        _ofile << ' ';
      }

      _convert_point (*begin, ix, iy);
      _ofile << ix << ' ' << iy;

      begin++;
    }
    _ofile << std::endl;

    // Write the shape factors: 0 for the endpoints and (factor) for each
    // of the midpoints.
    int     i;

    _ofile << '\t' << "0.000";
    for (i = 0; i < n_points - 2; i++)
      _ofile << ' ' << factor;
    _ofile << '\t' << "0.000" << std::endl;

    // Write the arrows, if necessary.
    if (forward_arrow)
        _write_arrow_line ();

    if (backward_arrow)
        _write_arrow_line ();

    return;
  }

  /*!
   * Write an arrow line.
   */
  void _write_arrow_line ()
  {
    int  width  = static_cast<int> (_scale * CGAL::to_double(_arrow_width));
    int  height = static_cast<int> (_scale * CGAL::to_double(_arrow_height));

    _ofile << _arrow_type << ' '
           << "0 "                        // Arrow style (always 0). 
           << _line_width << ' ' 
           << width << ' '
           << height << std::endl;
    return;
  }

  /*!
   * Write a text box.
   */
  void _write_text (const Point_2& pos,
		    const unsigned char *text,
		    const int&          len_text,
		    const double& angle,
		    const Fig_color& font_color,
		    const Fig_font&  font,
		    const int&       font_size)
  {
    // Compute the text-box dimensions.
    const int    text_height = font_size * 1200 / 80;
    const int    text_width = len_text * font_size * 1200 / 160;

    // Write the text properties.
    _ofile << "4 0 "                      // Desginate left-justified text.
	   << font_color << ' '
	   << _depth << ' '
	   << "0 "                        // Pen style (not in use, always 0).
	   << font << ' '
	   << font_size << ' '
	   << angle << ' '
	   << "2 "                        // Indicates a special LaTeX font.
	   << text_height << ' '
	   << text_width;
    
    // Write the position coordinates.
    int     ix, iy;

    _convert_point (pos, ix, iy);
    _ofile << ' ' << ix << ' ' << iy << ' ';

    // Write the text.
    char    oct[10];
    int     i;

    for (i = 0; i < len_text; i++)
    {
      if (text[i] >= ' ' && text[i] < 128 && text[i] != '\\')
      {
	// If the current character is printable, just write it.
	_ofile << static_cast<char>(text[i]);
      }
      else
      {
	// Convert the current character to an octal string and write it.
	sprintf (oct, "\\%03o", text[i]);
	_ofile << oct;
      }
    }
   
    // Write the end-of-string sequence.
    _ofile << "\\001" << std::endl;

    return;	
  }

  /*!
   * Reset all user-defined colors.
   */
  void _reset_colors ()
  {
    int     i;

    for (i = 0; i < FIG_FIRST_USER_DEFINED_COLOR; i++)
        colors[i] = true;

    for (i = FIG_FIRST_USER_DEFINED_COLOR;
         i < FIG_LAST_USER_DEFINED_COLOR; i++)
    {
        colors[i] = false;
    }

    return;
  }
};

} //namespace CGAL

#endif
