// Copyright (c) 2007  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Fernando Cacciola
//
// Descriptions of the file format can be found at
// http://www.autodesk.com/techpubs/autocad/acad2000/dxf/

#ifndef CGAL_DXF_STREAM_H
#define CGAL_DXF_STREAM_H

#include <CGAL/license/Straight_skeleton_2.h>


#include <CGAL/basic.h>
#include <CGAL/Polygon_2.h>

#include <CGAL/IO/Dxf_writer.h>

#include <vector>
#include <fstream>
#include <stdio.h>

namespace CGAL {

class Dxf_layer
{
public:

  Dxf_layer( std::string aStr ) : mStr(aStr) {}
  
  std::string str() const { return mStr ; }
  
private:
  
  std::string mStr ;
} ;
 
template <class Kernel_>
class Dxf_stream
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
  typedef CGAL::Polygon_2<Kernel>                Polygon_2;
  typedef typename Kernel::Circle_2              Circle_2;

protected:

  // Data members:
  Dxf_writer          mWriter ;
  int                 mDefaultDxfColor;
  int                 mDxfColor;
  Color               mCgalColor ; 
  std::string         mLayer ;
  
  struct Color_less
  {
    bool operator() ( Color const& a, Color const& b ) const
    {
      return Color_value(a) < Color_value(b);        
    }
    
    static int Color_value ( Color const& c ) 
    {
      return ( int(c.r()) << 16 ) + ( int(c.g()) << 8  ) + ( int(c.b()) ) ;  
    }
  } ;
  
  typedef std::map<Color,int,Color_less>  Color_table ;
  typedef typename Color_table::const_iterator Color_table_iterator ;  
  Color_table mColorTable ;
   
private:

  // Copy constructor and assignment operator - not supported.
  Dxf_stream (const Dxf_stream<Kernel>& );
  const Dxf_stream<Kernel>& operator= (const Dxf_stream<Kernel>& );

public:

  /// \name Constructors and destructor.
  //@{

  /*!
   * Constructor.
   * \param filename The name of the output FIG file.
   */
  Dxf_stream ( std::ostream& out )
    :
    mWriter          (out)
   ,mDefaultDxfColor (255)
   ,mDxfColor        (255)
   ,mCgalColor       (WHITE)
   ,mLayer           ("0")
  {
    setup_initial_color_table();
  }  

  /*!
   * Destructor.
   */
  virtual ~Dxf_stream () {}
  //@}

  /// \name Accessing drawing properties.
  //@{
  
  /*!
   * Get the current layer.
   */
  std::string layer() const { return mLayer ; }
  
  /*!
   * Get the current CGAL color.
   */
  Color color () const { return mCgalColor ; }

  /*!
   * Get the current DXF color.
   */
  int dxf_color () const { return mDxfColor ; }

  /*!
   * Get the current DXF color.
   */
  int default_dxf_color () const { return mDefaultDxfColor ; }
  
  /// \name Set the drawing properties.
  //@{

  /*!
   * Set the current layer.
   */
  void set_layer ( std::string aLayer ) { mLayer = aLayer ; }
  
  /*!
   * Set the current color.
   * \pre The color must be defined.
   */
  void set_color ( Color aColor )
  {
    mCgalColor = aColor ;
    
    Color_table_iterator f = mColorTable.find(aColor);
    if ( f != mColorTable.end() )
         mDxfColor = f->second ;
    else mDxfColor = mDefaultDxfColor ;
  }


  /*!
   * Sets the default DXF color in case a CGAL color is unmapped.
   * \param aDxfColor  The default DXF color.
   */
  void define_default_dxf_color ( int aDxfColor )
  {
    mDefaultDxfColor = aDxfColor ;
  }
  
  /*!
   * Adds a mapping between a CGAL Color and a DXF color.
   * \param aCgalColor The CGAL color.
   * \param aDxfColor  The DXF color.
   */
  void define_color ( Color const& aCgalColor, int aDxfColor )
  {
    mColorTable.insert( std::make_pair(aCgalColor,aDxfColor) ) ;
  }

  //@}

  /// \name Writing objects.
  //@{

  /*!
   * Write a 2D segment.
   */
  void write_segment_2 (const Segment_2& seg)
  {
    mWriter.add_segment_2( seg.source(), seg.target(), mLayer, mDxfColor ) ;
  }


  /*!
   * Write a 2D polyline.
   * \param begin An iterator of the control points (of type Point_2).
   * \param end A past-the-end iterator for the control points.
   */
  template <class Input_iterator>
  void write_polyline_2 (const Input_iterator& begin, const Input_iterator& end)
  {
    mWriter.add_polyline_2( begin, end, false, mLayer, mDxfColor ) ;
  }

  /*!
   * Write a 2D polygon (there is an added segment between the last vertex and the first)
   * \param begin An iterator of the control points (of type Point_2).
   * \param end A past-the-end iterator for the control points.
   */
  template <class Input_iterator>
  void write_polygon_2 (const Input_iterator& begin, const Input_iterator& end)
  {
    mWriter.add_polyline_2( begin, end, true, mLayer, mDxfColor ) ;
  }
  
  /*!
   * Write a 2D polyline but as a sequence of line segments
   * \param begin An iterator of the control points (of type Point_2).
   * \param end A past-the-end iterator for the control points.
   */
  template <class Input_iterator>
  void write_open_segment_chain_2 (const Input_iterator& begin, const Input_iterator& end)
  {
    mWriter.add_segments_2( begin, end, false, mLayer, mDxfColor ) ;
  }

  /*!
   * Write a 2D closed polyline but as a sequence of line segments
   * \param begin An iterator of the control points (of type Point_2).
   * \param end A past-the-end iterator for the control points.
   */
  template <class Input_iterator>
  void write_closed_segment_chain_2 (const Input_iterator& begin, const Input_iterator& end)
  {
    mWriter.add_segments_2( begin, end, true, mLayer, mDxfColor ) ;
  }
  
  /*!
   * Write a 2D (closed) polygon.
   */
  void write_polygon (const Polygon_2& pgn)
  {
    mWriter.add_polyline_2( pgn.begin(), pgn.end(), true, mLayer, mDxfColor ) ;
  }


  /// \name Setting the draw properties via the << operator.
  //@{

  /*!
   * Set the current layer.
   */
  Dxf_stream& operator<< ( Dxf_layer const& aLayer )
  {
    set_layer ( aLayer.str() );
    return (*this);
  }
  
  /*!
   * Set the current color.
   */
  Dxf_stream& operator<< ( Color const& aColor )
  {
    set_color (aColor);
    return (*this);
  }


  /// \name Drawing objects via the << operator.
  //@{

  /*!
   * Write a line segment.
   */
  Dxf_stream& operator<< (const Segment_2& seg)
  {
    write_segment_2 (seg);
    return (*this);
  }

  /*!
   * Write a polygon.
   */
  Dxf_stream& operator<< (const Polygon_2& pgn)
  {
    write_polygon_2 (pgn);
    return (*this);
  }

  //@}

protected:

  void setup_initial_color_table()
  {
    define_color(BLACK,0);
    define_color(RED,1);
    define_color(YELLOW,2);
    define_color(GREEN,3);
    define_color(PURPLE,4);
    define_color(BLUE,5);
    define_color(VIOLET,6);
    define_color(WHITE,7);
    define_color(GRAY,8);
  }
  
};

} // end namespace CGAL

#endif
