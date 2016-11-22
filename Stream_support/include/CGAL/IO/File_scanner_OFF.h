// Copyright (c) 1997  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
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
// Author(s)     : Lutz Kettner  <kettner@mpi-sb.mpg.de>

#ifndef CGAL_IO_FILE_SCANNER_OFF_H
#define CGAL_IO_FILE_SCANNER_OFF_H 1

#include <CGAL/config.h>
#include <cstddef>
#include <CGAL/IO/binary_file_io.h>
#include <CGAL/IO/File_header_OFF.h>
#include <iostream>
#include <sstream>
#include <boost/cstdint.hpp>

#include <CGAL/Point_3.h>
#include <CGAL/Vector_3.h>

namespace CGAL {

class CGAL_EXPORT File_scanner_OFF : public File_header_OFF {
    std::istream&  m_in;
    bool           normals_read;
    void skip_comment() { m_in >> skip_comment_OFF; }
public:
    File_scanner_OFF( std::istream& in, bool verbose = false)
      : File_header_OFF(verbose), m_in(in), normals_read(false) {
        in >> static_cast<File_header_OFF&>( *this);
    }
    File_scanner_OFF( std::istream& in, const File_header_OFF& header)
      : File_header_OFF(header), m_in(in), normals_read(false) {}

    std::istream& in() { return m_in; }

    // The scan_vertex() routine is provided for multiple
    // coordinate types to support parameterized polytopes.
    void scan_vertex( float&  x, float&  y, float&  z) {
        if ( binary()) {
            I_Binary_read_big_endian_float32( m_in, x);
            I_Binary_read_big_endian_float32( m_in, y);
            I_Binary_read_big_endian_float32( m_in, z);
            if ( is_homogeneous()) {
                float w;
                I_Binary_read_big_endian_float32( m_in, w);
                x /= w;
                y /= w;
                z /= w;
            }
        } else {
            skip_comment();
            m_in >> iformat(x) >> iformat(y) >> iformat(z);
            if ( is_homogeneous()) {
                float w;
                m_in >> iformat(w);
                x /= w;
                y /= w;
                z /= w;
            }
        }
    }
    void scan_vertex( double& x, double& y, double& z) {
        if ( binary()) {
            float f;
            I_Binary_read_big_endian_float32( m_in, f);
            x = f;
            I_Binary_read_big_endian_float32( m_in, f);
            y = f;
            I_Binary_read_big_endian_float32( m_in, f);
            z = f;
            if ( is_homogeneous()) {
                I_Binary_read_big_endian_float32( m_in, f);
                x /= f;
                y /= f;
                z /= f;
            }
        } else {
            skip_comment();
            m_in >> iformat(x) >> iformat(y) >> iformat(z);
            if ( is_homogeneous()) {
                double w;
                m_in >> iformat(w);
                x /= w;
                y /= w;
                z /= w;
            }
        }
    }
    void scan_vertex( int& x, int& y, int& z) {
        if ( binary()) {
            float fx, fy, fz;
            I_Binary_read_big_endian_float32( m_in, fx);
            I_Binary_read_big_endian_float32( m_in, fy);
            I_Binary_read_big_endian_float32( m_in, fz);
            if ( is_homogeneous()) {
                float fw;
                I_Binary_read_big_endian_float32( m_in, fw);
                x = int( fx / fw);
                y = int( fy / fw);
                z = int( fz / fw);
            } else {
                x = int(fx);
                y = int(fy);
                z = int(fz);
            }
        } else {
            skip_comment();
            if ( is_homogeneous()) {
                double fx, fy, fz, fw;
                m_in >> iformat(fx) >> iformat(fy) >> iformat(fz) >> iformat(fw);
                x = int( fx / fw);
                y = int( fy / fw);
                z = int( fz / fw);
            } else {
                double d;
                m_in >> iformat(d);
                x = int(d);
                m_in >> iformat(d);
                y = int(d);
                m_in >> iformat(d);
                z = int(d);
            }
        }
    }

    void scan_vertex( float&  x, float&  y, float&  z, float&  w) {
        w = 1;
        if ( binary()) {
            I_Binary_read_big_endian_float32( m_in, x);
            I_Binary_read_big_endian_float32( m_in, y);
            I_Binary_read_big_endian_float32( m_in, z);
            if ( is_homogeneous())
                I_Binary_read_big_endian_float32( m_in, w);
        } else {
            skip_comment();
            m_in >> iformat(x) >> iformat(y) >> iformat(z);
            if ( is_homogeneous())
              m_in >> iformat(w);
        }
    }
    void scan_vertex( double& x, double& y, double& z, double& w) {
        w = 1;
        if ( binary()) {
            float f;
            I_Binary_read_big_endian_float32( m_in, f);
            x = f;
            I_Binary_read_big_endian_float32( m_in, f);
            y = f;
            I_Binary_read_big_endian_float32( m_in, f);
            z = f;
            if ( is_homogeneous()) {
                I_Binary_read_big_endian_float32( m_in, f);
                w = f;
            }
        } else {
            skip_comment();
            m_in >> iformat(x);
            m_in >> iformat(y);
            m_in >> iformat(z);
            if ( is_homogeneous())
              m_in >> iformat(w);
        }
    }
    void scan_vertex( int& x, int& y, int& z, int& w) {
        w = 1;
        if ( binary()) {
            float f;
            I_Binary_read_big_endian_float32( m_in, f);
            x = int(f);
            I_Binary_read_big_endian_float32( m_in, f);
            y = int(f);
            I_Binary_read_big_endian_float32( m_in, f);
            z = int(f);
            if ( is_homogeneous()) {
                I_Binary_read_big_endian_float32( m_in, f);
                w = int(f);
            }
        } else {
            skip_comment();
            double d;
            m_in >> iformat(d);
            x = int(d);
            m_in >> iformat(d);
            y = int(d);
            m_in >> iformat(d);
            z = int(d);
            if ( is_homogeneous()) {
                m_in >> iformat(d);
                w = int(d);
            }
        }
    }

    void scan_normal( float&  x, float&  y, float&  z) {
        if ( has_normals()) {
            normals_read = true;
            if ( binary()) {
                I_Binary_read_big_endian_float32( m_in, x);
                I_Binary_read_big_endian_float32( m_in, y);
                I_Binary_read_big_endian_float32( m_in, z);
                if ( is_homogeneous()) {
                    float w;
                    I_Binary_read_big_endian_float32( m_in, w);
                    x /= w;
                    y /= w;
                    z /= w;
                }
            } else {
                m_in >> iformat(x) >> iformat(y) >> iformat(z);
                if ( is_homogeneous()) {
                    float w;
                    m_in >> iformat(w);
                    x /= w;
                    y /= w;
                    z /= w;
                }
            }
        }
    }
    void scan_normal( double& x, double& y, double& z) {
        if ( has_normals()) {
            normals_read = true;
            if ( binary()) {
                float fx, fy, fz;
                I_Binary_read_big_endian_float32( m_in, fx);
                I_Binary_read_big_endian_float32( m_in, fy);
                I_Binary_read_big_endian_float32( m_in, fz);
                if ( is_homogeneous()) {
                    float fw;
                    I_Binary_read_big_endian_float32( m_in, fw);
                    x = fx / fw;
                    y = fy / fw;
                    z = fz / fw;
                } else {
                    x = fx;
                    y = fy;
                    z = fz;
                }
            } else {
                if ( is_homogeneous()) {
                    float fx, fy, fz, fw;
                    m_in >> iformat(fx) >> iformat(fy) >> iformat(fz) >> iformat(fw);
                    x = fx / fw;
                    y = fy / fw;
                    z = fz / fw;
                } else
                    m_in >> iformat(x) >> iformat(y) >> iformat(z);
            }
        }
    }
    void scan_normal( int& x, int& y, int& z) {
        if ( has_normals()) {
            normals_read = true;
            if ( binary()) {
                float fx, fy, fz;
                I_Binary_read_big_endian_float32( m_in, fx);
                I_Binary_read_big_endian_float32( m_in, fy);
                I_Binary_read_big_endian_float32( m_in, fz);
                if ( is_homogeneous()) {
                    float fw;
                    I_Binary_read_big_endian_float32( m_in, fw);
                    x = int( fx / fw);
                    y = int( fy / fw);
                    z = int( fz / fw);
                } else {
                    x = int(fx);
                    y = int(fy);
                    z = int(fz);
                }
            } else {
                if ( is_homogeneous()) {
                    float fx, fy, fz, fw;
                    m_in >> iformat(fx) >> iformat(fy) >> iformat(fz) >> iformat(fw);
                    x = int( fx / fw);
                    y = int( fy / fw);
                    z = int( fz / fw);
                } else {
                    double d;
                    m_in >> iformat(d);
                    x = int(d);
                    m_in >> iformat(d);
                    y = int(d);
                    m_in >> iformat(d);
                    z = int(d);
                }
            }
        }
    }

    void scan_normal( float&  x, float&  y, float&  z, float&  w) {
        w = 1;
        if ( has_normals()) {
            normals_read = true;
            if ( binary()) {
                I_Binary_read_big_endian_float32( m_in, x);
                I_Binary_read_big_endian_float32( m_in, y);
                I_Binary_read_big_endian_float32( m_in, z);
                if ( is_homogeneous())
                    I_Binary_read_big_endian_float32( m_in, w);
            } else {
                m_in >> iformat(x) >> iformat(y) >> iformat(z);
                if ( is_homogeneous())
                    m_in >> iformat(w);
            }
        }
    }
    void scan_normal( double& x, double& y, double& z, double& w) {
        w = 1;
        if ( has_normals()) {
            normals_read = true;
            if ( binary()) {
                float f;
                I_Binary_read_big_endian_float32( m_in, f);
                x = f;
                I_Binary_read_big_endian_float32( m_in, f);
                y = f;
                I_Binary_read_big_endian_float32( m_in, f);
                z = f;
                if ( is_homogeneous()) {
                    I_Binary_read_big_endian_float32( m_in, f);
                    w = f;
                }
            } else {
                m_in >> iformat(x) >> iformat(y) >> iformat(z);
                if ( is_homogeneous())
                    m_in >> iformat(w);
            }
        }
    }
    void scan_normal( int& x, int& y, int& z, int& w) {
        w = 1;
        if ( has_normals()) {
            normals_read = true;
            if ( binary()) {
                float f;
                I_Binary_read_big_endian_float32( m_in, f);
                x = int(f);
                I_Binary_read_big_endian_float32( m_in, f);
                y = int(f);
                I_Binary_read_big_endian_float32( m_in, f);
                z = int(f);
                if ( is_homogeneous()) {
                    I_Binary_read_big_endian_float32( m_in, f);
                    w = int(f);
                }
            } else {
                double d;
                m_in >> iformat(d);
                x = int(d);
                m_in >> iformat(d);
                y = int(d);
                m_in >> iformat(d);
                z = int(d);
                if ( is_homogeneous()) {
                    m_in >> iformat(d);
                    w = int(d);
                }
            }
        }
    }

    static const Color& get_indexed_color(int id)
    {
      static const Color color[149] = {
        Color(255, 255, 255, 191),
        Color(255, 255, 255, 191),
        Color(255, 255, 255, 191),
        Color(255, 255, 255, 191),
        Color(255, 255, 255, 191),
        Color(255, 255, 255, 191),
        Color(178, 38, 25, 191),
        Color(51, 51, 204, 191),
        Color(229, 153, 5, 191),
        Color(25, 76, 204, 191),
        Color(25, 178, 51, 191),
        Color(204, 204, 102, 191),
        Color(178, 178, 0, 191),
        Color(178, 0, 178, 191),
        Color(0, 178, 178, 191),
        Color(229, 0, 51, 191),
        Color(51, 229, 0, 191),
        Color(0, 51, 229, 191),
        Color(191, 191, 191, 191),
        Color(204, 102, 0, 191),
        Color(204, 102, 0, 191),
        Color(0, 102, 204, 191),
        Color(0, 102, 204, 191),
        Color(0, 204, 102, 191),
        Color(0, 204, 102, 191),
        Color(102, 0, 204, 191),
        Color(102, 0, 204, 191),
        Color(204, 0, 102, 191),
        Color(204, 0, 102, 191),
        Color(178, 127, 51, 191),
        Color(178, 127, 51, 191),
        Color(178, 178, 0, 191),
        Color(178, 0, 178, 191),
        Color(0, 178, 178, 191),
        Color(229, 0, 0, 191),
        Color(0, 229, 0, 191),
        Color(0, 0, 229, 191),
        Color(191, 191, 191, 191),
        Color(204, 102, 0, 191),
        Color(102, 204, 0, 191),
        Color(0, 102, 204, 191),
        Color(0, 204, 102, 191),
        Color(102, 0, 204, 191),
        Color(204, 0, 102, 191),
        Color(178, 178, 0, 191),
        Color(178, 0, 178, 191),
        Color(0, 178, 178, 191),
        Color(229, 0, 0, 191),
        Color(0, 229, 0, 191),
        Color(0, 0, 229, 191),
        Color(191, 191, 191, 191),
        Color(204, 102, 0, 191),
        Color(102, 204, 0, 191),
        Color(0, 102, 204, 191),
        Color(0, 204, 102, 191),
        Color(102, 0, 204, 191),
        Color(204, 0, 102, 191),
        Color(178, 178, 0, 191),
        Color(178, 0, 178, 191),
        Color(0, 178, 178, 191),
        Color(229, 0, 0, 191),
        Color(0, 229, 0, 191),
        Color(0, 0, 229, 191),
        Color(191, 191, 191, 191),
        Color(204, 102, 0, 191),
        Color(102, 204, 0, 191),
        Color(0, 102, 204, 191),
        Color(0, 204, 102, 191),
        Color(102, 0, 204, 191),
        Color(204, 0, 102, 191),
        Color(255, 255, 255, 191),
        Color(255, 255, 255, 191),
        Color(255, 255, 255, 191),
        Color(255, 255, 255, 191),
        Color(255, 255, 255, 191),
        Color(255, 255, 255, 191),
        Color(12, 76, 25, 191),
        Color(178, 2, 25, 191),
        Color(51, 12, 153, 191),
        Color(229, 229, 5, 191),
        Color(0, 51, 102, 191),
        Color(25, 102, 102, 191),
        Color(204, 204, 204, 191),
        Color(178, 178, 0, 191),
        Color(178, 178, 0, 191),
        Color(178, 0, 178, 191),
        Color(178, 0, 178, 191),
        Color(0, 178, 178, 191),
        Color(0, 178, 178, 191),
        Color(229, 0, 0, 191),
        Color(229, 0, 0, 191),
        Color(0, 229, 0, 191),
        Color(0, 229, 0, 191),
        Color(0, 0, 229, 191),
        Color(0, 0, 229, 191),
        Color(191, 191, 191, 191),
        Color(191, 191, 191, 191),
        Color(204, 102, 0, 191),
        Color(204, 102, 0, 191),
        Color(0, 102, 204, 191),
        Color(0, 102, 204, 191),
        Color(0, 204, 102, 191),
        Color(0, 204, 102, 191),
        Color(102, 0, 204, 191),
        Color(102, 0, 204, 191),
        Color(204, 0, 102, 191),
        Color(204, 0, 102, 191),
        Color(178, 127, 51, 191),
        Color(178, 127, 51, 191),
        Color(178, 178, 0, 191),
        Color(178, 0, 178, 191),
        Color(0, 178, 178, 191),
        Color(229, 0, 0, 191),
        Color(0, 229, 0, 191),
        Color(0, 0, 229, 191),
        Color(191, 191, 191, 191),
        Color(204, 102, 0, 191),
        Color(102, 204, 0, 191),
        Color(0, 102, 204, 191),
        Color(0, 204, 102, 191),
        Color(102, 0, 204, 191),
        Color(204, 0, 102, 191),
        Color(178, 178, 0, 191),
        Color(178, 0, 178, 191),
        Color(0, 178, 178, 191),
        Color(229, 0, 0, 191),
        Color(0, 229, 0, 191),
        Color(0, 0, 229, 191),
        Color(191, 191, 191, 191),
        Color(204, 102, 0, 191),
        Color(102, 204, 0, 191),
        Color(0, 102, 204, 191),
        Color(0, 204, 102, 191),
        Color(102, 0, 204, 191),
        Color(204, 0, 102, 191),
        Color(178, 178, 0, 191),
        Color(178, 0, 178, 191),
        Color(0, 178, 178, 191),
        Color(229, 0, 0, 191),
        Color(0, 229, 0, 191),
        Color(0, 0, 229, 191),
        Color(191, 191, 191, 191),
        Color(204, 102, 0, 191),
        Color(102, 204, 0, 191),
        Color(0, 102, 204, 191),
        Color(0, 204, 102, 191),
        Color(102, 0, 204, 191),
        Color(204, 0, 102, 191),
        Color(120, 120, 120, 120) };
      if(id > 148) id =148;
      return color[id];
    }

    static CGAL::Color get_color_from_line(std::istream &is)
    {

        std::string color_info;
        bool is_float = false;

        std::string col;
        //get the line content
        std::getline(is, col);
        //split it into strings
        std::istringstream iss(col);
        //holds the rgb values
        unsigned char rgb[3];
        int index =0;
        //split the string into numbers
        while(iss>>color_info){
            //stop if comment is read
            if(color_info.at(0) == '#')
                break;
            //detect if the value is float
            for(int c = 0; c<static_cast<int>(color_info.length()); c++)
            {
                if(color_info.at(c) == '.')
                {
                    is_float = true;
                    break;
                }
            }
            //if the value is of float type, convert it into an int
            if(is_float)
                rgb[index] = static_cast<unsigned char>(atof(color_info.c_str())*255);
            //else stores the value
            else
                rgb[index] = static_cast<unsigned char>(atoi(color_info.c_str()));
            index++;
            if(index == 3)
             break;
        }
        CGAL::Color color;
        //if there were only one number, fetch the color in the color map
        if(index<2)
            color = get_indexed_color(rgb[0]);
        //else create the coor with the 3 values;
        else
            color = CGAL::Color(rgb[0], rgb[1], rgb[2]);
        return color;
    }

    void scan_color( unsigned char& r, unsigned char& g, unsigned char& b) {
        if ( binary()) {
            float fr, fg, fb;
            I_Binary_read_big_endian_float32( m_in, fr);
            I_Binary_read_big_endian_float32( m_in, fg);
            I_Binary_read_big_endian_float32( m_in, fb);
            r = (unsigned char)(fr);
            g = (unsigned char)(fg);
            b = (unsigned char)(fb);

        } else {
            CGAL::Color color = get_color_from_line(m_in);
            r = color.red();
            g = color.green();
            b = color.blue();
        }
    }

  void skip_to_next_vertex( std::size_t current_vertex);

  void scan_facet( std::size_t& size, std::size_t CGAL_assertion_code(current_facet)) {
        CGAL_assertion( current_facet < size_of_facets());
        if ( binary()){
            boost::int32_t i32;
            I_Binary_read_big_endian_integer32( m_in, i32);
            size = i32;
        } else {
            skip_comment();
            m_in >> size;
        }
    }

  void scan_facet_vertex_index( std::size_t& index,
                                std::size_t current_facet) {
    if ( binary()){
      boost::int32_t i32;
      I_Binary_read_big_endian_integer32( m_in, i32);
      index = i32;
    } else
      m_in >> index;

    if( m_in.fail()) {
      if ( verbose()) {
        std::cerr << " " << std::endl;
        std::cerr << "File_scanner_OFF::" << std::endl;
        std::cerr << "scan_facet_vertex_index(): input error:  "
          "cannot read OFF file beyond facet "
                  << current_facet << "." << std::endl;
      }
      set_off_header( false);
      return;
    }
    bool error  = index < index_offset();
    index -= index_offset();

    if(error || (index >= size_of_vertices())) {
      m_in.clear( std::ios::failbit);
      if ( verbose()) {
        std::cerr << " " << std::endl;
        std::cerr << "File_scanner_OFF::" << std::endl;
        std::cerr << "scan_facet_vertex_index(): input error: "
          "facet " << current_facet << ": vertex index "
                  << index + index_offset() << ": is out of range."
                  << std::endl;
      }
      set_off_header( false);
      return;
    }
  }

  void skip_to_next_facet( std::size_t current_facet);
};

template < class Point> inline
Point&
file_scan_vertex( File_scanner_OFF& scanner, Point& p) {
    typedef typename Point::R R;
    typedef typename R::RT    RT;
    double x, y, z, w;
    scanner.scan_vertex( x, y, z, w);
    if ( w == 1)
        p = Point( RT(x), RT(y), RT(z));
    else
        p = Point( RT(x), RT(y), RT(z), RT(w));
    return p;
}

template < class T_Color> inline
T_Color&
file_scan_color( File_scanner_OFF& scanner, T_Color& c) {
    unsigned char r, g, b;
    scanner.scan_color(r,g,b);
        c = T_Color(r,g,b);
    return c;
}

template < class Vector> inline
Vector&
file_scan_normal( File_scanner_OFF& scanner, Vector& v) {
    typedef typename Vector::R R;
    typedef typename R::RT     RT;
    double x, y, z, w;
    scanner.scan_normal( x, y, z, w);
    if ( w == 1)
        v = Vector( RT(x), RT(y), RT(z));
    else
        v = Vector( RT(x), RT(y), RT(z), RT(w));
    return v;
}

} //namespace CGAL

#ifdef CGAL_HEADER_ONLY
#include <CGAL/IO/File_scanner_OFF_impl.h>
#endif // CGAL_HEADER_ONLY

#endif // CGAL_IO_FILE_SCANNER_OFF_H //
// EOF //
