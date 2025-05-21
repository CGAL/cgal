// Copyright (c) 1997
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org);
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Lutz Kettner  <kettner@mpi-sb.mpg.de>

#ifndef CGAL_IO_OFF_FILE_SCANNER_OFF_H
#define CGAL_IO_OFF_FILE_SCANNER_OFF_H

#include <CGAL/config.h>

#include <CGAL/IO/binary_file_io.h>
#include <CGAL/IO/OFF/File_header_OFF.h>
#include <CGAL/IO/io.h>

#include <boost/cstdint.hpp>

#include <vector>
#include <cstddef>
#include <iostream>
#include <sstream>
#include <vector>

namespace CGAL {

class File_scanner_OFF
  : public File_header_OFF
{
  std::vector<double> entries;
  std::size_t color_entries;
  std::size_t first_color_index;
  std::string line;
  std::istream& m_in;
  bool normals_read;

  void skip_comment() { m_in >> skip_comment_OFF; }

public:
  File_scanner_OFF(std::istream& in, bool verbose = false)
    : File_header_OFF(verbose), m_in(in), normals_read(false)
  {
    in >> static_cast<File_header_OFF&>(*this);
  }

  File_scanner_OFF(std::istream& in, const File_header_OFF& header)
    :
      File_header_OFF(header), m_in(in), normals_read(false)
  { }

  std::istream& in() { return m_in; }

// Coordinates
  void scan_vertex(double& x, double& y, double& z, double& w)
  {
    w = 1;
    if(binary())
    {
      float f;
      I_Binary_read_big_endian_float32(m_in, f);
      x = f;
      I_Binary_read_big_endian_float32(m_in, f);
      y = f;
      I_Binary_read_big_endian_float32(m_in, f);
      z = f;
      if(is_homogeneous())
      {
        I_Binary_read_big_endian_float32(m_in, f);
        w = f;
      }
    }
    else
    {
      skip_comment();
      line.clear();
      std::getline(m_in, line);
      // First remove the comment if there is one
      std::size_t pos = line.find('#');
      if(pos != std::string::npos){
        line = line.substr(0,pos);
      }

      // Read all numbers in the line
      std::istringstream issline(line);
      entries.clear();
      double d;
      while(issline >> IO::iformat(d)){
        entries.push_back(d);
      }

      if(has_colors()){
        // Compute how many entries are there for the color
        int H = (is_homogeneous())? 1:0;
        first_color_index = 3 + H;

        color_entries = entries.size();
        color_entries -= 3 + H; // coordinates
        if(has_normals()){
          first_color_index += 3 + H;
          color_entries -= 3 + H;
        }
        if(has_textures()){
          color_entries -= 2 + H;
        }
        // now color_entries should be 0, 1, 3, or 4 for the color
        if(color_entries != 0)
          m_has_vcolors = true;
      }

      if(entries.size() < 3)
      {
        m_in.clear(std::ios::badbit);
        if(verbose())
          std::cerr<<"error while reading vertex."<<std::endl;
        return;
      }
      x = entries[0];
      y = entries[1];
      z = entries[2];

      if(is_homogeneous()){
        if(entries.size() < 4){
          m_in.clear(std::ios::badbit);
          if(verbose())
            std::cerr<<"error while reading vertex."<<std::endl;
          return;
        }
        w = entries[3];
      }
    }
  }

  void scan_vertex(float& x, float& y, float& z, float& w)
  {
    double dx(0), dy(0), dz(0), dw(0);
    scan_vertex(dx, dy, dz, dw);
    x = static_cast<float>(dx);
    y = static_cast<float>(dy);
    z = static_cast<float>(dz);
    w = static_cast<float>(dw);
  }

  void scan_vertex(int& x, int& y, int& z, int& w)
  {
    double dx(0), dy(0), dz(0), dw(0);
    scan_vertex(dx, dy, dz, dw);
    x = static_cast<int>(dx);
    y = static_cast<int>(dy);
    z = static_cast<int>(dz);
    w = static_cast<int>(dw);
  }

  void scan_vertex(double& x, double& y, double& z)
  {
    double dx(0), dy(0), dz(0), dw(0);
    scan_vertex(dx, dy, dz, dw);
    x = dx / dw;
    y = dy / dw;
    z = dz / dw;
  }

  void scan_vertex(float& x, float& y, float& z)
  {
    double dx(0), dy(0), dz(0);
    scan_vertex(dx, dy, dz);
    x = static_cast<float>(dx);
    y = static_cast<float>(dy);
    z = static_cast<float>(dz);
  }

  void scan_vertex(int& x, int& y, int& z)
  {
    double dx(0), dy(0), dz(0);
    scan_vertex(dx, dy, dz);
    x = static_cast<int>(dx);
    y = static_cast<int>(dy);
    z = static_cast<int>(dz);
  }

// Textures
  void scan_texture(double& x, double& y, double& w)
  {
    w=1;
    if(has_textures())
    {
      if(binary())
      {
        float fx, fy;
        I_Binary_read_big_endian_float32(m_in, fx);
        I_Binary_read_big_endian_float32(m_in, fy);
        if(is_homogeneous())
        {
          float fw;
          I_Binary_read_big_endian_float32(m_in, fw);
          x = fx / fw;
          y = fy / fw;
        } else
        {
          x = fx;
          y = fy;
        }
      }
      else
      {
        std::size_t first_texture_index = first_color_index + color_entries;
        x = entries[first_texture_index];
        y = entries[first_texture_index + 1];
        if(is_homogeneous())
        {
          if(entries.size() <= first_texture_index + 2)
          {
            m_in.clear(std::ios::badbit);
            if(verbose())
              std::cerr<<"error while reading texture."<<std::endl;
            return;
          }
          w = entries[first_texture_index + 2];
        }
      }

    }
    else
    {
      x=0;
      y=0;
    }
  }

  void scan_texture(float& x, float& y, float& w)
  {
    double dx, dy, dw;
    scan_texture(dx, dy, dw);
    x = static_cast<float>(dx);
    y = static_cast<float>(dy);
    w = static_cast<float>(dw);
  }

  void scan_texture(int& x, int& y, int& w)
  {
    double dx, dy, dw;
    scan_texture(dx, dy, dw);
    x = static_cast<int>(dx);
    y = static_cast<int>(dy);
    w = static_cast<int>(dw);
  }

  void scan_texture(double& x, double& y)
  {
    double dx, dy, dw;
    scan_texture(dx, dy, dw);
    x = dx / dw;
    y = dy / dw;
  }

  void scan_texture(float& x, float& y)
  {
    double dx, dy;
    scan_texture(dx, dy);
    x = static_cast<float>(dx);
    y = static_cast<float>(dy);
  }

  void scan_texture(int& x, int& y)
  {
    double dx, dy;
    scan_texture(dx, dy);
    x = static_cast<int>(dx);
    y = static_cast<int>(dy);
  }

  // Normals

  void scan_normal(double& x, double& y, double& z, double& w)
  {
    w = 1;
    if(has_normals())
    {
      normals_read = true;
      if(binary())
      {
        float f;
        I_Binary_read_big_endian_float32(m_in, f);
        x = f;
        I_Binary_read_big_endian_float32(m_in, f);
        y = f;
        I_Binary_read_big_endian_float32(m_in, f);
        z = f;
        if(is_homogeneous())
        {
          I_Binary_read_big_endian_float32(m_in, f);
          w = f;
        }
      }
      else
      {
        std::size_t first_normal_index = (is_homogeneous())? 4:3;
        if(entries.size() <= first_normal_index + 2)
        {
          m_in.clear(std::ios::badbit);
          if(verbose())
            std::cerr<<"error while reading normal."<<std::endl;
          return;
        }
        x = entries[first_normal_index];
        y = entries[first_normal_index + 1];
        z = entries[first_normal_index + 2];

        if(is_homogeneous())
        {
          if(entries.size() <= first_normal_index + 3){
            m_in.clear(std::ios::badbit);
            if(verbose())
              std::cerr<<"error while reading normal."<<std::endl;
            return;
          }
          w = entries[first_normal_index + 3];
        }
      }
    }
  }

  void scan_normal(float& x, float& y, float& z, float& w)
  {
    double dx(0), dy(0), dz(0), dw(0);
    scan_normal(dx, dy, dz, dw);
    x = static_cast<float>(dx);
    y = static_cast<float>(dy);
    z = static_cast<float>(dz);
    w = static_cast<float>(dw);
  }

  void scan_normal(int& x, int& y, int& z, int& w)
  {
    double dx(0), dy(0), dz(0), dw(0);
    scan_normal(dx, dy, dz, dw);
    x = static_cast<int>(dx);
    y = static_cast<int>(dy);
    z = static_cast<int>(dz);
    w = static_cast<int>(dw);
  }

  void scan_normal(double& x, double& y, double& z)
  {
    double dx(0), dy(0), dz(0), dw(0);
    scan_normal(dx, dy, dz, dw);
    x = dx / dw;
    y = dy / dw;
    z = dz / dw;
  }

  void scan_normal(float& x, float& y, float& z)
  {
    double dx(0), dy(0), dz(0);
    scan_normal(dx, dy, dz);
    x = static_cast<float>(dx);
    y = static_cast<float>(dy);
    z = static_cast<float>(dz);
  }

  void scan_normal(int& x, int& y, int& z)
  {
    double dx(0), dy(0), dz(0);
    scan_normal(dx, dy, dz);
    x = static_cast<int>(dx);
    y = static_cast<int>(dy);
    z = static_cast<int>(dz);
  }

  static const IO::Color& get_indexed_color(int id)
  {
    static const IO::Color color[149] = {
      IO::Color(255, 255, 255, 191),
      IO::Color(255, 255, 255, 191),
      IO::Color(255, 255, 255, 191),
      IO::Color(255, 255, 255, 191),
      IO::Color(255, 255, 255, 191),
      IO::Color(255, 255, 255, 191),
      IO::Color(178, 38, 25, 191),
      IO::Color(51, 51, 204, 191),
      IO::Color(229, 153, 5, 191),
      IO::Color(25, 76, 204, 191),
      IO::Color(25, 178, 51, 191),
      IO::Color(204, 204, 102, 191),
      IO::Color(178, 178, 0, 191),
      IO::Color(178, 0, 178, 191),
      IO::Color(0, 178, 178, 191),
      IO::Color(229, 0, 51, 191),
      IO::Color(51, 229, 0, 191),
      IO::Color(0, 51, 229, 191),
      IO::Color(191, 191, 191, 191),
      IO::Color(204, 102, 0, 191),
      IO::Color(204, 102, 0, 191),
      IO::Color(0, 102, 204, 191),
      IO::Color(0, 102, 204, 191),
      IO::Color(0, 204, 102, 191),
      IO::Color(0, 204, 102, 191),
      IO::Color(102, 0, 204, 191),
      IO::Color(102, 0, 204, 191),
      IO::Color(204, 0, 102, 191),
      IO::Color(204, 0, 102, 191),
      IO::Color(178, 127, 51, 191),
      IO::Color(178, 127, 51, 191),
      IO::Color(178, 178, 0, 191),
      IO::Color(178, 0, 178, 191),
      IO::Color(0, 178, 178, 191),
      IO::Color(229, 0, 0, 191),
      IO::Color(0, 229, 0, 191),
      IO::Color(0, 0, 229, 191),
      IO::Color(191, 191, 191, 191),
      IO::Color(204, 102, 0, 191),
      IO::Color(102, 204, 0, 191),
      IO::Color(0, 102, 204, 191),
      IO::Color(0, 204, 102, 191),
      IO::Color(102, 0, 204, 191),
      IO::Color(204, 0, 102, 191),
      IO::Color(178, 178, 0, 191),
      IO::Color(178, 0, 178, 191),
      IO::Color(0, 178, 178, 191),
      IO::Color(229, 0, 0, 191),
      IO::Color(0, 229, 0, 191),
      IO::Color(0, 0, 229, 191),
      IO::Color(191, 191, 191, 191),
      IO::Color(204, 102, 0, 191),
      IO::Color(102, 204, 0, 191),
      IO::Color(0, 102, 204, 191),
      IO::Color(0, 204, 102, 191),
      IO::Color(102, 0, 204, 191),
      IO::Color(204, 0, 102, 191),
      IO::Color(178, 178, 0, 191),
      IO::Color(178, 0, 178, 191),
      IO::Color(0, 178, 178, 191),
      IO::Color(229, 0, 0, 191),
      IO::Color(0, 229, 0, 191),
      IO::Color(0, 0, 229, 191),
      IO::Color(191, 191, 191, 191),
      IO::Color(204, 102, 0, 191),
      IO::Color(102, 204, 0, 191),
      IO::Color(0, 102, 204, 191),
      IO::Color(0, 204, 102, 191),
      IO::Color(102, 0, 204, 191),
      IO::Color(204, 0, 102, 191),
      IO::Color(255, 255, 255, 191),
      IO::Color(255, 255, 255, 191),
      IO::Color(255, 255, 255, 191),
      IO::Color(255, 255, 255, 191),
      IO::Color(255, 255, 255, 191),
      IO::Color(255, 255, 255, 191),
      IO::Color(12, 76, 25, 191),
      IO::Color(178, 2, 25, 191),
      IO::Color(51, 12, 153, 191),
      IO::Color(229, 229, 5, 191),
      IO::Color(0, 51, 102, 191),
      IO::Color(25, 102, 102, 191),
      IO::Color(204, 204, 204, 191),
      IO::Color(178, 178, 0, 191),
      IO::Color(178, 178, 0, 191),
      IO::Color(178, 0, 178, 191),
      IO::Color(178, 0, 178, 191),
      IO::Color(0, 178, 178, 191),
      IO::Color(0, 178, 178, 191),
      IO::Color(229, 0, 0, 191),
      IO::Color(229, 0, 0, 191),
      IO::Color(0, 229, 0, 191),
      IO::Color(0, 229, 0, 191),
      IO::Color(0, 0, 229, 191),
      IO::Color(0, 0, 229, 191),
      IO::Color(191, 191, 191, 191),
      IO::Color(191, 191, 191, 191),
      IO::Color(204, 102, 0, 191),
      IO::Color(204, 102, 0, 191),
      IO::Color(0, 102, 204, 191),
      IO::Color(0, 102, 204, 191),
      IO::Color(0, 204, 102, 191),
      IO::Color(0, 204, 102, 191),
      IO::Color(102, 0, 204, 191),
      IO::Color(102, 0, 204, 191),
      IO::Color(204, 0, 102, 191),
      IO::Color(204, 0, 102, 191),
      IO::Color(178, 127, 51, 191),
      IO::Color(178, 127, 51, 191),
      IO::Color(178, 178, 0, 191),
      IO::Color(178, 0, 178, 191),
      IO::Color(0, 178, 178, 191),
      IO::Color(229, 0, 0, 191),
      IO::Color(0, 229, 0, 191),
      IO::Color(0, 0, 229, 191),
      IO::Color(191, 191, 191, 191),
      IO::Color(204, 102, 0, 191),
      IO::Color(102, 204, 0, 191),
      IO::Color(0, 102, 204, 191),
      IO::Color(0, 204, 102, 191),
      IO::Color(102, 0, 204, 191),
      IO::Color(204, 0, 102, 191),
      IO::Color(178, 178, 0, 191),
      IO::Color(178, 0, 178, 191),
      IO::Color(0, 178, 178, 191),
      IO::Color(229, 0, 0, 191),
      IO::Color(0, 229, 0, 191),
      IO::Color(0, 0, 229, 191),
      IO::Color(191, 191, 191, 191),
      IO::Color(204, 102, 0, 191),
      IO::Color(102, 204, 0, 191),
      IO::Color(0, 102, 204, 191),
      IO::Color(0, 204, 102, 191),
      IO::Color(102, 0, 204, 191),
      IO::Color(204, 0, 102, 191),
      IO::Color(178, 178, 0, 191),
      IO::Color(178, 0, 178, 191),
      IO::Color(0, 178, 178, 191),
      IO::Color(229, 0, 0, 191),
      IO::Color(0, 229, 0, 191),
      IO::Color(0, 0, 229, 191),
      IO::Color(191, 191, 191, 191),
      IO::Color(204, 102, 0, 191),
      IO::Color(102, 204, 0, 191),
      IO::Color(0, 102, 204, 191),
      IO::Color(0, 204, 102, 191),
      IO::Color(102, 0, 204, 191),
      IO::Color(204, 0, 102, 191),
      IO::Color(120, 120, 120, 120) };
    if(id > 148) id =148;
    return color[id];
  }

  static CGAL::IO::Color get_color_from_line(std::istream &is)
  {
    std::string color_info;
    bool is_float = false;

    std::string col;
    //get the line content
    std::streampos position = is.tellg();
    std::getline(is, col);
    //split it into strings
    std::istringstream iss(col);
    //holds the rgb values
    unsigned char rgb[3]{};
    int index =0;
    //split the string into numbers
    while(iss>>color_info){
      //stop if comment is read
      if(color_info.at(0) == '#')
        break;
      //detect if the value is float
      for(int c = 0; c<static_cast<int>(color_info.length()); ++c)
      {
        if(color_info.at(c) == '.')
        {
          is_float = true;
       //   break;
        }
        if(color_info.at(c) == '#')
        {
          color_info.resize(c);
          break;
        }
      }

      //if the value is of float type, convert it into an int
      if(is_float)
        rgb[index] = static_cast<unsigned char>(atof(color_info.c_str())*255);

      //else stores the value
      else
        rgb[index] = static_cast<unsigned char>(atoi(color_info.c_str()));
      ++index;
      if(index == 3)
        break;
    }
    CGAL::IO::Color color;
    //if there were only one number, fetch the color in the color map
    if(index < 2)
    {
      color = get_indexed_color(rgb[0]);
      //else create the color with the 3 values;
    }
    else{
      color = CGAL::IO::Color(rgb[0], rgb[1], rgb[2]);
    }
    std::iostream::pos_type ss_pos = iss.tellg();
    if(ss_pos != std::iostream::pos_type(-1))
    {
      position +=ss_pos;
      is.seekg(position);
    }
    return color;
  }

  void scan_color(unsigned char& r, unsigned char& g, unsigned char& b)
  {
    if(binary())
    {
      float fr, fg, fb;
      I_Binary_read_big_endian_float32(m_in, fr);
      I_Binary_read_big_endian_float32(m_in, fg);
      I_Binary_read_big_endian_float32(m_in, fb);
      r = (unsigned char)(fr);
      g = (unsigned char)(fg);
      b = (unsigned char)(fb);

    }
    else
    {
      CGAL::IO::Color color;
      if(color_entries == 1){
        color = get_indexed_color(static_cast<int>(entries[first_color_index])); // the index in the color map
        r = color.red();
        g = color.green();
        b = color.blue();
        return;
      }
      double rd = entries[first_color_index];
      double gd = entries[first_color_index + 1];
      double bd = entries[first_color_index + 2];

      if( (floor(rd) == rd) &&  (floor(gd) == gd) && (floor(bd) == bd)){
        // we have to do with integers
        r = static_cast<unsigned char>(rd);
        g = static_cast<unsigned char>(gd);
        b = static_cast<unsigned char>(bd);
      }else{
        // we have to do with floats
        r = static_cast<unsigned char>(rd*255);
        g = static_cast<unsigned char>(gd*255);
        b = static_cast<unsigned char>(bd*255);
      }
      if(color_entries == 4){
        //double alphad = entries[first_color_index + 3];
        // it seems that we ignore it.
      }
    }
  }

  void skip_to_next_vertex(std::size_t current_vertex)
  {
    CGAL_assertion(current_vertex < size_of_vertices());
    if(binary())
    {
      float f;
      if(has_normals() && ! normals_read) {
        I_Binary_read_big_endian_float32(m_in, f);
        I_Binary_read_big_endian_float32(m_in, f);
        I_Binary_read_big_endian_float32(m_in, f);
        if(is_homogeneous())
          I_Binary_read_big_endian_float32(m_in, f);
      }

      if(has_colors())
      {
        std::int32_t k;
        I_Binary_read_big_endian_integer32(m_in, k);
        if(k<0 || k>4)
        {
          m_in.clear(std::ios::badbit);
          if(verbose())
          {
            std::cerr << " " << std::endl;
            std::cerr << "File_scanner_OFF::" << std::endl;
            std::cerr << "skip_to_next_vertex(): input error: bad "
                         " number of color indices at vertex "
                      << current_vertex << "." << std::endl;
          }

          set_off_header(false);
          return;
        }

        while(k--)
        {
          float dummy;
          I_Binary_read_big_endian_float32(m_in, dummy);
        }
      }
    }
  }

  void scan_facet(std::size_t& size, std::size_t CGAL_assertion_code(current_facet))
  {
    CGAL_assertion(current_facet < size_of_facets());
    if(binary())
    {
      std::int32_t i32;
      I_Binary_read_big_endian_integer32(m_in, i32);
      size = i32;
    }
    else
    {
      skip_comment();
      line.clear();
      std::getline(m_in, line);
      // First remove the comment if there is one
      std::size_t pos = line.find('#');
      if(pos != std::string::npos){
        line = line.substr(0,pos);
      }

      // Read all numbers in the line
      std::istringstream issline(line);
      entries.clear();
      double d;
      while(issline >> IO::iformat(d)){
        entries.push_back(d);
      }
      if(entries.empty())
      {
        m_in.clear(std::ios::badbit);
        size = 0;
        return;
      }
      size = static_cast<std::size_t>(entries[0]);
      if(has_colors()){
        // Compute how many entries are there for the color
        first_color_index = size + 1;

        color_entries = entries.size();
        color_entries -= size  + 1; // coordinates
        // now color_entries should be 0, 1, 3, or 4 for the color
        if(color_entries > 0)
          m_has_fcolors = true;
      }
    }
  }

  void scan_facet_vertex_index(std::size_t& index,
                               const std::size_t& current_entry,
                               std::size_t current_facet)
  {
    if(binary()){
      std::int32_t i32;
      I_Binary_read_big_endian_integer32(m_in, i32);
      index = i32;
    }
    else
    {
      if(entries.size() <= current_entry )
      {
        m_in.clear(std::ios::badbit);
        if(verbose())
          std::cerr<<"error while reading facet. Missing index."<<std::endl;
        index=0;
        return;
      }
      index = static_cast<std::size_t>(entries[current_entry]);
    }

    if(m_in.fail())
    {
      if(verbose())
      {
        std::cerr << " " << std::endl;
        std::cerr << "File_scanner_OFF::" << std::endl;
        std::cerr << "scan_facet_vertex_index(): input error:  "
                     "cannot read OFF file beyond facet "
                  << current_facet << "." << std::endl;
      }
      index=0;
      set_off_header(false);
      return;
    }

    bool error  = index < index_offset();
    index -= index_offset();

    if(error || (index >= size_of_vertices()))
    {
      m_in.clear(std::ios::failbit);
      if(verbose())
      {
        std::cerr << " " << std::endl;
        std::cerr << "File_scanner_OFF::" << std::endl;
        std::cerr << "scan_facet_vertex_index(): input error: "
                     "facet " << current_facet << ": vertex index "
                  << index + index_offset() << ": is out of range."
                  << std::endl;
      }
      index = 0;
      set_off_header(false);
      return;
    }
  }

  void skip_to_next_facet(std::size_t current_facet)
  {
    // Take care of trailing information like color triples.
    if(binary())
    {
      std::int32_t k;
      I_Binary_read_big_endian_integer32(m_in, k);
      if(k<0 || k>4)
      {
        m_in.clear(std::ios::badbit);
        if(verbose())
        {
          std::cerr << " " << std::endl;
          std::cerr << "File_scanner_OFF::" << std::endl;
          std::cerr << "skip_to_next_facet(): input error: bad "
                       "number of color indices at vertex "
                    << current_facet << "." << std::endl;
        }

        set_off_header(false);
        return;
      }

      while (k--)
      {
        float dummy;
        I_Binary_read_big_endian_float32(m_in, dummy);
      }
    }
  }
};

template < class Point> inline
Point& file_scan_vertex(File_scanner_OFF& scanner, Point& p)
{
  typedef typename Point::R R;
  typedef typename R::RT    RT;
  double x(0), y(0), z(0), w(0);

  scanner.scan_vertex(x, y, z, w);

  if(w == 1)
    p = Point(RT(x), RT(y), RT(z));
  else
    p = Point(RT(x), RT(y), RT(z), RT(w));

  return p;
}

template < class T_Color> inline
T_Color& file_scan_color(File_scanner_OFF& scanner, T_Color& c)
{
  unsigned char r, g, b;
  scanner.scan_color(r,g,b);
  c = T_Color(r,g,b);

  return c;
}

template < class Vector> inline
Vector&
file_scan_normal(File_scanner_OFF& scanner, Vector& v)
{
  typedef typename Vector::R R;
  typedef typename R::RT     RT;

  double x(0), y(0), z(0), w(0);
  scanner.scan_normal(x, y, z, w);
  if(w == 1)
    v = Vector(RT(x), RT(y), RT(z));
  else
    v = Vector(RT(x), RT(y), RT(z), RT(w));

  return v;
}

} //namespace CGAL

#endif // CGAL_IO_OFF_FILE_SCANNER_OFF_H
