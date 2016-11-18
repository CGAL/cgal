// Copyright (c) 2016 GeometryFactory
// All rights reserved.
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
// Author(s)     : Laurent Rineau

#ifndef CGAL_READ_SEP_IMAGE_DATA_H
#define CGAL_READ_SEP_IMAGE_DATA_H

#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>

#include <limits>

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>

#ifndef BOOST_FILESYSTEM_VERSION
// That macro was not defined in previous releases of Boost.
#  define BOOST_FILESYSTEM_VERSION 2
#endif

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>

#include <CGAL/ImageIO.h>
#include <CGAL/Image_3.h>

namespace read_file {

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;

template <typename Iterator>
bool parse_comment_or_empty(Iterator first, Iterator last) {
  using qi::phrase_parse;
  using ascii::space;
  using qi::char_;
  using qi::string;
  using qi::eoi;
  
  bool r = phrase_parse(
                        first,                           /*< start iterator >*/
                        last,                            /*< end iterator >*/
                        string("#") >> *(char_)|eoi,     /*< the parser >*/
                        space                            /*< the skip-parser >*/
                        );
  if (first != last) // fail if we did not get a full match
    return false;
  return r;
}

template <typename Iterator, typename Parser, typename Value>
bool parse_assignement(Iterator first, Iterator last, 
                       Parser parser, 
                       std::string& parameter_name, Value& value)
{
  using qi::phrase_parse;
  using ascii::space;
  using qi::char_;
  using boost::spirit::qi::_1;
  using boost::phoenix::ref;
  using boost::phoenix::push_back;

  bool r = phrase_parse(first,                               /*< start iterator >*/
                        last,                                /*< end iterator >*/
                        *(char_ - '=')[push_back( boost::phoenix::ref(parameter_name), _1)] 
                        >> '=' 
                        >> parser[ boost::phoenix::ref(value) = _1],         /*< the parser >*/
                        space                                /*< the skip-parser >*/
                        );
  if (first != last) // fail if we did not get a full match
    return false;
  return r;
}

template <typename Iterator>
bool parse_string_assignement(Iterator first, Iterator last, 
                              std::string& parameter_name,
                              std::string& value)
{
  using qi::phrase_parse;
  using ascii::space;
  using qi::char_;
  using boost::spirit::qi::_1;
  using boost::phoenix::ref;
  using boost::phoenix::push_back;

  bool r = phrase_parse(first,                               /*< start iterator >*/
                        last,                                /*< end iterator >*/
                        *(char_ - '=')[push_back( boost::phoenix::ref(parameter_name), _1)] 
                        >> '=' 
                        >> *char_[push_back( boost::phoenix::ref(value), _1)],
                                                             /*< the parser >*/
                        space                                /*< the skip-parser >*/
                        );
  if (first != last) // fail if we did not get a full match
    return false;
  return r;
}

} //  end namespace read_file

template <typename T>
class Sep_reader 
{
  
protected:
  
  std::string _fileName;

  boost::array<std::size_t, 3> _n;
  boost::array<double, 3> _d;
  boost::array<double, 3> _o;
  
  T* _data;
  _image* _im;
  CGAL::Image_3* _cgal_image;

  bool valid;
  std::string err_msg;

  std::map<std::string, int> int_dict;
  std::map<std::string, double> double_dict;
  std::map<std::string, std::string> string_dict;

  int dimension;

  const T special_value;

  bool parseHeader(std::string fileName) {
    using boost::spirit::qi::int_;
    using boost::spirit::qi::double_;

    std::ifstream input(fileName.c_str());
    if(!input) {
      std::cerr << "Error: cannot open the header file \"" 
                << fileName << "\"!\n";
      return false;
    }
    std::string line;
    while(std::getline(input, line)) {
      if(read_file::parse_comment_or_empty(line.begin(), line.end())) continue;
      std::string parameter;
      int i;
      if(read_file::parse_assignement(line.begin(), line.end(),
                                      int_,
                                      parameter,
                                      i))
      {
        int_dict[parameter] = i;
        continue;
      }
      parameter.clear();
      double d;
      if(read_file::parse_assignement(line.begin(), line.end(),
                                      double_,
                                      parameter,
                                      d))
      {
        // std::cerr << parameter << "=" << d << std::endl;
        double_dict[parameter] = d;
        continue;
      }
      parameter.clear();
      std::string s;
      if(read_file::parse_string_assignement(line.begin(), line.end(),
                                             parameter,
                                             s))
      {
        // std::cerr << parameter << "=" << s << std::endl;
        string_dict[parameter] = s;
        continue;
      }
      return false;
    }
    return true;
  }

  double& o(int i) { return _o[i-1]; }
  double& d(int i) { return _d[i-1]; }
  std::size_t& n(int i) { return _n[i-1]; }
public:
  const double& o(int i) const { return _o[i-1]; }
  const double& d(int i) const { return _d[i-1]; }
  const std::size_t& n(int i) const { return _n[i-1]; }

  bool is_valid() const { return valid; }
  std::string error_msg() const { return err_msg; }

  Sep_reader(std::string fileName, 
             const T special_value = std::numeric_limits<T>::infinity())
    : _fileName(fileName), _im(0), _cgal_image(0),
      valid(false), err_msg()
    , dimension(0), special_value(special_value)
  {
    if(!parseHeader(fileName)) {
      err_msg = "Invalid header file \"";
      err_msg += fileName;
      err_msg += "\"";
      return;
    }
    n(1) = int_dict["n1"];
    n(2) = int_dict["n2"];
    n(3) = int_dict["n3"];
    if(n(3)==0) {
      n(3) = 1;
      dimension = 2;
    }
    else {
      dimension = 3;
    }
    d(1) = double_dict["d1"];
    d(2) = double_dict["d2"];
    d(3) = double_dict["d3"];
    o(1) = double_dict["o1"];
    o(2) = double_dict["o2"];
    o(3) = double_dict["o3"];

    std::cout << "Header file: " << fileName << std::endl;
    std::cout << "Parameters:\n";
    std::cout << "dimension = " << dimension << "\n";
    for(int i = 1; i <= 3 ; ++i) {
      std::cout << "d" << i << " = " << d(i) << "\n";
      std::cout << "n" << i << " = " << n(i) << "\n";
      std::cout << "o" << i << " = " << o(i) << "\n";
    }
    std::cout << "in = " << string_dict["in"] << std::endl;

    boost::filesystem::path headerFile(fileName);
    boost::filesystem::path dataFile(string_dict["in"]);
#if BOOST_FILESYSTEM_VERSION == 2
    dataFile = boost::filesystem::complete(dataFile, 
                                           boost::filesystem::complete(headerFile.parent_path()));
#else
    dataFile = boost::filesystem::absolute(dataFile, 
                                           boost::filesystem::absolute(headerFile.parent_path()));
#endif
    if(!load_data(dataFile.string())) {
      return;
      err_msg = "Invalid data file \"";
      err_msg += dataFile.string();
      err_msg += "\"";
    }
    std::cout << "min = " << (min)() << std::endl;
    std::cout << "special value  = " << special_value << std::endl;
    valid = true;
  }

  ~Sep_reader()
  {
    if(_cgal_image != 0) delete _cgal_image;
    std::cerr << "Sep_reader destroyed\n";
  }

  CGAL::Image_3* cgal_image() { return _cgal_image; }
  const CGAL::Image_3* cgal_image() const { return _cgal_image; }

  const T& raw_value(std::size_t i, std::size_t j, std::size_t k = 0) const {
    const T* result = &_data[i + n(1) * ( j + n(2)*k )];
    // if(*result <= special_value) result = &raw_value(i, j, k + 1);
    return *result;
  }
  
  T interpolation(const double& x, const double& y) const
  {
    const double lx = (x - o(1)) / d(1);
    const double ly = (y - o(2)) / d(2);

    const std::size_t i1 = convert2index<1>(x);
    const std::size_t j1 = convert2index<2>(y);
    std::size_t i2 = i1 + 1;
    std::size_t j2 = j1 + 1;

    if(i2 == n(1)) i2 = n(1) - 1;
    if(j2 == n(2)) j2 = n(2) - 1;

    const double di2 = i2 - lx;
    const double di1 = lx - i1;
    const double dj2 = j2 - ly;
    const double dj1 = ly - j1;

    const double v11 = raw_value(i1, j1);
    const double v12 = raw_value(i1, j2);

    return 
      (v11 * di2 + v12 * di1) * dj2 + 
      (v11 * di2 + v12 * di1) * dj1;
  }

  const T& operator()(const double& x, const double& y) const
  {
    const std::size_t i = convert2index<1>(x);
    const std::size_t j = convert2index<2>(y);
    return raw_value(i, j, 0);
  }


  const T& operator()(const double& x, const double& y, const double& z) const
  {
    const std::size_t i = convert2index<1>(x);
    const std::size_t j = convert2index<2>(y);
    const std::size_t k = convert2index<3>(z);
    return raw_value(i, j, k);
  }

  template <int dim>
  std::size_t convert2index(const double& x) const {
    const double k = (x - o(dim)) / d(dim);
    if (k < 0) return 0;
    else if (k >= n(dim)) return n(dim) - 1;
    else return static_cast<std::size_t>(k);
  }

  template <int dim>
  double index2coord(const std::size_t& i) const {
    return o(dim) + d(dim) * i;
  }


  T min BOOST_PREVENT_MACRO_SUBSTITUTION() const {
    return *std::min_element(_data, _data + n(1) * n(2) * n(3));
  }
  T max BOOST_PREVENT_MACRO_SUBSTITUTION() const {
    return *std::max_element(_data, _data + n(1) * n(2) * n(3));
  }

protected :

  bool load_data(std::string dataFilename) 
  {
    if(_im) delete _im;
    _im = new _image;
    _im->xdim = static_cast<unsigned int>(n(1));
    _im->ydim = static_cast<unsigned int>(n(2));
    _im->zdim = static_cast<unsigned int>(n(3));
    _im->vdim = 1;

    // spacing
    _im->vx = d(1);
    _im->vy = d(2);
    _im->vz = d(3);
    if(_im->vz == 0.f) _im->vz=1.f;

    // image offset
    _im->tx = static_cast<float>(o(1));
    _im->ty = static_cast<float>(o(2));
    _im->tz = static_cast<float>(o(3));

    // image center
    _im->cx = _im->cy = _im->cz = 0;
    // image rotation
    _im->rx = _im->ry = _im->rz = 0.0;

    _im->fd = NULL;
    _im->openMode = OM_CLOSE;
    if(string_dict["data_format"] == "native_float" ||
       string_dict["data_format"] == "\"native_float\"")
    {
      std::cerr << "little endian\n";
      _im->endianness =        END_LITTLE;
    } else {
      std::cerr << "big endian\n";
      _im->endianness =        END_BIG;
    }

    _im->dataMode = DM_BINARY;

    // no user string
    _im->user = NULL;
    _im->nuser = 0;

    // word type (unsigned byte)
    _im->wdim = sizeof(T);
    _im->wordKind = WK_FLOAT;
    _im->vectMode = VM_SCALAR;
    _im->sign = SGN_SIGNED;
    _im->imageFormat = NULL;

    ::_openReadImage(_im, dataFilename.c_str());
    if(!_im->fd) return false;

    // Compute number of element
    const std::size_t size = n(1) * n(2) * n(3);
    
    // Allocate array
    _data = new T[size];
    _im->data = (void*)_data;
    
    // // Read file
    if(_readImageData(_im) < 0) return false;
   // char* buffer = reinterpret_cast<char*>(_data);
    // in.read (buffer, size * sizeof(T));

    ImageIO_close(_im);
    _cgal_image = new CGAL::Image_3(_im);
    return true;
  }
};

#endif // CGAL_READ_SEP_IMAGE_DATA_H
