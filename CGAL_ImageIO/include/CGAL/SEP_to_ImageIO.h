// Copyright (c) 2016, 2017 GeometryFactory
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
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_SEP_READER_IMAGEIO_HPP
#define CGAL_SEP_READER_IMAGEIO_HPP

#include "SEP_header.h"
#include <CGAL/ImageIO.h>
#include <CGAL/Image_3.h>

#include <fstream>
#include <algorithm>
#include <string>

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>

#ifndef BOOST_FILESYSTEM_VERSION
// That macro was not defined in previous releases of Boost.
#  define BOOST_FILESYSTEM_VERSION 2
#endif

#include <CGAL/IO/binary_file_io.h>

namespace CGAL {

template <typename T>
class SEP_to_ImageIO : public SEP_header
{
  T* _data;
  _image* _im;
  CGAL::Image_3* _cgal_image;

  bool valid;
  std::string err_msg;

public:
  bool is_valid() const { return valid; }
  std::string error_msg() const { return err_msg; }

  SEP_to_ImageIO(std::string fileName)
    : SEP_header(fileName),
      _data(0),
      _im(0), _cgal_image(0),
      valid(false), err_msg()
  {
    if(dimension() < 0) {
      err_msg = "Invalid header file \"";
      err_msg += fileName;
      err_msg += "\"";
      return;
    }
    display_information(fileName, std::cout);

    boost::filesystem::path headerFile(fileName);
    boost::filesystem::path dataFile(string_field("in"));
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
    valid = true;
  }

  ~SEP_to_ImageIO()
  {
    if(_cgal_image != 0) delete _cgal_image;
  }

  CGAL::Image_3* cgal_image() { return _cgal_image; }
  const CGAL::Image_3* cgal_image() const { return _cgal_image; }

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
    if(string_field("data_format") == "native_float" ||
       string_field("data_format") == "\"native_float\"")
    {
      _im->endianness =        END_LITTLE;
    } else {
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

} // end namespace CGAL

#endif // CGAL_SEP_READER_IMAGEIO_HPP
