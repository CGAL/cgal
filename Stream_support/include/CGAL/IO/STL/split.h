// Copyright (c) 2015 GeometryFactory
//
// This file is part of CGAL (www.cgal.org);
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Andreas Fabri,
//                 Mael Rouxel-Labbé

#ifndef CGAL_IO_STL_SPLIT_H
#define CGAL_IO_STL_SPLIT_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <boost/cstdint.hpp>

#include <CGAL/Bbox_3.h>


namespace CGAL {
namespace IO {


inline
bool split_binary_STL(std::istream& is,
                      std::string stem,
                      Bbox_3& bb,
                      int ncells,
                      const bool verbose = false)
{
  double xmin = bb.xmin()-0.1;
  double xmax = bb.xmax()+0.1;
  double ymin = bb.ymin()-0.1;
  double ymax = bb.ymax()+0.1;
  double zmin = bb.zmin()-0.1;
  double zmax = bb.zmax()+0.1;

  double cwx = (xmax - xmin) / ncells;
  double cwy = (ymax - ymin) / ncells;
  double cwz = (zmax - zmin) / ncells;

  std::vector<std::ofstream> streams;
  streams.reserve(ncells*ncells*ncells);

  std::vector<int> ntriangles(ncells*ncells*ncells,0);
  for(int i = 0; i < ncells; i++)
    for(int j = 0; j < ncells; j++)
      for(int k = 0; k < ncells; k++){
        int index = i*ncells*ncells + j*ncells +k;
        std::string fileName = stem + "_" + std:: to_string(i) + "_"  + std:: to_string(j) + "_" + std:: to_string(k) +   ".stl";
        streams.emplace_back(std::ofstream {fileName.c_str(), std::ios::binary} );
        streams[index] << "FileType: Binary                                                                ";
        boost::uint32_t N32 = 0; // has to be rewritten later
        streams[index].write(reinterpret_cast<const char *>(&N32), sizeof(N32));
      }


  // Start from the beginning again to simplify things
  is.clear();
  is.seekg(0, std::ios::beg);

  if(!is.good())
    return false;

  // Discard the first 80 chars (unused header)
  int pos = 0;
  char c;

  if(verbose)
    std::cout << "header: ";

  while(pos < 80)
  {
    is.read(reinterpret_cast<char*>(&c), sizeof(c));
    if(!is.good())
      break;

    if(verbose)
      std::cout << c;

    ++pos;
  }

  if(verbose)
    std::cout << std::endl;

  if(pos != 80)
    return true; // empty file

  int index = 0;

  boost::uint32_t N32;
  if(!(is.read(reinterpret_cast<char*>(&N32), sizeof(N32))))
  {
    if(verbose)
      std::cerr << "Error while reading number of facets" << std::endl;

    return false;
  }

  unsigned int N = N32;
  if(verbose)
    std::cout << N << " facets to read" << std::endl;

  for(unsigned int ti=0; ti<N; ++ti)
  {
    float normal[3];
    if(!(is.read(reinterpret_cast<char*>(&normal[0]), sizeof(normal[0]))) ||
       !(is.read(reinterpret_cast<char*>(&normal[1]), sizeof(normal[1]))) ||
       !(is.read(reinterpret_cast<char*>(&normal[2]), sizeof(normal[2]))))
    {
      if(verbose)
        std::cerr << "Error while reading normal coordinates (premature end of file)" << std::endl;

      return false;
    }

    int xrange[2], yrange[2], zrange[2];

    float f=0, x[3], y[3], z[3];
    for(int j=0; j<3; ++j)
    {
      if(!(is.read(reinterpret_cast<char*>(&x[j]), sizeof(f))) ||
         !(is.read(reinterpret_cast<char*>(&y[j]), sizeof(f))) ||
         !(is.read(reinterpret_cast<char*>(&z[j]), sizeof(f))))
      {
        if(verbose)
          std::cerr << "Error while reading vertex coordinates (premature end of file)" << std::endl;

        return false;
      }

      int xr = std::floor( (x[j]-xmin)/cwx );
      int yr = std::floor( (y[j]-ymin)/cwy );
      int zr = std::floor( (z[j]-zmin)/cwz );
      if(j ==0){
        xrange[0] = xr; xrange[1] = xr;
        yrange[0] = yr; yrange[1] = yr;
        zrange[0] = zr; zrange[1] = zr;
      }else{
        xrange[0] = (std::min)(xrange[0],xr); xrange[1] = (std::max)(xrange[1],xr);
        yrange[0] = (std::min)(yrange[0],yr); yrange[1] = (std::max)(yrange[1],yr);
        zrange[0] = (std::min)(zrange[0],zr); zrange[1] = (std::max)(zrange[1],zr);
      }
    }

    for(int i = xrange[0]; i <= xrange[1]; ++i)
      for(int j = yrange[0]; j <= yrange[1]; ++j)
        for(int k = zrange[0]; k <= zrange[1]; ++k){
          index = i * ncells * ncells + j * ncells + k;
          ++ntriangles[index];
          streams[index].write(reinterpret_cast<const char *>(&normal[0]), sizeof(f));
          streams[index].write(reinterpret_cast<const char *>(&normal[1]), sizeof(f));
          streams[index].write(reinterpret_cast<const char *>(&normal[2]), sizeof(f));
          streams[index].write(reinterpret_cast<const char *>(&x[0]), sizeof(f));
          streams[index].write(reinterpret_cast<const char *>(&y[0]), sizeof(f));
          streams[index].write(reinterpret_cast<const char *>(&z[0]), sizeof(f));
          streams[index].write(reinterpret_cast<const char *>(&x[1]), sizeof(f));
          streams[index].write(reinterpret_cast<const char *>(&y[1]), sizeof(f));
          streams[index].write(reinterpret_cast<const char *>(&z[1]), sizeof(f));
          streams[index].write(reinterpret_cast<const char *>(&x[2]), sizeof(f));
          streams[index].write(reinterpret_cast<const char *>(&y[2]), sizeof(f));
          streams[index].write(reinterpret_cast<const char *>(&z[2]), sizeof(f));
          streams[index] << "  ";
        }

    // Read so-called attribute byte count and ignore it
    char c;
    if(!(is.read(reinterpret_cast<char*>(&c), sizeof(c))) ||
       !(is.read(reinterpret_cast<char*>(&c), sizeof(c))))
    {
      if(verbose)
        std::cerr << "Error while reading attribute byte count (premature end of file)" << std::endl;

      return false;
    }
  }

  for(int i = 0; i < streams.size(); ++i){
     streams[i].seekp(80);
     N32 =  static_cast<boost::uint32_t>(ntriangles[i]);
     streams[i].write(reinterpret_cast<const char *>(&N32), sizeof(N32));
     std::cout << "# triangles : " << N32 << std::endl;
     streams[i].close();
  }

  return !is.fail();
}

}
}


#endif // CGAL_IO_STL_SPLIT_H
