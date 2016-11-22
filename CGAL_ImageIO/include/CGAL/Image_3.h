// Copyright (c) 2005-2008  INRIA Sophia-Antipolis (France).
//               2008 GeometryFactory
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
// Author(s)     : Laurent Rineau, Pierre Alliez

#ifndef CGAL_IMAGE_3_H
#define CGAL_IMAGE_3_H

#include <CGAL/basic.h>
#include <CGAL/array.h>

#include <boost/shared_ptr.hpp>
#include <boost/format.hpp>
#include <boost/unordered_set.hpp>
#include <CGAL/ImageIO.h>
#include <CGAL/function_objects.h>

#include <limits>
#include <set>
#include <cstdlib>

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4244 4251) // double float conversion loss of data and dll linkage
#endif

class vtkImageData;

namespace CGAL {

namespace ImageIO {

template <typename T>
struct Indicator_factory
{
  class Indicator : public std::unary_function<T, double>
  {
    const T label;
  public:
    Indicator(T i) : label(i) {}

    double operator()(T x) const
    {
      return (x == label) ? 1. : 0.;
    }
  }; // end nested class Indicator

  Indicator indicator(T i) const {
    return Indicator(i);
  }
};

} // end namespace CGAL::ImageIO

class CGAL_IMAGEIO_EXPORT Image_3
{

  struct Image_deleter {
    void operator()(_image* image)
    {
      ::_freeImage(image);
    }
  };
public:
  typedef boost::shared_ptr<_image> Image_shared_ptr;
  typedef Image_shared_ptr Pointer;

protected:
  Image_shared_ptr image_ptr;

   // implementation in src/CGAL_ImageIO/Image_3.cpp
  bool private_read(_image* im);

public:
  Image_3()
    : image_ptr()
  {
  }

  Image_3(const Image_3& bi)
    : image_ptr(bi.image_ptr)
  {
//     std::cerr << "Image_3::copy_constructor\n";
  }

  Image_3(_image* im) 
  {
    private_read(im);
  }

  ~Image_3()
  {
  }

  const _image* image() const
  {
    return image_ptr.get();
  }

  _image* image()
  {
    return image_ptr.get();
  }

  const void* data() const
  {
    return image()->data;
  }

  void* data()
  {
    return image()->data;
  }

  void set_data(void* d)
  {
    image()->data = d;
  }

  std::size_t xdim() const { return image_ptr->xdim; }
  std::size_t ydim() const { return image_ptr->ydim; }
  std::size_t zdim() const { return image_ptr->zdim; }

  std::size_t size() const { return xdim() * ydim() * zdim(); }

  double vx() const { return image_ptr->vx; }
  double vy() const { return image_ptr->vy; }
  double vz() const { return image_ptr->vz; }

  float value(const std::size_t i,
              const std::size_t j,
              const std::size_t k) const
  {
    return ::evaluate(image(),i,j,k);
  }

public:

  bool read(const char* file)
  {
    return private_read(::_readImage(file));
  }

  bool read_raw(const char* file,
                const unsigned int rx,
                const unsigned int ry,
                const unsigned int rz,
                const double vx = 1,
                const double vy = 1,
                const double vz = 1,
                const unsigned int offset = 0,
                std::size_t wdim=1,
                WORD_KIND wk = WK_FIXED,
                SIGN sign = SGN_UNSIGNED
                )
  {
    return private_read(::_readImage_raw(file,
                                         rx,ry,rz,
                                         vx,vy,vz,offset,
                                         wdim, wk, sign));
  }

public:
  template <typename Image_word_type,
	    typename Target_type,
	    typename Coord_type,
	    class Image_transform>
  Target_type
  trilinear_interpolation(const Coord_type&x, 
			  const Coord_type&y, 
			  const Coord_type&z,
			  const Target_type& value_outside =
			    Target_type(),
			  Image_transform transform = 
			    Image_transform() ) const;

  // default Image_transform = CGAL::Identity
  template <typename Image_word_type,
	    typename Target_type,
	    typename Coord_type>
  Target_type
  trilinear_interpolation(const Coord_type&x, 
			  const Coord_type&y, 
			  const Coord_type&z,
			  const Target_type& value_outside =
			  Target_type()) const
  {
    return trilinear_interpolation<
      Image_word_type,
      Target_type>(x, y, z, value_outside,
			  CGAL::Identity<Image_word_type>());
  }

  template <typename Image_word_type,
	    typename Coord_type,
            typename Target_type>
  Target_type
  labellized_trilinear_interpolation(const Coord_type&x,
				     const Coord_type&y,
				     const Coord_type&z,
				     const Target_type& value_outside =
  				       Target_type()) const
  {
    CGAL::ImageIO::Indicator_factory<Image_word_type> indicator_factory;
    return labellized_trilinear_interpolation<Image_word_type>
      (x, y, z, value_outside, indicator_factory);
  }

  template <typename Image_word_type,
	    typename Coord_type,
            typename Target_type,
            typename Indicator_factory>
  Target_type
  labellized_trilinear_interpolation
    (const Coord_type&x,
     const Coord_type&y,
     const Coord_type&z,
     const Target_type& value_outside,
     const Indicator_factory indicator_factory) const;

}; // end Image_3

template <typename Image_word_type,
	  typename Target_type,
	  typename Coord_type,
	  class Image_transform>
Target_type
Image_3::trilinear_interpolation(const Coord_type& x, 
				 const Coord_type& y, 
				 const Coord_type& z,
				 const Target_type& value_outside,
				 Image_transform transform) const 
{
  // Check on double/float coordinates, because (int)-0.1 gives 0
  if ( x < 0 || y < 0 || z < 0 )
    return value_outside;
  
  const Coord_type lx = x / image()->vx;
  const Coord_type ly = y / image()->vy;
  const Coord_type lz = z / image()->vz;
  const std::size_t dimx = xdim();
  const std::size_t dimy = ydim();
  const std::size_t dimz = zdim();
  const std::size_t dimxy = dimx*dimy;
  
  if(lx < 0 ||
     ly < 0 ||
     lz < 0 ||
     lz >= Coord_type(dimz-1) ||
     ly >= Coord_type(dimy-1) ||
     lx >= Coord_type(dimx-1))
  {
    return value_outside;
  }  

  // images are indexed by (z,y,x)
  const int i1 = (int)(lz); 
  const int j1 = (int)(ly);
  const int k1 = (int)(lx);
  const int i2 = i1 + 1;
  const int j2 = j1 + 1;
  const int k2 = k1 + 1;

  /*   We assume (x,y,z) lies in the following cube.
   *   a, b, c, d, e, f, g, h are the value of the image at the corresponding
   *   voxels:
   *
   *
   *     x        z
   *     |       /
   *        f___ __ g
   *       /|      /|
   *     e/_|____h/ |
   *     |  |    |  |
   *     |  |b___|_c|
   *     | /     | /
   *    a|/_____d|/  _y
   *    
   *
   * a = val(i1, j1, k1)
   * b = val(i2, j1, k1)
   * c = val(i2, j2, k1)
   * d = val(i1, j2, k1)
   * e = val(i1, j1, k2)
   * f = val(i2, j1, k2)
   * g = val(i2, j2, k2)
   * h = val(i1, j2, k2)
   */

  Image_word_type* ptr = (Image_word_type*)image()->data;
  ptr += i1 * dimxy + j1 * dimx + k1;
  const Target_type a = Target_type(transform(*ptr));
  const Target_type e = Target_type(transform(*(ptr+1)));
  ptr += dimxy; // i2 * dimxy + j1 * dimx + k1;
  const Target_type b = Target_type(transform(*ptr));
  const Target_type f = Target_type(transform(*(ptr+1)));
  ptr += dimx; // i2 * dimxy + j2 * dimx + k1
  const Target_type c = Target_type(transform(*ptr));
  const Target_type g = Target_type(transform(*(ptr+1)));
  ptr -= dimxy; // i1 * dimxy + j2 * dimx + k1
  const Target_type d = Target_type(transform(*ptr));
  const Target_type h = Target_type(transform(*(ptr+1)));
  

//   const Target_type a = ((Image_word_type*)image()->data)[i1 * dimxy + j1 * dimx + k1];
//   const Target_type b = ((Image_word_type*)image()->data)[i2 * dimxy + j1 * dimx + k1];
//   const Target_type c = ((Image_word_type*)image()->data)[i2 * dimxy + j2 * dimx + k1];
//   const Target_type d = ((Image_word_type*)image()->data)[i1 * dimxy + j2 * dimx + k1];
//   const Target_type e = ((Image_word_type*)image()->data)[i1 * dimxy + j1 * dimx + k2];
//   const Target_type f = ((Image_word_type*)image()->data)[i2 * dimxy + j1 * dimx + k2];
//   const Target_type g = ((Image_word_type*)image()->data)[i2 * dimxy + j2 * dimx + k2];
//   const Target_type h = ((Image_word_type*)image()->data)[i1 * dimxy + j2 * dimx + k2];

//   const Target_type outside = Target_type(transform(value_outside);

//   if(x < 0.f ||
//      y < 0.f ||
//      z < 0.f ||
//      i1 >= dimz ||
//      j1 >= dimy ||
//      k1 >= dimx)
//   {
//     return outside;
//   }

//   Target_type a, b, c, d, e, f, g, h;

//   if(k1 < 0) {
//     a = b = c = d = outside;
//   }
//   else {
//     if(j1 < 0) {
//       a = b = outside;
//     }
//     else {
//       if(i1 < 0)
// 	a = outside;
//       else
// 	a = ((Image_word_type*)image()->data)[i1 * dimxy + j1 * dimx + k1];

//       if(i2 >= dimz)
// 	b = outside;
//       else 
// 	b = ((Image_word_type*)image()->data)[i2 * dimxy + j1 * dimx + k1];
//     }

//     if(j2 >= dimy) {
//       c = d = outside;
//     }
//     else {
//       if(i1 < 0)
// 	d = outside;
//       else
// 	d = ((Image_word_type*)image()->data)[i1 * dimxy + j2 * dimx + k1];

//       if(i2 >= dimz)
// 	c = outside;
//       else
// 	c = ((Image_word_type*)image()->data)[i2 * dimxy + j2 * dimx + k1];
//     }
//   }

//   if(k2 >= dimx) {
//     e = f = g = h = outside;
//   }
//   else {
//     if(j1 < 0) {
//       e = f = outside;
//     }
//     else {
//       if(i1 < 0)
// 	e = outside;
//       else
// 	e = ((Image_word_type*)image()->data)[i1 * dimxy + j1 * dimx + k2];

//       if(i2 >= dimz)
// 	f = outside;
//       else 
// 	f = ((Image_word_type*)image()->data)[i2 * dimxy + j1 * dimx + k2];
//     }

//     if(j2 >= dimy) {
//       g = h = outside;
//     }
//     else {
//       if(i1 < 0)
// 	h = outside;
//       else
// 	h = ((Image_word_type*)image()->data)[i1 * dimxy + j2 * dimx + k2];

//       if(i2 >= dimz)
// 	g = outside;
//       else
// 	g = ((Image_word_type*)image()->data)[i2 * dimxy + j2 * dimx + k2];
//     }
//   }

  const Target_type di2 = i2 - lz;
  const Target_type di1 = lz - i1;
  const Target_type dj2 = j2 - ly;
  const Target_type dj1 = ly - j1;
  const Target_type dk2 = k2 - lx;
  const Target_type dk1 = lx - k1;
//   std::cerr << di2 << " " << di1 << "\n";
//   std::cerr << dj2 << " " << dj1 << "\n";
//   std::cerr << dk2 << " " << dk1 << "\n";

  return ( (  ( a * di2 + b * di1 ) * dj2 + 
	      ( d * di2 + c * di1 ) * dj1   ) * dk2 +
	   (  ( e * di2 + f * di1 ) * dj2 + 
	      ( h * di2 + g * di1 ) * dj1   ) * dk1 );
} // end trilinear_interpolation


template <typename Image_word_type,
	  typename Coord_type,
          typename Target_type,
          typename Indicator_factory>
Target_type
Image_3::labellized_trilinear_interpolation
  (const Coord_type& x,
   const Coord_type& y,
   const Coord_type& z,
   const Target_type& value_outside,
   const Indicator_factory indicator_factory) const
{
  // Check on double/float coordinates, because (int)-0.1 gives 0
  if ( x < 0 || y < 0 || z < 0 ) return value_outside;
  
  Coord_type lx = x / image()->vx;
  Coord_type ly = y / image()->vy;
  Coord_type lz = z / image()->vz;
  const std::size_t dimx = xdim();
  const std::size_t dimy = ydim();
  const std::size_t dimz = zdim();
  
  if( lx < 0 ||
      ly < 0 ||
      lz < 0 ||
     lz >= Coord_type(dimz-1) ||
     ly >= Coord_type(dimy-1) ||
     lx >= Coord_type(dimx-1))
  {
    return value_outside;
  }  

  // images are indexed by (z,y,x)
  const int i1 = (int)(lz); 
  const int j1 = (int)(ly);
  const int k1 = (int)(lx);
  const int i2 = i1 + 1;
  const int j2 = j1 + 1;

  CGAL::cpp11::array<std::size_t,8> index;
  index[0] = (i1 * dimy + j1) * dimx + k1;
  index[1] = index[0] + 1;
  index[2] = (i1 * dimy + j2) * dimx + k1;
  index[3] = index[2] + 1;
  index[4] = (i2 * dimy + j1) * dimx + k1;
  index[5] = index[4] + 1;
  index[6] = (i2 * dimy + j2) * dimx + k1;
  index[7] = index[6] + 1;

  CGAL::cpp11::array<Image_word_type,8> labels;
  
  labels[0] = ((Image_word_type*)image()->data)[index[0]];
  int lc = 1;
  for(int lci=1; lci<8; ++lci){
    bool found = false;
    Image_word_type iwt = ((Image_word_type*)image()->data)[index[lci]];
    for(int lcj=0; lcj < lc; ++lcj){
      if(iwt == labels[lcj]){
        found = true;
        break;
      }
    }
    if(found) continue;
    labels[lc] = iwt;
    ++lc;
  }

  CGAL_HISTOGRAM_PROFILER(
    "Number of labels around a vertex, Image_3::labellized_trilinear_interpolation()", 
    static_cast<unsigned int>(lc));

  if(lc == 1) {
    return labels[0];
  }
    
  double best_value = 0.;
  Image_word_type best = 0;
  for(int i = 0; i < lc; ++i)
  {
    Image_word_type iwt = labels[i];
    const double r = 
      trilinear_interpolation<Image_word_type,double,Coord_type>(
        x, y, z, value_outside, indicator_factory.indicator(iwt));

    if(r > best_value) {
      best = iwt;
      best_value = r;
    }
  }
//   CGAL_assertion(best_value > 0.5);
  return best;
}

} // end namespace CGAL

#ifdef CGAL_HEADER_ONLY
#include <CGAL/Image_3_impl.h>
#endif // CGAL_HEADER_ONLY

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif
 
#endif // CGAL_IMAGE_3_H
