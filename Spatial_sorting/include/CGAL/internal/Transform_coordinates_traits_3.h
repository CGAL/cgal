// Copyright (c) 2013  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Olivier Devillers
//                 Pedro Machado Manhaes de Castro

#ifndef CGAL_INTERNAL_SCALE_COORDINATES_ADAPTOR_TRAITS_3_H
#define CGAL_INTERNAL_SCALE_COORDINATES_ADAPTOR_TRAITS_3_H

namespace CGAL { 

namespace internal {
	
template <int x, int y, int z, int ord>
struct Transform_constant_struct;

template <> struct Transform_constant_struct<1,1,0,0> { enum {value = 0}; };
template <> struct Transform_constant_struct<-1,1,0,0> { enum {value = 1}; };
template <> struct Transform_constant_struct<1,-1,0,0> { enum {value = 2}; };
template <> struct Transform_constant_struct<-1,-1,0,0> { enum {value = 3}; };
template <> struct Transform_constant_struct<1,1,0,1> { enum {value = 4}; };
template <> struct Transform_constant_struct<-1,1,0,1> { enum {value = 5}; };
template <> struct Transform_constant_struct<1,-1,0,1> { enum {value = 6}; };
template <> struct Transform_constant_struct<-1,-1,0,1> { enum {value = 7}; };

template <> struct Transform_constant_struct<1,0,1,0> { enum {value = 8}; };
template <> struct Transform_constant_struct<-1,0,1,0> { enum {value = 9}; };
template <> struct Transform_constant_struct<1,0,-1,0> { enum {value = 10}; };
template <> struct Transform_constant_struct<-1,0,-1,0> { enum {value = 11}; };
template <> struct Transform_constant_struct<1,0,1,1> { enum {value = 12}; };
template <> struct Transform_constant_struct<-1,0,1,1> { enum {value = 13}; };
template <> struct Transform_constant_struct<1,0,-1,1> { enum {value = 14}; };
template <> struct Transform_constant_struct<-1,0,-1,1> { enum {value = 15}; };

template <> struct Transform_constant_struct<0,1,1,0> { enum {value = 16}; };
template <> struct Transform_constant_struct<0,-1,1,0> { enum {value = 17}; };
template <> struct Transform_constant_struct<0,1,-1,0> { enum {value = 18}; };
template <> struct Transform_constant_struct<0,-1,-1,0> { enum {value = 19}; };
template <> struct Transform_constant_struct<0,1,1,1> { enum {value = 20}; };
template <> struct Transform_constant_struct<0,-1,1,1> { enum {value = 21}; };
template <> struct Transform_constant_struct<0,1,-1,1> { enum {value = 22}; };
template <> struct Transform_constant_struct<0,-1,-1,1> { enum {value = 23}; };

template <class R, int opt>
struct Coordinate_value_adaptor;

template <class R>
struct Coordinate_value_adaptor<R,0> {
	static typename R::FT x(const typename R::Point_3& p) {return p.x();}
	static typename R::FT y(const typename R::Point_3& p) {return p.y();}
};

template <class R>
struct Coordinate_value_adaptor<R,1> {
	static typename R::FT x(const typename R::Point_3& p) {return -p.x();}
	static typename R::FT y(const typename R::Point_3& p) {return p.y();}
};

template <class R>
struct Coordinate_value_adaptor<R,2> {
	static typename R::FT x(const typename R::Point_3& p) {return p.x();}
	static typename R::FT y(const typename R::Point_3& p) {return -p.y();}
};

template <class R>
struct Coordinate_value_adaptor<R,3> {
	static typename R::FT x(const typename R::Point_3& p) {return -p.x();}
	static typename R::FT y(const typename R::Point_3& p) {return -p.y();}
};

template <class R>
struct Coordinate_value_adaptor<R,4> {
	static typename R::FT x(const typename R::Point_3& p) {return p.y();}
	static typename R::FT y(const typename R::Point_3& p) {return p.x();}
};

template <class R>
struct Coordinate_value_adaptor<R,5> {
	static typename R::FT x(const typename R::Point_3& p) {return -p.y();}
	static typename R::FT y(const typename R::Point_3& p) {return p.x();}
};

template <class R>
struct Coordinate_value_adaptor<R,6> {
	static typename R::FT x(const typename R::Point_3& p) {return p.y();}
	static typename R::FT y(const typename R::Point_3& p) {return -p.x();}
};

template <class R>
struct Coordinate_value_adaptor<R,7> {
	static typename R::FT x(const typename R::Point_3& p) {return -p.y();}
	static typename R::FT y(const typename R::Point_3& p) {return -p.x();}
};

template <class R>
struct Coordinate_value_adaptor<R,8> {
	static typename R::FT x(const typename R::Point_3& p) {return p.x();}
	static typename R::FT y(const typename R::Point_3& p) {return p.z();}
};

template <class R>
struct Coordinate_value_adaptor<R,9> {
	static typename R::FT x(const typename R::Point_3& p) {return -p.x();}
	static typename R::FT y(const typename R::Point_3& p) {return p.z();}
};

template <class R>
struct Coordinate_value_adaptor<R,10> {
	static typename R::FT x(const typename R::Point_3& p) {return p.x();}
	static typename R::FT y(const typename R::Point_3& p) {return -p.z();}
};

template <class R>
struct Coordinate_value_adaptor<R,11> {
	static typename R::FT x(const typename R::Point_3& p) {return -p.x();}
	static typename R::FT y(const typename R::Point_3& p) {return -p.z();}
};

template <class R>
struct Coordinate_value_adaptor<R,12> {
	static typename R::FT x(const typename R::Point_3& p) {return p.z();}
	static typename R::FT y(const typename R::Point_3& p) {return p.x();}
};

template <class R>
struct Coordinate_value_adaptor<R,13> {
	static typename R::FT x(const typename R::Point_3& p) {return -p.z();}
	static typename R::FT y(const typename R::Point_3& p) {return p.x();}
};

template <class R>
struct Coordinate_value_adaptor<R,14> {
	static typename R::FT x(const typename R::Point_3& p) {return p.z();}
	static typename R::FT y(const typename R::Point_3& p) {return -p.x();}
};

template <class R>
struct Coordinate_value_adaptor<R,15> {
	static typename R::FT x(const typename R::Point_3& p) {return -p.z();}
	static typename R::FT y(const typename R::Point_3& p) {return -p.x();}
};

template <class R>
struct Coordinate_value_adaptor<R,16> {
	static typename R::FT x(const typename R::Point_3& p) {return p.y();}
	static typename R::FT y(const typename R::Point_3& p) {return p.z();}
};

template <class R>
struct Coordinate_value_adaptor<R,17> {
	static typename R::FT x(const typename R::Point_3& p) {return -p.y();}
	static typename R::FT y(const typename R::Point_3& p) {return p.z();}
};

template <class R>
struct Coordinate_value_adaptor<R,18> {
	static typename R::FT x(const typename R::Point_3& p) {return p.y();}
	static typename R::FT y(const typename R::Point_3& p) {return -p.z();}
};

template <class R>
struct Coordinate_value_adaptor<R,19> {
	static typename R::FT x(const typename R::Point_3& p) {return -p.y();}
	static typename R::FT y(const typename R::Point_3& p) {return -p.z();}
};

template <class R>
struct Coordinate_value_adaptor<R,20> {
	static typename R::FT x(const typename R::Point_3& p) {return p.z();}
	static typename R::FT y(const typename R::Point_3& p) {return p.y();}
};

template <class R>
struct Coordinate_value_adaptor<R,21> {
	static typename R::FT x(const typename R::Point_3& p) {return -p.z();}
	static typename R::FT y(const typename R::Point_3& p) {return p.y();}
};

template <class R>
struct Coordinate_value_adaptor<R,22> {
	static typename R::FT x(const typename R::Point_3& p) {return p.z();}
	static typename R::FT y(const typename R::Point_3& p) {return -p.y();}
};

template <class R>
struct Coordinate_value_adaptor<R,23> {
	static typename R::FT x(const typename R::Point_3& p) {return -p.z();}
	static typename R::FT y(const typename R::Point_3& p) {return -p.y();}
};

template <class R, int opt>
class Compute_x_2
{
public:
  typedef typename R::Point_3     Point; 
  typename R::FT x(const Point &p) const { return Coordinate_value_adaptor<R,opt>::x(p); }
  typename R::FT operator()(const Point& p) const { return x(p); }
};

template <class R, int opt>
class Compute_y_2 
{
public:
  typedef typename R::Point_3     Point; 
  typename R::FT y(const Point &p) const { return Coordinate_value_adaptor<R,opt>::y(p); }
  typename R::FT operator()(const Point& p) const { return y(p); }
};

template <class R, int opt>
class Less_x_2
{
public:
  typedef typename R::Point_3     Point; 
  typename R::FT x(const Point &p) const { return Coordinate_value_adaptor<R,opt>::x(p); }
  bool operator()(const Point& p, const Point& q) const { return x(p) < x(q); }
};

template <class R, int opt>
class Less_y_2
{
public:
  typedef typename R::Point_3     Point; 
  typename R::FT y(const Point &p) const { return Coordinate_value_adaptor<R,opt>::y(p); }
  bool operator()(const Point& p, const Point& q) const { return y(p) < y(q); }
};

template <class R, int x, int y, int z, int ord>
struct Transform_coordinates_traits_3 {
	private:
		enum {opt = Transform_constant_struct<x,y,z,ord>::value};
	
	public:
		typedef Transform_coordinates_traits_3<R,x,y,z,ord>  	  		Traits;
		typedef R                                               		Rp;
		typedef typename Rp::Point_3                            		Point_2;
		typedef Less_x_2<R,opt>        									Less_x;
		typedef Less_y_2<R,opt>       									Less_y;
		typedef Compute_x_2<R,opt>     									Compute_x;
		typedef Compute_y_2<R,opt>     									Compute_y;
		
		Transform_coordinates_traits_3(){}
		Transform_coordinates_traits_3(const Transform_coordinates_traits_3&){}
				
		Less_x less_x_2_object() const { return Less_x(); }
		Less_y less_y_2_object() const { return Less_y(); }
		Compute_x compute_x_2_object() const { return Compute_x(); }
		Compute_y compute_y_2_object() const { return Compute_y(); }
};
  

} } //namespace CGAL::internal

#endif // CGAL_INTERNAL_SCALE_COORDINATES_ADAPTOR_TRAITS_3_H
