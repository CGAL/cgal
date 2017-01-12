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
// Author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//                 Pedro Machado Manhaes de Castro  <pmmc@cin.ufpe.br>
//                 Alexandru Tifrea
//                 Maxime Gimeno

#ifndef CGAL_POINT_GENERATORS_3_H
#define CGAL_POINT_GENERATORS_3_H 1
#include <CGAL/generators.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/number_type_basic.h>
#include <CGAL/internal/Generic_random_point_generator.h>
#include <CGAL/boost/graph/property_maps.h>

namespace CGAL {

template < class P, class Creator = 
                   Creator_uniform_3<typename Kernel_traits<P>::Kernel::RT,P> >
class Random_points_in_sphere_3 : public Random_generator_base<P> {
    void generate_point();
public:
    typedef Random_points_in_sphere_3<P,Creator> This;
    Random_points_in_sphere_3( double r = 1, Random& rnd = CGAL::get_default_random())
        // g is an input iterator creating points of type `P' uniformly
        // distributed in the open sphere with radius r, i.e. |`*g'| < r .
        // Three random numbers are needed from `rnd' for each point
    : Random_generator_base<P>( r, rnd) { generate_point(); }
    This& operator++()    {
        generate_point();
        return *this;
    }
    This  operator++(int) {
        This tmp = *this;
        ++(*this);
        return tmp;
    }
};

template < class P, class Creator >
void
Random_points_in_sphere_3<P,Creator>::
generate_point() {
  // A strip between z and z+dz has an area independant of z
    typedef typename Creator::argument_type T;
    double alpha = this->_rnd.get_double() * 2.0 * CGAL_PI;
    double z     = 2 * this->_rnd.get_double() - 1.0;
    double r     = std::sqrt( 1 - z * z);
    r *= std::pow( this->_rnd.get_double() , 1.0/3.0 );  
    Creator creator;
    this->d_item = creator( T(this->d_range * r * std::cos(alpha)),
                            T(this->d_range * r * std::sin(alpha)),
                            T(this->d_range * z));
}


template < class P, class Creator = 
                   Creator_uniform_3<typename Kernel_traits<P>::Kernel::RT,P> >
class Random_points_on_sphere_3 : public Random_generator_base<P> {
    void generate_point();
public:
    typedef Random_points_on_sphere_3<P,Creator> This;
    Random_points_on_sphere_3( double r = 1, Random& rnd = CGAL::get_default_random())
        // g is an input iterator creating points of type `P' uniformly
        // distributed on the sphere with radius r, i.e. |`*g'| == r . A
        // two random numbers are needed from `rnd' for each point.
    : Random_generator_base<P>( r, rnd) { generate_point(); }
    This& operator++()    {
        generate_point();
        return *this;
    }
    This  operator++(int) {
        This tmp = *this;
        ++(*this);
        return tmp;
    }
};

template < class P, class Creator >
void
Random_points_on_sphere_3<P,Creator>::
generate_point() {
  // A strip between z and z+dz has an area independant of z
    typedef typename Creator::argument_type T;
    double alpha = this->_rnd.get_double() * 2.0 * CGAL_PI;
    double z     = 2 * this->_rnd.get_double() - 1.0;
    double r     = std::sqrt( 1 - z * z);
    Creator creator;
    this->d_item = creator( T(this->d_range * r * std::cos(alpha)),
                            T(this->d_range * r * std::sin(alpha)),
                            T(this->d_range * z));
}


template < class P, class Creator = 
                   Creator_uniform_3<typename Kernel_traits<P>::Kernel::RT,P> >
class Random_points_in_cube_3 : public Random_generator_base<P>{
    void generate_point();
public:
    typedef Random_points_in_cube_3<P,Creator> This;
    Random_points_in_cube_3( double a = 1, Random& rnd = CGAL::get_default_random())
    : Random_generator_base<P>( a, rnd) { generate_point(); }
    This& operator++()    {
        generate_point();
        return *this;
    }
    This  operator++(int) {
        This tmp = *this;
        ++(*this);
        return tmp;
    }
};

template < class P, class Creator >
void
Random_points_in_cube_3<P,Creator>::
generate_point() {
    typedef typename Creator::argument_type T;
    Creator creator;
    this->d_item =
	     creator( T(this->d_range * ( 2 * this->_rnd.get_double() - 1.0)),
                      T(this->d_range * ( 2 * this->_rnd.get_double() - 1.0)),
                      T(this->d_range * ( 2 * this->_rnd.get_double() - 1.0)));
}


template <class OutputIterator, class Creator>
OutputIterator
points_on_cube_grid_3( double a, std::size_t n, 
                         OutputIterator o, Creator creator)
{
    if  (n == 0)
        return o;

    int m = int(std::ceil(
                  std::sqrt(std::sqrt(static_cast<double>(n)))));

    while (m*m*m < int(n)) m++;

    double base = -a;  // Left and bottom boundary.
    double step = 2*a/(m-1);
    int j = 0;
    int k = 0;
    double px = base;
    double py = base;
    double pz = base;
    *o++ = creator( px, py, pz);
    for (std::size_t i = 1; i < n; i++) {
        j++;
        if ( j == m) {
           k++;
           if ( k == m) {
              py = base;
              px = base;
              pz = pz + step;
              k = 0;
           }
           else {
              px = base;
              py = py + step;
           }
           j = 0;
        } else {
           px = px + step;
        }
        *o++ = creator( px, py, pz);
    }
    return o;
}

template <class OutputIterator>
OutputIterator
points_on_cube_grid_3( double a, std::size_t n, OutputIterator o)
{
    typedef std::iterator_traits<OutputIterator> ITraits;
    typedef typename ITraits::value_type         P;
    return points_on_square_grid_3(a, n, o, 
                 Creator_uniform_3<typename Kernel_traits<P>::Kernel::RT,P>());
}

template < class P, class Creator = 
Creator_uniform_3<typename Kernel_traits<P>::Kernel::RT,P> >
class Random_points_in_triangle_3 : public Random_generator_base<P> {
	P _p,_q,_r;
	void generate_point();
public:
	typedef P result_type;
	typedef Random_points_in_triangle_3<P, Creator> This;
	typedef typename Kernel_traits<P>::Kernel::Triangle_3 Triangle_3;
	Random_points_in_triangle_3() {}
	Random_points_in_triangle_3( const This& x,Random& rnd = get_default_random())
	: Random_generator_base<P>( 1, rnd ),_p(x._p),_q(x._q),_r(x._r) {
		generate_point();
	}
	Random_points_in_triangle_3( const P& p, const P& q, const P& r, Random& rnd = get_default_random())
	: Random_generator_base<P>( 1, rnd ),_p(p),_q(q),_r(r) {
		generate_point();
	}
	Random_points_in_triangle_3( const Triangle_3& triangle,Random& rnd = get_default_random())
	: Random_generator_base<P>( 1,
			rnd),_p(triangle[0]),_q(triangle[1]),_r(triangle[2]) {
		generate_point();
	}
	This& operator++() {
		generate_point();
		return *this;
	}
	This operator++(int) {
		This tmp = *this;
		++(*this);
		return tmp;
	}
};


template<class P, class Creator >
void Random_points_in_triangle_3<P, Creator>::generate_point() {
	typedef typename Creator::argument_type T;
	Creator creator;
	double a1 = this->_rnd.get_double(0,1);
	double a2 = this->_rnd.get_double(0,1);
	if(a1>a2) std::swap(a1,a2);
	double b1 = a1;
	double b2 = a2-a1;
	double b3 = 1.0-a2;
	T ret[3];
	for(int i = 0; i < 3; ++i) {
	    ret[i] = T(to_double(_p[i])*b1+to_double(_q[i])*b2+to_double(_r[i])*b3);
	}
	this->d_item = creator(ret[0],ret[1],ret[2]);
}

template < class P, class Creator = 
                   Creator_uniform_3<typename Kernel_traits<P>::Kernel::RT,P> >
class Random_points_on_segment_3 : public Random_generator_base<P> {
    P _p;
    P _q;
    void generate_point();
public:
    typedef Random_points_on_segment_3<P,Creator> This;
    Random_points_on_segment_3() {}
    Random_points_on_segment_3( const P& p,
                                const P& q,
                                Random& rnd = CGAL::get_default_random())
        // g is an input iterator creating points of type `P' uniformly
        // distributed on the segment from p to q except q, i.e. `*g' ==
        // \lambda p + (1-\lambda)\, q where 0 <= \lambda < 1 . A single
        // random number is needed from `rnd' for each point.
      : Random_generator_base<P>( (std::max)( (std::max)( (std::max)(to_double(p.x()), to_double(q.x())),
                                                         (std::max)(to_double(p.y()), to_double(q.y()))),
                                                         (std::max)(to_double(p.z()), to_double(q.z()))),
                                  rnd) , _p(p), _q(q)
    {
      generate_point();
    }

    template <class Segment_3>
    Random_points_on_segment_3( const Segment_3& s,
                                Random& rnd = CGAL::get_default_random())
        // g is an input iterator creating points of type `P' uniformly
        // distributed on the segment from p to q except q, i.e. `*g' ==
        // \lambda p + (1-\lambda)\, q where 0 <= \lambda < 1 . A single
        // random number is needed from `rnd' for each point.
      : Random_generator_base<P>( (std::max)( (std::max)( (std::max)(to_double(s[0].x()), to_double(s[1].x())),
                                                         (std::max)(to_double(s[0].y()), to_double(s[1].y()))),
                                                         (std::max)(to_double(s[0].z()), to_double(s[1].z()))),
                                  rnd) , _p(s[0]), _q(s[1])
    {
      generate_point();
    }

    const P&  source() const { return _p; }
    const P&  target() const { return _q; }
    This& operator++()    {
        generate_point();
        return *this;
    }
    This  operator++(int) {
        This tmp = *this;
        ++(*this);
        return tmp;
    }
};

template < class P, class Creator >
void
Random_points_on_segment_3<P,Creator>::
generate_point() {
    typedef typename Creator::argument_type  T;
    double la = this->_rnd.get_double();
    double mu = 1.0 - la;
    Creator creator;
    this->d_item = creator(T(mu * to_double(_p.x()) + la * to_double(_q.x())),
                           T(mu * to_double(_p.y()) + la * to_double(_q.y())),
                           T(mu * to_double(_p.z()) + la * to_double(_q.z())));
}

template < class P, class Creator =
Creator_uniform_3<typename Kernel_traits<P>::Kernel::RT,P> >
class Random_points_in_tetrahedron_3 : public Random_generator_base<P> {
	P _p,_q,_r,_s;
	void generate_point();
public:
	typedef P result_type;
	typedef Random_points_in_tetrahedron_3<P, Creator> This;
	typedef typename Kernel_traits<P>::Kernel::Tetrahedron_3 Tetrahedron_3;
	Random_points_in_tetrahedron_3() {}
	Random_points_in_tetrahedron_3( const This& x,Random& rnd = get_default_random())
	: Random_generator_base<P>( 1, rnd ),_p(x._p),_q(x._q),_r(x._r),_s(x._s) {
		generate_point();
	}
	Random_points_in_tetrahedron_3( const P& p, const P& q, const P& r, const P& s,Random& rnd = get_default_random())
	: Random_generator_base<P>( 1, rnd ),_p(p),_q(q),_r(r),_s(s) {
		generate_point();
	}
	Random_points_in_tetrahedron_3( const Tetrahedron_3& tetrahedron,Random& rnd = get_default_random())
	: Random_generator_base<P>( 1, rnd),_p(tetrahedron[0]),_q(tetrahedron[1]),_r(tetrahedron[2]),_s(tetrahedron[3]) {
		generate_point();
	}
	This& operator++() {
		generate_point();
		return *this;
	}
	This operator++(int) {
		This tmp = *this;
		++(*this);
		return tmp;
	}
};

template<class P, class Creator >
void Random_points_in_tetrahedron_3<P, Creator>::generate_point() {
	typedef typename Creator::argument_type T;
	Creator creator;
	double a[3];
	for(int i = 0; i < 3; ++i) {
		a[i]=this->_rnd.get_double(0,1);
	}
	std::sort(a,a+3);
	double b[4];
	b[0]=a[0];
	b[1]=a[1]-a[0];
	b[2]=a[2]-a[1];
	b[3]=1.0-a[2];
	T ret[3];
	for(int i = 0; i < 3; ++i) {
	    ret[i] = T(to_double(_p[i])*b[0]+to_double(_q[i])*b[1]+to_double(_r[i])*b[2]+to_double(_s[i])*b[3]);
	}
	this->d_item = creator(ret[0],ret[1],ret[2]);
}



template <class TriangleMesh,
          class VertexPointMap = typename boost::property_map<TriangleMesh,
                                                              CGAL::vertex_point_t>::const_type,
          class Creator = Creator_uniform_3<
                            typename Kernel_traits< typename boost::property_traits<VertexPointMap>::value_type >::Kernel::RT,
                            typename boost::property_traits<VertexPointMap>::value_type >
>
struct Random_points_in_triangle_mesh_3
  : public Generic_random_point_generator<
             typename boost::graph_traits <TriangleMesh>::face_descriptor ,
             CGAL::Property_map_to_unary_function<CGAL::Triangle_from_face_descriptor_map<
                                                    TriangleMesh, VertexPointMap > >,
             Random_points_in_triangle_3<typename boost::property_traits<VertexPointMap>::value_type, Creator>,
             typename boost::property_traits<VertexPointMap>::value_type>
{
  typedef typename boost::property_traits<VertexPointMap>::value_type  P;
  typedef Generic_random_point_generator<
            typename boost::graph_traits <TriangleMesh>::face_descriptor ,
            CGAL::Property_map_to_unary_function<CGAL::Triangle_from_face_descriptor_map<
            TriangleMesh,VertexPointMap> >,
            Random_points_in_triangle_3<P, Creator> , P>                   Base;
  typedef typename CGAL::Triangle_from_face_descriptor_map<
                     TriangleMesh,VertexPointMap>                          Object_from_id;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor      Id;
  typedef P result_type;
  typedef Random_points_in_triangle_mesh_3< TriangleMesh, VertexPointMap, Creator>  This;


  Random_points_in_triangle_mesh_3( const TriangleMesh& mesh,Random& rnd = get_default_random())
    : Base( faces(mesh),
            CGAL::Property_map_to_unary_function<Object_from_id>(Object_from_id(&mesh, get(vertex_point, mesh))),
            internal::Apply_approx_sqrt<typename Kernel_traits<P>::Kernel::Compute_squared_area_3>(),
            rnd )
  {
  }
  Random_points_in_triangle_mesh_3( const TriangleMesh& mesh, VertexPointMap vpm, Random& rnd = get_default_random())
    : Base( faces(mesh),
            CGAL::Property_map_to_unary_function<Object_from_id>(Object_from_id(&mesh, vpm)),
            internal::Apply_approx_sqrt<typename Kernel_traits<P>::Kernel::Compute_squared_area_3>(),
            rnd )
  {
  }
  This& operator++() {
    Base::generate_point();
    return *this;
  }
  This operator++(int) {
    This tmp = *this;
    ++(*this);
    return tmp;
  }
  double mesh_area() const
  {
    return this->sum_of_weights();
  }
};

template <class EdgeListGraph,
          class VertexPointMap = typename boost::property_map<EdgeListGraph,
                                                              CGAL::vertex_point_t>::const_type,
          class Creator = Creator_uniform_3<
                            typename Kernel_traits< typename boost::property_traits<VertexPointMap>::value_type >::Kernel::RT,
                            typename boost::property_traits<VertexPointMap>::value_type >
>
struct Random_points_on_edge_list_graph_3
  : public Generic_random_point_generator<
             typename boost::graph_traits <EdgeListGraph>::edge_descriptor,
             CGAL::Property_map_to_unary_function<CGAL::Segment_from_edge_descriptor_map<
                                                  EdgeListGraph, VertexPointMap > >,
             Random_points_on_segment_3<typename boost::property_traits<VertexPointMap>::value_type, Creator>,
             typename boost::property_traits<VertexPointMap>::value_type>
{
  typedef typename boost::property_traits<VertexPointMap>::value_type  P;
  typedef Generic_random_point_generator<
            typename boost::graph_traits <EdgeListGraph>::edge_descriptor,
            CGAL::Property_map_to_unary_function<CGAL::Segment_from_edge_descriptor_map<
                                                 EdgeListGraph, VertexPointMap > >,
            Random_points_on_segment_3<P, Creator> , P>                    Base;
  typedef typename CGAL::Segment_from_edge_descriptor_map<
                     EdgeListGraph,VertexPointMap>                          Object_from_id;
  typedef typename boost::graph_traits<EdgeListGraph>::edge_descriptor      Id;
  typedef P result_type;
  typedef Random_points_on_edge_list_graph_3< EdgeListGraph, VertexPointMap, Creator>  This;

  Random_points_on_edge_list_graph_3( const EdgeListGraph& mesh,Random& rnd = get_default_random())
    : Base( edges(mesh),
            CGAL::Property_map_to_unary_function<Object_from_id>(Object_from_id(&mesh, get(vertex_point, mesh))),
            internal::Apply_approx_sqrt<typename Kernel_traits<P>::Kernel::Compute_squared_length_3>(),
            rnd )
  {
  }
  Random_points_on_edge_list_graph_3( const EdgeListGraph& mesh, VertexPointMap vpm, Random& rnd = get_default_random())
    : Base( edges(mesh),
            CGAL::Property_map_to_unary_function<Object_from_id>(Object_from_id(&mesh, vpm)),
            internal::Apply_approx_sqrt<typename Kernel_traits<P>::Kernel::Compute_squared_length_3>(),
            rnd )
  {
  }
  This& operator++() {
    Base::generate_point();
    return *this;
  }
  This operator++(int) {
    This tmp = *this;
    ++(*this);
    return tmp;
  }
  double mesh_length() const
  {
    return this->sum_of_weights();
  }
};

namespace internal
{

template<class T>
class Triangle_from_face_C3t3
{
  typedef typename T::Triangle                           Triangle;
  typedef typename T::Point                              Point;
  typedef std::pair<typename T::Cell_handle, int> Face;
public:
  typedef Triangle result_type;

  Triangle_from_face_C3t3()
  {}
  Triangle operator()(Face face)const
  {
    typename T::Cell_handle cell = face.first;
    int index = face.second;
    const Point& pa = cell->vertex((index+1)&3)->point();
    const Point& pb = cell->vertex((index+2)&3)->point();
    const Point& pc = cell->vertex((index+3)&3)->point();
    return Triangle(pa, pb, pc);
  }
};

template<class T>
class Tetrahedron_from_cell_C3t3
{
  typedef typename T::Cell_handle                              Cell;
  typedef typename T::Point                                    Point;
  typedef typename Kernel_traits<Point>::Kernel::Tetrahedron_3 Tetrahedron;
public:
  typedef Tetrahedron result_type;

  Tetrahedron_from_cell_C3t3()
  {}
  Tetrahedron operator()(Cell cell)const
  {
    const Point&  p0 = cell->vertex(0)->point();
    const Point&  p1 = cell->vertex(1)->point();
    const Point&  p2 = cell->vertex(2)->point();
    const Point&  p3 = cell->vertex(3)->point();
    return Tetrahedron(p0,p1,p2,p3);
  }
};
}//end namespace internal

template <class C3t3,
          class Creator = Creator_uniform_3<
                            typename Kernel_traits< typename C3t3::Point >::Kernel::RT,
                            typename C3t3::Point >
>
struct Random_points_in_tetrahedral_mesh_boundary_3
  : public Generic_random_point_generator<
             std::pair<typename C3t3::Triangulation::Cell_handle, int>,
             internal::Triangle_from_face_C3t3<typename C3t3::Triangulation>,
             Random_points_in_triangle_3<typename C3t3::Point>,
             typename C3t3::Point>
{
  typedef Generic_random_point_generator<
            std::pair<typename C3t3::Triangulation::Cell_handle, int>,
            internal::Triangle_from_face_C3t3<typename C3t3::Triangulation>,
            Random_points_in_triangle_3<typename C3t3::Point, Creator>,
            typename C3t3::Point>                                    Base;
  typedef std::pair<typename C3t3::Triangulation::Cell_handle, int>  Id;
  typedef typename C3t3::Point                                       result_type;
  typedef Random_points_in_tetrahedral_mesh_boundary_3<C3t3, Creator> This;


  Random_points_in_tetrahedral_mesh_boundary_3( const C3t3& c3t3,Random& rnd = get_default_random())
    : Base( make_range( c3t3.facets_in_complex_begin(),
                        c3t3.facets_in_complex_end()),
            internal::Triangle_from_face_C3t3<typename C3t3::Triangulation>(),
            internal::Apply_approx_sqrt<typename Kernel_traits<typename C3t3::Point>::Kernel::Compute_squared_area_3>(),
            rnd )
  {
  }
  This& operator++() {
    Base::generate_point();
    return *this;
  }
  This operator++(int) {
    This tmp = *this;
    ++(*this);
    return tmp;
  }
};

template <class C3t3,
          class Creator = Creator_uniform_3<
                            typename Kernel_traits< typename C3t3::Point >::Kernel::RT,
                            typename C3t3::Point >
>
struct Random_points_in_tetrahedral_mesh_3
  : public Generic_random_point_generator<
             typename C3t3::Triangulation::Cell_handle,
             internal::Tetrahedron_from_cell_C3t3<typename C3t3::Triangulation>,
             Random_points_in_tetrahedron_3<typename C3t3::Point>,
             typename C3t3::Point>
{
  typedef Generic_random_point_generator<
            typename C3t3::Triangulation::Cell_handle,
            internal::Tetrahedron_from_cell_C3t3<typename C3t3::Triangulation>,
            Random_points_in_tetrahedron_3<typename C3t3::Point, Creator>,
            typename C3t3::Point>                                              Base;
  typedef typename C3t3::Triangulation::Cell_handle                            Id;
  typedef typename C3t3::Point                                                 result_type;
  typedef Random_points_in_tetrahedral_mesh_3<C3t3, Creator>                   This;


  Random_points_in_tetrahedral_mesh_3( const C3t3& c3t3,Random& rnd = get_default_random())
    : Base( CGAL::make_prevent_deref_range(c3t3.cells_in_complex_begin(),
                                           c3t3.cells_in_complex_end()),
            internal::Tetrahedron_from_cell_C3t3<typename C3t3::Triangulation>(),
            typename Kernel_traits<typename C3t3::Point>::Kernel::Compute_volume_3(),
            rnd )
  {
  }
  This& operator++() {
    Base::generate_point();
    return *this;
  }
  This operator++(int) {
    This tmp = *this;
    ++(*this);
    return tmp;
  }
};


template <class Point_3,
          class Triangle_3=typename Kernel_traits<Point_3>::Kernel::Triangle_3,
          class Creator = Creator_uniform_3<
                            typename Kernel_traits< Point_3 >::Kernel::RT,
                            Point_3 >
         >
struct Random_points_in_triangles_3
    : public Generic_random_point_generator<const Triangle_3*,
                                            internal::Deref<Triangle_3>,
                                            Random_points_in_triangle_3<Point_3>,
                                            Point_3>
{
  typedef Generic_random_point_generator<const Triangle_3*,
                                         internal::Deref<Triangle_3>,
                                         Random_points_in_triangle_3<Point_3, Creator>,
                                         Point_3>            Base;
  typedef const Triangle_3*                                         Id;
  typedef Point_3                                                   result_type;
  typedef Random_points_in_triangles_3<Point_3, Triangle_3, Creator>  This;

  template<typename TriangleRange>
  Random_points_in_triangles_3( const TriangleRange& triangles, Random& rnd = get_default_random())
    : Base(make_range( boost::make_transform_iterator(triangles.begin(), internal::Address_of<Triangle_3>()),
                       boost::make_transform_iterator(triangles.end(), internal::Address_of<Triangle_3>()) ),
           internal::Deref<Triangle_3>(),
           internal::Apply_approx_sqrt<typename Kernel_traits<Point_3>::Kernel::Compute_squared_area_3>()
           ,rnd )
  {
  }
  This& operator++() {
    Base::generate_point();
    return *this;
  }
  This operator++(int) {
    This tmp = *this;
    ++(*this);
    return tmp;
  }
};

} //namespace CGAL


#endif // CGAL_POINT_GENERATORS_3_H //
// EOF //
