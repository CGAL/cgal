// Copyright (c) 2011  GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Sebastien Loriot
// 


#ifndef CGAL_RECONSTRUCTION_FROM_PARALLEL_SLICES_3_CONTOUR_PROVIDERS_H
#define CGAL_RECONSTRUCTION_FROM_PARALLEL_SLICES_3_CONTOUR_PROVIDERS_H

namespace CGAL{

template <class Point_2,unsigned int constant_coordinate>
struct All_polygons_in_one_file_axis_aligned_planes{
  CGAL_static_assertion(constant_coordinate<3);
  
  std::list< std::pair<std::vector<Point_2>,double> > polygons;
  typename std::list< std::pair<std::vector<Point_2>,double> >::iterator current;
  typename std::vector<Point_2>::iterator itpt;
  All_polygons_in_one_file_axis_aligned_planes(char* fname)
  {
    std::ifstream input(fname);
    int nbpt;
    double last_elevation;
    double coords[3];
    int layer_nb=-1;
    bool order_check;
    
    while (input >> nbpt && input){
      input >> coords[0] >> coords[1] >> coords[2];
      
      if (layer_nb==-1) ++layer_nb;
      else
        if( last_elevation!= coords[constant_coordinate] ){
          if (layer_nb==0) order_check = (last_elevation < coords[constant_coordinate]);
          else assert( order_check == (last_elevation < coords[constant_coordinate]) );
          ++layer_nb;
        }
      last_elevation=coords[constant_coordinate];
      
      polygons.push_back(std::make_pair(std::vector<Point_2>(),coords[constant_coordinate]));
      std::vector<Point_2>& points=polygons.back().first;
      points.reserve(nbpt);
  
      points.push_back(Point_2(coords[(constant_coordinate+1)%3],coords[(constant_coordinate+2)%3]));
      for (int i=1;i<nbpt;++i){
        input >> coords[0] >> coords[1] >> coords[2];
        points.push_back(Point_2(coords[(constant_coordinate+1)%3],coords[(constant_coordinate+2)%3]));  
      }
    }
    current=polygons.begin();
    itpt=current->first.begin();
  }
  
  bool empty() const{
    return current==polygons.end();
  }

  bool has_next_planar_contour() const {
    return CGAL::cpp0x::next(current)!=polygons.end();
  }
  
  void next_polygon(){
    assert(current!=polygons.end());
    ++current;
    assert(current!=polygons.end());
    itpt=current->first.begin();
  }
  
  std::size_t number_of_points() const {
    return current->first.size();
  }
  
  bool has_another_component() {
    typename std::list< std::pair<std::vector<Point_2>,double > >::iterator next_one=CGAL::cpp0x::next(current);  
    return next_one!=polygons.end() && next_one->second==current->second;
  }
  
  Point_2 get_point(){ return *itpt++; }

  Point_2 get_point(double& z){
    z=current->second;
    return *itpt++;
  }
};
  
template <class PolygonIterator,class CoordIterator>
class Polygon_range_in_axis_aligned_planes{
  typedef typename std::iterator_traits<PolygonIterator>::value_type Polygon_2;
  typedef typename Polygon_2::Traits::Point_2 Point_2;
  typedef typename Polygon_2::Vertex_circulator Vertex_circulator;
  
  PolygonIterator m_next,m_end;
  CoordIterator m_z_next,m_z_end;
  double m_z;
  Polygon_2 m_current_polygon;
  Vertex_circulator itpt;


public:
  void next_polygon(){
    assert(m_next!=m_end);
    assert(m_z_next!=m_z_end);
    m_current_polygon=*m_next;
    ++m_next;
    m_z=*m_z_next;
    ++m_z_next;
    itpt=m_current_polygon.vertices_circulator();
  }
  
  Polygon_range_in_axis_aligned_planes(PolygonIterator begin,PolygonIterator end,CoordIterator z_begin,CoordIterator z_end):
    m_next(begin),m_end(end),m_z_next(z_begin),m_z_end(z_end)
  {
    if (m_next!=m_end)
      next_polygon();
  }
  
  bool empty() const{
    return m_next==m_end && m_current_polygon.is_empty();
  }

  bool has_next_planar_contour() const {
    return m_next!=m_end;
  }
  
  std::size_t number_of_points() const {
    return m_current_polygon.size()+1;
  }
  
  bool has_another_component() {
    return m_next!=m_end && m_z==*m_z_next;
  }
  
  Point_2 get_point(){ return *itpt++; }

  Point_2 get_point(double& z){
    z=m_z;
    return *itpt++;
  }
};


template <class Kernel>
class Polygon_as_vector_of_Point_3_in_axis_aligned_planes{
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Point_2 Point_2;
  
  const std::vector< std::vector<Point_3> >& m_contours;
  
  unsigned int m_current_polygon,m_current_point,m_cst_coord;
  
public:
  void next_polygon(){
    ++m_current_polygon;
    m_current_point=0;
  }
  
  Polygon_as_vector_of_Point_3_in_axis_aligned_planes(const std::vector< std::vector<Point_3> >& contours,int cst_coord):
    m_contours(contours),m_current_polygon(0),m_current_point(0),m_cst_coord(cst_coord)
  {}
  
  bool empty() const{
    return m_contours.empty();
  }

  bool has_next_planar_contour() const {
    return m_current_polygon!=m_contours.size()-1;
  }
  
  std::size_t number_of_points() const {
    return m_contours[m_current_polygon].size();
  }
  
  bool has_another_component() {
    return m_current_polygon!=m_contours.size()-1 && 
           m_contours[m_current_polygon+1][0][m_cst_coord]==m_contours[m_current_polygon][0][m_cst_coord];
  }
  
  Point_2 get_point(){ 
    double x = m_contours[m_current_polygon][m_current_point][(m_cst_coord+1)%3];
    double y = m_contours[m_current_polygon][m_current_point][(m_cst_coord+2)%3];
    ++m_current_point;
    return Point_2(x,y); 
  }

  Point_2 get_point(double& z){
    z=m_contours[m_current_polygon][0][m_cst_coord];
    return get_point();
  }
};




} //namespace CGAL

#endif //CGAL_RECONSTRUCTION_FROM_PARALLEL_SLICES_3_CONTOUR_PROVIDERS_H
