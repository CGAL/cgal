// Copyright (c) 2014  GeometryFactory (France).
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
// Author(s)     : Andreas Fabri and Sebastien Loriot
//


#ifndef CGAL_RECONSTRUCTION_FROM_PARALLEL_SLICES_3_CHECK_AND_FIX_INPUT_H
#define CGAL_RECONSTRUCTION_FROM_PARALLEL_SLICES_3_CHECK_AND_FIX_INPUT_H

#include <CGAL/Polygon_2.h>
#include <CGAL/intersections.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/function_objects.h>
#include <boost/shared_ptr.hpp>

#include <set>
#include <vector>

namespace CGAL{

namespace Reconstruction_from_parallel_slices {

enum Error_code {
  VALID_OR_TRIVIALLY_FIXED_POLYGON,
  INCONSISTENT_CST_COORD,
  DEGENERATE_POLYGON,
  POLYGON_NOT_SIMPLE,
  LOGIC_ERROR };

template <class Kernel>
Error_code
create_polygon(std::vector<typename Kernel::Point_3>& points_3d,
               int constant_coordinate, bool verbose)
{
  std::size_t nb_pts=points_3d.size();

  //remove consecutive duplicated points
  for (std::size_t i=1; i< nb_pts; ++i)
  {
    if ( points_3d[i][constant_coordinate]!=
         points_3d.front()[constant_coordinate] )
    {
      return INCONSISTENT_CST_COORD;
    }
    if ( points_3d[i]==points_3d[i-1] )
    {
      if (verbose) std::cerr << "Warning: removing duplicated point...\n";
      points_3d.erase(cpp11::next(points_3d.begin(),i));
      --i;
      --nb_pts;
      continue;
    }
  }

  if ( points_3d.front()!=points_3d.back() )
  {
    if (verbose) std::cerr << "Warning: closing polygon...\n";
    points_3d.push_back(points_3d.front());
    ++nb_pts;
  }

  //handle collinear points not correctly ordered
  /// \todo copy the code of Andreas
  /// \todo if a polygon contains only collinear points handle it as degenerate

  //now create the polygon
  Polygon_2<Kernel> polygon;

  for (std::size_t i=0; i< nb_pts-1; ++i)
  {
    typename Kernel::Point_2 pt( points_3d[i][(constant_coordinate+1)%3],
                                 points_3d[i][(constant_coordinate+2)%3] );
    polygon.push_back(pt);
  }

  if (polygon.size() < 3)
  {
    return DEGENERATE_POLYGON;
  }

  if ( !polygon.is_simple() )
  {
    return POLYGON_NOT_SIMPLE;
  }

  return VALID_OR_TRIVIALLY_FIXED_POLYGON;
}


////////////////////////////////

namespace internal{

template<class Kernel>
struct Box_with_segment_and_polygon_id
  : public Box_intersection_d::Box_d< typename Kernel::FT, 2, Box_intersection_d::ID_EXPLICIT>
{
  std::size_t polygon_index;
  typename Kernel::Segment_2 segment;

  typedef Box_intersection_d::Box_d< typename Kernel::FT, 2, Box_intersection_d::ID_EXPLICIT> Base;
  Box_with_segment_and_polygon_id(const typename Kernel::Segment_2& s, std::size_t i)
    :Base(s.bbox()),polygon_index(i), segment(s)
  {}

  Box_with_segment_and_polygon_id(const typename Kernel::Point_2& p1,
                                  const typename Kernel::Point_2& p2,
                                  std::size_t i)
    :Base(p1.bbox()+p2.bbox()),polygon_index(i), segment(p1,p2)
  {}
};

template <class Kernel>
struct Report_inters{
  std::set<std::pair<std::size_t,std::size_t> >& m_intersecting_polygons;
  typedef Box_with_segment_and_polygon_id<Kernel> Box;

  Report_inters(std::set<std::pair<std::size_t,std::size_t> >& intersecting_polygons)
    : m_intersecting_polygons(intersecting_polygons){}

  void operator() ( const Box& a, const Box& b) {
    if (a.polygon_index==b.polygon_index) return;

    if( do_intersect(a.segment,b.segment) ){
        m_intersecting_polygons.insert(
          a.polygon_index < b.polygon_index ?
          std::make_pair(a.polygon_index, b.polygon_index) :
          std::make_pair(b.polygon_index, a.polygon_index) );
    }
  }
};

} //end of namespace internal

template <class Kernel>
bool check_intersection_in_slice(
  const std::vector< boost::shared_ptr< std::vector< typename Kernel::Point_3 > > >& slice,
  int constant_coordinate,
  std::set<std::pair<std::size_t, std::size_t> >& intersecting_polygons)
{
  using namespace internal;

  typedef Box_with_segment_and_polygon_id<Kernel> Box;
  Report_inters<Kernel> report_inters(intersecting_polygons);

  std::vector<Box> boxes;
  std::size_t nb_contours=slice.size();
  for (std::size_t k=0;k<nb_contours; ++k){
    const std::vector< typename Kernel::Point_3 >& contour=*slice[k];
    std::size_t nb_edges=contour.size();
    typename Kernel::Point_2 prev(
      contour[0][(constant_coordinate+1)%3],
      contour[0][(constant_coordinate+2)%3]
    );
    for (std::size_t j=1;j<nb_edges;++j)
    {
      typename Kernel::Point_2 point(
        contour[j][(constant_coordinate+1)%3],
        contour[j][(constant_coordinate+2)%3]
      );
      boxes.push_back(Box(prev,point,k));
      prev=point;
    }
  }

  box_self_intersection_d( boxes.begin(), boxes.end(), report_inters);

  return intersecting_polygons.empty();
}

template <class Kernel, bool verbose=false>
class Contour_checker_and_fixer{
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Point_2 Point_2;
  typedef boost::shared_ptr< std::vector<typename Kernel::Point_3> > Point_container; /// \todo this should be a template
  typedef std::vector< boost::shared_ptr< std::vector<Point_3> > > Slice;
  enum State{SLICE_BEGAN, SLICE_ENDED};
///data members
  std::vector< Slice > m_slices;
  int m_constant_coordinate;
  //variables used for the input checking
  std::set< std::pair<std::size_t,std::size_t> > m_intersecting_polygons;
  bool m_last_contour_valid;
  State m_state;
  //variables used for the contour provider
  std::size_t m_current_slice;
  std::size_t m_current_polygon;
  std::size_t m_current_point;
public:
///Constructor
  Contour_checker_and_fixer()
    :m_constant_coordinate(-1)
    ,m_last_contour_valid(true)
    ,m_state(SLICE_ENDED)
    ,m_current_slice(0)
    ,m_current_polygon(0)
    ,m_current_point(0)
  {}

  Contour_checker_and_fixer(int constant_coordinate)
    :m_constant_coordinate(constant_coordinate)
    ,m_last_contour_valid(true)
    ,m_state(SLICE_ENDED)
    ,m_current_slice(0)
    ,m_current_polygon(0)
    ,m_current_point(0)
  {}

  bool begin_slice()
  {
    if ( m_state==SLICE_BEGAN ) return false;
    if ( !m_intersecting_polygons.empty() ) return false;
    if ( !m_last_contour_valid ) return false;
    m_slices.push_back( Slice() );
    m_state=SLICE_BEGAN;
    return true;
  }

  Error_code
  add_contour_to_slice(
    boost::shared_ptr< std::vector<typename Kernel::Point_3> > contour)
  {
    if (m_constant_coordinate==-1) return LOGIC_ERROR;
    if (!m_last_contour_valid) return LOGIC_ERROR;
    if (m_state==SLICE_ENDED) return LOGIC_ERROR;
    if ( contour->empty() ) return DEGENERATE_POLYGON;

    Error_code errcode =
      create_polygon<Kernel>( *contour, m_constant_coordinate,verbose );

    if (errcode==POLYGON_NOT_SIMPLE)
    {
      m_last_contour_valid=false;
      m_slices.back().push_back(contour);
    }
    else
      if (errcode==VALID_OR_TRIVIALLY_FIXED_POLYGON)
        m_slices.back().push_back(contour);
    return errcode;
  }

  //returns true if not empty and if the polygon does not self intersect
  bool end_slice()
  {
    if (m_state==SLICE_ENDED) return false;
    if ( m_slices.back().empty() ) return false;
    if (!m_last_contour_valid) return false;

    m_intersecting_polygons.clear();

    if (m_slices.back().size()<2)
      check_intersection_in_slice<Kernel>(m_slices.back(),
                                          m_constant_coordinate,
                                          m_intersecting_polygons);

    if (m_intersecting_polygons.empty())
    {
      m_state=SLICE_ENDED;
      return true;
    }

    return false;
  }

  //remove the last contour added
  void pop_slice_back(){
    m_slices.pop_back();
    m_intersecting_polygons.clear();
    m_last_contour_valid=true;
    m_state=SLICE_ENDED;
  }

  void pop_contour_back(){
    if (m_slices.empty()) return;
    m_slices.back().pop_back();
    m_last_contour_valid=true;
  }

  bool last_slice_is_valid() const { return m_state==SLICE_ENDED; }
  bool last_contour_is_valid() const{ return m_last_contour_valid; }

  int constant_coordinate() const{
    return m_constant_coordinate;
  }


  void set_constant_coordinate(int i)
  {
    m_constant_coordinate=i;
  }

//  void fix_last_contour();
//  void last_contour_is_valid();
//  void sort_slices();

  /// public members required to be a contour provider
  bool empty() const{
    return m_slices.empty();
  }

  #ifndef SWIG
  void next_polygon(){
    ++m_current_polygon;
    if ( m_current_polygon == m_slices[ m_current_slice ].size() )
    {
      ++m_current_slice;
      m_current_polygon=0;
    }
    m_current_point=0;
  }

  bool has_next_planar_contour() const {
    return m_slices.size()-1 != m_current_slice ||
           m_slices.back().size()-1!=m_current_polygon;
  }

  std::size_t number_of_points() const {
    return m_slices[m_current_slice][m_current_polygon]->size();
  }

  bool has_another_component() const {
    return m_slices[m_current_slice].size()-1 != m_current_polygon;
  }

  Point_2 get_point(){
    CGAL_assertion( m_constant_coordinate!=-1 );
    double x = m_slices[m_current_slice][m_current_polygon]->operator[](m_current_point)[(m_constant_coordinate+1)%3];
    double y = m_slices[m_current_slice][m_current_polygon]->operator[](m_current_point)[(m_constant_coordinate+2)%3];
    ++m_current_point;
    return Point_2(x,y);
  }

  Point_2 get_point(double& z){
    z=m_slices[m_current_slice][m_current_polygon]->operator[](m_current_point)[m_constant_coordinate];
    return get_point();
  }

  /// Debug functions
  void output_slices(std::ostream& output)
  {
    for (std::size_t s=0; s< m_slices.size(); ++s)
    {
      Slice& slice = m_slices[s];
      CGAL_assertion( !slice.empty() );
      for (std::size_t p=0; p<slice.size(); ++p)
      {
        output << slice[p]->size() << "\n";
        std::copy( slice[p]->begin(), slice[p]->end(),
                   std::ostream_iterator<Point_3>(output, "\n") );
      }
    }
  }
  #endif
};

} } //end of namespace CGAL::Reconstruction_from_parallel_slices

#endif // CGAL_RECONSTRUCTION_FROM_PARALLEL_SLICES_3_CHECK_AND_FIX_INPUT_H