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
#include <CGAL/property_map.h>
#include <boost/iterator/counting_iterator.hpp>
#include <CGAL/Kernel_traits.h>

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

namespace internal{
  template <class Ppmap>
  struct Identical_points{
    typedef typename boost::property_traits<Ppmap>::key_type Key;

    Ppmap ppmap;
    Key ref;
    Identical_points(Ppmap ppmap_, Key ref_)
    : ppmap(ppmap_)
    , ref(ref_)
    {}
    bool operator()(const Key& k) const
    {
      return get(ppmap,k) == get(ppmap,ref);
    }
  };
} //end of namespace internal


template <class KeyType, class PointPmap>
bool
find_safe_start_to_fix_consecutive_overlapping_segments(
  std::list< KeyType >& points, PointPmap ppmap)
{
  typedef typename std::list< KeyType >::iterator iterator;

  iterator pit = points.begin();
  iterator qit = cpp11::next(pit);
  iterator rit = cpp11::next(qit);
  std::size_t count = points.size();

  points.pop_back();
  internal::Identical_points<PointPmap> same_pt(ppmap, *pit);
  while(std::find_if(qit, points.end(),same_pt) != points.end()){
   //change the start point (degree > 2)
    points.splice(points.end(), points, pit);
    if(--count == 0){
      points.push_back(points.front());
      return false;
    }
    pit = qit;
    ++qit;
    ++rit;
  }

  while(collinear( get(ppmap, *pit),
                   get(ppmap, *qit),
                   get(ppmap, *rit) ) &&
        collinear_are_strictly_ordered_along_line(get(ppmap, *pit),
                                                  get(ppmap, *qit),
                                                  get(ppmap, *rit))
  ){
    // change the start point (collinear overlapping)
    points.splice(points.end(), points, pit);
    if(--count == 0){
      points.push_back(points.front());
      return false;
    }
    pit = qit;
    ++qit;
    ++rit;
  }
  points.push_back(points.front());
  return true;
}

/// Given a list of points representing a polygon, remove whether
/// consecutive collinear segments overlapping
/// find_safe_start_to_fix_consecutive_overlapping_segments is supposed to be
/// called prior to this function to prepare the data.
/// ex: if my segments are ab and bc, then ab is kept
///  --a----c----b--  --> --a------b--
/// \return the number of points deleted
template <class KeyType, class PointPmap>
std::size_t fix_consecutive_overlapping_segments(
  std::list< KeyType >& points, PointPmap ppmap)
{
  typedef typename boost::property_traits<PointPmap>::reference Point_2_ref;
  typedef typename std::list< KeyType >::iterator iterator;

  std::size_t nb_points_erased=0;

  iterator pit = points.begin();
  while(true){
    if(points.size() == 4){
      break;
    }
    iterator qit = pit;
    ++qit;
    if(qit == points.end()){
      break;
    }
    Point_2_ref p = get(ppmap, *pit);
    Point_2_ref q = get(ppmap, *qit);

    iterator rit = qit;
    ++rit;
    if(rit == points.end()){
      Point_2_ref r = get(ppmap, *(++points.begin()));
      if(CGAL::collinear(p,q,r) &&
         (! CGAL::collinear_are_strictly_ordered_along_line(p,q,r))){
        points.pop_back();
        points.pop_front();
        points.push_back(points.front()); // duplicate the first point
      }
      break;
    }
    Point_2_ref r = get(ppmap, *rit);

    if(CGAL::collinear(p,q,r)){
      if(! CGAL::collinear_are_strictly_ordered_along_line(p,q,r)){
        points.erase(qit);
        ++nb_points_erased;

        if(get(ppmap, *pit) == get(ppmap, *rit)){
          // identical consecutive points
          points.erase(rit);
          ++nb_points_erased;

        }
        if(pit != points.begin()){
            --pit;
          }
        continue;  // and do not advance pit, because its two successors are different
      }
    }
    ++pit;
  }
  return nb_points_erased;
}

namespace internal{
  template <class Kernel>
  struct Get_2D_point_from_3D_point_id{
    typedef std::size_t key_type;
    typedef typename Kernel::Point_2 value_type;
    typedef value_type reference;
    typedef boost::readable_property_map_tag category;

    const std::vector< typename Kernel::Point_3 >& points_3d;
    int constant_coordinate;

    Get_2D_point_from_3D_point_id(
      const std::vector< typename Kernel::Point_3 >& points_3d_,
      int constant_coordinate_)
    : points_3d(points_3d_)
    , constant_coordinate(constant_coordinate_)
    {}

    friend reference get(Get_2D_point_from_3D_point_id map, key_type k) {
      return value_type( map.points_3d[k][(map.constant_coordinate+1)%3],
                         map.points_3d[k][(map.constant_coordinate+2)%3] );
    }
  };
}

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
  /// \todo if a polygon contains only collinear points handle it as degenerate
  ///       but it is only while the code does not handle it correctly because
  ///       it should be allowed as well as isolated points

  std::list<std::size_t> point_indices(
    boost::counting_iterator<std::size_t>(0),
    boost::counting_iterator<std::size_t>(points_3d.size()) );
  internal::Get_2D_point_from_3D_point_id<Kernel>
    ppmap(points_3d, constant_coordinate);
  if( find_safe_start_to_fix_consecutive_overlapping_segments
        (point_indices, ppmap) )
  {
    std::size_t nb_pt_erased =
      fix_consecutive_overlapping_segments(point_indices, ppmap);
    if (nb_pt_erased!=0)
    {
      if (verbose)
        std::cerr << "Warning: removed some overlapping segments...\n";
      //we need to rewrite the vector of points
      std::vector<typename Kernel::Point_3> tmp;
      tmp.reserve( points_3d.size()-nb_pt_erased );
      for (std::list<std::size_t>::iterator it=point_indices.begin(),
            end=point_indices.end(); it!=end; ++it)
      {
        tmp.push_back( points_3d[*it] );
      }
      points_3d.swap(tmp);
    }
  }
  else
    //only collinear points with overlapping segments have been found
    return DEGENERATE_POLYGON;

  //check that the polygon is of dimension 2
  nb_pts=points_3d.size();
  if (nb_pts < 4) return DEGENERATE_POLYGON;
  for (std::size_t i=0, j=1, k=2;;)
  {
    if ( CGAL::coplanar_orientation(points_3d[i],points_3d[j],points_3d[k]) !=
         CGAL::COLLINEAR ) break;
    if (++k==nb_pts) return DEGENERATE_POLYGON;
  }

  //now create the polygon
  Polygon_2<Kernel> polygon;

  for (std::size_t i=0, end=points_3d.size()-1;i!=end;++i)
  {
    polygon.push_back(get(ppmap, i));
  }

  if ( !polygon.is_simple() )
  {
    return POLYGON_NOT_SIMPLE;
  }

  return VALID_OR_TRIVIALLY_FIXED_POLYGON;
}


namespace internal {
template <class Kernel>
struct Less {
  bool operator()(const std::pair<typename Kernel::Direction_2, int> d1,
                  const std::pair<typename Kernel::Direction_2, int> d2) const
  {
    return CGAL::compare_angle_with_x_axis(d1.first,d2.first) == CGAL::SMALLER;
  }
};

template <class Kernel>
int
crossing( const typename Kernel::Point_2& p,
          const typename Kernel::Point_2& q,
          const typename Kernel::Point_2& q2,
          const typename Kernel::Point_2& r,
          const typename Kernel::Point_2& s2)
{
  //std::cerr << "crossing test:\n" << p << std::endl;
  //std::cerr << q2 << std::endl;
  //std::cerr << r << std::endl;
  //std::cerr << s2 << std::endl;
  typedef typename Kernel::Direction_2 Direction_2;

  std::vector<std::pair<Direction_2, int> > directions;
  directions.push_back(std::make_pair(Direction_2(p-q),0));
  directions.push_back(std::make_pair(Direction_2(q2-q),0));
  directions.push_back(std::make_pair(Direction_2(r-q),1));
  directions.push_back(std::make_pair(Direction_2(s2-q),1));

  Less<Kernel> less;
  sort(directions.begin(), directions.end(),less);
  if((directions[0].first == directions[1].first) ||
     (directions[1].first == directions[2].first) ||
     (directions[2].first == directions[3].first)){
    return 0; // collinear segments
  }
  if((directions[0].second == directions[1].second) ||
     (directions[1].second == directions[2].second) ||
     (directions[2].second == directions[3].second)){
    return -1; // no crossing
  }
  return 1; // a crossing
}

} //end of namespace internal

template <class KeyType, class PointPmap>
bool fix_self_intersections(std::list< KeyType >& points, PointPmap ppmap)
{
  typedef typename boost::property_traits<PointPmap>::reference Point_2_ref;
  typedef typename boost::property_traits<PointPmap>::value_type Point_2;
  typedef typename std::list< KeyType >::iterator iterator;
  typedef typename CGAL::Kernel_traits<Point_2>::Kernel Kernel;
  typedef typename Kernel::Segment_2 Segment_2;

  iterator pit = points.begin();
  iterator end = points.end();

  while(true){
    bool swapped = false;
    const Point_2& p = get(ppmap, *pit);
    iterator qit = pit;
    ++qit;
    if(qit == end){
      break;
    }

    const Point_2& q = get(ppmap, *qit);
    //std::cerr << "\nouter loop: " << p << "  " << q << std::endl;
    iterator rit = qit;
    ++rit;
    while(true){
      if(rit == end){
        break;
      }
      const Point_2& r = get(ppmap, *rit);
      iterator sit = rit;
      ++sit;
      if(sit == end){
        break;
      }

      const Point_2& s = get(ppmap, *sit);
      //std::cerr << "  inner loop: " << r << "  " << s << std::endl;
      Segment_2 segA(p,q), segB(r,s);
      if(CGAL::do_intersect(segA,segB)){
        //std::cerr << "intersection" << std::endl;
        if(p!=r && p!=s && q!=r && q!=s){
          //std::cerr << "real intersection" << std::endl;
          //std::cerr << segA << std::endl;
          //std::cerr << segB << std::endl;
          std::list<KeyType> tmp;
          tmp.splice(tmp.begin(),
                     points,
                     qit,
                     sit);
          tmp.reverse();
          points.splice(sit,tmp);
          swapped = true;
          break;
        }
        else if (q == s){
          iterator q2it = qit; ++q2it;
          const Point_2& q2 = get(ppmap, *q2it);
          iterator s2it = sit; ++s2it;
          const Point_2& s2 = get(ppmap, *s2it);
          int sign = internal::crossing<Kernel>(p,q,q2,r,s2);
          if(sign == 1){
            //std::cerr << "crossing--------------------" << std::endl;
            std::list<KeyType> tmp;
            tmp.splice(tmp.begin(),
                       points,
                       q2it,
                       sit);
            tmp.reverse();
            points.splice(s2it,tmp);
            if((! CGAL::do_intersect(Segment_2(p,r),
                                     Segment_2(s,s2))) &&
               (! CGAL::do_intersect(Segment_2(p,r),
                                     Segment_2(q,q2)))){
              points.erase(qit);
            } else if((! CGAL::do_intersect(Segment_2(q2,s2),
                                           Segment_2(p,q))) &&
                      (! CGAL::do_intersect(Segment_2(q2,s2),
                                            Segment_2(q,r)))
                      ) {
              points.erase(sit);

            } else {
              //std::cerr << "shortcut of pqq2 or rss2 would introduce an intersection" << std::endl;
              return false;
            }
            points.splice(qit,tmp);
          } else if(sign != 1){
            //std::cerr << "touching--------------------" << std::endl;
            //std::cerr << "2\n" << p << std::endl << q2 << std::endl;
            //std::cerr << "2\n" << r << std::endl << s2 << std::endl;
            if((! CGAL::do_intersect(Segment_2(p,q2),
                                     Segment_2(r,s))) &&
               (! CGAL::do_intersect(Segment_2(p,q2),
                                     Segment_2(s2,s)))){
              points.erase(qit);
            } else if((! CGAL::do_intersect(Segment_2(r,s2),
                                           Segment_2(p,q))) &&
                      (! CGAL::do_intersect(Segment_2(r,s2),
                                            Segment_2(q2,q)))
                      ) {
              points.erase(sit);

            } else {
              return false;
              //std::cerr << "shortcut of pqq2 or rss2 would introduce an intersection" << std::endl;
            }
          } else {
            //std::cerr << "todo: treat overlapping segments" << std::endl;
            return false;
          }
        }
      }
      ++rit;
        }
    if(! swapped){
      ++pit;
    }
  }
  return true;
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

    /// \todo allow polygons to share a vertex?
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
  typedef boost::shared_ptr< std::vector<Point_3> > Point_container; /// \todo this should be a template
  typedef std::vector< Point_container > Slice;
  enum State{SLICE_BEGUN, SLICE_ENDED};
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
    if ( m_state==SLICE_BEGUN ) return false;
    if ( !m_intersecting_polygons.empty() ) return false;
    if ( !m_last_contour_valid ) return false;
    m_slices.push_back( Slice() );
    m_state=SLICE_BEGUN;
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

    if (m_slices.back().size()>1)
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

  bool fix_self_intersections(){
    if (m_last_contour_valid) return true;

    std::vector<typename Kernel::Point_3>& points_3d=*m_slices.back().back();
    std::list<std::size_t> point_indices(
      boost::counting_iterator<std::size_t>(0),
      boost::counting_iterator<std::size_t>(points_3d.size()) );

    internal::Get_2D_point_from_3D_point_id<Kernel>
      ppmap(points_3d, m_constant_coordinate);

    if ( ::CGAL::Reconstruction_from_parallel_slices::
          fix_self_intersections(point_indices, ppmap) )
    {
      //we need to rewrite the vector of points
      std::vector<typename Kernel::Point_3> tmp;
      tmp.reserve( points_3d.size() ); // the bound is not tight
      for (std::list<std::size_t>::iterator it=point_indices.begin(),
            end=point_indices.end(); it!=end; ++it)
      {
        tmp.push_back( points_3d[*it] );
      }
      points_3d.swap(tmp);

      //make sure the polygon is now simple!
      Polygon_2<Kernel> polygon;

      for (std::size_t i=0, end=points_3d.size()-1;i!=end;++i)
      {
        polygon.push_back(get(ppmap, i));
      }

      if ( !polygon.is_simple() ) return false;

      m_last_contour_valid=true;
      return true;
    }
    else
      return false;
  }

  std::size_t size_of_contours() const {
    return m_slices.empty() ? 0 : m_slices.back().size();
  }

  boost::shared_ptr< std::vector<typename Kernel::Point_3> >
  contour(std::size_t i){
    if (m_slices.empty() || m_slices.back().size() <= i )
      return boost::shared_ptr< std::vector<typename Kernel::Point_3> >(
          new std::vector<typename Kernel::Point_3>()
        );
    return m_slices.back()[i];
  }

  boost::shared_ptr< std::vector<typename Kernel::Point_3> >
  contours_back()
  {
    if ( m_slices.empty() || m_slices.back().empty() )
      return boost::shared_ptr< std::vector<typename Kernel::Point_3> >(
          new std::vector<typename Kernel::Point_3>()
        );
    return m_slices.back().back();
  }

  std::set< std::pair<std::size_t,std::size_t> >
  intersection_contours(){
    return m_intersecting_polygons;
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