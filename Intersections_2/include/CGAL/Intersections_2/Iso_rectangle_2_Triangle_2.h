// Copyright (c) 2002-2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Radu Ursu

#ifndef CGAL_INTERSECTIONS_2_ISO_RECTANGLE_2_TRIANGLE_2_H
#define CGAL_INTERSECTIONS_2_ISO_RECTANGLE_2_TRIANGLE_2_H

#include <CGAL/Triangle_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Intersections_2/Segment_2_Segment_2.h>
#include <CGAL/Intersection_traits_2.h>
#include <CGAL/Intersections_2/Iso_rectangle_2_Segment_2.h>

#include <vector>
#include <list>


namespace CGAL{

namespace Intersections {

namespace internal {

  template <class K>
  typename Intersection_traits<K, typename K::Triangle_2, typename K::Iso_rectangle_2>::result_type
  intersection(const typename K::Triangle_2 &t, const typename K::Iso_rectangle_2 &r, const K& kk)
  {
    typedef typename K::FT FT;
    typedef typename K::Segment_2 Segment;
    typedef typename K::Triangle_2 Triangle;
    typedef typename K::Point_2   Point;

    FT xr1, yr1, xr2, yr2;
    bool position[3][4] = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
    bool is_inside[3] = {false, false, false}; //true when a vertex is inside rectangle

    xr1 = r.xmin(); xr2 = r.xmax();
    yr1 = r.ymax(); yr2 = r.ymin();

    Point upper_left, lower_right;
    Point p[3]; //the vertices of triangle

    upper_left = Point(xr1, yr1); //upper left
    lower_right = Point(xr2, yr2); //lower right

    p[0] = t.vertex(0);
    p[1] = t.vertex(1);
    p[2] = t.vertex(2);

    //check the points of the triangle
    for(int i=0; i<3; i++){
      if(!(compare_x(p[i], upper_left) == SMALLER))
        if(!(compare_x(p[i], lower_right) == LARGER))
          if(!(compare_y(p[i], upper_left) == LARGER))
            if(!(compare_y(p[i], lower_right) == SMALLER))
              is_inside[i] = true; //the point is inside
            else
              position[i][2] = true; //South
          else
            position[i][0] = true; //North
        else
        {
          position[i][3] = true; //East
          if(compare_y(p[i], upper_left) == LARGER)
            position[i][0] = true; //North
          else if(compare_y(p[i], lower_right) == SMALLER)
            position[i][2] = true; //South
        }
      else
      {
        position[i][1] = true;  //West
        if(compare_y(p[i], upper_left) == LARGER)
          position[i][0] = true; //North
        else if(compare_y(p[i], lower_right) == SMALLER)
          position[i][2] = true; //South

      }
    }

    //test if the triangle is completely to the left, right, below or above the rectangle
    bool intersection = true; //true if it could be a intersection with the rectangle
    for(int j=0; j<4; j++)
      if(position[0][j] && position[1][j] && position[2][j] ){
        intersection = false;
        break;
      }

    if(intersection){
      Segment s[4]; //the segments that identify the N, W, S, E
      bool outside = false;
      bool status_in_list[3] = {false, false, false}; //true if the triangle's points are in the result vector
      std::list<int> last_intersected;
      int last_intersected_segment = 5; //could be 0=N, 1=W, 2=S, 3=E
      last_intersected.push_back(5);
      int status_intersected[4] = {0, 0, 0, 0}; //the number of intersections for each segment
      std::vector<Point> result; //the vector containing the result vertices
      int next; //the index of the next vertex

      s[0] = Segment(Point(xr2, yr1), Point(xr1, yr1)); //N
      s[1] = Segment(Point(xr1, yr1), Point(xr1, yr2)); //W
      s[2] = Segment(Point(xr1, yr2), Point(xr2, yr2)); //S
      s[3] = Segment(Point(xr2, yr2), Point(xr2, yr1)); //E

      //assign to p[i] the vertices of the triangle in ccw order
      if(t.orientation() == CGAL::CLOCKWISE)
      {
        p[0] = t.vertex(2);
        p[2] = t.vertex(0);

        is_inside[0] = is_inside[2] ^ is_inside[0];
        is_inside[2] = is_inside[2] ^ is_inside[0];
        is_inside[0] = is_inside[0] ^ is_inside[2];

        for(int i=0; i<4; i++){
          position[0][i] = position[2][i] ^ position[0][i];
          position[2][i] = position[2][i] ^ position[0][i];
          position[0][i] = position[0][i] ^ position[2][i];
        }
      }

      for(int index=0; index<3; index++) //for each vertex
      {
        next=(index+1)%3;
        if(is_inside[index]){ // true if first is inside
          if(!status_in_list[index]){  //if is not yet in the list
            result.push_back(p[index]);
            status_in_list[index] = true;
          }
          if(is_inside[next]){ //true if second is inside
            //add the points in the vector
            if(!status_in_list[next]){ // if is not yet in the list
              result.push_back(p[next]);
              status_in_list[next] = true;
            }
          } else { //I'm going out, the second is outside
            for(int j=0; j<4; j++) // for all directions
              if(position[next][j]) // if it's a second point direction
              {
                //test for intersection
                typename Intersection_traits<K, Segment, Segment>::result_type
                  v = internal::intersection(Segment(p[index], p[next]), s[j], kk);
                if(v) {
                  if(const Point *p_obj = intersect_get<Point>(v))
                {
                  //intersection found
                  outside = true;
                  result.push_back(*p_obj); //add the intersection point
                  if(last_intersected.back()!=j)
                    last_intersected.push_back(j);
                  status_intersected[j]++;
                }
              }
          }
          }
        } else { //the first point is outside
          for(int j=0; j<4; j++) // for all directions
            if(position[index][j]) //watch only the first point directions
            {
              //test for intersection
              typename Intersection_traits<K, Segment, Segment>::result_type
                v = internal::intersection(Segment(p[index], p[next]), s[j], kk);
              if(v) {
                if(const Point *p_obj = intersect_get<Point>(v))
              {
                //intersection found
                outside = false;
                last_intersected_segment = last_intersected.back();
                if(last_intersected_segment != 5 && last_intersected_segment != j && status_intersected[j] == 0){
                  //add the target of each rectangle segment in the list
                  if(last_intersected_segment < j)
                    while(last_intersected_segment < j){
                      result.push_back(s[last_intersected_segment].target());
                      last_intersected_segment++;
                    }
                  else{
                    while(last_intersected_segment < 4){
                      result.push_back(s[last_intersected_segment].target());
                      last_intersected_segment++;
                    }
                    last_intersected_segment = 0;
                    while(last_intersected_segment < j){
                      result.push_back(s[last_intersected_segment].target());
                      last_intersected_segment++;
                    }
                  }
                }
                result.push_back(*p_obj); //add the intersection point in the list
                if(last_intersected.back()!=j)
                  last_intersected.push_back(j);
                status_intersected[j]++;
                if(!is_inside[next]){ //if the second point is not inside search a second intersection point
                  for(j=0; j<4; j++) // for all directions
                    if(position[next][j])
                    {
                      //test for intersection
                        typename Intersection_traits<K, Segment, Segment>
                          ::result_type
                          v = internal::intersection(Segment(p[index], p[next]), s[j], kk);
                        if(v) {
                          if(const Point *p_obj = intersect_get<Point>(v))
                           //found the second intersection
                      {
                        outside = true;
                        result.push_back(*p_obj);
                        if(last_intersected.back()!=j)
                          last_intersected.push_back(j);
                        status_intersected[j]++;
                      }
                    }
                      }
                }//end if the second point is not inside
              }
              } // end v
            }
        }//end else (the first point is outside)
      }//endfor
      if(outside){
        std::list<int>::const_iterator it = last_intersected.begin();
        while(*it == 5)
          it++;
        last_intersected_segment = *it;
        int j = last_intersected.back();
        if(last_intersected_segment != 5 && last_intersected_segment != j){
          //add the target of each rectangle segment in the list
          if(last_intersected_segment > j)
            while(last_intersected_segment > j){
              result.push_back(s[j].target());
              j++;
            }
          else{
            while(j<4){
              result.push_back(s[j].target());
              j++;
            }
            j = 0;
            while(j<last_intersected_segment){
              result.push_back(s[j].target());
              j++;
            }
          }
        }
      }//end if(outside)

      //test if were not intersections
      //between the triangle and the rectangle
      if(status_intersected[0] == 0 && status_intersected[1] == 0 &&
        status_intersected[2] == 0 && status_intersected[3] == 0)
      {
        //should test if the rectangle is inside the triangle
        if(t.bounded_side(upper_left) == CGAL::ON_BOUNDED_SIDE){
          for(int k=0; k<4; k++)
            result.push_back(s[k].source());
        }
      }
      //remove duplicated consecutive points
      typename std::vector<Point>::iterator last = std::unique(result.begin(),result.end());
      result.erase(last,result.end());

      while(result.size() > 1 && result.back() == result.front())
        result.pop_back();

      switch(result.size()){
        case 0:
          return intersection_return<typename K::Intersect_2, Triangle, typename K::Iso_rectangle_2>();
        case 1:
          return intersection_return<typename K::Intersect_2, Triangle, typename K::Iso_rectangle_2>(result[0]);
        case 2:
          return intersection_return<typename K::Intersect_2, Triangle, typename K::Iso_rectangle_2>(Segment(result[0], result[1]));
        case 3:
          return intersection_return<typename K::Intersect_2, Triangle, typename K::Iso_rectangle_2>(Triangle(result[0], result[1], result[2]));
        default:
          return intersection_return<typename K::Intersect_2, Triangle, typename K::Iso_rectangle_2>(result);
      }

    }//end if(intersection)
    return intersection_return<typename K::Intersect_2, Triangle, typename K::Iso_rectangle_2>();
  }//end intersection

  template <class K>
  typename Intersection_traits<K, typename K::Triangle_2, typename K::Iso_rectangle_2>::result_type
  inline intersection(const typename K::Iso_rectangle_2 &r, const typename K::Triangle_2 &t, const K& k)
  {
    return intersection(t,r,k);
  }

  template <class K>
  bool do_intersect(
    const typename K::Triangle_2 &tr,
    const typename K::Iso_rectangle_2 &ir,
    const K& k)
  {
    //1) check if at least one vertex of tr is not outside ir
    //2) if not, check if at least on vertex of tr is not outside tr

    typename K::Has_on_unbounded_side_2 unbounded_side=k.has_on_unbounded_side_2_object();
    typename K::Construct_vertex_2 vertex=k.construct_vertex_2_object();

    for (int i=0;i<3;++i)
      if ( !unbounded_side( ir,vertex(tr,i) ) ) return true;
    for (int i=0;i<4;++i)
      if ( !unbounded_side( tr,vertex(ir,i) ) ) return true;

    typename K::Construct_segment_2 segment=k.construct_segment_2_object();
    for (int i=0;i<3;++i)
      if ( do_intersect(
            segment(vertex(tr,i),vertex(tr,(i+1)%3)),
            ir, k) ) return true;

    return false;
  }

  template <class K>
  inline bool do_intersect(
    const typename K::Iso_rectangle_2 &ir,
    const typename K::Triangle_2 &tr,
    const K& k)
  {
    return do_intersect(tr,ir,k);
  }

} //namespace internal
} // namespace Intersections

CGAL_INTERSECTION_FUNCTION(Triangle_2, Iso_rectangle_2, 2)
CGAL_DO_INTERSECT_FUNCTION(Triangle_2, Iso_rectangle_2, 2)

}//end namespace

#endif
