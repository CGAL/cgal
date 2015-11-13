// Copyright (c) 2015 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : 
//

#ifndef CGAL_STRUCTURE_POINT_SET_3_H
#define CGAL_STRUCTURE_POINT_SET_3_H

#include <CGAL/property_map.h>
#include <CGAL/point_set_processing_assertions.h>
#include <CGAL/assertions.h>

#include <CGAL/centroid.h>

#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Fuzzy_iso_box.h>
#include <CGAL/Search_traits_d.h>

#include <iterator>
#include <list>


namespace CGAL {


// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
/// \cond SKIP_IN_MANUAL
namespace internal {

  template <typename Traits>
  class Point_set_structuring
  {
  public:

    typedef Point_set_structuring<Traits> Self;

    typedef typename Traits::FT FT;
    typedef typename Traits::Point_3 Point;
    typedef typename Traits::Vector_3 Vector;
    typedef typename Traits::Segment_3 Segment;
    typedef typename Traits::Line_3 Line;

    typedef typename Traits::Plane_3 Plane;

    typedef typename Traits::Point_map Point_map;
    typedef typename Traits::Normal_map Normal_map;
    typedef typename Traits::Input_range Input_range;

    typedef typename Input_range::iterator Input_iterator;

    typedef Shape_detection_3::Shape_base<Traits> Shape; 
    typedef Shape_detection_3::Plane<Traits> Plane_shape;


  private:

    const std::size_t minus1;
    
    class My_point_property_map{
      const std::vector<Point>& points;
    public:
      typedef Point value_type;
      typedef const value_type& reference;
      typedef std::size_t key_type;
      typedef boost::lvalue_property_map_tag category;  
      My_point_property_map (const std::vector<Point>& pts) : points (pts) {}
      reference operator[] (key_type k) const { return points[k]; }
      friend inline reference get (const My_point_property_map& ppmap, key_type i) 
      { return ppmap[i]; }
    };

    enum Point_status { POINT, EDGE, CORNER, SKIPPED };
    
    struct Edge
    {
      CGAL::cpp11::array<std::size_t, 2> planes;
      std::vector<std::size_t> indices; // Points belonging to intersection
      Line support;
      bool active;

      Edge (std::size_t a, std::size_t b)
      { planes[0] = a; planes[1] = b; active = true; }
    };
    struct Corner
    {
      std::vector<std::size_t> planes;
      std::vector<std::size_t> edges;
      Point support;
      bool active;

      Corner (std::size_t p1, std::size_t p2, std::size_t p3,
              std::size_t e1, std::size_t e2, std::size_t e3)
      {
        planes.resize (3); planes[0] = p1; planes[1] = p2; planes[2] = p3;
        edges.resize (3); edges[0] = e1; edges[1] = e2; edges[2] = e3;
        active = true;
      }
    };
      

    Traits m_traits;

    std::vector<Point> m_points;
    std::vector<std::size_t> m_indices;
    std::vector<Point_status> m_status;
    Point_map m_point_pmap;
    Normal_map m_normal_pmap;
    
    std::vector<boost::shared_ptr<Plane_shape> > m_planes;
    std::vector<Edge> m_edges;
    std::vector<Corner> m_corners;
    
  public:

    Point_set_structuring (Traits t = Traits ())
      : minus1 (static_cast<std::size_t>(-1)), m_traits (t)
    {
    }

    
    Point_set_structuring (Input_iterator begin, Input_iterator end,
                           const Shape_detection_3::Efficient_RANSAC<Traits>& shape_detection)
      : minus1 (static_cast<std::size_t>(-1)), m_traits (shape_detection.traits())
    {
      for (Input_iterator it = begin; it != end; ++ it)
        m_points.push_back (get(m_point_pmap, *it));

      m_indices = std::vector<std::size_t> (m_points.size (), minus1);
      m_status = std::vector<Point_status> (m_points.size (), POINT);

      BOOST_FOREACH (boost::shared_ptr<Shape> shape, shape_detection.shapes())
        {
          boost::shared_ptr<Plane_shape> pshape
            = boost::dynamic_pointer_cast<Plane_shape>(shape);
        
          // Ignore all shapes other than plane
          if (pshape == boost::shared_ptr<Plane_shape>())
            continue;
          m_planes.push_back (pshape);

          for (std::size_t i = 0; i < pshape->indices_of_assigned_points().size (); ++ i)
            m_indices[pshape->indices_of_assigned_points()[i]] = m_planes.size () - 1;
        }

    }

    
    virtual ~Point_set_structuring ()
    {
      clear ();
    }

    void clear ()
    {

    }

    void run (double epsilon, double attraction_factor = 3.)
    {


      double radius = epsilon * attraction_factor;
      
      std::cerr << "Finding adjacent primitives... " << std::endl;
      find_pairs_of_adjacent_primitives (radius);
      std::cerr << " -> Found " << m_edges.size () << " pair(s) of adjacent primitives." << std::endl;

      std::cerr << "Computing edges... " << std::endl;
      compute_edges (epsilon);
      std::cerr << " -> Done" << std::endl;

      std::cerr << "Creating edge-anchor points... " << std::endl;
      {
        std::size_t size_before = m_points.size ();
        create_edge_anchor_points (radius, epsilon);
        std::cerr << " -> " << m_points.size () - size_before << " anchor point(s) created." << std::endl;
      }

      std::cerr << "Computating first set of corners... " << std::endl;
      compute_corners (radius);
      std::cerr << " -> Found " << m_corners.size () << " triple(s) of adjacent primitives/edges." << std::endl;

      std::cerr << "Merging corners... " << std::endl;
      {
        std::size_t size_before = m_points.size ();
        merge_corners (radius);
        std::cerr << " -> " << m_points.size () - size_before << " corner point(s) created." << std::endl;
      }

    }

  private:

    void find_pairs_of_adjacent_primitives (double radius)
    {
      typedef typename Traits::Search_traits Search_traits_base;
      typedef Search_traits_adapter <std::size_t, My_point_property_map, Search_traits_base> Search_traits;
      typedef CGAL::Kd_tree<Search_traits> Tree;
      typedef CGAL::Fuzzy_sphere<Search_traits> Fuzzy_sphere;

      My_point_property_map pmap (m_points);

      Tree tree (boost::counting_iterator<std::size_t> (0),
                 boost::counting_iterator<std::size_t> (m_points.size()),
                 typename Tree::Splitter(),
                 Search_traits (pmap));

      std::vector<std::vector<bool> > adjacency_table (m_planes.size (),
                                                       std::vector<bool> (m_planes.size (), false));

      //compute a basic adjacency relation (two primitives are neighbors
      //if at least one point of the primitive 1 is a k-nearest neighbor
      //of a point of the primitive 2 and vice versa)
      for (std::size_t i = 0; i < m_points.size (); ++ i)
        {
          std::size_t ind_i = m_indices[i];

          if (ind_i == minus1)
            continue;

          Fuzzy_sphere query (i, radius, 0., tree.traits());
          
          std::vector<std::size_t> neighbors;
          tree.search (std::back_inserter (neighbors), query); // WIP: SegFaults so far...

          
          for (std::size_t k = 0; k < neighbors.size(); ++ k)
            {
              std::size_t ind_k = m_indices[neighbors[k]];
              if (ind_k != minus1 && ind_k != ind_i)
                adjacency_table[ind_i][ind_k] = true;
            }
        }

      //verify the symmetry and store the pairs of primitives in
      //m_edges
      for (std::size_t i = 0; i < adjacency_table.size() - 1; ++ i)
        for (std::size_t j = i + 1; j < adjacency_table[i].size(); ++ j)
          if ((adjacency_table[i][j]) && (adjacency_table[j][i]))
            m_edges.push_back (Edge (i, j));

    }

    void compute_edges (double epsilon)
    {
      for (std::size_t i = 0; i < m_edges.size(); ++ i)
        {
          boost::shared_ptr<Plane_shape> plane1 = m_planes[m_edges[i].planes[0]];
          boost::shared_ptr<Plane_shape> plane2 = m_planes[m_edges[i].planes[1]];       

          double angle_A = std::acos (std::abs (plane1->plane_normal() * plane2->plane_normal()));
          double angle_B = CGAL_PI - angle_A;

          CGAL::Object ob_temp = CGAL::intersection (static_cast<Plane>(*plane1),
                                                     static_cast<Plane>(*plane2));
          if (!assign (m_edges[i].support, ob_temp))
            {
              std::cerr << "Warning: bad plane/plane intersection" << std::endl;
              continue;
            }

          Vector direction_p1 (0., 0., 0.);
          for (std::size_t k = 0; k < plane1->indices_of_assigned_points ().size(); ++ k)
            {
              std::size_t index_point = plane1->indices_of_assigned_points ()[k];
              
              const Point& point = m_points[index_point];
              Point projected = m_edges[i].support.projection (point);
              if (std::sqrt (CGAL::squared_distance (point, projected))
                  < 2 * std::min (4., 1 / std::sin (angle_A)) * epsilon
                  && m_status[index_point] != SKIPPED)
                direction_p1 = direction_p1 + Vector (projected, point);
            }
          if (direction_p1.squared_length() > 0)
            direction_p1 = direction_p1 / std::sqrt (direction_p1 * direction_p1);

          Vector direction_p2 (0., 0., 0.);
          for (std::size_t k = 0; k < plane2->indices_of_assigned_points ().size(); ++ k)
            {
              std::size_t index_point = plane2->indices_of_assigned_points ()[k];
              
              const Point& point = m_points[index_point];
              Point projected = m_edges[i].support.projection (point);
              if (std::sqrt (CGAL::squared_distance (point, projected))
                  < 2 * std::min (4., 1 / std::sin (angle_A)) * epsilon
                  && m_status[index_point] != SKIPPED)
                direction_p2 = direction_p2 + Vector (projected, point);
            }
          if (direction_p2.squared_length() > 0)
            direction_p2 = direction_p2 / std::sqrt (direction_p2 * direction_p2);

          double angle = std::acos (direction_p1 * direction_p2);
      
          if (direction_p1.squared_length() == 0
              || direction_p2.squared_length() == 0
              || (std::fabs (angle - angle_A) > 1e-2
                  && std::fabs (angle - angle_B) > 1e-2 ))
            {
              m_edges[i].active = false;
            }
        }
    }

    void create_edge_anchor_points (double radius, double epsilon)
    {
      double d_DeltaEdge = std::sqrt (2.) * epsilon;
      double r_edge = d_DeltaEdge / 2.;
      
      for (std::size_t i = 0; i < m_edges.size(); ++ i)
        {
          boost::shared_ptr<Plane_shape> plane1 = m_planes[m_edges[i].planes[0]];
          boost::shared_ptr<Plane_shape> plane2 = m_planes[m_edges[i].planes[1]];       

          const Line& line = m_edges[i].support;

          if (!(m_edges[i].active))
            {
              continue;
            }
							
          //find set of points close (<attraction_radius) to the edge and store in intersection_points
          std::vector<std::size_t> intersection_points;
          for (std::size_t k = 0; k < plane1->indices_of_assigned_points().size(); ++ k)
            {
              std::size_t index_point = plane1->indices_of_assigned_points()[k];
              Point point = m_points[index_point];
              Point projected = line.projection (point);
              if (CGAL::squared_distance (point, projected) < radius * radius)
                intersection_points.push_back (index_point);
            }
          for (std::size_t k = 0; k < plane2->indices_of_assigned_points().size(); ++ k)
            {
              std::size_t index_point = plane2->indices_of_assigned_points()[k];
              Point point = m_points[index_point];
              Point projected = line.projection (point);
              if (CGAL::squared_distance (point, projected) < radius * radius)
                intersection_points.push_back (index_point);
            }

          if (intersection_points.empty ())
            {
              continue;
            }

          const Point& t0 = m_points[intersection_points[0]];
          Point t0p = line.projection (t0);
          double dmin = 0.;
          double dmax = 0.;
          Point Pmin = t0p;
          Point Pmax = t0p;
          Vector dir = line.to_vector ();
          
          //compute the segment of the edge
          for (std::size_t k = 0; k < intersection_points.size(); ++ k)
            {
              std::size_t ind = intersection_points[k];
              const Point& point = m_points[ind];
              Point projected = line.projection (point);
              double d = Vector (t0p, projected) * dir;
                  
              if (d < dmin)
                {
                  dmin = d;
                  Pmin = projected;
                }
              else if (d > dmax)
                {
                  dmax = d;
                  Pmax = projected;
                }
            }

          //faire un partitionnement ds une image 1D en votant si
          //a la fois au moins un point de plan1 et aussi de plan
          //2 tombent dans une case (meme pas que pour les plans).
          Segment seg (Pmin,Pmax);
          int number_of_division = std::sqrt (seg.squared_length ()) / d_DeltaEdge + 1;
          std::vector<std::vector<std::size_t> > division_tab (number_of_division);

          for (std::size_t k = 0; k < intersection_points.size(); ++ k)
            {
              std::size_t ind = intersection_points[k];
              const Point& point = m_points[ind];
              Point projected = line.projection (point);

              std::size_t tab_index = std::sqrt (CGAL::squared_distance (seg[0], projected)) / d_DeltaEdge;

              division_tab[tab_index].push_back (ind);
            }

          //C1-CREATE the EDGE
          std::vector<int> index_of_edge_points;
          for (std::size_t j = 0; j < division_tab.size(); ++ j)
            {
              bool p1found = false, p2found = false;
              for (std::size_t k = 0; k < division_tab[j].size () && !(p1found && p2found); ++ k)
                {
                  if (m_indices[division_tab[j][k]] == m_edges[i].planes[0])
                    p1found = true;
                  if (m_indices[division_tab[j][k]] == m_edges[i].planes[1])
                    p2found = true;
                }

              if (!(p1found && p2found))
                {
                  division_tab[j].clear();
                  continue;
                }
              
              Point perfect (seg[0].x() + (seg[1].x() - seg[0].x()) * (j + 0.5) / (double)number_of_division,
                             seg[0].y() + (seg[1].y() - seg[0].y()) * (j + 0.5) / (double)number_of_division,
                             seg[0].z() + (seg[1].z() - seg[0].z()) * (j + 0.5) / (double)number_of_division);

              // keep closest point, replace it by perfect one and skip the others
              double dist_min = (std::numeric_limits<double>::max)();
              std::size_t index_best = 0;

              for (std::size_t k = 0; k < division_tab[j].size(); ++ k)
                {
                  std::size_t inde = division_tab[j][k];
                  m_status[inde] = SKIPPED; // Deactive all points except best (below)
                  double distance = CGAL::squared_distance (perfect, m_points[inde]);
                  if (distance < dist_min)
                    {
                      dist_min = distance;
                      index_best = inde;
                    }
                }

              m_points[index_best] = perfect;
              m_status[index_best] = EDGE;
              m_edges[i].indices.push_back (index_best);

            }

          //C2-CREATE the ANCHOR
          Vector direction_p1(0,0,0);
          Vector direction_p2(0,0,0);

          for (std::size_t j = 0; j < division_tab.size() - 1; ++ j)
            {
              if (division_tab[j].empty () || division_tab[j+1].empty ())
                continue;
              Point anchor (seg[0].x() + (seg[1].x() - seg[0].x()) * (j + 1) / (double)number_of_division,
                            seg[0].y() + (seg[1].y() - seg[0].y()) * (j + 1) / (double)number_of_division,
                            seg[0].z() + (seg[1].z() - seg[0].z()) * (j + 1) / (double)number_of_division);
              
              Plane ortho = seg.supporting_line().perpendicular_plane(anchor); 

              std::vector<Point> pts1, pts2;
              //Computation of the permanent angle and directions
              for (std::size_t k = 0; k < division_tab[j].size(); ++ k)
                { 
                  std::size_t inde = division_tab[j][k];
                  std::size_t plane = m_indices[inde];
                  if (plane == m_edges[i].planes[0])
                    pts1.push_back (m_points[inde]);
                  else if (plane == m_edges[i].planes[1])
                    pts2.push_back (m_points[inde]);
                }

              Point centroid1 = CGAL::centroid (pts1.begin (), pts1.end ());
              Point centroid2 = CGAL::centroid (pts2.begin (), pts2.end ());

              Line line_p1;
              CGAL::Object ob_temp1 = CGAL::intersection (static_cast<Plane> (*plane1), ortho);
              if (!assign(line_p1, ob_temp1))
                std::cout<<"Warning: bad plane/plane intersection"<<std::endl;
              else
                {
                  Vector vecp1 = line_p1.to_vector();
                  vecp1 = vecp1/ std::sqrt (vecp1 * vecp1);
                  Vector vtest1 (anchor, centroid1);
                  if (vtest1 * vecp1<0)
                    vecp1 = -vecp1;

                  direction_p1 = direction_p1+vecp1;

                  Point anchor1 = anchor + vecp1 * r_edge;
                  m_points.push_back (anchor1);
                  m_indices.push_back (m_edges[i].planes[0]);
                  m_status.push_back (POINT);
                }

              Line line_p2;
              CGAL::Object ob_temp2 = CGAL::intersection (static_cast<Plane> (*plane2),ortho);
              if (!assign(line_p2, ob_temp2))
                std::cout<<"Warning: bad plane/plane intersection"<<std::endl;
              else
                {
                  Vector vecp2 = line_p2.to_vector();
                  vecp2 = vecp2 / std::sqrt (vecp2 * vecp2);
                  Vector vtest2 (anchor, centroid2);
                  if (vtest2 * vecp2 < 0)
                    vecp2 =- vecp2;

                  direction_p2 = direction_p2+vecp2;

                  Point anchor2 = anchor + vecp2 * r_edge;
                  m_points.push_back (anchor2);
                  m_indices.push_back (m_edges[i].planes[1]);
                  m_status.push_back (POINT);
                }

            }

          //if not information enough (not enough edges to create
          //anchor) we unactivate the edge, else we update the angle
          //and directions
          if ( !(direction_p1.squared_length()>0 || direction_p2.squared_length()>0) )
            {
              m_edges[i].active = false;
              for (std::size_t j = 0; j < m_edges[i].indices.size (); ++ j)
                m_status[m_edges[i].indices[j]] = SKIPPED;
            }
        }
    }

    void compute_corners (double radius)
    {

      for (std::size_t i = 0; i < m_edges.size () - 2; ++ i)
        {
          if (!(m_edges[i].active))
            continue;

          for (std::size_t j = i + 1; j < m_edges.size () - 1; ++ j)
            {
              if (!(m_edges[j].active))
                continue;

              for (std::size_t k = j + 1; k < m_edges.size (); ++ k)
                {
                  if (!(m_edges[k].active))
                    continue;

                  std::set<std::size_t> planes;
                  planes.insert (m_edges[i].planes[0]);
                  planes.insert (m_edges[i].planes[1]);
                  planes.insert (m_edges[j].planes[0]);
                  planes.insert (m_edges[j].planes[1]);
                  planes.insert (m_edges[k].planes[0]);
                  planes.insert (m_edges[k].planes[1]);

                  if (planes.size () == 3) // Triple found
                    {
                      std::vector<std::size_t> vecplanes (planes.begin (), planes.end ());
                      m_corners.push_back (Corner (vecplanes[0], vecplanes[1], vecplanes[2],
                                                   i, j, k));
                    }
                }
            }
        }

      for (std::size_t i = 0; i < m_corners.size (); ++ i)
        {
          //calcul pt d'intersection des 3 plans
          Plane plane1 = static_cast<Plane> (*(m_planes[m_corners[i].planes[0]]));
          Plane plane2 = static_cast<Plane> (*(m_planes[m_corners[i].planes[1]]));
          Plane plane3 = static_cast<Plane> (*(m_planes[m_corners[i].planes[2]]));
          Line line;

          CGAL::Object ob_temp = CGAL::intersection(plane1, plane2);
          if (!assign(line, ob_temp))
            {
              std::cerr << "Warning: bad plane/plane intersection" << std::endl;
              continue;
            }
          else
            {
              CGAL::Object ob_temp2 = CGAL::intersection (line, plane3);
              if (!assign (m_corners[i].support, ob_temp2))
                {
                  std::cerr << "Warning: bad plane/line intersection" << std::endl;
                  continue;
                }
            }

          // test if point is in bbox + delta
          CGAL::Bbox_3 bbox = CGAL::bbox_3 (m_points.begin (), m_points.end ());
          
          double margin_x = 0.1 * (bbox.xmax() - bbox.xmin());
          double X_min = bbox.xmin() - margin_x;
          double X_max = bbox.xmax() + margin_x; 
          double margin_y = 0.1 * (bbox.ymax() - bbox.ymin());
          double Y_min = bbox.ymin() - margin_y;
          double Y_max = bbox.ymax() + margin_y; 
          double margin_z = 0.1* (bbox.zmax() - bbox.zmin());
          double Z_min = bbox.zmin() - margin_z;
          double Z_max = bbox.zmax() + margin_z;
          
          if ((m_corners[i].support.x() < X_min) || (m_corners[i].support.x() > X_max)
              || (m_corners[i].support.y() < Y_min) || (m_corners[i].support.y() > Y_max)
              || (m_corners[i].support.z() < Z_min) || (m_corners[i].support.z() > Z_max))
            break;

          // test if corner is in neighborhood of at least one point each of the 3 planes
          std::vector<bool> neighborhood (3, false);

          for (std::size_t k = 0; k < 3; ++ k)
            {
              for (std::size_t j = 0; j < m_edges[m_corners[i].edges[k]].indices.size(); ++ j)
                {
                  const Point& p = m_points[m_edges[m_corners[i].edges[k]].indices[j]];

                  if (CGAL::squared_distance (m_corners[i].support, p) < radius * radius)
                    {
                      neighborhood[k] = true;
                      break;
                    }
                }
            }

          if ( !(neighborhood[0] && neighborhood[1] && neighborhood[2]) )
            m_corners[i].active = false;
        }
    }

    void merge_corners (double radius)
    {
      for (std::size_t k = 0; k < m_corners.size(); ++ k)
        {
          if (!(m_corners[k].active))
            continue;

          int count_plane_number=3;
          
          for (std::size_t kb = k + 1; kb < m_corners.size(); ++ kb)
            {
              if (!(m_corners[kb].active))
                continue;

              int count_new_plane = 0;

              if (CGAL::squared_distance (m_corners[kb].support, m_corners[k].support) >= radius * radius)
                continue;

              for (std::size_t i = 0; i < m_corners[kb].planes.size (); ++ i)
                {
                  bool testtt = true; 
                  for (std::size_t l = 0; l < m_corners[k].planes.size(); ++ l)
                    if (m_corners[kb].planes[i] == m_corners[k].planes[l])
                      {
                        testtt = false;
                        break;
                      }
                  if (!testtt)
                    continue;

                  m_corners[k].planes.push_back (m_corners[kb].planes[i]);
                  ++ count_new_plane;
                  m_corners[kb].active = false;

                  std::vector<bool> is_edge_in (3, false);
                  for (std::size_t l = 0; l < m_corners[k].edges.size(); ++ l)
                    {
                      for (std::size_t j = 0; j < 3; ++ i)
                        if (m_corners[k].edges[l] == m_corners[kb].edges[j])
                          is_edge_in[j] = true;
                    }
                  for (std::size_t j = 0; j < 3; ++ i)
                    if (!(is_edge_in[j]))
                      m_corners[k].edges.push_back (m_corners[kb].edges[j]);

                }
              
              //update barycenter
              m_corners[k].support = CGAL::barycenter (m_corners[k].support, count_plane_number,
                                                       m_corners[kb].support, count_new_plane);
              count_plane_number += count_new_plane;
            }

          m_points.push_back (m_corners[k].support);
          m_indices.push_back (minus1);
          m_status.push_back (CORNER);
        }
    }
    
  };
  
} /* namespace internal */
/// \endcond



// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------

/// \ingroup PkgPointSetProcessing
/// TODO documentation

// This variant requires the kernel.
template <typename InputIterator,
          typename PointPMap,
          typename EfficientRANSACTraits,
          typename Kernel
>
void
structure_point_set (InputIterator first,  ///< iterator over the first input point.
                     InputIterator beyond, ///< past-the-end iterator over the input points.
                     PointPMap point_pmap, ///< property map: value_type of InputIterator -> Point_3
                     Shape_detection_3::Efficient_RANSAC<EfficientRANSACTraits>&
                     shape_detection, ///< shape detection engine
                     double radius, ///< attraction radius
                     const Kernel& /*kernel*/) ///< geometric traits.
{
  internal::Point_set_structuring<EfficientRANSACTraits> pss
    (first, beyond, shape_detection);
  pss.run (radius);
}

/// @cond SKIP_IN_MANUAL
// This variant deduces the kernel from the iterator type.
template <typename InputIterator,
          typename PointPMap,
          typename EfficientRANSACTraits
>
void
structure_point_set (InputIterator first,    ///< iterator over the first input point.
                     InputIterator beyond,   ///< past-the-end iterator over the input points.
                     PointPMap point_pmap, ///< property map: value_type of InputIterator -> Point_3
                     Shape_detection_3::Efficient_RANSAC<EfficientRANSACTraits>&
                     shape_detection, ///< shape detection engine
                     double radius) ///< attraction radius
{
  typedef typename boost::property_traits<PointPMap>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel Kernel;
  return structure_point_set (
    first,beyond,
    point_pmap,
    shape_detection,
    radius,
    Kernel());
}
/// @endcond

/// @cond SKIP_IN_MANUAL
// This variant creates a default point property map = Identity_property_map.
template < typename InputIterator, typename EfficientRANSACTraits >
void
structure_point_set (InputIterator first,    ///< iterator over the first input point.
                     InputIterator beyond,   ///< past-the-end iterator over the input points.
                     Shape_detection_3::Efficient_RANSAC<EfficientRANSACTraits>&
                     shape_detection, ///< shape detection engine
                     double radius) ///< attraction radius
{
  return structure_point_set (
    first,beyond,
#ifdef CGAL_USE_PROPERTY_MAPS_API_V1
    make_dereference_property_map(first),
#else
    make_identity_property_map(
    typename std::iterator_traits<InputIterator>::value_type()),
#endif
    shape_detection,
    radius);
}
/// @endcond


} //namespace CGAL

#endif // CGAL_STRUCTURE_POINT_SET_3_H

