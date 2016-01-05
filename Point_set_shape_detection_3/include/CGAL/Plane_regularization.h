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
// Author(s)     : Florent Lafarge, Simon Giraudot
//

/**
 * \ingroup PkgPointSetShapeDetection3
 * \file CGAL/Plane_regularization.h
 *
 */


#ifndef CGAL_PLANE_REGULARIZATION_H
#define CGAL_PLANE_REGULARIZATION_H

#include <CGAL/Shape_detection_3.h>
#include <CGAL/centroid.h>

#include <boost/foreach.hpp>


namespace CGAL {

  
/*!
\ingroup PkgPointSetShapeDetection3
\brief A plane regularization algorithm applied as a post-processing
to a shape detection algorithm.

Given a set of detected planes with their respective inlier sets, this
class enables to regularize the planes: planes almost parallel are
made exactly parallel. In addition, some additional regularization can
be performed:

- Plane clusters that are almost orthogonal can be made exactly
  orthogonal.

- Planes that are parallel and almost coplanar can be made exactly
  coplanar.

- Planes that are almost symmetrical with a user-defined axis can be
  made exactly symmetrical.

Planes are directly modified. Points are left unaltered, as well as
their relationships to planes (no transfer of point from a primitive
plane to another).

The implementation follows \cgalCite{cgal:vla-lod-15}.

\tparam Traits a model of `EfficientRANSACTraits`

*/
  template <typename Traits>
  class Plane_regularization
  {
  public:

    /// \cond SKIP_IN_MANUAL
    typedef Plane_regularization<Traits> Self;
    /// \endcond

    /// \name Types 
    /// @{
    /// \cond SKIP_IN_MANUAL
    typedef typename Traits::FT FT;
    typedef typename Traits::Point_3 Point;
    typedef typename Traits::Vector_3 Vector;
    typedef typename Traits::Line_3 Line;

    /// \endcond
    typedef typename Traits::Plane_3 Plane; ///< Raw plane type

    typedef typename Traits::Point_map Point_map;
    ///< property map to access the location of an input point.
    typedef typename Traits::Normal_map Normal_map;
    ///< property map to access the unoriented normal of an input point
    typedef typename Traits::Input_range Input_range;
    ///< Model of the concept `Range` with random access iterators, providing input points and normals
    /// through the following two property maps.

    typedef typename Input_range::iterator Input_iterator; ///< Iterator on input data

    typedef Shape_detection_3::Shape_base<Traits> Shape; ///< Shape type.
    typedef Shape_detection_3::Plane<Traits> Plane_shape; ///< Plane type.
    /// @}

  private:

    struct Plane_cluster
    {
      bool is_free;
      std::vector<std::size_t> planes;
      std::vector<std::size_t> coplanar_group;
      std::vector<std::size_t> orthogonal_clusters;
      Vector normal;
      FT cosangle_symmetry;
      FT area;
      FT cosangle_centroid;
    };

    Traits m_traits;

    Input_iterator m_input_begin;
    Input_iterator m_input_end;
    Point_map m_point_pmap;
    Normal_map m_normal_pmap;

    std::vector<boost::shared_ptr<Plane_shape> > m_planes;
    std::vector<Point> m_centroids;
    std::vector<FT> m_areas;

  public:

    /// \name Initialization 
    /// @{
    /*! 
      Constructs an empty plane regularization engine.
    */ 
    Plane_regularization (Traits t = Traits ())
      : m_traits (t)
    {

    }

    /*! 

      Constructs a plane regularization engine base on an input range
      of points with its related shape detection engine.

      \param input_range Range of input data.

      \param shape_detection Shape detection engine used to detect
      shapes from the input data. This engine may handle any types of
      primitive shapes but only planes will be regularized.

      \warning The `shape_detection` parameter must have already
      detected shapes and must have been using `input_range` as input.

    */ 
    Plane_regularization (Input_range& input_range,
                          const Shape_detection_3::Efficient_RANSAC<Traits>& shape_detection)
      : m_traits (shape_detection.traits())
    {
      m_input_begin = input_range.begin ();
      m_input_end = input_range.end ();

      BOOST_FOREACH (boost::shared_ptr<Shape> shape, shape_detection.shapes())
        {
          boost::shared_ptr<Plane_shape> pshape
            = boost::dynamic_pointer_cast<Plane_shape>(shape);
        
          // Ignore all shapes other than plane
          if (pshape == boost::shared_ptr<Plane_shape>())
            continue;
          m_planes.push_back (pshape);
        }

    }

    /*! 
      Releases all memory allocated by this instance.
    */ 
    virtual ~Plane_regularization ()
    {
      clear ();
    }
    /// @}

    /// \name Memory Management
    /// @{
    /*!
      Clear all internal structures.
     */ 
    void clear ()
    {
      std::vector<boost::shared_ptr<Plane_shape> > ().swap (m_planes);
      std::vector<Point> ().swap (m_centroids);
      std::vector<FT> ().swap (m_areas);
    }
    /// @}

    /// \name Regularization 
    /// @{
    /*! 

      Performs the plane regularization. Planes are directly modified.

      \param tolerance_angle Tolerance of deviation between normal
      vectors of planes so that they are considered parallel (in
      degrees).

      \param tolerance_coplanarity Maximal distance between two
      parallel planes such that they are considered coplanar. The
      default value is 0, meaning that coplanarity is not taken into
      account for regularization.

      \param regularize_orthogonality Make almost orthogonal clusters
      of plane exactly orthogonal.

      \param symmetry_direction Make clusters that are almost
      symmetrical in the symmetry direction exactly symmetrical. This
      parameter is ignored if it is equal to `CGAL::NULL_VECTOR`
      (default value).

      \return The number of clusters of parallel planes found.
    */ 

    std::size_t run (FT tolerance_angle = (FT)25.0,
                     FT tolerance_coplanarity = (FT)0.0,
                     bool regularize_orthogonality = true,
                     Vector symmetry_direction = CGAL::NULL_VECTOR)
    {
      compute_centroids_and_areas ();

      FT tolerance_cosangle = (FT)1. - std::cos (tolerance_angle);
      
      // clustering the parallel primitives and store them in clusters
      // & compute the normal, size and cos angle to the symmetry
      // direction of each cluster
      std::vector<Plane_cluster> clusters;
      compute_parallel_clusters (clusters, tolerance_cosangle, symmetry_direction);

      if (regularize_orthogonality)
        {
          //discovery orthogonal relationship between clusters 
          for (std::size_t i = 0; i < clusters.size(); ++ i)
            {
              for (std::size_t j = i + 1; j < clusters.size(); ++ j)
                {
              
                  if (std::fabs (clusters[i].normal * clusters[j].normal) < tolerance_cosangle)
                    {
                      clusters[i].orthogonal_clusters.push_back (j);
                      clusters[j].orthogonal_clusters.push_back (i);
                    }
                }
            }
        }
      
      //clustering the symmetry cosangle and store their centroids in
      //cosangle_centroids and the centroid index of each cluster in
      //list_cluster_index
      if (symmetry_direction != CGAL::NULL_VECTOR)
        cluster_symmetric_cosangles (clusters, tolerance_cosangle);

      //find subgraphs of mutually orthogonal clusters (store index of
      //clusters in subgraph_clusters), and select the cluster of
      //largest area
      if (regularize_orthogonality)
        subgraph_mutually_orthogonal_clusters (clusters, symmetry_direction);
      
      //recompute optimal plane for each primitive after normal regularization
      for (std::size_t i=0; i < clusters.size(); ++ i)
        {

          Vector vec_reg = clusters[i].normal;

          for (std::size_t j = 0; j < clusters[i].planes.size(); ++ j)
            {
              std::size_t index_prim = clusters[i].planes[j];
              Point pt_reg = m_planes[index_prim]->projection (m_centroids[index_prim]);
              if( m_planes[index_prim]->plane_normal () * vec_reg < 0)
                vec_reg=-vec_reg;
              Plane plane_reg(pt_reg,vec_reg);
              
              if( std::fabs(m_planes[index_prim]->plane_normal () * plane_reg.orthogonal_vector ()) > 1. - tolerance_cosangle)
                m_planes[index_prim]->update (plane_reg);
            }
        }


      //detecting co-planarity and store in list_coplanar_prim
      for (std::size_t i = 0; i < clusters.size(); ++ i)
        {
          Vector vec_reg = clusters[i].normal;

          for (std::size_t ip = 0; ip < clusters[i].planes.size(); ++ ip)
            clusters[i].coplanar_group.push_back (static_cast<std::size_t>(-1));

          std::size_t cop_index=0;

          for (std::size_t j = 0; j < clusters[i].planes.size(); ++ j)
            {
              std::size_t index_prim = clusters[i].planes[j];

              if (clusters[i].coplanar_group[j] == static_cast<std::size_t>(-1))
                {
                  clusters[i].coplanar_group[j] = cop_index;
			
                  Point pt_reg = m_planes[index_prim]->projection(m_centroids[index_prim]);
                  Plane plan_reg(pt_reg,vec_reg);

                  for (std::size_t k = j + 1; k < clusters[i].planes.size(); ++ k)
                    {
                      if (clusters[i].coplanar_group[k] == static_cast<std::size_t>(-1))
                        {
                          std::size_t index_prim_next = clusters[i].planes[k];
                          Point pt_reg_next = m_planes[index_prim_next]->projection(m_centroids[index_prim_next]);
                          Point pt_proj=plan_reg.projection(pt_reg_next);
                          FT distance=distance_Point(pt_reg_next,pt_proj);
					
                          if (distance < tolerance_coplanarity)
                            clusters[i].coplanar_group[k] = cop_index;
                        }
                    }
                  cop_index++; 
                }
            }

          //regularize primitive position by computing barycenter of cplanar planes
          std::vector<Point> pt_bary (cop_index, Point (0., 0., 0.));
          std::vector<FT> area (cop_index, 0.);
      
          for (std::size_t j = 0; j < clusters[i].planes.size (); ++ j)
            {
              std::size_t index_prim = clusters[i].planes[j];
              std::size_t group = clusters[i].coplanar_group[j];
              
              Point pt_reg = m_planes[index_prim]->projection(m_centroids[index_prim]);

              pt_bary[group] = CGAL::barycenter (pt_bary[group], area[group], pt_reg, m_areas[index_prim]); 
              area[group] += m_areas[index_prim];
            }


          for (std::size_t j = 0; j < clusters[i].planes.size (); ++ j)
            {
              std::size_t index_prim = clusters[i].planes[j];
              std::size_t group = clusters[i].coplanar_group[j];
              
              Plane plane_reg (pt_bary[group], vec_reg);
              
              if (m_planes[index_prim]->plane_normal ()
                  * plane_reg.orthogonal_vector() < 0)
                m_planes[index_prim]->update (plane_reg.opposite());
              else
                m_planes[index_prim]->update (plane_reg);
            }
        }
      
      return clusters.size ();
    }
    /// @}


  private:
    
    void compute_centroids_and_areas ()
    {
      for (std::size_t i = 0; i < m_planes.size (); ++ i)
        {
          std::vector < Point > listp;
          for (std::size_t j = 0; j < m_planes[i]->indices_of_assigned_points ().size (); ++ j)
            {
              std::size_t yy = m_planes[i]->indices_of_assigned_points()[j];
              Point pt = get (m_point_pmap, *(m_input_begin + yy));
              listp.push_back(pt);
            }
          m_centroids.push_back (CGAL::centroid (listp.begin (), listp.end ()));
          m_areas.push_back ((FT)(m_planes[i]->indices_of_assigned_points().size()) / (FT)100.);
        }
    }

    void compute_parallel_clusters (std::vector<Plane_cluster>& clusters, FT tolerance_cosangle,
                                    const Vector& symmetry_direction)
    {
      // find pairs of epsilon-parallel primitives and store them in parallel_planes
      std::vector<std::vector<std::size_t> > parallel_planes (m_planes.size ());
      for (std::size_t i = 0; i < m_planes.size (); ++ i)
        {
          Vector v1 = m_planes[i]->plane_normal ();
          
          for (std::size_t j = 0; j < m_planes.size(); ++ j)
            {
              if (i == j)
                continue;
              
              Vector v2 = m_planes[i]->plane_normal ();
              
              if (std::fabs (v1 * v2) > 1. - tolerance_cosangle)
                parallel_planes[i].push_back (j);
            }
        }


      std::vector<bool> is_available (m_planes.size (), true);
      
      for (std::size_t i = 0; i < m_planes.size(); ++ i)
        {

          if(is_available[i])
            {
              is_available[i] = false;

              clusters.push_back (Plane_cluster());
              Plane_cluster& clu = clusters.back ();

              //initialization containers
              clu.planes.push_back (i);
              
              std::vector<std::size_t> index_container_former_ring_parallel;
              index_container_former_ring_parallel.push_back(i);
              
              std::list<std::size_t> index_container_current_ring_parallel;

              //propagation over the pairs of epsilon-parallel primitives
              bool propagation=true;
              clu.normal = m_planes[i]->plane_normal ();
              clu.area = m_areas[i];
			
              do
                {
                  propagation = false;

                  for (std::size_t k = 0; k < index_container_former_ring_parallel.size(); ++ k)
                    {

                      std::size_t plane_index = index_container_former_ring_parallel[k];

                      for (std::size_t l = 0; l < parallel_planes[plane_index].size(); ++ l)
                        {
                          std::size_t it = parallel_planes[plane_index][l];
                          
                          Vector normal_it =  m_planes[it]->plane_normal ();

                          if(is_available[it]
                             && std::fabs (normal_it*clu.normal) > 1. - tolerance_cosangle )
                            {	
                              propagation = true;
                              index_container_current_ring_parallel.push_back(it);
                              is_available[it]=false;
                              
                              if(clu.normal * normal_it <0)
                                normal_it = -normal_it;

                              clu.normal = (FT)clu.area * clu.normal
                                + (FT)m_areas[it] * normal_it;
                              FT norm = (FT)1. / std::sqrt (clu.normal.squared_length()); 
                              clu.normal = norm * clu.normal;
                              clu.area += m_areas[it];
                            }	
                        }
                    }

                  //update containers
                  index_container_former_ring_parallel.clear();
                  for (std::list<std::size_t>::iterator it = index_container_current_ring_parallel.begin();
                       it != index_container_current_ring_parallel.end(); ++it)
                    {
                      index_container_former_ring_parallel.push_back(*it);
                      clu.planes.push_back(*it);
                    }
                  index_container_current_ring_parallel.clear();

                }
              while(propagation);

              if (symmetry_direction != CGAL::NULL_VECTOR)
                clu.cosangle_symmetry = std::fabs(symmetry_direction * clu.normal);
            }
        }
      is_available.clear();
    }

    void cluster_symmetric_cosangles (std::vector<Plane_cluster>& clusters, FT tolerance_cosangle)
    {
      std::vector < FT > cosangle_centroids;
      std::vector < std::size_t> list_cluster_index;
      for( std::size_t i = 0; i < clusters.size(); ++ i)
        list_cluster_index.push_back(static_cast<std::size_t>(-1));
      
      std::size_t mean_index = 0;
      for (std::size_t i = 0; i < clusters.size(); ++ i)
        {
          if(list_cluster_index[i] == static_cast<std::size_t>(-1))
            {
              list_cluster_index[i] = mean_index;
              FT mean = clusters[i].area * clusters[i].cosangle_symmetry;
              FT mean_area = clusters[i].area;
              
              for (std::size_t j = i+1; j < clusters.size(); ++ j)
                {
                  if (list_cluster_index[j] == static_cast<std::size_t>(-1)
                      && std::fabs (clusters[j].cosangle_symmetry -
                                    mean / mean_area) < tolerance_cosangle)
                    {
                      list_cluster_index[j] = mean_index;
                      mean_area += clusters[j].area;
                      mean += clusters[j].area * clusters[j].cosangle_symmetry;
                    }
                }
              ++ mean_index;
              mean /= mean_area;
              cosangle_centroids.push_back (mean);
            }
        }

      for (std::size_t i = 0; i < cosangle_centroids.size(); ++ i)
        {
          if (cosangle_centroids[i] < tolerance_cosangle)
            cosangle_centroids[i] = 0;
          else if (cosangle_centroids[i] > 1. - tolerance_cosangle)
            cosangle_centroids[i] = 1;
        }
      for (std::size_t i = 0; i < clusters.size(); ++ i)
        clusters[i].cosangle_symmetry = cosangle_centroids[list_cluster_index[i]];
    }

    void subgraph_mutually_orthogonal_clusters (std::vector<Plane_cluster>& clusters,
                                                const Vector& symmetry_direction)
    {
      std::vector < std::vector < std::size_t> > subgraph_clusters;
      std::vector < std::size_t> subgraph_clusters_max_area_index;

      for (std::size_t i = 0; i < clusters.size(); ++ i)
        clusters[i].is_free = true;

      for (std::size_t i = 0; i < clusters.size(); ++ i)
        {
          if(clusters[i].is_free)
            {
              clusters[i].is_free = false;
              FT max_area = clusters[i].area;
              std::size_t index_max_area = i;

              //initialization containers
              std::vector < std::size_t > index_container;
              index_container.push_back(i);
              std::vector < std::size_t > index_container_former_ring;
              index_container_former_ring.push_back(i);
              std::list < std::size_t > index_container_current_ring;

              //propagation
              bool propagation=true;
              do
                {
                  propagation=false;

                  //neighbors
                  for (std::size_t k=0;k<index_container_former_ring.size();k++)
                    {

                      std::size_t cluster_index=index_container_former_ring[k];

                      for (std::size_t j = 0; j < clusters[cluster_index].orthogonal_clusters.size(); ++ j)
                        {
                          if(clusters[j].is_free)
                            { 	
                              propagation = true;
                              index_container_current_ring.push_back(j);
                              clusters[j].is_free = false;

                              if(max_area < clusters[j].area)
                                {
                                  max_area = clusters[j].area;
                                  index_max_area = j;
                                }
                            }	
                        }
                    }

                  //update containers
                  index_container_former_ring.clear();
                  for(std::list < std::size_t>::iterator it = index_container_current_ring.begin();
                      it != index_container_current_ring.end(); ++it)
                    {
                      index_container_former_ring.push_back(*it);
                      index_container.push_back(*it);
                    }
                  index_container_current_ring.clear();

                }
              while(propagation);
              subgraph_clusters.push_back(index_container);
              subgraph_clusters_max_area_index.push_back(index_max_area);
            }
        }


      //create subgraphs of mutually orthogonal clusters in which the
      //largest cluster is excluded and store in
      //subgraph_clusters_prop
      std::vector < std::vector < std::size_t> > subgraph_clusters_prop;
      for (std::size_t i=0;i<subgraph_clusters.size(); i++)
        {
          std::size_t index=subgraph_clusters_max_area_index[i];
          std::vector < std::size_t> subgraph_clusters_prop_temp;
          for (std::size_t j=0;j<subgraph_clusters[i].size(); j++)
            if(subgraph_clusters[i][j]!=index)
              subgraph_clusters_prop_temp.push_back(subgraph_clusters[i][j]);

          subgraph_clusters_prop.push_back(subgraph_clusters_prop_temp);
        }



      //regularization of cluster normals : in eachsubgraph, we start
      //from the largest area cluster and we propage over the subgraph
      //by regularizing the normals of the clusters accorting to
      //orthogonality and cosangle to symmetry direction

      for (std::size_t i = 0; i < clusters.size(); ++ i)
        clusters[i].is_free = true;

      for (std::size_t i = 0; i < subgraph_clusters_prop.size(); ++ i)
        {
	
          std::size_t index_current=subgraph_clusters_max_area_index[i];
          Vector vec_current=regularize_normal(clusters[index_current].normal,
                                               symmetry_direction,
                                               clusters[index_current].cosangle_symmetry);
          clusters[index_current].normal = vec_current;
          clusters[index_current].is_free = false;

          //initialization containers
          std::vector < std::size_t> index_container;
          index_container.push_back(index_current);
          std::vector < std::size_t> index_container_former_ring;
          index_container_former_ring.push_back(index_current);
          std::list < std::size_t> index_container_current_ring;

          //propagation
          bool propagation=true;
          do
            {
              propagation=false;

              //neighbors
              for (std::size_t k=0;k<index_container_former_ring.size();k++)
                {

                  std::size_t cluster_index=index_container_former_ring[k];

                  for (std::size_t j = 0; j < clusters[cluster_index].orthogonal_clusters.size(); ++ j)
                    {
						
                      if(clusters[j].is_free)
                        { 	
							
                          propagation = true;
                          index_container_current_ring.push_back(j);
                          clusters[j].is_free = false;

                          Vector new_vect=regularize_normals_from_prior(clusters[cluster_index].normal,
                                                                        clusters[j].normal,
                                                                        symmetry_direction,
                                                                        clusters[j].cosangle_symmetry);
                          clusters[j].normal = new_vect;
                        }
                    }	
                }
			
              //update containers
              index_container_former_ring.clear();
              for(std::list < std::size_t>::iterator it = index_container_current_ring.begin();
                  it != index_container_current_ring.end(); ++it)
                {
                  index_container_former_ring.push_back(*it);
                  index_container.push_back(*it);
                }
              index_container_current_ring.clear();
            }while(propagation);
        }
    }
                                    

    FT distance_Point (const Point& a, const Point& b)
    {
      return std::sqrt (CGAL::squared_distance (a, b));
    }
  
    Vector regularize_normal (const Vector& n, const Vector& symmetry_direction,
                              FT cos_symmetry)
    {
      if (symmetry_direction == CGAL::NULL_VECTOR)
        return n;
      
      Point pt_symmetry = CGAL::ORIGIN + cos_symmetry* symmetry_direction;

      Plane plane_symmetry (pt_symmetry, symmetry_direction);
      Point pt_normal = CGAL::ORIGIN + n;

      if (n != symmetry_direction || n != -symmetry_direction)
        {
          Plane plane_cut (CGAL::ORIGIN, pt_normal, CGAL::ORIGIN + symmetry_direction);
          Line line;
          CGAL::Object ob_1 = CGAL::intersection(plane_cut, plane_symmetry);
          if (!assign(line, ob_1))
            return n;

          FT delta = std::sqrt ((FT)1. - cos_symmetry * cos_symmetry);

          Point projected_origin = line.projection (CGAL::ORIGIN);
          Vector line_vector (line);
          line_vector = line_vector / std::sqrt (line_vector * line_vector);
          Point pt1 = projected_origin + delta * line_vector;
          Point pt2 = projected_origin - delta * line_vector;

          if (CGAL::squared_distance (pt_normal, pt1) <= CGAL::squared_distance (pt_normal, pt2))
            return Vector (CGAL::ORIGIN, pt1);
          else
            return Vector (CGAL::ORIGIN, pt2);

      }
      else
        return n;
    }

  
    Vector regularize_normals_from_prior (const Vector& np,
                                          const Vector& n,
                                          const Vector& symmetry_direction,
                                          FT cos_symmetry)
    {
      if (symmetry_direction == CGAL::NULL_VECTOR)
        return n;

      Plane plane_orthogonality (CGAL::ORIGIN, np);
      Point pt_symmetry = CGAL::ORIGIN + cos_symmetry* symmetry_direction;

      Plane plane_symmetry (pt_symmetry, symmetry_direction);
		
      Line line;
      CGAL::Object ob_1 = CGAL::intersection (plane_orthogonality, plane_symmetry);
      if (!assign(line, ob_1))
        return regularize_normal (n, symmetry_direction, cos_symmetry);

      Point projected_origin = line.projection (CGAL::ORIGIN);
      FT R = CGAL::squared_distance (Point (CGAL::ORIGIN), projected_origin);

      if (R <= 1)  // 2 (or 1) possible points intersecting the unit sphere and line
        {
          FT delta = std::sqrt ((FT)1. - R);
          Vector line_vector(line); 
          line_vector = line_vector / std::sqrt (line_vector * line_vector);
          Point pt1 = projected_origin + delta * line_vector;
          Point pt2 = projected_origin - delta * line_vector;
			
          Point pt_n = CGAL::ORIGIN + n;
          if (CGAL::squared_distance (pt_n, pt1) <= CGAL::squared_distance (pt_n, pt2))
            return Vector (CGAL::ORIGIN, pt1);
          else
            return Vector (CGAL::ORIGIN, pt2);
        }
      else //no point intersecting the unit sphere and line
        return regularize_normal (n,symmetry_direction, cos_symmetry);

    }

  

  };


} // namespace CGAL

#endif // CGAL_PLANE_REGULARIZATION_H
