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
 * \file CGAL/regularize_planes.h
 *
 */


#ifndef CGAL_REGULARIZE_PLANES_H
#define CGAL_REGULARIZE_PLANES_H

#include <CGAL/license/Point_set_shape_detection_3.h>


#include <CGAL/Shape_detection_3.h>
#include <CGAL/centroid.h>
#include <CGAL/squared_distance_3.h>

#include <boost/foreach.hpp>


namespace CGAL {
  
// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
/// \cond SKIP_IN_MANUAL
namespace internal {
namespace PlaneRegularization {

template <typename Traits>
struct Plane_cluster
{
  bool is_free;
  std::vector<std::size_t> planes;
  std::vector<std::size_t> coplanar_group;
  std::vector<std::size_t> orthogonal_clusters;
  typename Traits::Vector_3 normal;
  typename Traits::FT cosangle_symmetry;
  typename Traits::FT area;
  typename Traits::FT cosangle_centroid;
};

  
template <typename Traits>
typename Traits::Vector_3 regularize_normal
  (const typename Traits::Vector_3& n,
   const typename Traits::Vector_3& symmetry_direction,
   typename Traits::FT cos_symmetry)
{
  typedef typename Traits::FT FT;
  typedef typename Traits::Point_3 Point;
  typedef typename Traits::Vector_3 Vector;
  typedef typename Traits::Line_3 Line;
  typedef typename Traits::Plane_3 Plane;

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

template <typename Traits>  
typename Traits::Vector_3 regularize_normals_from_prior
  (const typename Traits::Vector_3& np,
   const typename Traits::Vector_3& n,
   const typename Traits::Vector_3& symmetry_direction,
   typename Traits::FT cos_symmetry)
{
  typedef typename Traits::FT FT;
  typedef typename Traits::Point_3 Point;
  typedef typename Traits::Vector_3 Vector;
  typedef typename Traits::Line_3 Line;
  typedef typename Traits::Plane_3 Plane;

  Plane plane_orthogonality (CGAL::ORIGIN, np);
  Point pt_symmetry = CGAL::ORIGIN + cos_symmetry* symmetry_direction;

  Plane plane_symmetry (pt_symmetry, symmetry_direction);
		
  Line line;
  CGAL::Object ob_1 = CGAL::intersection (plane_orthogonality, plane_symmetry);
  if (!assign(line, ob_1))
    return regularize_normal<Traits> (n, symmetry_direction, cos_symmetry);

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
    return regularize_normal<Traits> (n,symmetry_direction, cos_symmetry);

}

template <typename Traits,
          typename RandomAccessIterator,
          typename PlaneContainer,
          typename PointPMap,
          typename CentroidContainer,
          typename AreaContainer>
void compute_centroids_and_areas (RandomAccessIterator input_begin,
                                  PlaneContainer& planes,
                                  PointPMap point_pmap,
                                  CentroidContainer& centroids,
                                  AreaContainer& areas)
{
  typedef typename Traits::FT FT;
  typedef typename Traits::Point_3 Point;
  
  for (std::size_t i = 0; i < planes.size (); ++ i)
    {
      std::vector < Point > listp;
      for (std::size_t j = 0; j < planes[i]->indices_of_assigned_points ().size (); ++ j)
        {
          std::size_t yy = planes[i]->indices_of_assigned_points()[j];
          Point pt = get (point_pmap, *(input_begin + yy));
          listp.push_back(pt);
        }
      centroids.push_back (CGAL::centroid (listp.begin (), listp.end ()));
      areas.push_back ((FT)(planes[i]->indices_of_assigned_points().size()) / (FT)100.);
    }
}


template <typename Traits,
          typename PlaneContainer,
          typename PlaneClusterContainer,
          typename AreaContainer>
void compute_parallel_clusters (PlaneContainer& planes,
                                PlaneClusterContainer& clusters,
                                AreaContainer& areas,
                                typename Traits::FT tolerance_cosangle,
                                const typename Traits::Vector_3& symmetry_direction)
{

  typedef typename Traits::FT FT;
  typedef typename Traits::Vector_3 Vector;
  
  typedef typename PlaneClusterContainer::value_type Plane_cluster;
  
  // find pairs of epsilon-parallel primitives and store them in parallel_planes
  std::vector<std::vector<std::size_t> > parallel_planes (planes.size ());
  for (std::size_t i = 0; i < planes.size (); ++ i)
    {
      Vector v1 = planes[i]->plane_normal ();
          
      for (std::size_t j = 0; j < planes.size(); ++ j)
        {
          if (i == j)
            continue;
              
          Vector v2 = planes[j]->plane_normal ();

          if (std::fabs (v1 * v2) > 1. - tolerance_cosangle)
            parallel_planes[i].push_back (j);
        }
    }


  std::vector<bool> is_available (planes.size (), true);
      
  for (std::size_t i = 0; i < planes.size(); ++ i)
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
          clu.normal = planes[i]->plane_normal ();
          clu.area = areas[i];
			
          do
            {
              propagation = false;

              for (std::size_t k = 0; k < index_container_former_ring_parallel.size(); ++ k)
                {

                  std::size_t plane_index = index_container_former_ring_parallel[k];

                  for (std::size_t l = 0; l < parallel_planes[plane_index].size(); ++ l)
                    {
                      std::size_t it = parallel_planes[plane_index][l];
                          
                      Vector normal_it =  planes[it]->plane_normal ();

                      if(is_available[it]
                         && std::fabs (normal_it*clu.normal) > 1. - tolerance_cosangle )
                        {	
                          propagation = true;
                          index_container_current_ring_parallel.push_back(it);
                          is_available[it]=false;
                              
                          if(clu.normal * normal_it <0)
                            normal_it = -normal_it;

                          clu.normal = (FT)clu.area * clu.normal
                            + (FT)areas[it] * normal_it;
                          FT norm = (FT)1. / std::sqrt (clu.normal.squared_length()); 
                          clu.normal = norm * clu.normal;
                          clu.area += areas[it];
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
            {
              clu.cosangle_symmetry = symmetry_direction * clu.normal;
              if (clu.cosangle_symmetry < 0.)
                {
                  clu.normal = -clu.normal;
                  clu.cosangle_symmetry = -clu.cosangle_symmetry;
                }
            }
        }
    }

  is_available.clear();
}

template <typename Traits,
          typename PlaneClusterContainer>
void cluster_symmetric_cosangles (PlaneClusterContainer& clusters,
                                  typename Traits::FT tolerance_cosangle,
                                  typename Traits::FT tolerance_cosangle_ortho)
{
  typedef typename Traits::FT FT;
  
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
                                mean / mean_area) < tolerance_cosangle_ortho)
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
      if (cosangle_centroids[i] < tolerance_cosangle_ortho)
        cosangle_centroids[i] = 0;
      else if (cosangle_centroids[i] > 1. - tolerance_cosangle)
        cosangle_centroids[i] = 1;
    }
  for (std::size_t i = 0; i < clusters.size(); ++ i)
    clusters[i].cosangle_symmetry = cosangle_centroids[list_cluster_index[i]];
}


template <typename Traits,
          typename PlaneClusterContainer>
void subgraph_mutually_orthogonal_clusters (PlaneClusterContainer& clusters,
                                            const typename Traits::Vector_3& symmetry_direction)
{
  typedef typename Traits::FT FT;
  typedef typename Traits::Vector_3 Vector;
  
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
                      std::size_t cluster_index_2 = clusters[cluster_index].orthogonal_clusters[j];
                      if(clusters[cluster_index_2].is_free)
                        {
                          propagation = true;
                          index_container_current_ring.push_back(cluster_index_2);
                          clusters[cluster_index_2].is_free = false;

                          if(max_area < clusters[cluster_index_2].area)
                            {
                              max_area = clusters[cluster_index_2].area;
                              index_max_area = cluster_index_2;
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

      Vector vec_current=regularize_normal<Traits>
        (clusters[index_current].normal,
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
                  std::size_t cluster_index_2 = clusters[cluster_index].orthogonal_clusters[j];						
                  if(clusters[cluster_index_2].is_free)
                    {
                      propagation = true;
                      index_container_current_ring.push_back(cluster_index_2);
                      clusters[cluster_index_2].is_free = false;

                      Vector new_vect=regularize_normals_from_prior<Traits>
                        (clusters[cluster_index].normal,
                         clusters[cluster_index_2].normal,
                         symmetry_direction,
                         clusters[cluster_index_2].cosangle_symmetry);
                      clusters[cluster_index_2].normal = new_vect;
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
                                    


} // namespace PlaneRegularization
} // namespace internal
/// \endcond


// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------

/// \ingroup PkgPointSetShapeDetection3
  
  /*! 

    Given a set of detected planes with their respective inlier sets,
    this function enables to regularize the planes: 

    - Planes near parallel can be made exactly parallel.

    - Planes near orthogonal can be made exactly orthogonal.

    - Planes parallel and near coplanar can be made exactly coplanar.

    - Planes near symmetrical with a user-defined axis can be made
    exactly symmetrical.

    Planes are directly modified. Points are left unaltered, as well as
    their relationships to planes (no transfer of point from a primitive
    plane to another).

    The implementation follows \cgalCite{cgal:vla-lod-15}.

    \tparam Traits a model of `EfficientRANSACTraits`

    \param shape_detection Shape detection object used to detect
    shapes from the input data. While the shape detection algorithm
    deals with several types of primitive shapes only planes can be
    regularized.

    \warning The `shape_detection` parameter must have already
    detected shapes. If no plane exists in it, the regularization
    function doesn't do anything.

    \param regularize_parallelism Select whether parallelism is
    regularized or not.

    \param regularize_orthogonality Select whether orthogonality is
    regularized or not.

    \param regularize_coplanarity Select whether coplanarity is
    regularized or not.

    \param regularize_axis_symmetry Select whether axis symmetry is
    regularized or not.

    \param tolerance_angle Tolerance of deviation between normal
    vectors of planes (in degrees) used for parallelism, orthogonality
    and axis symmetry. Default value is 25 degrees.

    \param tolerance_coplanarity Maximal distance between two parallel
    planes such that they are considered coplanar. Default value is
    0.01.

    \param symmetry_direction Chosen axis for symmetry
    regularization. Default value is the Z axis.
*/ 

template <typename EfficientRANSACTraits>
void regularize_planes (const Shape_detection_3::Efficient_RANSAC<EfficientRANSACTraits>& shape_detection,
                        bool regularize_parallelism,
                        bool regularize_orthogonality,
                        bool regularize_coplanarity,
                        bool regularize_axis_symmetry,
                        typename EfficientRANSACTraits::FT tolerance_angle
                        = (typename EfficientRANSACTraits::FT)25.0,
                        typename EfficientRANSACTraits::FT tolerance_coplanarity
                        = (typename EfficientRANSACTraits::FT)0.01,
                        typename EfficientRANSACTraits::Vector_3 symmetry_direction
                        = typename EfficientRANSACTraits::Vector_3
                        ((typename EfficientRANSACTraits::FT)0.,
                         (typename EfficientRANSACTraits::FT)0.,
                         (typename EfficientRANSACTraits::FT)1.))
{
  typedef typename EfficientRANSACTraits::FT FT;
  typedef typename EfficientRANSACTraits::Point_3 Point;
  typedef typename EfficientRANSACTraits::Vector_3 Vector;
  typedef typename EfficientRANSACTraits::Plane_3 Plane;

  typedef Shape_detection_3::Shape_base<EfficientRANSACTraits> Shape;
  typedef Shape_detection_3::Plane<EfficientRANSACTraits> Plane_shape;

  typedef typename internal::PlaneRegularization::Plane_cluster<EfficientRANSACTraits>
    Plane_cluster;

  typename EfficientRANSACTraits::Input_range::iterator input_begin = shape_detection.input_iterator_first();

  std::vector<boost::shared_ptr<Plane_shape> > planes;
    
  BOOST_FOREACH (boost::shared_ptr<Shape> shape, shape_detection.shapes())
    {
      boost::shared_ptr<Plane_shape> pshape
        = boost::dynamic_pointer_cast<Plane_shape>(shape);
        
      // Ignore all shapes other than plane
      if (pshape == boost::shared_ptr<Plane_shape>())
        continue;
      planes.push_back (pshape);
    }


  /*
   * Compute centroids and areas
   */
  std::vector<Point> centroids;
  std::vector<FT> areas;
  internal::PlaneRegularization::compute_centroids_and_areas<EfficientRANSACTraits>
    (input_begin, planes, shape_detection.point_map(), centroids, areas);

  tolerance_angle = tolerance_angle * (FT)CGAL_PI / (FT)(180);
  FT tolerance_cosangle = (FT)1. - std::cos (tolerance_angle);
  FT tolerance_cosangle_ortho = std::cos ((FT)0.5 * (FT)CGAL_PI - tolerance_angle);
      
  // clustering the parallel primitives and store them in clusters
  // & compute the normal, size and cos angle to the symmetry
  // direction of each cluster
  std::vector<Plane_cluster> clusters;
  internal::PlaneRegularization::compute_parallel_clusters<EfficientRANSACTraits>
    (planes, clusters, areas,
     (regularize_parallelism ? tolerance_cosangle : (FT)0.0),
     (regularize_axis_symmetry ? symmetry_direction : CGAL::NULL_VECTOR));

  if (regularize_orthogonality)
    {
      //discovery orthogonal relationship between clusters 
      for (std::size_t i = 0; i < clusters.size(); ++ i)
        {
          for (std::size_t j = i + 1; j < clusters.size(); ++ j)
            {
              if (std::fabs (clusters[i].normal * clusters[j].normal) < tolerance_cosangle_ortho)
                {
                  clusters[i].orthogonal_clusters.push_back (j);
                  clusters[j].orthogonal_clusters.push_back (i);
                }
            }
        }
    }
      
  if (regularize_axis_symmetry)
    {
      //clustering the symmetry cosangle and store their centroids in
      //cosangle_centroids and the centroid index of each cluster in
      //list_cluster_index
      internal::PlaneRegularization::cluster_symmetric_cosangles<EfficientRANSACTraits>
        (clusters, tolerance_cosangle, tolerance_cosangle_ortho);
    }
  
  //find subgraphs of mutually orthogonal clusters (store index of
  //clusters in subgraph_clusters), and select the cluster of
  //largest area
  if (regularize_orthogonality || regularize_axis_symmetry)
    internal::PlaneRegularization::subgraph_mutually_orthogonal_clusters<EfficientRANSACTraits>
      (clusters, (regularize_axis_symmetry ? symmetry_direction : CGAL::NULL_VECTOR));
      
  //recompute optimal plane for each primitive after normal regularization
  for (std::size_t i=0; i < clusters.size(); ++ i)
    {

      Vector vec_reg = clusters[i].normal;
      for (std::size_t j = 0; j < clusters[i].planes.size(); ++ j)
        {
          std::size_t index_prim = clusters[i].planes[j];
          Point pt_reg = planes[index_prim]->projection (centroids[index_prim]);
          if( planes[index_prim]->plane_normal () * vec_reg < 0)
            vec_reg=-vec_reg;
          Plane plane_reg(pt_reg,vec_reg);

          if( std::fabs(planes[index_prim]->plane_normal () * vec_reg) > 1. - tolerance_cosangle)
            planes[index_prim]->update (plane_reg);
        }
    }


  if (regularize_coplanarity)
    {
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
			
                  Point pt_reg = planes[index_prim]->projection(centroids[index_prim]);
                  Plane plan_reg(pt_reg,vec_reg);

                  for (std::size_t k = j + 1; k < clusters[i].planes.size(); ++ k)
                    {
                      if (clusters[i].coplanar_group[k] == static_cast<std::size_t>(-1))
                        {
                          std::size_t index_prim_next = clusters[i].planes[k];
                          Point pt_reg_next = planes[index_prim_next]->projection(centroids[index_prim_next]);
                          Point pt_proj=plan_reg.projection(pt_reg_next);
                          FT distance = std::sqrt (CGAL::squared_distance(pt_reg_next,pt_proj));

                          if (distance < tolerance_coplanarity)
                            clusters[i].coplanar_group[k] = cop_index;
                        }
                    }
                  cop_index++; 
                }
            }
          //regularize primitive position by computing barycenter of cplanar planes
          std::vector<Point> pt_bary (cop_index, Point ((FT)0., (FT)0., (FT)0.));
          std::vector<FT> area (cop_index, 0.);
      
          for (std::size_t j = 0; j < clusters[i].planes.size (); ++ j)
            {
              std::size_t index_prim = clusters[i].planes[j];
              std::size_t group = clusters[i].coplanar_group[j];
              
              Point pt_reg = planes[index_prim]->projection(centroids[index_prim]);

              pt_bary[group] = CGAL::barycenter (pt_bary[group], area[group], pt_reg, areas[index_prim]); 
              area[group] += areas[index_prim];
            }


          for (std::size_t j = 0; j < clusters[i].planes.size (); ++ j)
            {
              std::size_t index_prim = clusters[i].planes[j];
              std::size_t group = clusters[i].coplanar_group[j];
              Plane plane_reg (pt_bary[group], vec_reg);

              if (planes[index_prim]->plane_normal ()
                  * plane_reg.orthogonal_vector() < 0)
                planes[index_prim]->update (plane_reg.opposite());
              else
                planes[index_prim]->update (plane_reg);
            }
        }
    } 
}


} // namespace CGAL

#endif // CGAL_REGULARIZE_PLANES_H
