// Copyright (c) 2015 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Florent Lafarge, Simon Giraudot
//

#ifndef CGAL_REGULARIZE_PLANES_H
#define CGAL_REGULARIZE_PLANES_H

#include <CGAL/license/Shape_detection.h>

#include <CGAL/centroid.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Shape_detection/Efficient_RANSAC.h>

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

  Plane_cluster()
    : is_free (true)
    , normal (typename Traits::FT(0.),
              typename Traits::FT(0.),
              typename Traits::FT(1.))
    , cosangle_symmetry (typename Traits::FT(0.))
    , area (typename Traits::FT(0.))
    , cosangle_centroid (typename Traits::FT(0.))
  { }
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
          typename PointRange,
          typename PointMap,
          typename IndexMap>
void compute_centroids_and_areas (const PointRange& points,
                                  PointMap point_map,
                                  std::size_t nb_planes,
                                  IndexMap index_map,
                                  std::vector<typename Traits::Point_3>& centroids,
                                  std::vector<typename Traits::FT>& areas)
{
  typedef typename Traits::FT FT;
  typedef typename Traits::Point_3 Point;

  std::vector<std::vector<Point> > listp (nb_planes);

  for (std::size_t i = 0; i < points.size(); ++ i)
  {
    int idx = get (index_map, i);
    if (idx != -1)
      listp[std::size_t(idx)].push_back (get(point_map, *(points.begin() + i)));
  }

  centroids.reserve (nb_planes);
  areas.reserve (nb_planes);
  for (std::size_t i = 0; i < nb_planes; ++ i)
  {
    centroids.push_back (CGAL::centroid (listp[i].begin (), listp[i].end ()));
    areas.push_back ((FT)(listp[i].size() / (FT)100.));
  }
}


template <typename Traits,
          typename PlaneRange,
          typename PlaneMap>
void compute_parallel_clusters (PlaneRange& planes,
                                PlaneMap plane_map,
                                std::vector<Plane_cluster<Traits> >& clusters,
                                std::vector<typename Traits::FT>& areas,
                                typename Traits::FT tolerance_cosangle,
                                const typename Traits::Vector_3& symmetry_direction)
{

  typedef typename Traits::FT FT;
  typedef typename Traits::Vector_3 Vector;

  // find pairs of epsilon-parallel primitives and store them in parallel_planes
  std::vector<std::vector<std::size_t> > parallel_planes (planes.size ());
  for (std::size_t i = 0; i < std::size_t(planes.size ()); ++ i)
    {
      typename PlaneRange::iterator it = planes.begin() + i;
      Vector v1 = get(plane_map, *it).orthogonal_vector();

      for (std::size_t j = 0; j < std::size_t(planes.size()); ++ j)
        {
          if (i == j)
            continue;

          typename PlaneRange::iterator it2 = planes.begin() + j;
          Vector v2 = get(plane_map, *it2).orthogonal_vector();

          if (std::fabs (v1 * v2) > 1. - tolerance_cosangle)
            parallel_planes[i].push_back (j);
        }
    }


  std::vector<bool> is_available (planes.size (), true);

  for (std::size_t i = 0; i < std::size_t(planes.size()); ++ i)
    {

      if(is_available[i])
        {
          const typename Traits::Plane_3& plane = get(plane_map, *(planes.begin() + i));

          is_available[i] = false;

          clusters.push_back (Plane_cluster<Traits>());
          Plane_cluster<Traits>& clu = clusters.back ();

          //initialization containers
          clu.planes.push_back (i);

          std::vector<std::size_t> index_container_former_ring_parallel;
          index_container_former_ring_parallel.push_back(i);

          std::list<std::size_t> index_container_current_ring_parallel;

          //propagation over the pairs of epsilon-parallel primitives
          bool propagation=true;
          clu.normal = plane.orthogonal_vector ();
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

                      Vector normal_it = get(plane_map, *(planes.begin() + it)).orthogonal_vector ();

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

template <typename Traits>
void cluster_symmetric_cosangles (std::vector<Plane_cluster<Traits> >& clusters,
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


template <typename Traits>
void subgraph_mutually_orthogonal_clusters (std::vector<Plane_cluster<Traits> >& clusters,
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

  //regularization of cluster normals : in each subgraph, we start
  //from the largest area cluster and we propagate over the subgraph
  //by regularizing the normals of the clusters according to
  //orthogonality and cos angle to symmetry direction

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

/// \ingroup PkgShapeDetectionRef

  /*!
    Given a set of detected planes with their corresponding inlier sets,
    this function enables to regularize the planes:

    - Planes near parallel can be made exactly parallel;

    - Planes near orthogonal can be made exactly orthogonal;

    - Planes parallel and near coplanar can be made exactly coplanar;

    - Planes near symmetrical with a user-defined axis can be made
    exactly symmetrical.

    Planes are directly modified. Points are left unaltered, as well as
    their relationships to planes (no transfer of point from a primitive
    plane to another).

    The implementation follows \cgalCite{cgal:vla-lod-15}.

    \tparam PointRange must be a model of `ConstRange` with points.

    \tparam PointPMap must be a model of `ReadablePropertyMap` with value type `Kernel::Point_3`.
    It can be omitted if the value type of the iterator of `PointRange` is convertible to `Point_3<Kernel>`.

    \tparam PlaneRange must be a model of `Range` with planes.

    \tparam PlaneMap must be a model of `WritablePropertyMap` with value type `Kernel::Plane_3`.
    It can be omitted if the value type of the iterator of `PlaneRange` is convertible to `Plane_3<Kernel>`.

    \tparam IndexMap must be a model of `ReadablePropertyMap` with value type `int`.

    \tparam Kernel must be a geometric traits class.
    It can be omitted and deduced automatically from the value type of `PointMap`.

    \param points `ConstRange` of points
    \param point_map property map: value_type of `typename PointRange::const_iterator` -> `Point_3`
    \param planes range of planes
    \param plane_map property map: value_type of `typename PlaneRange::iterator` -> `Plane_3`
    \param index_map property map: index of point `std::size_t` -> index of plane `int` (-1 if the point is not assigned to a plane)

    \param regularize_parallelism select whether parallelism is
    regularized or not

    \param regularize_orthogonality select whether orthogonality is
    regularized or not

    \param regularize_coplanarity select whether coplanarity is
    regularized or not

    \param regularize_axis_symmetry select whether axis symmetry is
    regularized or not

    \param tolerance_angle tolerance of deviation between normal
    vectors of planes (in degrees) used for parallelism, orthogonality,
    and axis symmetry. %Default value is 25 degrees.

    \param tolerance_coplanarity maximal distance between two parallel
    planes such that they are considered coplanar. %Default value is
    0.01.

    \param symmetry_direction chosen axis for symmetry
    regularization. %Default value is the Z axis.
*/

// This variant requires all parameters
template <typename PointRange,
          typename PointMap,
          typename PlaneRange,
          typename PlaneMap,
          typename IndexMap,
          typename Kernel>
void regularize_planes(const PointRange& points,
                        PointMap point_map,
                        PlaneRange& planes,
                        PlaneMap plane_map,
                        IndexMap index_map,
                        const Kernel&,
                        bool regularize_parallelism,
                        bool regularize_orthogonality,
                        bool regularize_coplanarity,
                        bool regularize_axis_symmetry,
                        double tolerance_angle = 25.0,
                        double tolerance_coplanarity = 0.01,
                        typename Kernel::Vector_3 symmetry_direction
                        = typename Kernel::Vector_3 (0.0, 0.0, 1.0))
{
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Kernel::Plane_3 Plane;

  typedef typename internal::PlaneRegularization::Plane_cluster<Kernel>
    Plane_cluster;

  /*
   * Compute centroids and areas
   */
  std::vector<Point> centroids;
  std::vector<FT> areas;
  internal::PlaneRegularization::compute_centroids_and_areas<Kernel>
    (points, point_map, planes.size(), index_map, centroids, areas);

  tolerance_angle = tolerance_angle * (FT)CGAL_PI / (FT)(180);
  FT tolerance_cosangle = (FT)(1. - std::cos (tolerance_angle));
  FT tolerance_cosangle_ortho = (FT)(std::cos ((FT)0.5 * (FT)CGAL_PI - (FT)tolerance_angle));

  // clustering the parallel primitives and store them in clusters
  // & compute the normal, size and cos angle to the symmetry
  // direction of each cluster
  std::vector<Plane_cluster> clusters;
  internal::PlaneRegularization::compute_parallel_clusters<Kernel>
    (planes, plane_map, clusters, areas,
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
      //clustering the symmetry cos angle and store their centroids in
      //cosangle_centroids and the centroid index of each cluster in
      //list_cluster_index
      internal::PlaneRegularization::cluster_symmetric_cosangles<Kernel>
        (clusters, tolerance_cosangle, tolerance_cosangle_ortho);
    }

  //find subgraphs of mutually orthogonal clusters (store index of
  //clusters in subgraph_clusters), and select the cluster of
  //largest area
  if (regularize_orthogonality || regularize_axis_symmetry)
    internal::PlaneRegularization::subgraph_mutually_orthogonal_clusters<Kernel>
      (clusters, (regularize_axis_symmetry ? symmetry_direction : CGAL::NULL_VECTOR));

  //recompute optimal plane for each primitive after normal regularization
  for (std::size_t i=0; i < clusters.size(); ++ i)
    {
      Vector vec_reg = clusters[i].normal;
      for (std::size_t j = 0; j < clusters[i].planes.size(); ++ j)
        {
          std::size_t index_prim = clusters[i].planes[j];
          const Plane& plane = get(plane_map, *(planes.begin() + index_prim));

          Point pt_reg = plane.projection (centroids[index_prim]);
          if(plane.orthogonal_vector () * vec_reg < 0)
            vec_reg=-vec_reg;
          Plane plane_reg(pt_reg,vec_reg);

          if( std::fabs(plane.orthogonal_vector () * vec_reg) > 1. - tolerance_cosangle)
            put(plane_map, *(planes.begin() + index_prim), plane_reg);
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
                  const Plane& plane = get(plane_map, *(planes.begin() + index_prim));

                  clusters[i].coplanar_group[j] = cop_index;

                  Point pt_reg = plane.projection(centroids[index_prim]);
                  Plane plan_reg(pt_reg,vec_reg);

                  for (std::size_t k = j + 1; k < clusters[i].planes.size(); ++ k)
                    {
                      if (clusters[i].coplanar_group[k] == static_cast<std::size_t>(-1))
                        {
                          std::size_t index_prim_next = clusters[i].planes[k];
                          const Plane& plane_next = get(plane_map, *(planes.begin() + index_prim_next));
                          Point pt_reg_next = plane_next.projection(centroids[index_prim_next]);
                          Point pt_proj=plan_reg.projection(pt_reg_next);
                          FT distance = std::sqrt (CGAL::squared_distance(pt_reg_next,pt_proj));

                          if (distance < tolerance_coplanarity)
                            clusters[i].coplanar_group[k] = cop_index;
                        }
                    }
                  cop_index++;
                }
            }
          //regularize primitive position by computing barycenter of coplanar planes
          std::vector<Point> pt_bary (cop_index, Point ((FT)0., (FT)0., (FT)0.));
          std::vector<FT> area (cop_index, 0.);

          for (std::size_t j = 0; j < clusters[i].planes.size (); ++ j)
            {
              std::size_t index_prim = clusters[i].planes[j];
              std::size_t group = clusters[i].coplanar_group[j];

              Point pt_reg = get(plane_map, *(planes.begin()+index_prim)).projection(centroids[index_prim]);

              pt_bary[group] = CGAL::barycenter (pt_bary[group], area[group], pt_reg, areas[index_prim]);
              area[group] += areas[index_prim];
            }


          for (std::size_t j = 0; j < clusters[i].planes.size (); ++ j)
            {
              std::size_t index_prim = clusters[i].planes[j];
              std::size_t group = clusters[i].coplanar_group[j];
              Plane plane_reg (pt_bary[group], vec_reg);

              if (get(plane_map, *(planes.begin() + index_prim)).orthogonal_vector ()
                  * plane_reg.orthogonal_vector() < 0)
                put(plane_map, *(planes.begin() + index_prim), plane_reg.opposite());
              else
                put(plane_map, *(planes.begin() + index_prim), plane_reg);
            }
        }
    }
}


/// \cond SKIP_IN_MANUAL

// Workaround for bug reported here:
// https://developercommunity.visualstudio.com/content/problem/340310/unaccepted-typename-that-other-compilers-require.html
#if _MSC_VER == 1915
#define CGAL_TYPENAME_FOR_MSC
#else
#define CGAL_TYPENAME_FOR_MSC typename
#endif

// This variant deduces the kernel from the point property map.
template <typename PointRange,
          typename PointMap,
          typename PlaneRange,
          typename PlaneMap,
          typename IndexMap>
void regularize_planes (const PointRange& points,
                        PointMap point_map,
                        PlaneRange& planes,
                        PlaneMap plane_map,
                        IndexMap index_map,
                        bool regularize_parallelism,
                        bool regularize_orthogonality,
                        bool regularize_coplanarity,
                        bool regularize_axis_symmetry,
                        double tolerance_angle = 25.0,
                        double tolerance_coplanarity = 0.01,
                        typename Kernel_traits
                        <typename boost::property_traits
                        <PointMap>::value_type>::Kernel::Vector_3 symmetry_direction
                        = CGAL_TYPENAME_FOR_MSC Kernel_traits
                          <typename boost::property_traits
                          <PointMap>::value_type>::Kernel::Vector_3
                        (
                         CGAL_TYPENAME_FOR_MSC Kernel_traits
                         <typename boost::property_traits
                         <PointMap>::value_type>::Kernel::FT(0.),
                         CGAL_TYPENAME_FOR_MSC Kernel_traits
                         <typename boost::property_traits
                         <PointMap>::value_type>::Kernel::FT(0.),
                         CGAL_TYPENAME_FOR_MSC Kernel_traits
                         <typename boost::property_traits
                         <PointMap>::value_type>::Kernel::FT(1.)))
{
  typedef typename boost::property_traits<PointMap>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel Kernel;

  regularize_planes (points, point_map, planes, plane_map, index_map, Kernel(),
                     regularize_parallelism, regularize_orthogonality,
                     regularize_coplanarity, regularize_axis_symmetry,
                     tolerance_angle, tolerance_coplanarity, symmetry_direction);
}

#ifdef CGAL_TYPENAME_FOR_MSC
#undef CGAL_TYPENAME_FOR_MSC
#endif


/// \endcond

} // namespace CGAL

#endif // CGAL_REGULARIZE_PLANES_H
