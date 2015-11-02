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

  template <typename Traits>
  class Plane_regularization
  {
  public:

    typedef Plane_regularization<Traits> Self;

    typedef typename Traits::FT FT;
    typedef typename Traits::Point_3 Point;
    typedef typename Traits::Vector_3 Vector;
    typedef typename Traits::Plane_3 Plane;

    typedef typename Traits::Point_map Point_map;
    typedef typename Traits::Normal_map Normal_map;
    typedef typename Traits::Input_range Input_range;
    typedef typename Input_range::iterator Input_iterator;

    typedef Shape_detection_3::Shape_base<Traits> Shape;
    typedef Shape_detection_3::Plane<Traits> Plane_shape;
  

  private:

    struct Plane_cluster
    {
      std::vector<std::size_t> planes;
      Vector normal;
      FT cosangle_vertical;
      FT area;
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

    Plane_regularization (Traits t = Traits ())
      : m_traits (t)
    {

    }

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

    virtual ~Plane_regularization ()
    {
      clear ();
    }

    void clear ()
    {
      std::vector<boost::shared_ptr<Plane_shape> > ().swap (m_planes);
      std::vector<Point> ().swap (m_centroids);
      std::vector<FT> ().swap (m_areas);
    }


  
    std::size_t run (FT epsilon, FT tolerance_coplanarity)
    {
      compute_centroids_and_areas ();
      
      // find pairs of epsilon-parallel primitives and store them in parallel_planes
      
      std::vector < std::vector < std::size_t > > parallel_planes (m_planes.size ());
      
      for (std::size_t i = 0; i < m_planes.size (); ++ i)
        {
          Vector v1 = m_planes[i]->plane_normal ();
          
          for (std::size_t j = 0; j < m_planes.size(); ++ j)
            {
              if (i == j)
                continue;
              
              Vector v2 = m_planes[i]->plane_normal ();
              
              if (std::fabs (v1 * v2) > 1. - epsilon)
                parallel_planes[i].push_back (j);
            }
        }
	


      // clustering the parallel primitives and store them in clusters
      // & compute the normal, size and cos angle to the vertical of each cluster
      
      std::vector<Plane_cluster> clusters;
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
                             && std::fabs (normal_it*clu.normal) > 1. - epsilon )
                            {	
                              propagation = true;
                              index_container_current_ring_parallel.push_back(it);
                              is_available[it]=false;
                              
                              if(clu.normal * normal_it <0)
                                normal_it = -normal_it;

                              clu.normal = (FT)clu.area * clu.normal
                                + (FT)m_areas[it] * normal_it;
                              FT norm = 1. / std::sqrt (clu.normal.squared_length()); 
                              clu.normal = norm * clu.normal;
                              clu.area += m_areas[it];
                            }	
                        }
                    }

                  //update containers
                  index_container_former_ring_parallel.clear();
                  for (std::list < int >::iterator it = index_container_current_ring_parallel.begin();
                       it != index_container_current_ring_parallel.end(); ++it)
                    {
                      index_container_former_ring_parallel.push_back(*it);
                      clu.planes.push_back(*it);
                    }
                  index_container_current_ring_parallel.clear();

                }
              while(propagation);

              Vector v_vertical(0.,0.,1.);
              clu.cosangle_vertical = std::fabs(v_vertical*clu.normal);
            }
        }
      is_available.clear();

      //discovery orthogonal relationship between clusters 
      std::vector < std::vector < bool > > group_planes_orthogonal;
      for (std::size_t i = 0; i < clusters.size(); ++ i)
        {
          std::vector<bool> gp_tmp;
          for (std::size_t j = 0; j < clusters.size(); ++ j)
            gp_tmp.push_back(false);
          group_planes_orthogonal.push_back (gp_tmp);
        }

      for (std::size_t i = 0; i < group_planes_orthogonal.size(); ++ i)
        {
          for (std::size_t j = 0; j < group_planes_orthogonal.size(); ++ j)
            {

              if (i != j && std::fabs (clusters[i].normal * clusters[j].normal) < epsilon)
                {
                  group_planes_orthogonal[i][j]=true; 
                  group_planes_orthogonal[j][i]=true;
                }
            }
        }





      //clustering the vertical cosangle and store their centroids in cosangle_centroids and the centroid index of each cluster in list_cluster_index 
      std::vector < double > cosangle_centroids;
      std::vector < int > list_cluster_index;
      for( std::size_t i = 0; i < clusters.size(); ++ i)
        list_cluster_index.push_back(-1);
      
      int mean_index = 0;
      for (std::size_t i = 0; i < clusters.size(); ++ i)
        {
          if(list_cluster_index[i]<0)
            {
              list_cluster_index[i] = mean_index;
              double mean = clusters[i].area * clusters[i].cosangle_vertical;
              double mean_area = clusters[i].area;
              
              for (std::size_t j = i+1; j < clusters.size(); ++ j)
                {
                  if (list_cluster_index[j] < 0 && std::fabs (clusters[j].cosangle_vertical -
                                                              mean / mean_area) < epsilon)
                    {
                      list_cluster_index[j] = mean_index;
                      mean_area += clusters[j].area;
                      mean += clusters[j].area * clusters[j].cosangle_vertical;
                    }
                }
              ++ mean_index;
              mean /= mean_area;
              cosangle_centroids.push_back (mean);
            }
        }
      
      //desactive Z-verticalité
      for (std::size_t i = 0; i < cosangle_centroids.size(); ++ i)
        {
          if (cosangle_centroids[i] < epsilon)
            cosangle_centroids[i] = 0;
          else if (cosangle_centroids[i] > 1. - epsilon)
            cosangle_centroids[i] = 1;
        }
      for (std::size_t i = 0; i < group_planes_orthogonal.size(); ++ i)
        clusters[i].cosangle_vertical = cosangle_centroids[list_cluster_index[i]];
	
      //display console
      /*
      std::cout<<std::endl<<std::endl<<"clusters of parallel primitives:";
      for (std::size_t i=0; i<list_parallel_planes.size();i++)
        {
          std::cout<<std::endl<<i<<" -> ";
          for (std::size_t j=0; j<list_parallel_planes[i].size();j++)
            std::cout<<list_parallel_planes[i][j]<<"  ";
        }

      std::cout<<std::endl<<std::endl<<"pairs of orthogonal clusters:";
      for (std::size_t i=0; i<group_planes_orthogonal.size();i++)
        {
          std::cout<<std::endl<<i<<" -> ";
          for (std::size_t j=0;j<group_planes_orthogonal[i].size();j++)
            {
              if(group_planes_orthogonal[i][j])
                std::cout<<j<<"  ";
            }
          std::cout<<"     -> "<<list_cluster_cosangle_vertical[i]<<"  -> "<<cosangle_centroids[list_cluster_index[i]];
        }
      */


      //find subgraphs of mutually orthogonal clusters (store index of clusters in subgraph_clusters), and select the cluster of largest area
      std::vector < std::vector < int > > subgraph_clusters;
      std::vector < int > subgraph_clusters_max_area_index;
      std::vector < bool > is_free;
      for (std::size_t i = 0; i < clusters.size(); ++ i)
        is_free.push_back(true);
      for (std::size_t i = 0; i < clusters.size(); ++ i)
        {
          if(is_free[i])
            {
              is_free[i]=false;
              double max_area = clusters[i].area;
              int index_max_area = i;

              //initialization containers
              std::vector < int > index_container;
              index_container.push_back(i);
              std::vector < int > index_container_former_ring;
              index_container_former_ring.push_back(i);
              std::list < int > index_container_current_ring;

              //propagation
              bool propagation=true;
              do
                {
                  propagation=false;

                  //neighbors
                  for (std::size_t k=0;k<index_container_former_ring.size();k++)
                    {

                      int cluster_index=index_container_former_ring[k];

                      for (std::size_t j=0;j<group_planes_orthogonal[cluster_index].size();j++)
                        {
                          if(group_planes_orthogonal[cluster_index][j] && is_free[j])
                            { 	
                              propagation=true;
                              index_container_current_ring.push_back(j);
                              is_free[j]=false;

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
                  for(std::list < int >::iterator it = index_container_current_ring.begin();
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
      is_free.clear();


      //create subgraphs of mutually orthogonal clusters in which the largest cluster is excluded and store in subgraph_clusters_prop
      std::vector < std::vector < int > > subgraph_clusters_prop;
      for (std::size_t i=0;i<subgraph_clusters.size(); i++)
        {
          int index=subgraph_clusters_max_area_index[i];
          std::vector < int > subgraph_clusters_prop_temp;
          for (std::size_t j=0;j<subgraph_clusters[i].size(); j++)
            if(subgraph_clusters[i][j]!=index)
              subgraph_clusters_prop_temp.push_back(subgraph_clusters[i][j]);

          subgraph_clusters_prop.push_back(subgraph_clusters_prop_temp);
        }



      //display console
      /*
        for (std::size_t i=0;i<subgraph_clusters_prop.size(); i++)
        {
        std::cout<<std::endl<<std::endl<<"subgraph "<<i<<" ("<<subgraph_clusters_max_area_index[i]<<" max area): ";
        for (std::size_t j=0;j<subgraph_clusters_prop[i].size(); j++) std::cout<<subgraph_clusters_prop[i][j]<<"  ";
        }
      */


      //regularization of cluster normals : in eachsubgraph, we start from the largest area cluster and we propage over the subgraph by regularizing the normals of the clusters accorting to orthogonality and cosangle to vertical
      std::vector< bool > cluster_is_available; 
      for (std::size_t i = 0; i < clusters.size(); ++ i)
        cluster_is_available.push_back(true);

      for (std::size_t i = 0; i < subgraph_clusters_prop.size(); ++ i)
        {
	
          int index_current=subgraph_clusters_max_area_index[i];
          Vector vec_current=regularize_normal(clusters[index_current].normal,
                                               clusters[index_current].cosangle_vertical);
          clusters[index_current].normal = vec_current;
          cluster_is_available[index_current] = false;

          //initialization containers
          std::vector < int > index_container;
          index_container.push_back(index_current);
          std::vector < int > index_container_former_ring;
          index_container_former_ring.push_back(index_current);
          std::list < int > index_container_current_ring;

          //propagation
          bool propagation=true;
          do
            {
              propagation=false;

              //neighbors
              for (std::size_t k=0;k<index_container_former_ring.size();k++)
                {

                  int cluster_index=index_container_former_ring[k];

                  for (std::size_t j=0;j<group_planes_orthogonal[cluster_index].size();j++)
                    {
						
                      if(group_planes_orthogonal[cluster_index][j] && cluster_is_available[j])
                        { 	
							
                          propagation=true;
                          index_container_current_ring.push_back(j);
                          cluster_is_available[j]=false;

                          Vector new_vect=regularize_normals_from_prior(clusters[cluster_index].normal,
                                                                        clusters[j].normal,
                                                                        clusters[j].cosangle_vertical);
                          clusters[j].normal = new_vect;
                        }
                    }	
                }
			
              //update containers
              index_container_former_ring.clear();
              for(std::list < int >::iterator it = index_container_current_ring.begin();
                  it != index_container_current_ring.end(); ++it)
                {
                  index_container_former_ring.push_back(*it);
                  index_container.push_back(*it);
                }
              index_container_current_ring.clear();
            }while(propagation);
        }



      //recompute optimal plane for each primitive after normal regularization
      for (std::size_t i=0; i < clusters.size(); ++ i)
        {

          Vector vec_reg = clusters[i].normal;

          for (std::size_t j = 0; j < clusters[i].planes.size(); ++ j)
            {
              int index_prim = clusters[i].planes[j];
              Point pt_reg = m_planes[index_prim]->projection (m_centroids[index_prim]);
              if( m_planes[index_prim]->plane_normal () * vec_reg < 0)
                vec_reg=-vec_reg;
              Plane plane_reg(pt_reg,vec_reg);
		
              if( std::fabs(m_planes[index_prim]->plane_normal () * plane_reg.orthogonal_vector ()) > 1. - epsilon)
                {
                  m_planes[index_prim]->update (plane_reg);
                }
            }
        }





      //detecting co-planarity and store in list_coplanar_prim
      std::vector< std::vector< std::vector < int > > > list_coplanar_prim;
      for (std::size_t i = 0; i < clusters.size(); ++ i)
        {

          std::vector< std::vector < int > > list_coplanar_prim_tmp;
          Vector vec_reg = clusters[i].normal;
          std::vector < int > list_cop_index;
          for (std::size_t ip = 0; ip < clusters[i].planes.size(); ++ ip)
            list_cop_index.push_back(-1);

          int cop_index=0;

          for (std::size_t j = 0; j < clusters[i].planes.size(); ++ j)
            {
              int index_prim = clusters[i].planes[j];

              if (list_cop_index[j] < 0)
                {
		
                  std::vector < int > list_coplanar_prim_tmp_tmp;
                  list_cop_index[j]=cop_index;
                  list_coplanar_prim_tmp_tmp.push_back(index_prim);
			
                  Point pt_reg = m_planes[index_prim]->projection(m_centroids[index_prim]);
                  Plane plan_reg(pt_reg,vec_reg);

                  for (std::size_t k = j+1; k < clusters[i].planes.size(); ++ k)
                    {
                      if (list_cop_index[k] < 0)
                        {
                          int index_prim_next = clusters[i].planes[k];
                          Point pt_reg_next = m_planes[index_prim_next]->projection(m_centroids[index_prim_next]);
                          Point pt_proj=plan_reg.projection(pt_reg_next);
                          double distance=distance_Point(pt_reg_next,pt_proj);
					
                          if(distance<tolerance_coplanarity )
                            {
                              list_cop_index[k]=cop_index;
                              list_coplanar_prim_tmp_tmp.push_back(index_prim_next);
                            }
                        }
                    }
                  list_coplanar_prim_tmp.push_back(list_coplanar_prim_tmp_tmp);
                  cop_index++; 
                }
            }
          list_coplanar_prim.push_back(list_coplanar_prim_tmp);
        }



      //regularize primitive position by computing barycenter of coplanar planes
      std::vector < std::vector < int > > list_primitive_reg_index_extracted_planes;
      std::vector < Plane > list_primitive_reg;

      for (std::size_t i=0;i<list_coplanar_prim.size();i++)
        {
          for (std::size_t j=0;j<list_coplanar_prim[i].size();j++)
            {

              Point pt_bary(0.,0.,0.);
              double area=0;

              for (std::size_t k=0; k<list_coplanar_prim[i][j].size();k++)
                {
                  int index_prim=list_coplanar_prim[i][j][k];
                  Point pt_reg = m_planes[index_prim]->projection(m_centroids[index_prim]);

                  pt_bary=barycenter(pt_bary, area,pt_reg,m_areas[index_prim]); 
                  area+=m_areas[index_prim];
                }
              Vector vec_reg = m_planes[list_coplanar_prim[i][j][0]]->plane_normal ();

              Plane plane_reg(pt_bary,vec_reg);

              bool is_reg_used=false;
              std::vector< int > list_primitive_reg_index_extracted_planes_tmp1;

              for (std::size_t k=0; k<list_coplanar_prim[i][j].size();k++)
                {
                  int index_prim=list_coplanar_prim[i][j][k];
                  if( std::fabs(m_planes[index_prim]->plane_normal () * plane_reg.orthogonal_vector()) > 1. - epsilon)
                    {
                      if(m_planes[index_prim]->plane_normal () * plane_reg.orthogonal_vector()<0)
                        m_planes[index_prim]->update (plane_reg.opposite());
                      else
                        m_planes[index_prim]->update (plane_reg);
                      is_reg_used=true;
                      list_primitive_reg_index_extracted_planes_tmp1.push_back(index_prim);
                    }
                  else{
                    list_primitive_reg.push_back(static_cast<Plane> (*(m_planes[index_prim])));
                    std::vector< int > list_primitive_reg_index_extracted_planes_tmp;
                    list_primitive_reg_index_extracted_planes_tmp.push_back(index_prim);
                    list_primitive_reg_index_extracted_planes.push_back(list_primitive_reg_index_extracted_planes_tmp);
                  }
                }
              if(is_reg_used) {
                list_primitive_reg.push_back(plane_reg);
                list_primitive_reg_index_extracted_planes.push_back(list_primitive_reg_index_extracted_planes_tmp1);
              }
            }
        }

      //      std::cout<<std::endl<<std::endl<<"NB planes final = "<<list_primitive_reg.size()<<std::endl<<std::endl;



      //merge similar planes in plane_point_index and extracted planes and HPS[i].primitive_index
      // std::vector < std::vector < int > > plane_point_index_temp;
      // std::vector < Plane > extracted_planes_temp;
      // std::vector < bool > has_been_merged;
      // for (std::size_t i=0; i< m_planes.size();i++)
      //   has_been_merged.push_back(false);

      // for (std::size_t i=0; i< m_planes.size();i++)
      //   {
	
      //     if (!has_been_merged[i])
      //       {
      //         extracted_planes_temp.push_back (m_planes[i]);
      //         int label_index=extracted_planes_temp.size()-1;
      //         plane_point_index_temp.push_back(m_planes[i]->indices_of_assigned_points[i]);
      //         for (std::size_t k=0; k< m_planes[i]->indices_of_assigned_points().size();k++)
      //           {
      //             int index_pt=m_planes[i]->indices_of_assigned_points[k];
      //             primitive_index[index_pt]=label_index;
      //             label_plane[index_pt]=label_index;
      //           }
	
      //         for (std::size_t j=i+1;j< m_planes.size();j++)
      //           {

      //             if(m_planes[i]==m_planes[j])
      //               { //if identical (do opposite plane too ?) then store the second in the first
			
      //                 has_been_merged[j]=true;

      //                 std::vector< int > plane_point_index_new
      //                   = plane_point_index_temp[m_planes.size()-1];
      //                 for (std::size_t k=0; k< m_planes[j]->indices_of_assigned_points().size();k++)
      //                   {
      //                     int ind=m_planes[j]->indices_of_assigned_points[k];
      //                     plane_point_index_new.push_back(ind);
      //                     primitive_index[ind]=label_index;
      //                     label_plane[ind]=label_index;
      //                   }
      //                 plane_point_index_temp[plane_point_index_temp.size()-1]=plane_point_index_new;

      //               }
      //           }
      //       }
      //   }

      // TODO
      // m_planes=extracted_planes_temp;
      // m_planes->indices_of_assigned_points=plane_point_index_temp;

      return list_primitive_reg.size ();
    }

    void compute_centroids_and_areas ()
    {
      for (std::size_t i = 0; i < m_planes.size (); ++ i)
        {
          std::vector < Point > listp;
          for (std::size_t j = 0; j < m_planes[i]->indices_of_assigned_points ().size (); ++ j)
            {
              int yy = m_planes[i]->indices_of_assigned_points()[j];
              Point pt = get (m_point_pmap, *(m_input_begin + yy));
              listp.push_back(pt);
            }
          m_centroids.push_back (CGAL::centroid (listp.begin (), listp.end ()));
          m_areas.push_back ((double)(m_planes[i]->indices_of_assigned_points().size()) / 100.);
        }
    }
    
    FT distance_Point (const Point& a, const Point& b)
    {
      return std::sqrt (CGAL::squared_distance (a, b));
    }
  
    Vector regularize_normal (const Vector& n, FT cos_vertical)
    {
      FT A = 1 - cos_vertical * cos_vertical;
      FT B = 1 + (n.y() * n.y()) / (n.x() * n.x());
    
      FT vx = std::sqrt (A/B);
    
      if (n.x() < 0)
        vx = -vx;
    
      FT vy = vx * (n.y() / n.x()); 

      Vector res (vx, vy, cos_vertical);

      return res / std::sqrt (res * res);
    }

  
    Vector regularize_normals_from_prior (const Vector& np,
                                          const Vector& n,
                                          FT cos_vertical)
    {
      FT vx, vy;

      if (np.x() != 0)
        { 
          FT a = (np.y() * np.y()) / (np.x() * np.x()) + 1;
          FT b = 2 * np.y() * np.z() * cos_vertical / np.x();
          FT c= cos_vertical * cos_vertical-1;

          if (4 * a * c > b * b)
            return regularize_normal (n, cos_vertical); 
          else
            {
              FT delta = std::sqrt (b * b-4 * a * c);
              FT vy1= (-b-delta) / (2 * a);
              FT vy2= (-b+delta) / (2 * a);

              vy = (std::fabs(n.y()-vy1) < std::fabs(n.y()-vy2))
                ? vy1 : vy2;

              vx = (-np.y() * vy-np.z() * cos_vertical) / np.x();
            }
        }
      else if (np.y() != 0)
        {
          vy = -np.z() * cos_vertical / np.y();
          vx = std::sqrt (1 - cos_vertical * cos_vertical - vy * vy);
        
          if (n.x() < 0)
            vx = -vx;
        }
      else
        return regularize_normal (n, cos_vertical); 

      Vector res (vx, vy, cos_vertical);
      FT norm = std::max(1e-5, 1. / sqrt(res.squared_length ()));

      return norm * res;
    }

  

  };


}; // namespace CGAL

#endif // CGAL_PLANE_REGULARIZATION_H
