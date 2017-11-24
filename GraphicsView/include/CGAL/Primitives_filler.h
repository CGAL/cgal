#ifndef PRIMITIVES_FILLER_H
#define PRIMITIVES_FILLER_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_2_projection_traits_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Cartesian_converter.h>
#include <boost/unordered_map.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Local_kernel;
typedef Local_kernel::Point_3  Local_point;
typedef Local_kernel::Vector_3 Local_vector;
struct Vertex_info
{
  Local_vector v;
};

struct Face_info
{
  bool exist_edge[3];
  bool is_external;
  bool is_process;
};

typedef CGAL::Triangulation_2_projection_traits_3<CGAL::Exact_predicates_inexact_constructions_kernel> P_traits;
typedef CGAL::Triangulation_2_projection_traits_3<CGAL::Exact_predicates_inexact_constructions_kernel> Local_traits;
typedef CGAL::Triangulation_vertex_base_with_info_2<Vertex_info, P_traits> Vb;

typedef CGAL::Triangulation_face_base_with_info_2<Face_info, P_traits> Fb1;

typedef CGAL::Constrained_triangulation_face_base_2<P_traits, Fb1>    Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>                   TDS;
typedef CGAL::Exact_predicates_tag                                    Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<P_traits, TDS,
Itag>              CDT;

namespace CGAL {
namespace internal {

template<int dim>
struct Geom_utils;

template<>
struct Geom_utils<3>
{
  template<typename KPoint>
  static Local_point get_local_point(const KPoint& p)
  {
    CGAL::Cartesian_converter<typename CGAL::Kernel_traits<KPoint>::Kernel, Local_kernel> converter;
    return converter(p);
  }

  template<typename KVector>
  static Local_vector get_local_vector(const KVector& v)
  {
    CGAL::Cartesian_converter<typename CGAL::Kernel_traits<KVector>::Kernel, Local_kernel> converter;
    return converter(v);
  }
};

template<>
struct Geom_utils<2>
{
  template<typename KPoint>
  static Local_point get_local_point(const KPoint& p)
  {
    CGAL::Cartesian_converter<typename CGAL::Kernel_traits<KPoint>::Kernel, Local_kernel> converter;
    return Local_point(converter(p.x()),0,converter(p.y()));
  }

  template<typename KVector>
  static Local_vector get_local_vector(const KVector& v)
  {
    CGAL::Cartesian_converter<typename CGAL::Kernel_traits<KVector>::Kernel, Local_kernel> converter;
    return Local_vector(converter(v.x()),0,converter(v.y()));
  }
};

template<typename KPoint>
Local_point get_local_point(const KPoint& p)
{
  return Geom_utils<CGAL::Ambient_dimension<KPoint>::value>::get_local_point(p);
}

template<typename KVector>
Local_vector get_local_vector(const KVector& v)
{
  return Geom_utils<CGAL::Ambient_dimension<KVector>::value>::get_local_vector(v);
}

} // End namespace internal

class Primitive_filler
{
public:
  Primitive_filler()
    : segments_vector(NULL),
      segment_color_vector(NULL),
      faces_vector(NULL),
      flat_normal_vector(NULL),
      smooth_normal_vector(NULL),
      face_color_vector(NULL),
      idx_segments_vector(NULL),
      idx_faces_vector(NULL)
  {
    m_face_started=false;
    m_data_is_indexed = false;
  }

  Primitive_filler(std::vector<float>* p_segments_vector,
                   std::vector<float>* p_segment_color_vector,
                   std::vector<float>* p_faces_vector,
                   std::vector<float>* p_flat_normal_vector,
                   std::vector<float>* p_smooth_normal_vector,
                   std::vector<float>* p_face_color_vector)
    : segments_vector(p_segments_vector),
      segment_color_vector(p_segment_color_vector),
      faces_vector(p_faces_vector),
      flat_normal_vector(p_flat_normal_vector),
      smooth_normal_vector(p_smooth_normal_vector),
      face_color_vector(p_face_color_vector),
      idx_segments_vector(NULL),
      idx_faces_vector(NULL)

  {
    m_face_started = false;
    m_data_is_indexed = false;
  }

  Primitive_filler(std::vector<unsigned int>* p_idx_segments_vector,
                   std::vector<unsigned int>* p_idx_faces_vector,
                   std::vector<float>* p_smooth_normal_vector,
                   std::vector<float>* p_vertex_color_vector)
    : segments_vector(NULL),
      segment_color_vector(NULL),
      faces_vector(NULL),
      flat_normal_vector(NULL),
      smooth_normal_vector(p_smooth_normal_vector),
      face_color_vector(p_vertex_color_vector),
      idx_segments_vector(p_idx_segments_vector),
      idx_faces_vector(p_idx_faces_vector)
  {
    m_face_started = false;
    m_data_is_indexed = true;
  }

  template<typename KPoint, typename KVector>
  void add_indexed_point(const KPoint& kp,
                         const unsigned int& index,
                         const KVector& p_normal)
  {
    Local_point p = internal::get_local_point(kp);
    points_vector.push_back(p);
    point_index_map[p] = index;
    add_normal(internal::get_local_vector(p_normal), *smooth_normal_vector);
  }

  template<typename KPoint, typename KVector>
  void add_indexed_point(const KPoint& kp,
                         const unsigned int& index,
                         const KVector& p_normal,
                         const CGAL::Color& color)
  {
    add_indexed_point(kp, index, p_normal);
    add_color(color, *face_color_vector);
  }

  template<typename KPoint>
  void add_mono_segment(const KPoint& p1, const KPoint& p2)
  {
    if(segments_vector)
    {
      add_point(p1, *segments_vector);
      add_point(p2, *segments_vector);
    }
    else if(idx_segments_vector)
    {
      idx_segments_vector->push_back(point_index_map[p1]);
      idx_segments_vector->push_back(point_index_map[p2]);
    }
  }


   template<typename KPoint>
   void add_colored_segment(const KPoint& p1, const KPoint& p2,
                            const CGAL::Color& acolor)
   {
     if(segments_vector)
     {
       add_point(p1, *segments_vector);
       add_point(p2, *segments_vector);
     }
     else if(idx_segments_vector)
     {
       idx_segments_vector->push_back(point_index_map[p1]);
       idx_segments_vector->push_back(point_index_map[p2]);
     }
     if(segment_color_vector)
     {
       add_color(acolor, *segment_color_vector);
       add_color(acolor, *segment_color_vector);
     }
   }


   static bool is_facet_convex(const std::vector<Local_point>& facet,
                        const Local_vector& N)
   {

     typedef Local_kernel::Orientation Orientation;
     Orientation orientation;
     Local_vector normal = N;
     bool normal_is_ok;
     std::size_t id = 0;
     do{
       normal_is_ok = true;
       //Initializes the facet orientation

       Local_point S,T;
       S = facet[id];
       T = facet[id+1];
       Local_vector V1 = Local_vector((T-S).x(), (T-S).y(), (T-S).z());
       S = T;
       T = facet[id+2];
       Local_vector V2 = Local_vector((T-S).x(), (T-S).y(), (T-S).z());

       if(normal == Local_vector(0,0,0))
         normal_is_ok = false;
       if(normal_is_ok)
       {
         orientation = Local_kernel::Orientation_3()(V1, V2, normal);
         if( orientation == CGAL::COPLANAR )
           normal_is_ok = false;
       }
       //Checks if the normal is good : if the normal is null
       // or if it is coplanar to the facet, we need another one.
       if(!normal_is_ok)
       {
         normal = CGAL::cross_product(V1,V2);
       }

     }while( ++id != facet.size() && !normal_is_ok);
     //if no good normal can be found, stop here.
     if (!normal_is_ok)
       return false;

     //computes convexness

     //re-initializes he_end;

     for(id=0; id<facet.size(); ++id)
     {
       Local_point S,T;
       S = facet[id%facet.size()];
       T = facet[(id+1)%facet.size()];
       Local_vector V1 = Local_vector((T-S).x(), (T-S).y(), (T-S).z());
       S = T;
       T = facet[(id+2)%facet.size()];
       Local_vector V2 = Local_vector((T-S).x(), (T-S).y(), (T-S).z());
       Orientation res = Local_kernel::Orientation_3()(V1, V2, normal) ;

       if(res!= orientation && res != CGAL::ZERO)
         return false;
     }
     return true;
   }

  template<typename KVector>
  void face_begin(const KVector& p_normal)
  {
    if (m_face_started)
    {
      std::cerr<<"You cannot start a new face before finishing the previous one."<<std::endl;
      return;
    }

    normal_for_face = internal::get_local_vector(p_normal);
    m_face_started=true;
  }
  template<typename KVector>
  void mono_face_begin(const KVector& p_normal)
  {
    m_started_face_is_colored=false;
    face_begin(p_normal);
  }

  /// Start a new face, with a given color.
  template<typename KVector>
  void colored_face_begin(const CGAL::Color& acolor,
                          const KVector& p_normal)
  {
    if(m_data_is_indexed)
    {
      std::cerr<<"You cannot have per face color with indexed data."<<std::endl;
      return;
    }
    color_of_face=acolor;
    m_started_face_is_colored=true;
    face_begin(p_normal);
  }

  /// Add a point at the end of the current face
  /// With this method, it is not possible to use the Gourod shading.
  /// @param p the point to add
  template<typename KPoint>
  bool add_point_in_face(const KPoint& kp)
  {
    if (!m_face_started) return false;

    Local_point p=internal::get_local_point(kp);
    if (points_of_face.empty() || points_of_face.back()!=p)
    {
      points_of_face.push_back(p);
      return true;
    }
    return false;
  }

  /// Add a point at the end of the current face
  /// @param p the point to add
  /// @p_normal the vertex normal in this point (for Gourod shading)
  template<typename KPoint, typename KVector>
  void add_point_in_face(const KPoint& p, const KVector& p_normal)
  {
    if(m_data_is_indexed)
    {
      std::cerr<<"You cannot have per face normals with indexed data."<<std::endl;
      return;
    }
    if (!m_face_started) return;

    if (add_point_in_face(p))
    {
      vertex_normals_for_face.push_back(internal::get_local_vector(p_normal));
    }
  }

  /// End the face: compute the triangulation.
  void face_end()
  {
    if (points_of_face.size()<3)
    {
      std::cerr<<"ERROR: trying to triangulate a face with "<<points_of_face.size()<<" vertices."
              <<std::endl;

      m_face_started=false;
      points_of_face.clear();
      vertex_normals_for_face.clear();

      return;
    }


    if (points_of_face.size()==3) // Triangle: no need to triangulate
    {
      for (int i=0; i<3; ++i)
      {
        // The point
        if(faces_vector)
          add_point(points_of_face[i], *faces_vector);
        else if(idx_faces_vector)
          idx_faces_vector->push_back(point_index_map[ points_of_face[i] ]);
        if(!m_data_is_indexed)
        {
          // Its color
          if (face_color_vector && m_started_face_is_colored)
          { add_color(color_of_face, *face_color_vector); }

          // Its flat normal
          if(flat_normal_vector)
            add_normal(normal_for_face, *flat_normal_vector);

          // Its smoth normal (if given by the user)
          if (smooth_normal_vector && vertex_normals_for_face.size()==3)
          { // Here we have 3 vertex normals; we can use Gourod
            add_normal(vertex_normals_for_face[i], *smooth_normal_vector);
          }
          else if(smooth_normal_vector)
          { // Here user does not provide all vertex normals: we use face normal istead
            // and thus we will not be able to use Gourod
            add_normal(normal_for_face, *smooth_normal_vector);
          }
        }
      }
    }
    else if(is_facet_convex(points_of_face, normal_for_face))
    {
      if (points_of_face.size()==4) // Quad: no need for complex triangulation
      {
        if(faces_vector)
        {
          add_point(points_of_face[0], *faces_vector);
          add_point(points_of_face[1], *faces_vector);
          add_point(points_of_face[2], *faces_vector);

          add_point(points_of_face[0], *faces_vector);
          add_point(points_of_face[2], *faces_vector);
          add_point(points_of_face[3], *faces_vector);
        }
        else if(idx_faces_vector)
        {
          idx_faces_vector->push_back(point_index_map[ points_of_face[0] ]);
          idx_faces_vector->push_back(point_index_map[ points_of_face[1] ]);
          idx_faces_vector->push_back(point_index_map[ points_of_face[2] ]);

          idx_faces_vector->push_back(point_index_map[ points_of_face[0] ]);
          idx_faces_vector->push_back(point_index_map[ points_of_face[2] ]);
          idx_faces_vector->push_back(point_index_map[ points_of_face[3] ]);
        }
        if(!m_data_is_indexed)
          for(int i=0; i<6; ++i)
          {
            if (face_color_vector && m_started_face_is_colored)
              add_color(color_of_face, *face_color_vector);
            if(flat_normal_vector)
              add_normal(normal_for_face, *flat_normal_vector);
            if (smooth_normal_vector
                && vertex_normals_for_face.size()==4)
              add_normal(vertex_normals_for_face[0], *smooth_normal_vector);
            else if (smooth_normal_vector)
              add_normal(normal_for_face, *smooth_normal_vector);
          }
      }
      else //convex face with more than 4 points
      {
        Local_point p0,p1,p2;
        std::size_t id;
        for(id = 1; id<points_of_face.size()-1; ++id)
        {

          p0 = points_of_face[0];
          p1 = points_of_face[id];
          p2 = points_of_face[id+1];
          if(faces_vector)
          {
            add_point(p0, *faces_vector);
            add_point(p1, *faces_vector);
            add_point(p2, *faces_vector);
          }
          else if(idx_faces_vector)
          {
            idx_faces_vector->push_back(point_index_map[p0]);
            idx_faces_vector->push_back(point_index_map[p1]);
            idx_faces_vector->push_back(point_index_map[p2]);
          }

          if(!m_data_is_indexed)
          {
            if (face_color_vector && m_started_face_is_colored)
              for (int i = 0; i<6; ++i)
              {
                add_color(color_of_face, *face_color_vector);
              }
            if(flat_normal_vector)
              add_normal(normal_for_face, *flat_normal_vector);
            if (smooth_normal_vector
                && vertex_normals_for_face.size()==4)
              add_normal(vertex_normals_for_face[0], *smooth_normal_vector);
            else if (smooth_normal_vector)
              add_normal(normal_for_face, *smooth_normal_vector);
          }
        }
      }
    }
    else
    {

      try
      {
        P_traits cdt_traits(normal_for_face);
        CDT cdt(cdt_traits);

        bool with_vertex_normal=(vertex_normals_for_face.size()==points_of_face.size());

        // (1) We insert all the edges as contraint in the CDT.
        typename CDT::Vertex_handle previous=NULL, first=NULL;
        for (std::size_t i=0; i<points_of_face.size(); ++i)
        {
          typename CDT::Vertex_handle vh = cdt.insert(points_of_face[i]);
          if(first==NULL)
          { first=vh; }

          if (with_vertex_normal)
          { vh->info().v=vertex_normals_for_face[i]; }
          else
          { vh->info().v=normal_for_face; }

          if(previous!=NULL && previous!=vh)
          { cdt.insert_constraint(previous, vh); }
          previous=vh;
        }

        if (previous!=NULL && previous!=first)
          cdt.insert_constraint(previous, first);

        // (2) We mark all external triangles
        // (2.1) We initialize is_external and is_process values
        for(typename CDT::All_faces_iterator fit = cdt.all_faces_begin(),
            fitend = cdt.all_faces_end(); fit!=fitend; ++fit)
        {
          fit->info().is_external = true;
          fit->info().is_process = false;
        }
        // (2.2) We check if the facet is external or internal
        std::queue<typename CDT::Face_handle> face_queue;
        typename CDT::Face_handle face_internal = NULL;
        if (cdt.infinite_vertex()->face()!=NULL)
          face_queue.push(cdt.infinite_vertex()->face());
        while(! face_queue.empty() )
        {
          typename CDT::Face_handle fh = face_queue.front();
          face_queue.pop();
          if(!fh->info().is_process)
          {
            fh->info().is_process = true;
            for(int i=0; i<3; ++i)
            {
              if(!cdt.is_constrained(std::make_pair(fh, i)))
              {
                if (fh->neighbor(i)!=NULL)
                  face_queue.push(fh->neighbor(i));
              }
              else if (face_internal==NULL)
              {
                face_internal = fh->neighbor(i);
              }
            }
          }
        }

        if ( face_internal!=NULL )
          face_queue.push(face_internal);

        while(! face_queue.empty() )
        {
          typename CDT::Face_handle fh = face_queue.front();
          face_queue.pop();
          if(!fh->info().is_process)
          {
            fh->info().is_process = true;
            fh->info().is_external = false;
            for(int i=0; i<3; ++i)
            {
              if(!cdt.is_constrained(std::make_pair(fh, i)))
              {
                if (fh->neighbor(i)!=NULL)
                  face_queue.push(fh->neighbor(i));
              }
            }
          }
        }

        // (3) Now we iterates on the internal faces to add the vertices to the
        //     positions and the normals to the appropriate vectors
        for(typename CDT::Finite_faces_iterator ffit=cdt.finite_faces_begin(),
            ffitend = cdt.finite_faces_end(); ffit!=ffitend; ++ffit)
        {
          if(!ffit->info().is_external)
          {
            for(int i=0; i<3; ++i)
            {
              // The point
              if(faces_vector)
                add_point(ffit->vertex(i)->point(), *faces_vector);
              else if (idx_faces_vector)
                idx_faces_vector->push_back(point_index_map[ffit->vertex(i)->point()]);
              if(!m_data_is_indexed)
              {
                // Its color
                if (face_color_vector && m_started_face_is_colored)
                { add_color(color_of_face, *face_color_vector); }

                // Its flat normal
                if(flat_normal_vector)
                  add_normal(normal_for_face, *flat_normal_vector);

                // Its smooth normal (if given by the user)
                if(smooth_normal_vector)
                  add_normal(ffit->vertex(i)->info().v, *smooth_normal_vector);
              }
            }
          }
        }
      }
      catch(...)
      { // Triangulation crash: the face is not filled
        std::cout<<"Catch: face not filled."<<std::endl;
      }
    }

    m_face_started=false;
    points_of_face.clear();
    vertex_normals_for_face.clear();
  }


  template<typename KPoint>
  static void add_point(const KPoint& kp, std::vector<float>& point_vector)
  {
    Local_point p=internal::get_local_point(kp);
    point_vector.push_back(p.x());
    point_vector.push_back(p.y());
    point_vector.push_back(p.z());
  }

  static void add_color(const CGAL::Color& acolor, std::vector<float>& face_color_vector)
  {
    face_color_vector.push_back((float)acolor.red()/(float)255);
    face_color_vector.push_back((float)acolor.green()/(float)255);
    face_color_vector.push_back((float)acolor.blue()/(float)255);
  }

  template<typename KVector>
  static void add_normal(const KVector& kv, std::vector<float>& normal_vector)
  {
    Local_vector n=internal::get_local_vector(kv);
    normal_vector.push_back(n.x());
    normal_vector.push_back(n.y());
    normal_vector.push_back(n.z());
  }

private:

  // Local variables, used when we started a new face.
  bool m_face_started;
  bool m_data_is_indexed;
  bool m_started_face_is_colored;
  std::vector<float>* segments_vector;
  std::vector<float>* segment_color_vector;
  std::vector<float>* faces_vector;
  std::vector<float>* flat_normal_vector;
  std::vector<float>* smooth_normal_vector;
  std::vector<float>* face_color_vector;
  std::vector<unsigned int>* idx_segments_vector;
  std::vector<unsigned int>* idx_faces_vector;
  std::vector<Local_point> points_vector;
  std::vector<Local_point> points_of_face;
  std::vector<Local_vector> vertex_normals_for_face;
  Local_vector normal_for_face;
  CGAL::Color color_of_face;
  std::map<Local_point, std::size_t> point_index_map;
};
}//end namespace CGAL
#endif // PRIMITIVES_FILLER_H

