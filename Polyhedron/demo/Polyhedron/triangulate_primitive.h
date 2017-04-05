#ifndef TRIANGULATE_PRIMITIVE
#define TRIANGULATE_PRIMITIVE
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_2_projection_traits_3.h>
#include <CGAL/Three/Scene_item.h>
#include <queue>

#include <QColor>

//Make sure all the facets are triangles
template<class Mesh, typename Kernel, typename Index_type>
class FacetTriangulator
{
 typedef Kernel Traits;

 typedef typename boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;

 typedef typename Kernel::Vector_3 Vector;


 typedef CGAL::Triangulation_2_projection_traits_3<Traits>   P_traits;

 typedef CGAL::Triangulation_vertex_base_with_info_2<halfedge_descriptor,
                                                     P_traits>        Vb;

 struct Face_info {
   typename boost::graph_traits<Mesh>::halfedge_descriptor e[3];
   bool is_external;
 };

 typedef CGAL::Triangulation_face_base_with_info_2<Face_info,
                                                   P_traits>          Fb1;
 typedef CGAL::Constrained_triangulation_face_base_2<P_traits, Fb1>   Fb;
 typedef CGAL::Triangulation_data_structure_2<Vb,Fb>                  TDS;
 typedef CGAL::Exact_predicates_tag                                   Itag;

public:
 struct PointAndId {
  typename Kernel::Point_3 point;
  Index_type id;
 };

 typedef CGAL::Constrained_Delaunay_triangulation_2<P_traits,
                                                     TDS,
                                                     Itag>             CDT;
  CDT *cdt;
  CGAL::Unique_hash_map<typename CDT::Vertex_handle, Index_type> v2v;

  //Constructor
  FacetTriangulator(typename boost::graph_traits<Mesh>::face_descriptor fd,
                  const Vector& normal,
                    Mesh *poly,
                  const double item_diag)
  {
    std::vector<PointAndId> idPoints;
    BOOST_FOREACH(halfedge_descriptor he_circ, halfedges_around_face( halfedge(fd, *poly), *poly))
    {
      PointAndId idPoint;
      idPoint.point = get(boost::vertex_point,*poly,source(he_circ, *poly));
      idPoint.id = source(he_circ, *poly);
      idPoints.push_back(idPoint);

    }
    triangulate(idPoints, normal, item_diag);
  }
  FacetTriangulator(typename boost::graph_traits<Mesh>::face_descriptor fd,
                    const std::vector<typename Kernel::Point_3>& more_points,
                    const Vector& normal,
                    Mesh *poly,
                    const double item_diag)
  {
   std::vector<PointAndId> idPoints;
   BOOST_FOREACH(halfedge_descriptor he_circ, halfedges_around_face( halfedge(fd, *poly), *poly))
   {
    PointAndId idPoint;
    idPoint.point = get(boost::vertex_point,*poly,source(he_circ, *poly));
    idPoint.id = source(he_circ, *poly);
    idPoints.push_back(idPoint);

   }
   triangulate_with_points(idPoints,more_points, normal, item_diag);
  }

  FacetTriangulator(std::vector<PointAndId > &idPoints,
                  const Vector& normal,
                  const double item_diag)
  {
    triangulate(idPoints, normal, item_diag);
  }
  FacetTriangulator(std::vector<PointAndId > &idPoints,
                    const std::vector<typename Kernel::Point_3>& more_points,
                    const Vector& normal,
                    const double item_diag)
  {
   triangulate_with_points(idPoints, more_points, normal, item_diag);
  }

private:
  void triangulate( std::vector<PointAndId > &idPoints,
              const Vector& normal,
              const double item_diag )
  {
    P_traits cdt_traits(normal);
    cdt = new CDT(cdt_traits);
    typename CDT::Vertex_handle previous, first, last_inserted;
    //Compute a reasonable precision level used to decide
    //if two consecutive points in a facet can be estimated
    //equal.

    double min_sq_dist = CGAL::square(0.0001*item_diag);
    // Iterate the points of the facet and decide if they must be inserted in the CDT
    typename Kernel::FT x(0), y(0), z(0);

    BOOST_FOREACH(PointAndId idPoint, idPoints)
    {
     x += idPoint.point.x();
     y += idPoint.point.y();
     z += idPoint.point.z();
     typename CDT::Vertex_handle vh;
     //Always insert the first point, then only insert
     // if the distance with the previous is reasonable.
     if(first == 0 || CGAL::squared_distance(idPoint.point, previous->point()) > min_sq_dist)
     {
       vh = cdt->insert(idPoint.point);
       v2v[vh] = idPoint.id;
       if(first == 0) {
         first = vh;
       }
       if(previous != 0 && previous != vh) {
         double sq_dist = CGAL::squared_distance(previous->point(), vh->point());
         if(sq_dist > min_sq_dist)
         {
           cdt->insert_constraint(previous, vh);
           sq_dist = CGAL::squared_distance(previous->point(), first->point());
           if(sq_dist > min_sq_dist)
           {
             last_inserted = previous;
           }
         }
       }
     previous = vh;
     }
    }
    double sq_dist = CGAL::squared_distance(previous->point(), first->point());

    if(sq_dist > min_sq_dist)
    {
      cdt->insert_constraint(previous, first);
    }
    else
    {
      cdt->insert_constraint(last_inserted, first);
    }
    // sets mark is_external
    for(typename CDT::All_faces_iterator
        fit2 = cdt->all_faces_begin(),
        end = cdt->all_faces_end();
        fit2 != end; ++fit2)
    {
      fit2->info().is_external = false;
    }
    //check if the facet is external or internal
    std::queue<typename CDT::Face_handle> face_queue;
    face_queue.push(cdt->infinite_vertex()->face());
    while(! face_queue.empty() ) {
      typename CDT::Face_handle fh = face_queue.front();
      face_queue.pop();
      if(fh->info().is_external) continue;
      fh->info().is_external = true;
      for(int i = 0; i <3; ++i) {
        if(!cdt->is_constrained(std::make_pair(fh, i)))
        {
          face_queue.push(fh->neighbor(i));
        }
      }
    }
  }

  void triangulate_with_points( std::vector<PointAndId > &idPoints,
               const std::vector<typename Kernel::Point_3>& more_points,
               const Vector& normal,
               const double item_diag )
   {
     P_traits cdt_traits(normal);
     cdt = new CDT(cdt_traits);
     typename CDT::Vertex_handle previous, first, last_inserted;
     //Compute a reasonable precision level used to decide
     //if two consecutive points in a facet can be estimated
     //equal.

     double min_sq_dist = CGAL::square(0.0001*item_diag);
     // Iterate the points of the facet and decide if they must be inserted in the CDT
     BOOST_FOREACH(PointAndId idPoint, idPoints)
     {
      typename CDT::Vertex_handle vh;
      //Always insert the first point, then only insert
      // if the distance with the previous is reasonable.
      if(first == 0 || CGAL::squared_distance(idPoint.point, previous->point()) > min_sq_dist)
      {
        vh = cdt->insert(idPoint.point);
        v2v[vh] = idPoint.id;
        if(first == 0) {
          first = vh;
        }
        if(previous != 0 && previous != vh) {
          double sq_dist = CGAL::squared_distance(previous->point(), vh->point());
          if(sq_dist > min_sq_dist)
          {
            cdt->insert_constraint(previous, vh);
            sq_dist = CGAL::squared_distance(previous->point(), first->point());
            if(sq_dist > min_sq_dist)
            {
              last_inserted = previous;
            }
          }
        }
      previous = vh;
      }
     }
     double sq_dist = CGAL::squared_distance(previous->point(), first->point());

     if(sq_dist > min_sq_dist)
     {
       cdt->insert_constraint(previous, first);
     }
     else
     {
       cdt->insert_constraint(last_inserted, first);
     }
     BOOST_FOREACH(typename Kernel::Point_3 point, more_points)
     {
       cdt->insert(point);
     }
     // sets mark is_external
     for(typename CDT::All_faces_iterator
         fit2 = cdt->all_faces_begin(),
         end = cdt->all_faces_end();
         fit2 != end; ++fit2)
     {
       fit2->info().is_external = false;
     }
     //check if the facet is external or internal
     std::queue<typename CDT::Face_handle> face_queue;
     face_queue.push(cdt->infinite_vertex()->face());
     while(! face_queue.empty() ) {
       typename CDT::Face_handle fh = face_queue.front();
       face_queue.pop();
       if(fh->info().is_external) continue;
       fh->info().is_external = true;
       for(int i = 0; i <3; ++i) {
         if(!cdt->is_constrained(std::make_pair(fh, i)))
         {
           face_queue.push(fh->neighbor(i));
         }
       }
     }
  }
};

#endif // TRIANGULATE_PRIMITIVE

