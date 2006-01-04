
#ifndef GPS_BFS_INTERSECTION_VISITOR_H
#define GPS_BFS_INTERSECTION_VISITOR_H

#include <CGAL/Boolean_set_operations_2/Gps_bfs_base_visitor.h>

CGAL_BEGIN_NAMESPACE

template <class Arrangement_>
class Gps_bfs_intersection_visitor : public Gps_bfs_base_visitor<Arrangement_>
{
  typedef  Arrangement_                                  Arrangement;
  typedef typename Arrangement::Face_iterator            Face_iterator;
  typedef typename Arrangement::Halfedge_iterator        Halfedge_iterator;
  typedef Gps_bfs_base_visitor<Arrangement>              Base;
  typedef typename Base::Edges_hash                      Edges_hash;
  typedef typename Base::Faces_hash                      Faces_hash;
 
protected:

  unsigned int m_num_of_polygons;
public:

  Gps_bfs_intersection_visitor(Edges_hash* edges_hash,
                               Faces_hash* faces_hash,
                               unsigned int n_polygons): 
    Base(edges_hash, faces_hash),
    m_num_of_polygons(n_polygons)
  {}


  void flip_face(Face_iterator f1, Face_iterator f2, Halfedge_iterator he)
  {
    unsigned int ic_f2;
    ic_f2 = compute_ic(f1, f2, he);
    (*m_faces_hash)[f2] = ic_f2;
      
    CGAL_assertion(ic_f2 <= m_num_of_polygons);

    // only faces that have inside counter equal to the number of polygons
    // which are intersectd, will be marked true (containted)
    if(ic_f2 == m_num_of_polygons)
      f2->set_contained(true);
  }


};
CGAL_END_NAMESPACE

#endif
