
#ifndef GPS_BFS_JOIN_VISITOR_H
#define GPS_BFS_JOIN_VISITOR_H

#include <CGAL/Boolean_set_operations_2/Gps_bfs_base_visitor.h>

CGAL_BEGIN_NAMESPACE

template <class Arrangement_>
class Gps_bfs_join_visitor : public Gps_bfs_base_visitor<Arrangement_>
{
  typedef  Arrangement_                                  Arrangement;
  typedef typename Arrangement::Face_iterator            Face_iterator;
  typedef typename Arrangement::Halfedge_iterator        Halfedge_iterator;
  typedef Gps_bfs_base_visitor<Arrangement>              Base;
  typedef typename Base::Edges_hash                      Edges_hash;
  typedef typename Base::Faces_hash                      Faces_hash;

public:

  Gps_bfs_join_visitor(Edges_hash* edges_hash, Faces_hash* faces_hash, unsigned int n_pgn): 
    Base(edges_hash, faces_hash, n_pgn)
  {}


  void flip_face(Face_iterator f1, Face_iterator f2, Halfedge_iterator he)
  {
    unsigned int ic_f2;
    ic_f2 = this->compute_ic(f1, f2, he);
    (*(this->m_faces_hash))[f2] = ic_f2;
      
    if(ic_f2 > 0)
      f2->set_contained(true);
  }

  // mark the unbounded_face 
  void visit_ubf(Face_iterator ubf, unsigned int ubf_ic)
  {
    CGAL_assertion(ubf->is_unbounded());
    if(ubf_ic > 0)
      ubf->set_contained(true);
  }

};

CGAL_END_NAMESPACE

#endif
