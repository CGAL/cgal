#ifndef PVD_INSERT_H
#define PVD_INSERT_H


int get_new_id()
{
  static int gen = 0;
  gen++;

  assert( gen != 0 ); // guards against unwanted phenomena due to
                      // overflow 
  return gen;
}



template<class PVD, class Point>
typename PVD::Vertex_handle
insert_point(PVD& pvd, const Point& p)
{
  typename PVD::Vertex_handle v = pvd.insert(p);
  v->set_info( get_new_id() );
  return v;
}

template<class PVD, class Point>
typename PVD::Vertex_handle
insert_segment(PVD& pvd, const Point& p1, const Point_2& p2, int id)
{
  typedef typename PVD::Vertex_handle Vertex_handle;

  Vertex_handle v1 = pvd.insert(p1);
  Vertex_handle v2 = pvd.insert(p2);
  Vertex_handle v3 = pvd.insert(p1,p2);

  if ( v3 == Vertex_handle() ) {
    return v3;
  }

  v1->set_info( id );
  v2->set_info( id );
  v3->set_info( id );
  return v3;
}

template<class PVD, class Point>
typename PVD::Vertex_handle
insert_segment(PVD& pvd, const Point& p1, const Point_2& p2, 
	       typename PVD::Vertex_handle v, int id)
{
  typedef typename PVD::Vertex_handle Vertex_handle;

  Vertex_handle v1 = pvd.insert(p1, v);
  Vertex_handle v2 = pvd.insert(p2, v2);
  Vertex_handle v3 = pvd.insert(p1, p2, v1);

  if ( v3 == Vertex_handle() ) {
    return v3;
  }

  v1->set_info( id );
  v2->set_info( id );
  v3->set_info( id );
  return v3;
}


template<class PVD, class Point>
typename PVD::Vertex_handle
insert_segment(PVD& pvd, const Point& p1, const Point& p2)
{
  int id = get_new_id();
  return insert_segment(pvd, p1, p2, id);
}

template<class PVD, class Point>
typename PVD::Vertex_handle
insert_segment(PVD& pvd, const Point& p1, const Point& p2,
	       typename PVD::Vertex_handle v)
{
  int id = get_new_id();
  return insert_segment(pvd, p1, p2, v, id);
}


template<class PVD, class Polygon>
typename PVD::Vertex_handle
insert_polygon(PVD& pvd, const Polygon& pgn)
{
  typedef typename PVD::Vertex_handle Vertex_handle;

  int id = get_new_id();

  Vertex_handle v;
  int psize = pgn.size();
  for (int i = 0; i < psize; i++ ) {
    v = insert_segment( pvd, pgn[i], pgn[(i+1)%psize], id );
    if ( v == Vertex_handle() ) { break; }
  }
  return v;
}



#endif // PVD_INSERT_H
