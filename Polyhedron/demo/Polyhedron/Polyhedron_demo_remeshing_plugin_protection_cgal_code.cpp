#include <CGAL/config.h>
#include "Polyhedron_type.h"
#include "C2t3_type.h"
#include <set>
#include <map>

typedef Tr::Geom_traits::FT FT;

struct Less {
  template <typename Handle>
  bool operator()(const Handle& va, const Handle& vb) const {
    return &*va < &*vb;
  }
};

struct Insert_spheres {
  typedef Polyhedron::Halfedge_const_handle Halfedge_const_handle;
  typedef Polyhedron::Vertex_const_handle Vertex_const_handle;
  typedef Polyhedron::size_type size_type;

  typedef Tr::Geom_traits::Point_3 Point_3;

  typedef std::set<Vertex_const_handle, Less> Vertices_set;
  typedef std::map<Vertex_const_handle, size_type, Less> Vertices_counter;

  typedef std::set<Halfedge_const_handle, Less> Border_edges_set;

  Border_edges_set edges_to_consider;
  Vertices_counter border_vertices;
  Vertices_set corner_vertices;

  typedef std::vector<Tr::Geom_traits::Point_3> Polyline;
  
  typedef std::vector<Polyline> Polylines;
  
  Polylines polylines;

  C2t3& c2t3;
  Polyhedron* pMesh;
  const FT size;

  Halfedge_const_handle canonical(Halfedge_const_handle he) 
  {
    const Halfedge_const_handle& op = he->opposite();
    if(Less()(he, op))
      return he;
    else 
      return op;
  }

  /** Follow a polyline or a polygon, from the halfedge he. */
  void follow_half_edge(const Halfedge_const_handle he)
  {
    Border_edges_set::iterator it = edges_to_consider.find(canonical(he));
    if(it == edges_to_consider.end()) {
      return;
    }

    Polyline polyline;
    polyline.push_back(he->opposite()->vertex()->point());
    // std::cerr << "Start: " << he->opposite()->vertex()->point() << std::endl;

    Halfedge_const_handle current_he = he;
    do {
      CGAL_assertion(current_he->is_border() || 
                     current_he->opposite()->is_border());
      const size_type n = edges_to_consider.erase(canonical(current_he));
      CGAL_assertion(n > 0);
      Vertex_const_handle v = current_he->vertex();
      polyline.push_back(v->point());
      // std::cerr << v->point() << std::endl;
      if(corner_vertices.count(v) > 0) break;
      Polyhedron::Halfedge_around_vertex_const_circulator 
        loop_he = v->vertex_begin(), end(loop_he);
      ++loop_he;
      // CGAL_assertion((&*loop_he) != (&*current_he) );
      while((&*loop_he) == (&*current_he) ||
            (!loop_he->is_border() && !loop_he->opposite()->is_border()) ) {
        ++loop_he;
        // CGAL_assertion((&*loop_he) != (&*current_he) );
      }
      current_he = loop_he->opposite();
    } while(current_he != he );

    if(current_he == he)
      std::cerr << "New polyline, of size " << polyline.size() << std::endl;
    else
      std::cerr << "New polygon (cycle), of size " << polyline.size() << std::endl;

    polylines.push_back(polyline);
  }

  /** Loop around a corner vertex, and try to follow a polyline of border
      edges, from each incident edge. */
  void loop_around_corner(const Vertex_const_handle v)
  {
    Polyhedron::Halfedge_around_vertex_const_circulator 
      he = v->vertex_begin(), end(he);
    do {
      CGAL_assertion(he->vertex() == v);
      follow_half_edge(he->opposite());
      ++he;
    } while(he != end);
  }

  /** For a non-corner vertex v (that is incident to two border edges),
      mesure the angle between the two edges, and mark the vertex as corner
      edge, if the angle is < 120°. **/
  void mesure_angle(const Vertex_const_handle v)
  {
    Halfedge_const_handle e1;
    Halfedge_const_handle e2;
    Polyhedron::Halfedge_around_vertex_const_circulator he = v->vertex_begin(), end(he);
    // std::cerr << "mesure_handle(" << (void*)(&*v)
    //           << " = " << v->point() << ")";
    bool first = true;
    do {
      CGAL_assertion(he->vertex() == v);
      // std::cerr << he->opposite()->vertex()->point() << std::endl;
      if(he->is_border() || he->opposite()->is_border()) {
        if(first) {
          e1 = he;
          first = false;
        }
        else {
          CGAL_assertion(e2 == Halfedge_const_handle());
          e2 = he;
        }
        std::cerr << "x";
      }
      else
        std::cerr << ".";
      ++he;
    } while(he != end);
      std::cerr << "\n";
    const Point_3 pv = v->point();
    const Point_3 pa = e1->opposite()->vertex()->point();
    const Point_3 pb = e2->opposite()->vertex()->point();
    const Tr::Geom_traits::Vector_3 av = pv - pa;
    const Tr::Geom_traits::Vector_3 bv = pv - pb;
    const FT sc_prod = av * bv;
    if( sc_prod >= 0 ||
        (sc_prod < 0 && 
         CGAL::square(sc_prod) < (av * av) * (bv * bv) / 4 ) )
    {
      std::cerr << "Corner (" << pa << ", " << pv
                << ", " << pb << ")\n";
      corner_vertices.insert(v);
    }
  }

  void protect(Polyline::const_iterator begin,
               Polyline::const_iterator end) {
    CGAL_assertion(begin+1 != end);
    // for(Polyline::const_iterator it
    Polyline::const_iterator end2 = end;
    --end2;
    FT distance = 0;
    for(Polyline::const_iterator it = begin;
        it != end2;) {
      const Point& a = *it;
      const Point& b = *++it;
      distance += CGAL_NTS sqrt( CGAL::squared_distance(a, b) );
    }
    std::cerr << "Distance:   " << distance << std::endl;
    std::cerr << "Size:       " << size << std::endl;
    const size_type n = static_cast<size_type>(std::ceil(distance / size) + 0.5);
    const FT local_size = distance / n;
    std::cerr << n << std::endl;
    std::cerr << "Local size: " << local_size << std::endl;
    CGAL_assertion(local_size < size);
    
    Point_3 a(begin->point(), local_size/1.5);
    Point_3 b(end2->point(), local_size/1.5);
    c2t3.triangulation().insert(a);
    Polyline::const_iterator it = begin;
    ++it;
    FT small_distance_to_go = local_size;
    while(it != end) {
      const Point& a = *it;
      const Point& b = *++it;
      const FT  d = CGAL_NTS squared_distance(a, b);
      unsigned i = 0;
      for(; small_distance_to_go + i * local_size >= d;
          ++i)
      {
        const Point p = a +
          (small_distance_to_go + i * local_size) * ( b - a ) / d;
        c2t3.triangulation().insert(Point_3(p, local_size / 1.5));
      }
      small_distance_to_go -= d;
      small_distance_to_go += i * local_size;
    }
    c2t3.triangulation().insert(b);
    std::cerr << "One polyline is protected!\n";
  }

  Insert_spheres(C2t3& c2t3, Polyhedron* pMesh, const FT size)
    : c2t3(c2t3), pMesh(pMesh), size(size)
  {
    // That call orders the set of edges of the polyhedron, so that the
    // border edges are at the end of the sequence of edges.
    pMesh->normalize_border();

    // Iterate over border edges, and find out which vertices are corner
    // vertices (more than two incident border edges).
    for(Polyhedron::Edge_const_iterator 
          eit = pMesh->border_edges_begin (),
          end = pMesh->edges_end();
        eit != end; ++eit)
    {
      edges_to_consider.insert(canonical(eit));
      Polyhedron::Vertex_const_handle v = eit->vertex();
      for(unsigned i = 0; i < 2; ++i) {
        if(++border_vertices[v] == 3)
          corner_vertices.insert(v);
        v = eit->opposite()->vertex();
      }

    }
    std::cerr << "Corner vertices: " << corner_vertices.size() << std::endl;
    std::cerr << "Border vertices: " << border_vertices.size() << std::endl;
    
    // Iterate over non-corner border vertices, and mesure the angle.
    for(Vertices_counter::iterator it = border_vertices.begin(),
        end = border_vertices.end(); it != end; ++it)
    {
      const Vertex_const_handle v = it->first;
      if(corner_vertices.count(v) == 0) {
        CGAL_assertion(it->second == 2);
        mesure_angle(v);
      }
    }
    std::cerr << "New corner vertices: " << corner_vertices.size() << std::endl;

    // Follow the polylines...
    for(Vertices_set::iterator it = corner_vertices.begin(),
        end = corner_vertices.end(); it != end; ++it)
    {
      loop_around_corner(*it);
    }

    // ... and the cycles.
    while(! edges_to_consider.empty() ) {
      follow_half_edge(*edges_to_consider.begin());
    }

    for(Polylines::const_iterator it = polylines.begin(),
          end = polylines.end(); it != end; ++it)
    {
      protect(it->begin(), it->end());
    }
  }
};

void insert_spheres(C2t3& c2t3, Polyhedron* pMesh, const FT size) {

  Insert_spheres go(c2t3, pMesh, size);
}
