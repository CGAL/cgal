#ifndef COMPUTE_ANCHOR_3_H
#define COMPUTE_ANCHOR_3_H

#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Simplex_3.h>

CGAL_BEGIN_NAMESPACE

template < class RegularTriangulation3>
class Compute_anchor_3
{
public:
  typedef RegularTriangulation3                        Regular;
  typedef typename Regular::Geom_traits                Geom_traits;
  typedef typename Geom_traits::Weighted_point         Weighted_point;
  
  typedef typename Regular::Vertex_handle Vertex_handle;
  typedef typename Regular::Cell_handle   Cell_handle;
  typedef typename Regular::Facet         Facet;
  typedef typename Regular::Edge          Edge;
  
  typedef typename Regular::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Regular::Finite_edges_iterator    Finite_edges_iterator;
  typedef typename Regular::Finite_facets_iterator   Finite_facets_iterator;
  typedef typename Regular::Finite_cells_iterator    Finite_cells_iterator;
  typedef typename Regular::Facet_circulator         Facet_circulator;
  typedef Simplex_3<Regular>                         Simplex;


  Compute_anchor_3(Regular &reg) : reg(reg) {
  }

  ///////////////////////////////
  // Anchor functions
  ///////////////////////////////

  // Thest whether the anchor of edge (wp1,wp2) and wp2 are equal
  Sign test_anchor(Weighted_point &wp1, Weighted_point &wp2) {
    return 
      reg.geom_traits().in_smallest_orthogonal_sphere_3_object()(wp1, wp2);
  }
  // Test whether the anchor of e and anchor of e.first->vertex(i) are equal
  Sign test_anchor(Edge e, int i) {
    CGAL_assertion(!reg.is_infinite(e));
    CGAL_assertion(e.second == i || e.third == i);
    Weighted_point wp1, wp2;
    Cell_handle ch = e.first;
    return test_anchor(
      ch->vertex(e.second+e.third-i)->point(),
      ch->vertex(i)->point());
  }
  // Test whether the anchor of f and anchor of the edge f - f.first->vertex(i)
  // are equal
  Sign test_anchor(Facet f, int i){
    CGAL_assertion(!reg.is_infinite(f));
    CGAL_assertion(f.second != i);
    CGAL_assertion(0<=f.second && f.second<4);
    CGAL_assertion(0<=i && i<4);

    Weighted_point wp1, wp2, wp3;
    Cell_handle ch = f.first;
    wp3 = ch->vertex(i)->point();
    switch ((4+i-f.second)&3) {
      case 1: CGAL_assertion (((f.second+1)&3) == i);
        wp1 = ch->vertex((f.second+2)&3)->point();
        wp2 = ch->vertex((f.second+3)&3)->point();
        break;
      case 2: CGAL_assertion (((f.second+2)&3) == i);
        wp1 = ch->vertex((f.second+1)&3)->point();
        wp2 = ch->vertex((f.second+3)&3)->point();
        break;
      case 3: CGAL_assertion (((f.second+3)&3) == i);
        wp1 = ch->vertex((f.second+1)&3)->point();
        wp2 = ch->vertex((f.second+2)&3)->point();
        break;
      default:
        CGAL_assertion(false);
    }

    return 
      reg.geom_traits().in_smallest_orthogonal_sphere_3_object()(wp1, wp2, wp3);
  }


  // Test whether the anchor of ch and anchor of the facet (ch,i) are equal
  Sign test_anchor(Cell_handle ch, int i) {
    CGAL_assertion(!reg.is_infinite(ch));

    return reg.geom_traits().in_smallest_orthogonal_sphere_3_object()(
      ch->vertex((i+1)&3)->point(),
      ch->vertex((i+2)&3)->point(),
      ch->vertex((i+3)&3)->point(),
      ch->vertex(i)->point());		
  }

  Simplex compute_anchor_del( Edge const &e );
  Simplex compute_anchor_del( Facet const &f );
  Simplex compute_anchor_del( Cell_handle const ch );
  
  Simplex compute_anchor_vor( Vertex_handle const v );
  Simplex compute_anchor_vor( Edge const &e );
  Simplex compute_anchor_vor( Facet const &f );

  Simplex anchor_del( Vertex_handle v ) {
    return v;
  }
  Simplex anchor_del( Edge const &e ) {
    return compute_anchor_del(e);
  }
  Simplex anchor_del( Facet const &f ) {
    return compute_anchor_del(f);
  }
  Simplex anchor_del( Cell_handle ch ) {
    return compute_anchor_del(ch);
  }
  Simplex anchor_vor( Vertex_handle v ) {
    return compute_anchor_vor(v);
  }
  Simplex anchor_vor( Edge const &e ) {
    return compute_anchor_vor(e);
  }
  Simplex anchor_vor( Facet const &f ) {
    return compute_anchor_vor(f);
  }
  Simplex anchor_vor( Cell_handle const ch ) {
    return ch;
  }

private:
  Regular &reg;
};


// compute_anchor_del

template < class RegularTriangulation3 >
typename Compute_anchor_3<RegularTriangulation3>::Simplex
Compute_anchor_3<RegularTriangulation3>::
compute_anchor_del(Edge const &e) {
  CGAL_assertion(!reg.is_infinite(e));

  if (test_anchor(e, e.second)==NEGATIVE) {
    return Simplex(e.first->vertex(e.third));
  }
  if (test_anchor(e, e.third)==NEGATIVE) {
    return Simplex(e.first->vertex(e.second));
  } 
  return Simplex(e);
}

template < class RegularTriangulation3 >
typename Compute_anchor_3<RegularTriangulation3>::Simplex
Compute_anchor_3<RegularTriangulation3>::
compute_anchor_del(Facet const &f) {
  CGAL_assertion(!reg.is_infinite(f));

  Simplex s;
  int i;
  Sign sign;
  bool contains_focus = true;
  for (i=1; (i<4) && contains_focus; i++) {
    sign = test_anchor(f, (f.second+i)&3);
    contains_focus = (sign != NEGATIVE);
  }
  
  if (contains_focus) {
    s = Simplex(f);
  } else {
    i--;
    Edge e; e.first = f.first; 
    if (i==1) e.second = ((f.second+2)&3); 
    else e.second = ((f.second+1)&3);
    if (i==3) e.third = ((f.second+2)&3); 
    else e.third = ((f.second+3)&3);

    s = anchor_del(e);
    if (s.dimension() == 1) {
      return s;
    } else {
      // The anchor is the anchor of the other edge adjacent to tmp
      CGAL_assertion(s.dimension() == 0);
      Vertex_handle vh=s;
      e.second = e.first->index(vh);
      e.third = (f.second+i)&3;
      return anchor_del(e);
    }
  }
  
  return s;
}

template < class RegularTriangulation3 >
typename Compute_anchor_3<RegularTriangulation3>::Simplex
Compute_anchor_3<RegularTriangulation3>::
compute_anchor_del(Cell_handle const ch) {
  CGAL_assertion(!reg.is_infinite(ch));
  Simplex s;

  bool contains_focus = true;
  Sign sign;
  for (int i=0; (i<4) && contains_focus; i++) {
    sign = test_anchor(ch, i);
    if (sign == NEGATIVE) {
      contains_focus = false;
      s = anchor_del(Facet(ch,i));
    }
  }

  if (contains_focus) {
    return Simplex(ch);
  } else {
    Simplex tmp;

    bool found=true;
    while (true) {
      if (s.dimension() == 2) {
        return s;
      } else if (s.dimension() == 1) {
        // Test two adjacent facets
        Edge e=s;
        int ind1 = ch->index(e.first->vertex(e.second));
        int ind2 = ch->index(e.first->vertex(e.third));
        for (int i=0; (i<4) && found; i++) {
          if ((i != ind1) && (i != ind2)) {
            tmp = anchor_del(Facet(ch, i));
            found = (s == tmp);
          }
        }
      } else if (s.dimension() == 0) {
        // Test adjacent edges
        Vertex_handle vh=s;
        int index = ch->index(vh);
        for (int i=0; (i<4) && found; i++) {
          if (i != index) {
            tmp = anchor_del(Edge(ch, index, i));
            found = (s == tmp);
          }
        }
      } else {
        CGAL_assertion(false);
      }
      if (found) return s;
      found = true;
      s = tmp;
    }
  }
}

template < class RegularTriangulation3 >
typename Compute_anchor_3<RegularTriangulation3>::Simplex
Compute_anchor_3<RegularTriangulation3>::
compute_anchor_vor (Vertex_handle const v) {
  CGAL_assertion(!reg.is_infinite(v));
  CGAL_assertion(reg.is_vertex(v));

  Sign side;
  bool contains_focus=true;
  Simplex s;

  std::list<Vertex_handle> adj_vertices;
  typename std::list<Vertex_handle>::iterator adj_vertex;
  reg.incident_vertices(v, std::back_inserter(adj_vertices));
  
  adj_vertex = adj_vertices.begin();
  do {
    if (!reg.is_infinite(*adj_vertex)) {
      side = test_anchor(v->point(),(*adj_vertex)->point());
      contains_focus = (side != NEGATIVE);
    }
    adj_vertex++;
  } while ((adj_vertex != adj_vertices.end()) && contains_focus);

  if (contains_focus) {
    s = Simplex(v);
  } else {
    adj_vertex--;
    Edge e;
    if (!reg.is_edge(v, *adj_vertex, e.first, e.second, e.third)) {
      CGAL_assertion(false);
    }
    s = anchor_vor(e);
    while (1) {
      if (s.dimension() == 1) {
        // s lies on a Voronoi facet
        return s; 
      } else if (s.dimension() == 2) {
        // s lies on a Voronoi edge
        Facet f=s;
        int index = f.first->index(v);
        bool found = false;
        for (int i=1; (i<4) && (!found); i++) {
          if (((f.second+i)&3) != index) {
            Simplex stmp;
            stmp = anchor_vor(Edge(f.first, index, (f.second+i)&3));
            if (!(stmp == s)) {
              s = stmp; 
              found = true;
            }
          }
        }
        if (!found) return s;
      } else {
        //  s lies on a Voronoi vertex
        CGAL_assertion(s.dimension() == 3);
        Cell_handle ch=s;
        int index = ch->index(v);
        bool found = false;
        for (int i=1; (i<4) && (!found); i++) {
          Simplex stmp;
          stmp = anchor_vor(Facet(ch, (index+i)&3));
          if (!(stmp == s)) {
            s = stmp;
            found = true;
          }
        }
        if (!found) return s;
      }
    }
  }
  return s;
}

template < class RegularTriangulation3 >
typename Compute_anchor_3<RegularTriangulation3>::Simplex
Compute_anchor_3<RegularTriangulation3>::compute_anchor_vor (Edge const &e) {
  CGAL_assertion(!reg.is_infinite(e));

  Vertex_handle v0 = e.first->vertex(e.second);
  Vertex_handle v1 = e.first->vertex(e.third);

  bool contains_focus = true;
  Sign side = CGAL::NEGATIVE;
  Simplex s;
  
  Facet_circulator fcir, fstart;
  fcir = fstart = reg.incident_facets(e);
  int i;
  do {
    if (!reg.is_infinite(*fcir)) {
      i = 6 - (*fcir).second -
	(*fcir).first->index(v0) - (*fcir).first->index(v1);
      side = test_anchor(*fcir, i);
      if (side == NEGATIVE) {
        contains_focus = false;
        s = anchor_vor(Facet(*fcir));
      }
    }
    fcir++;
  } while (fcir != fstart && contains_focus);
  
  if (contains_focus) {
    s = Simplex(e);
    return s;
  } else {
    Simplex tmp;
    while (1) {
      CGAL_assertion(s.dimension() > 1);
      if (s.dimension() == 2) {
        return s;
      } else {
        CGAL_assertion(s.dimension() == 3);
        bool found = true;
        Cell_handle ch=s;
        int index0 = ch->index(v0);
        int index1 = ch->index(v1);
        for (int i=0; (i<4) && found; i++) {
          if ((i != index0) && (i != index1)) {
	    if (!reg.is_infinite(Facet(ch,i))) {
	      side = test_anchor(ch,6-index0-index1-i);
	      // Cannot be ZERO: checked earlier
	      if (side == NEGATIVE) {
		tmp = anchor_vor(Facet(ch,i));
		if (!(tmp == s)) {
		  s = tmp;
		  found = false;
		}
	      }
            }
          }
        }
        if (found) return s;
      }
    }
  }
}

template < class RegularTriangulation3 >
typename Compute_anchor_3<RegularTriangulation3>::Simplex
Compute_anchor_3<RegularTriangulation3>::compute_anchor_vor (Facet const &f) {
  CGAL_assertion(!reg.is_infinite(f));

  if (!reg.is_infinite(f.first)) {
    if (test_anchor(f.first, f.second)==NEGATIVE) {
      return anchor_vor(f.first);
    }
  }

  Cell_handle neighbor = f.first->neighbor(f.second);
  if (!reg.is_infinite(neighbor)) {
    int n_index = neighbor->index(f.first);
    if (test_anchor(neighbor, n_index)==NEGATIVE) {
      return anchor_vor(neighbor);
    }
  }

  return Simplex(f);
}


CGAL_END_NAMESPACE

#endif // COMPUTE_ANCHOR_3_H
