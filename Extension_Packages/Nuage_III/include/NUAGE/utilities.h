//=====================================================================
#ifndef UTILITIES_H
#define UTILITIES_H
//=====================================================================

#include <CGAL/Unique_hash_map.h>

void
construct_delaunay(const std::vector<Point> &V_p,
		   Triangulation_3& A)
{
  std::cout << "   Compute Delaunay Tetrahedrization" << std::endl; 
  t1.start();
  {
    for(std::vector<Point>::const_iterator v_it = V_p.begin();
	v_it != V_p.end(); ++v_it)
      {
	A.insert(*v_it);
      }
  }
  t1.stop();
  std::cout << "   Inserted " << A.number_of_vertices() << " points, "
	    <<  A.number_of_cells() << " cells computed in "
	    << t1.time() << " secondes." << std::endl;
  std::cout << "   Number of filter failures : " << 
    CGAL::Interval_base::number_of_failures << std::endl;
  if (A.dimension() < 3)
    {
      std::cout << "-- 2D sample of points ???" 
		<< std::endl;
      std::exit(0);
    }
  t1.reset();
}

//=====================================================================
//=====================================================================

//private
inline Next_border_elt* get_border_elt(const Vertex_handle& v1, const Vertex_handle& v2)
{
  return v1->get_border_elt((void*) &(*v2));
}

//public

inline IO_edge_type* get_border_IO_elt(const Vertex_handle& v1, const Vertex_handle& v2)
{
  return &get_border_elt(v1,v2)->second.first.second;
}

inline IO_edge_type* set_border_elt(const Vertex_handle& v1, const Vertex_handle& v2,
				    const Border_elt& e)
{
  v1->set_next_border_elt(Next_border_elt ((void*) &(*v2), e));
  return get_border_IO_elt(v1, v2);
}


inline IO_edge_type* set_again_border_elt(const Vertex_handle& v1, const Vertex_handle& v2,
					  const Border_elt& e)
{
  get_border_elt(v1,v2)->second = e;
  return get_border_IO_elt(v1, v2);
}

//---------------------------------------------------------------------

inline bool is_border_elt(Edge_like& key, Border_elt& result)
{
  Next_border_elt* it12 =  key.first->get_border_elt((void*) &(*key.second));
  if (it12 != NULL)
    {    
      result = it12->second;
      key = Edge_like(key.first, key.second);
      return true;
    }

  Next_border_elt* it21 =  key.second->get_border_elt((void*) &(*key.first));
  if (it21 != NULL)
    {    
      result = it21->second;
      key = Edge_like(key.second, key.first);
      return true;
    }
  return false;
}

//---------------------------------------------------------------------

inline bool is_ordered_border_elt(const Edge_like& key, Border_elt& result)
{
  Next_border_elt* it12 =  key.first->get_border_elt((void*) &(*key.second));
  if (it12 != NULL)
    {    
      result = it12->second;
      return true;
    }
  return false;
}

//---------------------------------------------------------------------

inline void
remove_border_elt(const Edge_like& ordered_key)
{
  ordered_key.first->remove_border_edge((void*) &(*ordered_key.second));
}

//---------------------------------------------------------------------

inline bool is_ordered_border_elt(const void_Edge_like& e, 
				  IO_edge_type* &ptr)
{
  Vertex_handle v1 = (Vertex*) e.first;

  Next_border_elt* it12 =  v1->get_border_elt(e.second);
  if (it12 != NULL)
    {   
      ptr = &it12->second.first.second;
      return true;
    }
  return false;
}

inline void set_incidence_request(const Vertex_handle& v,
				  const criteria& value,
				  const Edge_like& e)
{
  void_Edge_like ve((void*) &*e.first, (void*) &*e.second);
  Incidence_request_elt incident_elt(value, ve);
  v->set_incidence_request(incident_elt);
}

//---------------------------------------------------------------------

inline bool is_interior_edge(const Edge_like& key)
  // pour gerer certaines aretes interieures: a savoir celle encore connectee au 
  // bord (en fait seule, les aretes interieures reliant 2 bords nous
  // interressent...)
{
  return (key.first->is_interior_edge((void*) &(*key.second))||
	  key.second->is_interior_edge((void*) &(*key.first)));
}

//---------------------------------------------------------------------

#ifndef NOLAZY

inline coord_type get_lazy_squared_radius(const Cell_handle& c)
{
  if (c->get_lazy_squared_radius() != NULL)
    return *(c->get_lazy_squared_radius());
  /*
  c->set_lazy_squared_radius
    (CGAL::squared_radius(convert()(c->vertex(0)->point()),
			  convert()(c->vertex(1)->point()),
			  convert()(c->vertex(2)->point()),
			  convert()(c->vertex(3)->point())));
  */
  c->set_lazy_squared_radius
    (CGAL::squared_radius(c->vertex(0)->point(),
			  c->vertex(1)->point(),
			  c->vertex(2)->point(),
			  c->vertex(3)->point()));
  return *(c->get_lazy_squared_radius());
}

inline D_Point get_lazy_circumcenter(const Cell_handle& c)
{
  if (c->get_lazy_circumcenter() != NULL)
    return *(c->get_lazy_circumcenter());
  /*
  c->set_lazy_circumcenter
    (CGAL::circumcenter(convert()(c->vertex(0)->point()),
			convert()(c->vertex(1)->point()),
			convert()(c->vertex(2)->point()),
			convert()(c->vertex(3)->point())));
  */
  c->set_lazy_circumcenter
    (CGAL::circumcenter(c->vertex(0)->point(),
			c->vertex(1)->point(),
			c->vertex(2)->point(),
			c->vertex(3)->point()));
  return *(c->get_lazy_circumcenter());
}

#endif //NOLAZY

//---------------------------------------------------------------------

inline Edge_IFacet inc_facet_circ(const Edge_IFacet& e)
{
  Cell_handle c = (Cell*) e.first.first;
  int i = e.second;
  int i1 = e.first.second, i2 = e.first.third;
  int i3 = (6 - e.second - i1 - i2);
  
  Cell_handle n = c->neighbor(i);
  int j1 = n->index(c->vertex(i1)), j2 = n->index(c->vertex(i2));
  int j =  n->index(c->vertex(i3));
  return Edge_IFacet(void_Edge((void*) &*n, j1, j2), j);  
}

//---------------------------------------------------------------------

inline Edge_IFacet dec_facet_circ(const Edge_IFacet& e)
{
  Cell_handle c = (Cell*) e.first.first;
  int i = e.second;
  int i1 = e.first.second, i2 = e.first.third;
  int i3 = (6 - e.second - i1 - i2);
  
  Cell_handle n = c->neighbor(i3);
  int j1 = n->index(c->vertex(i1)), j2 = n->index(c->vertex(i2));
  int j =  n->index(c->vertex(i));
  return Edge_IFacet(void_Edge((void*) &*n, j1, j2), j);  
}

//---------------------------------------------------------------------

inline bool
my_coplanar(const D_Point& p, const D_Point& q, 
	    const D_Point& r, const D_Point& s)
{
  coord_type qpx = q.x()-p.x();
  coord_type qpy = q.y()-p.y();
  coord_type qpz = q.z()-p.z();
  coord_type rpx = r.x()-p.x();
  coord_type rpy = r.y()-p.y();
  coord_type rpz = r.z()-p.z();
  coord_type spx = s.x()-p.x();
  coord_type spy = s.y()-p.y();
  coord_type spz = s.z()-p.z();

  coord_type den = CGAL::det3x3_by_formula(qpx,qpy,qpz,
					   rpx,rpy,rpz,
					   spx,spy,spz);
  return (CGAL_NTS abs(den) < eps_3);
}

//---------------------------------------------------------------------


inline bool
my_collinear(const D_Point& p, const D_Point& q, const D_Point& s)
{
  coord_type psx = p.x()-s.x();
  coord_type psy = p.y()-s.y();
  coord_type psz = p.z()-s.z();
  coord_type qsx = q.x()-s.x();
  coord_type qsy = q.y()-s.y();
  coord_type qsz = q.z()-s.z();
  coord_type rsx = psy*qsz-psz*qsy;
  coord_type rsy = psz*qsx-psx*qsz;
  coord_type rsz = psx*qsy-psy*qsx;

  coord_type den = CGAL::det3x3_by_formula(psx,psy,psz,
					   qsx,qsy,qsz,
					   rsx,rsy,rsz);

  return (CGAL_NTS abs(den) < eps_3);
}

//---------------------------------------------------------------------

inline void
visu_facet(const Cell_handle& c, const int& i)
{
  c->select_facet(i);
  _facet_number++;
}

//=====================================================================
//=====================================================================

void
show_selected_facets(CGAL::Geomview_stream &gv, const Triangulation_3& A)
{ 

  // Header.
  bool ascii_bak = gv.get_ascii_mode();
  bool raw_bak = gv.set_raw(true);
  _vh_number = 0;
  for (Finite_vertices_iterator v_it = A.finite_vertices_begin();
       v_it != A.finite_vertices_end();
       v_it++)
   if (!v_it->is_exterior())
     {
       _vh_number++;
     }

  gv.set_binary_mode();
  gv << "(geometry " << gv.get_new_id("object")
     << " {appearance {}{ OFF BINARY\n"
     << _vh_number << _facet_number << 0;

  CGAL::Unique_hash_map<Vertex*, int> vertex_index_map(-1, A.number_of_vertices());

  int count(0);
  for (Finite_vertices_iterator v_it = A.finite_vertices_begin();
       v_it != A.finite_vertices_end();
       v_it++){
    CGAL::Unique_hash_map<Vertex*, int>::Data& d = vertex_index_map[&(*v_it)];
    if ((!v_it->is_exterior()) && d == -1){
      d = count;
      count++;
      gv << convert()(v_it->point())  << " \n";
    }
  }

  for(Finite_facets_iterator f_it = A.finite_facets_begin(); 
      f_it != A.finite_facets_end(); 
      f_it++)
    {
      Cell_handle n, c = (*f_it).first;
      int ni, ci = (*f_it).second;
      n = c->neighbor(ci);
      ni = n->index(c);
      int i1, i2 ,i3;

      if (c->is_selected_facet(ci))
	{
	  i1 = (ci+1) & 3;
	  i2 = (ci+2) & 3;
	  i3 = (ci+3) & 3;
	  gv << 3;
	  gv << vertex_index_map[&(*c->vertex(i1))];
	  gv << vertex_index_map[&(*c->vertex(i2))];
	  gv << vertex_index_map[&(*c->vertex(i3))];
	  gv << 0; // without color.
	  // gv << 4 << drand48() << drand48() << drand48() << 1.0; // random
	  // color
	}

       if (n->is_selected_facet(ni))
	{
	  i1 = (ni+1) & 3;
	  i2 = (ni+2) & 3;
	  i3 = (ni+3) & 3;
	  gv << 3;
	  gv << vertex_index_map[&(*n->vertex(i1))];
	  gv << vertex_index_map[&(*n->vertex(i2))];
	  gv << vertex_index_map[&(*n->vertex(i3))];
	  gv << 0; // without color.
	  // gv << 4 << drand48() << drand48() << drand48() << 1.0; // random
	  // color 
	}
    }
  // Footer.
  gv << "}})";
      
  gv.set_raw(raw_bak);
  gv.set_ascii_mode(ascii_bak);
}

//=====================================================================
//=====================================================================

int border_counter(const Triangulation_3& A)
{
  int _border_count(0);
  for(Finite_edges_iterator e_it=A.finite_edges_begin();
      e_it!=A.finite_edges_end();
      e_it++)
    {
      Cell_handle c = (*e_it).first;
      int i1 = (*e_it).second, i2 = (*e_it).third;
      Edge_like key(c->vertex(i1), c->vertex(i2));
      Border_elt result;

      if (is_border_elt(key, result))
	_border_count++;
    }
  return _border_count;
}

//=====================================================================
#endif //UTILITIES_H
//=====================================================================

