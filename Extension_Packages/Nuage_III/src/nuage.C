
//=====================================================================
// selection de facettes dans Delaunay....
//=====================================================================
#include <CGAL/basic.h>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <strstream>
#include <cassert>
#include <vector>
#include <list>


// Kernel
//#include <CGAL/Kernel_checker.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Static_filters.h>
#include <CGAL/Cartesian_converter.h>

#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Triangulation_geom_traits_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_hierarchy_3.h>

#include <CGAL/squared_distance_3.h>

#include <CGAL/IO/Geomview_stream.h>
#include <CGAL/IO/Triangulation_geomview_ostream_3.h>

#include <CGAL/Timer.h>

#include <NUAGE/Local_selection_vertex_base_3.h>
#include <NUAGE/Local_selection_cell_base_3.h>

#include "Parse.C"

//=====================================================================
//=====================================================================

typedef double coord_type;
typedef double NT;

//struct Rep : public CGAL::Filtered_kernel<CGAL::Simple_cartesian<NT> > {};
struct Rep : public CGAL::Static_filters<CGAL::Simple_cartesian<NT> > {};

//---
// struct leda_real_NT_converter
// {
//     leda_real
//     operator()(const NT &a) const
//     {
//         return a.exact();
//     }
// };

// typedef leda_real ENT;
// typedef CGAL::Simple_cartesian<NT> F_Rep;
// typedef CGAL::Simple_cartesian<ENT> E_Rep;
// typedef CGAL::Cartesian_converter<F_Rep, E_Rep, leda_real_NT_converter > TMP_Conv;
// typedef CGAL::Kernel_checker<F_Rep, E_Rep, TMP_Conv > Rep;
//---


struct D_Rep : public CGAL::Simple_cartesian<coord_type> {};

//conversion de coord_type en Filtered_exact
struct coord_type_NT_converter
{
    coord_type
    operator()(const NT &a) const
    {
        return CGAL::to_double(a);
    }
};

//==========================================

typedef CGAL::Cartesian_converter<Rep, D_Rep, coord_type_NT_converter > convert;

// af: does this make sense when the two kernels are the same:
//struct convert {
//  Rep::Point_3& operator()(Rep::Point_3& p) {return p;}
//};


typedef Rep::Point_3  Point;
typedef Rep::Vector_3 Vector;
typedef Rep::Sphere_3  Sphere;
typedef Rep::Segment_3  Segment;
typedef Rep::Ray_3  Ray;
typedef Rep::Line_3  Line;
typedef Rep::Triangle_3  Triangle;

typedef D_Rep::Point_3  D_Point;
typedef D_Rep::Vector_3 D_Vector;
typedef D_Rep::Sphere_3  D_Sphere;
typedef D_Rep::Segment_3  D_Segment;
typedef D_Rep::Ray_3  D_Ray;
typedef D_Rep::Line_3  D_Line;
typedef D_Rep::Triangle_3  D_Triangle;


typedef Rep Gt; // af
typedef Local_selection_vertex_base_3<Gt> my_Vb;
typedef CGAL::Triangulation_hierarchy_vertex_base_3<my_Vb> Vb;
typedef Local_selection_cell_base_3<Gt> Fb;


typedef CGAL::Triangulation_data_structure_3<Vb,Fb> Tds;
typedef CGAL::Delaunay_triangulation_3<Gt,Tds> Delaunay_Triangulation_3;
typedef CGAL::Triangulation_hierarchy_3<Delaunay_Triangulation_3> Triangulation_3;

typedef Triangulation_3::Cell  Cell;
typedef Triangulation_3::Vertex Vertex;
typedef Triangulation_3::Edge Edge;
typedef Triangulation_3::Facet Facet;
typedef Triangulation_3::Cell_handle  Cell_handle;
typedef Triangulation_3::Vertex_handle Vertex_handle;

typedef Triangulation_3::Cell_circulator  Cell_circulator;
typedef Triangulation_3::Facet_circulator Facet_circulator;

typedef Triangulation_3::Locate_type Locate_type;

typedef Triangulation_3::Finite_cells_iterator  Finite_cells_iterator;
typedef Triangulation_3::Finite_facets_iterator Finite_facets_iterator;
typedef Triangulation_3::Finite_vertices_iterator  Finite_vertices_iterator;
typedef Triangulation_3::Finite_edges_iterator  Finite_edges_iterator;

typedef Triangulation_3::All_cells_iterator  All_cells_iterator;
typedef Triangulation_3::All_facets_iterator All_facets_iterator;
typedef Triangulation_3::All_vertices_iterator  All_vertices_iterator;
typedef Triangulation_3::All_edges_iterator  All_edges_iterator;

typedef Vb::void_Edge void_Edge;
typedef Vb::Edge_IFacet Edge_IFacet;
typedef Vb::IO_edge_type IO_edge_type;
typedef Vb::criteria criteria;
typedef Vb::Radius_edge_type Radius_edge_type;
typedef Vb::Border_elt Border_elt;
typedef Vb::Next_border_elt Next_border_elt;
typedef Vb::Radius_ptr_type Radius_ptr_type;

typedef Vb::Incidence_request_iterator Incidence_request_iterator;
typedef Vb::void_Edge_like void_Edge_like;
typedef Vb::Incidence_request_elt Incidence_request_elt;

typedef std::pair< Vertex_handle, Vertex_handle > Edge_like;
typedef CGAL::Triple< Vertex_handle, Vertex_handle, Vertex_handle > Facet_like;

typedef std::list< Facet_like > Additional_facets_list_type;
typedef Additional_facets_list_type::iterator Additional_facets_list_iterator;

typedef std::multimap< criteria, IO_edge_type*, 
                       std::less<criteria> > Ordered_border_type;
typedef Ordered_border_type::iterator Ordered_border_iterator;

enum Validation_case {not_valid, not_valid_connecting_case, final_case,
		      ear_case, exterior_case, connecting_case};

//=====================================================================
//=====================================================================

Ordered_border_type* _ordered_border = new Ordered_border_type();
Additional_facets_list_type* 
       _additional_facets_list = new Additional_facets_list_type();
int _number_of_border(1);
//SLIVER_ANGULUS represente la qualite d'echantillonnage de la surface 
const coord_type SLIVER_ANGULUS = .86;
// DELTA represente la qualite d'echantillonnage du bord
coord_type DELTA;
coord_type K, min_K;
const coord_type eps(1e-7);
const coord_type inv_eps_2(coord_type(1)/(eps*eps)); // 1/(eps^2)
const coord_type eps_3(eps*eps*eps); // test de ^3 donc points tel 1e-7 soit petit
const criteria STANDBY_CANDIDATE(3);
const criteria STANDBY_CANDIDATE_BIS(STANDBY_CANDIDATE+1);
const criteria NOT_VALID_CANDIDATE(STANDBY_CANDIDATE+2);

//---------------------------------------------------------------------
// pour desactiver le dump dans geomview :
#define BLIND
// pour desactiver l'evaluation paresseuse...
// #define NOLAZY

CGAL::Timer t1;

//---------------------------------------------------------------------
//Pour une visu correcte
//pour retenir les facettes selectionnees
int _vh_number(0);
int _facet_number(0);

//---------------------------------------------------------------------
//Pour le post traitement
int _postprocessing_cont(0);
int _A_size_before_postprocessing(0);

//=====================================================================
//=====================================================================

#include <NUAGE/utilities.h>
#include <NUAGE/iofile_manipulator.h>

//=====================================================================
//=====================================================================
inline
coord_type get_smallest_radius_delaunay_sphere(const Triangulation_3 & A,
					       const Cell_handle& c,
					       const int& index)
{
  int i1, i2, i3;

  Cell_handle n = c->neighbor(index);
  // evaluation paresseuse...
  coord_type value = c->get_smallest_radius(index);
  if ((value >= 0)&&(n->get_smallest_radius(n->index(c)) == value))
    return value;
  /*
  D_Point cp0 = convert()(c->vertex(index)->point());
  D_Point cp1 = convert()(c->vertex((index+1) & 3)->point());
  D_Point cp2 = convert()(c->vertex((index+2) & 3)->point());
  D_Point cp3 = convert()(c->vertex((index+3) & 3)->point());

  D_Point np0=  convert()(n->vertex(0)->point());
  D_Point np1 = convert()(n->vertex(1)->point());
  D_Point np2 = convert()(n->vertex(2)->point());
  D_Point np3 = convert()(n->vertex(3)->point());
  */
  const D_Point& cp0 = c->vertex(index)->point();
  const D_Point& cp1 = c->vertex((index+1) & 3)->point();
  const D_Point& cp2 = c->vertex((index+2) & 3)->point();
  const D_Point& cp3 = c->vertex((index+3) & 3)->point();

  const D_Point& np0 = n->vertex(0)->point();
  const D_Point& np1 = n->vertex(1)->point();
  const D_Point& np2 = n->vertex(2)->point();
  const D_Point& np3 = n->vertex(3)->point();

  bool c_is_plane(my_coplanar(cp0, cp1, cp2, cp3));
  bool n_is_plane(my_coplanar(np0, np1, np2, np3));

  bool c_is_infinite(A.is_infinite(c));
  bool n_is_infinite(A.is_infinite(n));
  if ((c_is_plane && n_is_plane)||
      (c_is_plane && n_is_infinite)||
      (n_is_plane && c_is_infinite)||
      my_collinear(cp1, cp2, cp3))
    value = HUGE_VAL;
  else
    {
      if (c_is_infinite||n_is_infinite||c_is_plane||n_is_plane)
	{ 
	  int ind;
	  Cell_handle cc;
	  if(c_is_infinite||c_is_plane)
	    {
	      cc = n;
	      ind = n->index(c);
	    }
	  else
	    {
	      cc = c;
	      ind = index;
	    }
	  i1 = (ind+1) & 3;
	  i2 = (ind+2) & 3;
	  i3 = (ind+3) & 3;
	  /*
	  D_Point pp0 = convert()(cc->vertex(ind)->point());
	  D_Point pp1 = convert()(cc->vertex(i1)->point());
	  D_Point pp2 = convert()(cc->vertex(i2)->point());
	  D_Point pp3 = convert()(cc->vertex(i3)->point());
	  */
	  const D_Point& pp0 = cc->vertex(ind)->point();
	  const D_Point& pp1 = cc->vertex(i1)->point();
	  const D_Point& pp2 = cc->vertex(i2)->point();
	  const D_Point& pp3 = cc->vertex(i3)->point();

	  D_Sphere facet_sphere(pp1, pp2, pp3);
	  if (CGAL::squared_distance(facet_sphere.center(), pp0) <
	      facet_sphere.squared_radius())
	    {
#ifndef NOLAZY
	      value = get_lazy_squared_radius(cc);
#endif //NOLAZY
#ifdef NOLAZY
	      value = CGAL::squared_radius(pp0, pp1, pp2, pp3);
#endif //NOLAZY
	    }
	  else
	    value = facet_sphere.squared_radius();
	}
      else
	{
	  D_Point cc, cn;
#ifndef NOLAZY
	  cc = get_lazy_circumcenter(c);
	  cn = get_lazy_circumcenter(n);
#endif //NOLAZY
#ifdef NOLAZY
	  cc = CGAL::circumcenter(cp0, cp1, cp2, cp3);
	  cn = CGAL::circumcenter(np0, np1, np2, np3);
#endif //NOLAZY
	  // calcul de la distance de cp1 au segment dual cc, cn...
	  D_Vector V(cc - cn), Vc(cc - cp1), Vn(cp1 - cn);
	  coord_type ac(V * Vc), an(V * Vn), norm_V(V * V);
	  if ((ac > 0) && (an > 0))
	    {
	      value = (Vc*Vc) - ac*ac/norm_V;
	      if ((value < 0)||(norm_V > inv_eps_2))
		value = CGAL::squared_radius(cp1, cp2, cp3);
	    }
	  else
	    {
	      if (ac <= 0)
		value = CGAL::squared_distance(cc, cp1);
	      else // (an <= 0)
		value = CGAL::squared_distance(cn, cp1);
	    }
	}
    }
  // stockage des valeurs deja calculee...
  c->set_smallest_radius(index, value);
  n->set_smallest_radius(n->index(c), value);

  return value;
}

//---------------------------------------------------------------------

Radius_edge_type compute_value(const Edge_IFacet& e, const Triangulation_3& A)
{
  Cell_handle c = (Cell*) e.first.first;
  int i = e.second;
  int i1 = e.first.second, i2 = e.first.third;
  int i3 = 6 - e.second - i1 - i2;

  Edge_IFacet e_it = e, predone = dec_facet_circ(e);
  Cell_handle c_predone = (Cell*) predone.first.first;

  coord_type min_valueP = NOT_VALID_CANDIDATE, min_valueA = HUGE_VAL;
  Facet min_facet, min_facetA;
  bool border_facet(false);

  coord_type pscal;//, prec_pliure = e.third;

/*
  D_Point p1 = convert()(c->vertex(i1)->point());
  D_Point p2 = convert()(c->vertex(i2)->point());
  D_Point pc = convert()(c->vertex(i3)->point());
*/
  const D_Point& p1 = c->vertex(i1)->point();
  const D_Point& p2 = c->vertex(i2)->point();
  const D_Point& pc = c->vertex(i3)->point();
  
  D_Vector P2P1 = p1-p2, P2Pn, PnP1;

  D_Vector v2, v1 = CGAL::cross_product(pc-p2, P2P1);

  coord_type norm, norm1 = v1*v1;
  coord_type norm12 = P2P1*P2P1;
  //int count(0); 

  e_it = inc_facet_circ(e_it);
  bool succ_start(true);
  
  do
    {
      Cell_handle neigh = (Cell*) e_it.first.first;
      Facet facet_it(neigh, e_it.second);

      if (!A.is_infinite(facet_it))
	// &&!CGAL::collinear(p1, p2, pc) en principe inutile car HUGE_VAL ???
	{
	  int n_ind = facet_it.second;
	  int n_i1 = e_it.first.second;
	  int n_i2 = e_it.first.third;
	  int n_i3 = 6 - n_ind - n_i1 - n_i2; 

	  coord_type tmp = get_smallest_radius_delaunay_sphere(A, neigh, n_ind);
	      
	  // 	  bool is_on_same_border
	  // 	    (neigh->vertex(n_i3)->is_on_border(result12.second));

	  Edge_like el1(neigh->vertex(n_i1),neigh->vertex(n_i3)),
	    el2(neigh->vertex(n_i2),neigh->vertex(n_i3));

	  // si on veut ne s'autoriser que le meme bord pour vni3
	  //         if ((neigh->vertex(n_i3)->is_exterior() || is_on_same_border)&&
	  // si on veut pouvoir connecter des bords differents
	  //         if (neigh->vertex(n_i3)->not_interior()&&
	  
	  if ((tmp != HUGE_VAL)&&
	      neigh->vertex(n_i3)->not_interior()&&
	      (!is_interior_edge(el1))&&(!is_interior_edge(el2)))
	    {
	      D_Point pn = convert()(neigh->vertex(n_i3)->point());
	 
	      P2Pn = pn-p2;
	      v2 = CGAL::cross_product(P2P1,P2Pn);

	      //pas necessaire de normer pour un bon echantillon: 
	      //            on peut alors tester v1*v2 >= 0
	      norm =  CGAL::sqrt(norm1 * (v2*v2));
	      pscal = v1*v2;
	      //SLIVER_ANGULUS represente la qualite d'echantillonnage de la
	      //surface 
	      bool sliver_facet = ((succ_start || (neigh == c_predone))&& 
				   (pscal <= -SLIVER_ANGULUS*norm));
	      
	      if (succ_start) succ_start = false;

	      if (!sliver_facet)
		{
		  if (tmp < min_valueA)
		    {
		      PnP1 = p1-pn;
		      // DELTA represente la qualite d'echantillonnage du bord
		      border_facet = !((P2P1*P2Pn >= 
					-DELTA*CGAL::sqrt(norm12*(P2Pn*P2Pn)))&&
				       (P2P1*PnP1 >= 
					-DELTA*CGAL::sqrt(norm12*(PnP1*PnP1))));

		      min_facetA = facet_it; 
		      min_valueA = tmp;
		      min_valueP = pscal/norm;
		    }
		}
	    }
	}
      //count++;
      e_it = inc_facet_circ(e_it);
    }
  while(Cell_handle((Cell*) e_it.first.first) != c);

  criteria value;

  if ((min_valueA == HUGE_VAL) || border_facet) // bad facets case
    { 
      //std::cout << "aucune facette candidate parmi " << count-1 << std::endl;
      min_facet = Facet(c, i); // !!! sans aucune signification....
      value = NOT_VALID_CANDIDATE; // Attention a ne pas inserer dans PQ
    }
  else
    {
      min_facet = min_facetA; 

      //si on considere seulement la pliure value appartient a [0, 2]
      //value = coord_type(1) - min_valueP;
 
      // si la pliure est bonne on note suivant le alpha sinon on prend en compte la 
      // pliure seule... pour discriminer entre les bons slivers...
      // si on veut discriminer les facettes de bonnes pliures plus finement
      // alors -(1+1/min_valueA) app a [-inf, -1]
      // -min_valueP app a [-1, 1]

      if (min_valueP > SLIVER_ANGULUS)
	value = -(coord_type(1) + coord_type(1)/min_valueA);
      else
	{
	  //on refuse une trop grande non-uniformite
	  if (min_valueA <= K*get_smallest_radius_delaunay_sphere(A, c, i))
	    value = - min_valueP;
	  else
	    {
	      value = STANDBY_CANDIDATE; // tres mauvais candidat mauvaise pliure
	      // + grand alpha... a traiter plus tard....
	      min_K = 
		std::min(min_K,
			 min_valueA/get_smallest_radius_delaunay_sphere(A, c, i));
	    }
	}
    }
  
  Cell_handle n = min_facet.first;
  int ni1 = n->index(c->vertex(i1)), ni2 = n->index(c->vertex(i2));

  return 
    Radius_edge_type(value, IO_edge_type(e, Edge_IFacet
					 (void_Edge((void*) &(*n), ni1, ni2),
					  min_facet.second)));
}

//=====================================================================
//=====================================================================

bool
init(const Triangulation_3& A, const bool& re_init)
{
  Facet min_facet;
  coord_type min_value = HUGE_VAL;
  int i1, i2, i3;

  if (!re_init)
    for(Finite_facets_iterator facet_it=A.finite_facets_begin(); 
	facet_it!=A.finite_facets_end(); 
	facet_it++)
      {
	coord_type value = get_smallest_radius_delaunay_sphere(A, (*facet_it).first,
							       (*facet_it).second);
	if (value < min_value)
	  {
	    min_facet = *facet_it;
	    min_value = value;
	  }
      }
  else //if (re_init)
     for(Finite_facets_iterator facet_it=A.finite_facets_begin(); 
	 facet_it!=A.finite_facets_end(); 
	 facet_it++)
       {
	 Cell_handle c = (*facet_it).first;
	 int index = (*facet_it).second;	 
	 if (c->vertex((index+1) & 3)->is_exterior())
	   if (c->vertex((index+2) & 3)->is_exterior())
	     if (c->vertex((index+3) & 3)->is_exterior())
	      {
		coord_type value = get_smallest_radius_delaunay_sphere(A, c, index);

		if (value < min_value)
		  {
		    min_facet = *facet_it;
		    min_value = value;
		  }
	      }
       }

  if (min_value != HUGE_VAL)
    {
      Cell_handle c_min = min_facet.first;
      int ind = min_facet.second;
      i1 = (ind+1) & 3;
      i2 = (ind+2) & 3;
      i3 = (ind+3) & 3;

      Radius_edge_type e12, e23, e31;

      e12 = compute_value(Edge_IFacet(void_Edge((void*) &(*c_min), i1, i2), ind), A);
      e23 = compute_value(Edge_IFacet(void_Edge((void*) &(*c_min), i2, i3), ind), A);
      e31 = compute_value(Edge_IFacet(void_Edge((void*) &(*c_min), i3, i1), ind), A);

      IO_edge_type* p12 = set_border_elt(c_min->vertex(i1), c_min->vertex(i2),
					 Border_elt(e12, _number_of_border));
      IO_edge_type* p23 = set_border_elt(c_min->vertex(i2), c_min->vertex(i3),
					 Border_elt(e23, _number_of_border));
      IO_edge_type* p31 = set_border_elt(c_min->vertex(i3), c_min->vertex(i1),
					 Border_elt(e31, _number_of_border));
  
      c_min->vertex(i1)->inc_mark();
      c_min->vertex(i2)->inc_mark();
      c_min->vertex(i3)->inc_mark();
//       if (e12.first < NOT_VALID_CANDIDATE)
	_ordered_border->insert(Radius_ptr_type (e12.first, p12));
//       if (e23.first < NOT_VALID_CANDIDATE)
	_ordered_border->insert(Radius_ptr_type (e23.first, p23));
//       if (e31.first < NOT_VALID_CANDIDATE)
       _ordered_border->insert(Radius_ptr_type (e31.first, p31));

      // Pour une visu correcte_e_it_bis
      visu_facet(c_min, ind);
      return true;
    }
  return false;
}

//---------------------------------------------------------------------
// test de reciprocite avant de recoller une oreille anti-singularite
inline int
test_merge(const Triangulation_3& A,
	   const Edge_like& ordered_key, const Border_elt& result, 
	   const Vertex_handle& v, const coord_type& ear_alpha)
{
  Edge_IFacet Ifacet = result.first.second.first;
  /*
  D_Point p1 = convert()((ordered_key.first)->point());
  D_Point p2 = convert()((ordered_key.second)->point());
  D_Point pc = convert()(v->point());
  */
  
  const D_Point& p1 = (ordered_key.first)->point();
  const D_Point& p2 = (ordered_key.second)->point();
  const D_Point& pc = v->point();

  Cell_handle neigh = (Cell*) Ifacet.first.first;
  int n_ind = Ifacet.second;
  int n_i1 = Ifacet.first.second;
  int n_i2 = Ifacet.first.third;
  int n_i3 = (6 - n_ind - n_i1 - n_i2);

  D_Point pn = convert()(neigh->vertex(n_i3)->point());
  D_Vector v1 = CGAL::cross_product(pc-p2,p1-p2),
    v2 = CGAL::cross_product(p1-p2,pn-p2);
  coord_type norm = CGAL::sqrt((v1*v1)*(v2*v2));

  if (v1*v2 > SLIVER_ANGULUS*norm)    
    return 1; // label bonne pliure sinon:

  if (ear_alpha <= K*get_smallest_radius_delaunay_sphere(A, neigh, n_ind))
    return 2; // label alpha coherent...
    
  return 0; //sinon oreille a rejeter...
}

//---------------------------------------------------------------------
// test de reciprocite avant de recoller une oreille 
// controle de la propagation par rampe
// inline bool
// test_merge_ear(const Edge_like& ordered_key, const Border_elt& result, 
// 	       const Vertex_handle& v, const coord_type& ear_alpha)
// {
//   Edge_IFacet Ifacet = result.first.second.first;
//   D_Point p1 = convert()((ordered_key.first)->point());
//   D_Point p2 = convert()((ordered_key.second)->point());
//   D_Point pc = convert()(v->point());
  
//   Cell_handle neigh = (Cell*) Ifacet.first.first;
//   int n_ind = Ifacet.second;
//   int n_i1 = Ifacet.first.second;
//   int n_i2 = Ifacet.first.third;
//   int n_i3 = (6 - n_ind - n_i1 - n_i2);

//   D_Point pn = convert()(neigh->vertex(n_i3)->point());
//   D_Vector v1 = CGAL::cross_product(pc-p2,p1-p2),
//     v2 = CGAL::cross_product(p1-p2,pn-p2);

//   if (v1*v2 > 0)      
//     return true;

//   // si ca courbe plus que pi/2 il faut un tres bon echantillon
//   if (ear_alpha <= K*neigh->get_smallest_radius(n_ind))
//     return true;
    
//   return false;
// }

//---------------------------------------------------------------------

void
_ordered_map_erase(const criteria& value, const IO_edge_type* pkey)
{
  int number_of_conflict = _ordered_border->count(value);  
  int verif(0); 
  if (number_of_conflict == 1)
    {
      _ordered_border->erase(_ordered_border->find(value));
      verif++;
    }

  if (number_of_conflict > 1)
    {
      Ordered_border_iterator elt_it =
	_ordered_border->find(value);
      // si ca foire jamais on peut s'areter des que l'elt 
      // est trouve!!! 
      for(int jj=0; (jj<number_of_conflict)&&(verif<1); jj++)
	{	  
	  if (((long) elt_it->second) == ((long) pkey))
	    {
	      _ordered_border->erase(elt_it);
	      verif++;
	    }
	  elt_it++;
	}
    }

  if (verif > 1)
    {
      std::cerr << "+++Problem with key: " << value << 
	" containing " << number_of_conflict << " elts." << std::endl; 
      std::cerr << "   " << verif << 
	" elts removed from _ordered_border" << std::endl;
    }
}

//---------------------------------------------------------------------

void
force_merge(const Edge_like& ordered_key, const Border_elt& result)
{
  //  Border_map_iterator bord_it = _border_map.find(key);

  criteria value = result.first.first;
  IO_edge_type* pkey = get_border_IO_elt(ordered_key.first, ordered_key.second);

  _ordered_map_erase(value, pkey);

  remove_border_elt(ordered_key);
}

//---------------------------------------------------------------------

void dequeue_incidence_request(const Vertex_handle& v)
{
  if (v->is_incidence_requested())
    {
      for(Incidence_request_iterator v_it = v->incidence_request_begin();
	  v_it != v->get_incidence_request_end();
	  v_it++)
	{
	  IO_edge_type* ptr;
	  
	  if (is_ordered_border_elt(v_it->second, ptr))
	    _ordered_border->insert(Radius_ptr_type(v_it->first, ptr));
	}
      v->erase_incidence_request();
    }
}

//---------------------------------------------------------------------

bool
try_to_close_border(IO_edge_type* /*p*/)
{
  //=================== Pas une mauvaise idee, MAIS =================
  // ATTENTION C'est incompatible avec le post-traitement !!!

  // pour l'instant on ferme juste les triangles qui ne sont pas dans
  // Delaunay, des que l'on en trouve un...


//   Edge_IFacet e_Ifacet = p->first;
//   Cell_handle c = (Cell*) e_Ifacet.first.first;
//   Vertex_handle 
//     v1 = c->vertex(e_Ifacet.first.second),
//     v2 = c->vertex(e_Ifacet.first.third);

//   Edge_like el(v1, v2);
//   Border_elt result;
//   is_border_elt(el, result);
  
//   Next_border_elt* succ1 = el.second->get_next_on_border(result.second);
//   Vertex_handle  v_succ1 = (Vertex*) succ1->first;
//   Next_border_elt* succ2 = v_succ1->get_next_on_border(result.second);
//   Vertex_handle  v_succ2 = (Vertex*) succ2->first;
  
//   if (v_succ2 == el.first)
//     {
//       // dans ce cas on a a faire a un contour a trois cote qui n'a pas ete
//       // trouve comme candidat a fermer... certainement pas dans Delaunay...
//       remove_border_elt(el);
//       force_merge(Edge_like(el.second, v_succ1), succ1->second);
//       force_merge(Edge_like(v_succ1, el.first), succ2->second);
//       el.first->dec_mark();
//       el.second->dec_mark();
//       v_succ1->dec_mark();
//       // marquer la facette en question pour l'affichage...
//       _facet_number++;
//       _additional_facets_list->push_back(Facet_like(el.second, v_succ1, el.first));
//       return true;
//     }
  return false;
}


//---------------------------------------------------------------------

void
merge_ear(const Edge_like& ordered_el1, const Border_elt& result1, 
	  const Edge_like& ordered_key,
	  const Vertex_handle& v1, const Vertex_handle& v2,
	  const Edge_IFacet& edge_Ifacet_2, const Triangulation_3& A)
{  
  remove_border_elt(ordered_key);
  force_merge(ordered_el1, result1);
  
  Radius_edge_type e2 = compute_value(edge_Ifacet_2, A);

  IO_edge_type* p2;
  if (ordered_el1.first == v1)
    p2 = set_border_elt(v2, ordered_el1.second,
			Border_elt(e2,result1.second));
  else
    p2 = set_border_elt(ordered_el1.first, v2,
			Border_elt(e2,result1.second));

  v1->dec_mark();

  // if e2 contain HUGE_VAL there is no candidates to
  // continue: compute_value is not valid...

//   if (e2.first < NOT_VALID_CANDIDATE)
    _ordered_border->insert(Radius_ptr_type(e2.first, p2));
//   else
//     try_to_close_border(p2);

  //depiler les eventuelles requettes de connections avortees... zones etoilees, 
  //en effet le bord a change donc on peut peut etre maintenant.
  dequeue_incidence_request(v2);

  if (ordered_el1.first == v1)
    dequeue_incidence_request(ordered_el1.second);
  else
    dequeue_incidence_request(ordered_el1.first);
}

//---------------------------------------------------------------------

void
border_extend(const Edge_like& ordered_key, const Border_elt& result12, 
	      const Vertex_handle& v1, const Vertex_handle& v2,
	      const Vertex_handle& v3,
	      const Radius_edge_type& e1, const Radius_edge_type& e2,
	      IO_edge_type* &p1, IO_edge_type* &p2)
{
  remove_border_elt(ordered_key);

  //depiler v3 avant de le mettre a jour... pour reperer s'il est sur un bord
  if (v3->number_of_incident_border() > 0)
    dequeue_incidence_request(v3);

  if (ordered_key.first == v1)
    {
      p1 = set_border_elt(v1, v3, Border_elt(e1,result12.second));
      p2 = set_border_elt(v3, v2, Border_elt(e2,result12.second));
    }
  else
    {
      p2 = set_border_elt(v2, v3, Border_elt(e2,result12.second));
      p1 = set_border_elt(v3, v1, Border_elt(e1,result12.second));
    }

  v3->inc_mark(); 


  //depiler les eventuelles requettes de connections avortees... zones etoilees, 
  //en effet le bord a change donc on peut peut etre maintenant.
  dequeue_incidence_request(v1);
  dequeue_incidence_request(v2);
}

//=====================================================================
//=====================================================================

Validation_case
validate(const Edge_IFacet& edge_Efacet, const Triangulation_3& A, 
	 const criteria& value)
{
  int i = (6 - edge_Efacet.second 
	   - edge_Efacet.first.second
	   - edge_Efacet.first.third);
  Cell_handle c = (Cell*) edge_Efacet.first.first;

//   coord_type candidate_alpha =  c->get_smallest_radius(edge_Efacet.second);
  //  coord_type pre_pliure = edge_Efacet.third;

//  if ((c->vertex(i)->not_interior() > 0)&&(value > K*alpha_max))
//     return not_valid;

  Vertex_handle 
    v1 = c->vertex(edge_Efacet.first.second),
    v2 = c->vertex(edge_Efacet.first.third);
      	      
  Edge_like ordered_el1(c->vertex(i), v1);
  Edge_like ordered_el2(c->vertex(i), v2);
  Border_elt result1, result2, result12;

  Edge_like ordered_key(v1,v2);

  if (!is_border_elt(ordered_key, result12))
    std::cerr << "+++probleme coherence bord <validate>" << std::endl;

  bool is_border_el1 = is_border_elt(ordered_el1, result1),
    is_border_el2 = is_border_elt(ordered_el2, result2);

  //  bool is_on_same_border (c->vertex(i)->is_on_border(result12.second));

  Radius_edge_type e1, e2;

  if (c->vertex(i)->not_interior())
    {
      if ((!is_interior_edge(ordered_el1))&&
	  (!is_interior_edge(ordered_el2)))
	{ 
	  //toujours utile meme avec l'essai de try_to_close_border avant
	  //validate pour la resolution de singularite par oreille qui elle
	  //doit etre dans Delaunay.
	  if (is_border_el1&&is_border_el2)
	    { 
	      remove_border_elt(ordered_key);			  
	      force_merge(ordered_el1, result1);
	      force_merge(ordered_el2, result2);

	      v1->dec_mark();
	      v2->dec_mark();
	      c->vertex(i)->dec_mark();

	      //Pour une visu correcte
	      visu_facet(c, edge_Efacet.second);	

	      return final_case;
	    }

	  //--------------------------------------------------------------------- 
	  //on peut alors marquer v1 et on pourrait essayer de merger 
	  //sans faire de calcul inutile???
	  if (is_border_el1)
	    {
// 	      if (test_merge_ear(ordered_el1, result1, v2, candidate_alpha)&&
// 		  (result12.second==result1.second))// force a merger
// 		{
		  Edge_IFacet edge_Ifacet_2(void_Edge((void*) &(*c), i, edge_Efacet.first.third), 
					    edge_Efacet.second);
		  merge_ear(ordered_el1, result1, 
			    ordered_key, v1, v2, edge_Ifacet_2, A);

		  //Pour une visu correcte
		  visu_facet(c, edge_Efacet.second);

		  return ear_case;
// 		}
// 	      return not_valid;
	    }

	  //---------------------------------------------------------------------
	  //idem pour v2
	  if (is_border_el2)
	    {
// 	      if (test_merge_ear(ordered_el2, result2, v1, candidate_alpha)&&
// 		  (result12.second==result2.second))// force a merger
// 		{
		  Edge_IFacet edge_Ifacet_1(void_Edge((void*) &(*c), i, edge_Efacet.first.second), 
					    edge_Efacet.second);
		  merge_ear(ordered_el2, result2, 
			    ordered_key, v2, v1, edge_Ifacet_1, A);

		  //Pour une visu correcte
		  visu_facet(c, edge_Efacet.second);

		  return ear_case;
// 		}
// 	      return not_valid;
	    }
	    

	  //---------------------------------------------------------------------
	  if ((!is_border_el1)&&(!is_border_el2))
	    { 
	      // si on veut s'interdir de spliter un bord (pelure d'orange....) 
	      // seulement c->vertex(i)->is_exterior()
	      // pour s'autoriser des split de bord surface a bord->sphere ou Moebius...
	      // alors || is_on_same_border:
	      //       if (c->vertex(i)->is_exterior() || is_on_same_border) 
	      // pour passer au tore (changementde type de topologie)
	      // recoller deux bord different...
	      //       if (c->vertex(i)->not_interior() deja teste en haut	      
	      
	      if(c->vertex(i)->is_exterior())
		{		  
		  Edge_IFacet edge_Ifacet_1(void_Edge((void*) &(*c), i, edge_Efacet.first.second), 
					    edge_Efacet.second);
		  
		  Edge_IFacet edge_Ifacet_2(void_Edge((void*) &(*c), i, edge_Efacet.first.third), 
					    edge_Efacet.second);
		  e1 = compute_value(edge_Ifacet_1, A);
		  e2 = compute_value(edge_Ifacet_2, A);
		  
		  IO_edge_type* p1;  
		  IO_edge_type* p2;		  

		  border_extend(ordered_key, result12, 
				v1, v2, c->vertex(i),
				e1, e2, p1, p2);
	     
		  // if e1 contain HUGE_VAL there is no candidates to
		  // continue: compute_value is not valid...

// 		   if (e1.first < NOT_VALID_CANDIDATE)
		     _ordered_border->insert(Radius_ptr_type(e1.first, p1));
// 		   else
// 		     try_to_close_border(p1);
		  // if e2 contain HUGE_VAL there is no candidates to
		  // continue: compute_value is not valid...

// 		   if (e2.first < NOT_VALID_CANDIDATE)
		     _ordered_border->insert(Radius_ptr_type(e2.first, p2));
// 		   else
// 		     try_to_close_border(p2);
		  //Pour une visu correcte
		  visu_facet(c, edge_Efacet.second);

		  return exterior_case;
		}
	      else // c->vertex(i) is a border point (and now there's only 1
		// border incident to a point... _mark<1 even if th orientation
		// may be such as one vh has 2 successorson the same border...
		{		 
		  // a ce niveau on peut tester si le recollement se fait en
		  // maintenant la compatibilite d'orientation des bords (pour
		  // surface orientable...) ou si elle est brisee...
		  Edge_IFacet edge_Ifacet_1(void_Edge((void*) &(*c), i, edge_Efacet.first.second), 
					    edge_Efacet.second);
		  Edge_IFacet edge_Ifacet_2(void_Edge((void*) &(*c), i, edge_Efacet.first.third), 
					    edge_Efacet.second);

		  e1 = compute_value(edge_Ifacet_1, A);
		  e2 = compute_value(edge_Ifacet_2, A);

		  if ((e1.first >= STANDBY_CANDIDATE)&&(e2.first >= STANDBY_CANDIDATE)) 
		    return not_valid_connecting_case;

		  // vu compute value: les candidats oreilles fournis sont sans
		  // aretes interieures et le sommet oppose n'est pas non plus interieur
		  Edge_IFacet ear1 = e1.second.second;
		  Edge_IFacet ear2 = e2.second.second;

		  int ear1_i = (6 - ear1.second 
				- ear1.first.second
				- ear1.first.third);
		  Cell_handle ear1_c = (Cell*) ear1.first.first;
		  Border_elt result_ear1;

		  int ear2_i = (6 - ear2.second 
				- ear2.first.second
				- ear2.first.third);
		  Cell_handle ear2_c = (Cell*) ear2.first.first;
		  Border_elt result_ear2;

		  Edge_like ear1_e, ear2_e;
		  // pour maintenir la reconstruction d'une surface orientable :
		  // on verifie que les bords se recollent dans des sens opposes
		  if (ordered_key.first==v1)
		    {
		     ear1_e = Edge_like(c->vertex(i), ear1_c ->vertex(ear1_i));
		     ear2_e = Edge_like(ear2_c ->vertex(ear2_i), c->vertex(i));
		    }
		  else 
		    {
		      ear1_e = Edge_like(ear1_c ->vertex(ear1_i), c->vertex(i));
		      ear2_e = Edge_like(c->vertex(i), ear2_c ->vertex(ear2_i));
		    }
		    
		  //maintient la surface orientable
		  bool is_border_ear1 = is_ordered_border_elt(ear1_e, result_ear1);		  
		  bool is_border_ear2 = is_ordered_border_elt(ear2_e, result_ear2);
		  bool ear1_valid(false), ear2_valid(false);
		  //version sans controle d'orientabilite
// 		  bool is_border_ear1 = is_border_elt(ear1_e, result_ear1);		  
// 		  bool is_border_ear2 = is_border_elt(ear2_e, result_ear2);
		  if (is_border_ear1&&(e1.first < STANDBY_CANDIDATE)&&
		      (e1.first <=  value)&&
		      (result12.second==result_ear1.second))
		    {
		      ear1_valid = test_merge(A, ear1_e, result_ear1, v1,
					      get_smallest_radius_delaunay_sphere(A, ear1_c, 
										  ear1.second));
		      // si on veut etre plus restrictif
		      // et exiger au moins une bonne pliure:
// 		      int test_merge_ear1 = 
// 			test_merge(ear1_e, result_ear1, v1,
// 				   ear1_c->get_smallest_radius(ear1.second));	      
// 		      ear1_valid = (test_merge_ear1&&(e1.first < -1)&&
// 				    ((value < -1)||(test_merge_ear1 == 1)));
		    }
		  
		  if (is_border_ear2&&(e2.first < STANDBY_CANDIDATE)&&
		      (e2.first <= value)&&
		      (result12.second==result_ear2.second))
		    {
		      ear2_valid = test_merge(A, ear2_e, result_ear2, v2,
					      get_smallest_radius_delaunay_sphere(A, ear2_c, 
										  ear2.second));
		      // si on veut etre plus restrictif
		      // et exiger au moins une bonne pliure:
// 		      int test_merge_ear2 = 
// 			test_merge(ear2_e, result_ear2, v2,
// 				   ear2_c->get_smallest_radius(ear2.second));		      
// 		      ear2_valid = (test_merge_ear2&&(e2.first < -1)&&
// 				    ((value < -1)||(test_merge_ear2 == 1)));
		    } 

		  if ((!ear1_valid)&&(!ear2_valid)) 
		    return not_valid_connecting_case;

		  IO_edge_type* p1;  
		  IO_edge_type* p2;		  

		  border_extend(ordered_key, result12, 
				v1, v2, c->vertex(i),
				e1, e2, p1, p2);

		  if (ear1_valid&&ear2_valid&&(ear1_e==ear2_e))
		    {
		      if (e1.first < e2.first)
			{ 
			  Validation_case res = validate(ear1, A, e1.first);
			  if (!((res == ear_case)||(res == final_case)))
			    std::cerr << "+++probleme de recollement : cas " 
				      << res << std::endl;
			  e2 = compute_value(edge_Ifacet_2, A);

			  if (ordered_key.first == v1)
			    p2 = set_again_border_elt(c->vertex(i), v2,
						      Border_elt(e2, result2.second));
			  else
			    p2 = set_again_border_elt(v2, c->vertex(i),
						      Border_elt(e2, result2.second));

// 			  if (e2.first < NOT_VALID_CANDIDATE)
			    _ordered_border->insert(Radius_ptr_type(e2.first, p2));
// 			  else 
// 			    try_to_close_border(p2);
			}
		      else
			{
			  Validation_case res = validate(ear2, A, e2.first);
			  if (!((res == ear_case)||(res == final_case)))
			    std::cerr << "+++probleme de recollement : cas " 
				      << res << std::endl;
			  e1 = compute_value(edge_Ifacet_1, A);

			  if (ordered_key.first == v1) 
			    p1 = set_again_border_elt(v1, c->vertex(i),
						      Border_elt(e1, result1.second));
			  else 
			    p1 = set_again_border_elt(c->vertex(i), v1,
						      Border_elt(e1, result1.second));

// 			  if (e1.first < NOT_VALID_CANDIDATE)
			    _ordered_border->insert(Radius_ptr_type(e1.first, p1));
// 			  else 
// 			    try_to_close_border(p1);
			}
		    }
		  else// les deux oreilles ne se recollent pas sur la meme arete...
		    {
		      // on resoud la singularite.
		      if (ear1_valid)
			{
			  Validation_case res = validate(ear1, A, e1.first);
			  if (!((res == ear_case)||(res == final_case)))
			    std::cerr << "+++probleme de recollement : cas " 
				      << res << std::endl;
			}
		      if (ear2_valid)
			{
			  Validation_case res = validate(ear2, A, e2.first);
			  if (!((res == ear_case)||(res == final_case)))
			    std::cerr << "+++probleme de recollement : cas " 
				      << res << std::endl;
			}
		      // on met a jour la PQ s'il y a lieu... mais surtout pas
		      // avant la resolution de la singularite
		      if (!ear1_valid)
			{
// 			  if (e1.first < NOT_VALID_CANDIDATE)
			    _ordered_border->insert(Radius_ptr_type(e1.first, p1));
// 			  else
// 			    try_to_close_border(p1);
			}
		      if (!ear2_valid)		
			{
// 			  if (e2.first < NOT_VALID_CANDIDATE)
			    _ordered_border->insert(Radius_ptr_type(e2.first, p2));
// 			  else
// 			    try_to_close_border(p2);
			}
		    }

		  //Pour une visu correcte
		  visu_facet(c, edge_Efacet.second);

		  return connecting_case;
		}

// 	      if (is_on_same_border)
// 		{
// 		  _number_of_border++;
// 		  Incident_border_iterator tmp;
// 		  std::cout << "En train de separer deux bords :"
// 			    << result12.second << " et " << _number_of_border 
// 			    << std::endl;
// 		  Vertex_handle circ =  c->vertex(i), done = circ;
// 		  do
// 		    {
// 		      tmp = circ->get_next_on_border(result12.second);
// 		      (*tmp).second.second = _number_of_border;
// 		      circ = (Vertex*) (*tmp).first;
// 		    }
// 		  while(circ != done);
// 		}
// 	      else
// 		{ 
// 		  //certainement un probleme avec le recollement de bords differents
// 		  if (c->vertex(i)->not_interior() > 1)
// 		    {
// 		      std::cout << "En train de recoller des bords au bord :"
// 				<< result12.second << std::endl;
// 		      Incident_border_iterator tmp;
// 		      for(Incident_border_iterator it =
// 			    c->vertex(i)->first_incident();
// 			  it != c->vertex(i)->not_incident();
// 			  it++)
// 			{
// 			  int current_index = (*it).second.second;
// 			  if (current_index != result12.second)
// 			    {
// 			      Vertex_handle circ =  c->vertex(i), done = circ;
// 			      do
// 				{
// 				  tmp = circ->get_next_on_border(current_index);
// 				  (*tmp).second.second = result12.second;
// 				  circ = (Vertex*) (*tmp).first;
// 				}
// 			      while(circ != done);
// 			    }
// 			}
// 		    }
// 		}
	    }
	}							  	  
    }
  return not_valid;
}

//=====================================================================
//=====================================================================

void re_compute_values(const Triangulation_3& A)
{
  if(!_ordered_border->empty())
    {
      Ordered_border_type* _ordered_border_tmp = new Ordered_border_type();
      do 
	{
	  Ordered_border_iterator e_it = _ordered_border->begin();
	  Edge_IFacet mem_Ifacet =  e_it->second->first;
	  Cell_handle c_tmp = (Cell*) mem_Ifacet.first.first;
	  _ordered_border->erase(e_it);
	  Vertex_handle v1 = c_tmp->vertex(mem_Ifacet.first.second);
	  Vertex_handle v2 = c_tmp->vertex(mem_Ifacet.first.third);

	  Radius_edge_type new_candidate;
	  new_candidate = compute_value(mem_Ifacet, A);

// 	  if (new_candidate.first < NOT_VALID_CANDIDATE)
	    {	
	      if (new_candidate.first == STANDBY_CANDIDATE)
		{
		  // a garder pour un K un peu plus grand...
		  new_candidate.first = STANDBY_CANDIDATE_BIS;
		}

	      Border_elt result;
	      Edge_like key_tmp(v1,v2);
	      is_border_elt(key_tmp, result);
	      IO_edge_type* pnew = 
		set_again_border_elt(key_tmp.first, key_tmp.second, 
				     Border_elt (new_candidate, result.second));
	      _ordered_border_tmp->insert(Radius_ptr_type(new_candidate.first, pnew));
	    }
	}
      while(!_ordered_border->empty());
      delete _ordered_border;
      _ordered_border = _ordered_border_tmp;
    }
}

//---------------------------------------------------------------------

void 
#ifndef BLIND
loop(const Triangulation_3& A, 
     const coord_type& K_init, const coord_type& K_step, const coord_type& K_max, 
     CGAL::Geomview_stream& gv)
#endif //BLIND
#ifdef BLIND
loop(const Triangulation_3& A, 
     const coord_type& K_init, const coord_type& K_step, const coord_type& K_max)
#endif //BLIND
{
  // initilisation de la variable globale K: qualite d'echantillonnage requise
  K = K_init; // valeur d'initialisation de K pour commencer prudemment...
  //-------------------------------------------------------------------
// modif
  //int _facet_number_test = _last_component_facet_number+3;
  //bool close_result(false);
// modif
  int count(1);
  int nnn(-1);
  Vertex_handle v1, v2;
#ifndef BLIND
  std::cout << "   entrer un entier=nbre d'it (si =0 quit, <0 termine la loop): ";
  std::cin >> nnn;
#endif //BLIND
  t1.start();
  if (nnn == 0) std::exit(0);
  if (_ordered_border->empty()) return;
  do
    {
      min_K = HUGE_VAL; // pour retenir le prochain K necessaire pour progresser...
      do
	{ 
	  if (nnn > 0)
	    {
	      if (count%nnn == 0)
		{    
		  t1.stop();
		  //Pour une visu correcte
#ifndef BLIND
		  gv.clear();
		  show_selected_facets(gv, A); 

		  std::cout << "   surface partielle reconstuite en " << 
		    t1.time() << " secondes." << std::endl;
		  std::cout << "   " << _facet_number << 
		    " facettes, " << _vh_number << " sommets." << std::endl;
	  
		  std::cout << "   entrer un entier=nbre d'it (si <=0 quit): ";

		  std::cin >> nnn;
		  if (nnn <= 0) exit(0);
#endif //BLIND
		  t1.start();
		}
	      count++;
	    }

	  Ordered_border_iterator e_it = _ordered_border->begin();

	  criteria value = e_it->first;
	  if (value >= STANDBY_CANDIDATE)
	    re_compute_values(A);
	  else
	    {
	      Edge_IFacet candidate = e_it->second->second;
	      Cell_handle c_ext = (Cell*) candidate.first.first;
	      int i1, i2 , i3;
	      i1 = candidate.first.second;
	      i2 = candidate.first.third;
	      i3 = (6 - i1- i2 - candidate.second);

	      Edge_IFacet mem_Ifacet =  e_it->second->first;
	      Cell_handle c_tmp = (Cell*) mem_Ifacet.first.first;

	      v1 = c_tmp->vertex(mem_Ifacet.first.second);
	      v2 = c_tmp->vertex(mem_Ifacet.first.third);

	      Radius_edge_type mem_e_it(e_it->first, *e_it->second);
	      // Radius_ptr_type mem_first_it(*e_it);

	      _ordered_border->erase(e_it);
	      
// modif: Pour boucher les trous triangulaires avant de faire des conneries???
	      //if (_facet_number > _facet_number_test)
	      //close_result = try_to_close_border(e_it->second);

	      //if (!close_result)
	      //{
// fin de la modif...   
		  Validation_case validate_result = validate(candidate, A, value);
		  //      Cell_handle ccc = (Cell*) candidate.first.first;
		  if ((validate_result == not_valid)||
		      (validate_result == not_valid_connecting_case))
		    { 
		      Radius_edge_type new_candidate;
		      Border_elt result;
		      Edge_like key_tmp(v1,v2);
		      is_border_elt(key_tmp, result);

		      if (validate_result == not_valid_connecting_case)
			set_incidence_request(c_ext->vertex(i3), value, key_tmp);

		      if (validate_result == not_valid)
			{ 
			  new_candidate = compute_value(mem_Ifacet, A);
			  if ((new_candidate != mem_e_it))
// 			      &&(new_candidate.first < NOT_VALID_CANDIDATE))
			    {
			      IO_edge_type* pnew = 
				set_again_border_elt(key_tmp.first, key_tmp.second, 
						     Border_elt (new_candidate, result.second));
			      _ordered_border->insert(Radius_ptr_type(new_candidate.first,
								      pnew));
			    }
// 			  else
// 			    try_to_close_border(e_it->second);
			}
		    }
		  else // valid candidate...
		    {
		      //	  alpha_max = std::max(alpha_max, value);
		      //  if (validate_result != final_case)
		      // 	    {
		      // 	      Radius_edge_type v1_candidate, v2_candidate, v3_candidate;
		      // 	      if (v1->not_interior())
		      // 		if (v1->is_incidence_requested())//le bord ayant change autant essayer
		      // 		  v1_candidate = v1->get_best_incidence_request();
		      // 	      if (v2->not_interior())
		      // 		if (v2->is_incidence_requested())
		      // 		  v2_candidate = v2->get_best_incidence_request();
		      // 	      if (validate_result != exterior_case)
		      // 		if (c_ext->vertex(i3)->is_incidence_requested())
		      // 		  v3_candidate =
		      // 		    c_ext->vertex(i3)->get_best_incidence_request();
		      // 	    }
		    }
// modif
		  //}
// modif
	    }
	}
      while((!_ordered_border->empty())&&
	    (_ordered_border->begin()->first < STANDBY_CANDIDATE_BIS));

      K += std::max(K_step, min_K-K+eps); 
      // on augmente progressivement le K mais on a deja rempli sans
      // faire des betises auparavant...
    }
  while((!_ordered_border->empty())&&(K <= K_max)&&(min_K != HUGE_VAL));
  t1.stop();
  if ((min_K < HUGE_VAL)&&(!_ordered_border->empty()))
    std::cout << "   [ next K required = " << min_K << " ]" << std::endl;
}


//================== POST TRAITEMENT ==================================
//=====================================================================
#include <NUAGE/postprocessing.h>
//=====================================================================
//=====================================================================


//=====================================================================
//------------------ main ---------------------------------------------
//=====================================================================
//=====================================================================

int main(int argc,  char* argv[])
{
  //parse command line
  Options opt;
  std::cout << ">> option line for this execution is :" << std::endl;
  if (!parse(argc, argv, opt))
    std::exit(0);
  std::cout << std::endl << std::endl;
  DELTA = opt.DELTA;
  
  min_K = HUGE_VAL; // pour retenir le prochain K necessaire pour progresser...
  //------------------
  
  Triangulation_3 A;

#ifndef BLIND  
  CGAL::Geomview_stream gv_main(CGAL::Bbox_3(0,0,0, 2, 2, 2));
  gv_main.set_line_width(4);
  gv_main.set_trace(false);
  gv_main.set_bg_color(CGAL::Color(0, 200, 200));
#endif //BLIND

  std::vector<Point> L;
  bool re_init(false);
  int number_of_connected_comp(0);
  coord_type total_time(0);

  if (!opt.Section_file)
    file_input(opt.finname, opt.number_of_points, L);
  else
    section_file_input(opt.finname, opt.number_of_points, L);

  construct_delaunay(L, A);
  //assert(A.is_valid());
  
  L.clear();
  _A_size_before_postprocessing = A.number_of_vertices();

  do
    {
      number_of_connected_comp++;
      coord_type sum_time;
      if (re_init)
	std::cout << ">> searching another grain [init " <<
	  number_of_connected_comp << "] : "
		  << std::endl;
      t1.start();
      bool result_init = init(A, re_init);
      t1.stop();     
      if (result_init)
	{
// modif
	  //_last_component_facet_number = _facet_number;
// modif
	  sum_time = t1.time();
	  std::cout << std::endl << ">> initialisation finished after " << t1.time() << " seconds." <<
	    std::endl << std::endl;
	  t1.reset();
	  std::cout << ">> selection loop : " << std::endl;
#ifndef BLIND
	  loop(A, opt.K_init, opt.K_step, opt.K, gv_main);
#endif //BLIND
#ifdef BLIND
	  loop(A, opt.K_init, opt.K_step, opt.K);
#endif //BLIND
	  std::cout << "-- selection loop finished after  " << t1.time() 
		    << " seconds." << std::endl;

	  sum_time += t1.time();
	  t1.reset();

	  // debut du processus loop + postprocessing
	  if ((_facet_number > _A_size_before_postprocessing)&&
	      (opt.NB_BORDER_MAX > 0))
	    // en principe 2*nb_sommets = nb_facettes: y a encore de la marge!!!
	    {
	      bool post_trait;
	 
	      do
		{
		  t1.reset();
		  post_trait = postprocessing(A, opt.NB_BORDER_MAX);

		  std::cout << "-- postprocessing finished after " << t1.time() 
			    << " seconds." << std::endl;
		  sum_time += t1.time();
		  t1.reset();
		  if (post_trait)
		    {
		      std::cout << std::endl << ">> selection loop : " << std::endl;
#ifndef BLIND
		      loop(A, opt.K_init, opt.K_step, opt.K, gv_main);
#endif //BLIND
#ifdef BLIND
		      loop(A, opt.K_init, opt.K_step, opt.K);
#endif //BLIND
		      std::cout  << "-- selection loop finished after  " 
				 << t1.time() << " seconds." << std::endl;
	      
		      sum_time += t1.time();
		    }
		  else
		    std::cout << "   postprocessing finished." 
			      << std::endl << std::endl;
		}
	      while(post_trait);
	    }
#ifndef BLIND
	    gv_main.clear();
	    show_selected_facets(gv_main, A);
#endif //BLIND

	    //   gv_main << A;
	    std::cout << "-- reconstruction (init+loop+post) finished after  " <<
	      sum_time << " seconds." << std::endl;
	    total_time += sum_time;

#ifndef BLIND
	    std::cout << "   "  << _A_size_before_postprocessing - A.number_of_vertices()
		      << " point(s) ignored." << std::endl;
	    std::cout << "   Reconstructed surface: " << _facet_number << 
	      " facets, " << _vh_number << " vertices." << std::endl << std::endl;
#endif //BLIND

#ifdef BLIND
	    std::cout << std::endl;
#endif //BLIND

	  if ((A.number_of_vertices()-_vh_number) > 4)
	    {
	      int cont(1);
#ifndef BLIND
	      std::cout << ">> continue ? (si <= 0 quit): ";
	      std::cin >> cont;
#endif //BLIND
	      if (cont > 0)
		{
		  re_init = true;
		  // pour continuer il suffit de virer les sommets de _vh_vect...
		  t1.reset();
		}
	      else
		re_init = false;
	    }
	  else
	    re_init = false;
	}
      else
	{
	  std::cout << "-- no grains...."
		    << std::endl << std::endl;
	  re_init = false;
	}
    }while(re_init&&
	   ((number_of_connected_comp < opt.max_connected_comp)||
	    (opt.max_connected_comp < 0)));

  std::cout << "-- totalite de la selection terminee en  " <<
    total_time << " secondes." << std::endl;

 
  dump_in_file_selected_facets(opt.foutname, A, opt.contour, opt.out_format);


#ifdef BLIND
  std::cout << "   "  << _A_size_before_postprocessing - _vh_number
	    << " points ignored." << std::endl; 
  std::cout << "   Reconstructed surface: " << _facet_number << 
    " facets, " << _vh_number << " vertices." << std::endl;
  std::cout << "   "  << border_counter(A) << 
    " border edges." << std::endl;
  std::cout << "   number of connected components <= " 
	    << std::max(1, number_of_connected_comp-1)
	    << std::endl << std::endl;
#endif //BLIND

#ifndef BLIND
  int end_int(1);
  while(end_int)
    std::cin >> end_int;
#endif //BLIND
 
  return 0;
}

//=====================================================================
//=====================================================================
