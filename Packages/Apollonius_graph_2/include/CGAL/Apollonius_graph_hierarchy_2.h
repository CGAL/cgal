#ifndef CGAL_APOLLONIUS_GRAPH_HIERARCHY_2_H
#define CGAL_APOLLONIUS_GRAPH_HIERARCHY_2_H

#include <CGAL/Random.h>
#include <map>
#include <CGAL/Triangulation_hierarchy_2.h>

#include <CGAL/Apollonius_graph_2.h>
#include <CGAL/Apollonius_graph_data_structure_2.h>
#include <CGAL/Apollonius_graph_vertex_base_2.h>
#include <CGAL/Apollonius_graph_face_base_2.h>

CGAL_BEGIN_NAMESPACE

#if 0
template < class Vbb>
class Apollonius_graph_hierarchy_vertex_base_2
 : public Triangulation_hierarchy_vertex_base_2<Vbb>
{
 public:
  typedef Vbb V_Base;
  typedef typename V_Base::Point            Point;
  typedef typename V_Base::Weighted_point   Weighted_point;
};
#endif

// parameterization of the  hierarchy
const int ad_hierarchy_2__ratio    = 30;
const int ad_hierarchy_2__minsize  = 20;
const int ad_hierarchy_2__maxlevel = 5;
// maximal number of points is 30^5 = 24 millions !

template < class Gt, bool StoreTrivial = true,
  class Tds = Apollonius_graph_data_structure_2<
    Triangulation_hierarchy_vertex_base_2<
       Apollonius_graph_vertex_base_2<Gt,StoreTrivial> >,
    Apollonius_graph_face_base_2<Gt> > >
class Apollonius_graph_hierarchy_2
  : public Apollonius_graph_2< Gt, StoreTrivial, Tds >
{
  // public:
private:
  typedef Apollonius_graph_2<Gt,StoreTrivial,Tds> Apollonius_graph_base;
  typedef Apollonius_graph_base                   Ag_base;
  typedef typename Gt::Weighted_point             Weighted_point;

#if 0
  typedef Tds                          Data_structure;
  typedef Gt                           Geom_traits;
  typedef typename Gt::Bare_point      Point;
  typedef typename Gt::Site            Site;
  typedef typename Gt::Weight          Weight;

  typedef typename Ag_base::Face_handle      Face_handle;
  typedef typename Ag_base::Vertex_handle    Vertex_handle;
  typedef typename Ag_base::Edge             Edge;

  typedef typename Ag_base::Face_circulator       Face_circulator;
  typedef typename Ag_base::Edge_circulator       Edge_circulator;
  typedef typename Ag_base::Vertex_circulator     Vertex_circulator;

  typedef typename Ag_base::All_faces_iterator    All_faces_iterator;
  typedef typename Ag_base::Finite_faces_iterator Finite_faces_iterator;

  typedef typename Ag_base::All_vertices_iterator All_vertices_iterator;
  typedef typename Ag_base::Finite_vertices_iterator 
                                                     Finite_vertices_iterator;

  typedef typename Ag_base::All_edges_iterator    All_edges_iterator;
  typedef typename Ag_base::Finite_edges_iterator Finite_edges_iterator;

  typedef Finite_vertices_iterator   Vertex_iterator;
  typedef Finite_faces_iterator      Face_iterator;
  typedef Finite_edges_iterator      Edge_iterator;

#endif

 private:
  // here is the stack of triangulations which form the hierarchy
  Apollonius_graph_base*   hierarchy[ad_hierarchy_2__maxlevel];
  Random random; // random generator

public:
  Apollonius_graph_hierarchy_2
  (const Geom_traits& gt = Geom_traits());

  template<class Input_iterator>
  Apollonius_graph_hierarchy_2(Input_iterator first,
			       Input_iterator beyond,
			       const Geom_traits& gt = Geom_traits())
    : Apollonius_graph_hierarchy_2(gt)
  {
    for (Input_iterator it = first; it != beyond; ++it) {
      insert(*it);
    }
  }

  Apollonius_graph_hierarchy_2
  (const Apollonius_graph_hierarchy_2& ad);

  Apollonius_graph_hierarchy_2&
  operator=(const Apollonius_graph_hierarchy_2& ad);
  ~Apollonius_graph_hierarchy_2();

  //Helping
  void copy(const Apollonius_graph_hierarchy_2 &ad);

  void clear();

  // CHECKING
  bool is_valid(bool verbose = false, int level = 1) const;


  // INSERT REMOVE
  Vertex_handle insert(const Weighted_point& p);
 
  template < class Input_iterator >
  void insert(Input_iterator first, Input_iterator beyond,
	      bool sort_first = false)
  {
    //    int i(0);
    for (Input_iterator it = first; it != beyond; ++it) {
      insert(*it);
      //      i++;
    }
    //    return i;
#if 0
    Input_iterator it;

    typename Apollonius_graph_base::Weighted_point_less_than_comparator
      less_than(geom_traits());
    typename Apollonius_graph_base::Weighted_point_list wp_list;
    for (it = first; it != beyond; ++it)  wp_list.push_back(*it);
    std::sort(wp_list.begin(), wp_list.end(), less_than);

    typename Apollonius_graph_base::Weighted_point_list_iterator lit;
    for (lit = wp_list.begin(); lit != wp_list.end(); ++lit) {
      insert(*lit);
    }
    wp_list.clear();

    int n = number_of_vertices();
#if 0
    while(first != last){
      insert(*first);
      ++first;
    }
#endif
    return number_of_vertices() - n;
#endif
  }

  inline void remove(Vertex_handle v) {
    remove(v, true);
  }

  void remove(Vertex_handle v, bool insert_trivial);

public:
  // find nearest neighbor
  Vertex_handle  nearest_neighbor(const Point& p) const;

private:
  void  nearest_neighbor(const Point& p,
			 Vertex_handle vnear[ad_hierarchy_2__maxlevel])
    const; 
  int random_level();
};


template < class Gt, bool StoreTrivial, class Tds>
Apollonius_graph_hierarchy_2<Gt,StoreTrivial,Tds>::
Apollonius_graph_hierarchy_2(const Geom_traits& gt)
  : Apollonius_graph_base(gt), random((long)0)
{ 
  hierarchy[0] = this; 
  for(int i = 1; i < ad_hierarchy_2__maxlevel; ++i) {
    hierarchy[i] = new Apollonius_graph_base(gt);
  }
}


// copy constructor duplicates vertices and faces
template <class Gt, bool StoreTrivial, class Tds>
Apollonius_graph_hierarchy_2<Gt,StoreTrivial,Tds>::
Apollonius_graph_hierarchy_2
(const Apollonius_graph_hierarchy_2<Gt,StoreTrivial,Tds> &ad)
    : Apollonius_graph_base(), random((long)0)
{ 
  // create an empty triangulation to be able to delete it !
  hierarchy[0] = this; 
  for(int i = 1; i < ad_hierarchy_2__maxlevel; ++i)
    hierarchy[i] = new Apollonius_graph_base(ad.geom_traits());
  copy(ad);
} 
 

//Assignement
template <class Gt, bool StoreTrivial, class Tds>
Apollonius_graph_hierarchy_2<Gt,StoreTrivial,Tds> &
Apollonius_graph_hierarchy_2<Gt,StoreTrivial,Tds>::
operator=(const Apollonius_graph_hierarchy_2<Gt,StoreTrivial,Tds> &ad)
{
  copy(ad);
  return *this;
}

template <class Gt, bool StoreTrivial, class Tds>
void
Apollonius_graph_hierarchy_2<Gt,StoreTrivial,Tds>::   
copy
(const Apollonius_graph_hierarchy_2<Gt,StoreTrivial,Tds> &ad)
{
  std::map< const void*, void*, std::less<const void*> > V;
  {
    for(int i = 0; i < ad_hierarchy_2__maxlevel; ++i) {
      //      hierarchy[i]->copy_triangulation(*awvd.hierarchy[i]);
      *(hierarchy[i]) = *ad.hierarchy[i];
    }
  }
  //up and down have been copied in straightforward way
  // compute a map at lower level
  {
    for( Vertex_iterator it = hierarchy[0]->vertices_begin(); 
	 it != hierarchy[0]->vertices_end(); ++it) {
      if (it->up()) V[ ((Vertex*)(it->up()))->down() ] = &(*it);
    }
  }
  {
    for(int i = 1; i < ad_hierarchy_2__maxlevel; ++i) {
      for( Vertex_iterator it = hierarchy[i]->vertices_begin(); 
	   it != hierarchy[i]->vertices_end(); ++it) {
	// down pointer goes in original instead in copied triangulation
	it->set_down(V[it->down()]);
	// make reverse link
	((Vertex*)(it->down()))->set_up( &(*it) );
	// make map for next level
	if (it->up()) V[ ((Vertex*)(it->up()))->down() ] = &(*it);
      }
    }
  }
}

template <class Gt, bool StoreTrivial, class Tds>
Apollonius_graph_hierarchy_2<Gt,StoreTrivial,Tds>:: 
~Apollonius_graph_hierarchy_2()
{
  clear();
  for(int i = 1; i < ad_hierarchy_2__maxlevel; ++i) {
    delete hierarchy[i];
  }
}

template <class Gt, bool StoreTrivial, class Tds>
void
Apollonius_graph_hierarchy_2<Gt,StoreTrivial,Tds>:: 
clear()
{
  for(int i = 0; i < ad_hierarchy_2__maxlevel; ++i) {
    hierarchy[i]->clear();
  }
}

template<class Gt, bool StoreTrivial, class Tds>
bool
Apollonius_graph_hierarchy_2<Gt,StoreTrivial,Tds>:: 
is_valid(bool verbose, int level) const
{
  bool result(true);

  //verify correctness of triangulation at all levels
  for(int i = 0; i < ad_hierarchy_2__maxlevel; ++i) {
    if ( verbose ) {
      std::cout << "Level " << i << ": " << std::flush;
    }
    result = result && hierarchy[i]->is_valid(verbose, level);
    if ( verbose ) {
      std::cout << std::endl;
    }
  }
  //verify that lower level has no down pointers
  for( Vertex_iterator it = hierarchy[0]->vertices_begin(); 
       it != hierarchy[0]->vertices_end(); ++it) {
    result = result && ( it->down() == 0 );
  }

  //verify that other levels has down pointer and reciprocal link is fine
  for(int i = 1; i < ad_hierarchy_2__maxlevel; ++i) {
    for( Vertex_iterator it = hierarchy[i]->vertices_begin(); 
	 it != hierarchy[i]->vertices_end(); ++it) {
      result = result && 
	( ((Vertex*)((Vertex*)it->down())->up()) ==  &(*it) );
    }
  }
  return result;
}


template <class Gt, bool StoreTrivial, class Tds>
Apollonius_graph_hierarchy_2<Gt,StoreTrivial,Tds>::Vertex_handle
Apollonius_graph_hierarchy_2<Gt,StoreTrivial,Tds>::
insert(const Weighted_point &p)
{
  int vertex_level = random_level();

  int n = number_of_vertices();
  Vertex_handle vertex;
  Vertex_handle vnear[ad_hierarchy_2__maxlevel];
  //  typename Apollonius_graph_base::List l(*this);
  typename Apollonius_graph_base::List l;
  typename Apollonius_graph_base::Face_map fm;
  typename Apollonius_graph_base::Vertex_map v_trivial;

  if ( n <= 2 ) {
    if ( n == 0 ) {
      vertex = insert_first(p);
    } else if ( n == 1 ) {
      vertex = insert_second(p);
    } else if ( n == 2 ) {
      vertex = insert_third(p);
    }

    Vertex_handle previous = vertex;
    Vertex_handle first = vertex;

    int level = 1;
    while (level <= vertex_level ){
      vertex = hierarchy[level]->insert(p, vnear[level]);
      vertex->set_down((void *) &*previous); // link with level above
      previous->set_up((void *) &*vertex);
      previous = vertex;
      level++;
    }
    return first;
  }

  int n_trivial = 0;

  // locate the nearest neighbor using hierarchy
  nearest_neighbor(p, vnear);

  CGAL_assertion( vnear[0].ptr() != NULL );

  // check if it is trivial
  Weighted_point wp_nearest = vnear[0]->point();
  if ( is_trivial_test(wp_nearest, p) ) {
    vnear[0]->add_weighted_point(p);
    return Vertex_handle(NULL);
  }

  // find the first conflict
  typename Apollonius_graph_base::Face_circulator fc_start =
    vnear[0]->incident_faces();
  typename Apollonius_graph_base::Face_circulator fc = fc_start;

  Face_handle start_f;
  Sign s;
  do {
    Face_handle f(fc);
    s = incircle_test(f, p);

    if ( s == NEGATIVE ) {
      start_f = f;
      break;
    }
    ++fc;
  } while ( fc != fc_start );

  if ( s != NEGATIVE ) {
    typename Apollonius_graph_base::Edge_circulator ec_start = 
      vnear[0]->incident_edges();
    typename Apollonius_graph_base::Edge_circulator ec = ec_start;

    bool interior_in_conflict(false);
    typename Apollonius_graph_base::Edge e;
    do {
      e = *ec;
      interior_in_conflict = edge_interior_test(e, p, false);
      
      if ( interior_in_conflict ) { break; }
      ++ec;
    } while ( ec != ec_start );

    CGAL_assertion( interior_in_conflict );

    vertex = insert_degree_2(e, p);

    // insert at other levels
    Vertex_handle previous = vertex;
    Vertex_handle first = vertex;

    int level = 1;
    while (level <= vertex_level ){
      vertex = hierarchy[level]->insert(p, vnear[level]);
      vertex->set_down((void *) &*previous); // link with level above
      previous->set_up((void *) &*vertex);
      previous = vertex;
      level++;
    }
    return first;
  }

  vertex = create_vertex();
  vertex->set_point(p);

  initialize_conflict_region(start_f, l);
  expand_conflict_region(start_f, p, l, fm, v_trivial, NULL);
  n_trivial = v_trivial.size();

  if ( n_trivial != 0 ) {
    typename Apollonius_graph_base::Vertex_map::iterator it;
    for (it = v_trivial.begin(); it != v_trivial.end(); it++) {
      Vertex_handle vh = (*it).first;
      Vertex* v = static_cast<Vertex*>(vh.ptr());
      void * u = v->up();
      if ( u != NULL ) {
	v = (Vertex*)u;
	u = v->up();
	int l = 1;
	while ( true ) {
	  hierarchy[l++]->remove(v, false);
	  if (!u) break; 
	  if(l >= ad_hierarchy_2__maxlevel) { break; }
	  v = (Vertex*)u;
	  u = v->up();
	}
      }
    }
  }

  // now really insert at level 0
  retriangulate_conflict_region(vertex, l, fm, v_trivial);

  fm.clear();
  v_trivial.clear();
  // end of insertion at level 0

  // insert at other levels
  Vertex_handle previous = vertex;
  Vertex_handle first = vertex;

  if ( n_trivial != 0 ) {
    nearest_neighbor(p.point(), vnear);
  }
      
  int level = 1;
  while (level <= vertex_level ){
    vertex = hierarchy[level]->insert(p, vnear[level]);
    vertex->set_down((void *) &*previous); // link with level above
    previous->set_up((void *) &*vertex);
    previous = vertex;
    level++;
  }
  return first;
}

template <class Gt, bool StoreTrivial, class Tds>
void 
Apollonius_graph_hierarchy_2<Gt,StoreTrivial,Tds>::
remove(Vertex_handle v, bool insert_trivial)
{
  CGAL_triangulation_precondition( v != Vertex_handle());
  CGAL_triangulation_precondition( !is_infinite(v));

  // get the trivial circles
  Weighted_point_list wp_list;
  typename Vertex_base::Weighted_point_list_iterator wpit;

  if ( insert_trivial ) {
    for (wpit = v->weighted_points_begin();
	 wpit != v->weighted_points_end(); ++wpit) {
      wp_list.push_back(*wpit);
    }
  }

  // do the actual removal
  void * u = v->up();
  int l = 0;
  while ( true ) {
    hierarchy[l++]->remove(v, false);
    if (!u) break; 
    if(l >= ad_hierarchy_2__maxlevel) break;
    v = (Vertex*)u;
    u = v->up();
  }

  for (unsigned int i = 0; i < wp_list.size(); i++) {
    insert(wp_list[i]);
  }
}


template <class Gt, bool StoreTrivial, class Tds>
Apollonius_graph_hierarchy_2<Gt,StoreTrivial,Tds>::Vertex_handle 
Apollonius_graph_hierarchy_2<Gt,StoreTrivial,Tds>::
nearest_neighbor(const Point& p) const
{
  Vertex_handle vnear[ad_hierarchy_2__maxlevel];
  nearest_neighbor(p, vnear);
  return vnear[0];
}

template <class Gt, bool StoreTrivial, class Tds>
void
Apollonius_graph_hierarchy_2<Gt,StoreTrivial,Tds>::
nearest_neighbor(const Point& p,
		 Vertex_handle vnear[ad_hierarchy_2__maxlevel])
  const
{
  Vertex_handle nearest = 0;
  int level  = ad_hierarchy_2__maxlevel;

  // find the highest level with enough vertices
  while ( hierarchy[--level]->number_of_vertices() 
	  < ad_hierarchy_2__minsize ) {
    if ( !level ) break;  // do not go below 0
  }
  for (int i = level+1; i < ad_hierarchy_2__maxlevel; ++i) {
    vnear[i] = 0;
  }

  while ( level > 0 ) {
    vnear[level] = nearest =
      hierarchy[level]->nearest_neighbor(p, nearest);  

    CGAL_assertion( !hierarchy[level]->is_infinite(vnear[level]) );
    // go at the same vertex on level below
    nearest = (Vertex*)( nearest->down() );
    --level;
  }
  vnear[0] =
    hierarchy[level]->nearest_neighbor(p, nearest);  // at level 0
}

template <class Gt, bool StoreTrivial, class Tds>
int
Apollonius_graph_hierarchy_2<Gt,StoreTrivial,Tds>::
random_level()
{
  int l = 0;
  while (1) {
    if ( random(ad_hierarchy_2__ratio) ) break;
    ++l;
  }
  if (l >= ad_hierarchy_2__maxlevel)
    l = ad_hierarchy_2__maxlevel -1;
  return l;
}

CGAL_END_NAMESPACE

#endif // CGAL_APOLLONIUS_GRAPH_HIERARCHY_2_H
