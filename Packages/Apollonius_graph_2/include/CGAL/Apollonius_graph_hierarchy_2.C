// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : include/CGAL/Apollonius_graph_hierarchy_2.C
// package       : Apollonius_graph_2
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ======================================================================



// class implementation
//---------------------

CGAL_BEGIN_NAMESPACE

template < class Gt, bool StoreHidden, class Agds>
Apollonius_graph_hierarchy_2<Gt,StoreHidden,Agds>::
Apollonius_graph_hierarchy_2(const Geom_traits& gt)
  : Apollonius_graph_base(gt), random((long)0)
{ 
  hierarchy[0] = this; 
  for(int i = 1; i < ag_hierarchy_2__maxlevel; ++i) {
    hierarchy[i] = new Apollonius_graph_base(gt);
  }
}


// copy constructor duplicates vertices and faces
template <class Gt, bool StoreHidden, class Agds>
Apollonius_graph_hierarchy_2<Gt,StoreHidden,Agds>::
Apollonius_graph_hierarchy_2
(const Apollonius_graph_hierarchy_2<Gt,StoreHidden,Agds>& agh)
    : Apollonius_graph_base(), random((long)0)
{ 
  // create an empty triangulation to be able to delete it !
  hierarchy[0] = this; 
  for(int i = 1; i < ag_hierarchy_2__maxlevel; ++i)
    hierarchy[i] = new Apollonius_graph_base(agh.geom_traits());
  copy(agh);
} 
 

//Assignement
template <class Gt, bool StoreHidden, class Agds>
Apollonius_graph_hierarchy_2<Gt,StoreHidden,Agds> &
Apollonius_graph_hierarchy_2<Gt,StoreHidden,Agds>::
operator=(const Apollonius_graph_hierarchy_2<Gt,StoreHidden,Agds> &agh)
{
  copy(agh);
  return *this;
}

template <class Gt, bool StoreHidden, class Agds>
void
Apollonius_graph_hierarchy_2<Gt,StoreHidden,Agds>::
copy
(const Apollonius_graph_hierarchy_2<Gt,StoreHidden,Agds> &agh)
{
  std::map< const void*, void*, std::less<const void*> > V;
  for(int i = 0; i < ag_hierarchy_2__maxlevel; ++i) {
    //      hierarchy[i]->copy_triangulation(*awvd.hierarchy[i]);
    *(hierarchy[i]) = *agh.hierarchy[i];
  }

  //up and down have been copied in straightforward way
  // compute a map at lower level
  for( Finite_vertices_iterator it = hierarchy[0]->finite_vertices_begin(); 
       it != hierarchy[0]->finite_vertices_end(); ++it) {
    if (it->up()) V[ ((Vertex*)(it->up()))->down() ] = &(*it);
  }

  for(int i = 1; i < ag_hierarchy_2__maxlevel; ++i) {
    for( Finite_vertices_iterator it = hierarchy[i]->finite_vertices_begin(); 
	 it != hierarchy[i]->finite_vertices_end(); ++it) {
      // down pointer goes in original instead in copied triangulation
      it->set_down(V[it->down()]);
      // make reverse link
      ((Vertex*)(it->down()))->set_up( &(*it) );
      // make map for next level
      if (it->up()) V[ ((Vertex*)(it->up()))->down() ] = &(*it);
    }
  }
}

template <class Gt, bool StoreHidden, class Agds>
Apollonius_graph_hierarchy_2<Gt,StoreHidden,Agds>:: 
~Apollonius_graph_hierarchy_2()
{
  clear();
  for(int i = 1; i < ag_hierarchy_2__maxlevel; ++i) {
    delete hierarchy[i];
  }
}

template <class Gt, bool StoreHidden, class Agds>
void
Apollonius_graph_hierarchy_2<Gt,StoreHidden,Agds>:: 
clear()
{
  for(int i = 0; i < ag_hierarchy_2__maxlevel; ++i) {
    hierarchy[i]->clear();
  }
}

template<class Gt, bool StoreHidden, class Agds>
bool
Apollonius_graph_hierarchy_2<Gt,StoreHidden,Agds>:: 
is_valid(bool verbose, int level) const
{
  bool result(true);

  //verify correctness of triangulation at all levels
  for(int i = 0; i < ag_hierarchy_2__maxlevel; ++i) {
    if ( verbose ) {
      std::cout << "Level " << i << ": " << std::flush;
    }
    result = result && hierarchy[i]->is_valid(verbose, level);
    if ( verbose ) {
      std::cout << std::endl;
    }
  }
  //verify that lower level has no down pointers
  for( Finite_vertices_iterator it = hierarchy[0]->finite_vertices_begin(); 
       it != hierarchy[0]->finite_vertices_end(); ++it) {
    result = result && ( it->down() == 0 );
  }

  //verify that other levels has down pointer and reciprocal link is fine
  for(int i = 1; i < ag_hierarchy_2__maxlevel; ++i) {
    for( Finite_vertices_iterator it = hierarchy[i]->finite_vertices_begin(); 
	 it != hierarchy[i]->finite_vertices_end(); ++it) {
      result = result && 
	( ((Vertex*)((Vertex*)it->down())->up()) ==  &(*it) );
    }
  }
  return result;
}


template <class Gt, bool StoreHidden, class Agds>
typename Apollonius_graph_hierarchy_2<Gt,StoreHidden,Agds>::Vertex_handle
Apollonius_graph_hierarchy_2<Gt,StoreHidden,Agds>::
insert(const Weighted_point &p)
{
  int vertex_level = random_level();

  int n = number_of_vertices();
  Vertex_handle vertex;
  Vertex_handle vnear[ag_hierarchy_2__maxlevel];

  typename Apollonius_graph_base::List l;
  typename Apollonius_graph_base::Face_map fm;
  typename Apollonius_graph_base::Vertex_map v_hidden;

  if ( n <= 2 ) {
    if ( n == 0 ) {
      vertex = insert_first(p);
    } else if ( n == 1 ) {
      vertex = insert_second(p);
    } else if ( n == 2 ) {
      vertex = insert_third(p);
    }

    // if it hidden just return it right away
    if ( &(*vertex) == NULL ) {
      return vertex;
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

  int n_hidden = 0;

  // locate the nearest neighbor using hierarchy
  nearest_neighbor(p, vnear);

  CGAL_assertion( vnear[0] != NULL );

  // check if it is hidden
  Weighted_point wp_nearest = vnear[0]->point();
  if ( is_hidden(wp_nearest, p) ) {
    vnear[0]->add_hidden_weighted_point(p);
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
    s = incircle(f, p);

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
      interior_in_conflict = edge_interior(e, p, false);
      
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

  initialize_conflict_region(start_f, l);
  expand_conflict_region(start_f, p, l, fm, v_hidden, NULL);
  n_hidden = v_hidden.size();

  if ( n_hidden != 0 ) {
    int n_non_hidden = number_of_vertices() - n_hidden;
    if ( n_non_hidden < 2 ) {
      for(int i = 1; i < ag_hierarchy_2__maxlevel; ++i) {
	hierarchy[i]->clear();
      }

      if ( n_non_hidden == 1 ) {
	Vertex_handle non_hidden;
	Finite_vertices_iterator vit = finite_vertices_begin();
	do {
	  non_hidden = Vertex_handle(vit);
	  ++vit;
	} while ( v_hidden.find(non_hidden) != v_hidden.end() );

	non_hidden->set_up(NULL);
      }
    } else {
      typename Apollonius_graph_base::Vertex_map::iterator it;
      for (it = v_hidden.begin(); it != v_hidden.end(); it++) {
	Vertex_handle vh = (*it).first;
	Vertex* v = static_cast<Vertex*>( &(*vh) );
	void * u = v->up();
	if ( u != NULL ) {
	  v = (Vertex*)u;
	  u = v->up();
	  int l = 1;
	  while ( true ) {
	    hierarchy[l++]->remove(v);
	    if (!u) break; 
	    if(l >= ag_hierarchy_2__maxlevel) { break; }
	    v = (Vertex*)u;
	    u = v->up();
	  }
	}
      }
    }
  }

  // now really insert at level 0
  vertex = retriangulate_conflict_region(p, l, fm, v_hidden);

  fm.clear();
  v_hidden.clear();
  // end of insertion at level 0

  // insert at other levels
  Vertex_handle previous = vertex;
  Vertex_handle first = vertex;

  if ( n_hidden != 0 ) {
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

template <class Gt, bool StoreHidden, class Agds>
void 
Apollonius_graph_hierarchy_2<Gt,StoreHidden,Agds>::
remove(Vertex_handle v)
{
  CGAL_triangulation_precondition( v != Vertex_handle());
  CGAL_triangulation_precondition( !is_infinite(v));

  // get the hidden circles
  Weighted_point_list wp_list;
  typename Vertex_base::Hidden_weighted_point_iterator wpit;

  for (wpit = v->hidden_weighted_points_begin();
       wpit != v->hidden_weighted_points_end(); ++wpit) {
    wp_list.push_back(*wpit);
  }
  v->clear_hidden_weighted_point_container();

  // do the actual removal
  void * u = v->up();
  int l = 0;
  while ( true ) {
    hierarchy[l++]->remove(v);
    if (!u) break; 
    if(l >= ag_hierarchy_2__maxlevel) break;
    v = (Vertex*)u;
    u = v->up();
  }

  insert(wp_list.begin(), wp_list.end());
}


template <class Gt, bool StoreHidden, class Agds>
typename Apollonius_graph_hierarchy_2<Gt,StoreHidden,Agds>::Vertex_handle 
Apollonius_graph_hierarchy_2<Gt,StoreHidden,Agds>::
nearest_neighbor(const Point& p) const
{
  Vertex_handle vnear[ag_hierarchy_2__maxlevel];
  nearest_neighbor(p, vnear);
  return vnear[0];
}



template <class Gt, bool StoreHidden, class Agds>
void
Apollonius_graph_hierarchy_2<Gt,StoreHidden,Agds>::
swap(Apollonius_graph_hierarchy_2<Gt,StoreHidden,Agds>& agh)
{
  Ag_base* temp;
  Ag_base::swap(tr);
  for(int i = 1; i < ag_hierarchy_2__maxlevel; ++i){
    temp = hierarchy[i];
    hierarchy[i] = agh.hierarchy[i];
    agh.hierarchy[i]= temp;
  }
}




template <class Gt, bool StoreHidden, class Agds>
void
Apollonius_graph_hierarchy_2<Gt,StoreHidden,Agds>::
nearest_neighbor(const Point& p,
		 Vertex_handle vnear[ag_hierarchy_2__maxlevel])
  const
{
  Vertex_handle nearest = 0;
  int level  = ag_hierarchy_2__maxlevel;

  // find the highest level with enough vertices
  while ( hierarchy[--level]->number_of_vertices() 
	  < ag_hierarchy_2__minsize ) {
    if ( !level ) break;  // do not go below 0
  }
  for (int i = level+1; i < ag_hierarchy_2__maxlevel; ++i) {
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

template <class Gt, bool StoreHidden, class Agds>
int
Apollonius_graph_hierarchy_2<Gt,StoreHidden,Agds>::
random_level()
{
  int l = 0;
  while (1) {
    if ( random(ag_hierarchy_2__ratio) ) break;
    ++l;
  }
  if (l >= ag_hierarchy_2__maxlevel)
    l = ag_hierarchy_2__maxlevel -1;
  return l;
}

CGAL_END_NAMESPACE
