// Copyright (c) 2003  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>

#ifndef CGAL_APOLLONIUS_GRAPH_HIERARCHY_2_C
#define CGAL_APOLLONIUS_GRAPH_HIERARCHY_2_C


// class implementation
//---------------------

CGAL_BEGIN_NAMESPACE

template<class Gt, class Agds, class LTag>
void
Apollonius_graph_hierarchy_2<Gt,Agds,LTag>::
init_hierarchy(const Geom_traits& gt)
{
  hierarchy[0] = this; 
  for(unsigned int i = 1; i < ag_hierarchy_2__maxlevel; ++i) {
    hierarchy[i] = new Apollonius_graph_base(gt);
  }
}

template<class Gt, class Agds, class LTag>
Apollonius_graph_hierarchy_2<Gt,Agds,LTag>::
Apollonius_graph_hierarchy_2(const Geom_traits& gt)
  : Apollonius_graph_base(gt), random((long)0)
{ 
  init_hierarchy(gt);
}


// copy constructor duplicates vertices and faces
template<class Gt, class Agds, class LTag>
Apollonius_graph_hierarchy_2<Gt,Agds,LTag>::
Apollonius_graph_hierarchy_2
(const Apollonius_graph_hierarchy_2<Gt,Agds,LTag>& agh)
    : Apollonius_graph_base(agh.geom_traits()), random((long)0)
{ 
  init_hierarchy(agh.geom_traits());
  copy(agh);
} 
 

//Assignement
template<class Gt, class Agds, class LTag>
Apollonius_graph_hierarchy_2<Gt,Agds,LTag> &
Apollonius_graph_hierarchy_2<Gt,Agds,LTag>::
operator=(const Apollonius_graph_hierarchy_2<Gt,Agds,LTag> &agh)
{
  copy(agh);
  return *this;
}

template<class Gt, class Agds, class LTag>
void
Apollonius_graph_hierarchy_2<Gt,Agds,LTag>::
copy
(const Apollonius_graph_hierarchy_2<Gt,Agds,LTag> &agh)
{
  std::map< Vertex_handle, Vertex_handle > V;
  for(unsigned int i = 0; i < ag_hierarchy_2__maxlevel; ++i) {
    //      hierarchy[i]->copy_triangulation(*awvd.hierarchy[i]);
    *(hierarchy[i]) = *agh.hierarchy[i];
  }

  //up and down have been copied in straightforward way
  // compute a map at lower level
  for( Finite_vertices_iterator it = hierarchy[0]->finite_vertices_begin(); 
       it != hierarchy[0]->finite_vertices_end(); ++it) {
    if ( it->up() != NULL ) V[ it->up()->down() ] = it;
  }

  for(unsigned int i = 1; i < ag_hierarchy_2__maxlevel; ++i) {
    for( Finite_vertices_iterator it = hierarchy[i]->finite_vertices_begin(); 
	 it != hierarchy[i]->finite_vertices_end(); ++it) {
      // down pointer goes in original instead in copied triangulation
      it->set_down(V[it->down()]);
      // make reverse link
      it->down()->set_up( it );
      // make map for next level
      if ( it->up() != NULL ) V[ it->up()->down() ] = it;
    }
  }
}

template<class Gt, class Agds, class LTag>
Apollonius_graph_hierarchy_2<Gt,Agds,LTag>:: 
~Apollonius_graph_hierarchy_2()
{
  clear();
  for(unsigned int i = 1; i < ag_hierarchy_2__maxlevel; ++i) {
    delete hierarchy[i];
  }
}

template<class Gt, class Agds, class LTag>
void
Apollonius_graph_hierarchy_2<Gt,Agds,LTag>:: 
clear()
{
  for(unsigned int i = 0; i < ag_hierarchy_2__maxlevel; ++i) {
    hierarchy[i]->clear();
  }
}

template<class Gt, class Agds, class LTag>
bool
Apollonius_graph_hierarchy_2<Gt,Agds,LTag>:: 
is_valid(bool verbose, int level) const
{
  bool result(true);

  //verify correctness of triangulation at all levels
  for(unsigned int i = 0; i < ag_hierarchy_2__maxlevel; ++i) {
    if ( verbose ) {
      std::cerr << "Level " << i << ": " << std::flush;
    }
    result = result && hierarchy[i]->is_valid(verbose, level);
    if ( verbose ) {
      std::cerr << std::endl;
    }
  }
  //verify that lower level has no down pointers
  for( Finite_vertices_iterator it = hierarchy[0]->finite_vertices_begin(); 
       it != hierarchy[0]->finite_vertices_end(); ++it) {
    result = result && ( it->down() == 0 );
  }

  //verify that other levels has down pointer and reciprocal link is fine
  for(unsigned int i = 1; i < ag_hierarchy_2__maxlevel; ++i) {
    for( Finite_vertices_iterator it = hierarchy[i]->finite_vertices_begin(); 
	 it != hierarchy[i]->finite_vertices_end(); ++it) {
#ifdef CGAL_T2_USE_ITERATOR_AS_HANDLE 
      result = result && ( &*it == &*(it->down()->up()) );
#else
      result = result && ( it->down()->up() ==  it );
#endif
    }
  }
  return result;
}


template<class Gt, class Agds, class LTag>
typename Apollonius_graph_hierarchy_2<Gt,Agds,LTag>::Vertex_handle
Apollonius_graph_hierarchy_2<Gt,Agds,LTag>::
insert(const Site_2 &p)
{
  int vertex_level = random_level();

  size_type n = this->number_of_vertices();
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
    if ( vertex == NULL ) {
      return vertex;
    }

    Vertex_handle previous = vertex;
    Vertex_handle first = vertex;

    int level = 1;
    while (level <= vertex_level ){
      vertex = hierarchy[level]->insert(p, vnear[level]);
      vertex->set_down(previous); // link with level above
      previous->set_up(vertex);
      previous = vertex;
      level++;
    }
    return first;
  }

  int n_hidden = 0;

  // locate the nearest neighbor using hierarchy
  nearest_neighbor(p.point(), vnear);

  CGAL_assertion( vnear[0] != NULL );

  // check if it is hidden
  Site_2 wp_nearest = vnear[0]->site();
  if ( is_hidden(wp_nearest, p) ) {
    vnear[0]->add_hidden_site(p);
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
      vertex->set_down(previous); // link with level above
      previous->set_up(vertex);
      previous = vertex;
      level++;
    }
    return first;
  }

  initialize_conflict_region(start_f, l);
  expand_conflict_region(start_f, p, l, fm, v_hidden, NULL);
  n_hidden = v_hidden.size();

  if ( n_hidden != 0 ) {
    int n_non_hidden = this->number_of_vertices() - n_hidden;
    if ( n_non_hidden < 2 ) {
      for(unsigned int i = 1; i < ag_hierarchy_2__maxlevel; ++i) {
	hierarchy[i]->clear();
      }

      if ( n_non_hidden == 1 ) {
	Vertex_handle non_hidden;
	Finite_vertices_iterator vit = this->finite_vertices_begin();
	do {
	  non_hidden = Vertex_handle(vit);
	  ++vit;
	} while ( v_hidden.find(non_hidden) != v_hidden.end() );

	non_hidden->set_up(NULL);
      }
    } else {
      typename Apollonius_graph_base::Vertex_map::iterator it;
      for (it = v_hidden.begin(); it != v_hidden.end(); it++) {
	Vertex_handle v = (*it).first;
	Vertex_handle u = v->up();
	if ( u != NULL ) {
	  v = u;
	  u = v->up();
	  unsigned int l = 1;
	  while ( true ) {
	    hierarchy[l++]->remove(v);
	    if (u == NULL) break; 
	    if(l >= ag_hierarchy_2__maxlevel) { break; }
	    v = u;
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
    vertex->set_down(previous); // link with level above
    previous->set_up(vertex);
    previous = vertex;
    level++;
  }
  return first;
}

template<class Gt, class Agds, class LTag>
void 
Apollonius_graph_hierarchy_2<Gt,Agds,LTag>::
remove(Vertex_handle v)
{
  CGAL_triangulation_precondition( v != Vertex_handle());
  CGAL_triangulation_precondition( !is_infinite(v));

  // get the hidden circles
  typename Apollonius_graph_base::Site_list wp_list;
  typename Vertex::Hidden_sites_iterator wpit;

  for (wpit = v->hidden_sites_begin();
       wpit != v->hidden_sites_end(); ++wpit) {
    wp_list.push_back(*wpit);
  }
  v->clear_hidden_sites_container();

  // do the actual removal
  Vertex_handle u = v->up();
  unsigned int l = 0;
  while ( true ) {
    hierarchy[l++]->remove(v);
    if (u == NULL) break; 
    if(l >= ag_hierarchy_2__maxlevel) break;
    v = u;
    u = v->up();
  }

  insert(wp_list.begin(), wp_list.end());
}


template<class Gt, class Agds, class LTag>
typename Apollonius_graph_hierarchy_2<Gt,Agds,LTag>::Vertex_handle 
Apollonius_graph_hierarchy_2<Gt,Agds,LTag>::
nearest_neighbor(const Point_2& p) const
{
  Vertex_handle vnear[ag_hierarchy_2__maxlevel];
  nearest_neighbor(p, vnear);
  return vnear[0];
}



template<class Gt, class Agds, class LTag>
void
Apollonius_graph_hierarchy_2<Gt,Agds,LTag>::
swap(Apollonius_graph_hierarchy_2<Gt,Agds,LTag>& agh)
{
  Ag_base* temp;
  Ag_base::swap(agh);
  for(unsigned int i = 1; i < ag_hierarchy_2__maxlevel; ++i){
    temp = hierarchy[i];
    hierarchy[i] = agh.hierarchy[i];
    agh.hierarchy[i]= temp;
  }
}




template<class Gt, class Agds, class LTag>
void
Apollonius_graph_hierarchy_2<Gt,Agds,LTag>::
nearest_neighbor(const Point_2& p,
		 Vertex_handle vnear[ag_hierarchy_2__maxlevel])
  const
{
  Vertex_handle nearest = 0;
  unsigned int level  = ag_hierarchy_2__maxlevel;

  // find the highest level with enough vertices
  while ( hierarchy[--level]->number_of_vertices() 
	  < ag_hierarchy_2__minsize ) {
    if ( !level ) break;  // do not go below 0
  }
  for (unsigned int i = level+1; i < ag_hierarchy_2__maxlevel; ++i) {
    vnear[i] = 0;
  }

  while ( level > 0 ) {
    vnear[level] = nearest =
      hierarchy[level]->nearest_neighbor(p, nearest);  

    CGAL_assertion( !hierarchy[level]->is_infinite(vnear[level]) );
    // go at the same vertex on level below
    nearest =  nearest->down();
    --level;
  }
  vnear[0] =
    hierarchy[level]->nearest_neighbor(p, nearest);  // at level 0
}

template<class Gt, class Agds, class LTag>
int
Apollonius_graph_hierarchy_2<Gt,Agds,LTag>::
random_level()
{
  unsigned int l = 0;
  while (1) {
    if ( random(ag_hierarchy_2__ratio) ) break;
    ++l;
  }
  if (l >= ag_hierarchy_2__maxlevel)
    l = ag_hierarchy_2__maxlevel -1;
  return l;
}

CGAL_END_NAMESPACE


#endif // CGAL_APOLLONIUS_GRAPH_HIERARCHY_2_C
