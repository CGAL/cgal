// ============================================================================
//
// Copyright (c) 1998, 2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Triangulation_hierarchy_3.h
// revision      : $Revision$
// revision_date : $Date$
// package       : Triangulation3
// author(s)     : Olivier Devillers, Sylvain Pion
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_HIERARCHY_3_H
#define CGAL_TRIANGULATION_HIERARCHY_3_H

#include <CGAL/Random.h>

CGAL_BEGIN_NAMESPACE

template < class Vbb>
class Triangulation_hierarchy_vertex_base_3
 : public Vbb
{
 public:
  typedef Vbb V_Base;
  typedef typename V_Base::Point   Point;

  Triangulation_hierarchy_vertex_base_3()
    : V_Base(), _up(0), _down(0)
    {}
  Triangulation_hierarchy_vertex_base_3(const Point & p, void* f)
    : V_Base(p,f), _up(0), _down(0)
    {}
  Triangulation_hierarchy_vertex_base_3(const Point & p)
    : V_Base(p), _up(0), _down(0)
    {}

 public:  // for use in Triangulation_hierarchy only
  //  friend class Triangulation_hierarchy_3;
  void* up() const {return _up;}
  void* down() const {return _down;}
  void set_up(void *u) {_up=u;}
  void set_down(void *d) {if (this) _down=d;}

 private:
  void* _up;    // same vertex one level above
  void* _down;  // same vertex one level below
};

// parameterization of the  hierarchy
//const float Triangulation_hierarchy_3__ratio    = 30.0;
const int Triangulation_hierarchy_3__ratio    = 30;
const int Triangulation_hierarchy_3__minsize  = 20;
const int Triangulation_hierarchy_3__maxlevel = 5;
// maximal number of points is 30^5 = 24 millions !

template < class Tr>
class Triangulation_hierarchy_3
: public Tr
{
 public:
  typedef Tr                        Tr_Base;
  typedef typename Tr_Base::Geom_traits        Geom_traits;
  typedef typename Geom_traits::Point_3        Point;
  typedef typename Tr_Base::Vertex_handle     Vertex_handle;
  typedef typename Tr_Base::Cell_handle       Cell_handle;
  typedef typename Tr_Base::Vertex_iterator   Vertex_iterator;
  typedef typename Tr_Base::Vertex           Vertex;
  typedef typename Tr_Base::Locate_type       Locate_type;

 private:
  // here is the stack of triangulations which form the hierarchy
  Tr_Base*   hierarchy[Triangulation_hierarchy_3__maxlevel];
  Random random; // random generator

public:
  Triangulation_hierarchy_3(const Geom_traits& traits = Geom_traits());
  Triangulation_hierarchy_3(const Triangulation_hierarchy_3& tr);

  Triangulation_hierarchy_3 &operator=(const  Triangulation_hierarchy_3& tr);
  ~Triangulation_hierarchy_3();

  //Helping
  void copy_triangulation(const Triangulation_hierarchy_3 &tr);
  void swap(Triangulation_hierarchy_3 &tr);
  void clear();

  // CHECKING
  bool is_valid(bool verbose = false, int level = 0) const;

  // INSERT REMOVE
  Vertex_handle insert(const Point &p);
  Vertex_handle push_back(const Point &p);
 
  template < class InputIterator >
  int insert(InputIterator first, InputIterator last)
    {
      int n = number_of_vertices();
      while(first != last){
	insert(*first);
	++first;
      }
      return number_of_vertices() - n;
    }

  void remove_degree_3(Vertex_handle  v);
  void remove_first(Vertex_handle  v);
  void remove_second(Vertex_handle v);
  void remove(Vertex_handle  v);

  //LOCATE
  Cell_handle locate(const Point& p, Locate_type& lt, int& li, int& lj) const;
  Cell_handle locate(const Point& p) const;

private:
  void locate(const Point& p, Locate_type& lt, int& li, int& lj,
	      Cell_handle pos[Triangulation_hierarchy_3__maxlevel]) const;
  int random_level();

  // added to make the test program of usual triangulations work
  // undocumented
public:
  
  Vertex_handle insert(const Point &p, Cell_handle start)
  {
    return Tr_Base::insert(p, start);
  }

  Vertex_handle insert(const Point& p, Locate_type lt, Cell_handle loc, int li)
  {
    return Tr_Base::insert(p);
  }

  Cell_handle locate(const Point& p, Locate_type& lt, int& li, int& lj,
		     Cell_handle start) const
  {
    return Tr_Base::locate(p, lt, li, lj, start);
  }
};


template <class Tr >
Triangulation_hierarchy_3<Tr>::
Triangulation_hierarchy_3(const Geom_traits& traits)
  : Tr_Base(traits), random((long)0)
{ 
  hierarchy[0] = this; 
  for(int i=1;i<Triangulation_hierarchy_3__maxlevel;++i)
    hierarchy[i] = new Tr_Base(traits);
}

// copy constructor duplicates vertices and cells
template <class Tr>
Triangulation_hierarchy_3<Tr>::
Triangulation_hierarchy_3(const Triangulation_hierarchy_3<Tr> &tr)
    : Tr_Base(), random((long)0)
{ 
  // create an empty triangulation to be able to delete it !
  hierarchy[0] = this; 
  for(int i=1;i<Triangulation_hierarchy_3__maxlevel;++i)
    hierarchy[i] = new Tr_Base(tr.geom_traits());
  copy_triangulation(tr);
} 
 

//Assignement
template <class Tr>
Triangulation_hierarchy_3<Tr> &
Triangulation_hierarchy_3<Tr>::
operator=(const Triangulation_hierarchy_3<Tr> &tr)
{
  copy_triangulation(tr);
  return *this;
}


template <class Tr>
void
Triangulation_hierarchy_3<Tr>::   
copy_triangulation(const Triangulation_hierarchy_3<Tr> &tr)
{
  std::map< const void*, void* > V;

  for(int i=0; i<Triangulation_hierarchy_3__maxlevel; ++i)
    hierarchy[i]->copy_triangulation(*tr.hierarchy[i]);

  // up and down have been copied in straightforward way
  // compute a map at lower level

  for( Vertex_iterator it=hierarchy[0]->finite_vertices_begin(); 
       it != hierarchy[0]->vertices_end(); ++it)
    if (it->up())
      V[ ((Vertex*)(it->up()))->down() ] = &(*it);

  for(int j=1; j<Triangulation_hierarchy_3__maxlevel; ++j) {
    for( Vertex_iterator it=hierarchy[j]->finite_vertices_begin();
	 it != hierarchy[j]->vertices_end(); ++it) {
	// down pointer goes in original instead in copied triangulation
	it->set_down(V[it->down()]);
	// make reverse link
	((Vertex*)(it->down()))->set_up( &(*it) );
	// make map for next level
	if (it->up())
	    V[ ((Vertex*)(it->up()))->down() ] = &(*it);
    }
  }
}

template <class Tr>
void
Triangulation_hierarchy_3<Tr>:: 
swap(Triangulation_hierarchy_3<Tr> &tr)
{
//   Tr_Base** h= hierarchy;
//   hierarchy = tr.hierarchy;
//   tr.hierarchy = h;
  Tr_Base* temp;
  Tr_Base::swap(tr);
  for(int i= 1; i<Triangulation_hierarchy_3__maxlevel; ++i){
    temp = hierarchy[i];
    hierarchy[i] = tr.hierarchy[i];
    tr.hierarchy[i]= temp;
  }
}

template <class Tr>
Triangulation_hierarchy_3<Tr>:: 
~Triangulation_hierarchy_3()
{
  clear();
  for(int i= 1; i<Triangulation_hierarchy_3__maxlevel; ++i){ 
    delete hierarchy[i];
  }
}

template <class Tr>
void
Triangulation_hierarchy_3<Tr>:: 
clear()
{
        for(int i=0;i<Triangulation_hierarchy_3__maxlevel;++i)
	hierarchy[i]->clear();
}


template <class Tr>
bool
Triangulation_hierarchy_3<Tr>:: 
is_valid(bool verbose, int level) const
{
  bool result = true;
  
  //verify correctness of triangulation at all levels
  for(int i=0; i<Triangulation_hierarchy_3__maxlevel; ++i)
	result = result && hierarchy[i]->is_valid(verbose, level);

  //verify that lower level has no down pointers
  for( Vertex_iterator it = hierarchy[0]->finite_vertices_begin(); 
       it != hierarchy[0]->vertices_end(); ++it) 
    result = result && ( it->down() == 0 );

  //verify that other levels has down pointer and reciprocal link is fine
  for(int i=1; i<Triangulation_hierarchy_3__maxlevel; ++i)
    for( Vertex_iterator it = hierarchy[i]->finite_vertices_begin(); 
	 it != hierarchy[i]->vertices_end(); ++it) 
      result = result && 
	       ( ((Vertex*)((Vertex*)it->down())->up()) ==  &(*it) );

  //verify that other levels has down pointer and reciprocal link is fine
  for(int i=0; i<Triangulation_hierarchy_3__maxlevel-1; ++i)
    for( Vertex_iterator it = hierarchy[i]->finite_vertices_begin(); 
	 it != hierarchy[i]->vertices_end(); ++it) 
      result = result && ( ((Vertex*)it->up() == NULL) ||
	       ( ((Vertex*)((Vertex*)it->up())->down()) ==  &(*it) ));

  return result;
}

  
template <class Tr>
Triangulation_hierarchy_3<Tr>::Vertex_handle
Triangulation_hierarchy_3<Tr>::
insert(const Point &p)
{
  int vertex_level = random_level();
  Locate_type lt;
  int i, j;
  // locate using hierarchy
  Cell_handle positions[Triangulation_hierarchy_3__maxlevel];
  locate(p,lt,i,j,positions);
  //insert at level 0
  Vertex_handle vertex = hierarchy[0]->insert(p,positions[0]);
  Vertex_handle previous = vertex;
  Vertex_handle first = vertex;

  int level = 1;
  while (level <= vertex_level ){
    if (positions[level] != NULL)
      vertex = hierarchy[level]->insert(p, positions[level]);
    else
      vertex = hierarchy[level]->insert(p);
    vertex->set_down((void *) &*previous);// link with level above
    previous->set_up((void *) &*vertex);
    previous=vertex;
    level++;
  }
  return first;
}

template <class Tr>
inline
Triangulation_hierarchy_3<Tr>::Vertex_handle
Triangulation_hierarchy_3<Tr>::
push_back(const Point &p)
{
  return insert(p);
}

template <class Tr>
void 
Triangulation_hierarchy_3<Tr>::
remove(Vertex_handle v )
{
  void * u=v->up();
  int l = 0;
  while(1){
    hierarchy[l++]->remove(v);
    if (!u) break;
    if (l>Triangulation_hierarchy_3__maxlevel) break;
    v=(Vertex*)u; u=v->up();
  }
}

template <class Tr>
inline void 
Triangulation_hierarchy_3<Tr>::
remove_degree_3(Vertex_handle v )
{
  remove(v);
}

template <class Tr>
inline void 
Triangulation_hierarchy_3<Tr>::
remove_first(Vertex_handle v )
{
  remove(v);
}

template <class Tr>
inline void 
Triangulation_hierarchy_3<Tr>::
remove_second(Vertex_handle v )
{
  remove(v);
}

template <class Tr>
inline
Triangulation_hierarchy_3<Tr>::Cell_handle 
Triangulation_hierarchy_3<Tr>::
locate(const Point& p, Locate_type& lt, int& li, int& lj) const
{
  Cell_handle positions[Triangulation_hierarchy_3__maxlevel];
  locate(p,lt,li,lj,positions);
  return positions[0];
}

template <class Tr>
inline
Triangulation_hierarchy_3<Tr>::Cell_handle 
Triangulation_hierarchy_3<Tr>::
locate(const Point& p) const
{
  Locate_type lt;
  int li, lj;
  return locate(p, lt, li, lj);
}

template <class Tr>
void
Triangulation_hierarchy_3<Tr>::
locate(const Point& p, Locate_type& lt, int& li, int& lj,
       Cell_handle pos[Triangulation_hierarchy_3__maxlevel]) const
{
  Cell_handle position;
  Vertex_handle nearest;
  int level  = Triangulation_hierarchy_3__maxlevel;
  typename Geom_traits::Less_distance_to_point_3 
    closer = geom_traits().less_distance_to_point_3_object(p);

  // find the highest level with enough vertices
  while (hierarchy[--level]->number_of_vertices() 
	 < Triangulation_hierarchy_3__minsize){
    if ( ! level) break;  // do not go below 0
  }
  for (int i=level+1; i<Triangulation_hierarchy_3__maxlevel;++i) pos[i]=0;
  while(level > 0) {
    pos[level]=position=hierarchy[level]->locate(p,position);
    // locate at that level from "position"
    // result is stored in "position" for the next level
    // find the nearest between vertices 0 and 1
    if (hierarchy[level]->is_infinite(position->vertex(0)))
      nearest = position->vertex(1);
    else if (hierarchy[level]->is_infinite(position->vertex(1)))
      nearest = position->vertex(0);
    else if ( closer(position->vertex(0)->point(),
		     position->vertex(1)->point()))
      nearest = position->vertex(0);
    else
      nearest = position->vertex(1);
    // compare to vertex 2
    if ( dimension() >= 2
      && ( ! hierarchy[level]->is_infinite(position->vertex(2)))
      && ( closer( position->vertex(2)->point(), nearest->point())))
	nearest = position->vertex(2);
    // compare to vertex 3
    if ( dimension() == 3
      && ( ! hierarchy[level]->is_infinite(position->vertex(3)))
      && ( closer( position->vertex(3)->point(), nearest->point())))
	nearest = position->vertex(3);
    // go at the same vertex on level below
    nearest = (Vertex*)( nearest->down() );
    position = nearest->cell();                // incident cell
    --level;
  }
  pos[0]=hierarchy[level]->locate(p,lt,li,lj,position);  // at level 0
}


template <class Tr>
int
Triangulation_hierarchy_3<Tr>::
random_level()
{
  int l = 0;
  while (1) {
    if ( random(Triangulation_hierarchy_3__ratio) ) break;
    ++l;
  }
  if (l >= Triangulation_hierarchy_3__maxlevel)
    l = Triangulation_hierarchy_3__maxlevel -1;
  return l;
}

CGAL_END_NAMESPACE

#endif // CGAL_TRIANGULATION_HIERARCHY_3_H
