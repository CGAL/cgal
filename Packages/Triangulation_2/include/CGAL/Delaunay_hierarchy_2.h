// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
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
// file          : include/CGAL/Delaunay_hierarchy_2.h
// package       : Triangulation
// source        : 
// revision      : 
// revision_date : 
// package       : Delaunay hierarchy
// author(s)     : Olivier Devillers <Olivivier.Devillers@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef DELAUNAY_HIERARCHY_2_H
#define DELAUNAY_HIERARCHY_2_H

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Random.h>

CGAL_BEGIN_NAMESPACE

template < class Traits >
class Delaunay_hierarchy_vertex_base_2
 : public Triangulation_vertex_base_2<Traits>
{
 typedef Triangulation_vertex_base_2<Traits> Base;
 public:
  Delaunay_hierarchy_vertex_base_2()
    : Base(), _up(0), _down(0)
    {}
  Delaunay_hierarchy_vertex_base_2(const Point & p, void* f)
    : Base(p,f), _up(0), _down(0)
    {}
  Delaunay_hierarchy_vertex_base_2(const Point & p)
    : Base(p), _up(0), _down(0)
    {}

 public:  // for use in Delaunay_hierarchy only
  //  friend class Delaunay_hierarchy_2;
  void* up() {return _up;}
  void* down() {return _down;}
  void set_up(void *u) {_up=u;}
  void set_down(void *d) {if (this) _down=d;}

 private:
  void* _up;    // same vertex one level above
  void* _down;  // same vertex one level below
};

// parameterization of the Delaunay hierarchy
const float Delaunay_hierarchy_2__ratio    = 30.0;
const int   Delaunay_hierarchy_2__minsize  = 20;
const int   Delaunay_hierarchy_2__maxlevel = 5;
// maximal number of points is 30^5 = 24 millions !

template < class Traits, class Tds>
class Delaunay_hierarchy_2
: public Delaunay_triangulation_2<Traits,Tds>
{
 public:
  typedef Delaunay_triangulation_2<Traits,Tds> Delaunay;
  typedef typename Delaunnay::Vertex_handle    Vertex_handle;
  typedef typename Delaunnay::Face_handle      Face_handle;

 private:
  // here is the stack of triangulations which form the hierarchy
  Delaunay*   hierarchy[Delaunay_hierarchy_2__maxlevel];
  Random random; // random generator

public:
  Delaunay_hierarchy_2(const Traits& traits = Traits());
  Delaunay_hierarchy_2(const Delaunay_hierarchy_2& tr);

  Delaunay_hierarchy_2 &operator=(const  Delaunay_hierarchy_2& tr);

  //Helping
  void copy_triangulation(const Delaunay_hierarchy_2 &tr);
  void swap(Delaunay_hierarchy_2 &tr);
  void clear();

  // CHECKING
  bool is_valid() const;

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
  Face_handle  locate(const Point& p, Locate_type& lt,int& li) const;
  Face_handle  locate(const Point& p);

private:
  void  locate(const Point& p,
	       Locate_type& lt,
	       int& li,
	       Face_handle pos[Delaunay_hierarchy_2__maxlevel]) const;
  int random_level();

};


template <class Traits, class Tds >
Delaunay_hierarchy_2<Traits,Tds>::
Delaunay_hierarchy_2(const Traits& traits)
  : Delaunay(traits), random((long)0)
{ 
  hierarchy[0] = this; 
  for(int i=1;i<Delaunay_hierarchy_2__maxlevel;++i)
    hierarchy[i] = new Delaunay();
}


// copy constructor duplicates vertices and faces
template <class Gt, class Tds >
Delaunay_hierarchy_2<Gt, Tds>::
Delaunay_hierarchy_2(const Delaunay_hierarchy_2<Gt,Tds> &tr)
    : Delaunay(), random((long)0)
{ 
  // create an empty triangulation to be able to delete it !
  hierarchy[0] = this; 
  for(int i=1;i<Delaunay_hierarchy_2__maxlevel;++i)
    hierarchy[i] = new Delaunay();
  copy_triangulation(tr);
} 
 

//Assignement
template <class Gt, class Tds >
Delaunay_hierarchy_2<Gt, Tds> &
Delaunay_hierarchy_2<Gt, Tds>::
operator=(const Delaunay_hierarchy_2<Gt,Tds> &tr)
{
  copy_triangulation(tr);
  return *this;
}

template <class Gt, class Tds >
void
Delaunay_hierarchy_2<Gt, Tds>::   
copy_triangulation(const Delaunay_hierarchy_2<Gt,Tds> &tr)
{
  std::map< const void*, void*, std::less<const void*> > V;

  for(int i=0;i<Delaunay_hierarchy_2__maxlevel;++i)
    hierarchy[i]->copy_triangulation(*tr.hierarchy[i]);
  //up and down have been copied in straightforward way
  // compute a map at lower level
  for( Vertex_iterator it=hierarchy[0]->vertices_begin(); 
       it != hierarchy[0]->vertices_end(); ++it) {
    if (it->up()) V[ ((Vertex*)(it->up()))->down() ] = &(*it);
  }
  for(int i=1;i<Delaunay_hierarchy_2__maxlevel;++i) {
    for( Vertex_iterator it=hierarchy[i]->vertices_begin(); 
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

template <class Gt, class Tds >
void
Delaunay_hierarchy_2<Gt, Tds>:: 
swap(Delaunay_hierarchy_2<Gt,Tds> &tr)
{
  Delaunay** h= hierarchy;
  hierarchy = tr.hierarchy;
  tr.hierarchy = h;
}

template <class Gt, class Tds >
void
Delaunay_hierarchy_2<Gt, Tds>:: 
clear()
{
        for(int i=0;i<Delaunay_hierarchy_2__maxlevel;++i)
	hierarchy[i]->clear();
}


template <class Gt, class Tds >
bool
Delaunay_hierarchy_2<Gt, Tds>:: 
is_valid()
{
  bool result = true;
  //verify correctness of triangulation at all levels
  for(int i=0;i<Delaunay_hierarchy_2__maxlevel;++i)
	result = result && hierarchy[i]->is_valid();
  //verify that lower level has no down pointers
  for( Vertex_iterator it=hierarchy[0]->vertices_begin(); 
       it != hierarchy[0]->vertices_end(); ++it) 
    result = result && ( it->down() == 0 );
  //verify that other levels has down pointer and reciprocal link is 
fine
  for(int i=1;i<Delaunay_hierarchy_2__maxlevel;++i)
    for( Vertex_iterator it=hierarchy[i]->vertices_begin(); 
	 it != hierarchy[i]->vertices_end(); ++it) 
      result = result && ( ((Vertex*)((Vertex*)it->down())->up()) ==
			   &(*it) );
  return result;
}

  
template <class Traits, class Tds >
Delaunay_hierarchy_2<Traits,Tds>::Vertex_handle
Delaunay_hierarchy_2<Traits,Tds>::
insert(const Point &p)
{
  int vertex_level = random_level();
  Locate_type lt;
  int i;
  // locate using hierarchy
  Face_handle positions[Delaunay_hierarchy_2__maxlevel];
  locate(p,lt,i,positions);
  //insert at level 0
  Vertex_handle vertex=hierarchy[0]->insert(p,positions[0]);
  Vertex_handle previous=vertex;
  Vertex_handle first = vertex;
      
  int level  = 1;
  while (level <= vertex_level ){
    vertex=hierarchy[level]->insert(p,positions[level]);
    vertex->set_down((void *) &*previous);// link with level above
    previous->set_up((void *) &*vertex);
    previous=vertex;
    level++;
  }
  return first;
}

template <class Traits, class Tds >
inline
Delaunay_hierarchy_2<Traits,Tds>::Vertex_handle
Delaunay_hierarchy_2<Traits,Tds>::
push_back(const Point &p)
{
  return insert(p);
}

template <class Traits, class Tds >
void 
Delaunay_hierarchy_2<Traits,Tds>::
remove(Vertex_handle v )
{
  void * u=v->up();
  int l = 0 ;
  while(1){
    hierarchy[l++]->remove(v);
    if (!u) break; 
    if(l>Delaunay_hierarchy_2__maxlevel) break;
    v=(Vertex*)u; u=v->up();
  }
}

template <class Traits, class Tds >
inline void 
Delaunay_hierarchy_2<Traits,Tds>::
remove_degree_3(Vertex_handle v )
{
  remove(v);
}

template <class Traits, class Tds >
inline void 
Delaunay_hierarchy_2<Traits,Tds>::
remove_first(Vertex_handle v )
{
  remove(v);
}

template <class Traits, class Tds >
inline void 
Delaunay_hierarchy_2<Traits,Tds>::
remove_second(Vertex_handle v )
{
  remove(v)
}

template <class Traits, class Tds >
Delaunay_hierarchy_2<Traits,Tds>::Face_handle 
Delaunay_hierarchy_2<Traits,Tds>::
locate(const Point& p, Locate_type& lt, int& li) const
{
  Face_handle positions[Delaunay_hierarchy_2__maxlevel];
  locate(p,lt,li,positions);
  return positions[0];
}

template <class Traits, class Tds >
Delaunay_hierarchy_2<Traits,Tds>::Face_handle 
Delaunay_hierarchy_2<Traits,Tds>::
locate(const Point& p) const
{
  Locate_type lt;
  int li;
  return locate(p, lt, li);
}

template <class Traits, class Tds >
void
Delaunay_hierarchy_2<Traits,Tds>::
locate(const Point& p,
       Locate_type& lt,
       int& li,
       Face_handle pos[Delaunay_hierarchy_2__maxlevel]) const
{
  std::cerr << "locate of Delaunay_hierarchy_2" << std::endl;
  Face_handle position;
  Vertex_handle nearest;
  int level  = Delaunay_hierarchy_2__maxlevel;

  // find the highest level with enough vertices
  while (hierarchy[--level]->number_of_vertices() 
	 < Delaunay_hierarchy_2__minsize){
    if ( ! level) break;  // do not go below 0
  }
  for (int i=level+1; i<Delaunay_hierarchy_2__maxlevel;++i) pos[i]=0;
  while(level > 0) {
    pos[level]=position=hierarchy[level]->locate(p,position);  
    // locate at that level from "position"
    // result is stored in "position" for the next level
    // find the nearest between vertices 0 and 1
    if (hierarchy[level]->is_infinite(position->vertex(0)))
      nearest = position->vertex(1);
    else if (hierarchy[level]->is_infinite(position->vertex(1)))
      nearest = position->vertex(0);
    else if (  //geom_traits().
	     cmp_dist_to_point(  p,
				 position->vertex(0)->point(),
				 position->vertex(1)->point())==SMALLER)
      nearest = position->vertex(0);
    else
      nearest = position->vertex(1);
    // compare to vertex 2
    if ( !  hierarchy[level]->is_infinite(position->vertex(2)))
      if (    //geom_traits().
	  cmp_dist_to_point( p,
			     position->vertex(2)->point(),
			     nearest->point())==SMALLER)
	nearest = position->vertex(2);
    // go at the same vertex on level below
    nearest = (Vertex*)( nearest->down() );
    position = nearest->face();                // incident face
    --level;
  }
  pos[0]=hierarchy[level]->locate(p,lt,li,position);  // at level 0
}


template <class Traits, class Tds >
int
Delaunay_hierarchy_2<Traits,Tds>::
random_level()
{
  int l = 0;
  while (1) {
    if ( random(Delaunay_hierarchy_2__ratio) ) break;
    ++l;
  }
  if (l >= Delaunay_hierarchy_2__maxlevel)
    l = Delaunay_hierarchy_2__maxlevel -1;
  return l;
}

CGAL_END_NAMESPACE
#endif //DELAUNAY_HIERARCHY_2_H
