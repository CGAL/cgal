// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
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
// file          : Triangulation/include/CGAL/Triangulation_ds_iterators_2.h
// source        : $Source$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Mariette Yvinec
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_DS_ITERATORS_2_H
#define CGAL_TRIANGULATION_DS_ITERATORS_2_H



#include <pair.h>

#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_short_names_2.h>

//template< class Gt , class Vb, class Fb> 
//class CGAL_Triangulation_default_data_structure_2;


//template < class Gt , class Vb, class Fb>
template <class Tds>
class CGAL_Triangulation_ds_iterator_base_2
{
public:
  typedef typename Tds::Geom_traits Geom_traits;
  typedef typename Tds::Vertex Vertex;
  typedef typename Tds::Face  Face;
  typedef typename Tds::Edge Edge;

  CGAL_Triangulation_ds_iterator_base_2()
     : _tds(NULL), pos(NULL)
        {}

  CGAL_Triangulation_ds_iterator_base_2(Tds* tds)
     : _tds(tds), pos(NULL)
  {
    if(_tds->number_of_vertices() < 2) {
      return;
    }
    pos = _tds->infinite_vertex()->face();
  }

 CGAL_Triangulation_ds_iterator_base_2(Tds* tds, int i)
     :  _tds(tds), pos(NULL)
 {}

  CGAL_Triangulation_ds_iterator_base_2(Tds* tds, Face* f)
    : _tds(tds), pos(f)
  {}
   
protected:
        Tds*  _tds;
        Face* pos;
        
public:       
        static
        int
        ccw(int i)
        {
            return (i+1) % 3;
        }
        
        
        static
        int
        cw(int i)
        {
            return (i+2) % 3;
        }
        

        void
        increment()
        {
	  if (_tds->dimension() == 1 || _tds->dimension() == 0){
	    pos = pos->neighbor(0);
	    return;
	  }

            int max = maximum(pos);
            Face* next=pos->neighbor(max);         // tentative first child
            Face* parent;
            int max2 = maximum(next);
            if ( next->neighbor(cw(max2))== pos){
              // next is the first child of pos
              pos = next;
              return;
            }
            // look for the second child of an ancestor of pos
            next=pos->neighbor(ccw(max));          // tentatative second child
            while (1){
                max2 = maximum(next);
                if ( next->neighbor(cw(max2))== pos) // next is the second child of pos
                    { pos = next; return;}
                while (1){
                    parent = pos->neighbor(cw(max));        // go to parent
                    max = maximum(parent);
                    next=parent->neighbor(ccw(max));       // tentatative second child
                    if (next==pos)              // already coming back from this child
                        { pos = parent; continue; }
                    else
                        { pos = parent; break; }
                }
            }
        }
        

        void
        decrement()
        {
	  if(_tds->dimension() == 0){
	    pos = pos->neighbor(0);
	    return;
	  }

	    if (_tds->dimension() == 1){
	    pos = pos->neighbor(1);
	    return;
	  }

	    // dimension() ==2
	    int max = maximum(pos);
            Face* next=pos->neighbor(cw(max));     // parent of pos
            int max2 = maximum(next);
            if ( next->neighbor(max2) == pos)      // pos is the first child of next
                { pos = next; return;}
            pos = next->neighbor(max2);         // tentative first child
            max = maximum(pos);
            if ( pos->neighbor(cw(max))!=next)   // pos is not the first child of next
                { pos = next; return;}
            // look for "last" node in first subtree of next
            while (1){
                next = pos->neighbor(ccw(max));       // tentatative second child
                max2 = maximum(next);
                if (next->neighbor(cw(max2))!= pos){
                    //next is not the second child of pos
                    next=pos->neighbor(max);         // tentative first child
                    max2 = maximum(next);
                    if ( next->neighbor(cw(max2))!= pos){
                        //next is not first child of pos
                        return;
                    }
                }
                pos=next;
                max=max2;
            }
        }

        int maximum(const Face* f) const
        {
            if ( _tds->is_infinite(f) ){
                return f->index( _tds->infinite_vertex() );
            }
            if(_tds->geom_traits().compare_y(f->vertex(0)->point(),
                                               f->vertex(1)->point()) == CGAL_SMALLER)
                //  v0 < v1
                if(_tds->geom_traits().compare_y(f->vertex(2)->point(),
                                                   f->vertex(1)->point())==CGAL_SMALLER)
                    //  v0,v2 < v1
                    { return 1; }
                else
                    //  v0 < v1 <= v2
                    { return 2; }
            else
                //  v1 <= v0
        
                if(_tds->geom_traits().compare_y(f->vertex(1)->point(),
                                                   f->vertex(2)->point())!=CGAL_SMALLER)
                    //  v2 <= v1 <= v0
                    if(_tds->geom_traits().compare_y(
                                                   f->vertex(0)->point(),
                                                   f->vertex(1)->point()) == CGAL_EQUAL)
                        //  v2 <= v1 == v0
                        { return 1; }
                    else
                        //  v2 <= v1 < v0
                        { return 0; }
                else
                    //  v1<=v0, v1<v2
                    if(_tds->geom_traits().compare_y(f->vertex(0)->point(),
                                                       f->vertex(2)->point())
                       ==CGAL_SMALLER)
                        //  v1 <= v0 < v2
                        { return 2; }
                    else
                        //  v1 < v2 <= v0
                        { return 0; }
        
        }

};

// the following iterator visit all the Tds faces
// whatever may be the dimensionality of those faces
template<class Tds>
class CGAL_Triangulation_ds_face_iterator_2
  : public CGAL_Triangulation_ds_iterator_base_2<Tds>,
    public bidirectional_iterator<typename Tds::Face, ptrdiff_t>
{
public:
  typedef typename Tds::Geom_traits Geom_traits;
  typedef typename Tds::Vertex Vertex;
  typedef typename Tds::Face  Face;
  typedef typename Tds::Edge Edge;

  typedef CGAL_Triangulation_ds_iterator_base_2<Tds> Iterator_base;
  typedef CGAL_Triangulation_ds_face_iterator_2<Tds> Face_iterator;

        CGAL_Triangulation_ds_face_iterator_2()
            : Iterator_base()
        {}
        CGAL_Triangulation_ds_face_iterator_2(Tds * tds)
            : Iterator_base(tds)
        {}

        CGAL_Triangulation_ds_face_iterator_2(Tds* tds, int i)
	  : Iterator_base(tds,i)
        {}

        CGAL_Triangulation_ds_face_iterator_2(const Face_iterator& fi)
          : Iterator_base(fi._tds, fi.pos)
        {}
        
        Face_iterator&
        operator=(const Face_iterator& fi)
        {
	  pos = fi.pos;
	  _tds = fi._tds;
	  return *this;
        }

        bool
        operator==(const Face_iterator& fi) const
        {
            return ((pos == fi.pos )&&(_tds==fi._tds));
        }
        
        
        bool
        operator!=(const Face_iterator& fi) const
        {
            return !(*this == fi);
        }

        Face_iterator&
        operator++()
        {
            if ( pos == NULL ){
                return *this;    //  past-the-end iterator cannot advance
            }                    // include the case dimension()==0

	    increment();
	    if ( pos == _tds->infinite_face() ){
	      pos = NULL;  // complete tour
	    }  
	    return *this;           // next finite triangle found
        }

        Face_iterator&
        operator--()
        {
	  if (dimension() == 0) {return this;}

	  if (pos ==  _tds->infinite_face()) {
	    pos == NULL; //first face, can't decrement
	    return *this;
	  }
	  if (pos == NULL) { //  past the end iterator, can decrease
	    *this = Face_iterator(_tds); //first face, decrement will give the last one if any
	  }
	  decrement(); 
	  return *this;
	}
        
        Face_iterator
        operator++(int)
        {
            Face_iterator tmp(*this);
            ++(*this);
            return tmp;
        }
        
        
        Face_iterator
        operator--(int)
        {
            Face_iterator tmp(*this);
            --(*this);
            return tmp;
        }
        
  inline Face& operator*() const
  {
        return *pos;
  }
    
  inline Face*  operator->() const
  {
        return pos;
  }


};



template < class Tds>
class CGAL_Triangulation_ds_vertex_iterator_2
: public CGAL_Triangulation_ds_iterator_base_2<Tds>,
  public bidirectional_iterator<typename Tds::Vertex, ptrdiff_t>
{
public:
  typedef typename Tds::Geom_traits Geom_traits;
  typedef typename Tds::Vertex Vertex;
  typedef typename Tds::Face  Face;
  typedef typename Tds::Edge Edge;

  typedef CGAL_Triangulation_ds_iterator_base_2<Tds> Iterator_base;
  typedef CGAL_Triangulation_ds_vertex_iterator_2<Tds> Vertex_iterator;


private :
  int index;

public:
     CGAL_Triangulation_ds_vertex_iterator_2()
            : Iterator_base(), index(0)
        {}
    
    
    CGAL_Triangulation_ds_vertex_iterator_2(Tds * tds)
        :  Iterator_base(tds), index(0)
    {
        switch( tds->number_of_vertices() ){
        case 0: // past-the-end
            pos = NULL;
            return;
        case 1:
	  pos = (Face*)1 ; // different from any pointer;
	  return;         // the iterator must "point" to the only vertex of the triangulation
	
        default:
	  pos = tds->infinite_face();
	  index = 0;
	  while ( ! associated_vertex()){
	    increment();
	  }
	  return;
	}
    }
    
    CGAL_Triangulation_ds_vertex_iterator_2(Tds* tds, int i)
            : Iterator_base(tds,i), index(0)
    {}

private:
  void   increment()
  {
    if ( index==_tds->dimension()) {Iterator_base::increment(); index = 0;}
    else { index +=1; }
    return;
  }

  void decrement()
  {
    if (index == 0) { 
      Iterator_base::decrement();
      index = _tds->dimension();
    }
    else {index -= 1;}
    return;
  }

bool associated_vertex()
  {
    return ( pos->vertex(index)->face() == pos ); // marche en toute dimension
  }

public:
    Vertex_iterator&
    operator++()
    {
        if (pos==NULL){
	  return *this;            // cannot advance past-the-end iterator
        }                          // include the case number_of_vertices() == 0
	if (pos==(Face*)1){
	  pos == NULL;             //number_of_vertices() == 1
	  return *this; 
	}

        do{
            increment();
            if ( pos == _tds->infinite_face() && index == 0){
                pos = NULL;   // complete tour
		return *this;
            }
        }while ( ! associated_vertex() );
        return *this;
    }
    
public:    
    Vertex_iterator&
    operator--()
    {
        switch(_tds->number_of_vertices()) {
        case 0:
	  return *this;
        case 1:
            if(pos == NULL){
                *this = Vertex_iterator(_tds);
            } else {
                pos = NULL;
            }
	    return *this;            // can decrease past-the-end iterator
        default:
	  if (pos == _tds->infinite_face() && index == 0){ //first position, cannot decrease
	    pos = NULL;
	    return *this;
	  }
	  if (pos==NULL){ // can decrease past-the-end iterator
	    *this = Vertex_iterator(_tds);
	  }
	    
	  do{
	    decrement();
	    if ( pos == _tds->infinite_face() && index == 0){
	      pos = NULL;   // complete tour
	      return *this;
	    }
	  }while ( ! associated_vertex());
	  return *this;
        }
    }
    
    
    Vertex_iterator operator++(int)
    {
        Vertex_iterator tmp(*this);
        ++(*this);
        return tmp;
    }
    
    
    Vertex_iterator operator--(int)
    {
        Vertex_iterator tmp(*this);
        --(*this);
        return tmp;
    }
    
    bool operator==(const Vertex_iterator& fi) const
   {
    if ((pos==fi.pos)&&(_tds==fi._tds)){
      if (pos==NULL) {
	return true;
      }else{
	return (index == fi.index);
      }
    }
    else{
      return false;
    }
  }
 
    
    bool operator!=(const Vertex_iterator& fi) const
    {
        return !(*this == fi);
    }
    
  inline Vertex& operator*() const
  {
        return *(pos->vertex(index));
  }
    
  inline Vertex*  operator->() const
  {
        return pos->vertex(index);
  }

    
   
};



template <class Tds>
class CGAL_Triangulation_ds_edge_iterator_2
 : public CGAL_Triangulation_ds_iterator_base_2<Tds>,
   public bidirectional_iterator<typename Tds::Edge, ptrdiff_t>
{
public:
  typedef typename Tds::Geom_traits Geom_traits;
  typedef typename Tds::Vertex Vertex;
  typedef typename Tds::Face  Face;
  typedef typename Tds::Edge Edge;

  typedef CGAL_Triangulation_ds_iterator_base_2<Tds> Iterator_base;
  typedef CGAL_Triangulation_ds_edge_iterator_2<Tds> Edge_iterator;

private:
int index;

public:
     CGAL_Triangulation_ds_edge_iterator_2()
            : Iterator_base(),index(0)
        {}
    
    
    CGAL_Triangulation_ds_edge_iterator_2(Tds * tds)
        :  Iterator_base(tds), index(0) 
    {
      if (_tds->dimension()== 0){
	pos = NULL;                  // there is no edge
	return;
      }
      pos = _tds->infinite_face();
      if (_tds->dimension() == 1) {index = 2;}
      else {index = 0;}
      while ( !associated_edge() ){
	increment();
      }
    }
     
    
    CGAL_Triangulation_ds_edge_iterator_2(Tds* tds, int i)
            : Iterator_base(tds,i),index(0)
    {
       if (_tds->dimension() == 1) {index = 2;}
    }

  bool
  operator==(const Edge_iterator& fi) const
  {
    if ((pos==fi.pos)&&(_tds==fi._tds)){
      if (pos==NULL) {
	return true;
      }else{
	return (index==fi.index);
      }
    }
    else{
      return false;
    }
  }
    
    
    bool
    operator!=(const Edge_iterator& fi) const
    {
        return ! (*this == fi);
    }

private: 
void   increment()
  {
    if (_tds->dimension() == 1) {Iterator_base::increment();}
    else {
      if (index == 2) {
      Iterator_base::increment();
      index =  0;
      }
      else{
	index +=1;
      }
    }
    return;
  }

void decrement()
  {
    if (_tds->dimension() == 1) {Iterator_base::increment();}
    else {
      if (index == 0) {
	Iterator_base::decrement();
	index = 2;
      }
      else {index -= 1;}
      return;
    }
  }

bool associated_edge()
  {
    if (_tds->dimension() == 1) {return true;}
    int max = maximum(pos);
    if (index == cw(max))  {return true; }
    if (index == ccw(max)) { return false;}
    // index = maximun(pos)
    if ( pos->vertex(max) != _tds->infinite_vertex()){
      return ( _tds->geom_traits().compare_y(pos->vertex(cw(max))->point(),
					     pos->vertex(ccw(max))->point())
               == CGAL_LARGER);
    }
    else{ //pos->vertex(max) == infinite_vertex(), convex hull edge
      return ( _tds->geom_traits().compare_y(pos->vertex(cw(max))->point(),
					     pos->vertex(ccw(max))->point())
               == CGAL_LARGER ||
	       (_tds->geom_traits().compare_y(pos->vertex(cw(max))->point(),
					     pos->vertex(ccw(max))->point())
		== CGAL_EQUAL  // convex hull horizntal edges
		&&
		_tds->geom_traits().compare_x(pos->vertex(cw(max))->point(),
					     pos->vertex(ccw(max))->point())
               == CGAL_SMALLER )
	       );
    }
  }
    
public:   
Edge_iterator&
    operator++()
    {
        if (pos==NULL){
            return *this;            // cannot advance past-the-end iterator
        }

	if(_tds->dimension()==1){
	  Iterator_base::increment();
	  if ( pos == _tds->infinite_face()) { pos = NULL ; return *this;}
	}
	else{
         do{
            increment();
            if ( pos == _tds->infinite_face() && index == 0){
	      pos = NULL;   // complete tour
	      return *this;
            }
	 }while ( ! associated_edge() );
	}
        return *this;
    }
    



 Edge_iterator&
    operator--()
    {
      switch(_tds->dimension()){
      case 0:
	return *this;
      case 1:
	if (pos == _tds->infinite_face()) { pos=NULL; return *this;}
	if (pos == NULL) { pos= _tds->infinite_face;} 
	decrement(); 
	return *this;
      case 2:
	if (pos == _tds->infinite_face() && index == 0){
	  pos == NULL; //first edge, cannot decrement
	  return this;
	}
	if (pos == NULL){ //past the end, can decrement;
	  *this = Edge_iterator(_tds);
	}
	do{
	  decrement();
	  if ( pos == _tds->infinite_face() && index == 0){
	    pos = NULL;   // complete tour
	    return *this;
	  }
	}while ( ! associated_edge());
	return *this;
      }
    }
    
    
    Edge_iterator
    operator++(int)
    {
        Edge_iterator tmp(*this);
        ++(*this);
        return tmp;
    }
    
    
    Edge_iterator
    operator--(int)
    {
        Edge_iterator tmp(*this);
        --(*this);
        return tmp;
    }
    
    Edge
    operator*() const
    {
      return make_pair(pos, index);
    }
};

#endif CGAL_TRIANGULATION_DS_ITERATORS_2_H

