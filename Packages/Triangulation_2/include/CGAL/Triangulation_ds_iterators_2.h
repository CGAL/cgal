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
//   typedef Gt Geom_traits;
//   
//   typedef CGAL_Triangulation_ds_vertex_2<Vb,Fb> Vertex;
//   typedef CGAL_Triangulation_ds_face_2<Vb,Fb> Face;
//   typedef pair<Face*, int>  Edge;
// 
//   typedef CGAL_Triangulation_default_data_structure_2<Gt,Vb,Fb> Tds;

  typedef typename Tds::Geom_traits Geom_traits;
  typedef typename Tds::Vertex Vertex;
  typedef typename Tds::Face  Face;
  typedef typename Tds::Edge Edge;

  CGAL_Triangulation_ds_iterator_base_2()
     : pos(NULL),  _tds(NULL)
        {}

  CGAL_Triangulation_ds_iterator_base_2(Tds* tds)
     : pos(NULL),  _tds(tds)
  {
    if(_tds->number_of_vertices() < 2) {
      return;
    }
    pos = _tds->infinite_face();
  }

 CGAL_Triangulation_ds_iterator_base_2(Tds* tds, int i)
     : pos(NULL),  _tds(tds)
 {}

  CGAL_Triangulation_ds_iterator_base_2(Tds* tds, Face* f)
    : pos(f), _tds(tds)
  {}
   
protected:
        Tds*  _tds;
        Face* pos;
        
        
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



//template < class Gt , class Vb, class Fb> 
template<class Tds>
class CGAL_Triangulation_ds_face_iterator_2
//  : public CGAL_Triangulation_ds_iterator_base_2<Gt,Vb,Fb>,
//    public bidirectional_iterator<CGAL_Triangulation_ds_face_2<Vb,Fb>, ptrdiff_t>
  : public CGAL_Triangulation_ds_iterator_base_2<Tds>,
    public bidirectional_iterator<typename Tds::Face, ptrdiff_t>
{
public:
//   typedef Gt Geom_traits;
//   
//   typedef CGAL_Triangulation_ds_vertex_2<Vb,Fb> Vertex;
//   typedef CGAL_Triangulation_ds_face_2<Vb,Fb> Face;
//   typedef pair<Face*, int>  Edge;
// 
//   typedef CGAL_Triangulation_default_data_structure_2<Gt,Vb,Fb> Tds;
//   typedef CGAL_Triangulation_ds_iterator_base_2<Gt,Vb,Fb> Iterator_base;
// 
//   typedef CGAL_Triangulation_ds_face_iterator_2<Gt,Vb,Fb> Face_iterator;
//   //  typedef CGAL_Triangulation_ds_vertex_iterator_2<Gt,Vb,Fb> Vertex_iterator;
//   // typedef CGAL_Triangulation_ds_edge_iterator_2<Gt,Vb,Fb> Edge_iterator;
//   

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
        {
	  if (tds->number_of_vertices()<2){
	    pos = NULL;                   // there is no faces
	    return;
	  }
	  Face* start = pos;
            while (_tds->is_infinite(pos)){
                increment();
                if(pos == start){
                    pos = NULL;                  // there is no finite triangle
                    return;
                }
            }
	}

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
        operator!=(const Face_iterator& fi)
        {
            return !(*this == fi);
        }

        Face_iterator&
        operator++()
        {
            if ( pos == NULL ){
                return *this;    //  past-the-end iterator cannot advance
            }
            do{
                increment();
                if ( pos == (_tds->infinite_face()) ){
                    pos = NULL;  // complete tour
                    return *this;
                }
            }while (_tds->is_infinite(pos));
	    return *this;           // next finite triangle found
        }

        Face_iterator&
        operator--()
        {
            if ( pos == NULL ) {//  past the end iterator can decrease
                *this = Face_iterator(_tds); // first finite triangle
            }              //next loop will go to last finite triangle
        
            do{
                decrement();
                if ( pos == _tds->infinite_face()){
                    pos = NULL;  // complete tour
		    return *this;
                }
           }while (_tds->is_infinite(pos));
	   return *this;           // next finite triangle found
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


//template < class Gt , class Vb, class Fb> 
template < class Tds>
class CGAL_Triangulation_ds_vertex_iterator_2
//  : public CGAL_Triangulation_ds_iterator_base_2<Gt,Vb,Fb>,
//    public bidirectional_iterator<CGAL_Triangulation_ds_vertex_2<Vb,Fb>, ptrdiff_t>
: public CGAL_Triangulation_ds_iterator_base_2<Tds>,
  public bidirectional_iterator<typename Tds::Vertex, ptrdiff_t>
{
public:
//   typedef Gt Geom_traits;
//   
//   typedef CGAL_Triangulation_ds_vertex_2<Vb,Fb> Vertex;
//   typedef CGAL_Triangulation_ds_face_2<Vb,Fb> Face;
//   typedef pair<Face*, int>  Edge;
// 
//   typedef CGAL_Triangulation_default_data_structure_2<Gt,Vb,Fb> Tds;
//   typedef CGAL_Triangulation_ds_iterator_base_2<Gt,Vb,Fb> Iterator_base;
// 
//   typedef CGAL_Triangulation_ds_face_iterator_2<Gt,Vb,Fb> Face_iterator;
//   typedef CGAL_Triangulation_ds_vertex_iterator_2<Gt,Vb,Fb> Vertex_iterator;
// 
  typedef typename Tds::Geom_traits Geom_traits;
  typedef typename Tds::Vertex Vertex;
  typedef typename Tds::Face  Face;
  typedef typename Tds::Edge Edge;

  typedef CGAL_Triangulation_ds_iterator_base_2<Tds> Iterator_base;
  typedef CGAL_Triangulation_ds_vertex_iterator_2<Tds> Vertex_iterator;


     CGAL_Triangulation_ds_vertex_iterator_2()
            : Iterator_base()
        {}
    
    
    CGAL_Triangulation_ds_vertex_iterator_2(Tds * tds)
        :  Iterator_base(tds)
    {
        switch( _tds->number_of_vertices() ){
        case 0: // past-the-end
            pos = NULL;
            return;
        case 1:
            pos = (Face*)1 ; // different from any pointer;
            return;         // "points" to the only vertex of the triangulation
        default:
            {
                pos = _tds->infinite_face();
                while ( associated_vertex() ==NULL){
                    increment();
                }
                return;
            }
        }
    }
    
    CGAL_Triangulation_ds_vertex_iterator_2(Tds* tds, int i)
            : Iterator_base(tds,i)
    {}

    Vertex_iterator&
    operator++()
    {
        if (pos==NULL){
	  return *this;            // cannot advance past-the-end iterator
        }
        if (_tds->number_of_vertices()==1){
            pos = NULL; // past-the-end
	    return *this;
        }
        do{
            increment();
            if ( pos == _tds->infinite_face()){
                pos = NULL;   // complete tour
		return *this;
            }
        }while ( associated_vertex() ==NULL);
        return *this;
    }
    
    
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
            if (pos==NULL){
                *this = Vertex_iterator(_tds);
                --*this;
		return *this;            // can decrease past-the-end iterator
            }
            do{
                decrement();
                if ( pos == _tds->infinite_face()){
                    pos = NULL;   // complete tour
		    return *this;
                }
            }while ( associated_vertex()  ==NULL);
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
        return ( (pos == fi.pos ) && (_tds ==fi._tds) );
    }
    
    bool operator!=(const Vertex_iterator& fi) const
    {
        return !(*this == fi);
    }
    
  inline Vertex& operator*() const
  {
        return *(associated_vertex());
  }
    
  inline Vertex*  operator->() const
  {
        return associated_vertex();
  }

    
    Vertex*
    associated_vertex() const
    {
        if(pos == NULL) {
            return NULL;
        }
        switch(_tds->number_of_vertices() ){
        case 0:
            return NULL;
        case 1:
            return _tds->finite_vertex();
        case 2:
             return pos->vertex(cw(maximum(pos)));
        default:
            {
                int i = cw(maximum(pos));            // candidate associate vertex
                if(_tds->geom_traits().compare_y(pos->vertex(i)->point(),
                                                   pos->vertex(cw(i))->point())
                   ==CGAL_LARGER){
                    //   vcw(i) < vi
                    return pos->vertex(i);
                }
                if ( _tds->is_infinite(pos)){
                    Face* f=pos->neighbor(cw(i));
                                                 // next edge on the CH in cw order
                    int   j=f->index(_tds->infinite_vertex() );
                    CGAL_Comparison_result comp =
                        _tds->geom_traits().compare_y(pos->vertex(i)->point(),
                                                        f->vertex(cw(j))->point());
                    if (comp == CGAL_SMALLER){         //  vi smaller vertex
                        return pos->vertex(i);
                    }
                    if(comp ==CGAL_EQUAL){            // horizontal part of CH
                        comp = _tds->geom_traits().compare_x(
                                                        pos->vertex(i)->point(),
                                                        f->vertex(cw(j))->point());
                        if(comp == CGAL_LARGER){        // lower hull
                            return pos->vertex(i);
                        }
                        if(pos->vertex(cw(i)) == f->vertex(cw(j))){ // 1 dim. horiz.
                            return pos->vertex(i);
                        }
                    }
                }
            }
        }
        return NULL;
    }
    
};


//template < class Gt , class Vb, class Fb> 
template <class Tds>
class CGAL_Triangulation_ds_edge_iterator_2
//  : public CGAL_Triangulation_ds_iterator_base_2<Gt,Vb,Fb>,
//    public bidirectional_iterator<pair<CGAL_Triangulation_ds_face_2<Vb,Fb>,int>, ptrdiff_t>
 : public CGAL_Triangulation_ds_iterator_base_2<Tds>,
   public bidirectional_iterator<typename Tds::Edge, ptrdiff_t>
{
public:
// //   typedef Gt Geom_traits;
// //   
// //   typedef CGAL_Triangulation_ds_vertex_2<Vb,Fb> Vertex;
// //   typedef CGAL_Triangulation_ds_face_2<Vb,Fb> Face;
// //   typedef pair<Face*, int>  Edge;
// // 
// //   typedef CGAL_Triangulation_default_data_structure_2<Gt,Vb,Fb> Tds;
//   typedef CGAL_Triangulation_ds_iterator_base_2<Gt,Vb,Fb> Iterator_base;
// 
// //   typedef CGAL_Triangulation_ds_face_iterator_2<Gt,Vb,Fb> Face_iterator;
// //   typedef CGAL_Triangulation_ds_vertex_iterator_2<Gt,Vb,Fb> Vertex_iterator;
//   typedef CGAL_Triangulation_ds_edge_iterator_2<Gt,Vb,Fb> Edge_iterator;
// 
  typedef typename Tds::Geom_traits Geom_traits;
  typedef typename Tds::Vertex Vertex;
  typedef typename Tds::Face  Face;
  typedef typename Tds::Edge Edge;

  typedef CGAL_Triangulation_ds_iterator_base_2<Tds> Iterator_base;
  typedef CGAL_Triangulation_ds_edge_iterator_2<Tds> Edge_iterator;

     CGAL_Triangulation_ds_edge_iterator_2()
            : Iterator_base(),status(0)
        {}
    
    
    CGAL_Triangulation_ds_edge_iterator_2(Tds * tds)
        :  Iterator_base(tds)
    {
      if (tds->number_of_vertices()<2){
	pos = NULL;                   // there is no finite edge
	return;
      }
      pos = tds->infinite_face();
      Face* start = pos;
      while ( ! compute_status(true) ){
	increment();
	if(pos == start){
	  pos=NULL;                   // there is no finite triangle
	  return;
	}
      }
    } 
    
    CGAL_Triangulation_ds_edge_iterator_2(Tds* tds, int i)
            : Iterator_base(tds,i),status(0)
    {}

  bool
  operator==(const Edge_iterator& fi) const
  {
    if ((pos==fi.pos)&&(_tds==fi._tds)){
      if (pos==NULL) {
	return true;
      }else{
	return (status==fi.status);
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
    
Edge_iterator&
    operator++()
    {
        if (pos==NULL){
            return *this;            // cannot advance past-the-end iterator
        }
        if ( (status >=3) && (status <=5) ) {
            // current status is 3 + ccw(i), next status is 6+i
            status = 6 + cw( status -3);
        } else do{
            increment();
            if ( pos == _tds->infinite_face()){
                pos = NULL;   // tour complete
                return *this;
            }
        }while (!compute_status(true));
        return *this;
    }
    



 Edge_iterator&
    operator--()
    {
        if (pos==NULL){
            *this = Edge_iterator(_tds);
            --*this;
            return *this;            // can decrease past-the=end iteartor
        }
        if ( (status >=3) && (status <=5) ) {
            // current status is 6+i, next status is 3 + ccw(i)
            status = 3 + ccw( status -6);
        } else do{
            decrement();
            if ( pos == triangulation->infinite_face()){
                pos = NULL;   // tour complete
                return *this;
            }
        }while (compute_status(false));
        return *this;
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
        int j = (status<3) ? status : (status<6) ? status-3 : status-6;
        return make_pair(pos, j);
    }
    
    
private:
        
    int status;

    bool
    compute_status(bool forward)
    {
        if (pos == NULL){
            return false;
        }
        int i = maximum(pos);            // higher point
        if ( ! _tds->is_infinite(pos)){
            if(_tds->geom_traits().compare_y(pos->vertex(cw(i))->point(),
                                               pos->vertex(ccw(i))->point())
               == CGAL_LARGER){
                //  vccw(i) < vcw(i)
                // right edges = i cw(i)  cw(i) ccw(i)
                status = (forward) ? 3+ccw(i) : 6+i;
            } else {
                // right edge = i cw(i)
                status = ccw(i);
            }
            return true;
        }
        CGAL_Comparison_result comp = _tds->geom_traits().compare_y(
                                                      pos->vertex(cw(i))->point(),
                                                      pos->vertex(ccw(i))->point());
        if( comp==CGAL_SMALLER){
            return false;
        }
        if( comp==CGAL_LARGER){
            //  vccw(i) < vcw(i)
            // right edge =  cw(i) ccw(i)
            status = i;
            return true;
        }
        // comp = CGAL_EQUAL
        if ( _tds->geom_traits().compare_x(pos->vertex(cw(i))->point(),
                                             pos->vertex(ccw(i))->point())
             == CGAL_SMALLER){
            //  vccw(i) < vcw(i)
            // upper horizontal edge =  cw(i) ccw(i)
            status = i;
            return true;
        }
    
        return false;
    }
};

#endif CGAL_TRIANGULATION_DS_ITERATORS_2_H

