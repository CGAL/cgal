#ifndef CGAL_TRIANGULATION_CIRCULATORS_2_H
#define CGAL_TRIANGULATION_CIRCULATORS_2_H

#include <pair.h>
#include <iterator.h>
#include <CGAL/circulator.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/Triangulation_ds_circulators_2.h>

template <class Gt, class Tds >
class CGAL_Triangulation_2;

template <class Gt,  class Tds >
class CGAL_Triangulation_face_2;

template <class Gt,  class Tds >
class CGAL_Triangulation_vertex_2;

template < class Gt, class Tds >
class CGAL_Triangulation_face_circulator_2;

template < class Gt, class Tds >
class CGAL_Triangulation_vertex_circulator_2;

template < class Gt, class Tds >
class CGAL_Triangulation_edge_circulator_2;

template < class Gt, class Tds>
class CGAL_Triangulation_face_circulator_2
    : public CGAL_Bidirectional_circulator_base<CGAL_Triangulation_face_2<Gt,Tds>,ptrdiff_t,size_t>
{
public:
  typedef typename Tds::Vertex Ve;
  typedef typename Tds::Face Fa;
  typedef typename Tds::Face_circulator  Base_face_circulator;

  typedef CGAL_Triangulation_face_2<Gt,Tds> Face;
  typedef CGAL_Triangulation_vertex_2<Gt, Tds> Vertex;
  typedef CGAL_Triangulation_vertex_handle_2<Gt, Tds> Vertex_handle;
  typedef CGAL_Triangulation_face_handle_2<Gt, Tds>  Face_handle;
  typedef pair<Face_handle, int>     Edge;

  typedef CGAL_Triangulation_face_circulator_2<Gt,Tds>      Face_circulator;
  typedef CGAL_Triangulation_edge_circulator_2<Gt,Tds>      Edge_circulator;
  typedef CGAL_Triangulation_vertex_circulator_2<Gt,Tds>    Vertex_circulator;


  static int ccw(int i)
  {
        return (i+1) % 3;
  }

  static int cw(int i)
  {
    return (i+2) % 3;
  }

  CGAL_Triangulation_face_circulator_2()
    : _bfc()
  {}

  CGAL_Triangulation_face_circulator_2(const Vertex_handle& v)
    : _bfc(&(*v))
  {}

  CGAL_Triangulation_face_circulator_2(const Vertex_handle& v,
				       const Face_handle& f)
    : _bfc( &(*v), &(*f))
  {}

 
   
   CGAL_Triangulation_face_circulator_2(const Face_circulator &fc)
    : _bfc(fc._bfc)
  {}
   
  CGAL_Triangulation_face_circulator_2 &operator=(const Face_circulator &fc)
  {
    _bfc = fc._bfc;
    return *this;
  }
  
  Face_circulator& operator++()
  {
    _bfc++;
      return *this;
  }
        
  Face_circulator operator++(int)
  {
    Face_circulator tmp(*this);
    ++(*this);
    return tmp;
  }
        
  Face_circulator& operator--()
  {
    _bfc--;
    return *this;
  }
        
  Face_circulator operator--(int)
  {
    Face_circulator tmp(*this);
    --(*this);
    return tmp;
  }
  
  bool operator==(const Face_circulator &fc) const
  {
            return (_bfc == fc._bfc) ;
  }
        
        
  bool operator!=(const Face_circulator &fc) const
  {
            return ! (*this == fc);
  }
        
  bool  operator==(CGAL_NULL_TYPE n) const
  {
    CGAL_triangulation_assertion( n == NULL);
    return (_bfc == NULL);
  }
        
  bool
  operator!=(CGAL_NULL_TYPE n) const
  {
    CGAL_triangulation_assertion( n == NULL);
    return ! (*this == NULL);
  }
        
  inline Face& operator*() const
  {
    return (Face&  )(*_bfc);
  }

  inline Face* operator->() const
  {
    return   (Face*)( & (*_bfc));
  }
  
  
private:
  Base_face_circulator  _bfc;
};


template < class Gt, class Tds>
class CGAL_Triangulation_vertex_circulator_2
    : public CGAL_Bidirectional_circulator_base<CGAL_Triangulation_vertex_2<Gt,Tds>,ptrdiff_t,size_t>
{
public:
  typedef typename Tds::Vertex Ve;
  typedef typename Tds::Face Fa;
  typedef typename Tds::Vertex_circulator  Base_vertex_circulator;

  typedef CGAL_Triangulation_face_2<Gt,Tds> Face;
  typedef CGAL_Triangulation_vertex_2<Gt,Tds> Vertex;
  typedef CGAL_Triangulation_vertex_handle_2<Gt,Tds> Vertex_handle;
  typedef CGAL_Triangulation_face_handle_2<Gt,Tds> Face_handle;
  typedef pair<Face_handle, int>     Edge;

  typedef CGAL_Triangulation_face_circulator_2<Gt,Tds>      Face_circulator;
  typedef CGAL_Triangulation_edge_circulator_2<Gt,Tds>      Edge_circulator;
  typedef CGAL_Triangulation_vertex_circulator_2<Gt,Tds>    Vertex_circulator;


  static int ccw(int i)
  {
        return (i+1) % 3;
  }

  static int cw(int i)
  {
    return (i+2) % 3;
  }

  CGAL_Triangulation_vertex_circulator_2()
    : _bvc()
  {}

  CGAL_Triangulation_vertex_circulator_2(const Vertex_handle& v)
    : _bvc(&(*v))
  {}

  CGAL_Triangulation_vertex_circulator_2(const Vertex_handle& v,
				       const Face_handle& f)
    : _bvc( &(*v), &(*f))
  {}
   
   CGAL_Triangulation_vertex_circulator_2(const Vertex_circulator &vc)
    : _bvc(vc._bvc)
  {}
   
  CGAL_Triangulation_vertex_circulator_2 &operator=(const Vertex_circulator &vc)
  {
    _bvc = vc._bvc;
    return *this;
  }
  
  Vertex_circulator& operator++()
  {
    _bvc++;
    return *this;
  }
        
  Vertex_circulator operator++(int)
  {
    Vertex_circulator tmp(*this);
    ++(*this);
    return tmp;
  }
        
  Vertex_circulator& operator--()
  {
    _bvc--;
    return *this;
  }
        
  Vertex_circulator operator--(int)
  {
    Vertex_circulator tmp(*this);
    --(*this);
    return tmp;
  }
  
  bool operator==(const Vertex_circulator &vc) const
  {
            return (_bvc == vc._bvc) ;
  }
        
        
  bool operator!=(const Vertex_circulator &vc) const
  {
            return ! (*this == vc);
  }
        
  bool  operator==(CGAL_NULL_TYPE n) const
  {
    CGAL_triangulation_assertion( n == NULL);
    return (_bvc == NULL);
  }
        
  bool
  operator!=(CGAL_NULL_TYPE n) const
  {
    CGAL_triangulation_assertion( n == NULL);
    return ! (*this == NULL);
  }
        
  inline Vertex& operator*() const
  {
    return (Vertex &)(*_bvc);
  }

  inline Vertex* operator->() const
  {
    return   (Vertex*)( & (*_bvc));
  }
  
  
private:
  Base_vertex_circulator  _bvc;
};


template < class Gt, class Tds>
class CGAL_Triangulation_edge_circulator_2
  : public CGAL_Bidirectional_circulator_base< typename CGAL_Triangulation_2<Gt, Tds>::Edge,ptrdiff_t,size_t>
{
public:
  typedef Tds Triangulation_data_structure;

  typedef typename Tds::Vertex Ve;
  typedef typename Tds::Face Fa;
  typedef typename Tds::Edge_circulator  Base_edge_circulator;

  typedef CGAL_Triangulation_face_2<Gt,Tds> Face;
  typedef CGAL_Triangulation_vertex_2<Gt,Tds> Vertex;
  typedef CGAL_Triangulation_vertex_handle_2<Gt,Tds> Vertex_handle;
  typedef CGAL_Triangulation_face_handle_2<Gt,Tds> Face_handle;
  typedef pair<Face_handle, int>     Edge;

  typedef CGAL_Triangulation_face_circulator_2<Gt,Tds>      Face_circulator;
  typedef CGAL_Triangulation_edge_circulator_2<Gt,Tds>      Edge_circulator;
  typedef CGAL_Triangulation_vertex_circulator_2<Gt,Tds>    Vertex_circulator;


  static int ccw(int i)
  {
        return (i+1) % 3;
  }

  static int cw(int i)
  {
    return (i+2) % 3;
  }

  CGAL_Triangulation_edge_circulator_2()
    : _bec()
  {}

  CGAL_Triangulation_edge_circulator_2(const Vertex_handle& v)
    : _bec(&(*v))
  {}

  CGAL_Triangulation_edge_circulator_2(const Vertex_handle& v,
				       const Face_handle& f)
    : _bec( &(*v), &(*f))
  {}
   
   CGAL_Triangulation_edge_circulator_2(const Vertex_circulator &ec)
    : _bec(ec._bec)
  {}
   
  CGAL_Triangulation_edge_circulator_2 &operator=(const Edge_circulator &ec)
  {
    _bec = ec._bec;
    return *this;
  }
  
  Edge_circulator& operator++()
  {
    _bec++;
    return *this;
  }
        
  Edge_circulator operator++(int)
  {
    Vertex_circulator tmp(*this);
    ++(*this);
    return tmp;
  }
        
  Edge_circulator& operator--()
  {
    _bec--;
    return *this;
  }
        
  Edge_circulator operator--(int)
  {
    Vertex_circulator tmp(*this);
    --(*this);
    return tmp;
  }
  
  bool operator==(const Edge_circulator &ec) const
  {
            return (_bec == ec._bec) ;
  }
        
        
  bool operator!=(const Edge_circulator &ec) const
  {
            return ! (*this == ec);
  }
        
  bool  operator==(CGAL_NULL_TYPE n) const
  {
    CGAL_triangulation_assertion( n == NULL);
    return (_bec == NULL);
  }
        
  bool
  operator!=(CGAL_NULL_TYPE n) const
  {
    CGAL_triangulation_assertion( n == NULL);
    return ! (*this == NULL);
  }
        
  inline Edge& operator*() const
  {
    Face_handle fh = (Face *) ((*_bfc).first);
    return make_pair( fh  , (*_bfc).second );
  }

    
private:
  Base_edge_circulator  _bec;
};


#endif CGAL_TRIANGULATION_CIRCULATORS_2_H
