template < class Tds >
class CGAL_Triangulation_face_2;

template < class Tds >
class CGAL_Triangulation_vertex_2;

template < class Tds>
class CGAL_Triangulation_face_iterator_2;

template < class Tds>
class CGAL_Triangulation_vertex_iterator_2;

template < class Tds>
class CGAL_Triangulation_face_circulator_2;

template < class Tds>
class CGAL_Triangulation_vertex_circulator_2;


template < class Tds>
class CGAL_Triangulation_face_handle_2
  :public CGAL_Pointer<CGAL_Triangulation_face_2<Tds> > 
{
public:
  typedef CGAL_Pointer<CGAL_Triangulation_face_2<Tds> > Pointer;
  typedef CGAL_Triangulation_face_2<Tds> Face;
  typedef CGAL_Triangulation_face_handle_2<Tds> Face_handle;
  
  typedef CGAL_Triangulation_face_iterator_2<Tds>      Face_iterator;
  typedef CGAL_Triangulation_face_circulator_2<Tds>    Face_circulator;
  
  inline 
  CGAL_Triangulation_face_handle_2()
    : Pointer(NULL)
  {}

  inline  
  CGAL_Triangulation_face_handle_2(const Face* p)
    : Pointer((Face*)p)
  {}

  inline Face_handle& operator=(const Face*& p)
    {
        ptr() = p ;
        return *this;
    }
    
    inline Face_handle& operator=(const Face_handle& p)
    {
        ptr() = p.ptr();
        return *this;
    }
  
   inline  
    CGAL_Triangulation_face_handle_2(const Face_iterator& fit)
        : Pointer(&(*fit))
    {}
  

  inline  
   CGAL_Triangulation_face_handle_2(const Face_circulator& fit)
        : Pointer(&(*fit))
    {}
};

template < class Tds>
class CGAL_Triangulation_vertex_handle_2
  :public CGAL_Pointer<CGAL_Triangulation_vertex_2<Tds> > 
{
public:
  typedef CGAL_Pointer<CGAL_Triangulation_vertex_2<Tds> > Pointer;
  typedef CGAL_Triangulation_vertex_2<Tds> Vertex;
  typedef CGAL_Triangulation_vertex_handle_2<Tds> Vertex_handle;
  
  typedef CGAL_Triangulation_vertex_iterator_2<Tds>      Vertex_iterator;
  typedef CGAL_Triangulation_vertex_circulator_2<Tds>    Vertex_circulator;
  
  inline 
  CGAL_Triangulation_vertex_handle_2()
    : Pointer(NULL)
  {}

  inline  
  CGAL_Triangulation_vertex_handle_2(const Vertex* p)
        : Pointer((Vertex*)p)
    {}

  inline Vertex_handle& operator=(const Vertex*& p)
    {
        ptr() = p ;
        return *this;
    }
    
    inline Vertex_handle& operator=(const Vertex_handle& p)
    {
        ptr() = p.ptr();
        return *this;
    }
  
   inline  
   CGAL_Triangulation_vertex_handle_2(const Vertex_iterator& fit)
        : Pointer(&(*fit))
    {}
  

  inline  
  CGAL_Triangulation_vertex_handle_2(const Vertex_circulator& fit)
        : Pointer(&(*fit))
    {}
};
