#ifndef CGAL_AFSR_CELL_BASE_3_H
#define CGAL_AFSR_CELL_BASE_3_H

namespace CGAL {

template < class CellBase >
class AFSR_cell_base_3 : public CellBase
{

public:
  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename CellBase::template Rebind_TDS<TDS2>::Other  Cb2;
    typedef AFSR_cell_base_3<Cb2>                    Other;
  };

  typedef typename CellBase::Vertex_handle Vertex_handle;
  typedef typename CellBase::Cell_handle Cell_handle;

private:
     
#ifdef AFSR_FACET_NUMBER
  int _facet_number[4];
#endif
  typedef double coord_type; 
#ifdef AFSR_LAZY       
  typedef typename CGAL::Simple_cartesian<coord_type>::Point_3 D_Point;
#endif 
  //-------------------- DATA MEMBERS ---------------------------------

  coord_type* _smallest_radius_facet_tab;
  unsigned char selected_facet;
#ifdef AFSR_LAZY
  D_Point* _circumcenter;
  coord_type* _squared_radius;
#endif 

  //-------------------- CONSTRUCTORS ----------------------------------

public:
  
  AFSR_cell_base_3() 
    : CellBase(),
      _smallest_radius_facet_tab(NULL), selected_facet(0)
#ifdef AFSR_LAZY
      , _circumcenter(NULL), _squared_radius(NULL)
#endif
    {
#ifdef AFSR_FACET_NUMBER
      for(int i = 0; i < 4; i++){
	_facet_number[i] = -1;
      }
#endif
    }
  
  AFSR_cell_base_3(Vertex_handle v0, Vertex_handle v1, Vertex_handle v2, Vertex_handle v3)
    : CellBase( v0, v1, v2, v3),
      _smallest_radius_facet_tab(NULL), selected_facet(0)
#ifdef AFSR_LAZY
      , _circumcenter(NULL), _squared_radius(NULL)
#endif
    {
#ifdef FACET_NUMBER
      for(int i = 0; i < 4; i++){
	_facet_number[i] = -1;
      }
#endif
    }
  
  AFSR_cell_base_3(Vertex_handle v0, Vertex_handle v1, Vertex_handle v2, Vertex_handle v3,
			      Cell_handle n0, Cell_handle n1, Cell_handle n2, Cell_handle n3)
    : CellBase(v0,  v1,  v2, v3,
	       n0,  n1,  n2, n3),
      _smallest_radius_facet_tab(NULL), selected_facet(0)
#ifdef AFSR_LAZY
      , _circumcenter(NULL), _squared_radius(NULL)
#endif
    {
#ifdef AFSR_FACET_NUMBER
      for(int i = 0; i < 4; i++){
	_facet_number[i] = -1;
      }
#endif
    }

  //-------------------- DESTRUCTOR -----------------------------------

  inline ~AFSR_cell_base_3()
    {
      if (_smallest_radius_facet_tab != NULL)
	 delete[] _smallest_radius_facet_tab;
#ifdef AFSR_LAZY
      if (_circumcenter != NULL)
	delete _circumcenter;
      if (_squared_radius != NULL)
	delete _squared_radius;
#endif
    }

  //--------------------  MEMBER FUNCTIONS ----------------------------

public:
  
  inline void clear()
    {
      if (_smallest_radius_facet_tab != NULL)
	delete[] _smallest_radius_facet_tab;
      _smallest_radius_facet_tab = NULL;
      selected_facet = 0;
#ifdef AFSR_LAZY
      if (_circumcenter != NULL)
	delete _circumcenter;
      _circumcenter = NULL;
      if (_squared_radius != NULL)
	delete _squared_radius;
      _squared_radius = NULL;
#endif
    }

  //-------------------------------------------------------------------

  inline coord_type get_smallest_radius(const int& i)
    { 
      if (_smallest_radius_facet_tab == NULL)
	return -1;
      return _smallest_radius_facet_tab[i];
    }
  
  inline void set_smallest_radius(const int& i, const coord_type& c)
    {
      if (_smallest_radius_facet_tab == NULL)
	{
	  _smallest_radius_facet_tab = new coord_type[4];
	  for(int i = 0; i < 4; i++)
	    _smallest_radius_facet_tab[i] = -1;
	}
      _smallest_radius_facet_tab[i] = c;
    }

  // pour un controle de l'allocation memoire... utile???
  inline bool alloc_smallest_radius_tab(coord_type* ptr)
    {
      if (_smallest_radius_facet_tab==NULL)
	{
	  _smallest_radius_facet_tab = ptr;
	  return true;
	}
      return false;
    }

  
  //-------------------------------------------------------------------

#ifdef FACET_NUMBER
  void set_facet_number(int i, int n){}
  {
    _facet_number[i] = n;
  }

  int facet_number(int i)
  {
    return _facet_number[i];
  }
#else
  void set_facet_number(int, int){}
  int facet_number(int){return 0;}
#endif


  //-------------------------------------------------------------------
  
  inline void select_facet(const int& i)
    {
      selected_facet |= (1 << i);
    }
  
  inline void unselect_facet(const int& i)
    { 
      selected_facet &= (15 - (1 << i));
    }

  inline bool is_selected_facet(const int& i)
    {
      return (selected_facet & (1 << i)) != 0;
    }

  inline bool has_facet_on_surface(const int& i)
    {
      return (selected_facet & (1 << i)) != 0;
    }
 
#ifdef AFSR_LAZY
 
  //-------------------------------------------------------------------

  inline D_Point* get_lazy_circumcenter()
    {
      return _circumcenter;
    }

  inline void set_lazy_circumcenter(const D_Point& p)
    {
      _circumcenter = new D_Point(p);
    }
  
  //-------------------------------------------------------------------

  inline coord_type* get_lazy_squared_radius()
    { 
      return _squared_radius;
    }

  inline void set_lazy_squared_radius(const coord_type& sr)
    { 
      _squared_radius = new coord_type(sr);
    }

#endif //AFSR_LAZY

};

} // namespace CGAL

#endif // CGAL_AFSR_CELL_BASE_3_H
