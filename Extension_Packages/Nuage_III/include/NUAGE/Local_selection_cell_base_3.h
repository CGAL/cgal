// ======================================================================
//
// Copyright (c) 2002 INRIA
//
// ----------------------------------------------------------------------
//
// file          : include/NUAGE/Local_selection_cell_base_3.h
// package       : Local_selections_3(1.0)
// revision_date : $Date$
// author(s)     : Tran Kai Frank DA <Frank.Da@sophia.inria.fr>
//
// ======================================================================

#ifndef LOCAL_SELECTION_CELL_BASE_3_H
#define LOCAL_SELECTION_CELL_BASE_3_H

#include <CGAL/Triangulation_cell_base_3.h>

template < class Gt >
class Local_selection_cell_base_3 : public CGAL::Triangulation_cell_base_3<Gt>
{
private:
               
  typedef double coord_type;        
  typedef typename CGAL::Simple_cartesian<coord_type>::Point_3 D_Point;

  //-------------------- DATA MEMBERS ---------------------------------

  coord_type* _smallest_radius_facet_tab;
  bool* _selected_facet;
#ifndef NOLAZY
  D_Point* _circumcenter;
  coord_type* _squared_radius;
#endif //NOLAZY 

  //-------------------- CONSTRUCTORS ----------------------------------

public:
  
  Local_selection_cell_base_3() 
    : CGAL::Triangulation_cell_base_3<Gt>(),
      _smallest_radius_facet_tab(NULL), _selected_facet(NULL)
#ifndef NOLAZY
      , _circumcenter(NULL), _squared_radius(NULL)
#endif //NOLAZY
    {}
  
  Local_selection_cell_base_3(void* v0, void* v1, void* v2, void* v3)
    : CGAL::Triangulation_cell_base_3<Gt>( v0, v1, v2, v3),
      _smallest_radius_facet_tab(NULL), _selected_facet(NULL)
#ifndef NOLAZY
      , _circumcenter(NULL), _squared_radius(NULL)
#endif //NOLAZY
    {}
  
  Local_selection_cell_base_3(void* v0, void* v1, void* v2, void* v3,
			      void* n0, void* n1, void* n2, void* n3)
    : CGAL::Triangulation_cell_base_3<Gt>(v0,  v1,  v2, v3,
					  n0,  n1,  n2, n3),
      _smallest_radius_facet_tab(NULL), _selected_facet(NULL)
#ifndef NOLAZY
      , _circumcenter(NULL), _squared_radius(NULL)
#endif //NOLAZY
    {}

  //-------------------- DESTRUCTOR -----------------------------------

  inline ~Local_selection_cell_base_3()
    {
      if (_smallest_radius_facet_tab != NULL)
	 delete[] _smallest_radius_facet_tab;
      if (_selected_facet != NULL)
	delete[] _selected_facet;
#ifndef NOLAZY
      if (_circumcenter != NULL)
	delete _circumcenter;
      if (_squared_radius != NULL)
	delete _squared_radius;
#endif //NOLAZY
    }

  //--------------------  MEMBER FUNCTIONS ----------------------------

public:
  
  inline void clear()
    {
      if (_smallest_radius_facet_tab != NULL)
	delete[] _smallest_radius_facet_tab;
      if (_selected_facet != NULL)
	delete[] _selected_facet;
      _smallest_radius_facet_tab = NULL;
      _selected_facet = NULL;
#ifndef NOLAZY
      if (_circumcenter != NULL)
	delete _circumcenter;
      _circumcenter = NULL;
      if (_squared_radius != NULL)
	delete _squared_radius;
      _squared_radius = NULL;
#endif //NOLAZY
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
  
  inline void select_facet(const int& i)
    {
      if (_selected_facet==NULL)
	{
	  _selected_facet = new bool[4];
	   for(int i = 0; i < 4; i++)
	     _selected_facet[i] = false;
	}
      _selected_facet[i] = true;
    }
  
  inline void unselect_facet(const int& i)
    { 
      _selected_facet[i] = false;
    }

  inline bool is_selected_facet(const int& i)
    {
      if (_selected_facet==NULL) return false;
      return _selected_facet[i];
    }
 
#ifndef NOLAZY
 
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

#endif //NOLAZY

};

#endif LOCAL_SELECTION_CELL_BASE_3_H
