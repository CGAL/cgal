// ======================================================================
//
// Copyright (c) 2002 INRIA
//
// ----------------------------------------------------------------------
//

// file          : include/CGAL/Local_selection_vertex_base_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Tran Kai Frank DA <Frank.Da@sophia.inria.fr>
//
// ======================================================================

#ifndef LOCAL_SELECTION_VERTEX_BASE_3_H
#define LOCAL_SELECTION_VERTEX_BASE_3_H

#include <utility>

#include <list>
#include <map>
#include <set>
#include <CGAL/Triangulation_vertex_base_3.h>

template <class Gt>
class Local_selection_vertex_base_3 : public CGAL::Triangulation_vertex_base_3<Gt>
{
public:

  typedef CGAL::Triangulation_vertex_base_3<Gt> Base; //af
  typedef Base::Point Point; //af
  typedef double coord_type;
  
  typedef CGAL::Triple< void*, int, int > void_Edge;
  typedef std::pair< void_Edge, int > Edge_IFacet;
  typedef std::pair< Edge_IFacet, Edge_IFacet > IO_edge_type;

  //  typedef std::pair< coord_type, coord_type > criteria;
  typedef coord_type criteria;

  typedef std::pair< criteria, IO_edge_type > Radius_edge_type;
  typedef std::pair< Radius_edge_type, int > Border_elt;
  typedef std::pair< void*, Border_elt > Next_border_elt;

private:

  //par convention je remplis d'abord first et si necessaire second...
  typedef std::pair< Next_border_elt*,  Next_border_elt*> Intern_successors_type;

public:

  typedef std::pair< criteria, IO_edge_type* > Radius_ptr_type;
  typedef std::pair< void*, void* > void_Edge_like;
  typedef std::pair< criteria,  void_Edge_like > Incidence_request_elt;
  typedef std::list< Incidence_request_elt > Incidence_request_type;
  typedef typename Incidence_request_type::iterator Incidence_request_iterator;
  
  typedef std::set< void* > Interior_edge_set_type;
  
  //-------------------- DATA MEMBERS ---------------------------------

private:

  int _mark;
  int _visu_index;
  int _post_mark;
  Intern_successors_type* _incident_border;
  Incidence_request_type* _incidence_request;
  Interior_edge_set_type* _interior_edge_set;

  //-------------------- CONSTRUCTORS ---------------------------------

public:

  Local_selection_vertex_base_3()
    : CGAL::Triangulation_vertex_base_3<Gt>(), _mark(-1), 
      _visu_index(-1), _post_mark(-1)
    {
      _incident_border = new Intern_successors_type(new Next_border_elt(),
						    new Next_border_elt());
      _incident_border->first->first = NULL;
      _incident_border->second->first = NULL;
      _incidence_request = new Incidence_request_type();
      _interior_edge_set = new Interior_edge_set_type();
    }
  
  Local_selection_vertex_base_3(const Point & p)
    : CGAL::Triangulation_vertex_base_3<Gt>(p), _mark(-1), 
      _visu_index(-1), _post_mark(-1)
    {
      _incident_border = new Intern_successors_type(new Next_border_elt(),
						    new Next_border_elt());
      _incident_border->first->first = NULL;
      _incident_border->second->first = NULL;
      _incidence_request = new Incidence_request_type();
      _interior_edge_set = new Interior_edge_set_type();
    }
  
  Local_selection_vertex_base_3(const Point & p, void* f)
    : CGAL::Triangulation_vertex_base_3<Gt>(p, f), _mark(-1), 
      _visu_index(-1), _post_mark(-1)
    {
      _incident_border = new Intern_successors_type(new Next_border_elt(),
						    new Next_border_elt());
      _incident_border->first->first = NULL;
      _incident_border->second->first = NULL;
      _incidence_request = new Incidence_request_type();
      _interior_edge_set = new Interior_edge_set_type();
    }

  Local_selection_vertex_base_3(void* f)
    : CGAL::Triangulation_vertex_base_3<Gt>(f), _mark(-1), 
      _visu_index(-1), _post_mark(-1)
    {
      _incident_border = new Intern_successors_type(new Next_border_elt(),
						    new Next_border_elt());
      _incident_border->first->first = NULL;
      _incident_border->second->first = NULL;
      _incidence_request = new Incidence_request_type();
      _interior_edge_set = new Interior_edge_set_type();
    }

  //-------------------- DESTRUCTOR -----------------------------------

  ~Local_selection_vertex_base_3()
    {
      if (_incident_border != NULL)
	{
	  delete _incident_border->first;
	  delete _incident_border->second;
	  delete _incident_border;
	}
      if (_incidence_request != NULL)
	delete _incidence_request;
      if (_interior_edge_set != NULL)
	delete _interior_edge_set;
    }

  //-------------------- MEMBER FUNCTIONS -----------------------------

  inline void re_init()
    {
      if (_incident_border != NULL)
	{
	  delete _incident_border->first;
	  delete _incident_border->second;
	  delete _incident_border;
	}
      if (_incidence_request != NULL)
	delete _incidence_request;
      //      if (_interior_edge_set != NULL)
      //	delete _interior_edge_set;

      _incident_border = new Intern_successors_type(new Next_border_elt(),
						    new Next_border_elt());
      _incident_border->first->first = NULL;
      _incident_border->second->first = NULL;
      _incidence_request = new Incidence_request_type();
      //      _interior_edge_set = new Interior_edge_set_type();
      _mark = -1;
      _post_mark = -1;
      // Attention ne pas toucher a visu index sous peine de tout casser dans la 
      // visu de la partie deja reconstruite
    }

  //-------------------------------------------------------------------

  inline bool is_on_border(const int& i)
    {
      if (_incident_border == NULL) return false; //vh is interior
      if (_incident_border->first->first != NULL)
	{ 
	  if (_incident_border->second->first != NULL)
	    return ((_incident_border->first->second.second == i)||
		    (_incident_border->second->second.second == i));
	  return (_incident_border->first->second.second == i);
	}
      return false; //vh is still exterior   
    }

  inline Next_border_elt* get_next_on_border(const int& i)
    { 
      if (_incident_border == NULL) return NULL; //vh is interior
      if (_incident_border->first->first != NULL)
	if (_incident_border->first->second.second == i)
	  return _incident_border->first;
      if (_incident_border->second->first != NULL)
	if (_incident_border->second->second.second == i)
	  return _incident_border->second;
      return NULL;
    }


  inline void remove_border_edge(void* v)
    {
      if (_incident_border != NULL)
	{
	  if (_incident_border->second->first == v)
	    {
	      _incident_border->second->first = NULL;
	      if (_interior_edge_set != NULL){
		_interior_edge_set->insert(v);
		std::cout << "|set| = " << _interior_edge_set->size() << std::endl;
	      }
	      return;
	    }
	  if (_incident_border->first->first == v)
	    {
	      if (_incident_border->second->first != NULL)
		{
		  Next_border_elt* tmp = _incident_border->first;
		  _incident_border->first = _incident_border->second;
		  _incident_border->second = tmp;
		  _incident_border->second->first = NULL;
		  if (_interior_edge_set != NULL){
		    _interior_edge_set->insert(v);
		    std::cout << "|set| = " << _interior_edge_set->size() << std::endl;
		  }
		  return;
		}
	      else
		{ 
		  _incident_border->first->first = NULL;
		  if (_interior_edge_set != NULL){
		    _interior_edge_set->insert(v);
		    std::cout << "|set| = " << _interior_edge_set->size() << std::endl;
		  }
		  return;
		}
	    }
	}
    }

  inline bool is_border_edge(void* v)
    { 
      if (_incident_border == NULL) return false;
      return ((_incident_border->first->first == v)||
	      (_incident_border->second->first == v));
    }

  inline Next_border_elt* get_border_elt(void* v)
    {
      if (_incident_border == NULL) return NULL;
      if (_incident_border->first->first == v) return _incident_border->first;
      if (_incident_border->second->first == v) return _incident_border->second;
      return NULL; 
    }

  inline Next_border_elt* first_incident()
    {
      if (_incident_border == NULL) return NULL;
      return _incident_border->first;
    }

  inline Next_border_elt* second_incident()
    {
      if (_incident_border == NULL) return NULL;
      return _incident_border->second;
    }


  inline  void set_next_border_elt(const Next_border_elt& elt)
    {
      if (_incident_border->first->first == NULL)
	*_incident_border->first = elt;
      else
	{
	  if (_incident_border->second->first != NULL)
	    std::cerr << "+++probleme de MAJ du bord <Vertex_base>" << std::endl;
	  *_incident_border->second = elt;
	}
    }

  //-------------------------------------------------------------------
  // pour gerer certaines aretes interieures: a savoir celle encore connectee au 
  // bord (en fait seule, les aretes interieures reliant 2 bords nous
  // interressent...)

  inline bool is_interior_edge(void* v)
    {
      if (_interior_edge_set == NULL) return true;
      return (_interior_edge_set->find(v) != _interior_edge_set->end());
    }

  inline void set_interior_edge(void* v)
    {
      _interior_edge_set->insert(v);  
      std::cout << "|set| = " << _interior_edge_set.size() << std::endl;    
    }

  inline void remove_interior_edge(void* v)
    {
      typename Interior_edge_set_type::iterator
	it_tmp = _interior_edge_set->find(v);
      if(it_tmp != _interior_edge_set->end())
	_interior_edge_set->erase(it_tmp);
    }

  //-------------------------------------------------------------------
  
  inline void set_incidence_request(const Incidence_request_elt& e)
    {
      _incidence_request->push_back(e);
    }

  inline bool is_incidence_requested()
    {
      if (_incidence_request == NULL) return false;
      return (!_incidence_request->empty());
    }
  
  inline Incidence_request_iterator incidence_request_begin()
    {
      return _incidence_request->begin();
    }

  inline Incidence_request_iterator get_incidence_request_end()
    {
      return _incidence_request->end();
    }

  inline void erase_incidence_request()
    { 
      if (_incidence_request != NULL)
	_incidence_request->clear();
    }
  

  //-------------------------------------------------------------------

  inline bool is_on_border()
    {
      return (_mark > 0);
    }

  inline bool not_interior()
    {
      return (_mark != 0);
    }

  inline int number_of_incident_border()
    {
      return _mark;
    }

  inline bool is_exterior()
    {
      return (_mark < 0);
    }
  
  //-------------------------------------------------------------------

  inline void inc_mark()
    {
      if (_mark==-1)
	_mark=1;
      else
	_mark++;
    }

  inline void dec_mark()
    {
      _mark--;
      if(_mark == 0)
	{
	  delete _incident_border->first;
	  delete _incident_border->second;
	  delete _incident_border;
	  _incident_border = NULL;
	  delete _incidence_request;
	  //	  delete _interior_edge_set;
	  _incidence_request = NULL;
	  //	  _interior_edge_set = NULL;
	}
    }

  //-------------------------------------------------------------------

  inline void set_visu_index(const int& i)
    {
      _visu_index = i;
    }
  
  inline int get_visu_index()
    {
      return _visu_index;
    }

  //-------------------------------------------------------------------

  inline void set_post_mark(const int& i)
    {
      _post_mark = i;
    }
  
  inline bool is_post_marked(const int& i)
    {
      return (_post_mark == i);
    }
};

//-------------------------------------------------------------------
#endif LOCAL_SELECTION_VERTEX_BASE_3_H
//-------------------------------------------------------------------
