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


template <class VertexBase>
class Local_selection_vertex_base_3 : public VertexBase
{
public:

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename VertexBase::template Rebind_TDS<TDS2>::Other  Vb2;
    typedef Local_selection_vertex_base_3<Vb2>                    Other;
  };


  typedef VertexBase Base;
  typedef typename Base::Vertex_handle Vertex_handle;
  typedef typename Base::Cell_handle Cell_handle;
  typedef typename VertexBase::Point Point;
  typedef double coord_type;
  
  typedef CGAL::Triple< Cell_handle, int, int > Edge;
  typedef std::pair< Edge, int > Edge_IFacet;
  typedef std::pair< Edge_IFacet, Edge_IFacet > IO_edge_type;

  typedef coord_type criteria;

  typedef std::pair< criteria, IO_edge_type > Radius_edge_type;
  typedef std::pair< Radius_edge_type, int > Border_elt;
  typedef std::pair< Vertex_handle, Border_elt > Next_border_elt;

private:

  //par convention je remplis d'abord first et si necessaire second...
  typedef std::pair< Next_border_elt*,  Next_border_elt*> Intern_successors_type;

public:

  typedef std::pair< criteria, IO_edge_type* > Radius_ptr_type;
  typedef std::pair< Vertex_handle, Vertex_handle > Edge_like;
  typedef std::pair< criteria,  Edge_like > Incidence_request_elt;
  typedef std::list< Incidence_request_elt > Incidence_request_type;
  typedef typename Incidence_request_type::iterator Incidence_request_iterator;
  
  //  typedef std::set< void* > Interior_edge_set_type;
  
  //-------------------- DATA MEMBERS ---------------------------------

private:

  int _mark;
  int _post_mark;
  Intern_successors_type* _incident_border;

  // Instead of having a set per vertex, there is a global list.
  static std::list<Vertex_handle> interior_edges;
  // and two iterators per vertex in this list
  // Note that ie_last is not past the end
  // ie_first == ie_last == interior_edge.end()  iff  the set is empty
  typename std::list<Vertex_handle>::iterator ie_first, ie_last;


  // We do the same for the incidence requests
  static std::list< Incidence_request_elt > incidence_requests;
  typename std::list< Incidence_request_elt >::iterator ir_first, ir_last;
  //-------------------- CONSTRUCTORS ---------------------------------

public:

  Local_selection_vertex_base_3()
    : VertexBase(), _mark(-1), 
    _post_mark(-1), 
    ie_first(interior_edges.end()), ie_last(interior_edges.end()),
    ir_first(incidence_requests.end()), ir_last(incidence_requests.end())
    {
      _incident_border = new Intern_successors_type(new Next_border_elt(),
						    new Next_border_elt());
      _incident_border->first->first = NULL; 
      _incident_border->second->first = NULL;
    }
  
  Local_selection_vertex_base_3(const Point & p)
    : VertexBase(p), _mark(-1), 
    _post_mark(-1), 
    ie_first(interior_edges.end()), ie_last(interior_edges.end()),
    ir_first(incidence_requests.end()), ir_last(incidence_requests.end())
    {
      _incident_border = new Intern_successors_type(new Next_border_elt(),
						    new Next_border_elt());
      _incident_border->first->first = NULL;
      _incident_border->second->first = NULL;
    }
  
  Local_selection_vertex_base_3(const Point & p, Cell_handle f)
    : VertexBase(p, f), _mark(-1), 
    _post_mark(-1), 
    ie_first(interior_edges.end()), ie_last(interior_edges.end()),
    ir_first(incidence_requests.end()), ir_last(incidence_requests.end())
    {
      _incident_border = new Intern_successors_type(new Next_border_elt(),
						    new Next_border_elt());
      _incident_border->first->first = NULL;
      _incident_border->second->first = NULL;
    }

  Local_selection_vertex_base_3(Cell_handle f)
    : VertexBase(f), _mark(-1), 
    _post_mark(-1), 
    ie_first(interior_edges.end()), ie_last(interior_edges.end()),
    ir_first(incidence_requests.end()), ir_last(incidence_requests.end())
    {
      _incident_border = new Intern_successors_type(new Next_border_elt(),
						    new Next_border_elt());
      _incident_border->first->first = NULL;
      _incident_border->second->first = NULL;
    }

	Local_selection_vertex_base_3(const Local_selection_vertex_base_3& other)
    : VertexBase(), _mark(-1), 
    _post_mark(-1), 
    ie_first(interior_edges.end()), ie_last(interior_edges.end()),
    ir_first(incidence_requests.end()), ir_last(incidence_requests.end())
	{
      _incident_border = new Intern_successors_type(new Next_border_elt(),
						    new Next_border_elt());
      _incident_border->first->first = NULL;
      _incident_border->second->first = NULL;
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
      if(ir_first != incidence_requests.end()){
	assert(ir_last != incidence_requests.end());
	typename std::list< Incidence_request_elt >::iterator b(ir_first), e(ir_last);
	e++;
	incidence_requests.erase(b, e);
      }

      if(ie_first != interior_edges.end()){
	assert(ie_last != interior_edges.end());
	typename std::list<Vertex_handle>::iterator b(ie_first), e(ie_last);
	e++;
	interior_edges.erase(b, e);
      }
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

      if(ir_first != incidence_requests.end()){
	assert(ir_last != incidence_requests.end());
	typename std::list< Incidence_request_elt >::iterator b(ir_first), e(ir_last);
	e++;
	incidence_requests.erase(b, e);
	ir_first = incidence_requests.end();
	ir_last = incidence_requests.end();
      }

      _incident_border = new Intern_successors_type(new Next_border_elt(),
						    new Next_border_elt());
      _incident_border->first->first = NULL;
      _incident_border->second->first = NULL;
      _mark = -1;
      _post_mark = -1;
    }

  //-------------------------------------------------------------------

  inline bool is_on_border(const int& i) const
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

  inline Next_border_elt* get_next_on_border(const int& i) const
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


  inline void remove_border_edge(Vertex_handle v)
    {
      if (_incident_border != NULL)
	{
	  if (_incident_border->second->first == v)
	    {
	      _incident_border->second->first = NULL;
	      set_interior_edge(v);
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
		  set_interior_edge(v);
		  return;
		}
	      else
		{ 
		  _incident_border->first->first = NULL;
		  set_interior_edge(v);
		  return;
		}
	    }
	}
    }

  inline bool is_border_edge(Vertex_handle v) const
    { 
      if (_incident_border == NULL) return false;
      return ((_incident_border->first->first == v)||
	      (_incident_border->second->first == v));
    }

  inline Next_border_elt* get_border_elt(Vertex_handle v) const
    {
      if (_incident_border == NULL) return NULL;
      if (_incident_border->first->first == v) return _incident_border->first;
      if (_incident_border->second->first == v) return _incident_border->second;
      return NULL; 
    }

  inline Next_border_elt* first_incident() const
    {
      if (_incident_border == NULL) return NULL;
      return _incident_border->first;
    }

  inline Next_border_elt* second_incident() const
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

  inline bool is_interior_edge(Vertex_handle v) const
    {

      bool r1;
      if(ie_first == interior_edges.end()){
	r1 = false;
      }else {
	typename std::list<Vertex_handle>::iterator b(ie_first), e(ie_last);
	e++;
	typename std::list<Vertex_handle>::iterator r = std::find(b, e, v);
	r1 = ( r != e);
      }

      return r1;
    }

  inline void set_interior_edge(Vertex_handle v)
    {
      if(ie_last == interior_edges.end()){ // empty set
	assert(ie_first == ie_last);
	ie_last = interior_edges.insert(ie_last, v);
	ie_first = ie_last;
      } else {
	typename std::list<Vertex_handle>::iterator e(ie_last);
	e++;
#ifdef DEBUG
	typename std::list<Vertex_handle>::iterator r = std::find(ie_first, e, v);
#endif
	assert(r == e);
	ie_last = interior_edges.insert(e, v);
      }
    }

  inline void remove_interior_edge(Vertex_handle v)
    {
      if(ie_first == interior_edges.end()){
	assert(ie_last == ie_first);
      } else if(ie_first == ie_last){ // there is only one element
	if(*ie_first == v){
	  interior_edges.erase(ie_first);
	  ie_last = interior_edges.end();
	  ie_first = ie_last;
	}
      } else {
	typename std::list<Vertex_handle>::iterator b(ie_first), e(ie_last);
	e++;
	typename std::list<Vertex_handle>::iterator r = std::find(b, e, v);
	if(r != e){
	  if(r == ie_first){
	    ie_first++;
	  }
	  if(r == ie_last){
	    ie_last--;
	  }
	  interior_edges.erase(r);
	}
      }
    }

  //-------------------------------------------------------------------
  
  inline void set_incidence_request(const Incidence_request_elt& ir)
    {
      if(ir_last == incidence_requests.end()){
	assert(ir_first == ir_last);
	ir_last = incidence_requests.insert(ir_last, ir);
	ir_first = ir_last;
      } else {
	typename std::list<Incidence_request_elt>::iterator e(ir_last);
	e++;
	ir_last = incidence_requests.insert(e, ir);
      }
    }

  inline bool is_incidence_requested() const
    {
      if(ir_last == incidence_requests.end()){
	assert(ir_first == incidence_requests.end());
      }
      return (ir_last != incidence_requests.end());
    }
  
  inline Incidence_request_iterator incidence_request_begin()
    {
      return ir_first;
    }

  inline Incidence_request_iterator get_incidence_request_end()
    {
      if(ir_last != incidence_requests.end()){
	assert(ir_first != incidence_requests.end());
	Incidence_request_iterator it(ir_last);
	it++;
	return it;
      }
      return ir_last;
    }

  inline void erase_incidence_request()
    { 
      if(ir_last != incidence_requests.end()){
	assert(ir_first != incidence_requests.end());
	ir_last++;
	incidence_requests.erase(ir_first, ir_last);
	ir_first = incidence_requests.end();
	ir_last = incidence_requests.end();
      }
    }
  

  //-------------------------------------------------------------------

  inline bool is_on_border() const
    {
      return (_mark > 0);
    }

  inline bool not_interior() const
    {
      return (_mark != 0);
    }

  inline int number_of_incident_border() const
    {
      return _mark;
    }

  inline bool is_exterior() const
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
	  erase_incidence_request();
	}
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

template <class Gt>
std::list<typename Local_selection_vertex_base_3<Gt>::Vertex_handle> Local_selection_vertex_base_3<Gt>::interior_edges;

template <class Gt>
std::list<typename Local_selection_vertex_base_3<Gt>::Incidence_request_elt> Local_selection_vertex_base_3<Gt>::incidence_requests;

//-------------------------------------------------------------------
#endif //LOCAL_SELECTION_VERTEX_BASE_3_H
//-------------------------------------------------------------------
