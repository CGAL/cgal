// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 1999, October 13
//
// file          : include/CGAL/Topological_map.h
// package       : pm (4.08)
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Iddo Hanniel <hanniel@math.tau.ac.il>
//                 Oren Nechushtan <theoren@math.tau.ac.il>
//
//
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter       : 
// ======================================================================
#ifndef  CGAL_TOPOLOGICAL_MAP_H
#define  CGAL_TOPOLOGICAL_MAP_H

#ifndef CGAL_PLANAR_MAP_MISC_H
#include  <CGAL/Planar_map_2/Planar_map_misc.h>
#endif


CGAL_BEGIN_NAMESPACE


template <class Dcel >
class Topological_map
{

protected:
  typedef typename Dcel::Vertex	         dvertex;
  typedef typename Dcel::Halfedge        dedge;
  typedef typename Dcel::Face	         dface;

public:
  typedef typename Dcel::Size            Size;
  typedef typename Dcel::size_type       size_type;
  typedef typename Dcel::difference_type difference_type;
  typedef typename Dcel::difference_type Difference;
  
  typedef typename Dcel::iterator_category      iterator_category;

  //These were taken from Lutz's polyhedron probably I'll need them too  
#ifndef __SUNPRO_CC 
//in SunPro protected members are not reachable from nested classes
protected:
#endif
  // These extra (internal) typedefs are necessary to make
  // SunPro CC 4.2 happy. (And they are used.)
  typedef typename  Dcel::Vertex_iterator          TR_VI;
  typedef typename  Dcel::Vertex_const_iterator    TR_C_VI;
  typedef typename  Dcel::Halfedge_iterator        TR_HI;
  typedef typename  Dcel::Halfedge_const_iterator  TR_C_HI;

  typedef typename  Dcel::Face_iterator           TR_FI;
  typedef typename  Dcel::Face_const_iterator     TR_C_FI;
    
  typedef typename Dcel::Face::Holes_iterator TR_HOI;
  typedef typename Dcel::Face::Holes_const_iterator TR_C_HOI;


public:
  class Vertex;
  class Halfedge;
  class Face;

//in VC class's private members aren't available in nested classes
#if _MSC_VER>=1100
  friend Face;
#endif

  typedef _Polyhedron_iterator<
  TR_VI,
  Vertex,
  Difference, iterator_category> Vertex_iterator;
  
  typedef _Polyhedron_const_iterator<
  TR_C_VI, TR_VI,
  Vertex,
  Difference, iterator_category>       Vertex_const_iterator;
  
  typedef _Polyhedron_iterator<
  TR_HI,
  Halfedge,
  Difference, iterator_category>       Halfedge_iterator;
  
#ifndef __SUNPRO_CC
  typedef _Polyhedron_const_iterator<
  TR_C_HI, TR_HI,
  Halfedge,
  Difference, iterator_category>       Halfedge_const_iterator;
#endif // __SUNPRO_CC //
  
  typedef _Polyhedron_iterator<
  TR_FI,
  Face,
  Difference, iterator_category>       Face_iterator;
  
  typedef _Polyhedron_const_iterator<
  TR_C_FI, TR_FI,
  Face,
  Difference, iterator_category>       Face_const_iterator;


  typedef _Polyhedron_facet_circ<
  Halfedge,
  Halfedge_iterator,
  Forward_circulator_tag>            Ccb_halfedge_circulator;

  typedef _Polyhedron_vertex_circ<
  Halfedge,
  Halfedge_iterator,
  Forward_circulator_tag>            Halfedge_around_vertex_circulator;


#ifdef __SUNPRO_CC
  typedef _Polyhedron_halfedge_const_iterator<
  TR_C_HI, TR_HI,
  Ccb_halfedge_circulator,
  Halfedge_around_vertex_circulator,
  Halfedge,
  Difference, iterator_category>       Halfedge_const_iterator;
#endif // __SUNPRO_CC //
 
  
  typedef _Polyhedron_facet_const_circ<
  Halfedge,
  Halfedge_const_iterator,
  Forward_circulator_tag>       Ccb_halfedge_const_circulator;
  
  
  typedef _Polyhedron_vertex_const_circ<
  Halfedge,
  Halfedge_const_iterator,
  Forward_circulator_tag>      Halfedge_around_vertex_const_circulator;




  typedef _Polyhedron_iterator<
  TR_HOI,
  Ccb_halfedge_circulator,
  Difference, std::bidirectional_iterator_tag>       Holes_iterator;
  
  typedef _Polyhedron_const_iterator<
  TR_C_HOI, TR_HOI,
  Ccb_halfedge_const_circulator,
  Difference, std::bidirectional_iterator_tag>       Holes_const_iterator;


  
  typedef Vertex_iterator Vertex_handle;
  typedef Halfedge_iterator Halfedge_handle;
  typedef Face_iterator Face_handle;

  typedef Vertex_const_iterator Vertex_const_handle;
  typedef Halfedge_const_iterator Halfedge_const_handle;
  typedef Face_const_iterator Face_const_handle;
  

  /**************************************************************/
  /********************* V e r t e x ****************************/
  /**************************************************************/
  /**************************************************************/
  class Vertex : public Dcel::Vertex
  {
  public:
    typedef typename Dcel::Vertex Base;

    Vertex()
    {}
    
  
    //    Vertex(dvertex *v) : dvertex(*v) 
    Vertex(Base *v) : Base(*v) 
    {}


    //    Vertex(dvertex& v) : dvertex(v) 
    Vertex(Base& v) : Base(v) 
    {}


    bool is_incident_edge(Halfedge_const_handle e) const
    {
      return (&(*(e->source())) == this) || (&(*(e->target())) == this);
    }

    bool is_incident_face(Face_const_handle f) const
    {
      //      const typename Dcel::Halfedge* e_ptr=(dvertex::halfedge());
      const typename Dcel::Halfedge* e_ptr=(Base::halfedge());
      
      do{
        if ( e_ptr->face() == &(*f) )
          return true;

        e_ptr=e_ptr->opposite()->next();
        //      } while (e_ptr != dvertex::halfedge());
      } while (e_ptr != Base::halfedge());
    
      return false;
    }
    
    unsigned int degree() const
    {
      //      const typename Dcel::Halfedge* e=dvertex::halfedge();
      const typename Dcel::Halfedge* e=Base::halfedge();
      const typename Dcel::Halfedge* e1=e;
      int n=0;
      if (e1!=NULL)
      {
        do {
          n++;
          e1=e1->next()->opposite();
        } while (e1!=e);
      }
    return n;
    }

    Halfedge_around_vertex_circulator incident_halfedges() 
    {
      //return Halfedge_around_vertex_circulator(TR_HI(dvertex::halfedge()));
      return Halfedge_around_vertex_circulator(TR_HI(Base::halfedge()));
    }

    Halfedge_around_vertex_const_circulator incident_halfedges() const 
    {
      //return 
      //Halfedge_around_vertex_const_circulator(TR_C_HI(dvertex::halfedge())); 
      return 
	Halfedge_around_vertex_const_circulator(TR_C_HI(Base::halfedge())); 
    }

  private:
    //    void set_halfedge(Halfedge* h) {dvertex::set_halfedge(h);}
    void set_halfedge(Halfedge* h) {Base::set_halfedge(h);}
    
  };
  
  /**************************************************************/
  /**************************************************************/
  /********************* Halfedge *******************************/
  /**************************************************************/
  /**************************************************************/
  class Halfedge : public Dcel::Halfedge
  {
  public:
    typedef typename Dcel::Halfedge Base;

    Halfedge()
    {}
    
    //    Halfedge(dedge *e) : dedge(*e) {}
    Halfedge(Base *e) : Base(*e) {}

    //    Halfedge(dedge& e) : dedge(e) {}
    Halfedge(Base& e) : Base(e) {}


    Vertex_handle source()
    {
      //      return TR_VI(dedge::opposite()->vertex()); 
      return TR_VI(Base::opposite()->vertex()); 
    }

    Vertex_const_handle source() const
    {
      //      return TR_C_VI(dedge::opposite()->vertex()); 
      return TR_C_VI(Base::opposite()->vertex()); 
    }
    
    Vertex_handle target()	
    {
      //return TR_VI(dedge::vertex());  
      return TR_VI(Base::vertex());  
    }

    Vertex_const_handle target() const	
    {
      //return TR_C_VI(dedge::vertex()); 
      return TR_C_VI(Base::vertex()); 
    }
    
    Face_handle face() 
    {
      //return TR_FI(dedge::face()); 
      return TR_FI(Base::face()); 
    }

    Face_const_handle face() const
    {
      //return TR_C_FI(dedge::face()); 
      return TR_C_FI(Base::face()); 
    }
    
    Halfedge_handle twin() 
    {
      //return TR_HI(dedge::opposite()); 
      return TR_HI(Base::opposite()); 
    }

    Halfedge_const_handle twin() const 
    {
      //return TR_C_HI(dedge::opposite()); 
      return TR_C_HI(Base::opposite()); 
    }
    
    Halfedge_handle next_halfedge() 
    { 
      //return TR_HI(dedge::next()); 
      return TR_HI(Base::next()); 
    }

    Halfedge_const_handle next_halfedge() const 
    { 
      //return TR_C_HI(dedge::next()); 
      return TR_C_HI(Base::next()); 
    }
    

    Ccb_halfedge_circulator ccb()
    { 
      return Ccb_halfedge_circulator(TR_HI(this)); 
    }

    Ccb_halfedge_const_circulator ccb() const
    { 
      return Ccb_halfedge_const_circulator(TR_C_HI(this)); 
    }


  private: //hide from users the underlying structure
    //void  set_next( Halfedge* h)     { dedge::set_next(h);}
    void  set_next( Halfedge* h)     { Base::set_next(h);}
    //void  set_vertex( Vertex* ve)    { dedge::set_vertex(ve);}
    void  set_vertex( Vertex* ve)    { Base::set_vertex(ve);}
    //void  set_face( Face* face)      { dedge::set_face(face);}
    void  set_face( Face* face)      { Base::set_face(face);}
     
    };
  
  /**************************************************************/
  /**************************************************************/
  /********************* F a c e ********************************/
  /**************************************************************/
  /**************************************************************/
  class Face : public Dcel::Face
  {
  public:
    typedef typename Dcel::Face Base;

#ifndef _MSC_VER
    // the following two typedefs are needed for compilation on irix64 (CC7.30)
    typedef Topological_map<Dcel>::Holes_iterator 
    Holes_iterator; 
    typedef Topological_map<Dcel>::Holes_const_iterator 
    Holes_const_iterator; 
#endif
    
    Face()
    {}
    
    Face(Base *f) : Base(*f) {}

    Face(Base& f) : Base(f) {} 


    bool is_unbounded() const  
    {
      // face is not bounded iff it has no outer boundary
      return (Base::halfedge() == NULL);
    }
    

    // must explicitly say Topological_map::Holes_iterator. Otherwise
    // the CC compiler thinks we mean Base::Holes_iterator.
    Holes_iterator holes_begin() 
    {
      return TR_HOI(Base::holes_begin());
    }
    Holes_const_iterator holes_begin() const
    {
      return TR_C_HOI(Base::holes_begin()); 
    }
    
    Holes_iterator holes_end() 
    {
      return TR_HOI(Base::holes_end()); 
    }

    Holes_const_iterator holes_end() const 
    {
      return TR_C_HOI(Base::holes_end()); 
    }

    bool is_halfedge_on_inner_ccb(Halfedge_const_handle e) 
#if !(_MSC_VER>=1100) 
		const
#endif

    {
      dedge* dummy;
      return (/*Topological_map<Dcel>::*/
	      is_halfedge_on_inner_ccb(&(*e),
				       (typename Dcel::Face*)this,
				       dummy)); 
// oren corrections (dummy parameter - for bug in MSC)
    }
    
    bool is_halfedge_on_outer_ccb(Halfedge_const_handle e) 
#if !(_MSC_VER>=1100) 
		const
#endif
    {
      dedge* dummy;
      return /*Topological_map<Dcel>::*/
	is_halfedge_on_outer_ccb(&(*e),this,dummy);
// oren corrections (dummy parameter - bug in MSC)
    }
	
    bool does_outer_ccb_exist() const
    {
      return Base::halfedge() != NULL;
    }
    
    Halfedge_handle halfedge_on_outer_ccb() 
    {
      CGAL_precondition(does_outer_ccb_exist());      
      return TR_HI(Base::halfedge()); 
    }

    Halfedge_const_handle halfedge_on_outer_ccb() const 
    {
      CGAL_precondition(does_outer_ccb_exist());      
      return TR_C_HI(Base::halfedge());  
    }

    
    Ccb_halfedge_circulator outer_ccb() 
    {
      CGAL_precondition(does_outer_ccb_exist()); 
      return (halfedge_on_outer_ccb())->ccb();
    }
    
    Ccb_halfedge_const_circulator outer_ccb() const
    {
      CGAL_precondition(does_outer_ccb_exist()); 
      return (halfedge_on_outer_ccb())->ccb();
    }

  private:
    void  set_halfedge( Halfedge* h) { Base::set_halfedge(h);}

    void add_hole(Halfedge* h) {Base::add_hole(h);}

  };
  

  /**************************************************************/
  /**************************************************************/
  /********************* Topological_map *******************/
  /**************************************************************/
  /**************************************************************/

  Topological_map()
    : d()
  {
    u_face = d.new_face();
    u_face->set_halfedge(NULL);
  }
  


//INSERTIONS


//insertion from vertices may invalidate the holes containers (since a new
//face can be created which contains some of the holes of the old face.
//to validate the containers the function move_hole should be called.

//if a new face is created, returns the halfedge on the new face.
//the new face will ALWAYS (by definition) be the one from v1->v2 (==e1) 
//this sets a requirement on the user!!! (since a topological pm can't 
//know if v1->v2 is on the new face or v2->v1)

//since there can be a number of edges from the vertices, and we can't
//know topologically which will be in the new face and which will be
//outside, this should be provided by the user - defining prev1 and
//prev2 (previous to v1,v2 resp.)  the planar map will define this
//geometrically using the one previous to the curve clock wise.
  
  Halfedge_handle insert_at_vertices( Halfedge_handle previous1, 
                                      Halfedge_handle previous2) ;
  //moves ONE hole from f1 to f2
  void move_hole(Holes_iterator e, Face_handle f1, Face_handle f2);


  Halfedge_handle insert_from_vertex(Halfedge_handle previous)   ;

  Halfedge_handle insert_in_face_interior(Face_handle f)   ;


  Halfedge_handle split_edge (Halfedge_handle e )   ;
  
  Halfedge_handle merge_edge (Halfedge_handle e1, Halfedge_handle e2) ;

//returns the merged face after the deletion
  Face_handle remove_edge(Halfedge_handle e) ;

    
  Face_handle unbounded_face()
  {
    return TR_FI(u_face); 
  }

  Face_const_handle unbounded_face() const {
    return TR_C_FI(u_face);
  }
  
  Face_iterator faces_begin()
  { 
    return TR_FI(d.faces_begin()); 
  }

  Face_const_iterator faces_begin() const
  { 
    return TR_C_FI(d.faces_begin()); 
  }
  
  Face_iterator faces_end() 
  { 
    return TR_FI(d.faces_end()); 
  }

  Face_const_iterator faces_end() const 
  { 
    return TR_C_FI(d.faces_end()); 
  }
  
  Halfedge_iterator halfedges_begin() 
  {
    return TR_HI(d.halfedges_begin()); 
  }
  
  Halfedge_const_iterator halfedges_begin() const
  {
    return TR_C_HI(d.halfedges_begin()); 
  }
  
  Halfedge_iterator halfedges_end() 
  { 
    return TR_HI(d.halfedges_end()); 
  }

  Halfedge_const_iterator halfedges_end() const
  { 
    return TR_C_HI(d.halfedges_end()); 
  }
  
  Vertex_iterator vertices_begin() 
  { 
    return TR_VI(d.vertices_begin()); 
  }

  Vertex_const_iterator vertices_begin() const
  { 
    return TR_C_VI(d.vertices_begin()); 
  }
  
  Vertex_iterator vertices_end()
  {
    return TR_VI(d.vertices_end()); 
  }

  Vertex_const_iterator vertices_end() const
  {
    return TR_C_VI(d.vertices_end()); 
  }

  void clear()
  {
    d.delete_all();
    u_face = NULL;
    u_face = d.new_face();
    u_face->set_halfedge(NULL);
  }
  
  void assign(const Topological_map<Dcel> &tpm)
  {
	  d.delete_all();
	  u_face = NULL;
	  u_face = (Face*)d.assign(tpm.d, tpm.u_face);
  }

  
  /**************************************************************/
  /**************************************************************/
  /*********** Topological_map_2  - validity checks ********/
  /**************************************************************/
  /**************************************************************/
  
  bool is_valid(Vertex_const_handle v 
#if _MSC_VER>=1100 // VC template bug
	  ,Vertex* dummy=NULL
#endif
	  ) const
  {
    bool valid = true;
    
    // check if every edge from v has v as its target
    Halfedge_around_vertex_const_circulator ec = v->incident_halfedges();
    Halfedge_around_vertex_const_circulator start = ec;
    do {
      if ((*ec).target() != v)
        {
          valid = false;
        }
      ++ec;
    } while (ec != start);
    
    return valid;
  }

  
  bool is_valid(Halfedge_const_handle e
#if _MSC_VER>=1100 // VC template bug
	  ,Halfedge* dummy=NULL
#endif
	  ) const
  {
    //check relations with next
    if (e->target() != e->next_halfedge()->source())
      return false;
    
    return true;
  }
  
  bool is_valid(const Ccb_halfedge_const_circulator& start, 
                Face_const_handle f) const
  {
    bool valid = true;
    Ccb_halfedge_const_circulator circ = start;
    
    do {
      if ((*circ).face() != f)
        valid = false;
      ++circ;
    } while (circ != start);
    
    return valid;
  }
  
  bool is_valid(Face_const_handle f
#if _MSC_VER>=1100 // VC template bug
	  ,Face* dummy=NULL
#endif
	  ) const

  {
    bool valid = true;
    
    // check if all edges of f (on all ccb's) refer to f (as their face)
    Holes_const_iterator iccbit;
    Ccb_halfedge_const_circulator ccb_circ;
    
    if (f->does_outer_ccb_exist())
      {
        ccb_circ = f->outer_ccb();
        if (!is_valid(ccb_circ, f))
          valid = false;
      }
    
    for (iccbit = f->holes_begin(); iccbit != f->holes_end(); ++iccbit)
      {
        ccb_circ = *iccbit;
        if (!is_valid(ccb_circ, f))
          valid = false;
      }
    
    return valid;
  }
  
bool is_valid() const
  {
    bool valid = true;
    
    Vertex_const_iterator vi;
    for (vi = vertices_begin(); vi != vertices_end(); ++vi)
      {
        if (!is_valid(vi)) 
          {//an iterator is used for a handle
            valid = false;
          }
      }
    
    Halfedge_const_iterator ei;
    for (ei = halfedges_begin(); ei != halfedges_end(); ++ei)
      {
        if (!is_valid(ei))
          valid = false;
      }	
    
    Face_const_iterator fi;
    for (fi = faces_begin(); fi != faces_end(); ++fi)
      {
        if (!is_valid(fi))
          valid = false;
      }
    
    return valid;
  }



  /******************************************************************/
  /*********************   counting functions ***********************/
  /******************************************************************/

  Size number_of_faces() const
  {
    return d.size_of_faces();
  }
   

  //counts every halfedge (i.e always even)
  Size number_of_halfedges() const
  {
    return d.size_of_halfedges();
  }
   
  Size number_of_vertices() const
  {
    return d.size_of_vertices();
  }
  
//--------------------------------------------------------------- 
//PRIVATE IMPLEMENTATIONS    
//---------------------------------------------------------------

private:
  //find previous halfedge of given halfedge
  typename Dcel::Halfedge* get_prev( typename  Dcel::Halfedge* de) 
  {
    //find previous halfedge of de
    typename Dcel::Halfedge* prev=de;
    do {
      if  (prev->next()==de) break;
      prev=prev->next();
    }while (prev!=de);  //assumes a loop to itself is unacceptable
    CGAL_assertion(prev!=de);
    
    return prev;
}


bool find_and_erase_hole(typename Dcel::Halfedge* e, typename Dcel::Face* f)
{
  typename Dcel::Face::Holes_iterator it=f->holes_begin();
  bool ret=false;
  for (; it!=f->holes_end(); ++it) {
      if ( (*it)==e ) {
        f->erase_hole(it);
        ret=true;
        break;
      }
    }
    return ret;
}

protected:
  Dcel d;
  dface* u_face;
};




/////////////////////////////////////////////////////////////////
//                  implementation
/////////////////////////////////////////////////////////////////

template <class Dcel>
bool is_halfedge_on_inner_ccb(const typename Dcel::Halfedge* e, 
                              typename Dcel::Face* f, 
                              typename Dcel::Halfedge* &iccb) 
{
  typedef typename Dcel::Vertex	         dvertex;
  typedef typename Dcel::Halfedge        dedge;
  typedef typename Dcel::Face	         dface;
  //move over all holes in face and check 
  for (typename dface::Holes_iterator dhi=f->holes_begin();
       dhi!=f->holes_end();
       ++dhi)
    {
      iccb=(*dhi);
      const typename Dcel::Halfedge* aux=iccb;
      do
        {
          if (aux==e) 
            return true;
          aux=aux->next() ;
        } while (aux!=(*dhi));
    }
  iccb=NULL;
  return false;
}

//returns the representative halfedge of outer ccb in occb, or false
template <class Dcel>
bool is_halfedge_on_outer_ccb(const typename Dcel::Halfedge* e, 
                              typename Dcel::Face* f, 
                              typename Dcel::Halfedge* &occb)
{
  typedef typename Dcel::Vertex	         dvertex;
  typedef typename Dcel::Halfedge        dedge;
  typedef typename Dcel::Face	         dface;
  occb=f->halfedge();
  if ( occb ) {  //if f isn't the unbounded face
  const typename Dcel::Halfedge* aux=occb;    
    do
      {
        if (aux==e)
          return true;
        aux=aux->next();
      }while (aux!=f->halfedge());
  }
  
  occb=NULL;  
  return false;
}





//insertion from vertices may invalidate the holes containers (since a new
//face can be created which contains some of the holes of the old face.
//to validate the containers the function move_hole should be called.


//if a new face is created, returns the halfedge on the new face.  the
//new face will ALWAYS (by definition) be the one from v1->v2
//(==previous1->e1->..)  this sets a requirement on the user!!! (since
//a topological pm can't know if v1->v2 is on the new face or v2->v1)

//since there can be a number of edges from the vertices, and we can't
//know topologically which will be in the new face and which will be
//outside, this should be provided by the user - defining prev1 and
//prev2 (previous to v1,v2 resp.)  the planar map will define this
//geometrically using the one previous to the curve clock wise.
//prev1->target() is v1

template<class Dcel>
Topological_map<Dcel>::Halfedge_handle
Topological_map<Dcel>::
insert_at_vertices(Halfedge_handle previous1, Halfedge_handle previous2)  
{
  // vertices should be distinct
  CGAL_precondition(previous1->target()!=previous2->target()); 

  typename Dcel::Halfedge *prev1=&(*previous1), *prev2=&(*previous2);

  typename Dcel::Vertex* dv1=prev1->vertex();
  typename Dcel::Vertex* dv2=prev2->vertex();

  CGAL_precondition(prev1->face()==prev2->face());

  typename Dcel::Face* df=prev1->face();

  //will hold represantative halfedge of prev1/2
  typename Dcel::Halfedge *ccb1,*ccb2; 
  if (!(is_halfedge_on_outer_ccb<Dcel>(prev1,df,ccb1)) ) 
    is_halfedge_on_inner_ccb<Dcel>(prev1,df,ccb1);
  if (!(is_halfedge_on_outer_ccb<Dcel>(prev2,df,ccb2)) )
    is_halfedge_on_inner_ccb<Dcel>(prev2,df,ccb2);


  bool ccb1_is_inner = !(ccb1==df->halfedge());
  bool ccbs_equal = (ccb1 == ccb2);

  typename Dcel::Halfedge* e1=d.new_edge();
  typename Dcel::Halfedge* e2=e1->opposite();
  e1->set_vertex(dv1);
  e2->set_vertex(dv2);

  e1->set_next(prev1->next());
  e2->set_next(prev2->next());

  prev1->set_next(e2); 
  prev2->set_next(e1);

  //Cases:
  //a. ccb1!=ccb2 (at least one is inner ccb) - no face added
  if (!ccbs_equal) {
      e1->set_face(df);
      e2->set_face(df);
      
      //remove ONE inner ccb from df
      if (ccb1_is_inner) {
        //ccb1 becomes part of ccb2 
        find_and_erase_hole(ccb1,df);
      }
      else {
        //ccb2 becomes part of ccb1
        find_and_erase_hole(ccb2,df);
      }

      return TR_HI(e2);
  }

  //b. ccb1==ccb2
  //new face is created
  //the new face will ALWAYS (by definition) be the one from v1->v2 (==e2) 
  //this sets a requirement on the user!!!
  
  //the planar map  will define this geometrically 
  //(by finding the leftmost halfedge
  //(or left-down) when going from v1 to v2 (visiting v1 only once!!) and 
  //checking if it is up or down, if it is up then send v1->v2 else 
  //send v2->v1)

  typename Dcel::Face* nf = d.new_face();
  
  //b.1 the vertices are on outer ccb
  if (! ccb1_is_inner) {
    e2->set_face(nf);
    e1->set_face(df);

    nf->set_halfedge(e2); //set the represantative edge to be e2
    df->set_halfedge(e1);

    //set the face pointer for all halfedges of e2 to nf 
    ccb2=e2->next();
    while (ccb2!=e2) 
      {
        ccb2->set_face(nf);
        ccb2=ccb2->next();
      }
   }

  //b.2 the vertices are on an inner ccb
  else {

    e1->set_face(df); 
    e2->set_face(nf);
    
    nf->set_halfedge(e2);

    ccb2=e2->next(); 
    while (ccb2!=e2) 
      {
        ccb2->set_face(nf);
        ccb2=ccb2->next();
      }
    
    //remove ccb1 from holes_container (maybe it is now inside the new face)
    find_and_erase_hole(ccb1,df);
    df->add_hole(e1);

  }
  
  return TR_HI(e2);
}


template<class Dcel> 
void
Topological_map<Dcel>::
move_hole(Holes_iterator e, Face_handle f1, Face_handle f2)
{
  //move hole from df1 to df2 and set face pointer of each edge on hole to df1
  typename Dcel::Face* df1=&(*f1);
  typename Dcel::Face* df2=&(*f2);
  
  //set face pointers to df1
  typename Dcel::Halfedge* last=&(*(*(e))); 
  typename Dcel::Halfedge* first=last;
  do  {
    first->set_face(df2);
    first=first->next();
  } while (first!=last) ;     
  
  //copy to df2
  df2->add_hole(first);
  df1->erase_hole(e.current_iterator());
  
}  


//returns halfedge which is previous->next
template<class Dcel>
Topological_map<Dcel>::Halfedge_handle
Topological_map<Dcel>::
insert_from_vertex(Halfedge_handle previous)
{
  typename Dcel::Halfedge* prev=&(*previous);
  typename Dcel::Face* df=prev->face();

  typename Dcel::Vertex* dv1=prev->vertex();  
  
  typename Dcel::Vertex* dv2=d.new_vertex();   //the new vertex
  typename Dcel::Halfedge* h1=d.new_edge();

  typename Dcel::Halfedge* h2=h1->opposite();
  
  h1->set_vertex(dv1);
  h2->set_vertex(dv2);

  h1->set_face(df);
  h2->set_face(df);
  
  dv2->set_halfedge(h2);

  h2->set_next(h1);       
  h1->set_next(prev->next());
  
  prev->set_next(h2);

  return TR_HI(h2);
}


template<class Dcel>
Topological_map<Dcel>::Halfedge_handle
Topological_map<Dcel>::
insert_in_face_interior(Face_handle f)
{
  typename Dcel::Vertex* v1=d.new_vertex(); 
  typename Dcel::Vertex* v2=d.new_vertex(); 
  
  typename Dcel::Halfedge* h1 = d.new_edge();
  typename Dcel::Halfedge* h2 = h1->opposite(); 
  
  typename Dcel::Face* fp =&(*f);  

  h1->set_next(h2);
  h1->set_vertex(v1);
  h1->set_face(fp);

  h2->set_next(h1);
  h2->set_vertex(v2);
  h2->set_face(fp);
  
  v1->set_halfedge(h1);
  v2->set_halfedge(h2);
  
  fp->add_hole(h1);
  
  return TR_HI(h1); 
}      


template<class Dcel>
Topological_map<Dcel>::Halfedge_handle
Topological_map<Dcel>::
split_edge (Halfedge_handle e)
{
  typename Dcel::Halfedge* e1=&(*e); 
  typename Dcel::Halfedge* e2=e1->opposite();
  
  //find previous halfedge of e2
  typename Dcel::Halfedge* prev2=get_prev(e2);
  
  typename Dcel::Vertex* v=d.new_vertex();  
  
  typename Dcel::Halfedge* e3=d.new_edge();
  typename Dcel::Halfedge* e4=e3->opposite();
  
  v->set_halfedge(e4);  

  if (e1->next()!=e2)
    e3->set_next(e1->next());
  else
    e3->set_next(e4); 

  e4->set_vertex(v);
  e3->set_face(e1->face());
  e4->set_next(e2);
  e3->set_vertex(e1->vertex()); 
  e4->set_face(e2->face());
  
  e1->vertex()->set_halfedge(e3); //makes the incident halfedge e3 
  //(so if e1 was the incident halfedge it comes instead ) - changed 

  if (prev2!=e1)
    prev2->set_next(e4);  
  //otherwise e3->set_next(e4) - taken care above
  
  e1->set_next(e3);
  e1->set_vertex(v);

  return TR_HI(e1); 
}

template<class Dcel>
Topological_map<Dcel>::Halfedge_handle
Topological_map<Dcel>::
merge_edge (Halfedge_handle e1, Halfedge_handle e2) 
{
    //check e1->e2 and that degree(e1.target)==2 (i.e no other edge connected)
    CGAL_assertion(e1->target()==e2->source());
    CGAL_assertion(e1->target()->degree()==2);

    typename Dcel::Halfedge* de1=&(*e1);
    typename Dcel::Halfedge* de1t=de1->opposite();
    typename Dcel::Halfedge* de2=&(*e2);
    typename Dcel::Halfedge* de2t=de2->opposite();

    typename Dcel::Vertex* v=de1->vertex();
    
    typename Dcel::Face* f=de2->face();
    typename Dcel::Face* ft=de1t->face();

   
    //at the end de1 and de1t will remain and de2t,de2 will be deleted
    //check if they are a "represantative" of a hole or outer ccb of face
    if (f->halfedge()==de2) 
      f->set_halfedge(de1);
    else {
      if (find_and_erase_hole(de2,f))
        f->add_hole(de1);
    }

    if (ft->halfedge()==de2t) 
      ft->set_halfedge(de1t);
    else {
      if (find_and_erase_hole(de2t,f))  
        f->add_hole(de1t);
    }
    
    //in case de2 is representative halfedge of the target vertex  
    de2->vertex()->set_halfedge(de1); 

    if (de2->next() != de2t) {  //de2 is not a tip of antenna
      //find previous halfedge of de2t
      typename Dcel::Halfedge* prev2=get_prev(de2t);
      de1->set_next(de2->next());
      prev2->set_next(de1t);      
    }
    else {
      de1->set_next(de1t);
    }

    de1->set_vertex(de2->vertex());  

    d.delete_edge(de2);
    d.delete_vertex(v);

    return TR_HI(de1); 
}


template<class Dcel>
Topological_map<Dcel>::Face_handle
Topological_map<Dcel>::
remove_edge(Halfedge_handle e)
{
  typename Dcel::Halfedge* de1=&(*e);
  typename Dcel::Halfedge* de2=de1->opposite();
  typename Dcel::Face* df1=de1->face();
  typename Dcel::Face* df2=de2->face();

  typename Dcel::Halfedge* prev1;
  typename Dcel::Halfedge* prev2;

  //CASES:

  if (df1==df2) { //case a. antenna - no face deleted

    if (de1->next()==de2 && de2->next()==de1) {
      //a.1 only one edge as a hole
      if (!find_and_erase_hole(de1,df1))
        find_and_erase_hole(de2,df1);
      d.delete_vertex(de2->vertex());
      d.delete_vertex(de1->vertex());
      d.delete_edge(de1);

      return TR_FI(df1); 
    }
    
    
    if (de1->next()==de2 || de2->next()==de1) {
      //a.2 a tip of the antenna
      if (de1->next()==de2) {  //de1 points at the tip
        prev1=get_prev(de1);
        prev1->set_next(de2->next());
        
        //check if they are a represantative halfedge if so exchange
        if (df1->halfedge()==de2) 
          df1->set_halfedge(prev1);
        else {
          if (find_and_erase_hole(de2,df1))
            df1->add_hole(prev1);
        }
    
        if (df1->halfedge()==de1) 
          df1->set_halfedge(prev1);
        else {
          if (find_and_erase_hole(de1,df1))  
            df1->add_hole(prev1);
        }
        
        //check if de2 is incident halfedge of de2->vertex and if so - 
        //change to prev1 (iddo 16/6)
        if (de2->vertex()->halfedge() == de2)
          de2->vertex()->set_halfedge(prev1);

        d.delete_vertex(de1->vertex());  //changed to conform with HDS
     }
      else { //de2 points at the tip - find previous halfedge of de2 
        prev2=get_prev(de2);
        prev2->set_next(de1->next());

       //check if they are represantative halfedge if so exchange
        if (df1->halfedge()==de2) 
          df1->set_halfedge(prev2);
        else {
          if (find_and_erase_hole(de2,df1))
            df1->add_hole(prev2);
        }
    
        if (df1->halfedge()==de1) 
          df1->set_halfedge(prev2);
        else {
          if (find_and_erase_hole(de1,df1))  
            df1->add_hole(prev2);
        }

        //check if de1 is incident halfedge of de1->vertex and if so - 
        //change to prev2 
        if (de1->vertex()->halfedge() == de1)
          de1->vertex()->set_halfedge(prev2) ;

        d.delete_vertex(de2->vertex());  
      }
      
      d.delete_edge(de1);

      return TR_FI(df1); 
    }      
    
    //a.3 a middle of an antenna - a new hole is created
    prev1=get_prev(de1);
    prev2=get_prev(de2);

    //find the ccb of de1/de2        
    typename Dcel::Halfedge *ccb1;
    if (is_halfedge_on_outer_ccb<Dcel>(prev1,df1,ccb1)) { 
      //antenna is connected to outer ccb, split a hole from it.
      //we assume in this case that de1->next is the hole 
      //(to be determined genoetrically by the user) 
      
      df1->add_hole(de1->next());
      df1->set_halfedge(prev1);
    }

    else { //antenna is a hole - split it into two holes
      is_halfedge_on_inner_ccb<Dcel>(prev1,df1,ccb1);
#ifndef CGAL_NO_ASSERTIONS // in order to avoid warnings
      bool hole_found = 
#endif
	find_and_erase_hole(ccb1,df1); 
      CGAL_assertion(hole_found) ;//ccb1 must be a hole in df1
      df1->add_hole(prev1);
      df1->add_hole(prev2);
    }
    
    prev1->set_next(de2->next());
    prev2->set_next(de1->next());    
    
    //check if de1 is incident halfedge of de1->vertex() and if so change to
    //prev2 same for de2
    if (de1->vertex()->halfedge() == de1)
      de1->vertex()->set_halfedge(prev2);

    if (de2->vertex()->halfedge() == de2)
      de2->vertex()->set_halfedge(prev1);

    d.delete_edge(de1);

    return TR_FI(df1); 
  }

  else {    //case b. edge between faces - merge faces
    typename Dcel::Halfedge* ccb;
    bool de1_is_on_hole=!(is_halfedge_on_outer_ccb<Dcel>(de1,df1,ccb));
    bool de2_is_on_hole=!(is_halfedge_on_outer_ccb<Dcel>(de2,df2,ccb));
    if ( !de1_is_on_hole && !de2_is_on_hole ) {
      //b.1 both halfedges are on outer boundary
      prev1=get_prev(de1);
      prev2=get_prev(de2);
      
      //arbitrarily df2 will be removed - change all edges so point to df1
      //set the face pointer for all halfedges of df2 to df1
      typename Dcel::Halfedge* aux=de2->next();
      while (aux!=de2) 
        {
          aux->set_face(df1);
          aux=aux->next();
        }
      
      //move holes from df2 to df1 and set face pointer of each edge on 
      //hole to df1  
      typename Dcel::Face::Holes_iterator hi = df2->holes_begin();
      for ( ; hi!=df2->holes_end(); ++hi) {
        //set face pointers to df1
        typename Dcel::Halfedge* aux1=(*hi);
        do  {
          aux1->set_face(df1);
          aux1=aux1->next();
        } while (aux1!=(*hi)) ;     
      }      
      //copy the holes to the holes list of df1
      for (hi = df2->holes_begin(); hi!=df2->holes_end(); ++hi) {
        df1->add_hole(*hi);
      }
      df2->erase_holes(df2->holes_begin(),df2->holes_end());

      //check if represantative edge
      if (de1==df1->halfedge()) 
        df1->set_halfedge(prev1);

      //check if de1 is incident halfedge of de1->vertex() and if so change to
      //prev2 same for de2 
      if (de1->vertex()->halfedge() == de1)
        de1->vertex()->set_halfedge(prev2);
      
      if (de2->vertex()->halfedge() == de2)
        de2->vertex()->set_halfedge(prev1);
      
      
      prev1->set_next(de2->next());
      prev2->set_next(de1->next());
      
      d.delete_face(df2);
      d.delete_edge(de1);
      
      return TR_FI(df1); 
    }
    else { //edge is on a hole
        prev1=get_prev(de1);
        prev2=get_prev(de2);
 
      if (de1_is_on_hole) { //df1 is the outer face df2 was the hole
         //df2 will be removed - change all edges so point to df1
        //set the face pointer for all halfedges of df2 to df1
        typename Dcel::Halfedge* aux=de2->next();
        while (aux!=de2) {
          aux->set_face(df1);
          aux=aux->next();
        }
        
        //move holes from df2 to df1 and set face pointer of each edge on 
	//hole to df1

      typename Dcel::Face::Holes_iterator hi = df2->holes_begin();
      for ( ; hi!=df2->holes_end(); ++hi) {
        //set face pointers to df1
        typename Dcel::Halfedge* aux1=(*hi);
        do  {
          aux1->set_face(df1);
          aux1=aux1->next();
        } while (aux1!=(*hi)) ;     
      }      
      //copy to the holes list of df1
      for (hi = df2->holes_begin(); hi!=df2->holes_end(); ++hi) {
        df1->add_hole(*hi);
      }
      df2->erase_holes(df2->holes_begin(),df2->holes_end());

      
      //check if de1 is a represantative hole of df1 (it is not on the outer) 
      if (find_and_erase_hole(de1,df1))
          df1->add_hole(prev1);

      //check if de1 is incident halfedge of de1->vertex() and if so change to
      //prev2 same for de2 
      if (de1->vertex()->halfedge() == de1)
        de1->vertex()->set_halfedge(prev2);
      
      if (de2->vertex()->halfedge() == de2)
        de2->vertex()->set_halfedge(prev1);
      


        prev1->set_next(de2->next());
        prev2->set_next(de1->next());
        
        d.delete_face(df2);
        d.delete_edge(de1);
        
        return TR_FI(df1);
      }
      else { //df2 is the outer face
        //erase df1
        //df1 will be removed - change all edges so point to df2
        typename Dcel::Halfedge* aux=de1->next();
        while (aux!=de1) {
          aux->set_face(df2);
          aux=aux->next();
        }
        
        //move holes from df1 to df2 and set face pointer of each edge on 
	//hole to df1  

      typename Dcel::Face::Holes_iterator hi = df1->holes_begin();
      for ( ; hi!=df1->holes_end(); ++hi) {
        //set face pointers to df2
        typename Dcel::Halfedge* aux1=(*hi);
        do  {
          aux1->set_face(df2);
          aux1=aux1->next();
        } while (aux1!=(*hi)) ;     
      }      
      //copy to the holes list of df2
      for (hi = df1->holes_begin(); hi!=df1->holes_end(); ++hi) {
        df2->add_hole(*hi);
      }
      df1->erase_holes(df1->holes_begin(),df1->holes_end());

      //check if de2 is a represantative hole of df2 
      if (find_and_erase_hole(de2,df2))
        df2->add_hole(prev2);
      
      //check if de1 is incident halfedge of de1->vertex() and if so change to
      //to prev2 same for de2 
      if (de1->vertex()->halfedge() == de1)
        de1->vertex()->set_halfedge(prev2);
        
      if (de2->vertex()->halfedge() == de2)
        de2->vertex()->set_halfedge(prev1);
      
      prev1->set_next(de2->next());
      prev2->set_next(de1->next());
        
      d.delete_face(df1);
      d.delete_edge(de1);
      
      return TR_FI(df2); 

      } //df2 is outer..
    } //edge is on a hole..
  } //case b. (face deleted)
}

 

CGAL_END_NAMESPACE

#else   
#error  Header file .h included twice
#endif  

/*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
 *     
 * 
\*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*/









