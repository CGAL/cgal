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
// release       : $CGAL_Revision: CGAL-2.3-I-44 $
// release_date  : $CGAL_Date: 2001/03/09 $
//
// file          : include/CGAL/Planar_map_2.h
// package       : pm (5.45)
// maintainer    : Eyal Flato <flato@math.tau.ac.il>
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Iddo Hanniel <hanniel@math.tau.ac.il>
//                 Eyal Flato
//                 Oren Nechushtan <theoren@math.tau.ac.il>
//                 Eti Ezra <estere@post.tau.ac.il>
//
//
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter       : 
// ======================================================================
#ifndef CGAL_PLANAR_MAP_2_H
#define CGAL_PLANAR_MAP_2_H

#ifndef CGAL_PLANAR_MAP_MISC_H
#include <CGAL/Planar_map_2/Planar_map_misc.h>
#endif

#ifndef CGAL_TOPOLOGICAL_MAP_H
#include <CGAL/Topological_map.h>
#endif

#ifndef CGAL_NO_PM_DEFAULT_POINT_LOCATION

#ifndef CGAL_PM_DEFAULT_POINT_LOCATION_H
#include <CGAL/Pm_default_point_location.h>
#endif

#ifndef CGAL_PM_WALK_ALONG_LINE_POINT_LOCATION_H
#include <CGAL/Pm_walk_along_line_point_location.h>
#endif

#ifndef CGAL_PM_NAIVE_POINT_LOCATION_H
#include <CGAL/Pm_naive_point_location.h>
#endif
//#else // CGAL_NO_PM_DEFAULT_POINT_LOCATION
#ifndef CGAL_PM_POINT_LOCATION_BASE_H
#include <CGAL/Pm_point_location_base.h>
#endif
#endif // CGAL_NO_PM_DEFAULT_POINT_LOCATION

// for solving the dynamic cast in the copy constructor, 
// these lines will be removed after writing 
// copy construtor for point location.
#ifndef CGAL_PM_WALK_ALONG_LINE_POINT_LOCATION_H
#include <CGAL/Pm_walk_along_line_point_location.h>
#endif

#ifndef CGAL_PM_NAIVE_POINT_LOCATION_H
#include <CGAL/Pm_naive_point_location.h>
#endif
// end.

/*#ifndef CGAL_NO_PM_DEFAULT_POINT_LOCATION
  #ifndef CGAL_PM_DEFAULT_POINT_LOCATION_H
  #include <CGAL/Pm_default_point_location.h>
  #endif
  #else // CGAL_NO_PM_DEFAULT_POINT_LOCATION
  #ifndef CGAL_PM_POINT_LOCATION_BASE_H
  #include <CGAL/Pm_point_location_base.h>
#endif*/

//#endif // CGAL_NO_PM_DEFAULT_POINT_LOCATION

// default bounding box for finite curves
#ifndef CGAL_PM_UNBOUNDING_BOX_H
#include <CGAL/Pm_unbounding_box.h>
#endif

// default bounding box for infinite curves
#ifndef CGAL_PM_DYNAMIC_CLOSED_BOUNDING_BOX_H
#include <CGAL/Pm_dynamic_open_bounding_box.h>
#endif

#ifndef CGAL_IO_PM_FILE_SCANNER_H
#include <CGAL/IO/Pm_file_scanner.h>
#endif // CGAL_IO_PM_FILE_SCANNER_H

#include <list>

CGAL_BEGIN_NAMESPACE

////////////////////////////////////////////////////////////////////////////
//      PLANAR_MAP_2   

template <class Dcel_, class Traits_> 
class Planar_map_2 : public Topological_map<Dcel_>{
public:
  typedef Dcel_ Dcel;
  typedef Traits_ Traits;
  typedef Planar_map_2<Dcel,Traits> Self;
  typedef Planar_map_traits_wrap<Traits> Traits_wrap;
  typedef typename Traits::X_curve X_curve;
  typedef typename Traits::Point Point;
  typedef std::list<X_curve> X_curve_container;
  
  typedef Topological_map<Dcel> TPM;
  typedef typename TPM::Vertex_iterator Vertex_iterator;
  typedef typename TPM::Halfedge_iterator Halfedge_iterator;
  typedef typename TPM::Face_iterator Face_iterator;
  typedef typename TPM::Vertex_const_iterator Vertex_const_iterator;
  typedef typename TPM::Halfedge_const_iterator Halfedge_const_iterator;
  typedef typename TPM::Face_const_iterator Face_const_iterator;
  typedef typename TPM::Vertex_handle Vertex_handle;
  typedef typename TPM::Vertex_const_handle Vertex_const_handle;
  typedef typename TPM::Halfedge_handle Halfedge_handle;
  typedef typename TPM::Face_handle Face_handle;
  typedef typename TPM::Halfedge_const_handle Halfedge_const_handle;
  typedef typename TPM::Face_const_handle Face_const_handle;
  typedef typename TPM::Halfedge_around_vertex_circulator
  Halfedge_around_vertex_circulator;
  typedef typename TPM::Halfedge_around_vertex_const_circulator 
  Halfedge_around_vertex_const_circulator;
  typedef typename TPM::Holes_iterator Holes_iterator;
  typedef typename TPM::Holes_const_iterator Holes_const_iterator;
  typedef typename TPM::Ccb_halfedge_const_circulator 
  Ccb_halfedge_const_circulator;
  typedef typename TPM::Ccb_halfedge_circulator Ccb_halfedge_circulator;
  
  typedef typename TPM::Size Size;
  
  typedef Pm_point_location_base<Self> Point_location_base;
  typedef Pm_bounding_box_base<Self> Bounding_box_base;
  
  typedef enum{ 
    VERTEX = 1, 
      EDGE, 
      FACE , 
      UNBOUNDED_VERTEX, 
      UNBOUNDED_EDGE, 
      UNBOUNDED_FACE } Locate_type ;
  
#ifndef CGAL_NO_PM_DEFAULT_POINT_LOCATION
  // constructor #1 - no parameters
  Planar_map_2 ()
  {
    traits = new Traits_wrap();
    use_delete_traits = true;
    
    pl = new Pm_default_point_location<Self>;
    use_delete_pl = true;
    pl->init(*this,*traits);
    
    bb=init_default_bounding_box((Traits*)traits);
    use_delete_bb=true;
    bb->init(*this,*traits);
  }    
  
#endif

  // constructor #2 - set only the PL
  Planar_map_2(Point_location_base *pl_ptr)
  {
    traits = new Traits_wrap();
    use_delete_traits = true;
    
    pl = pl_ptr;
    use_delete_pl = false;
    pl->init(*this,*traits);
    
    bb = init_default_bounding_box((Traits*)traits);
    use_delete_bb = true;
    bb->init(*this,*traits);
  }
  
  // constructor #3 - copy traits, set pl and bb
  // set NULLs for defaults
  Planar_map_2(const Traits &tr_, Point_location_base *pl_ptr, 
	       Bounding_box_base *bb_ptr)
  {
    traits = new Traits_wrap(tr_);
    use_delete_traits = true;
    
    if (pl_ptr == NULL)
      {
#ifndef CGAL_NO_PM_DEFAULT_POINT_LOCATION
        pl = new Pm_default_point_location<Self>;
        use_delete_pl = true;
        pl->init(*this,*traits);
#else
        assert(0); 
	// if no default PL is defined you must supply a pl.
#endif
      }
    else
      {
        pl = pl_ptr;
        use_delete_pl = false;
        pl->init(*this,*traits);
      }
    
    if (bb_ptr == NULL)
      {
        bb=init_default_bounding_box((Traits*)traits);
        use_delete_bb=true;
        bb->init(*this,*traits);
      }
    else
      {
        bb = bb_ptr;
        use_delete_bb = false;
        bb->init(*this,*traits);
      }
  }
  
  // constructor #4 - set traits, pl and bb
  // set NULLs for defaults
  Planar_map_2(Traits_wrap *tr_ptr, Point_location_base *pl_ptr, 
	       Bounding_box_base *bb_ptr)
  {
    if (tr_ptr == NULL)
      {
        traits = new Traits_wrap();
        use_delete_traits = true;
      }
    else
      {
        traits = tr_ptr;
        use_delete_traits = false;
      }
    
    if (pl_ptr == NULL)
      {
#ifndef CGAL_NO_PM_DEFAULT_POINT_LOCATION
        pl = new Pm_default_point_location<Self>;
        use_delete_pl = true;
        pl->init(*this,*traits);
#else
        assert(0); 
        // if no default PL is defined you must supply a pl.
#endif
      }
    else
      {
        pl = pl_ptr;
        use_delete_pl = false;
        pl->init(*this,*traits);
      }
    
    if (bb_ptr == NULL)
      {
        bb=init_default_bounding_box((Traits*)traits);
        use_delete_bb=true;
        bb->init(*this,*traits);
      }
    else
      {
        bb = bb_ptr;
        use_delete_bb = false;
        bb->init(*this,*traits);
      }
  }

///////////////////////////////////////////////////////////////////////////////
//                                 Copy constructor.
///////////////////////////////////////////////////////////////////////////////
  Planar_map_2(const Self& pm){
    // doing the same as Planar_map_2(pm.get_traits(),pm.get_point_location(),
    //                                pm.get_point_bbox());
    
    typedef Pm_naive_point_location<Planar_map_2<Dcel,Traits> >  Pm_naive;
    typedef Pm_naive*                                    Pm_naive_pointer;

    traits = new Traits_wrap();
    use_delete_traits = true;

    if (Pm_naive_pointer tmp_pl = dynamic_cast<Pm_naive_pointer>(pm.pl) ){
      //cout<<"Naive"<<std::endl;
      pl = new Pm_naive_point_location<Self>;
    }
    else if (Pm_walk_along_line_point_location<Self>* tmp_pl = 
	     dynamic_cast<Pm_walk_along_line_point_location<Self>*>(pm.pl) ){
      pl = new Pm_walk_along_line_point_location<Self>;
      //cout<<"Walk"<<std::endl;
    }
    else{
      //cout<<"Default"<<std::endl;
#ifndef CGAL_NO_PM_DEFAULT_POINT_LOCATION
      pl = new Pm_default_point_location<Self>;
#else
      assert(0); 
      // if no default PL is defined you must supply a pl.
#endif
    }
    use_delete_pl = true;
    pl->init(*this,*traits);
    
    bb=init_default_bounding_box((Traits*)traits);
    use_delete_bb=true;
    bb->init(*this,*traits);
    
    assign(pm);
    
    Halfedge_iterator h_iter;
    for (h_iter = halfedges_begin(); 
	 h_iter != halfedges_end(); 
	 h_iter++, h_iter++)
      pl->insert(h_iter, h_iter->curve());
    
    for (Vertex_iterator v_iter = vertices_begin(); 
	 v_iter != vertices_end(); 
	 v_iter++)
      bb->insert(v_iter->point());
    
    for (h_iter = halfedges_begin(); 
	 h_iter !=  halfedges_end(); 
	 h_iter++, h_iter++)
      bb->insert(h_iter->curve());
  }
  
  /////////////////////////////////////////////////////////////////////////////
  //                  Reading Planar map functions. 
  ////////////////////////////////////////////////////////////////////////////
  bool read (std::istream &in)
  {
    clear();
    
    Pm_file_scanner<Self>  scanner(in); 
    
    return scan_planar_map(scanner);
  }    

  template <class Scanner>
  bool read (std::istream &in, Scanner& scanner)
  {
    clear(); 
    
    return scan_planar_map(scanner);
  } 
  
  /*
#ifndef CGAL_NO_PM_DEFAULT_POINT_LOCATION
    
  Planar_map_2 (const Traits& tr_=Traits()) :
    pl(new Pm_default_point_location<Self>),use_delete_pl(true)
  {
    traits = new Traits_wrap(tr_);
    use_delete_traits = true;
    pl->init(*this,*traits);
    bb=init_default_bounding_box(traits);
    use_delete_bb=true;
    bb->init(*this,*traits);
  }    
  
  Planar_map_2( Bounding_box_base *bb_ptr, const Traits& tr_): 
    pl(new Pm_default_point_location<Self>),bb(bb_ptr),
    use_delete_pl(true),use_delete_bb(false)
  {
    traits = new Traits_wrap(tr_);
    use_delete_traits = true;
    if (!bb) 
      {				
        bb=init_default_bounding_box(&tr_);
        use_delete_bb=true;
      }
    pl->init(*this,*traits);
    bb->init(*this,*traits);
  }    
  
  Planar_map_2( Bounding_box_base *bb_ptr): 
    pl(new Pm_default_point_location<Self>),bb(bb_ptr), 
    traits(new Traits_wrap),
    use_delete_pl(true),use_delete_bb(false),use_delete_traits(true)
  {
    if (!bb) 
      {				
        bb=init_default_bounding_box((Traits*)0);
        use_delete_bb=true;
      }
    pl->init(*this,*traits);
    bb->init(*this,*traits);
  }    
#endif // CGAL_NO_PM_DEFAULT_POINT_LOCATION
  
  Planar_map_2 (Point_location_base *pl_ptr, Traits_wrap *tr_ptr) : 
    pl(pl_ptr),use_delete_pl(false), traits(tr_ptr),use_delete_traits(false)
  {
    pl->init(*this,*traits);
    bb=init_default_bounding_box(traits);
    use_delete_bb=true;
    bb->init(*this,*traits);
  }
  
  Planar_map_2(Point_location_base *pl_ptr,const Traits& tr_) : 
    pl(pl_ptr), use_delete_pl(false)
  {
    traits = new Traits_wrap(tr_);
    use_delete_traits = true;
    bb = init_default_bounding_box(traits);
    use_delete_bb = true;
    // initialize the bounding box.
    pl->init(*this,*traits);
    bb->init(*this,*traits);
  }
  
  Planar_map_2(Point_location_base *pl_ptr,Bounding_box_base *bb_ptr, 
               const Traits& tr_) : 
    pl(pl_ptr),bb(bb_ptr),
    use_delete_pl(false),use_delete_bb(false),
  {
    traits = new Traits_wrap(tr_);
    use_delete_traits = true;
    pl->init(*this,traits);
    bb->init(*this,traits);
  }
  
  Planar_map_2(Point_location_base *pl_ptr):
    pl(pl_ptr),bb(init_default_bounding_box((Traits*)0)),
    traits(new Traits_wrap),
    use_delete_pl(false),use_delete_bb(true),use_delete_traits(true)
  {
    if (!traits) 
      {				
        traits=new Traits_wrap;
        use_delete_traits=true;
      }
    // initialize the bounding box.
    pl->init(*this,*traits);
    bb->init(*this,*traits);
  }
  
  Planar_map_2(
               Point_location_base *pl_ptr,Bounding_box_base *bb_ptr):
    pl(pl_ptr),bb(bb_ptr),
    traits(new Traits_wrap),			
    use_delete_pl(false),
    use_delete_bb(false),
    use_delete_traits(true)
  {
    pl->init(*this,*traits);
    bb->init(*this,*traits);
  }
*/
  
  virtual ~Planar_map_2 () {
    if (use_delete_pl) delete pl;
    if (use_delete_bb) delete bb;
    if (use_delete_traits) delete traits;
  }

  // Inserts a new curve cv in the interior of the face f.
  // Returns the halfedge which is directed in the same way as the curve cv
  // (i.e., traits->curve_source(cv) == h.source()).
  Halfedge_handle insert_in_face_interior(const X_curve& cv, Face_handle f) 
  {
    Halfedge_handle h = Topological_map<Dcel>::insert_in_face_interior(f);
    h->set_curve(cv);  //should set the curve of the twin as well but for now
    h->twin()->set_curve(cv);
			
    //pl->insert(h);  //maybe should be above
    //iddo - for arrangement
    pl->insert(h,cv);

    h->source()->set_point(traits->curve_source(cv));
    h->target()->set_point(traits->curve_target(cv));
			
    return h;
  }
		
  Halfedge_handle insert_from_vertex(const X_curve& cv, 
				     Vertex_handle v1, bool source) 
  {
    //find the previous of cv.
    Halfedge_around_vertex_circulator
      previous=v1->incident_halfedges(),
      after=previous,
      infinite_loop=previous;
    ++after;
			
    if (after!=previous) {
      while (!(traits->curve_is_between_cw(cv,previous->curve(),
					   after->curve(),v1->point()))) {
	previous=after;
	++after;
	if (previous==infinite_loop)  // infinite loop indication
	  {
	    std::cerr << std::endl << "Planar_map_2::insert_from_vertex(" <<
	      "const X_curve& cv, Vertex_handle v1, " << 
	      "bool source) called with previously " <<
	      "inserted curve " << std::endl;
	    return Halfedge_handle();
	  }
      }
    }
			
    Halfedge_handle h = Topological_map<Dcel>::insert_from_vertex(previous);  
    h->set_curve(cv);  
    h->twin()->set_curve(cv);
			
    //pl->insert(h);  //maybe should be above
    //iddo - for arrangement
    pl->insert(h,cv);
			
    //h is now pointing from v1
    if (source)
      h->target()->set_point(traits->curve_target(cv));
    else
      h->target()->set_point(traits->curve_source(cv));

    return h;
  }
		
  Halfedge_handle insert_at_vertices(const X_curve& cv, 
				     Vertex_handle v1, Vertex_handle v2) 
  {
    Size num_before=number_of_faces();
			
    Halfedge_around_vertex_circulator 
      previous1=v1->incident_halfedges(),
      previous2=v2->incident_halfedges(),
      after=previous1,
      infinite_loop=previous1;
    ++after;
			
    if (after!=previous1) {
      while (!(traits->curve_is_between_cw(cv,previous1->curve(),
					   after->curve(),v1->point()))) {
	previous1=after;
	++after;
	if (previous1==infinite_loop)  // infinite loop indication
	  {
	    std::cerr << std::endl << "Planar_map_2::insert_at_vertices(" <<
	      "const X_curve& cv, Vertex_const_handle v1, " <<
	      "Vertex_const_handle v2) called with previously " <<
	      "inserted curve " << std::endl;
	    return Halfedge_handle();
	  }
      }
    }    
			
    after=previous2;
    infinite_loop=previous2;
    ++after;
			
    if (after!=previous2) {
      while (!(traits->curve_is_between_cw(cv,previous2->curve(),
					   after->curve(),v2->point()))) {
	previous2=after;
	++after;
	if (previous2==infinite_loop) // infinite loop indication
	  {
	    std::cerr << std::endl << "Planar_map_2::insert_at_vertices(" <<
	      "const X_curve& cv, Vertex_const_handle v1," <<
	      "Vertex_const_handle v2) called with previously " <<
	      "inserted curve " << std::endl;
	    return Halfedge_handle();
	  }
      }
    }    
			
			
    bool prev1_before_prev2 = prev1_inside_hole(previous1,previous2,cv);
    Halfedge_handle h;
    if (prev1_before_prev2)
      h = Topological_map<Dcel>::insert_at_vertices(previous1,previous2); 
    else
      h = Topological_map<Dcel>::insert_at_vertices(previous2,previous1);
			
    h->set_curve(cv); 
    h->twin()->set_curve(cv);
			
    Size num_after=number_of_faces();
    if (num_after-num_before) { //if additional face was added - move holes
      Face_handle nf = h->face(); //the new face will always be the one 
                                  // pointed at by h
      Face_handle of = h->twin()->face(); //old face
				
      Holes_iterator it=of->holes_begin();
      while (it!=of->holes_end()) {
					
	//check if the hole is inside new face
	//        if ( point_is_in((*it)->target()->point(),nf) ) {
	//new for arrangement
	if ( point_is_in((*it)->target()->point(),h,cv) ) {
	  Holes_iterator tmp=it;  //deletion invalidates iterators so... 
	  ++it;   //assumes only the erased iterator is invalidated (like stl
	  //list) 
						
	  move_hole( tmp,of,nf); 
	}
	else
	  ++it;
      }
				
    }
			
    if (!prev1_before_prev2) h=h->twin();
			
			
    //pl->insert(h);
    //iddo - for arrangement
    pl->insert(h,cv);
			
    return h;
  }   
		
  bool is_empty() const {return halfedges_begin()==halfedges_end();}
  
  // Note that the return types are abstract base classes
  const Point_location_base* get_point_location() const {return pl;}
  const Bounding_box_base* get_bounding_box() const {return bb;}
  
  const Traits_wrap& get_traits() { return *traits;}

		
private:
	
  //a private implementation which defines if previous1 is on an outer ccb of 
  //the new face (returns true) or on an inner ccb (returns false)
  bool prev1_inside_hole(Halfedge_const_handle previous1,
			 Halfedge_const_handle previous2,
			 const X_curve& cv) 
  {
    /*
      Defining geometrically whether there is a new face an if 
      there is find if previous1 is on the outside of the new 
      face (send previous1,previous2) or on the inside of the new 
      face (send previous2,previous1)
		
      The algorithm: 
      1. go over all the halfedges of the face which will hold 
      previous1 (since the new face is not constructed yet, 
      this is modeled by going from previous2->next to previous1 and then
      over the new curve)
		
      2. find if the left-most-lower 
      halfedge  in the path (i.e, the one with the leftmost down target and 
      is the lowest to the right among the incident edges of this vertex)
      is directed left (we are on the outside) or right (we are inside )
      (if not on same ccb then it doesn't matter and return true)
    */		
    Ccb_halfedge_const_circulator left_edge(previous2);
    ++left_edge;
    Ccb_halfedge_const_circulator first(previous2),curr(left_edge),
      last(previous1);
    ++last; //we want the previous1 to be checked as well 
		
    Point left = previous2->target()->point();
    bool b;
		
    do {
      //source
      b=false;
      if (traits->point_is_left( curr->source()->point(),left)) 
	b=true;
      else
	if (traits->point_is_same(curr->source()->point(),left)) {
	  if (traits->curve_is_vertical(curr->curve()) &&
	      traits->point_is_lower(curr->target()->point(),left) )
	    b=true;
	  else
	    if (traits->curve_compare_at_x_right(curr->curve(),
						 left_edge->curve(),
						 left)==SMALLER ) 
	      b=true;
	}
				
      if (b) {
	left=curr->source()->point();
	left_edge=curr;
      }
				
				//target
      b=false;
      if (traits->point_is_left( curr->target()->point(),left))
	b=true;
      if (traits->point_is_same(curr->target()->point(),left)) {
	if (traits->curve_is_vertical(curr->curve()) &&
	    traits->point_is_lower(curr->source()->point(),left) )
	  b=true;
	else
	  if (traits->curve_compare_at_x_right(curr->curve(),
					       left_edge->curve(),
					       left)==SMALLER ) 
	    b=true;
						
	//we want in the degenerate case to return the halfedge 
	//pointing _at_ the left point 
	  else
	    if ( (curr)==(left_edge->twin()) )
	      b=true;
							
      }
				
      if (b) {
	left=curr->target()->point();
	left_edge=curr;
      }
				
      ++curr;
    } while ( (curr != first) && (curr != last) );
		
		
    //test the new curve against left_edge
    if (traits->point_is_same(traits->curve_target(cv),left)||
	traits->point_is_same(traits->curve_source(cv),left)) {
      if (traits->curve_is_vertical(cv)) {
	return (traits->point_is_lower(previous2->target()->point(),
				       previous1->target()->point()));
      }
      else
	if (traits->curve_compare_at_x_right(cv,left_edge->curve(),
					     left)==SMALLER ) {  
	  return (traits->point_is_left(previous1->target()->point(),
					previous2->target()->point()));
	}
    }
		
    //check if left_edge is from left to right
    if (traits->curve_is_vertical(left_edge->curve())) {
      if (traits->point_is_lower(left_edge->source()->point(),
				 left_edge->target()->point()))
	return false;
      else
	return true;
    }
		
    return (traits->point_is_left(left_edge->source()->point(),
				  left_edge->target()->point()));
		
  }
  
  
public:
  Halfedge_handle insert(const X_curve& cv) {
    CGAL_assertion(bb);
    bb->insert(cv);
    
    if (traits->curve_is_degenerate(cv)) 
      {
        std::cerr << "\nPlanar_map_2::insert(const X_curve& cv) " << 
          "called with a null length curve " << cv << std::flush;
        return Halfedge_handle();
      }
    
    Locate_type lt1,lt2;
    Point p1=traits->curve_source(cv);
    Point p2=traits->curve_target(cv);
    
    // The point location may not change the bounding box.
    Halfedge_handle h1=((const Point_location_base*)pl)->locate(p1,lt1);
    Halfedge_handle h2=((const Point_location_base*)pl)->locate(p2,lt2);
    
    if (lt1==EDGE || lt1==UNBOUNDED_EDGE) 
      // the curve intersects the bounding box.
      {
        Halfedge_handle h=h1,h2;
        bb->split_boundary_edge(h,h1,h2,p1);
        // make sure the intersection point is in the map, 
        // i.e. split the halfedge that contains its.
        lt1=VERTEX; 
      }
    if (lt2==EDGE || lt2==UNBOUNDED_EDGE) 
      {
        Halfedge_handle h1,h=h2;
        bb->split_boundary_edge(h,h1,h2,p2);
        // make sure the intersection point is in the map, 
        // i.e. split the halfedge that contains its.
        lt1=VERTEX; 
      }		
    if (lt1==VERTEX && lt2==VERTEX) 
      return insert_at_vertices(cv,h1->target(),h2->target()); 
    
    if (lt1==VERTEX && lt2!=VERTEX)
      return insert_from_vertex(cv,h1->target(),true); 
    if (lt1!=VERTEX && lt2==VERTEX)
      return insert_from_vertex(cv,h2->target(),false)->twin();
    
    if (lt1==UNBOUNDED_FACE)
      return insert_in_face_interior(cv,unbounded_face());
    
    if (lt1==FACE)
      return insert_in_face_interior(cv,h1->face());
    
    CGAL_postcondition(lt1==VERTEX||lt1==UNBOUNDED_FACE||lt1==FACE);
    return Halfedge_handle();
    
  }
  
  template <class X_curve_iterator>
  Halfedge_iterator insert(const X_curve_iterator& begin,
			   const X_curve_iterator& end) {
    X_curve_iterator it=begin;
    Halfedge_iterator out;
    if (it!=end) 
      {
        out=insert(*it);
        it++;
      }
    while (it!=end)
      {
        insert(*it);
        it++;
      }
    return out;
  }
  
  Halfedge_handle split_edge(Halfedge_handle e, 
                             const X_curve& c1, const X_curve& c2)
  {
    
    CGAL_precondition(traits->point_is_same(traits->curve_source(c2),
					    traits->curve_target(c1)));
		
    CGAL_precondition(
		      traits->point_is_same(
					    traits->curve_source(c1),
					    e->source()->point()) &&
		      traits->point_is_same(	
					    traits->curve_target(c2),
					    e->target()->point()) ||
		      traits->point_is_same(
					    traits->curve_source(c1),
					    e->target()->point()) &&
		      traits->point_is_same(
					    traits->curve_target(c2),
					    e->source()->point()) );
		
    X_curve cv(e->curve());
		
    Halfedge_handle h = Topological_map<Dcel>::split_edge(e);
		
    if (traits->point_is_same(traits->curve_source(c1),h->source()->point())) {
      h->set_curve(c1);
      h->twin()->set_curve(c1);
      h->next_halfedge()->set_curve(c2);
      h->next_halfedge()->twin()->set_curve(c2);
      h->target()->set_point(traits->curve_target(c1));
      pl->split_edge(cv,h,h->next_halfedge(),c1,c2);
    }
    else {
      h->set_curve(c2);
      h->twin()->set_curve(c2);
      h->next_halfedge()->set_curve(c1);
      h->next_halfedge()->twin()->set_curve(c1);
      h->target()->set_point(traits->curve_target(c1));
      pl->split_edge(cv,h,h->next_halfedge(),c2,c1);
    }
		
    return h;
  }
	
	
  Halfedge_handle merge_edge(Halfedge_handle e1, Halfedge_handle e2, 
			     const X_curve& cv) {
    CGAL_precondition( (traits->point_is_same(traits->curve_source(cv),
					      e1->source()->point() )&&
			traits->point_is_same(traits->curve_target(cv),
					      e2->target()->point())) || 
		       (traits->point_is_same(traits->curve_target(cv),
					      e1->source()->point() )&&
			traits->point_is_same(traits->curve_source(cv),
					      e2->target()->point())) );
		
    X_curve c1(e1->curve()), c2(e2->curve());
		
    Halfedge_handle h = Topological_map<Dcel>::merge_edge(e1,e2); 
    h->set_curve(cv);
    h->twin()->set_curve(cv);
		
    //pl->merge_edge(c1,c2,h);
    //iddo - for arrangement
    pl->merge_edge(c1,c2,h,cv);
		
		
    return h;
  }
	
  Face_handle remove_edge(Halfedge_handle e) {
    pl->remove_edge(e);
		
    //if a new hole can be created define geometrically the 
    //halfedge (e or e->twin) that points at the new hole.
    //if the leftmost point in the path e...e->twin
    //is left of the leftmost point in the path e->twin ... e
    //then e->twin  points at the hole created.
		
    if (e->face() == e->twin()->face() ) {
      Ccb_halfedge_circulator ccb_e=e->ccb() ;
      Ccb_halfedge_circulator ccb_t=e->twin()->ccb();
			
      Point e_left=e->target()->point();
      Point t_left=ccb_t->target()->point();
			
      //find the leftmost point in the path from e to its twin
      Ccb_halfedge_circulator aux=ccb_e;
      do {
	if (traits->compare_x(aux->target()->point(),e_left)==SMALLER) {
	  e_left=aux->target()->point();
	}
      } while (++aux!=ccb_t);
			
      //find the leftmost point in the path from the twin to e
      aux=ccb_t;
      do {
	if (traits->compare_x(aux->target()->point(),t_left)==SMALLER) {
	  t_left=aux->target()->point();
	}        
      } while (++aux!=ccb_e);
			
      //compare the two left points
      if (traits->compare_x(t_left,e_left) == SMALLER) //e points at hole 
	return Topological_map<Dcel>::remove_edge(e);
      else
	return Topological_map<Dcel>::remove_edge(e->twin());
			
    }
    else {
			
      return Topological_map<Dcel>::remove_edge(e);
    }
		
  }
	
	
  Halfedge_handle vertical_ray_shoot(const Point& p,
				     Locate_type &lt, bool up)
  {
    CGAL_precondition(pl);
    return pl->vertical_ray_shoot(p,lt,up);
  }
	
  
  Halfedge_const_handle vertical_ray_shoot(const Point& p,
					   Locate_type &lt, bool up) const
  {
    CGAL_precondition(pl);
    return ((const Point_location_base*)pl)->vertical_ray_shoot(p,lt,up);
    // The type cast to const is there to ensure that the planar map 
    // is not changed.
  }
	
	
	
  Halfedge_handle locate(const Point& p, Locate_type &lt) {
    CGAL_precondition(pl);
    return pl->locate(p,lt);
  }
	
  Halfedge_const_handle locate(const Point& p, Locate_type &lt) const {
    CGAL_precondition(pl);
    return ((const Point_location_base*)pl)->locate(p,lt);
    // The type cast to const is there to ensure that the planar map 
    // is not changed.
  }
  
protected:  //private implementation
  //returns true if the point  is inside (in the strict sense) of nf
  //algorithm: count the intersections of a vertical ray shoot - are they odd ?
  //assumes outer ccb exists
  
  //CHANGE - for this to work in the arrangement, it should not pass over
  //the curve just inserted (i.e the halfedge)
  
  /*
    bool point_is_in(const Point& p, Face_const_handle nf) const
    {
    int count = 0;
    
    Ccb_halfedge_const_circulator circ = nf->outer_ccb();
    do {
    ++circ;
    } while ((traits->curve_is_vertical(circ->curve()))&&
              circ!=nf->outer_ccb());
    if (circ==nf->outer_ccb() && traits->curve_is_vertical(circ->curve()) )
    return false; 
    //if the whole ccb is vertical then the point is out.
    //else advance to a non vertical curve 
    
    Ccb_halfedge_const_circulator last = circ;
    
    
    do {
    if (traits->point_is_same(circ->target()->point(), p)) 
    //point is on outer ccb 
    return false;
    if (!traits->curve_is_vertical(circ->curve())) {
    
    if ( (traits->curve_get_point_status(circ->curve(),p) == 
    Traits::UNDER_CURVE) && 
    !(traits->point_is_same_x(circ->source()->point(),p)) ) {  
    //point is under curve in the range (source,target]
    
    if (traits->point_is_same_x(circ->target()->point(),p)) {
    //p is exactly under a vertex of the ccb - if next is not on the 
    //same side of the vertical line from p as circ is, 
    //we count one more intersection
    
    Ccb_halfedge_const_circulator next=circ;
    ++next;
    if (traits->curve_is_vertical(next->curve())) {
    //advance to non-vertical edge
    while (traits->curve_is_vertical(next->curve())) {
    ++next;
    }
    }
    if ( (traits->point_is_right(circ->source()->point(),p)&&
    traits->point_is_left(next->target()->point(),p)) ||
    (traits->point_is_left(circ->source()->point(),p)&&
    traits->point_is_right(next->target()->point(),p)) ) {
    
    ++count;
    }
    }
    else {
    ++count;
    }
    }
    }
    } while (++circ!=last);
    
    return (count%2 != 0);  //if count is odd return true
    
    }
  */
  
  bool point_is_in(const Point& p, Halfedge_const_handle ne,
		   const X_curve& ncv) const
    // returns true if the points is inside the planar map.
  {
    int count = 0;
    
    //    Ccb_halfedge_const_circulator circ = nf->outer_ccb();
    Ccb_halfedge_const_circulator circ = ne;
    do {
      ++circ;
    } while (circ!=ne && traits->curve_is_vertical(circ->curve()));
    if (circ==ne && traits->curve_is_vertical(ncv) )
      return false; 
    //if the whole ccb is vertical then the point is out.
    //else advance to a non vertical curve 
    
    Ccb_halfedge_const_circulator last = circ;
    
    
    do {
      X_curve circv;
      if (circ!=ne) { //not on the new halfedge (circ has a curve
	circv=circ->curve();
      }
      else { //maybe doesn't have a curve yet (e.g in arrangement)
	circv=ncv;
      }
      
      if (traits->point_is_same(circ->target()->point(), p)) 
	//point is on outer ccb 
	return false;
      if (!traits->curve_is_vertical(circv)) {
	
	if ( (traits->curve_get_point_status(circv,p) == 
	      Traits::UNDER_CURVE) && 
	     !(traits->point_is_same_x(circ->source()->point(),p)) ) {  
	  //point is under curve in the range (source,target]
	  
	  if (traits->point_is_same_x(circ->target()->point(),p)) {
	    //p is exactly under a vertex of the ccb - if next is not on the 
	    //same side of the vertical line from p as circ is, 
	    //we count one more intersection
	    
	    Ccb_halfedge_const_circulator next=circ;
	    ++next;
	    X_curve nextcv;
	    if (next!=ne) {
	      nextcv=next->curve();
	    }
	    else {
	      nextcv=ncv;
	    }
	    if (traits->curve_is_vertical(nextcv)) {
	      //advance to non-vertical edge
	      while (traits->curve_is_vertical(nextcv)) {
		if (next!=ne) {
		  nextcv=next->curve();
		}
		else {
		  nextcv=ncv;
		}
		++next;
		
	      }
	    }
	    if ( (traits->point_is_right(circ->source()->point(),p)&&
		  traits->point_is_left(next->target()->point(),p)) ||
		 (traits->point_is_left(circ->source()->point(),p)&&
		  traits->point_is_right(next->target()->point(),p)) ) {
	      
	      ++count;
	    }
	  }
	  else {
	    ++count;
	  }
	}
      }
      
      
    } while (++circ!=last);
    
    return (count%2 != 0);  //if count is odd return true
    
  }
  
  
  /////////////////////////////////////////////////////////
  // Assignment functions 
  // Those are temporary functions and it should not be used
public:
  void assign(const Self &pm)
  {
    TPM::assign(pm);
    // traits->assign(pm->traits);
    // bb->assign(pm->bb);
    // pl->assign(pm->pl);
  }

  /*Self& operator=(const Self& pm)
    {
    if (this != &pm)
    {
    X_curve_container l;
    pm.x_curve_container(l);
    clear();
    insert(l.begin(), l.end());
    }
    return *this;
    }
  */
  
  Self& operator=(const Self& pm)
  {
    if (this != &pm) {
      clear();
      assign(pm);
                 
      Halfedge_iterator h_iter;
      for (h_iter = halfedges_begin(); 
	   h_iter != halfedges_end(); 
	   h_iter++, h_iter++)
	pl->insert(h_iter, h_iter->curve());
                 
      for (Vertex_iterator v_iter = vertices_begin(); 
	   v_iter != vertices_end(); 
	   v_iter++)
	bb->insert(v_iter->point());
                 
      for (h_iter = halfedges_begin(); 
	   h_iter !=  halfedges_end(); 
	   h_iter++, h_iter++)
	bb->insert(h_iter->curve());
    }
    return *this;
  }
	
  // used in implementation of operator=(
  void clear() 
  {
    pl->clear();
    TPM::clear();
    /*
      Halfedge_iterator it=halfedges_begin(),prev=it,it_e=halfedges_end();
      while (it!=it_e) {++it;++it;remove_edge(prev);prev=it;}
    */
    bb->clear();
  }
  
protected:
	// used in implementation of operator=(
  void x_curve_container(X_curve_container &l) const
  {
    Halfedge_const_iterator it=halfedges_begin(),it_e=halfedges_end();
    while (it!=it_e) 
      {
	l.push_back(it->curve());
	++it;++it;
      }
  }
  //
  /////////////////////////////////////////////////////////



protected:
	// default initializer for the bounding box.
#include <CGAL/Planar_map_2/Bounding_box_special_initializer.h>

#ifdef CGAL_PM_DEBUG // for private debugging use
	
public:
  void debug()
  {
    if (pl) (pl->debug());
  }

#endif	
	
private:
///////////////////////////////////////////////////////////////////////////
//                 Scanning Arrangement.
///////////////////////////////////////////////////////////////////////// 
template <class Scanner>
bool  scan_planar_map (Scanner& scanner){
  if (!build_dcel(scanner))
    return false;

  return true;
}
  
template <class Scanner>
bool  build_dcel (Scanner& scanner) {
  typedef typename Dcel::Vertex	                  D_vertex;
  typedef typename Dcel::Halfedge                 D_halfedge;
  typedef typename Dcel::Face	                  D_face;
  
  typedef typename  Dcel::Vertex_iterator          D_vetrex_iterator;
  typedef typename  Dcel::Vertex_const_iterator    D_vetrex_const_iterator;
  typedef typename  Dcel::Halfedge_iterator        D_halfedge_iterator;
  typedef typename  Dcel::Halfedge_const_iterator  D_halfedge_const_iterator;
  typedef typename  Dcel::Face_iterator            D_face_iterator;
  typedef typename  Dcel::Face_const_iterator      D_face_const_iterator;

  // keeping a vector of halfedges (to access them easily by their indices).
  std::vector<D_halfedge* >  halfedges_vec;  

  std::vector<D_vertex* >    vertices_vec; 
 
  if ( ! scanner.in())
    return 0;

  scanner.scan_pm_vhf_sizes();
  if ( ! scanner.in()){
    std::cerr << "can't read vhf values"<<std::endl;
    clear();
    return false;
  }

  // read in all vertices
  unsigned int  i;
  for (i = 0; i < scanner.number_of_vertices(); i++) {
    D_vertex* nv = d.new_vertex();
    Point p;

    /*scanner.scan_vertex_attributes (nv);
    if ( ! scanner.in()){
      std::cerr << "can't read vertex attributes"<<std::endl;
      clear();
      return false;
      }*/

    scanner.scan_vertex (nv);
    if ( ! scanner.in()){
      std::cerr << "can't read vertex"<<std::endl;
      clear();
      return false;
    }
    //nv->set_point(p);

    // for debug.
    //std::cout<<"Reading vertex no " <<i<<" point is ";
    //std::cout<<nv->point()<<std::endl;

    vertices_vec.push_back(nv);

    bb->insert(nv->point());
    
    //scanner.skip_to_next_vertex();
    //if ( ! scanner.in()){
    //  std::cerr << "can't skip to next vertex"<<std::endl;
    //  scanner.in().clear( std::ios::badbit);
    //  clear();
    //  return false;
    // }
  }
  
  for (i = 0; i < scanner.number_of_halfedges(); i++, i++){
    D_halfedge *nh = NULL;
    void  *nv1, *nv2;
    std::size_t index1, index2;
    X_curve cv;

    //std::cout<<"Reading Edge no " <<i<<std::endl;

    nh = d.new_edge();
    
    /*scanner.scan_halfedge_attributes (nh);
    if ( ! scanner.in()){
      std::cerr << "can't read halfedge attributes"<<std::endl;
      clear();
      return false;
      }*/

    scanner.scan_index(index1);
    if ( ! scanner.in()){
      std::cerr << "can't read source of halfedge"<<std::endl;
      clear();
      return false;
    }
    cv = scanner.scan_halfedge(nh);
    if ( ! scanner.in()){
      std::cerr << "can't read halfedge"<<std::endl;
      clear();
      return false;
    }
    //nh->set_curve(cv);

    /*scanner.scan_halfedge_attributes (nh->opposite());
    if ( ! scanner.in()){
      std::cerr << "can't read halfedge attributes"<<std::endl;
      clear();
      return false;
      }*/

    scanner.scan_index (index2);
    if ( ! scanner.in()){
      std::cerr << "can't read source of halfedge"<<std::endl;
      clear();
      return false;
    }
    scanner.scan_halfedge(nh->opposite());
    if ( ! scanner.in()){
      std::cerr << "can't read halfedge"<<std::endl;
      clear();
      return false;
    }
    //nh->opposite()->set_curve(cv);

    nv1 = vertices_vec[index1];
    ((D_vertex*) nv1)->set_halfedge(nh); 
    nh->set_vertex((D_vertex*) nv1);
    //for debug
    //std::cout<<((D_vertex*) nv1)->point()<<std::endl;
    
    nv2 = vertices_vec[index2];
    ((D_vertex*) nv2)->set_halfedge(nh->opposite()); 
    nh->opposite()->set_vertex((D_vertex*) nv2);
    //for debug
    //std::cout<<((D_vertex*) nv2)->point()<<std::endl;

    pl->insert(D_halfedge_iterator(nh->opposite()), cv);
    bb->insert(cv);
    
    halfedges_vec.push_back(nh);
    halfedges_vec.push_back(nh->opposite());

    //scanner.skip_to_next_halfedge();
    //if ( ! scanner.in()){
    //  std::cerr << "can't skip to next halfedge"<<std::endl;
    //  scanner.in().clear( std::ios::badbit);
    //  clear();
    //  return false;
    //} 
  }
  
  // read in all facets
  for (i = 0; i < scanner.number_of_faces(); i++) {
    //std::size_t  num_of_holes, num_halfedges_on_outer_ccb;
    
    //std::cout<<"Reading Face no " <<i<<std::endl;
   
    D_face* nf = u_face; //this is the unbounded face.
    if (i > 0)  // else - allocate the bounded face.
      nf = d.new_face();

    scanner.scan_face(nf);
    if ( ! scanner.in()){
      std::cerr << "can't read face"<<std::endl;
      clear();
      return false;
    }
    
    /*
    //if (i > 0){ // not an unbounded face. Scanning the outer ccb.
    scanner.scan_face_number(num_halfedges_on_outer_ccb, i);
    if ( ! scanner.in()){
    std::cerr << "can't read face number"<<std::endl;
    scanner.in().clear( std::ios::badbit);
    clear();
    return false;
    }
    
    // not an unbounded face. Scanning the outer ccb.
    if (num_halfedges_on_outer_ccb > 0) {
    std::size_t  index, prev_index, first_index;
    for (unsigned int j = 0; j < num_halfedges_on_outer_ccb; j++) {
    
    scanner.scan_index(index);
    if ( ! scanner.in()){
    std::cerr << "can't read halfedge's index on face"<<std::endl;
    scanner.in().clear( std::ios::badbit);
    clear();
    return false;
    }
    
    D_halfedge* nh = halfedges_vec[index];
    
    // for debugging.
    //std::cout<<"source of haledge : "<<nh->vertex()->point()<<std::endl;
        
    if (j > 0) {
    D_halfedge* prev_nh = halfedges_vec[prev_index];
    prev_nh->set_next(nh);
    }
    else {
    nf->set_halfedge(nh);
    first_index = index;
    }
    
    nh->set_face(nf); 
    
    prev_index = index;
    }
    
    // making the last halfedge point to the first one (cyclic order).
    D_halfedge* nh = halfedges_vec[first_index];
    D_halfedge* prev_nh = halfedges_vec[prev_index];
      prev_nh->set_next(nh);
      }
      
    scanner.scan_face_number(num_of_holes, i);
    if ( ! scanner.in()){
    std::cerr << "can't read number holes in face"<<std::endl;
    scanner.in().clear( std::ios::badbit);
    clear();
    return false;
    }
    
    // take care the hols.
    for (unsigned int k = 0; k < num_of_holes; k++){
    std::size_t  num_halfedges_on_inner_ccb;
      
    scanner.scan_face_number(num_halfedges_on_inner_ccb, i);
    if ( ! scanner.in()){
    std::cerr << "can't read number of halfedges in hole"<<std::endl;
    scanner.in().clear( std::ios::badbit);
    clear();
    return false;
    }
    
      std::size_t  index, prev_index, first_index;
      for (unsigned int j = 0; j < num_halfedges_on_inner_ccb; j++) {
      scanner.scan_index(index);
      if ( ! scanner.in()){
      std::cerr << "can't read halfedge's index on hole"<<std::endl;
      scanner.in().clear( std::ios::badbit);
      clear();
      return false;
      }
      
      D_halfedge* nh = halfedges_vec[index];
        
      // for debugging.
      //std::cout<<"source of haledge : "<<nh->vertex()->point()<<std::endl;
      
      if (j > 0) {
      D_halfedge* prev_nh = halfedges_vec[prev_index];
      prev_nh->set_next(nh);
      }
      else {
      nf->add_hole(nh);
      first_index = index;
      }
      
      nh->set_face(nf); 
        
      prev_index = index;
      }
      
      // making the last halfedge point to the first one (cyclic order).
      D_halfedge* nh = halfedges_vec[first_index];
      D_halfedge* prev_nh = halfedges_vec[prev_index];
      prev_nh->set_next(nh);
      }
      scanner.skip_to_next_face(i);
      if ( ! scanner.in()){
      std::cerr << "can't skip to next face"<<std::endl;
      scanner.in().clear( std::ios::badbit);
      clear();
      return false;
      }*/
  }
  
  if ( ! scanner.in() ) {
    scanner.in().clear( std::ios::badbit);
    return false;
  } 
  
  return true;
}	
  
  // #else was removed by eti.
protected:
	
  Point_location_base *pl;
  Bounding_box_base *bb;
  Traits_wrap*  traits;
private:
  bool use_delete_pl;
  bool use_delete_bb;
  bool use_delete_traits;
};

CGAL_END_NAMESPACE


#endif








