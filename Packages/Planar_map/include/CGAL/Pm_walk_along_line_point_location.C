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
// file          : include/CGAL/Pm_walk_along_line_point_location.C
// package       : pm (4.08)
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Oren Nechushtan <theoren@math.tau.ac.il>
//                 Iddo Hanniel <hanniel@math.tau.ac.il>
//
//
// coordinator   : Tel-Aviv University (Dan Halperin halperin<@math.tau.ac.il>)
//
// Chapter       : 
// ======================================================================
#ifndef CGAL_PM_WALK_ALONG_LINE_POINT_LOCATION_C
#define CGAL_PM_WALK_ALONG_LINE_POINT_LOCATION_C

#ifndef CGAL_PM_WALK_ALONG_LINE_POINT_LOCATION_H
#include <CGAL/Pm_walk_along_line_point_location.h>
#endif // CGAL_PM_WALK_ALONG_LINE_POINT_LOCATION_H

CGAL_BEGIN_NAMESPACE

//if unbounded face - returns NULL or some edge on unbounded face 
//if its a vertex returns a halfedge pointing _at_ it
template <class Planar_map>
Pm_walk_along_line_point_location<Planar_map>::Halfedge_handle
Pm_walk_along_line_point_location<Planar_map>::locate(
                                                      const Point& p, 
                                                      Locate_type& lt) const
{
  Face_handle f=pm->unbounded_face(),last=pm->faces_end();  
  // invariant: p always lies in f's interior or holes
  
  Halfedge_handle e=pm->halfedges_end(); // closest halfedge so far
  lt = Planar_map::UNBOUNDED_FACE;
  while(f!=last) // stop at innermost level
    {
      last = f;
      Holes_iterator it=f->holes_begin(),end=f->holes_end();
      while(it!=end && last==f)     // handle holes
        {
          if (find_closest(p,*it,true,true,e,lt))
            switch (lt)
              {
              case Planar_map::UNBOUNDED_FACE:
                break;
              case Planar_map::FACE:
                walk_along_line(p,true,true,e,lt);
                switch(lt)
                  {
                  case Planar_map::VERTEX:
                    
#ifdef CGAL_PM_DEBUG
                    
                    CGAL_assertion(
                      traits->point_is_same(e->target()->point(),p)
                      );
                    
#endif
                    
                  case Planar_map::EDGE:
                    return e;
                  case Planar_map::FACE:
                    f=e->face();
                    break;
                  default:
                    CGAL_assertion(lt==Planar_map::FACE||
					lt==Planar_map::EDGE||
                                   lt==Planar_map::VERTEX);
                    break;
                  }
                break;
              case Planar_map::VERTEX:
                
#ifdef CGAL_PM_DEBUG
                
                CGAL_assertion(traits->point_is_same(e->target()->point(),p));
                
#endif
                
              case Planar_map::EDGE:
                return e;
              default:
                CGAL_assertion(
                               lt==Planar_map::UNBOUNDED_FACE||
                               lt==Planar_map::FACE||
                               lt==Planar_map::EDGE||
                               lt==Planar_map::VERTEX);
              }
          ++it;
        }
    }
  if (lt==Planar_map::UNBOUNDED_FACE && f!=pm->unbounded_face()) 
    lt=Planar_map::FACE;
  if (e==pm->halfedges_end() && pm->number_of_halfedges()>0) 
    return Halfedge_handle(*(pm->unbounded_face()->holes_begin()));
  
  return e;
}

template <class Planar_map>
Pm_walk_along_line_point_location<Planar_map>::Halfedge_handle
Pm_walk_along_line_point_location<Planar_map>::locate(
						      const Point& p, 
						      Locate_type& lt){
  ((Bounding_box*)get_bounding_box())->insert(p);
  Halfedge_handle h=((cPLp)this)->locate(p,lt);
  if (!((Bounding_box*)get_bounding_box())->locate(p,lt,h))
    h=((cPLp)this)->locate(p,lt);
  return h;
}

template <class Planar_map>
Pm_walk_along_line_point_location<Planar_map>::Halfedge_handle
Pm_walk_along_line_point_location<Planar_map>::vertical_ray_shoot(
                                                 const Point& p, 
                                                 Locate_type& lt, 
                                                 bool up) const
{
  Face_handle f=pm->unbounded_face(),last=pm->faces_end();  
  // p always lies in f's interior
  Halfedge_handle e=pm->halfedges_end(); // closest halfedge so far
  lt = Planar_map::UNBOUNDED_FACE;
  while(f!=last) // stop at innermost level
    {
      last = f;
      Holes_iterator it=f->holes_begin(),end=f->holes_end();
      while(it!=end && last==f)     // handle holes
        {
          if (find_closest(p,*it,up,false,e,lt))
            {
              switch (lt)
                {
                case Planar_map::VERTEX:
                  
#ifdef CGAL_PM_DEBUG
                  
                  CGAL_assertion(
                                 traits->point_is_same_x(
                                                         e->target()->point(),
                                                         p)
                                 );
                  
#endif
                  
                case Planar_map::EDGE:
                  
#ifdef CGAL_PM_WALK_DEBUG
                  std::cerr << "\ncalling walk_along_line(" << p << "," 
                            << up << ",false," << e->curve() << "," << lt 
                            << ")";
#endif
                  
                  walk_along_line(p,up,false,e,lt);
                  switch(lt)
                    {
                    case Planar_map::VERTEX:
                      f=e->twin()->face();
                      break;
                    case Planar_map::EDGE:
                      
#ifdef CGAL_PM_DEBUG
                      
                      CGAL_assertion(
                                     up == traits->point_is_left_low(
                                       e->target()->point(),
                                       e->source()->point()
                                       )
                                     );
                      
#endif                      
                    case Planar_map::UNBOUNDED_FACE:
                      break;
                    default:
                      CGAL_assertion(
                                     lt==Planar_map::UNBOUNDED_FACE||
                                     lt==Planar_map::EDGE||
                                     lt==Planar_map::VERTEX
                                     );
                      break;
                    }
                  break;
                default:
                  CGAL_assertion(
                                 lt==Planar_map::EDGE||
                                 lt==Planar_map::VERTEX);
                  break;
                }
            }
          ++it;
        }
    }
  if (e==pm->halfedges_end()) 
    lt=Planar_map::UNBOUNDED_FACE;
  
  /* symmetric diagrams for downward ray shoot
     x
     |\
     | x
     p    => VERTEX
     
     x
     \
     |\  => EDGE
     | x
     p
     
     x        x
     |        |
     p=x   or  p   => EDGE
     |
     x
     
     x
     |
     x
     => VERTEX
     p
  */
  return e;
}


template <class Planar_map>
Pm_walk_along_line_point_location<Planar_map>::Halfedge_handle
Pm_walk_along_line_point_location<Planar_map>::vertical_ray_shoot(
                                                 const Point& p, 
                                                 Locate_type& lt, bool up){
  /* Make sure the source point is in the bounding box on the output */
  ((Bounding_box*)get_bounding_box())->insert(p);
  Halfedge_handle h=((cPLp)this)->vertical_ray_shoot(p,lt,up);
  /* Apply the bounding box on the output */
	if (!((Bounding_box*)get_bounding_box())->vertical_ray_shoot(p,lt,up,h))
	{
		h=((cPLp)this)->vertical_ray_shoot(p,lt,up);
		CGAL_assertion(lt!=Planar_map::UNBOUNDED_FACE);
	}
	return h;
}

///IMPLEMENTATION ////////////////////////////////////////

/*
 Point p -  the source of the vertical ray shoot
 bool up -  true if the ray shoot is upwards, otherwise false (downwards)
 bool including - true if the ray includes its source, false otherwise

 The intersection of the vertical ray shoot is represented using:
 Locate_type lt - the type of the intersection
 Halfedge_handle e - an halfedge representing the intersection depending on
{
 lt == point, ( a vertex in the planar map )
      find the first halfedge pointing at p, when going clockwise if 
      highest==true - start from 12 oclock, else start from 6 oclock 
      DEBUG: highest = up or !up ???
 lt == curve, ( an halfedge in the planar map )
      DEBUG: ?
 lt == face, ( a face in the planar map )
      DEBUG: ?
 lt == unbounded face, ( an unbounded face in the planar map)
      DEBUG: ?
}

*/


template <class Planar_map>
void Pm_walk_along_line_point_location<Planar_map>::walk_along_line(
						 const Point& p,
						 bool up,
						 bool including, 
					      // bool type,
						 Halfedge_handle& e,
						 Locate_type& lt) const 
{
  bool type = including;
  Face_handle f=(type || lt!=Planar_map::VERTEX) ? e->face() : 
    e->twin()->face()    ,last;
  do
    {
      last=f;
      /*
        
          x
         /
        x
        
        p

        A situation where CGAL_assertion(f!=pm->unbounded_face()) doesn't hold.
     */

#ifdef CGAL_PM_WALK_DEBUG

      std::cerr << "\n pre find_closest(" << p << ", , " << up 
		<< "," << including << ",(" << e->source()->point() 
		<< "," << e->target()->point()  << ")," << lt << ");";
#endif

        if (f!=pm->unbounded_face()) 
          {

#ifdef CGAL_PM_DEBUG

      bool found = 

#endif
        
        find_closest(p,f->outer_ccb(),up,including,e,lt);

#ifdef CGAL_PM_DEBUG
      
      CGAL_assertion(found);
      
#endif
          }
      
        f=(type || lt!=Planar_map::VERTEX) ? e->face() : e->twin()->face();
    }
  while((type == (lt==Planar_map::UNBOUNDED_FACE)) && last!=f);
}

template <class Planar_map>
bool Pm_walk_along_line_point_location<Planar_map>::find_closest(
	     const Point& p,
	     const Ccb_halfedge_circulator& c,
	     bool up,bool including, // bool type,
	     Halfedge_handle& e,
	     Locate_type& lt) const
{
  bool type = including; // for possible future implementation
  bool intersection = e!=pm->halfedges_end(); 
  // used to answer is it known that ray shoot intersects curves?
  bool inside = false; // used to calculate if point is inside ccb
  Ccb_halfedge_circulator curr=c;
  do
    {
      const X_curve& cv=curr->curve(),&ecv=e->curve();
      const Point& p1=traits->curve_source(cv),&p2 = traits->curve_target(cv);
      Curve_point_status s = traits->curve_get_point_status(cv, p);
      if ( s == (up ? Traits::UNDER_CURVE : Traits::ABOVE_CURVE)
// && !traits->curve_is_vertical(cv)
	) 
        /* cv is a non vertical curve intersecting the vertical ray shoot 
               x
             / 
           x |     ( for a vertical ray shoot upwards )
             |
             p

        The vertical case is excluded for lexicographically the curve is not intersecting with the ray:

           x        |  x
           |        | /
           x   =>   |x 
           |        |
           p        p
        */

	{
	  if (traits->point_is_left_low(p,p1)!=traits->point_is_left_low(p,p2))
	    {
	      inside = !inside; // count parity of curves intersecting ray in their interior
	    }
	  if (!intersection || 
	      traits->curve_compare_at_x(ecv, cv , p) == (up ? LARGER : SMALLER)
	      )
	    {
              // orient e leftlow
              if ( up == traits->point_is_left_low(
						   curr->target()->point(),
						   curr->source()->point())
		   ) 
                e = curr ;
              else
                e = curr->twin();
              if (!type)
                if (!traits->curve_is_vertical(cv))
                  {
                    if (traits->point_is_same_x(e->source()->point(),p)) 
		      { e=e->twin(); lt=Planar_map::VERTEX;}
                    else if (traits->point_is_same_x(e->target()->point(),p)) 
		      { lt=Planar_map::VERTEX;}
                    else lt=Planar_map::EDGE;
                  }
                else // for p is within e'th x range
                  lt=Planar_map::VERTEX;
              intersection = true;
	    }
	  else if (e!=curr && e!=curr->twin() && 
		   traits->curve_compare_at_x(ecv, cv , p) == EQUAL)
	    /* first intersection point of ray and curve is an end point like

	                   x x
	   x--x---x        |/ 
	      |            x
	      |            |
	      p       or   p

	    */
	    {
	      Point q = traits->curve_source(cv);
	      if (q!=traits->curve_source(ecv) &&
		  q!=traits->curve_target(ecv))
	      q=traits->curve_target(cv);
	      if ((up ? traits->curve_compare_at_x_from_bottom(ecv,cv,q) :
		  traits->curve_compare_at_x_from_top(ecv,cv,q)) == LARGER)
                if (type != (traits->point_is_same(curr->target()->point(),q))) 
                  e = curr;
                else
                  e= curr->twin();
              if (!type) lt=Planar_map::VERTEX; 
	      // lt should be already Planar_map::VERTEX
            }

#ifdef CGAL_PM_DEBUG

	  else
	    {
	      CGAL_assertion(
			 traits->curve_compare_at_x(ecv, cv , p) == 
			 (!up ? LARGER : SMALLER) ||
			 e==curr || e==curr->twin()
	      );
	    }

#endif // CGAL_PM_DEBUG

        }
      else if ( s == Traits::ON_CURVE )
	{
	  if (!including)
	  /* The vertical ray shoot is not including p itself,
	     thus we are interested only in vertical curves that
	     extend upwards
	     Remark:
	     The Locate type is always EDGE
	  */
	    {
	      if (traits->curve_is_vertical(cv) && traits->point_is_higher(traits->curve_highest(cv),p))

		/*
		  x       x
		  |       |
		 p=x  or  p
		          |
			  x
		*/

		{
		  lt = Planar_map::EDGE;
		  if (up==traits->point_is_left_low(curr->target()->point(),curr->source()->point()))
                    e = curr;
                  else 
                    e = curr->twin();

#ifdef CGAL_PM_WALK_DEBUG

                  std::cerr << "\n find_closest(" 
			    << p << ", , " << up << "," << including 
			    << ",(" << "," << e->target()->point() 
			    << ")," << lt << ");";

#endif

		  return true;
		}
	    }
	  else // including
	    {

	      // p is in interior of curr->curve();
	      if ( !traits->point_is_same(p,traits->curve_source(cv)) && 
		   !traits->point_is_same(p,traits->curve_target(cv)))
		{
		  lt = Planar_map::EDGE;
		  if (up==traits->point_is_left_low(curr->target()->point(),curr->source()->point()))
                    e = curr;
                  else 
                    e = curr->twin();
		}
	      else // end point
		{
		  lt = Planar_map::VERTEX;

#ifdef CGAL_PM_DEBUG
		  
		  CGAL_assertion(curr!=pm->halfedges_end());

#endif
                  if (traits->point_is_same(curr->target()->point(),p))
                    e = find_vertex_representation(curr,p,up);	
                  else
                    e = find_vertex_representation(curr->twin(),p,up);

#ifdef CGAL_PM_DEBUG

		  CGAL_assertion(traits->point_is_same(e->target()->point(),p));

#endif

		}

#ifdef CGAL_PM_WALK_DEBUG

              std::cerr << "\n find_closest(" 
			<< p << ", , " << up << "," << including 
			<< ",(" << e->source()->point() 
			<< "," << e->target()->point()  
			<< ")," << lt 
			<< ");";

#endif

	      return true;
	    }
	}
      ++curr;
    }
  while (curr!=c);

  if (!intersection) {
    lt=Planar_map::UNBOUNDED_FACE;
    return false;
  }
  if (type) lt = (inside ? Planar_map::FACE : Planar_map::UNBOUNDED_FACE);

#ifdef CGAL_PM_WALK_DEBUG

  std::cerr << "\n find_closest(" << p << ", , " << up << "," << including 
	    << ",(" << e->source()->point() << "," << e->target()->point()  
	    << ")," << lt << ");";

#endif

  return true;
}

CGAL_END_NAMESPACE

#endif
