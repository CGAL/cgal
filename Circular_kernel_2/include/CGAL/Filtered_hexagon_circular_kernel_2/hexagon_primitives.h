// Copyright (c) 2003-2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s)     : Monique Teillaud, Sylvain Pion, Constantinos Tsirogiannis , Pedro Machado

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_HEXAGON_PRIMITIVES_H
#define CGAL_HEXAGON_PRIMITIVES_H

#include <vector>
#include <stack>
#include <CGAL/Polygon_2.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/number_utils.h>
#include <CGAL/Line_2.h>
#include <CGAL/Object.h>
#include <CGAL/Lazy.h>
#include <CGAL/Circular_kernel_2/internal_functions_on_circular_arc_2.h>
#include <CGAL/Circular_kernel_converter.h>

namespace CGAL{
 namespace internal{

using CGAL::CircularFunctors::advanced_make_xy_monotone;

CGAL::Polygon_2 < CGAL::Simple_cartesian< double > >
  construct_polygon_from_bbox(const CGAL::Bbox_2& bb )
  {
    typedef CGAL::Simple_cartesian< double > dK;
    typedef dK::Point_2                      Point_2;
    
    CGAL::Polygon_2 < dK > pgn;
  
    
    pgn.insert(pgn.vertices_end(),Point_2(bb.xmin(), bb.ymin()));
    pgn.insert(pgn.vertices_end(),Point_2(bb.xmax(), bb.ymin()));
    pgn.insert(pgn.vertices_end(),Point_2(bb.xmax(), bb.ymax()));
    pgn.insert(pgn.vertices_end(),Point_2(bb.xmin(), bb.ymax()));
    
    return pgn;    
  }


template <typename CK>
CGAL::Polygon_2 < CGAL::Simple_cartesian< double > >
 construct_bounding_hexagon_for_line_arc_2( const typename CK::Line_arc_2 &a) 
   {
	typedef  typename CK::Line_arc_2             Line_arc_2;
	typedef  CGAL::Simple_cartesian<double>      dK;
	typedef  dK::Point_2                         Point_2;
	typedef  CGAL::Polygon_2<dK>                 Polygon_2;
	   
	CGAL::Bbox_2 src_bb= a.source().bbox(),
   		     trgt_bb= a.target().bbox();
		     	  
	if( (src_bb.xmin()<=trgt_bb.xmin() && src_bb.xmax()>=trgt_bb.xmax()) ||
	    (src_bb.ymin()<=trgt_bb.ymin() && src_bb.ymax()>=trgt_bb.ymax()) ||
	    (src_bb.xmin()>=trgt_bb.xmin() && src_bb.xmax()<=trgt_bb.xmax()) ||	  
	    (src_bb.ymin()>=trgt_bb.ymin() && src_bb.ymax()<=trgt_bb.ymax())  )
	    
	    return construct_polygon_from_bbox(src_bb+trgt_bb);
	    	
	Polygon_2 pgn;
	
	bool tmp;
	
	if( (tmp=(src_bb.xmin()<trgt_bb.xmin() && src_bb.ymax()>trgt_bb.ymax())) ||
	              (trgt_bb.xmin()<src_bb.xmin() && trgt_bb.ymax()>src_bb.ymax()))
	{
	  CGAL::Bbox_2 bb1=(tmp)? src_bb : trgt_bb,
	               bb2=(tmp)? trgt_bb : src_bb;
	
	  pgn.push_back(Point_2(bb1.xmin(),bb1.ymin()));
	  pgn.push_back(Point_2(bb2.xmin(),bb2.ymin()));
	  pgn.push_back(Point_2(bb2.xmax(),bb2.ymin()));
	  pgn.push_back(Point_2(bb2.xmax(),bb2.ymax()));
	  pgn.push_back(Point_2(bb1.xmax(),bb1.ymax()));
	  pgn.push_back(Point_2(bb1.xmin(),bb1.ymax()));
	}
	else if( ( tmp=(src_bb.xmin()<trgt_bb.xmin() && src_bb.ymin()<trgt_bb.ymin())) ||
	              (trgt_bb.xmin()<src_bb.xmin() && trgt_bb.ymin()<src_bb.ymin()))
	{
	  CGAL::Bbox_2 bb1=(tmp)? src_bb : trgt_bb,
	               bb2=(tmp)? trgt_bb : src_bb;
	
	  pgn.push_back(Point_2(bb1.xmax(),bb1.ymin()));
	  pgn.push_back(Point_2(bb2.xmax(),bb2.ymin()));
	  pgn.push_back(Point_2(bb2.xmax(),bb2.ymax()));
	  pgn.push_back(Point_2(bb2.xmin(),bb2.ymax()));
	  pgn.push_back(Point_2(bb1.xmin(),bb1.ymax()));
	  pgn.push_back(Point_2(bb1.xmin(),bb1.ymin()));	       
	}
	
	
	return pgn;
   }



template <typename CK, typename Nested_pair>
CGAL::Polygon_2 < CGAL::Simple_cartesian< double > >
 construct_bounding_hexagon_2( const Nested_pair &n) 
   {
      typedef typename CK::Circular_arc_2      Circular_arc_2;
      typedef  CGAL::Simple_cartesian<double>  dK;
      typedef  dK::Point_2                     Point_2;
      typedef  CGAL::Polygon_2<dK>             Polygon_2;
      typedef  dK::Line_2                      Line_2;
      typedef  dK::Direction_2                 Direction_2;
      typedef  dK::Vector_2                    Vector_2;
      typedef  typename CK::Root_of_2          Root_of_2;
      typedef  typename CK::FT 	               FT;
		
		
      Polygon_2 plgn;

      const Circular_arc_2 *a=CGAL::object_cast<Circular_arc_2>(&n.first);
	
				
      // For these  and the is_on_ booleans it is required to use a lazy nt
      CGAL_assertion(a->is_x_monotone());
      CGAL_assertion(a->is_y_monotone());

		
      CGAL::Bbox_2 left_bb(a->left().bbox()),
      right_bb(a->right().bbox());


      double ymin=(left_bb.ymin()<right_bb.ymin())?left_bb.ymin():right_bb.ymin();
      double ymax=(left_bb.ymax()>right_bb.ymax())?left_bb.ymax():right_bb.ymax();
	

      bool is_on_upper_part=n.second.first;
      bool is_on_left_side =n.second.second ;
	
		
      // Define a hypothetical worst case center (_a,_b)
      double _a= (is_on_left_side)? right_bb.xmin() : left_bb.xmax() , _b;
		
      if(is_on_upper_part)
	_b= (is_on_left_side)? left_bb.ymax(): right_bb.ymax();
      else
	_b= (is_on_left_side)? left_bb.ymin(): right_bb.ymin();


      //Define a worst case radius
      double  r1= (is_on_left_side)? right_bb.xmin()-left_bb.xmin():
	right_bb.xmax()-left_bb.xmax();
      double  r2= (is_on_upper_part)? right_bb.ymax()-left_bb.ymax():
	left_bb.ymin()-right_bb.ymin();
			    
      if(!is_on_left_side)
	r2=-r2;
		
      Interval_nt<> rad= CGAL::max(r1,r2);
      Interval_nt<> sqr2=sqrt(Interval_nt<>(2));
      Interval_nt<> _tngx= (is_on_left_side)? Interval_nt<>(_a)-rad/sqr2:
					      Interval_nt<>(_a)+rad/sqr2;
      Interval_nt<> _tngy= (is_on_upper_part)? Interval_nt<>(_b)+rad/sqr2:
					       Interval_nt<>(_b)-rad/sqr2;		

      double tngx= (is_on_left_side)? to_interval(_tngx).first:
				      to_interval(_tngx).second;
      double tngy= (is_on_upper_part)? to_interval(_tngy).second:
				       to_interval(_tngy).first;

      // Indicate direction
      double dir= ((is_on_upper_part==is_on_left_side)? 1:-1);
      CGAL::Interval_nt<> space(to_interval(tngy-dir*tngx).first,to_interval(tngy-dir*tngx).second);
		
      double x1= (is_on_left_side)?left_bb.xmin():right_bb.xmax();
      double y2;
		
      if(is_on_upper_part)
        y2= (is_on_left_side)? right_bb.ymax(): left_bb.ymax();
      else
	y2= (is_on_left_side)? right_bb.ymin(): left_bb.ymin();
		
   
      CGAL::Interval_nt<> _y1= dir*x1+space;
      CGAL::Interval_nt<> _x2= (y2-space)*dir;
		
      double y1=(is_on_upper_part)? to_interval(_y1).second: to_interval(_y1).first;
      double x2=(is_on_left_side)? to_interval(_x2).first: to_interval(_x2).second;
		
				
      //Points 1 & 2 or just the 1st one
		
      if (x2<left_bb.xmin()||x2>right_bb.xmax()||y1<ymin||y1>ymax)
        plgn.push_back(Point_2(((is_on_left_side)?left_bb.xmin():right_bb.xmax()),
				                  ((is_on_upper_part)? ymax : ymin)));
      else if (is_on_left_side){
	plgn.push_back(Point_2(x1,y1));
	plgn.push_back(Point_2(x2,y2));
      } else{
	plgn.push_back(Point_2(x2,y2));
	plgn.push_back(Point_2(x1,y1));
      }
		
      //Point 3 -- we get it through the righmost endpoint
      if(is_on_upper_part==is_on_left_side)
	plgn.push_back(Point_2(right_bb.xmax(),right_bb.ymax()));
      else 
	plgn.push_back(Point_2(right_bb.xmax(),right_bb.ymin()));

		
	if(is_on_upper_part && is_on_left_side && 
	   (right_bb.ymin()<=left_bb.ymin() || left_bb.xmax()>=right_bb.xmax() ) )
	      plgn.push_back(Point_2(right_bb.xmax(), left_bb.ymin()));
	else if(is_on_upper_part && !is_on_left_side && 
	   (left_bb.ymin()<=right_bb.ymin() || right_bb.xmin()<= left_bb.xmin()))
	      plgn.push_back(Point_2(left_bb.xmin(), right_bb.ymin()));
	else if(!is_on_upper_part && is_on_left_side && 
	   (right_bb.ymax()>=left_bb.ymax() || left_bb.xmax()>=right_bb.xmax()))
	      plgn.push_back(Point_2(right_bb.xmax(), left_bb.ymax()));
	else if(!is_on_upper_part && !is_on_left_side && 
           (left_bb.ymax()>=right_bb.ymax() || right_bb.xmin()<= left_bb.xmin()))
	      plgn.push_back(Point_2(left_bb.xmin(), right_bb.ymax()));	
	else
	{	
		
	  // Point 4 
	  double _x1,_y1;
	  if(is_on_left_side)
	    _x1=right_bb.xmax();
	  else
	    _x1=right_bb.xmin();
		
	  if(is_on_upper_part)
	    _y1=right_bb.ymin();
	  else
	    _y1=right_bb.ymax();
		
	  plgn.push_back(Point_2(_x1,_y1));
	
		
	  //Point 5	
			
	  if(is_on_left_side)
	    _x1=left_bb.xmax();
	  else
	    _x1=left_bb.xmin();
		
	  if(is_on_upper_part)
	    _y1=left_bb.ymin();
	  else
	    _y1=left_bb.ymax();
			
	  plgn.push_back(Point_2(_x1,_y1));
	}
		
		
	//Point 6
	if(is_on_upper_part==is_on_left_side)
	  plgn.push_back(Point_2(left_bb.xmin(),left_bb.ymin()));
	else
	  plgn.push_back(Point_2(left_bb.xmin(),left_bb.ymax()));

		
        //Reverse (except the last one) the order of the vertices if the hexagon
        //is on the upper part , so that all the hexagons
        //appear with their vertices in counterclockwise order		
		
        if(is_on_upper_part)
        {
          std::stack<Point_2> stck;
          int siz=plgn.size();
        
          for(int i=0;i<siz-1;++i)
            stck.push(plgn[i]);

          Point_2 last_one=plgn[siz-1];
          plgn.clear();

          

          for(int i=0;i<siz-1;++i)
          {
            plgn.push_back(stck.top());
            stck.pop();
          }

          plgn.push_back(last_one);
        }
          
		
		
	return plgn;       
	}


template<typename CK, typename Output_iterator>
Output_iterator construct_bounding_hexagons_2(const typename CK::Circular_arc_2 &a, Output_iterator res)		
  {

    typedef std::pair< bool,bool >  Position;
    typedef std::pair<CGAL::Object,Position>           Verbose_type; 
    typedef std::vector<Verbose_type >                 Xy_mon_vector;	
    typedef CGAL::Simple_cartesian<double>             K;
    typedef std::vector< CGAL::Polygon_2<K> >          Polygon_vector;

    Xy_mon_vector xy_arcs;
    Polygon_vector polygons;

    advanced_make_xy_monotone<CK>(a,std::back_inserter(xy_arcs));


    for(unsigned int i=0;i<xy_arcs.size();i++)
      *res++=construct_bounding_hexagon_2<CK>(xy_arcs.at(i));

    return res;

  }	

 bool  _inner_intersect_hexagon_predicate(const CGAL::Polygon_2<Simple_cartesian< double > > &a,
                                          const CGAL::Polygon_2<Simple_cartesian<double> > &b)
 {

    typedef Simple_cartesian<double> CK;
    typedef CK::RT      RT;
    typedef CK::Point_2 Point_2;
    typedef CGAL::Polygon_2<CK> Polygon_2;
    typedef Simple_cartesian<Interval_nt<> > intrv_CK;
    typedef intrv_CK::Line_2  intrv_Line_2;
    typedef intrv_CK::Point_2  intrv_Point_2;

    
  bool pred=(a[0].x()!=a[1].x() && a[0].y()!=a[1].y());
      
          Interval_nt<> d1=a[0].x(),d2=a[0].y();
          Interval_nt<> e1=a[1].x(),e2=a[1].y();
          intrv_Point_2 a0( d1, d2 );
          intrv_Point_2 a1( e1, e2 );

          intrv_Line_2 tmp_ln(a0,a1);

   	
      if(pred && !tmp_ln.is_vertical())
        {
	  int i;
	
	  for(i=0;i<b.size();i++)
	    {
	      
                //if(b[i]==a[0] || b[i]==a[1])
                //  return true;	     

	      std::pair<double,double> app_y=tmp_ln.y_at_x(b[i].x()).pair();

              if(b[i].y()>=app_y.first )
	        if(i==0)
	          break;
	        else 
	          return true;
	    }
	
	  if(i!=0)
	    return false;

	 }
	

	  //Condition implied: the segment used for testing is vertical
	 if(a[2+pred].x()==a[3+pred].x())
	   return true;
	 
          Interval_nt<> d12=a[2+pred].x(),d22=a[2+pred].y();
          Interval_nt<> e12=a[3+pred].x(),e22=a[3+pred].y();
          intrv_Point_2 a2( d12, d22 );
          intrv_Point_2 a3( e12, e22 );

          intrv_Line_2 tmp_ln2(a2,a3);
	 
	 for(int i=0;i<b.size();i++) 
	 {
	   std::pair<double,double> app_y=tmp_ln2.y_at_x(b[i].x()).pair();

           if(b[i].y()<=app_y.second ) 
	     return true;	
	 }
			
	 return false;
			
}


bool _do_intersect_hexagon_2(const CGAL::Polygon_2<Simple_cartesian< double > > &p1,
			     const CGAL::Polygon_2<Simple_cartesian<double> > &p2)
  {
    typedef Simple_cartesian<double> CK;
    typedef CK::RT      RT;
    typedef CK::Point_2 Point_2;
    typedef CGAL::Polygon_2<CK> Polygon_2;
    typedef Simple_cartesian<Interval_nt<> > intrv_CK;
    typedef intrv_CK::Line_2  intrv_Line_2;
    typedef intrv_CK::Point_2  intrv_Point_2;
		
    Polygon_2 a,b;
    bool intersect,tmp=0;
    Point_2 frst,scnd;
	
    int size1=p1.size(),
    	size2=p2.size();

    if(size1==4 && size2==4)
      return true;

    if(size1==4)
      {	
	a=p2;
	b=p1;
      } else if (size2==4)
      {
	a=p1;
	b=p2;
      } else 
      {
		 
	CGAL::Bbox_2 test1(p1.left_vertex()->x(),p1.bottom_vertex()->y(), 
			           p1.right_vertex()->x() ,p1.top_vertex()->y());
	
	for(int i=0;i<size2;i++)
	  if( tmp=(p2[i].x()>=test1.xmin() && p2[i].x()<=test1.xmax() &&
              p2[i].y()>=test1.ymin() && p2[i].y()<=test1.ymax()))
			break;
		
	if(tmp)
	{
	  a=p1;
	  b=p2;
	}else
	{
	  a=p2;
	  b=p1;
	}

      }

       intersect=_inner_intersect_hexagon_predicate(a,b);

       if(intersect==false || b.size()==4)
         return intersect;
       else
         return _inner_intersect_hexagon_predicate(b,a);

   }

// This is meant for x_monotone and for non-x_monotone cases
template < typename Hex_iterator1, typename Hex_iterator2 >
bool do_intersect_hexagons_2(Hex_iterator1 a_begin, Hex_iterator1 a_end, Hex_iterator2 b_begin, Hex_iterator2 b_end)
  {

    for(Hex_iterator1 it1=a_begin;it1!=a_end;it1++)
      for(Hex_iterator2 it2=b_begin;it2!=b_end;it2++) 
	if( !(it2->top_vertex()->y() < it1->bottom_vertex()->y() || 
	      it2->bottom_vertex()->y() > it1->top_vertex()->y() ||     
	      it2->right_vertex()->x() < it1->left_vertex()->x() || 
	      it2->left_vertex()->x() > it1->right_vertex()->x() ))     
	{
	  bool tmp = _do_intersect_hexagon_2(*it1,*it2);
          if (tmp)
	    return true;
            
        }

	return false;

  }// do_intersect_hexagons_2


 }//namespace internal


//Lazy like functors that there is no use to be included in the kernel

template < class CK, class Hexagon>
  class Hexagon_construction_with_interval_2 {

  typedef typename CK::Circular_arc_2  Circular_arc_2;
  typedef typename CK::Line_arc_2  Line_arc_2;
  typedef CGAL::Simple_cartesian<CGAL::Interval_nt<> >                   FK;
  typedef CGAL::Circular_kernel_2< FK,CGAL::Algebraic_kernel_for_circles_2_2<FK::RT> >   CK2;
  typedef CGAL::Circular_kernel_converter<CK,CK2>                          Conv;

  public:

  template < class OutputIterator>
  OutputIterator  
  operator()(const Circular_arc_2 &a, OutputIterator res) const
    {
      Conv cnv;
      static const bool Protection = true;
      try{return internal::construct_bounding_hexagons_2<CK2>(cnv(a),res);}
      catch (Uncertain_conversion_exception)
      {
         CGAL::Protect_FPU_rounding<!Protection> P(CGAL_FE_TONEAREST);
         return internal::construct_bounding_hexagons_2<CK>(a,res);
      }

   }

  Hexagon operator()(const Line_arc_2 &a) const
    {
      Conv cnv;
      static const bool Protection = true;
      try{return internal::construct_bounding_hexagon_for_line_arc_2<CK2>(cnv(a));}
      catch (Uncertain_conversion_exception)
      {
         CGAL::Protect_FPU_rounding<!Protection> P(CGAL_FE_TONEAREST);
         return internal::construct_bounding_hexagon_for_line_arc_2<CK>(a);
      }

   }

};


template < class CK, class Hexagon>
 class Hexagon_construction_on_lazy_kernel_2 {

  typedef typename CK::Circular_arc_2  Circular_arc_2;
  typedef typename CK::Line_arc_2  Line_arc_2;

  public:

  template < class OutputIterator>
  OutputIterator  
  operator()(const Circular_arc_2 &a, OutputIterator res) const
    { 
      static const bool Protection = true;      
      try{return internal::construct_bounding_hexagons_2<typename CK::AK>(a.approx(),res);}
      catch (Uncertain_conversion_exception)
      {
         CGAL::Protect_FPU_rounding<!Protection> P(CGAL_FE_TONEAREST);
         return internal::construct_bounding_hexagons_2<typename CK::EK>(a.exact(),res);
      }

   }

  Hexagon  operator()(const Line_arc_2 &a) const
    { 
      static const bool Protection = true;
      try{return internal::construct_bounding_hexagon_for_line_arc_2<typename CK::AK>(a.approx());}
      catch (Uncertain_conversion_exception)
      {
         CGAL::Protect_FPU_rounding<!Protection> P(CGAL_FE_TONEAREST);
         return internal::construct_bounding_hexagon_for_line_arc_2<typename CK::EK>(a.exact());
      }

    }
};


}// namespace CGAL

#endif // CGAL_HEXAGON_PRIMITIVES_H
