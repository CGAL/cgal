#ifndef CGAL_HEXAGON_PRIMITIVES_H
#define CGAL_HEXAGON_PRIMITIVES_H

#include <vector>
#include <CGAL/Polygon_2.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/number_utils.h>
#include <CGAL/Line_2.h>
#include <CGAL/Interval_arithmetic.h>
#include <CGAL/Object.h>
#include <CGAL/Curved_kernel/internal_functions_on_circular_arc_2.h>

namespace CGAL{
 namespace CGALi{

using CGAL::CircularFunctors::advanced_make_xy_monotone;


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
      assert(a->is_x_monotone());
      assert(a->is_y_monotone());

		
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
	
	
// Special orientation predicate used in to check to which side of a diagonal line ,of 
// degree 45 or 135 passing through the point a ,lies point b. The slope of the line
// (45 or 135) is indicated by the boolean argument. For now,is supposed only to serve
// the do_intersect predicate for the exagons, thus it serves a kernel encapsulating
// a lazy nt. 

CGAL::Comparison_result diagonal_orientation(Simple_cartesian<double>::Point_2 a, 
				             Simple_cartesian<double>::Point_2 b, bool dir)
  {	
    int i= ((dir)? -1 :1 );
    return CGAL::compare( to_interval(b.y()+i*b.x()),to_interval(a.y()+i*a.x()));
  }	
	







bool _do_intersect_hexagon_2(const CGAL::Polygon_2<Simple_cartesian< double > > &p1,
			     const CGAL::Polygon_2<Simple_cartesian<double> > &p2)
  {
    typedef Simple_cartesian<double> CK;
    typedef CK::RT      RT;
    typedef CK::Point_2 Point_2;
    typedef CK::Line_2      Line_2;
    typedef CGAL::Polygon_2<CK> Polygon_2;
		
    Polygon_2 a,b;
    bool pred,direction,tmp;
    CGAL::Comparison_result side;
    Point_2 frst,scnd;
    CGAL::Bbox_2 bb;

	
    int size1=p1.size(),
    	size2=p2.size();

    if(size1==4 && size2==4)
      return true;

    if(size1==4)
      {	
	a=p2;
	b=p1;
	bb=CGAL::Bbox_2(p2.left_vertex()->x(),p2.bottom_vertex()->y(), 
			p2.right_vertex()->x() ,p2.top_vertex()->y());

      } else if (size2==4)
      {
	a=p1;
	b=p2;
	bb=CGAL::Bbox_2(p1.left_vertex()->x(),p1.bottom_vertex()->y(), 
			p1.right_vertex()->x() ,p1.top_vertex()->y());
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
	  bb=test1;
	}else
	{
	  a=p2;
	  b=p1;
	  bb=CGAL::Bbox_2(p2.left_vertex()->x(),p2.bottom_vertex()->y(), 
			  p2.right_vertex()->x() ,p2.top_vertex()->y());
	}

      }

      pred=(a[0].x()!=a[1].x() && a[0].y()!=a[1].y());
	
      if(pred)
        {
	  bool temp;

          direction=(a[0].y()<a[1].y());
				
	  // Condition implied : is_on_upper_part
	  side= ( (  (a[1].y()==a[2].y())==direction )? LARGER : SMALLER);
	
	  for(int i=0;i<b.size();i++)
	    {
	      bool res=(diagonal_orientation(a[0],b[i],direction)==side);
	      if(temp=(!res && i==0))
	      break;
	      else if(!res)
	      return true;
	    }
	
	  if(!temp)
	    return false;

	 }
	
	 if(a[pred].y()==bb.ymax() || (pred && a[0].y()==bb.ymax()))
	   side=LARGER;
         else 
	   side=SMALLER;

	 for(int i=0;i<b.size();i++) 
	 {
	   Line_2 tmp_ln(a[2+pred],a[3+pred]);
	   std::pair<double,double> app_y=to_interval(tmp_ln.y_at_x(b[i].x()));
		
           if((side==SMALLER && b[i].y()<=app_y.second) ||(side==LARGER && b[i].y()>=app_y.first))
	     return true;	
	 }
			
	 return false;
			
}



// This is meant for x_monotone and for non-x_monotone cases
template < typename Hex_iterator>
bool do_intersect_hexagons_2(Hex_iterator a_begin, Hex_iterator a_end, Hex_iterator b_begin, Hex_iterator b_end)
  {

    for(Hex_iterator it1=a_begin;it1!=a_end;it1++)
      for(Hex_iterator it2=b_begin;it2!=b_end;it2++) 
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


 }//namespace CGALi

}// namespace CGAL

#endif // CGAL_HEXAGON_PRIMITIVES_H
