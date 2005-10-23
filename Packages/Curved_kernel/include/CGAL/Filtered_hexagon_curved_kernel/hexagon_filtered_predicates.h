// Author : Constantinos Tsirogiannis

#ifndef CGAL_HEXAGON_FILTERED_PREDICATES_H  
#define CGAL_HEXAGON_FILTERED_PREDICATES_H  


#include <cassert>
#include <fstream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h> 
#include <CGAL/enum.h>
#include <CGAL/Object.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Interval_arithmetic.h>
#include <CGAL/Filtered_hexagon_curved_kernel/hexagon_primitives.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {


template <class HK>
class In_range_2
  {
    typedef typename HK::Curved_kernel                           CK;
    typedef typename HK::Circular_arc_point_2                    Circular_arc_point_2;
    typedef typename HK::Circular_arc_2                          Circular_arc_2;
    typedef typename HK::Line_arc_2                              Line_arc_2;

   public:

    typedef bool result_type;

   private:

    template <class Arc_2>
    result_type
    _in_range_2(const Arc_2 &a, const Circular_arc_point_2 &p) const 
    {
      std::pair<double,double> pr;
        
          if(a.has_no_hexagons())
           a.construct_hexagons();
      
      typename Arc_2::Hexagon_const_iterator hips;

         pr = to_interval(p.x());

      hips=a.hexagons_begin();

      while(hips!=a.hexagons_end())
	{ 
	  if( (pr.second >= (*(*hips).left_vertex()).x() &&  
               pr.second <= (*(*hips).right_vertex()).x() ) ||
	      (pr.first >= (*(*hips).left_vertex()).x() &&
               pr.first <= (*(*hips).right_vertex()).x() )  )	
                return CK().in_range_2_object()( a.arc() ,p);

	  hips++;			
	}

      return false;
     }

   public:

    result_type
    operator()( const Circular_arc_2 &a, const Circular_arc_point_2 &p) const
    { 
      CGAL_precondition( a.arc().is_x_monotone());
      return _in_range_2(a,p);
   }

    result_type
    operator()( const Line_arc_2 &a, const Circular_arc_point_2 &p) const
    { return _in_range_2(a,p);}
    

  };


template <class HK>
class Construct_Circular_min_vertex_2
  {
    typedef typename HK::Curved_kernel                       CK;
    typedef typename HK::Circular_arc_point_2                Circular_arc_point_2;
    typedef typename HK::Circular_arc_2                      Circular_arc_2;

   public:

    typedef Circular_arc_point_2 result_type;
    
    template <typename T>
    result_type
      operator()(const T& a) const
    { return a.left();}
    
  };

template <class HK>
class Construct_Circular_max_vertex_2
  {
    typedef typename HK::Curved_kernel                           CK;
    typedef typename HK::Circular_arc_point_2                 Circular_arc_point_2;
    typedef typename HK::Circular_arc_2                          Circular_arc_2;

   public:

    typedef Circular_arc_point_2 result_type;

    template <typename T>
     result_type
      operator()(const T& a) const
    { return a.right();  }

  };

template <class HK>
class Is_vertical_2
  {
    typedef typename HK::Curved_kernel                           CK;
    typedef typename HK::Circular_arc_point_2                 Circular_arc_point_2;
    typedef typename HK::Circular_arc_2                          Circular_arc_2;

   public:

    typedef bool result_type;
    
    template <typename T>
    result_type
      operator()(const T& a) const
    { return CK().is_vertical_2_object()(a.arc()); }
  };

 


template <class HK>
class Compare_y_at_x_2
  {
    typedef typename HK::Curved_kernel                                    CK;
    typedef typename HK::Circular_arc_2                                   Circular_arc_2;
    typedef typename HK::Circular_arc_point_2                             Circular_arc_point_2;
    typedef typename HK::Line_arc_2                                       Line_arc_2;
    typedef typename Circular_arc_2::Hexagon                              Hexagon;
    typedef typename Simple_cartesian<double>::Point_2		          Point_2; //Attention!!!
    typedef Exact_predicates_inexact_constructions_kernel                 EK;

  public:

    typedef Comparison_result result_type;

  private:

    template <class Arc_2>
    result_type
    _compare_y_at_x_2(const Circular_arc_point_2 &p,const Arc_2 &a) const
    {
      CGAL_precondition_code(bool tmp=In_range_2<HK>()(a,p));
      CGAL_precondition(tmp );

      if(a.has_no_hexagons())
        a.construct_hexagons();
	
      Comparison_result rel_pos=EQUAL;
      Bbox_2 bb = p.bbox();

      Hexagon     pnt_vec;

      pnt_vec.push_back(Point_2(bb.xmin(),bb.ymin()));
      pnt_vec.push_back(Point_2(bb.xmin(),bb.ymax()));
      pnt_vec.push_back(Point_2(bb.xmax(),bb.ymax()));
      pnt_vec.push_back(Point_2(bb.xmax(),bb.ymin()));

      typename Arc_2::Hexagon_const_iterator hips=a.hexagons_begin();

      while(hips!=a.hexagons_end())
	{
	  if( !(bb.xmin() > hips->right_vertex()->x() || bb.xmax() < hips->left_vertex()->x())  )
	    if(bb.ymin() > hips->top_vertex()->y()) 
	      rel_pos=LARGER;
	    else if(bb.ymax() < hips->bottom_vertex()->y()) 
	      rel_pos=SMALLER;
	    else
              {
		rel_pos=EQUAL;
		Hexagon hxgn= *hips;
		if(_do_intersect_hexagon_2(hxgn,pnt_vec))
                   return CK().compare_y_at_x_2_object()(p,a.arc());
                  

	        EK::Point_2 a1(hxgn[0].x(),hxgn[0].y()),
			    a2(hxgn[1].x(),hxgn[1].y()),
			    a3(pnt_vec[0].x(),pnt_vec[0].y());
		
		if ( EK().orientation_2_object()(a1,a2,a3)==RIGHT_TURN)
		  rel_pos=SMALLER;
		else 
		  rel_pos=LARGER;
	
	      }

	      hips++;
	    }

	    return rel_pos;
    }

  public:

    result_type
    operator()( const Circular_arc_point_2 &p,const Circular_arc_2 &a ) const
    {
      CGAL_kernel_precondition( a.arc().is_x_monotone());
      return _compare_y_at_x_2(p,a);
    }

    result_type
    operator()( const Circular_arc_point_2 &p,const Line_arc_2 &a ) const
     {return _compare_y_at_x_2(p,a);}

};



template <class HK>
class Equal_2
  {
    typedef typename HK::Curved_kernel                                    CK;
    typedef typename HK::Circular_arc_2                                   Circular_arc_2;
    typedef typename HK::Circular_arc_point_2                             Circular_arc_point_2;
    typedef typename HK::Line_arc_2                                       Line_arc_2;

  public:

    typedef bool result_type;

  private:

    template <class Arc_2>
    result_type
    _equal_2(const Arc_2 &a,const Arc_2 &b) const
    {
      if(a.has_no_hexagons())
        a.construct_hexagons();

	if(b.has_no_hexagons())
        b.construct_hexagons();
	
      if(a.hexagon_no() != b.hexagon_no())
	  return false;

      typename Arc_2::Hexagon_const_iterator hips=a.hexagons_begin(),
						  bips=b.hexagons_begin();

      while(hips!=a.hexagons_end())
      {	
	if(hips->size()!=bips->size())
	  return false;

	for(int j=0;j< hips->size() ;j++)
	  if( (*hips)[j]!=(*bips)[j] )
	    return false;

	hips++;
	bips++;
      }

      return CK().equal_2_object()( a.arc(),b.arc() );

    }

  public:

    result_type
    operator()( const Circular_arc_point_2 &a ,
                const Circular_arc_point_2 &b ) const
	 {  return CK().equal_2_object()(a,b);}

    result_type
    operator()( const Circular_arc_2 &a , const Circular_arc_2 &b ) const
    {
      CGAL_precondition( a.arc().is_x_monotone());
      CGAL_precondition( b.arc().is_x_monotone());
      return _equal_2(a,b);      
    }

    result_type
    operator()( const Line_arc_2 &a ,
                const Line_arc_2 &b ) const
    {  return _equal_2(a,b);}

    result_type
    operator()( const Circular_arc_2 &a ,
                const Line_arc_2 &b ) const
	 {  return false;}

    result_type
    operator()( const Line_arc_2 &a ,
                const Circular_arc_2 &b ) const
	 {  return false;}

};


template <class HK>
class Do_overlap_2
  {
    typedef typename HK::Curved_kernel                                  CK;
    typedef typename HK::Circular_arc_2                                 Circular_arc_2;
    typedef typename HK::Line_arc_2                                     Line_arc_2;

  public:

    typedef bool result_type;

  private:

    template <class Arc_2>
    result_type
    _do_overlap_2(const Arc_2 &a, const Arc_2 &b) const
    {
      if(a.has_no_hexagons())
        a.construct_hexagons();
	
      if(b.has_no_hexagons())
        b.construct_hexagons();

      if(do_intersect_hexagons_2( a.hexagons_begin(),a.hexagons_end(),b.hexagons_begin(),b.hexagons_end() ) )
        return CK().do_overlap_2_object()(a.arc(),b.arc());

      return false;
    }


  public:
    
    result_type
    operator()( const Circular_arc_2 &a , const Circular_arc_2 &b ) const
    {
      CGAL_precondition( a.arc().is_x_monotone());
      CGAL_precondition( b.arc().is_x_monotone());
      return _do_overlap_2(a,b); 
    }

    result_type
    operator()( const Line_arc_2 &a ,
                const Line_arc_2 &b ) const
    {  return _do_overlap_2(a,b);}

    result_type
    operator()( const Circular_arc_2 &a ,
                const Line_arc_2 &b ) const
	 {  return false;}

    result_type
    operator()( const Line_arc_2 &a ,
                const Circular_arc_2 &b ) const
	 {  return false;}

};


// This predicate cannot be filtered

 template < class HK >
  class Compare_y_to_right_2
  {
    typedef typename HK::Curved_kernel            CK;
    typedef typename HK::Circular_arc_2           Circular_arc_2;
    typedef typename HK::Circular_arc_point_2  Circular_arc_point_2;

  public:
    typedef Comparison_result result_type;
    
    template <typename T1, typename T2>
    result_type
    operator()(const T1 &a1,
               const T2 &a2,
               const Circular_arc_point_2 &p) const
    { return CK().compare_y_to_right_2_object()(a1.arc(), a2.arc(), p); }
    
  };





  template < class HK >
  class Make_x_monotone_2
  {
    typedef typename HK::Curved_kernel            CK;
    typedef typename HK::Rcirc_arc_2              Rcirc_arc_2;
    typedef typename HK::Circular_arc_2       	  Circular_arc_2;
    typedef typename HK::Line_arc_2               Line_arc_2;

  public:
    // typedef OutputIterator result_type;

    template < class OutputIterator >
    OutputIterator
    operator()(const Circular_arc_2 &A, OutputIterator res)
      { 
	std::vector<Object> vec;
	
	CK().make_x_monotone_2_object()(A.arc(), std::back_inserter(vec) );
 
	for(unsigned i=0; i<vec.size() ; i++)
	{
	  const Rcirc_arc_2 *tmp_arc;

	  tmp_arc=object_cast<Rcirc_arc_2>(&vec.at(i));
		
	  assert(tmp_arc!=NULL);
		
	  *res++ = make_object( Circular_arc_2(*tmp_arc) );
	}

	return res;
      }
    
    template < class OutputIterator >
    OutputIterator
    operator()(const Line_arc_2 &A, OutputIterator res)
      { 
	*res++ = make_object(A);
	return res;
      }
  };


  template < class HK >
  class Intersect_2
  {
    public:

    typedef typename HK::Curved_kernel            CK;
    typedef typename HK::Circular_arc_2           Circular_arc_2;
    typedef typename HK::Rcirc_arc_2              Rcirc_arc_2;
    typedef typename HK::Rline_arc_2              Rline_arc_2;
    typedef typename HK::Circle_2                 Circle;
    typedef typename HK::Line_arc_2               Line_arc_2;

    template < class OutputIterator >
    OutputIterator
    operator()(const Circle & c1, const Circle & c2, OutputIterator res)
      { return CK().intersect_2_object()(c1,c2,res); }

     template < class OutputIterator >
    OutputIterator
    operator()(const Circular_arc_2 & c1, const Circular_arc_2 & c2, 
	       OutputIterator res)
      { 
        if(c1.has_no_hexagons())
          c1.construct_hexagons();

        if(c2.has_no_hexagons())
          c2.construct_hexagons();

        if(!do_intersect_hexagons_2(c1.hexagons_begin(),c1.hexagons_end(),c2.hexagons_begin(),c2.hexagons_end()) )
         return res;

	std::vector<Object> vec;
	
	CK().intersect_2_object()(c1.arc(),c2.arc(),std::back_inserter(vec)); 

 	for(unsigned i=0; i<vec.size() ; i++)
	{
	  const Rcirc_arc_2 *tmp_arc;

	  if ( (tmp_arc=object_cast<Rcirc_arc_2>(&vec.at(i)) )!=NULL )
	    *res++ = make_object( Circular_arc_2(*tmp_arc) );
	  else
	    *res++ = vec.at(i);
	}
	
	return res;
	

      }

     template < class OutputIterator >
    OutputIterator
    operator()(const Line_arc_2 & c1, const Line_arc_2 & c2, 
	       OutputIterator res)
      { 
        if(c1.has_no_hexagons())
          c1.construct_hexagons();

        if(c2.has_no_hexagons())
          c2.construct_hexagons();

        if(!do_intersect_hexagons_2(c1.hexagons_begin(),c1.hexagons_end(),c2.hexagons_begin(),c2.hexagons_end()) )
         return res;


	std::vector<Object> vec;
	
	CK().intersect_2_object()(c1.arc(),c2.arc(),std::back_inserter(vec)); 

 	for(unsigned i=0; i<vec.size() ; i++)
	{
	  const Rline_arc_2 *tmp_arc;

	  if ( (tmp_arc=object_cast<Rline_arc_2>(&vec.at(i)) )!=NULL )
	    *res++ = make_object( Line_arc_2(*tmp_arc) );
	  else
	    *res++ = vec.at(i);
	}
	
	return res;	
      }

     template < class OutputIterator >
    OutputIterator
    operator()(const Circular_arc_2 & c1, const Line_arc_2 & c2, 
	       OutputIterator res)
      { 
        if(c1.has_no_hexagons())
          c1.construct_hexagons();

        if(c2.has_no_hexagons())
          c2.construct_hexagons();

        if(!do_intersect_hexagons_2(c1.hexagons_begin(),c1.hexagons_end(),c2.hexagons_begin(),c2.hexagons_end()) )
         return res;


	std::vector<Object> vec;
	
	CK().intersect_2_object()(c1.arc(),c2.arc(),std::back_inserter(vec)); 

 	for(unsigned i=0; i<vec.size() ; i++)
	{
	  const Rcirc_arc_2 *tmp_arc;
	  const Rline_arc_2 *tmp_line;

	  if ( (tmp_arc=object_cast<Rcirc_arc_2>(&vec.at(i)) )!=NULL )
	    *res++ = make_object( Circular_arc_2(*tmp_arc) );
	  else if ( (tmp_line=object_cast<Rline_arc_2>(&vec.at(i)) )!=NULL )
	    *res++ = make_object( Line_arc_2(*tmp_line) );
	  else
	    *res++ = vec.at(i);
	}

	return res;	

      }

      template < class OutputIterator >
    OutputIterator
    operator()(const Line_arc_2 & c1, const Circular_arc_2 & c2, 
	       OutputIterator res)
      {
	return operator()(c2,c1,res);
      }
     
    
  };


  template < class HK >
  class Split_2
  {

    typedef typename HK::Curved_kernel            CK;
    typedef typename HK::Rcirc_arc_2              Rcirc_arc_2;
    typedef typename HK::Rline_arc_2              Rline_arc_2;
    typedef typename HK::Circular_arc_2           Circular_arc_2;
    typedef typename HK::Line_arc_2               Line_arc_2;
    typedef typename HK::Circular_arc_point_2  Circular_arc_point_2;

  public:
    typedef void result_type;

    result_type
    operator()(const Circular_arc_2 &A, 
	       const Circular_arc_point_2 &p,
	       Circular_arc_2 &ha1, Circular_arc_2 &ha2) const
    {  
      Rcirc_arc_2 ca1 , ca2;

      CK().split_2_object()(A.arc(), p, ca1, ca2);

      ha1=Circular_arc_2(ca1); 
      ha2=Circular_arc_2(ca2); 
    }
    
    result_type
    operator()(const Line_arc_2 &A, 
	       const Circular_arc_point_2 &p,
	       Line_arc_2 &ha1, Line_arc_2 &ha2) const
    {  
      Rline_arc_2 ca1 , ca2;

      CK().split_2_object()(A.arc(), p, ca1, ca2);

      ha1=Line_arc_2(ca1); 
      ha2=Line_arc_2(ca2); 
    }
    
  };


  
} //CGALi

CGAL_END_NAMESPACE  

#endif // CGAL_HEXAGON_FILTERED_PREDICATES_H  
