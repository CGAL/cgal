// Copyright (c) 2003  INRIA Sophia-Antipolis (France) and
//                     Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// Authors : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//           Sylvain Pion     <Sylvain.Pion@sophia.inria.fr>
//           Julien Hazebrouck
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (CGAL - Effective Computational Geometry for Curves and Surfaces) 

// file : include/CGAL/Line_arc_traits.h

#ifndef CGAL_VARIANT_TRAITS_H
#define CGAL_VARIANT_TRAITS_H

#include <CGAL/basic.h>
#include <cassert>
#include <boost/variant.hpp>
#include <CGAL/Curved_kernel/internal_functions_on_line_2.h>
#include <CGAL/Curved_kernel/internal_functions_on_circle_2.h>
#include <CGAL/Curved_kernel/internal_functions_on_circular_arc_2.h>
#include <CGAL/Curved_kernel/internal_functions_on_line_arc_2.h>


namespace CGAL {
  namespace VariantFunctors{

    // Takes an iterator range of Object(Line/Circular_arc),
    // returns an Object(Variant(Line/Circular_arc)).
    // Do nothing for Object(Endpoint).
    template <class CK, class Arc1, class Arc2, class OutputIterator>
      OutputIterator object_to_object_variant(const std::vector<CGAL::Object>& res1, OutputIterator res2){
      
      typedef typename CK::Circular_arc_point_2      Circular_arc_point_2;
      for(std::vector<CGAL::Object>::const_iterator it = res1.begin(); 
	  it != res1.end(); ++it ){
	if(const Arc1 *arc = CGAL::object_cast< Arc1 >(&*it)){
	  boost::variant< Arc1, Arc2 > v =  *arc;
	  *res2++ = make_object(v);
	}
	else if (const Arc2 *line = CGAL::object_cast< Arc2 >(&*it)){
	  boost::variant< Arc1, Arc2 > v =  *line;
	  *res2++ = make_object(v);
	}
	else{
	  *res2++ = *it;
	}
      }
      return res2;	  
    }


    template <class CurvedKernel, class Arc1, class Arc2>
    class In_range_2
    {
    public:
      typedef typename CurvedKernel::Circular_arc_point_2      Circular_arc_point_2;
      typedef bool result_type;
      
      result_type
	operator()(const boost::variant< Arc1, Arc2 > &a,
		   const Circular_arc_point_2 &p) const
      { 
	if ( const Arc1* arc1 = boost::get<Arc1>( &a ) ){
	  return CGAL::CircularFunctors::point_in_range<CurvedKernel>(*arc1, p);
	}
	else {
	  const Arc2* arc2 = boost::get<Arc2>( &a );
	  return CGAL::CircularFunctors::point_in_range<CurvedKernel>(*arc2, p);
	}
      }
    };
  

    template <class CurvedKernel, class Arc1, class Arc2>
    class Compare_y_to_right_2
    {
    public:
      typedef CGAL::Comparison_result result_type;
      typedef typename CurvedKernel::Circular_arc_point_2      Circular_arc_point_2;
    
      result_type
	operator()(const boost::variant< Arc1, Arc2 > &a1,
		   const boost::variant< Arc1, Arc2 > &a2,
		   const Circular_arc_point_2 &p) const
      { 
	if ( const Arc1* arc1 = boost::get<Arc1>( &a1 ) ){
	  if ( const Arc1* arc2 = boost::get<Arc1>( &a2 ) ){
	    return CGAL::CircularFunctors::compare_y_to_right<CurvedKernel>(*arc1, *arc2, p);
	  }
	  else {
	    const Arc2* arc2 = boost::get<Arc2>( &a2 );
	    return CGAL::CircularFunctors::compare_y_to_right<CurvedKernel>(*arc1, *arc2, p);
	  }
	}
	const Arc2* arc1 = boost::get<Arc2>( &a1 );
	if ( const Arc1* arc2 = boost::get<Arc1>( &a2 ) ){
	  return CGAL::CircularFunctors::compare_y_to_right<CurvedKernel>(*arc1, *arc2, p);
	}
	//else {
	const Arc2* arc2 = boost::get<Arc2>( &a2 );
	return CGAL::CircularFunctors::compare_y_to_right<CurvedKernel>(*arc1, *arc2, p);
	//}
      }  
    };


    template <class CurvedKernel>
    class Variant_Equal_2
      : public boost::static_visitor<bool>
    {
    public :

      template < typename T >
      bool
      operator()(const T &a0, const T &a1) const
      {
	return typename CurvedKernel::Equal_2()(a0, a1);
      }

      template < typename T1, typename T2 >      
      bool
      operator()(const T1 &, const T2 &) const
      {
	return false;
      }
    };

      

    template <class CurvedKernel, class Arc1, class Arc2>
    class Equal_2
      : public 
      CurvedKernel::Equal_2
    {
    public:
      typedef boost::variant< Arc1, Arc2 >  Curve_2;
      typedef bool result_type;

      using CurvedKernel::Equal_2::operator();

      result_type
	operator()(const Curve_2 &a0, const Curve_2 &a1) const
      {
	return boost::apply_visitor( Variant_Equal_2<CurvedKernel>(), a0, a1 );
      }
      
    };

    


    template <class CurvedKernel, class Arc1, class Arc2>
    class Compare_y_at_x_2
    {
    public:
      typedef typename CurvedKernel::Circular_arc_point_2      Circular_arc_point_2;
      typedef CGAL::Comparison_result result_type;

      result_type
	operator() (const Circular_arc_point_2 &p,
		    const boost::variant< Arc1, Arc2 > &A1) const
      { 
	if ( const Arc1* arc1 = boost::get<Arc1>( &A1 ) ){
	  return CGAL::CircularFunctors::compare_y_at_x<CurvedKernel>(p, *arc1); 
	}
	else {
	  const Arc2* arc2 = boost::get<Arc2>( &A1 );
	  return CGAL::CircularFunctors::compare_y_at_x<CurvedKernel>(p, *arc2);
	}
      }
    };

     template <class CurvedKernel>
    class Variant_Do_overlap_2
      : public boost::static_visitor<bool>
    {
    public :

      template < typename T >
      bool
      operator()(const T &a0, const T &a1) const
      {
	return typename CurvedKernel::Do_overlap_2()(a0, a1);
      }

      template < typename T1, typename T2 >      
      bool
      operator()(const T1 &, const T2 &) const
      {
	return false;
      }
    };
    

    template <class CurvedKernel, class Arc1, class Arc2>
    class Do_overlap_2
    {
    public:
      typedef typename CurvedKernel::Circular_arc_point_2      Circular_arc_point_2;
      typedef bool result_type;

      result_type
	operator()(const boost::variant< Arc1, Arc2 > &A0,
		   const boost::variant< Arc1, Arc2 > &A1) const
      { 
	return boost::apply_visitor( Variant_Do_overlap_2<CurvedKernel>(), A0, A1 );
      }    
    };


    template <class CurvedKernel, class Arc1, class Arc2>
    class Make_x_monotone_2
    {
    public:
      typedef typename CurvedKernel::Circular_arc_point_2      Circular_arc_point_2;

      template < class OutputIterator >
	OutputIterator
	operator()(const boost::variant< Arc1, Arc2 > &A, OutputIterator res)
      { if ( const Arc1* arc1 = boost::get<Arc1>( &A ) ){
	  std::vector<CGAL::Object> container;
	  CGAL::CircularFunctors::make_x_monotone<CurvedKernel>(*arc1,  std::back_inserter(container));
	  object_to_object_variant<CurvedKernel, Arc1, Arc2>(container, res);
	  return res;
	}
	else {
	  const Arc2* arc2 = boost::get<Arc2>( &A );
	  std::vector<CGAL::Object> container;
	  CGAL::CircularFunctors::make_x_monotone<CurvedKernel>(*arc2,  std::back_inserter(container));
	  object_to_object_variant<CurvedKernel, Arc1, Arc2>(container, res);
	  return res;
	}
      }
    };


    
    template <class CurvedKernel, class Arc1, class Arc2>
    class Intersect_2
    {
    public:
    typedef typename CurvedKernel::Circular_arc_point_2      Circular_arc_point_2;
      
      template < class OutputIterator >
	OutputIterator
	operator()(const boost::variant< Arc1, Arc2 > &c1,
		   const boost::variant< Arc1, Arc2 > &c2,
		   OutputIterator res) const
      { 
	if ( const Arc1* arc1 = boost::get<Arc1>( &c1 ) ){
	  if ( const Arc1* arc2 = boost::get<Arc1>( &c2 ) ){
	    std::vector<CGAL::Object> container;
	    CGAL::CircularFunctors::intersect_2<CurvedKernel>(*arc1,*arc2,  std::back_inserter(container));
	    object_to_object_variant<CurvedKernel, Arc1, Arc2>(container, res);
	    return res;
	  }
	  else if ( const Arc2* arc2 = boost::get<Arc2>( &c2 ) ){
	    std::vector<CGAL::Object> container;
	    CGAL::CircularFunctors::intersect_2<CurvedKernel>(*arc1, *arc2, std::back_inserter(container));
	    object_to_object_variant<CurvedKernel, Arc1, Arc2>(container, res);
	    return res;
	  }
	}
	else {
	  const Arc2* arc1 = boost::get<Arc2>( &c1 );
	  if ( const Arc1* arc2 = boost::get<Arc1>( &c2 ) ){
	    std::vector<CGAL::Object> container;
	    CGAL::CircularFunctors::intersect_2<CurvedKernel>(*arc1, *arc2,  std::back_inserter(container));
	    object_to_object_variant<CurvedKernel, Arc1, Arc2>(container, res);
	    return res;
	  }
	  const Arc2* arc2 = boost::get<Arc2>( &c2 );
	  std::vector<CGAL::Object> container;
	  CGAL::CircularFunctors::intersect_2<CurvedKernel>(*arc1, *arc2,  std::back_inserter(container));
	  object_to_object_variant<CurvedKernel, Arc1, Arc2>(container, res);
	  return res;
	}
	return res;//for no warning
      }
    
    };

    
    template <class CurvedKernel, class Arc1, class Arc2>
    class Split_2
    {

    public:
    typedef typename CurvedKernel::Circular_arc_point_2      Circular_arc_point_2;
      typedef void result_type;
      result_type
	operator()(const boost::variant< Arc1, Arc2 > &A, 
		   const Circular_arc_point_2 &p,
		   boost::variant< Arc1, Arc2 > &ca1,
		   boost::variant< Arc1, Arc2 > &ca2) const
      { 
	if ( const Arc1* arc1 = boost::get<Arc1>( &A ) ){
	  Arc1 carc1;
	  Arc1 carc2;
	  CGAL::CircularFunctors::split<CurvedKernel>(*arc1, p, carc1, carc2); 
	  ca1 = carc1;
	  ca2 = carc2;
	  return ;
	
	}
	else{
	  const Arc2* arc2 = boost::get<Arc2>( &A );
	  Arc2 cline1;
	  Arc2 cline2;
	  CGAL::CircularFunctors::split<CurvedKernel>(*arc2, p, cline1, cline2); 
	  ca1 = cline1;
	  ca2 = cline2;
	  return ;
	
	}
      }
    };


     template <class CurvedKernel>
    class Variant_Construct_min_vertex_2
      : public boost::static_visitor<const typename CurvedKernel::Circular_arc_point_2&>
    {
    public :

      template < typename T >
      const typename CurvedKernel::Circular_arc_point_2&
      operator()(const T &a) const
      {
	bool precondition = typename CurvedKernel::Compare_xy_2()(a.left(),a.right())==CGAL::SMALLER;
	CGAL_kernel_precondition(precondition);
	return (a.left());
      }
    };

    template <class CurvedKernel, class Arc1, class Arc2>
    class Construct_min_vertex_2
    {
      typedef typename CurvedKernel::Circular_arc_point_2      Point_2;
    public:
      const Point_2& operator() (const boost::variant< Arc1, Arc2 > & cv) const
      {
	return boost::apply_visitor( Variant_Construct_min_vertex_2<CurvedKernel>(), cv );
      }
    };





     template <class CurvedKernel>
    class Variant_Construct_max_vertex_2
       : public boost::static_visitor<const typename CurvedKernel::Circular_arc_point_2&>
    {
    public :

      template < typename T >
      const typename CurvedKernel::Circular_arc_point_2&
      operator()(const T &a) const
      {
	bool precondition = typename CurvedKernel::Compare_xy_2()(a.left(),a.right())==CGAL::SMALLER;
	CGAL_kernel_precondition(precondition); 
	return (a.right());
      }
    };


    
    template <class CurvedKernel, class Arc1, class Arc2>
    class Construct_max_vertex_2
    {
      typedef typename CurvedKernel::Circular_arc_point_2      Point_2;
    public:
      /*!
       * Get the right endpoint of the x-monotone curve (segment).
       * \param cv The curve.
       * \return The right endpoint.
       */
      const Point_2& operator() (const boost::variant< Arc1, Arc2 > & cv) const
      {
	return boost::apply_visitor( Variant_Construct_max_vertex_2<CurvedKernel>(), cv );
      }
    };

      template <class CurvedKernel>
    class Variant_Is_vertical_2
      : public boost::static_visitor<bool>
    {
    public :

      template < typename T >
      bool
      operator()(const T &a) const
      {
	return CGAL::CircularFunctors::is_vertical<CurvedKernel>(a);
      }
    };

    template <class CurvedKernel, class Arc1, class Arc2>
    class Is_vertical_2
    {
    public:
      typedef bool result_type;

      bool operator() (const boost::variant< Arc1, Arc2 >& cv) const
      {
	return boost::apply_visitor( Variant_Is_vertical_2<CurvedKernel>(), cv );
      }
    };

  }
  /// Traits class for CGAL::Arrangement_2 (and similar) based on a CurvedKernel.

    template < typename CurvedKernel, typename Arc1, typename Arc2>
    class Variant_traits {

    typedef Variant_traits< CurvedKernel, Arc1, Arc2 >   Self;

  public:
  
    typedef CurvedKernel Kernel;
    typedef typename CurvedKernel::Circular_arc_point_2      Circular_arc_point_2;

    typedef typename CurvedKernel::Circular_arc_point_2      Point;
    typedef typename CurvedKernel::Circular_arc_point_2      Point_2;
  
    typedef CGAL::Tag_false                        Has_left_category;
    typedef CGAL::Tag_false 			 Has_merge_category;
  
    typedef boost::variant< Arc1, Arc2 > Curve_2;
    typedef boost::variant< Arc1, Arc2 > X_monotone_curve_2;

  private:
    CurvedKernel ck;
  public:
  
    Variant_traits(const CurvedKernel &k = CurvedKernel())
      : ck(k) {}

    typedef typename CurvedKernel::Compare_x_2           Compare_x_2;
    typedef typename CurvedKernel::Compare_xy_2          Compare_xy_2;
    typedef VariantFunctors::Construct_min_vertex_2<CurvedKernel, Arc1, Arc2>  Construct_min_vertex_2;
    typedef VariantFunctors::Construct_max_vertex_2<CurvedKernel, Arc1, Arc2>  Construct_max_vertex_2;
    typedef VariantFunctors::Is_vertical_2<CurvedKernel, Arc1, Arc2>           Is_vertical_2;
    typedef VariantFunctors::Compare_y_at_x_2<CurvedKernel, Arc1, Arc2>      Compare_y_at_x_2;
    typedef VariantFunctors::Compare_y_to_right_2<CurvedKernel, Arc1, Arc2>  Compare_y_at_x_right_2; 
    typedef VariantFunctors::Equal_2<CurvedKernel, Arc1, Arc2>               Equal_2;
    typedef VariantFunctors::Make_x_monotone_2<CurvedKernel, Arc1, Arc2>     Make_x_monotone_2;
    typedef VariantFunctors::Split_2<CurvedKernel, Arc1, Arc2>               Split_2;
    typedef VariantFunctors::Intersect_2<CurvedKernel, Arc1, Arc2> Intersect_2;

  
 Compare_x_2 compare_x_2_object() const
  { return ck.compare_x_2_object(); }

  Compare_xy_2 compare_xy_2_object() const
  { return ck.compare_xy_2_object(); }

  Compare_y_at_x_2 compare_y_at_x_2_object() const 
  { return Compare_y_at_x_2(); }

  Compare_y_at_x_right_2 compare_y_at_x_right_2_object() const 
  { return Compare_y_at_x_right_2(); }

  Equal_2 equal_2_object() const
  { return Equal_2(); }

  Make_x_monotone_2 make_x_monotone_2_object() const
  { return Make_x_monotone_2(); }

  Split_2 split_2_object() const
  { return Split_2(); }

  Intersect_2 intersect_2_object() const
    { return Intersect_2(); }

  Construct_min_vertex_2 construct_min_vertex_2_object() const
    { return Construct_min_vertex_2(); }

  Construct_max_vertex_2 construct_max_vertex_2_object() const
    { return Construct_max_vertex_2(); }

  Is_vertical_2 is_vertical_2_object() const
    { return Is_vertical_2();}


};

} // namespace CGAL

#endif // CGAL_VARIANT_TRAITS_H
