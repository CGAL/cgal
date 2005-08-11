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
    template <class CK,class OutputIterator >
      OutputIterator object_to_object_variant(const std::vector<CGAL::Object>& res1, OutputIterator res2){
      typedef typename CK::Line_arc_2  Line_arc_2;
      typedef typename CK::Circular_arc_2  Circular_arc_2;
      typedef typename CK::Circular_arc_endpoint_2      Circular_arc_endpoint_2;
      for(std::vector<CGAL::Object>::const_iterator it = res1.begin(); 
	  it != res1.end(); ++it ){
	if(const Circular_arc_2 *arc = CGAL::object_cast< Circular_arc_2 >(&*it)){
	  boost::variant< Circular_arc_2, Line_arc_2 > v =  *arc;
	  *res2++ = make_object(v);
	}
	else if (const Line_arc_2 *line = CGAL::object_cast< Line_arc_2 >(&*it)){
	  boost::variant< Circular_arc_2, Line_arc_2 > v =  *line;
	  *res2++ = make_object(v);
	}
	else{
	  *res2++ = *it;
	}
      }
      return res2;	  
    }

    template <class CurvedKernel>
    class In_range_2
    {
    public:
      typedef typename CurvedKernel::Line_arc_2  Line_arc_2;
      typedef typename CurvedKernel::Circular_arc_2  Circular_arc_2;
      typedef typename CurvedKernel::Circular_arc_endpoint_2      Circular_arc_endpoint_2;
      typedef bool result_type;
      
      result_type
	operator()(const boost::variant< Circular_arc_2, Line_arc_2 > &a,
		   const Circular_arc_endpoint_2 &p) const
      { 
	if ( const Circular_arc_2* arc1 = boost::get<Circular_arc_2>( &a ) ){
	  return CGAL::CircularFunctors::point_in_range<CurvedKernel>(*arc1, p);
	}
	else {
	  const Line_arc_2* line1 = boost::get<Line_arc_2>( &a );
	  return CGAL::CircularFunctors::point_in_range<CurvedKernel>(*line1, p);
	}
      }
    };
  

    template <class CurvedKernel>
    class Compare_y_to_right_2
    {
    public:
      typedef CGAL::Comparison_result result_type;
      typedef typename CurvedKernel::Line_arc_2  Line_arc_2;
      typedef typename CurvedKernel::Circular_arc_2  Circular_arc_2;
      typedef typename CurvedKernel::Circular_arc_endpoint_2      Circular_arc_endpoint_2;
    
      result_type
	operator()(const boost::variant< Circular_arc_2, Line_arc_2 > &a1,
		   const boost::variant< Circular_arc_2, Line_arc_2 > &a2,
		   const Circular_arc_endpoint_2 &p) const
      { 
	if ( const Circular_arc_2* arc1 = boost::get<Circular_arc_2>( &a1 ) ){
	  if ( const Circular_arc_2* arc2 = boost::get<Circular_arc_2>( &a2 ) ){
	    return CGAL::CircularFunctors::compare_y_to_right<CurvedKernel>(*arc1, *arc2, p);
	  }
	  else {
	    const Line_arc_2* line2 = boost::get<Line_arc_2>( &a2 );
	    if (CGAL::CircularFunctors::compare_y_to_right<CurvedKernel>(*line2, *arc1, p) == CGAL::LARGER)
	      return CGAL::SMALLER;
	    return CGAL::LARGER;
	  }
	}
	const Line_arc_2* line1 = boost::get<Line_arc_2>( &a1 );
	if ( const Circular_arc_2* arc2 = boost::get<Circular_arc_2>( &a2 ) ){
	  return CGAL::CircularFunctors::compare_y_to_right<CurvedKernel>(*line1, *arc2, p);
	}
	//else {
	const Line_arc_2* line2 = boost::get<Line_arc_2>( &a2 );
	return CGAL::CircularFunctors::compare_y_to_right<CurvedKernel>(*line1, *line2, p);
	//}
      }  
    };


    template <class CurvedKernel>
    class Equal_2
      : public //CurvedKernel::Linear_kernel::Equal_2,
      CurvedKernel::Equal_2
    {
    public:
      //       typedef typename VariantTraits::Curve_2  Curve_2;
      typedef typename CurvedKernel::Line_arc_2   C1;
      typedef typename CurvedKernel::Circular_arc_2  C2;
      typedef boost::variant< C2, C1 >  Curve_2;
      typedef bool result_type;

      //using CurvedKernel::Linear_kernel::Equal_2::operator();
      using CurvedKernel::Equal_2::operator();

      result_type
	operator()(const Curve_2 &a0, const Curve_2 &a1) const
      { 
	if ( const C1* arc1 = boost::get<C1>( &a0 ) ){
	  if ( const C1* arc2 = boost::get<C1>( &a1 ) ){
	    return this->operator()(*arc1, *arc2);
	  }
	  CGAL_assertion(boost::get<C2>( &a1 ) );
	  return false;
	}
	else {
	  const C2* line1 = boost::get<C2>( &a0 );
	  if ( boost::get<C1>( &a1 ) ){
	    return false;
	  }
	  //else {
	  const C2* line2 = boost::get<C2>( &a1 );
	  return this->operator()(*line1, *line2);
	  //}
	}
	return false;//for no warning
      }
    };


    template <class CurvedKernel>
    class Compare_y_at_x_2
    {
    public:
      typedef typename CurvedKernel::Line_arc_2  Line_arc_2;
      typedef typename CurvedKernel::Circular_arc_2  Circular_arc_2;
      typedef typename CurvedKernel::Circular_arc_endpoint_2      Circular_arc_endpoint_2;
      typedef CGAL::Comparison_result result_type;

      result_type
	operator() (const Circular_arc_endpoint_2 &p,
		    const boost::variant< Circular_arc_2, Line_arc_2 > &A1) const
      { 
	if ( const Circular_arc_2* arc1 = boost::get<Circular_arc_2>( &A1 ) ){
	  return CGAL::CircularFunctors::compare_y_at_x<CurvedKernel>(p, *arc1); 
	}
	else {
	  const Line_arc_2* line1 = boost::get<Line_arc_2>( &A1 );
	  return CGAL::CircularFunctors::compare_y_at_x<CurvedKernel>(p, *line1);
	}
      }
    };

    template <class CurvedKernel>
    class Do_overlap_2
    {
    public:
      typedef typename CurvedKernel::Line_arc_2  Line_arc_2;
      typedef typename CurvedKernel::Circular_arc_2  Circular_arc_2;
      typedef typename CurvedKernel::Circular_arc_endpoint_2      Circular_arc_endpoint_2;
      typedef bool result_type;

      result_type
	operator()(const boost::variant< Circular_arc_2, Line_arc_2 > &A0,
		   const boost::variant< Circular_arc_2, Line_arc_2 > &A1) const
      { 
	if ( const Circular_arc_2* arc1 = boost::get<Circular_arc_2>( &A0 ) ){
	  if ( const Circular_arc_2* arc2 = boost::get<Circular_arc_2>( &A1 ) ){
	    return CGAL::CircularFunctors::do_overlap<CurvedKernel>(*arc1, *arc2);
	  }
	  else if ( const Line_arc_2* line2 = boost::get<Line_arc_2>( &A1 ) ){
	    return false;
	  }
	}
	else {
	  const Line_arc_2* line1 = boost::get<Line_arc_2>( &A0 );
	  if ( const Circular_arc_2* arc2 = boost::get<Circular_arc_2>( &A1 ) ){
	    return false;
	  }
	  //else {
	  const Line_arc_2* line2 = boost::get<Line_arc_2>( &A1 );
	  return CGAL::CircularFunctors::do_overlap<CurvedKernel>(*line1, *line2);
	  //}
	}
      }    
    };

    template <class CurvedKernel>
    class Make_x_monotone_2
    {
    public:
      typedef typename CurvedKernel::Line_arc_2  Line_arc_2;
      typedef typename CurvedKernel::Circular_arc_2  Circular_arc_2;
      typedef typename CurvedKernel::Circular_arc_endpoint_2      Circular_arc_endpoint_2;

      template < class OutputIterator >
	OutputIterator
	operator()(const boost::variant< Circular_arc_2, Line_arc_2 > &A, OutputIterator res)
      { if ( const Circular_arc_2* arc1 = boost::get<Circular_arc_2>( &A ) ){
	  std::vector<CGAL::Object> container;
	  CGAL::CircularFunctors::make_x_monotone<CurvedKernel>(*arc1,  std::back_inserter(container));
	  object_to_object_variant<CurvedKernel>(container, res);
	  return res;
	}
	else {
	  *res++ = make_object(A);
	  return res;
	}
      }
    };


    
    template <class CurvedKernel>
    class Construct_intersections_2
    {
    public:
      typedef typename CurvedKernel::Line_arc_2  Line_arc_2;
    typedef typename CurvedKernel::Circular_arc_2  Circular_arc_2;
    typedef typename CurvedKernel::Circular_arc_endpoint_2      Circular_arc_endpoint_2;
      
      template < class OutputIterator >
	OutputIterator
	operator()(const boost::variant< Circular_arc_2, Line_arc_2 > &c1,
		   const boost::variant< Circular_arc_2, Line_arc_2 > &c2,
		   OutputIterator res) const
      { 
	if ( const Circular_arc_2* arc1 = boost::get<Circular_arc_2>( &c1 ) ){
	  if ( const Circular_arc_2* arc2 = boost::get<Circular_arc_2>( &c2 ) ){
	    std::vector<CGAL::Object> container;
	    CGAL::CircularFunctors::construct_intersections_2<CurvedKernel>(*arc1,*arc2,  std::back_inserter(container));
	    object_to_object_variant<CurvedKernel>(container, res);
	    return res;
	  }
	  else if ( const Line_arc_2* line2 = boost::get<Line_arc_2>( &c2 ) ){
	    std::vector<CGAL::Object> container;
	    CGAL::CircularFunctors::construct_intersections_2<CurvedKernel>(*line2, *arc1,  std::back_inserter(container));
	    object_to_object_variant<CurvedKernel>(container, res);
	    return res;
	  }
	}
	else {
	  const Line_arc_2* line1 = boost::get<Line_arc_2>( &c1 );
	  if ( const Circular_arc_2* arc2 = boost::get<Circular_arc_2>( &c2 ) ){
	    std::vector<CGAL::Object> container;
	    CGAL::CircularFunctors::construct_intersections_2<CurvedKernel>(*line1, *arc2,  std::back_inserter(container));
	    object_to_object_variant<CurvedKernel>(container, res);
	    return res;
	  }
	  //else {
	  const Line_arc_2* line2 = boost::get<Line_arc_2>( &c2 );
	  std::vector<CGAL::Object> container;
	  CGAL::CircularFunctors::construct_intersections_2<CurvedKernel>(*line1, *line2,  std::back_inserter(container));
	   object_to_object_variant<CurvedKernel>(container, res);
	  return res;
	  //}
	}
	return res;//for no warning
      }
    
    };

    
    template <class CurvedKernel>
    class Split_2
    {

    public:
      typedef typename CurvedKernel::Line_arc_2  Line_arc_2;
    typedef typename CurvedKernel::Circular_arc_2  Circular_arc_2;
    typedef typename CurvedKernel::Circular_arc_endpoint_2      Circular_arc_endpoint_2;
      typedef void result_type;
      result_type
	operator()(const boost::variant< Circular_arc_2, Line_arc_2 > &A, 
		   const Circular_arc_endpoint_2 &p,
		   boost::variant< Circular_arc_2, Line_arc_2 > &ca1,
		   boost::variant< Circular_arc_2, Line_arc_2 > &ca2) const
      { 
	if ( const Circular_arc_2* arc1 = boost::get<Circular_arc_2>( &A ) ){
	  Circular_arc_2 carc1;
	  Circular_arc_2 carc2;
	  CGAL::CircularFunctors::split<CurvedKernel>(*arc1, p, carc1, carc2); 
	  ca1 = carc1;
	  ca2 = carc2;
	  return ;
	
	}
	else{
	  const Line_arc_2* line1 = boost::get<Line_arc_2>( &A );
	  Line_arc_2 cline1;
	  Line_arc_2 cline2;
	  CGAL::CircularFunctors::split<CurvedKernel>(*line1, p, cline1, cline2); 
	  ca1 = cline1;
	  ca2 = cline2;
	  return ;
	
	}
      }
    };



  }
  /// Traits class for CGAL::Arrangement_2 (and similar) based on a CurvedKernel.

  template < typename CurvedKernel >
    class Variant_traits {

    typedef Variant_traits< CurvedKernel >   Self;

  public:
  
    typedef CurvedKernel Kernel;
    typedef typename CurvedKernel::Line_arc_2  Line_arc_2;
    typedef typename CurvedKernel::Circular_arc_2  Circular_arc_2;
    typedef typename CurvedKernel::Circular_arc_endpoint_2      Circular_arc_endpoint_2;

    typedef typename CurvedKernel::Circular_arc_endpoint_2      Point;
    typedef typename CurvedKernel::Circular_arc_endpoint_2      Point_2;
  
    typedef CGAL::Tag_false                        Has_left_category;
    typedef CGAL::Tag_false 			 Has_merge_category;
  
    typedef boost::variant< Circular_arc_2, Line_arc_2 > Curve_2;
    typedef boost::variant< Circular_arc_2, Line_arc_2 > X_monotone_curve_2;

  private:
    CurvedKernel ck;
  public:
  
    Variant_traits(const CurvedKernel &k = CurvedKernel())
      : ck(k) {}

    typedef typename CurvedKernel::Compare_x_2           Compare_x_2;
    typedef typename CurvedKernel::Compare_xy_2          Compare_xy_2;
    typedef VariantFunctors::Compare_y_at_x_2<CurvedKernel>      Compare_y_at_x_2;
    typedef VariantFunctors::Compare_y_to_right_2<CurvedKernel>  Compare_y_at_x_right_2; 
    typedef VariantFunctors::Equal_2<CurvedKernel>               Equal_2;
    typedef VariantFunctors::Make_x_monotone_2<CurvedKernel>     Make_x_monotone_2;
    typedef VariantFunctors::Split_2<CurvedKernel>               Split_2;
    typedef VariantFunctors::Construct_intersections_2<CurvedKernel> Intersect_2;

    class Construct_min_vertex_2
    {
    public:
      const Point_2& operator() (const X_monotone_curve_2 & cv) const
      {
	if ( const Circular_arc_2* arc = boost::get<Circular_arc_2>( &cv ) ){
	  CGAL_kernel_precondition( CGAL::CircularFunctors::compare_xy<CurvedKernel>(arc->left(),arc->right())==CGAL::SMALLER);
	  return (arc->left());
	}  
	else {
	  const Line_arc_2* line = boost::get<Line_arc_2>( &cv );
	  CGAL_kernel_precondition( CGAL::CircularFunctors::compare_xy<CurvedKernel>(line->left(),line->right())==CGAL::SMALLER);
	  return (line->left());
	}
      }
    };

    class Construct_max_vertex_2
    {
    public:
      /*!
       * Get the right endpoint of the x-monotone curve (segment).
       * \param cv The curve.
       * \return The right endpoint.
       */
      const Point_2& operator() (const X_monotone_curve_2 & cv) const
      {
	if ( const Circular_arc_2* arc = boost::get<Circular_arc_2>( &cv ) ){
	  CGAL_kernel_precondition( CGAL::CircularFunctors::compare_xy<CurvedKernel>
				    (arc->left(),arc->right())==CGAL::SMALLER);
	  return (arc->right());
	}
	else {
	  const Line_arc_2* line = boost::get<Line_arc_2>( &cv );
	  CGAL_kernel_precondition( CGAL::CircularFunctors::compare_xy<CurvedKernel>
				    (line->left(),line->right())==CGAL::SMALLER);
	  return (line->right());
	}
      }
    };

    class Is_vertical_2
    {
    public:
      typedef bool result_type;

      bool operator() (const X_monotone_curve_2& cv) const
      {
	if ( const Line_arc_2* line = boost::get<Line_arc_2>( &cv ) ){
	  return line->supporting_line().is_vertical();
	}
	else {
	  return false;
	}
      }
    };

    

  
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
