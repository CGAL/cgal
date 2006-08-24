#ifndef AH_FILTERED_PREDICATES_H_
#define AH_FILTERED_PREDICATES_H_

#include <CGAL/enum.h>
#include <CGAL/Object.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Filtered_interval_circular_kernel/ah_primitives.h>
#include <CGAL/MP_Float.h>

CGAL_BEGIN_NAMESPACE

// I think this name is better then Interval_functors
namespace AH_functors { 

template <class AHK>
class Compare_x_2  {
  typedef typename AHK::Circular_kernel       CK;
  typedef typename AHK::Circular_arc_point_2  Circular_arc_point_2;

public:
  typedef Comparison_result result_type;

public:
  result_type operator()(const Circular_arc_point_2 &a, 
    const Circular_arc_point_2 &b) const {
    Bbox_2 bb1=a.bbox(),bb2=b.bbox();
    if( bb1.xmin()>bb2.xmax() ) return LARGER;
    if( bb1.xmax()<bb2.xmin() ) return SMALLER;
    return CK().compare_x_2_object()(a,b);
  }
};

template <class AHK>
class Compare_y_2  {
  typedef typename AHK::Circular_kernel         CK;
  typedef typename AHK::Circular_arc_point_2    Circular_arc_point_2;

public:
  typedef Comparison_result result_type;

public:
  result_type operator()(const Circular_arc_point_2 &a, 
    const Circular_arc_point_2 &b) const {	
    Bbox_2 bb1=a.bbox(),bb2=b.bbox();
    if( bb1.ymin()>bb2.ymax() )  return LARGER;
    if( bb1.ymax()<bb2.ymin() )  return SMALLER;
    return CK().compare_y_2_object()(a,b);
  }
};

template <class AHK>
class Compare_xy_2 {
  typedef typename AHK::Circular_kernel        CK;
  typedef typename AHK::Circular_arc_point_2   Circular_arc_point_2;

public:
  typedef Comparison_result result_type;

public:
  result_type operator()(const Circular_arc_point_2 &a, 
    const Circular_arc_point_2 &b) const {
    typename AHK::Compare_x_2 compx;
    typename AHK::Compare_y_2 compy;
    Comparison_result tmp;
    if( (tmp=compx(a,b))!=EQUAL) return tmp;
    return compy(a,b);
  }
};

// RETR
template <class AHK>
class In_x_range_2 {
  typedef typename AHK::Circular_kernel       CK;
  typedef typename AHK::Circular_arc_point_2  Circular_arc_point_2;
  typedef typename AHK::Circular_arc_2        Circular_arc_2;
  typedef typename AHK::Line_arc_2            Line_arc_2;

public:
  typedef bool result_type;

private:
  template <class Arc_2>
  result_type _in_x_range_2(const Arc_2 &a, 
    const Circular_arc_point_2 &p) const {
    Bbox_2 bb11 = a.source().bbox(),
           bb12 = a.target().bbox(),
           bb2=p.bbox();      
    if(bb11.xmin() > bb12.xmax()) {
      if(bb2.xmax() < bb12.xmin()) return false;
      else if(bb2.xmin() > bb11.xmax()) return false;
      else if(bb12.xmax() < bb2.xmin() &&
       	      bb2.xmax() < bb11.xmin()) return true;
    } else if(bb11.xmax() < bb12.xmin()) {
      if(bb2.xmax() < bb11.xmin()) return false;
      else if(bb2.xmin() > bb12.xmax()) return false;
      else if(bb11.xmax() < bb2.xmin() &&
              bb2.xmax() < bb12.xmin()) return true;
    } else {
      if(bb2.xmin() > std::max(bb11.xmax(),bb12.xmax())) return false;
      if(bb2.xmax() < std::min(bb11.xmin(),bb12.xmin())) return false;
    }
    
    typename CK::In_x_range_2 Range;
    return Range(a,p);
  }

public:
  result_type operator()(const Circular_arc_2 &a, 
    const Circular_arc_point_2 &p) const { 
    CGAL_precondition(a.is_x_monotone());
    return _in_x_range_2(a,p);
  }
   
  result_type operator()(const Line_arc_2 &a, 
    const Circular_arc_point_2 &p) const { 
    return _in_x_range_2(a,p);
  }
};

template <class AHK>
class Construct_circular_source_vertex_2  {
  typedef typename AHK::Circular_kernel        CK;
  typedef typename AHK::Circular_arc_point_2   Circular_arc_point_2;
  typedef typename AHK::Circular_arc_2         Circular_arc_2;

public:
  typedef Circular_arc_point_2    result_type;
  typedef const result_type &     qualified_result_type;

  template <typename T>
  result_type operator()(const T& a) const {
    return CK().construct_circular_source_vertex_2_object()(a);
  }
};

template <class AHK>
class Construct_circular_target_vertex_2 {
  typedef typename AHK::Circular_kernel        CK;
  typedef typename AHK::Circular_arc_point_2   Circular_arc_point_2;
  typedef typename AHK::Circular_arc_2         Circular_arc_2;

public:
  typedef Circular_arc_point_2 result_type;
  typedef const result_type &     qualified_result_type;

  template <typename T>
  result_type operator()(const T& a) const {
  	return CK().construct_circular_target_vertex_2_object()(a);
  }
};

template <class AHK>
class Construct_circular_min_vertex_2 {
  typedef typename AHK::Circular_kernel         CK;
  typedef typename AHK::Circular_arc_point_2    Circular_arc_point_2;
  typedef typename AHK::Circular_arc_2          Circular_arc_2;

public:
  typedef Circular_arc_point_2 result_type;
    
  template <typename T>
  result_type operator()(const T& a) const {
    return CK().construct_circular_min_vertex_2_object()(a);
  }
};

template <class AHK>
class Construct_circular_max_vertex_2 {
  typedef typename AHK::Circular_kernel          CK;
  typedef typename AHK::Circular_arc_point_2     Circular_arc_point_2;
  typedef typename AHK::Circular_arc_2           Circular_arc_2;

public:

  typedef Circular_arc_point_2 result_type;

  template <typename T>
  result_type operator()(const T& a) const { 
  	return CK().construct_circular_max_vertex_2_object()(a); 
  }
};

template <class AHK>
class Is_vertical_2 {
  typedef typename AHK::Circular_kernel            CK;
  typedef typename AHK::Circular_arc_point_2       Circular_arc_point_2;
  typedef typename AHK::Circular_arc_2             Circular_arc_2;

public:
  typedef bool result_type;
    
  template <typename T>
  result_type operator()(const T& a) const { 
  	return CK().is_vertical_2_object()(a); 
  }
};

// RETR
template <class AHK>
class Compare_y_at_x_2 {
  typedef typename AHK::Circular_kernel          CK;
  typedef typename AHK::Circular_arc_2           Circular_arc_2;
  typedef typename AHK::Circular_arc_point_2     Circular_arc_point_2;
  typedef typename AHK::Line_arc_2               Line_arc_2;

public:
  typedef Comparison_result result_type;

private:
  template <class Arc_2>
  result_type _compare_y_at_x_2(const Circular_arc_point_2 &p,
    const Arc_2 &a) const {
    CGAL_precondition_code(bool tmp=In_x_range_2<AHK>()(a,p));
    CGAL_precondition(tmp);
    Bbox_2 bb1=a.bbox(),bb2=p.bbox();
    if(bb1.ymin()>bb2.ymax()) return SMALLER;
    else if(bb1.ymax()<bb2.ymin()) return LARGER;
    return CK().compare_y_at_x_2_object()(p,a);
  }

public:
  result_type operator()(const Circular_arc_point_2 &p,
    const Circular_arc_2 &a ) const {   
    CGAL_precondition( a.is_x_monotone());
    return _compare_y_at_x_2(p,a);
  }

  result_type operator()(const Circular_arc_point_2 &p,
    const Line_arc_2 &a ) const {
    return _compare_y_at_x_2(p,a);
  }
};

template <class AHK>
class Has_on_2 {
  typedef typename AHK::Circular_kernel        CK;
  typedef typename AHK::Circular_arc_2         Circular_arc_2;
  typedef typename AHK::Circular_arc_point_2   Circular_arc_point_2;
  typedef typename AHK::Line_arc_2             Line_arc_2;

public:
  typedef bool result_type;

private:
  template <class Arc_2>
  result_type _has_on_2(const Arc_2 &a, 
    const Circular_arc_point_2 &p) const {
    Bbox_2 bb1=a.bbox(),bb2=p.bbox();
    if(do_overlap(bb1,bb2))
      return CK().has_on_2_object()(a,p);
    return false;
  }

public:
  result_type
  operator()(const Circular_arc_2 &a,
    const Circular_arc_point_2 &p ) const {     
    CGAL_precondition( a.is_x_monotone());
    return _has_on_2(a,p);
  }

  result_type operator()(const Line_arc_2 &a, 
    const Circular_arc_point_2 &p ) const {
    return _has_on_2(a,p);
  }
};

template <class AHK>
class Equal_2 {
  typedef typename AHK::Circular_kernel           CK;
  typedef typename AHK::Circular_arc_2            Circular_arc_2;
  typedef typename AHK::Circular_arc_point_2      Circular_arc_point_2;
  typedef typename AHK::Line_arc_2                Line_arc_2;

public:
  typedef bool result_type;

private:
  template <class Arc_2>
  result_type _equal_2(const Arc_2 &a,const Arc_2 &b) const {
    Bbox_2 bb11=a.source().bbox(),
           bb12=a.target().bbox(),
           bb21=b.source().bbox(),
           bb22=b.target().bbox();
    if(bb11.xmin() > bb21.xmax()) return false;
    if(bb11.xmax() < bb21.xmin()) return false;
    if(bb11.ymin() > bb21.ymax()) return false;
    if(bb11.ymax() < bb21.ymin()) return false;
    if(bb12.xmin() > bb22.xmax()) return false;
    if(bb12.xmax() < bb22.xmin()) return false;
    if(bb12.ymin() > bb22.ymax()) return false;
    if(bb12.ymax() < bb22.ymin()) return false;
    return CK().equal_2_object()(a,b);
  }

public:
  result_type operator()(const Circular_arc_point_2 &a,
    const Circular_arc_point_2 &b ) const { 
    Bbox_2 bb1=a.bbox(),bb2=b.bbox();
    if(bb1.xmin() > bb2.xmax()) return false;
    if(bb1.xmax() < bb2.xmin()) return false;
    if(bb1.ymin() > bb2.ymax()) return false;
    if(bb1.ymax() < bb2.ymin()) return false;
    return CK().equal_2_object()(a,b);
  }

  result_type operator()(const Circular_arc_2 &a, 
    const Circular_arc_2 &b ) const {
    CGAL_precondition( a.is_x_monotone());
    CGAL_precondition( b.is_x_monotone());
    return _equal_2(a,b);      
  }

  result_type operator()(const Line_arc_2 &a,
    const Line_arc_2 &b ) const {  
    return _equal_2(a,b);
  }

  result_type operator()(const Circular_arc_2 &a,
    const Line_arc_2 &b) const {
    return false;
  }

  result_type operator()(const Line_arc_2 &a, 
    const Circular_arc_2 &b) const { 
  	return false;
  }
};

template <class AHK>
class Do_overlap_2 {
  typedef typename AHK::Circular_kernel     CK;
  typedef typename AHK::Circular_arc_2      Circular_arc_2;
  typedef typename AHK::Line_arc_2          Line_arc_2;

public:
  typedef bool result_type;

private:
  template <class Arc_2>
  result_type _do_overlap_2(const Arc_2 &a, const Arc_2 &b) const {
    Bbox_2 bb1=a.bbox(),bb2=b.bbox();  
    if(do_overlap(bb1,bb2))
      return CK().do_overlap_2_object()(a,b);
    return false;        
  }

public:
  result_type operator()(const Circular_arc_2 &a, 
    const Circular_arc_2 &b ) const {
    CGAL_precondition( a.is_x_monotone());
    CGAL_precondition( b.is_x_monotone());
    return _do_overlap_2(a,b); 
  }

  result_type operator()( const Line_arc_2 &a ,
    const Line_arc_2 &b ) const {
    return _do_overlap_2(a,b);
  }

  result_type operator()(const Circular_arc_2 &a,
    const Line_arc_2 &b ) const {
    return false;
  }

  result_type operator()(const Line_arc_2 &a,
    const Circular_arc_2 &b) const {
    return false;
  }
};


// This predicate cannot be filtered
template < class AHK >
class Compare_y_to_right_2 {
  typedef typename AHK::Circular_kernel       CK;
  typedef typename AHK::Circular_arc_2        Circular_arc_2;
  typedef typename AHK::Circular_arc_point_2  Circular_arc_point_2;

public:
  typedef Comparison_result result_type;
    
  template <typename T1, typename T2>
  result_type operator()(const T1 &a1, const T2 &a2,
    const Circular_arc_point_2 &p) const { 
    return CK().compare_y_to_right_2_object()(a1, a2, p); 
  }
};

template < class AHK >
class Make_x_monotone_2 {
  typedef typename AHK::Circular_kernel       CK;
  typedef typename AHK::Rcirc_arc_2           Rcirc_arc_2;
  typedef typename AHK::Circular_arc_2        Circular_arc_2;
  typedef typename AHK::Line_arc_2            Line_arc_2;

public:
  template < class OutputIterator >
  OutputIterator operator()(const Circular_arc_2 &A, 
    OutputIterator res) { 
    std::vector<Object> vec;
	CK().make_x_monotone_2_object()(A, std::back_inserter(vec));
  	for(unsigned i=0; i<vec.size() ; i++) {
	  const Rcirc_arc_2 *tmp_arc;
	  tmp_arc=object_cast<Rcirc_arc_2>(&vec.at(i));
	  assert(tmp_arc!=NULL);
	  *res++ = make_object( Circular_arc_2(*tmp_arc) );
	} return res;
  }
    
  template < class OutputIterator >
  OutputIterator operator()(const Line_arc_2 &A, OutputIterator res) { 
	*res++ = make_object(A);
	return res;
  }
};


template < class AHK >
class Intersect_2 {
public:

  typedef typename AHK::Circular_kernel          CK;
  typedef typename AHK::Circular_arc_2           Circular_arc_2;
  typedef typename AHK::Circular_arc_point_2     Circular_arc_point_2;
  typedef typename AHK::Line_arc_2               Line_arc_2;
  typedef typename AHK::Rcirc_arc_2              Rcirc_arc_2;
  typedef typename AHK::Rline_arc_2              Rline_arc_2;
  typedef typename AHK::Rcirc_arc_point_2        Rcirc_arc_point_2;
  typedef typename AHK::Circle_2                 Circle;

  template < class OutputIterator >
  OutputIterator operator()(const Circle & c1, 
    const Circle & c2, OutputIterator res) { 
    // Filtering
    Bbox_2 bb1=c1.bbox(),bb2=c2.bbox();
    if(!do_overlap(bb1,bb2)) return res;
    if(CGAL::CGALi::intersect<CK>(c1,c2) == CGAL::CGALi::DONT_INTERSECT) {
      return res;
    }
    std::vector<Object> vec;
    CK().intersect_2_object()(c1,c2,std::back_inserter(vec));
    for(unsigned i=0; i<vec.size() ; i++) {
      const std::pair<Rcirc_arc_point_2, unsigned> *tmp_point;
      if ((tmp_point=object_cast<std::pair<Rcirc_arc_point_2, unsigned> >(
            &vec.at(i)))!=NULL ) {
        *res++ = make_object(std::make_pair(
          Circular_arc_point_2(tmp_point->first),tmp_point->second));
      } else *res++=vec.at(i);
    }  return res;
  }

  template < class OutputIterator >
  OutputIterator operator()(const Circular_arc_2 & c1, 
    const Circular_arc_2 & c2, OutputIterator res) { 
    
    // Filtering
    Bbox_2 bb1=c1.bbox(),bb2=c2.bbox();
    if(!do_overlap(bb1,bb2)) return res;
    if(CGAL::CGALi::intersect<CK>(c1,c2) == CGAL::CGALi::DONT_INTERSECT) {
      return res;
    } 
    
    std::vector<Object> vec;
	CK().intersect_2_object()(c1,c2,std::back_inserter(vec)); 
   	for(unsigned i=0; i<vec.size() ; i++) {
	  const Rcirc_arc_2 *tmp_arc;
  	  if ( (tmp_arc=object_cast<Rcirc_arc_2>(&vec.at(i)) )!=NULL )
	    *res++ = make_object( Circular_arc_2(*tmp_arc) );
	  else {
        const std::pair<Rcirc_arc_point_2, unsigned> *tmp_point;
        tmp_point=object_cast<std::pair<Rcirc_arc_point_2, unsigned> >(
          &vec.at(i));
        assert(tmp_point!=NULL);
        *res++ = make_object( std::make_pair(
          Circular_arc_point_2(tmp_point->first),tmp_point->second));   
	  }
    } return res;
  }

  template < class OutputIterator >
  OutputIterator operator()(const Line_arc_2 & c1, 
    const Line_arc_2 & c2, OutputIterator res) { 
    
    // Filtering
    Bbox_2 bb1=c1.bbox(),bb2=c2.bbox();
    if(!do_overlap(bb1,bb2)) return res;
    if(CGAL::CGALi::intersect<CK>(c1,c2) == CGAL::CGALi::DONT_INTERSECT) {
      return res;
    }
    
    std::vector<Object> vec;
	CK().intersect_2_object()(c1,c2,std::back_inserter(vec)); 
    for(unsigned i=0; i<vec.size() ; i++) {
	  const Rline_arc_2 *tmp_arc;
  	  if ((tmp_arc=object_cast<Rline_arc_2>(&vec.at(i)) )!=NULL)
	    *res++ = make_object( Line_arc_2(*tmp_arc) );
	  else {
        const std::pair<Rcirc_arc_point_2, unsigned> *tmp_point;
        tmp_point=object_cast<std::pair<Rcirc_arc_point_2, unsigned> >(
          &vec.at(i));
        assert(tmp_point!=NULL);
	    *res++ = make_object( std::make_pair(Circular_arc_point_2(
	      tmp_point->first),tmp_point->second));
      }
	} return res;	
  }

  template < class OutputIterator >
  OutputIterator operator()(const Circular_arc_2 & c1, 
    const Line_arc_2 & c2, OutputIterator res) { 
    // Filtering
    Bbox_2 bb1=c1.bbox(),bb2=c2.bbox();
    if(!do_overlap(bb1,bb2)) return res;
    if(CGAL::CGALi::intersect<CK>(c2,c1) == CGAL::CGALi::DONT_INTERSECT) {
      return res;
    }
    
	std::vector<Object> vec;
	CK().intersect_2_object()(c1,c2,std::back_inserter(vec)); 
 	for(unsigned i=0; i<vec.size() ; i++) {
	  const Rcirc_arc_2 *tmp_arc;
	  const Rline_arc_2 *tmp_line;
	  // Can this happen?
	  if ( (tmp_arc=object_cast<Rcirc_arc_2>(&vec.at(i)) )!=NULL )  
	    *res++ = make_object( Circular_arc_2(*tmp_arc) );
      //Can this also happen?
	  else if ( (tmp_line=object_cast<Rline_arc_2>(&vec.at(i)) )!=NULL ) 
	    *res++ = make_object( Line_arc_2(*tmp_line) );
	  else {
        const std::pair<Rcirc_arc_point_2, unsigned> *tmp_point;
        tmp_point=object_cast<std::pair<Rcirc_arc_point_2, unsigned> >(
          &vec.at(i));
        assert(tmp_point!=NULL);
        *res++ = make_object( std::make_pair(
          Circular_arc_point_2(tmp_point->first),tmp_point->second));
      }
	} return res;	
  }

  template < class OutputIterator >
  OutputIterator operator()(const Line_arc_2 & c1, 
    const Circular_arc_2 & c2, OutputIterator res) {
    return operator()(c2,c1,res);
  }
};

template < class AHK >
class Split_2 {

  typedef typename AHK::Circular_kernel            CK;
  typedef typename AHK::Rcirc_arc_2              Rcirc_arc_2;
  typedef typename AHK::Rline_arc_2              Rline_arc_2;
  typedef typename AHK::Circular_arc_2           Circular_arc_2;
  typedef typename AHK::Line_arc_2               Line_arc_2;
  typedef typename AHK::Circular_arc_point_2  Circular_arc_point_2;

public:

  typedef void result_type;

  result_type operator()(const Circular_arc_2 &A, 
    const Circular_arc_point_2 &p,
    Circular_arc_2 &ha1, Circular_arc_2 &ha2) const {  
    CK().split_2_object()(A, p, ha1, ha2);
  }
    
  result_type operator()(const Line_arc_2 &A,
    const Circular_arc_point_2 &p,
    Line_arc_2 &ha1, Line_arc_2 &ha2) const {  
      CK().split_2_object()(A, p, ha1, ha2);
  }
};

} //AH_functors

CGAL_END_NAMESPACE   

#endif /*AH_FILTERED_PREDICATES_H_*/
