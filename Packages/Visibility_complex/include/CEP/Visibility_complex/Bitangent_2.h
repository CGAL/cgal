#ifndef CGAL_BITANGENT_2_H
#define CGAL_BITANGENT_2_H

#include <cmath>
#include <list>

#ifndef CGAL_CONVEX_ARC_2_H
#include <CEP/Visibility_complex/Arc_2.h>
#endif 

CGAL_BEGIN_NAMESPACE

// -----------------------------------------------------------------------------
// -------------------- Base Bitangent class -----------------------------------

template< class D_ >
class Bitangent_base 
{
public:
    // -------------------------------------------------------------------------
    typedef D_                                        Disk;
    typedef const Disk*                           Disk_handle;
    // -------------------------------------------------------------------------
    enum Type { LL , RR , LR , RL };
    class Type_util {
	public:
	Type operator()(bool b1,bool b2) const {
	    if ( b1 &&  b2) return LL;
	    if (!b1 &&  b2) return RL;
	    if (!b1 && !b2) return RR;
	    return LR;
	}
    };
public:
    // Constructeurs -----------------------------------------------------------
    Bitangent_base() 
	: _type(LL) , _source_object(0)     , _target_object(0) { }
    Bitangent_base(Type t, Disk_handle o1, Disk_handle o2)  
	: _type(t)  , _source_object(o1) , _target_object(o2)   { }
    // ----------------------- Operators ---------------------------------------
    bool operator==(const Bitangent_base& b) const{
	return (b.type() == type() && source_object() == b.source_object() && 
				      target_object() == b.target_object());
    }
    bool operator!=(const Bitangent_base& b) const{ return !(*this == b); }
    // -------- return the opposite oriented bitangent -------------------------
    Type type()                    const { return _type;          }
    // ---- accesing the two objects defining the bitangent --------------------
    Disk_handle source_object() const { return _source_object; }
    Disk_handle target_object() const { return _target_object; }
    // ---------- informations on the type of the bitangent -------------------- 
    bool is_left_right()   const { return (type() == LR);                 }
    bool is_left_left()    const { return (type() == LL);                 }
    bool is_right_right()  const { return (type() == RR);                 }
    bool is_right_left()   const { return (type() == RL);                 }
    bool is_left_xx()      const { return (type() == LL || type() == LR); }
    bool is_xx_left()      const { return (type() == LL || type() == RL); }
    bool is_right_xx()     const { return (type() == RR || type() == RL); }
    bool is_xx_right()     const { return (type() == RR || type() == LR); }
    bool is_internal()     const { return (type() == RL || type() == LR); }
    bool is_external()     const { return (type() == LL || type() == RR); }
    // -------------------------------------------------------------------------
private :
    Type              _type;
    Disk_handle   _source_object , _target_object;
};

//------------------------------------------------------------------------------
//------------------ General definition of Bitangent_2 -------------------------
//------------------------------------------------------------------------------

template< class D_ >
struct Bitangent_2 
    : public Bitangent_base<D_> 
{
    //--------------------------------------------------------------------------
    typedef D_                          Disk;
    typedef Arc_2<Disk>      Arc_2;
    typedef Bitangent_base<Disk>    Base;
    typedef typename Base::Disk_handle        Disk_handle;
    typedef typename Base::Type                  Type;
    typedef typename D_::Point_2        Point_2;
    typedef typename D_::Segment_2      Segment_2;
    //--------------------------------------------------------------------------
    Bitangent_2() : Base() { }
    Bitangent_2(const Point_2& v1 , const Point_2& v2 , Type t ,
		Disk_handle start, Disk_handle finish)
	: Base(v1,v2,t,start,finish) { }
    Bitangent_2(Type t, const Arc_2& source,const Arc_2& target) 
	: Base(t,source.object(),target.object()) { }
    Bitangent_2(Type t ,  Disk_handle o1 , Disk_handle o2) 
	: Base(t,o1,o2) { }
    //--------------------------------------------------------------------------
    bool operator==(const Bitangent_2& b) const 
    { return Base::operator==(b); }
    bool operator!=(const Bitangent_2& b) const 
    { return Base::operator!=(b); }
    //--------------------------------------------------------------------------
};

//------------------------------------------------------------------------------
//------------------ Partial Specialization for CGAL::Polygon_2 ----------------
//------------------------------------------------------------------------------

#ifdef CGAL_POLYGON_2_H
template< class R_ , class C_ >
struct Bitangent_2 < Polygon_2<R_,C_> >
    : public R_::Segment_2 , 
      public Bitangent_base< Polygon_2<R_,C_> >
{
    // -------------------------------------------------------------------------
    typedef R_                                      R;
    typedef typename R::FT                          FT;
    typedef Polygon_2<R_,C_>                        Disk;
    typedef Bitangent_base<Disk>                Base;
    typedef CGAL::Arc_2<Disk>                  Arc_2;
    typedef CGAL::Arc_2<Disk>                  CCC;
    typedef typename Base::Disk_handle                    Disk_handle;
    typedef typename R_::Segment_2                  Segment_2;
    typedef typename R_::Point_2                    Point_2;
    typedef typename Base::Type                              Type;
    // -------------------------------------------------------------------------
    typedef typename CCC::Vertex_const_iterator Vertex_iterator;
    typedef typename CCC::Vertex_const_iterator Vertex_const_iterator;
    typedef typename Disk::Vertex_const_circulator Vertex_circulator;
    typedef typename Disk::Vertex_const_circulator Vertex_const_circulator;
    // Constructeurs -----------------------------------------------------------
    Bitangent_2() : Base() { }
    Bitangent_2(const Point_2& v1 , const Point_2& v2 , Type t ,
		Disk_handle start, Disk_handle finish)
	: Segment_2(v1,v2) , Base(t,start,finish) { }
    Bitangent_2(Type t, const Arc_2& source, const Arc_2& target) {
	if (source.begin() == source.end() || target.begin() == target.end()) {
	    *this = Bitangent_2(t,source.object(),target.object());
	    return;
	}
	bool t1,t2;

	Vertex_const_iterator it_source = source.begin();
	Vertex_const_iterator it_source_succ;
	Vertex_const_iterator it_target = target.begin();
	Vertex_const_iterator it_target_succ;

	Point_2 p_source(*it_source);
	Point_2 p_target(*it_target);
	Point_2 p_target_succ,p_source_succ;

	bool finished_source = false;
	bool finished_target = false;

	do {
	    if (it_source == source.end()) finished_source = true;
	    if (it_target == target.end()) finished_target = true;
	    if (finished_target && finished_source) break;

	    p_source = *(it_source); 
	    p_target = *(it_target); 
	    
	    it_source_succ = it_source; ++it_source_succ; 
	    if (it_source_succ != source.end()) 
		 p_source_succ = *(it_source_succ);
	    else p_source_succ = p_source;

	    it_target_succ = it_target; ++it_target_succ; 
	    if (it_target_succ != target.end()) 
		 p_target_succ = *(it_target_succ);
	    else p_target_succ = p_target;
	    
	    t1 = ( (t == RR || t == RL) && !leftturn (p_source,p_target,p_source_succ))
	      || ( (t == LL || t == LR) && !rightturn(p_source,p_target,p_source_succ));
	    if ( collinear(p_source,p_target,p_source_succ) && p_source != p_source_succ)
	      t1 = are_ordered_along_line(p_source_succ,p_source,p_target);
	    
	    t2 = ( (t == LL || t == RL) && !rightturn(p_source,p_target,p_target_succ))
	      || ( (t == RR || t == LR) && !leftturn (p_source,p_target,p_target_succ));
	    if ( collinear(p_source,p_target,p_target_succ) && p_target != p_target_succ)
	      t2 = are_ordered_along_line(p_source,p_target,p_target_succ);

	    if (finished_target)         ++it_source;
	    else if (finished_source)    ++it_target;
	    else if (t1 == 1 && t2 == 0) ++it_target;      
	    else if (t1 == 0 && t2 == 1) ++it_source;
	    else if (t1 == 0 && t2 == 0) {
		// Le test d'orientation est different suivant que 
		// source_is_left et target_is_left sont de meme 
		// signe ou non 
		Point_2 tmp(p_target_succ + (p_source_succ - p_target));
		if ( (    (t == LR || t == RL)    
		       && rightturn(p_source_succ, p_source, tmp))
		     || ( (t == LL || t == RR)
		       && leftturn (p_source_succ, p_source, tmp)))
		     ++it_target;
		else ++it_source;
	    }
	} while (t1 == 0 || t2 == 0);
	*this = Bitangent_2(p_source,p_target,t,source_obj,target_obj);
    }
    Bitangent_2(Type t ,  Disk_handle o1 , Disk_handle o2){ 
	bool exists = true;
	bool found = false;
	bool found_start = false;
	bool found_finish = false;
	bool found_start_succ = false;
	bool found_start_pred = false;
	bool found_finish_succ = false;
	bool found_finish_pred = false;
	int  loop_counter = o1->size() + o2->size();

	Vertex_circulator start  = o1->vertices_circulator();
	Vertex_circulator finish = o2->vertices_circulator();
	do {
	Vertex_circulator start_succ  = start;  ++start_succ;
	Vertex_circulator finish_succ = finish; ++finish_succ; 
	Vertex_circulator start_pred  = start;  --start_pred;
	Vertex_circulator finish_pred = finish; --finish_pred;

	if (t == LL || t == LR) { 
	    // --------------------------------------------------------------------
	    found_start_succ = 
	    (leftturn(*start,*finish,*start_succ)                      || 
	     *start == *start_succ                                     ||
	     (collinear(*start,*finish,*start_succ) &&
	      are_ordered_along_line(*start_succ,*start,*finish)));
	    // --------------------------------------------------------------------
	    found_start_pred = 
	    (leftturn(*start,*finish,*start_pred)                      || 
	     *start == *start_pred                                     || 
	     (collinear(*start,*finish,*start_pred) &&
	      are_ordered_along_line(*start_pred,*start,*finish)));
	    // --------------------------------------------------------------------
	}
	else   {
	    // --------------------------------------------------------------------
	    found_start_succ = 
	    (rightturn(*start,*finish,*start_succ)                     || 
	     *start == *start_succ                                     ||
	     (collinear(*start,*finish,*start_succ) &&
	      are_ordered_along_line(*start_succ,*start,*finish)));
	    // --------------------------------------------------------------------
	    found_start_pred = 
	    (rightturn(*start,*finish,*start_pred)                     || 
	     *start == *start_pred                                     || 
	     (collinear(*start,*finish,*start_pred) &&
	      are_ordered_along_line(*start_pred,*start,*finish)));
	    // --------------------------------------------------------------------
	}
	found_start = (found_start_pred && found_start_succ);

	if (t == LL || t == RL) {
	    // --------------------------------------------------------------------
	    found_finish_succ = 
	    (rightturn(*finish,*start,*finish_succ)                    || 
	     *finish == *finish_succ                                   ||
	     (collinear(*finish,*start,*finish_succ) &&
	      are_ordered_along_line(*start,*finish,*finish_succ)));
	    // --------------------------------------------------------------------
	    found_finish_pred = 
	    (rightturn(*finish,*start,*finish_pred)                    || 
	     *finish == *finish_pred                                   ||
	     (collinear(*finish,*start,*finish_pred) && 
	      are_ordered_along_line(*start,*finish,*finish_pred)));
	    // --------------------------------------------------------------------
	}
	else  {
	    // --------------------------------------------------------------------
	    found_finish_succ = 
	    (leftturn(*finish,*start,*finish_succ)                     || 
	     *finish == *finish_succ                                   ||
	     (collinear(*finish,*start,*finish_succ) &&
	      are_ordered_along_line(*start,*finish,*finish_succ)));
	    // --------------------------------------------------------------------
	    found_finish_pred = 
	    (leftturn(*finish,*start,*finish_pred)                     || 
	     *finish == *finish_pred                                   ||
	     (collinear(*finish,*start,*finish_pred) && 
	      are_ordered_along_line(*start,*finish,*finish_pred)));
	    // --------------------------------------------------------------------
	}
	found_finish = (found_finish_pred && found_finish_succ);

	/* If found == true we have found the bitangent */
	found = (found_start && found_finish);

	/* Otherwise we compute one more determinant to decide on which */
	/* object to advance */

	if (found == false) 
	  {
	    if (( (t == LR || t == RL) 
		      && rightturn(*start_succ,*start,*finish_succ + (*start_succ - *finish)))||
		 ((t == LL || t == RR) 
		      && leftturn (*start_succ,*start,*finish_succ + (*start_succ - *finish))))
	      {
		++finish;
	      }
	    else if (((t == LR || t == RL) 
		      && leftturn (*start_succ,*start,*finish_succ + (*start_succ - *finish)))||
		     ((t == LL || t == RR) 
		      && rightturn(*start_succ,*start,*finish_succ + (*start_succ - *finish))))
	      {
		++start;
	      }
	    else if (found_start == true && found_finish == false) ++finish;
	    else ++start;
	  }
	--loop_counter; 
	if (loop_counter < -100) {exists = false; break; }
	} while (found == false);

	if (exists){ *this = Bitangent_2(*start,*finish,t,o1,o2); }
	else { *this = Bitangent_2(); }
    }
    //--------------------------------------------------------------------------
    bool operator==(const Bitangent_2& b) const 
    { return Base::operator==(b); }
    bool operator!=(const Bitangent_2& b) const 
    { return Base::operator!=(b); }
    //--------------------------------------------------------------------------
};
#endif // CGAL_POLYGON_2_H

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//------------------ Partial Specialization for CGAL::Point_2 ------------------
//------------------------------------------------------------------------------

template< class R_>
struct Bitangent_2 < Point_2<R_> >
    : public R_::Segment_2 , 
      public Bitangent_base< Point_2<R_> >
{
    // -------------------------------------------------------------------------
    typedef R_                                      R;
    typedef typename R::FT                          FT;
    typedef Point_2<R_>                             Disk;
    typedef Bitangent_base<Disk>                Base;
    typedef Arc_2<Disk>                  Arc_2;
    typedef typename R_::Segment_2                  Segment_2;

    typedef typename R_::Point_2                    Point_2;
    typedef typename Base::Type                              Type;
    typedef typename Base::Disk_handle                    Disk_handle;
    // Constructeurs -----------------------------------------------------------
    Bitangent_2() : Base() { }
    Bitangent_2(const Point_2& v1 , const Point_2& v2 , 
		Type t, Disk_handle start, Disk_handle finish)
	: Segment_2(v1,v2) , Base(t,start,finish) { }
    Bitangent_2(Type t, const Arc_2& source, const Arc_2& target) 
	: Segment_2(*source.object(),*target.object()) ,
	  Base(t,source.object(),target.object()) { }
    Bitangent_2(Type t , Disk_handle o1 , Disk_handle o2) 
	: R_::Segment_2(*o1,*o2) , Base(t,o1,o2) { }
    //--------------------------------------------------------------------------
    bool operator==(const Bitangent_2& b) const 
    { return Base::operator==(b); }
    bool operator!=(const Bitangent_2& b) const 
    { return Base::operator!=(b); }
    //--------------------------------------------------------------------------
};

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//---------- Partial Specialization for CGAL::Circle_by_radius -----------------
//------------------------------------------------------------------------------

#ifdef CGAL_CIRCLE_BY_RADIUS_2_H
template< class R_>
class Bitangent_2 < Circle_by_radius_2<R_> >
    : public R_::Segment_2 , 
      public Bitangent_base< Circle_by_radius_2<R_> >
{
public:
    // -------------------------------------------------------------------------
    typedef R_                                     R;
    typedef typename R::FT                         FT;
    typedef Circle_by_radius_2<R_>                 Disk;
    typedef Bitangent_base<Disk>               Base;
    typedef Base::Disk_handle                   Disk_handle;
    typedef Arc_2<Disk>                 Arc_2;
    typedef typename R::Segment_2                  Segment_2;
    typedef typename R::Point_2                    Point_2;
    typedef Base::Type                             Type;
    // -------------------------------------------------------------------------
public:
    // Constructeurs -----------------------------------------------------------
    Bitangent_2() : Base() { }
    Bitangent_2(const Point_2& v1 , const Point_2& v2 , 
		Type t , Disk_handle start, Disk_handle finish)
	: Segment_2(v1,v2) , Base(t,start,finish) { }
    Bitangent_2(Type t ,  Disk_handle o1 , Disk_handle o2) : Base(t,o1,o2) 
    { 
	compute();
    }
    Bitangent_2(Type t, const Arc_2& source, const Arc_2& target) 
    { 
	*this = Bitangent_2(t,source.object(),target.object()); 
	compute();
    }
    //--------------------------------------------------------------------------
    Point_2 source() const 
    {
	return Point_2(source_object()->center().x() + R1*pbra,
		       source_object()->center().y() - R1*parb);
    }
    Point_2 target() const
    {
	return Point_2(target_object()->center().x() + R2*pbra,
		       target_object()->center().y() - R2*parb);
    }
    //--------------------------------------------------------------------------
    bool operator==(const Bitangent_2& b) const 
    { return Base::operator==(b); }
    bool operator!=(const Bitangent_2& b) const 
    { return Base::operator!=(b); }
    //--------------------------------------------------------------------------
private:
    void compute() {
	R1 = (is_left_xx())?  source_object()->radius(): 
			    - source_object()->radius();
	R2 = (is_xx_left())?  target_object()->radius(): 
			    - target_object()->radius();

	a  = target_object()->center().x() - source_object()->center().x();
	b  = target_object()->center().y() - source_object()->center().y();
	r  = R2 - R1;
	FT aabb = a*a + b*b;
	p  = CGAL::sqrt(aabb - r*r);
	pbra = (p*b - r*a) / aabb;
	parb = (p*a + r*b) / aabb;
    }
    FT R1,R2,a,b,r,p,pbra,parb;
};
#endif  // CGAL_CIRCLE_BY_RADIUS_2_H
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//------------------ Partial Specialization for CGAL::Segment_2 ----------------
//------------------------------------------------------------------------------

#ifdef CGAL_SEGMENT_2_H
template< class R_>
class Bitangent_2 < Segment_2<R_> >
    : public R_::Segment_2 , 
      public Bitangent_base< Segment_2<R_> >
{
public:
    // -------------------------------------------------------------------------
    typedef R_                          R;
    typedef typename R::FT              FT;
    typedef typename R::Segment_2       Segment_2;
    typedef typename R::Point_2         Point_2;
    typedef Bitangent_base<Segment_2>   Base;
    typedef typename Base::Disk                  Disk;
    typedef typename Base::Disk_handle           Disk_handle;
    typedef Arc_2<Disk>                 Arc_2;
    typedef typename Base::Type                  Type;
    // -------------------------------------------------------------------------
public:
    // Constructeurs -----------------------------------------------------------
    Bitangent_2() : Base() { }
    Bitangent_2(const Point_2& v1 , const Point_2& v2 , 
		Type t, Disk_handle start, Disk_handle finish)
	: Segment_2(v1,v2) , Base(t,start,finish) { }
    Bitangent_2(Type t ,  Disk_handle o1 , Disk_handle o2) {
	if (is_bitangent(t,o1->source(),o2->source(),
			   o1->target(),o2->target()))
	     *this = Bitangent_2(o1->source(),o2->source(),t,o1,o2);
	else if (is_bitangent(t,o1->source(),o2->target(),
				o1->target(),o2->source()))
	     *this = Bitangent_2(o1->source(),o2->target(),t,o1,o2);
	else if (is_bitangent(t,o1->target(),o2->source(),
				o1->source(),o2->target()))
	     *this = Bitangent_2(o1->target(),o2->source(),t,o1,o2);
	else *this = Bitangent_2(o1->target(),o2->target(),t,o1,o2);
    }
    Bitangent_2(Type t, const Arc_2& source, const Arc_2& target)
    { *this = Bitangent_2(t,source.object(),target.object()); }
    //--------------------------------------------------------------------------
    bool operator==(const Bitangent_2& b) const 
    { return Base::operator==(b); }
    bool operator!=(const Bitangent_2& b) const 
    { return Base::operator!=(b); }
    // -------------------------------------------------------------------------
private:
    // b = (p1 , p2) and q1 and q2 are the two other points respectively on
    // source and target object. Returns true if the bitangent is valid.
    bool is_bitangent(Type t , const Point_2& p1, const Point_2& p2,
			       const Point_2& q1, const Point_2& q2)
    {
	return ((t == LL 
		 && (leftturn (p1,p2,q1) ||
		     (collinear(p1,p2,q1) && are_ordered_along_line(q1,p1,p2)))
		 && (leftturn (p1,p2,q2) ||
		     (collinear(p1,p2,q2) && are_ordered_along_line(p1,p2,q2))))
	     ||	(t == LR 
		 && (leftturn (p1,p2,q1) ||
		     (collinear(p1,p2,q1) && are_ordered_along_line(q1,p1,p2)))
		 && (rightturn(p1,p2,q2) ||
		     (collinear(p1,p2,q2) && are_ordered_along_line(p1,p2,q2))))
             ||	(t == RR 
		 && (rightturn(p1,p2,q1) ||
		     (collinear(p1,p2,q1) && are_ordered_along_line(q1,p1,p2)))
		 && (rightturn(p1,p2,q2) ||
		     (collinear(p1,p2,q2) && are_ordered_along_line(p1,p2,q2))))
             ||	(t == RL 
		 && (rightturn(p1,p2,q1) ||
		     (collinear(p1,p2,q1) && are_ordered_along_line(q1,p1,p2)))
		 && (leftturn (p1,p2,q2) ||
		     (collinear(p1,p2,q2) && are_ordered_along_line(p1,p2,q2)))));
    }
    //--------------------------------------------------------------------------
};
#endif // CGAL_SEGMENT_2_H

//------------------------------------------------------------------------------

template < class D_ >
std::ostream &
operator<<(std::ostream &os, const Bitangent_2<D_> &b)
{
    switch (b.type()) {
	case Bitangent_2<D_>::LL : os<<"LL "; break;
	case Bitangent_2<D_>::LR : os<<"LR "; break;
	case Bitangent_2<D_>::RL : os<<"RL "; break;
	case Bitangent_2<D_>::RR : os<<"RR "; break;
	default: os<<"Unknown type ";
    }
    typedef Bitangent_2<D_> Bitangent_2;
    os << static_cast<const typename Bitangent_2::Segment_2&>(b);
    return os;
}

//------------------------------------------------------------------------------

CGAL_END_NAMESPACE

#endif
