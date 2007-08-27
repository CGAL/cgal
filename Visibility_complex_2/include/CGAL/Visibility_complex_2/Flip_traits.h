// Copyright (c) 2001-2004  ENS of Paris (France).
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
// Author(s)     : Pierre Angelier, Michel Pocchiola

#ifndef CGAL_VISIBILITY_COMPLEX_2_FLIP_TRAITS_H
#define CGAL_VISIBILITY_COMPLEX_2_FLIP_TRAITS_H

CGAL_BEGIN_NAMESPACE
namespace Visibility_complex_2_details {

// -----------------------------------------------------------------------------

#ifdef FLIP_STATISTICS
static int chi1_calls = 0;
static int chi2_calls = 0;
static int chi1_naive_calls = 0;
static int chi2_naive_calls = 0;
static int chi3_calls = 0;
static int chi3_used  = 0;
static int phiR_equal_phi  = 0;
static int phiL_equal_phi  = 0;
static int psi_computed = 0;
#endif

// -----------------------------------------------------------------------------

template < class A_ >
class Chi2_flip_traits 
{
public:
    // -------------------------------------------------------------------------
    typedef A_ Antichain;
    typedef typename Antichain::Disk_handle Disk_handle;
    typedef typename Antichain::Vertex_handle Vertex_handle;
    typedef typename Antichain::Edge_handle   Edge_handle;
    typedef typename Antichain::Face_handle   Face_handle;
    // -------------------------------------------------------------------------
    // CONSTRUCTOR
    Chi2_flip_traits(/*const*/ Antichain* an) : a(an) { }
    // -------------------------------------------------------------------------
  template <class Tr> struct Skip {
    typename Tr::Sup sup;
    typename Tr::Ccw_edge ccw_edge;
    typename Tr::Inf inf;
    typename Tr::Cw_edge cw_edge;
    typename Tr::Source_object source_object;
    typename Tr::Target_object target_object;
    typename Tr::Cw_source_edge cw_source_edge;
    typename Tr::Cw_target_edge cw_target_edge;
    typename Tr::Ccw_source_edge ccw_source_edge;
    typename Tr::Ccw_target_edge ccw_target_edge;

    
    bool is_in_pt(Vertex_handle v) const {
      if (v->is_constraint()) return true;
      if (inf(v)) {
        if (inf(v)->top_edge()) return inf(v)->top_edge()->is_in_antichain();
        if (inf(v)->bottom_edge()) return inf(v)->bottom_edge()->is_in_antichain();
      }
      return false;
    }

    Edge_handle next_edge(Edge_handle e) const {
      Disk_handle o1=e->object();
      Disk_handle o2;
      if (sup(e)->target_object()==o1) o2=sup(e)->source_object(); else
        o2=sup(e)->target_object(); 
      if (is_in_pt(sup(e))) return ccw_edge(sup(e),o2); else 
        return ccw_edge(sup(e),o1);
    }

    Edge_handle forward(Edge_handle e) {
      while (e&&!is_in_pt(sup(e))) e=ccw_edge(sup(e),e->object());
      return e;
    }

    Edge_handle backward(Edge_handle e) {
      while (e&&!is_in_pt(inf(e))) e=cw_edge(inf(e),e->object());
      return e;
    }

    Vertex_handle cwL(Vertex_handle v) {
      Disk_handle d=target_object(v);
      do v=inf(cw_edge(v,d)); while (v&&!is_in_pt(v));
      return v;
    }

    Vertex_handle cwR(Vertex_handle v) {
      Disk_handle d=source_object(v);
      do v=inf(cw_edge(v,d)); while (v&&!is_in_pt(v));
      return v;
    }

    Vertex_handle ccR(Vertex_handle v) {
      Disk_handle d=target_object(v);
      do v=sup(ccw_edge(v,d)); while (v&&!is_in_pt(v));
      return v;
    }

    Vertex_handle ccL(Vertex_handle v) {
      Disk_handle d=source_object(v);
      do v=sup(ccw_edge(v,d)); while (v&&!is_in_pt(v));
      return v;
    }
  };



    template < class Tr >
    Edge_handle   first_cusp_edge(Vertex_handle v, Tr)    const;
    template < class Tr >
    Edge_handle   second_cusp_edge(Vertex_handle v, Tr)   const;
    template < class Tr >
    Edge_handle   third_cusp_edge(Vertex_handle v, Tr)    const;
    template < class Tr >
    Edge_handle   start_edge(Vertex_handle v, Tr)         const;
    template < class Tr >
    Vertex_handle walk(Edge_handle r, Edge_handle l, Tr)  const;

    // -------------------------------------------------------------------------
    template < class Tr >
    Vertex_handle operator()(const Vertex_handle& v , Tr tr) const
    {
	typename Tr::Right_traits right;
	typename Tr::Left_traits  left;
	return walk(start_edge(v,right),start_edge(v,left),tr);
    }
    // -------------------------------------------------------------------------
public:
    // -------------------------------------------------------------------------
    // Returns false if phiv is tangent to r
    template < class Tr >
    bool advance(Edge_handle r, Vertex_handle phiv, Tr) const;
    // -------------------------------------------------------------------------
    // Returns the vertex leaving the support disk of r and entering the support
    // disk of l. Check the equality to sup(r),inf(r),sup(l),inf(l) to avoid
    // duplicates
    template < class Tr >
    Vertex_handle bit(Edge_handle r, Edge_handle l, Tr) const;
  //template < class Tr >
  //Vertex_handle bit(Edge_handle r, Edge_handle l, Tr) const
  //{ return bit(r,l, Tr()); }
    // -------------------------------------------------------------------------
protected:
    /*const*/ Antichain* a; 
};

// -----------------------------------------------------------------------------

template < class A_ >
template < class Tr >
inline bool
Chi2_flip_traits<A_>::advance(Edge_handle r, Vertex_handle phiv, Tr) const
{
    typename Tr::Sup sup; typename Tr::Sign sign;
    typename Tr::Less_vertex_handle chi1; 
    Skip<Tr> skip;
    if (sup(r) == 0 || (sign(r) && !sup(r)->is_constraint()&&skip.is_in_pt(sup(r)))
        ) return false;
    //typename Tr::Target_object target_object; 
    //if (target_object(sup(r)) != r->object()) return false;
#ifdef FLIP_STATISTICS
    ++chi1_calls;
    if (chi1(phiv,sup(r))) ++chi1_naive_calls;
#endif
    return chi1(sup(r),phiv);
}

// -----------------------------------------------------------------------------

template < class A_ >
template < class Tr >
inline 
typename Chi2_flip_traits<A_>::Vertex_handle
Chi2_flip_traits<A_>::bit(Edge_handle r, Edge_handle l, Tr) const
{
    typename Tr::Sup sup; typename Tr::Inf inf;
    typename Tr::Vertex_creator vertex;
    Vertex_handle phi = vertex(r,l);
    if (sup(r) != 0 && *phi == *sup(r)) { delete phi; return sup(r); }
    if (inf(r) != 0 && *phi == *inf(r)) { delete phi; return inf(r); }
    if (sup(l) != 0 && *phi == *sup(l)) { delete phi; return sup(l); }
    if (inf(l) != 0 && *phi == *inf(l)) { delete phi; return inf(l); }
    return phi;
}

// -----------------------------------------------------------------------------
// Returns the first arc of Rleave(b). 
// Constraint case not treated FIXME

template < class A_ >
template < class Tr >
typename Chi2_flip_traits<A_>::Edge_handle
Chi2_flip_traits<A_>::first_cusp_edge(Vertex_handle v, Tr ) const 
{

    // Operators used by this method
    typename Tr::Dr dr; typename Tr::Ul ul; 
    typename Tr::Sup sup; 
    typename Tr::CcR ccR; typename Tr::CwR cwR; 
    typename Tr::Ccw_source_edge ccw_source_edge;
    typename Tr::Ccw_target_edge ccw_target_edge;
    typename Tr::Cw_target_edge  cw_target_edge;
    typename Tr::Cw_source_edge  cw_source_edge;
    typename Tr::Is_left_xx is_left_xx;
    typename Tr::Is_xx_left is_xx_left;
    Skip<Tr> skip;
    // ------------------------------------------------------------------------- 
    // xx-right bitangents
    if (!is_xx_left(v)) return ccw_target_edge(v);
    // ------------------------------------------------------------------------- 
    Vertex_handle cusp1 = sup(dr(cw_source_edge(v)));
    CGAL_precondition(cusp1 != 0);
    // ------------------------------------------------------------------------- 
    // right-left bitangents
    if (!is_left_xx(v)) return (is_xx_left(cusp1)) ? ccw_source_edge(cusp1) :
						      skip.backward(cw_target_edge(cusp1));
    // ------------------------------------------------------------------------- 
    // left-left bitangents
    Vertex_handle cusp2 = sup(ul(cw_source_edge(v)));
    CGAL_precondition(cusp2 != 0);
    // Rtwo(v) is an arc ?
    if (skip.ccR(skip.cwR(cusp2)) == cusp2) return skip.backward(cw_source_edge(cusp2)); 
    // Otherwise return edge on cusp1
    return (is_xx_left(cusp1)) ? ccw_source_edge(cusp1) :
				  skip.backward(cw_target_edge(cusp1));
    // ------------------------------------------------------------------------- 
}

// -----------------------------------------------------------------------------
// Returns the last arc of Rleave(b). 
// Constraint case not treated FIXME

template < class A_ >
template < class Tr >
typename Chi2_flip_traits<A_>::Edge_handle
Chi2_flip_traits<A_>::second_cusp_edge(Vertex_handle v, Tr ) const 
{
    // ------------------------------------------------------------------------- 
    // Operators used by this method
    typename Tr::Ul ul; 
    typename Tr::Sup sup; typename Tr::Inf inf; 
    typename Tr::CwR cwR; 
    typename Tr::Cw_edge         cw_edge;
    typename Tr::Cw_source_edge  cw_source_edge;
    typename Tr::Cw_target_edge  cw_target_edge;
    typename Tr::Is_left_xx      is_left_xx;
    typename Tr::Source_object   source_object;
    Skip<Tr> skip;
    // ------------------------------------------------------------------------- 
    // left-xx bitangents
    if (is_left_xx(v)) {
	Vertex_handle cusp2 = sup(ul(cw_source_edge(v)));
	CGAL_precondition(cusp2 != 0);
	return skip.backward(cw_source_edge(cusp2));
    }
    // ------------------------------------------------------------------------- 
    // right-xx bitangents
    Edge_handle e = cw_edge(cwR(v)->pi(),source_object(v));
    return e;
//     return (e == 0 || inf(e) == 0) ? cw_target_edge(v->pi()) : e;
    // ------------------------------------------------------------------------- 
}

// -----------------------------------------------------------------------------
// Returns the last arc of Rleave(b). 
// Constraint case not treated FIXME

template < class A_ >
template < class Tr >
typename Chi2_flip_traits<A_>::Edge_handle
Chi2_flip_traits<A_>::third_cusp_edge(Vertex_handle v, Tr ) const 
{
    // ------------------------------------------------------------------------- 
    // Operators used by this method
    typename Tr::Inf inf; 
    typename Tr::CwR cwR;
    typename Tr::Cw_edge cw_edge;
    typename Tr::Cw_target_edge  cw_target_edge;
    typename Tr::Source_object source_object;
    // ------------------------------------------------------------------------- 
    Edge_handle e = cw_edge(cwR(v)->pi(),source_object(v));
    return (e == 0 || inf(e) == 0) ? cw_target_edge(v->pi()) : e;
    // ------------------------------------------------------------------------- 
}

// -----------------------------------------------------------------------------

template < class A_ >
template < class Tr >
typename Chi2_flip_traits<A_>::Edge_handle
Chi2_flip_traits<A_>::start_edge(Vertex_handle v, Tr) const 
{
    // ------------------------------------------------------------------------- 
    // Operators used by this method
    typename Tr::Dr dr; typename Tr::Ul ul; typename Tr::Ur ur;
    typename Tr::Inf inf; typename Tr::Sup sup; 
//     typename Tr::CcR ccR; typename Tr::CcL ccL;
//     typename Tr::CwR cwR; typename Tr::CwL cwL;
    typename Tr::Ccw_target_edge ccw_target_edge;
    typename Tr::Ccw_source_edge ccw_source_edge;
    typename Tr::Cw_target_edge cw_target_edge;
    typename Tr::Cw_source_edge cw_source_edge;
    typename Tr::Is_left_xx is_left_xx;
    typename Tr::Is_xx_left is_xx_left;
    typename Tr::Less_vertex_handle chi2;
    typename Tr::Source_object source_object;
    typename Tr::Target_object target_object;
    Skip<Tr> skip;
    // ------------------------------------------------------------------------- 
//     while (is_xx_left(v)&&source_object(v)==target_object(cwR(v))) v=ccR(v); return ccw_source_edge(v);

    // If v is xx-right we have access to the first cusp point directly
    if (!is_xx_left(v)) return ccw_target_edge(v);
    // ------------------------------------------------------------------------- 
    // The minimal bitangent of the second side of R(b) is the sink of 
    // the face face dr(cw_source_edge(v)).
    Vertex_handle cusp1 = sup(dr(cw_source_edge(v)));
    CGAL_precondition(cusp1 != 0);
    
//      return ccw_source_edge(cusp1);
    Vertex_handle tmp=0;


    // ------------------------------------------------------------------------- 
    // For a left-xx bitangent we have access to the first bitangent of the
    // third side of R(b). we use this information to see if we can jump
    // directly to the second cusp point.
    if (is_left_xx(v)) {
	Vertex_handle cusp2 = sup(ul(cw_source_edge(v)));
	if (!cusp2->is_constraint()) {
	    if (!is_left_xx(cusp2) && (a->is_on_convex_hull(cusp1) || 
				       chi2(cusp2,cusp1))) {
              while (tmp=skip.cwR(cusp2),!tmp->is_constraint() &&
                      (sup(tmp) != 0 || inf(tmp) == 0)) 
// 		    cusp2 = cwR(cusp2);
                cusp2=tmp;
		// First case  : Rtwo(b) is an arc
		if (skip.ccR(tmp) == cusp2) 
		    return skip.backward(cw_source_edge(cusp2));
		// Second case : Rtwo(b) contains only constraints 
		cusp1 = v;	
		while (tmp=skip.ccR(cusp1),(is_xx_left(cusp1) && 
                                            skip.cwR(tmp) == cusp1))
		    cusp1 = tmp;
		return (is_xx_left(cusp1)) ? ccw_source_edge(tmp):
					     ccw_target_edge(cusp1);
	    }
	    else if (check_tag(typename Tr::Chi3::is_valid())&&is_left_xx(cusp2)) {
		// Compute the maximal non constraint bitangent of Rtwo(b)
		// to apply the chi3 predicate.
              while (tmp=skip.cwR(cusp2),!skip.is_in_pt(cusp2)||
                     (tmp->is_constraint() && skip.ccR(tmp) == cusp2))
		    cusp2 = tmp; // constraints of Rthr(b)
//                 typename A_::Disk_handle d=source_object(cusp2);
// 		cusp2 = cwR(cusp2); // maximal bitangent of Rtwo(b)
                cusp2=tmp;
		while (cusp2->is_constraint() && 
                       (tmp=skip.cwL(cusp2),(skip.ccL(tmp) == cusp2))
                      //   : target_object(cusp2)==d
                       )
                  cusp2 = tmp;// cwL(cusp2); // constraints of Rtwo(b)
		typename Tr::Chi3 chi3;
		if (chi3(cusp2,v)) return ccw_source_edge(cusp2);
	    }
	}
    }
    // ------------------------------------------------------------------------- 
    // We have no clue where to start the walk of R(b), so we simply walk 
    // on Rone(v) till the first cusp point. This is NOT efficient !
    // I should probably maintain a pointer to this cusp point
    if (cusp1->is_constraint()) { 
	cusp1 = v;	
	while (is_xx_left(cusp1) && (tmp=skip.ccR(cusp1),skip.cwR(tmp) == cusp1))
          cusp1=tmp;
// 	    cusp1 = ccR(cusp1);
	return (is_xx_left(tmp)) ? ccw_source_edge((tmp)):
				     ccw_target_edge(tmp);
    }
    // ------------------------------------------------------------------------- 
    // if cusp1 is xx-left, jump directly to its tail
    if (is_xx_left(cusp1)) return ccw_source_edge(cusp1);
    // ------------------------------------------------------------------------- 
    // if the arc preceeding cusp1 is a cusp arc return this arc
    if (skip.ccR(skip.cwL(cusp1)) == cusp1)   return skip.backward(cw_target_edge(cusp1));
    // ------------------------------------------------------------------------- 
    // At this point we are left with the case where the arc preceeding cusp1 is
    // regular. 
    // Let vR be the successor of v in R(v) amongst the non-constraint vertices.
    Vertex_handle vR = sup(ur(cw_target_edge(v)));
    CGAL_precondition(vR != 0);
    // If vR is a constraint then vL defines the second cusp point.
    if (vR->is_constraint()) return (!is_xx_left(vR)) ? ccw_target_edge(vR):
						        ccw_source_edge(vR);
    // Otherwise if pi(phi_*(vR)) does not leave the first atom of Rtwo(v) then
    // phi(v) doesn't either. Jump to the source edge of cusp1.
    CGAL_precondition(inf(vR) != 0);
    Vertex_handle phisvR = inf(inf(vR));
    CGAL_precondition(phisvR != 0);
    if (sup(cusp1) == 0 && chi2(phisvR,cusp1)) 
	return ccw_source_edge(cusp1);
    // ------------------------------------------------------------------------- 
    // Now we know that pi(phi_*(vR)) enters the first atom of Rtwo(v). 
    // The first atom of Rtwo(v) might be spliced by several bitangents in
    // pi(G_*).  Traverse all these bitangents in G_* (this will be done only
    // once).
    Vertex_handle w = skip.cwL(cusp1);
    while (!w->is_constraint() && source_object(w) == target_object(cusp1)) 
	w = skip.cwR(w);
    // Case 1: w is regular
    if (!w->is_constraint()) return ccw_target_edge(w);
    // Case 2: w is a constraint. Walk back on constraints till we find the cusp
    // point.
    while (!is_xx_left(w) && (tmp=skip.cwL(w),skip.ccL(tmp) == w)) w = tmp;
    return (!is_xx_left(w)) ? skip.backward(cw_target_edge(w)) : ccw_source_edge(w);
    // ------------------------------------------------------------------------- 
}

// -----------------------------------------------------------------------------

template < class A_ >
template < class Tr >
typename Chi2_flip_traits<A_>::Vertex_handle
Chi2_flip_traits<A_>::walk(Edge_handle r , Edge_handle l, Tr) const 
{
    // ------------------------------------------------------------------------- 
    // Operators used by this method
    typename Tr::Less_vertex_handle chi2; 
    typename Tr::Sup sup; 
    typename Tr::Inf inf; 
    typename Tr::Splice splice; 


    // ------------------------------------------------------------------------- 
    typename Tr::Right_traits right_traits=typename Tr::Right_traits();
    typename Tr::Left_traits  left_traits=typename Tr::Left_traits();
    Skip<typename Tr::Right_traits> skipr;
    Skip<typename Tr::Left_traits> skipl;
    // ------------------------------------------------------------------------- 
    Vertex_handle phiv = bit(r,l, Tr());
    bool bbb=false;
    bbb=sup(r)==phiv||inf(r)==phiv||sup(l)==phiv||inf(l)==phiv;
    bool advance_r = advance(r,phiv, right_traits);
    bool advance_l = advance(l,phiv, left_traits);
    while (advance_r || advance_l) {
	// If we have to advance on both sides we use the chi2 predicate to
	// choose.
	if (advance_r && advance_l) advance_r = chi2(sup(r),sup(l));
	if (advance_r) {
          r=skipr.next_edge(r);

// 	    if (r->sign()||source_object(sup(r)) == r->object())
// 		 r = ccw_target_edge(sup(r));
// 	    else r = ccw_source_edge(sup(r));
	}
	else {
// 	    if (!l->sign()||target_object(sup(l)) == l->object())
// 		 l = ccw_source_edge(sup(l));
// 	    else l = ccw_target_edge(sup(l));
          l=skipl.next_edge(l);
	}
        if (!bbb) delete phiv;
	phiv = bit(r,l, Tr());
        bbb=sup(r)==phiv||inf(r)==phiv||sup(l)==phiv||inf(l)==phiv;
	advance_r = advance(r,phiv, right_traits);
	advance_l = advance(l,phiv, left_traits);
    }
    // ------------------------------------------------------------------------- 
    splice(r,phiv); splice(l,phiv); 
    return phiv;
    // ------------------------------------------------------------------------- 
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

template < class A_ >
class PhiR_flip_traits : public Chi2_flip_traits<A_>
{
public:
    // -------------------------------------------------------------------------
    typedef A_ Antichain;
    typedef typename Antichain::Vertex_handle Vertex_handle;
    typedef typename Antichain::Edge_handle   Edge_handle;
    typedef typename Antichain::Edge          Edge;
    typedef typename Antichain::Face_handle   Face_handle;
    typedef Chi2_flip_traits<A_>          Base;
    // -------------------------------------------------------------------------
    // CONSTRUCTOR
    PhiR_flip_traits(/*const*/ Antichain* an) : Base(an) { }
    // -------------------------------------------------------------------------
    template < class Tr >
    Vertex_handle operator()(const Vertex_handle& v , Tr) const
    {
	CGAL_precondition(!this->a->is_on_convex_hull(v));
	typename Tr::Right_traits right;
	typename Tr::Left_traits  left;
	typename Tr::Ccw_source_edge  ccw_source_edge;
	typename Tr::Ccw_target_edge  ccw_target_edge;
	return walk(ccw_source_edge(phiR(v,right)),
		    ccw_target_edge(phiR(v,left)),Tr());
    }
    // -------------------------------------------------------------------------
protected:
    // -------------------------------------------------------------------------
    template < class Tr >
    Vertex_handle phiR(const Vertex_handle& v , Tr tr) const;
    template < class Tr >
    Vertex_handle phiR_walk(Edge_handle r , Edge_handle l , 
			    Vertex_handle phis,Tr tr) const;
    // ------------------------------------------------------------------------
    template < class Tr >
    Vertex_handle bit_with_pi(Edge_handle r, Edge_handle l,Tr tr)  const
    {
	typename Tr::Vertex_creator vertex;
	CGAL_precondition(l != 0);
	Edge_handle lpi = new Edge(*l); 
	lpi->set_sign(!l->sign());
	Vertex_handle v = vertex(r,lpi); 
        lpi->set_adjacent_faces(0,0,0);
        lpi->set_sup(0);
        lpi->set_inf(0);
        delete lpi;
	return v;
    }
    // -------------------------------------------------------------------------
};

// -----------------------------------------------------------------------------

template < class A_ >
template < class Tr >
inline typename PhiR_flip_traits<A_>::Vertex_handle
PhiR_flip_traits<A_>::
phiR_walk(Edge_handle r , Edge_handle l , Vertex_handle phis,  Tr tr) const
{
    // ------------------------------------------------------------------------- 
    // Operators used by this method
    typename Tr::Sup sup;
    typename Tr::Less_vertex_handle chi1;
    typename Tr::Ccw_target_edge ccw_target_edge;
    typename Tr::Ccw_source_edge ccw_source_edge;
    typename Tr::Set_ccw_edge    set_ccw_edge;
    typename Tr::Right_traits right_traits;
    typename Base::template Skip<typename Tr::Right_traits> skipr;
    typename Base::template Skip<typename Tr::Left_traits> skipl;
    // ------------------------------------------------------------------------- 
    Vertex_handle phiR = bit_with_pi(r,l,Tr());
    bool advance_r = advance(r,phiR, right_traits);
    bool advance_l = (!advance_r && sup(l) != phis && chi1(sup(l)->pi(),phiR));
#ifdef FLIP_STATISTICS
    if (!advance_r && sup(l) != phis) ++chi1_calls;
#endif
    while (advance_r || advance_l) {
      if (advance_r) r = skipr.next_edge(r);
      else           l = skipl.next_edge(l);
// 	if (advance_r) r = ccw_source_edge(sup(r));
// 	else           l = ccw_target_edge(sup(l));
	delete phiR; phiR = bit_with_pi(r,l,Tr());
	advance_r = advance(r,phiR, right_traits);
	advance_l = (!advance_r && sup(l) != phis && chi1(sup(l)->pi(),phiR));
#ifdef FLIP_STATISTICS
	if (!advance_r && sup(l) != phis) ++chi1_calls;
#endif
    }
    // ------------------------------------------------------------------------- 
    set_ccw_edge(phiR,r);
    return phiR;
    // ------------------------------------------------------------------------- 
}

// -----------------------------------------------------------------------------

template < class A_ >
template < class Tr >
inline typename PhiR_flip_traits<A_>::Vertex_handle
PhiR_flip_traits<A_>::phiR(const Vertex_handle& v , Tr tr) const
{
    // ------------------------------------------------------------------------- 
    // Operators used by this method
    typename Tr::Inf inf;
    typename Tr::CcR ccR; 
    typename Tr::Ccw_source_edge ccw_source_edge;
    typename Base::template Skip<Tr> skip;
    // ------------------------------------------------------------------------- 
    Edge_handle r = first_cusp_edge(v,Tr());
    Edge_handle l = ccw_source_edge(v->pi());
    Vertex_handle phis = inf(inf(skip.ccR(v)));
    // ------------------------------------------------------------------------- 
    return phiR_walk(r,l,phis,tr);
    // ------------------------------------------------------------------------- 
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

template < class A_ >
class Chi1_flip_traits : public PhiR_flip_traits<A_>
{
public:
    // -------------------------------------------------------------------------
    typedef A_ Antichain;
    typedef typename Antichain::Vertex_handle Vertex_handle;
    typedef typename Antichain::Edge_handle   Edge_handle;
    typedef typename Antichain::Edge          Edge;
    typedef typename Antichain::Face_handle   Face_handle;
    typedef typename Antichain::Right_ccw_traits Right_ccw_traits;
    typedef typename Antichain::Left_ccw_traits  Left_ccw_traits;
    typedef typename Antichain::Right_cw_traits  Right_cw_traits;
    typedef typename Antichain::Left_cw_traits   Left_cw_traits;
    typedef PhiR_flip_traits<A_>          Base;
    // -------------------------------------------------------------------------
    // CONSTRUCTOR
    Chi1_flip_traits(/*const*/ Antichain* an) : Base(an) { }
    // -------------------------------------------------------------------------
    template < class Tr >
    Vertex_handle operator()(const Vertex_handle& v , Tr) const;
    // -------------------------------------------------------------------------
private:
    // -------------------------------------------------------------------------
    template < class Tr >
    Vertex_handle chi1_walk(Edge_handle l, Edge_handle r, Tr tr,
			    bool end_known) const;
    // -------------------------------------------------------------------------
};

// -----------------------------------------------------------------------------

template < class A_ >
template < class Tr >
inline typename Chi1_flip_traits<A_>::Vertex_handle
Chi1_flip_traits<A_>::
chi1_walk(Edge_handle l , Edge_handle r, Tr tr,bool end_known = true) const
{
    // -------------------------------------------------------------------------
    // Operators used in this method
    typename Tr::Sign            sign;
    typename Tr::Splice          splice;
    typename Tr::Sup sup; typename Tr::Inf inf;
    typename Tr::Target_object target_object;
    typename Tr::Less_vertex_handle chi1;
    typename Tr::Cw_target_edge  cw_target_edge;
    typename Tr::Cw_source_edge  cw_source_edge;
    typename Tr::Ccw_target_edge ccw_target_edge;
    typename Tr::Ccw_source_edge ccw_source_edge;
    typename Tr::Bottom_edge     bottom_edge;
    typename Tr::Right_traits right_traits;
    typename Tr::Left_cw_traits left_cw_traits;
    typename Tr::Right_cw_traits right_cw_traits;
    typename Base::template Skip<typename Tr::Right_traits> skipr;
    typename Base::template Skip<typename Tr::Left_traits> skipl;
    // -------------------------------------------------------------------------
    // If l is negative then l is the first atom of Lthr(v)
    if ( !sign(l) || l == cw_target_edge(sup(l)) ) return walk(r,l,Tr());
    // -------------------------------------------------------------------------
    Edge_handle lmin = l; // save initial value of l
    // -------------------------------------------------------------------------
    // Else compute psi(l)
    Face_handle   f    = inf(sup(l));
    Vertex_handle phis = inf(f); // sentinel when walking on G_*
    Edge_handle lmax = (end_known) ?  cw_source_edge(phis) : 
				     ccw_target_edge(phis);
    Vertex_handle min  = inf(bottom_edge(f));
    if (target_object(min) == l->object())  
	 l = third_cusp_edge(min,right_cw_traits);
    else l = second_cusp_edge(min,left_cw_traits);
    // -------------------------------------------------------------------------
    Vertex_handle phiR = bit_with_pi(r,l,Tr());
    bool advance_r = advance(r,phiR,right_traits);
    bool advance_l = (!advance_r && l != lmax && chi1(sup(l)->pi(),phiR));
#ifdef FLIP_STATISTICS
    if (!advance_r && l != lmax) ++chi1_calls;
#endif
    while (advance_r || advance_l) {
// 	if (advance_r) r = ccw_source_edge(sup(r));
// 	else           l = ccw_target_edge(sup(l));
      if (advance_r) r = skipr.next_edge(r);
      else           l = skipr.next_edge(l);
	if (l == ccw_target_edge(phis)) break;
	delete phiR; phiR = bit_with_pi(r,l,Tr());
	advance_r = advance(r,phiR,right_traits);
	advance_l = (!advance_r && l != lmax && chi1(sup(l)->pi(),phiR));
#ifdef FLIP_STATISTICS
	if (!advance_r && l != lmax) ++chi1_calls;
#endif
    }
    // -------------------------------------------------------------------------
    // Special case when phiR(sup(l)) leaves first atom of Lthr(v)
    if (l == skipr.forward(ccw_target_edge(phis))) {
	CGAL_precondition(end_known == false);
	l = lmin;
	while (sign(l)) l = skipr.next_edge(l);//ccw_target_edge(sup(l));
	sup(lmin)->set_phiR(l);
	return walk(r,l,Tr()); // phi(v) is xx-right enters Lthr(v)
    }
    // -------------------------------------------------------------------------
    if (!sign(l)) { 
	CGAL_precondition(l->object() == lmin->object());
	Vertex_handle phiv = bit(r,lmin,Tr());
#ifdef FLIP_STATISTICS
	++chi1_calls;
#endif
	if (chi1(phiv,sup(lmin))) {
	    splice(lmin,phiv); splice(r,phiv);
	    return phiv;
	}
	else delete phiv;
    }
    sup(lmin)->set_phiR(r);
    return 0;
    // -------------------------------------------------------------------------
}

// -----------------------------------------------------------------------------

template < class A_ >
template < class Tr >
inline typename Chi1_flip_traits<A_>::Vertex_handle
Chi1_flip_traits<A_>::operator()(const Vertex_handle& v , Tr tr) const
{
    // -------------------------------------------------------------------------
    // Operators used in this method
    typename Tr::Is_xx_left      is_xx_left;
    typename Tr::Is_left_xx      is_left_xx;
    typename Tr::Sup sup; typename Tr::Inf inf;
    typename Tr::Less_vertex_handle chi1;
    typename Tr::Vertex_creator  vertex;
    typename Tr::Splice          splice;
    typename Tr::Ccw_target_edge ccw_target_edge;
    typename Tr::Ccw_source_edge ccw_source_edge;
    typename Tr::Cw_target_edge  cw_target_edge;
    typename Tr::Cw_source_edge  cw_source_edge;
    typename Tr::CcR ccR; typename Tr::CcL ccL;
    typename Tr::CwR cwR; typename Tr::CwL cwL;
    typename Tr::Right_traits right_traits;
    typename Tr::Left_traits  left_traits;
    typename Tr::Right_cw_traits  right_cw_traits;
    typename Tr::Left_cw_traits   left_cw_traits;
    typename Base::template Skip<typename Tr::Right_traits> skipr;
    typename Base::template Skip<typename Tr::Left_traits> skipl;
    // -------------------------------------------------------------------------
    // Starting point on Rtwo(v)
    Edge_handle r = 0;
    Vertex_handle tmp=0;
    if (is_xx_left(v) && (tmp=skipr.ccR(v),skipr.cwR(tmp) == v)) {
	if (tmp->phiR() == 0) {
	    tmp->set_phiR_vertex(phiR(v,right_traits));
	    tmp->set_phiR(ccw_source_edge(tmp->phiR_vertex()));
	}
	r = tmp->phiR();
    }
    else r = first_cusp_edge(v,right_traits);
    CGAL_precondition(r != 0);
    // -------------------------------------------------------------------------
    // Starting point on Rtwo(v)
    Edge_handle l = 0;
    if (!is_left_xx(v) && (tmp=skipr.ccL(v),skipr.cwL(tmp) == v)) { 
	if (tmp->phiL() == 0) {
	    tmp->set_phiL_vertex(phiR(v,left_traits));
	    tmp->set_phiL(ccw_target_edge(tmp->phiL_vertex()));
	}
	l = tmp->phiL();
    }
    else l = first_cusp_edge(v,left_traits);
    CGAL_precondition(l != 0);
    // -------------------------------------------------------------------------
    Vertex_handle phiv = 0; // the vertex we return
    // -------------------------------------------------------------------------
    // If l is the first atom of Ltwo(v)
    if (ccL(inf(l)) == sup(l)) {
	if (v->phiR() == 0) {
	    Edge_handle l;
            tmp=cwR(v);
	    if (!is_left_xx(v) && (skipr.ccR(tmp) == v))
		 l = skipr.backward(cw_source_edge(tmp));
	    else if (skipr.ccR(tmp) == v)
		 l = third_cusp_edge(tmp,right_cw_traits);
	    else l = second_cusp_edge(tmp,left_cw_traits);
	    v->set_phiR_vertex(phiR_walk(r,l,inf(inf(v)),right_traits));
	    v->set_phiR(ccw_source_edge(v->phiR_vertex()));
	    CGAL_precondition(v->phiR() != 0);
	}
	CGAL_precondition(v->phiR() != 0);
	phiv = vertex(v->phiR(),l);
#ifdef FLIP_STATISTICS
	++chi1_calls;
#endif
	if (chi1(phiv,sup(l))) { 
#ifdef FLIP_STATISTICS
	    ++phiR_equal_phi;
#endif
	    splice(l,phiv); splice(v->phiR(),phiv); 
	}
	else {
	    delete phiv; phiv = 0;
	    l = skipr.next_edge(l);//ccw_target_edge(sup(l));
	}
    }
    // -------------------------------------------------------------------------
    // If phiv was not found at the previuos step walk by computing successive
    // phiR.
    if (phiv == 0) {
	Edge_handle lmin = l;
	phiv = chi1_walk(l,r,tr,false);
	while ( phiv == 0 ) {
          l = skipr.next_edge(l);//ccw_target_edge(sup(l));
	    phiv = chi1_walk(l,r,tr,false);
#ifdef FLIP_STATISTICS
	    ++psi_computed;
#endif
	}
	while (lmin != l) {
	    if (sup(lmin)->phiR() == cw_source_edge(phiv))
		sup(lmin)->set_phiR(ccw_source_edge(phiv));
	    lmin = skipr.next_edge(lmin);//ccw_target_edge(sup(lmin));
	}
    }
    // -------------------------------------------------------------------------
    // Update of the phiR operator.
    if (phiv == ccR(v)) {
	phiv->set_phiR(cw_source_edge(v->pi()));
	phiv->set_phiR_vertex(v->pi());
    }
    else if (!is_xx_left(v) || cwR(ccR(v)) != v) { // regular arc already done
	Vertex_handle phisvR = inf(inf(ccR(v)));
	CGAL_precondition(phisvR != 0);
#ifdef FLIP_STATISTICS
	if (v->phiL() == 0 && is_left_xx(phiv)) ++chi1_calls;
#endif
	if (is_left_xx(phiv) && chi1(phisvR->pi(),phiv)){
#ifdef FLIP_STATISTICS
	     ++phiL_equal_phi;
#endif
	     ccR(v)->set_phiL(cw_source_edge(phiv));
	     ccR(v)->set_phiL_vertex(phisvR->pi());
	}
	else if (v->phiL() == cw_target_edge(phiv))
	     ccR(v)->set_phiL(ccw_target_edge(phiv));
	else {
	    ccR(v)->set_phiL(v->phiL());
	    ccR(v)->set_phiL_vertex(v->phiL_vertex());
	}
    }
    // -------------------------------------------------------------------------
    // Update of the phiL operator.
    if (phiv == ccL(v)) {
	phiv->set_phiL(cw_target_edge(v->pi()));
	phiv->set_phiL_vertex(v->pi());
    }
    else if (is_left_xx(v) || cwL(ccL(v)) != v) {
	Vertex_handle phisvL = inf(inf(ccL(v)));
	CGAL_precondition(phisvL != 0);
#ifdef FLIP_STATISTICS
	if (v->phiR() == 0 && !is_xx_left(phiv)) ++chi1_calls;
#endif
	if (!is_xx_left(phiv) && chi1(phisvL->pi(),phiv)) {
#ifdef FLIP_STATISTICS
	     ++phiR_equal_phi;
#endif
	     ccL(v)->set_phiR(cw_target_edge(phiv));
	     ccL(v)->set_phiR_vertex(phisvL->pi());
	}
	else if (v->phiR() == cw_source_edge(phiv))
	     ccL(v)->set_phiR(ccw_source_edge(phiv));
	else {
	    ccL(v)->set_phiR(v->phiR());
	    ccL(v)->set_phiR_vertex(v->phiR_vertex());
	}
    }
    // -------------------------------------------------------------------------
    return phiv;
    // -------------------------------------------------------------------------
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

template < class A_ >
class Bw_flip_traits : public Chi2_flip_traits<A_>
{
public:
    // -------------------------------------------------------------------------
    typedef A_ Antichain;
    typedef typename Antichain::Vertex_handle Vertex_handle;
    typedef typename Antichain::Edge_handle   Edge_handle;
    typedef typename Antichain::Edge          Edge;
    typedef typename Antichain::Face_handle   Face_handle;
    typedef Chi2_flip_traits<A_>          Base;
    // -------------------------------------------------------------------------
    // CONSTRUCTOR
    Bw_flip_traits(/*const*/ Antichain* an) : Base(an) { }
    // -------------------------------------------------------------------------
    template < class Tr >
    Vertex_handle operator()(const Vertex_handle& v , Tr) const
    {
	CGAL_precondition(!this->a->is_on_convex_hull(v));
	typename Tr::Right_traits right;
	typename Tr::Left_traits  left;
	return walk_backward<Tr>(second_cusp_edge(v,right),
				  second_cusp_edge(v,left));
			    
    }
    // -------------------------------------------------------------------------
private:
    template < class Tr >
    Vertex_handle walk_backward(const Edge_handle& right,
			        const Edge_handle& left) const;
    template < class Tr >
    bool advance(Edge_handle r, Vertex_handle phiv) const;
};

// -----------------------------------------------------------------------------
// Returns false if phiv is tangent to r

template < class A_ >
template < class Tr >
inline bool
Bw_flip_traits<A_>::advance(Edge_handle r, Vertex_handle phiv) const
{
    typename Tr::Inf inf; 
    typename Tr::Less_vertex_handle chi1; 
    typename Tr::Target_object target_object; 
    typename Tr::Is_xx_left is_xx_left;
    // Cases when r is first edge
    if (is_xx_left(inf(r))) return false; 
    if (target_object(inf(r)) == r->object()) return false;
    // Else compute chi1
#ifdef FLIP_STATISTICS
    ++chi1_calls;
#endif
    return chi1(phiv,inf(r));
}

// -----------------------------------------------------------------------------

template < class A_ >
template < class Tr >
typename Bw_flip_traits<A_>::Vertex_handle
Bw_flip_traits<A_>::walk_backward(const Edge_handle& right , 
					const Edge_handle& left) const 
{
    // ------------------------------------------------------------------------- 
    // Operators used by this method
    typename Tr::Cw_target_edge cw_target_edge;
    typename Tr::Cw_source_edge cw_source_edge;
    typename Tr::Inf inf; 
    typename Tr::Splice splice; 
    // ------------------------------------------------------------------------- 
    typename Tr::Right_traits right_traits;
    typename Tr::Left_traits  left_traits;
    // ------------------------------------------------------------------------- 
    Edge_handle r = right; Edge_handle l = left;
    Vertex_handle phiv = bit(r,l,Tr());
    bool advance_r = advance(r,phiv,right_traits);
    bool advance_l = (!advance_r && advance(l,phiv,left_traits));
    while (advance_r || advance_l) {
	if (advance_r) r = cw_target_edge(inf(r));
	else           l = cw_source_edge(inf(l));
	delete phiv; phiv = bit(r,l,Tr());
	advance_r = advance(r,phiv,right_traits);
	advance_l = (!advance_r && advance(l,phiv,left_traits));
    }
    // ------------------------------------------------------------------------- 
    splice(r,phiv); splice(l,phiv); 
    return phiv;
    // ------------------------------------------------------------------------- 
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------


// -----------------------------------------------------------------------------

struct Flip_traits {
    template < class A_ >
    struct Flip_wrapper {
	typedef Chi2_flip_traits<A_> Flip_traits;
    };
};

// -----------------------------------------------------------------------------

struct Phir_flip_traits {
    template < class A_ >
    struct Flip_wrapper {
	typedef PhiR_flip_traits<A_> Flip_traits;
    };
};

// -----------------------------------------------------------------------------

struct Matroidal_flip_traits {
    template < class A_ >
    struct Flip_wrapper {
	typedef Chi1_flip_traits<A_> Flip_traits;
    };
};

// -----------------------------------------------------------------------------

struct Backward_flip_traits {
    template < class A_ >
    struct Flip_wrapper {
	typedef Bw_flip_traits<A_> Flip_traits;
    };
};

// -----------------------------------------------------------------------------

}
CGAL_END_NAMESPACE

#endif // VISIBILITY_COMPLEX_2_FLIP_TRAITS_H
