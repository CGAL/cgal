#ifndef VISIBILITY_COMPLEX_FLIP_TRAITS_H
#define VISIBILITY_COMPLEX_FLIP_TRAITS_H

CGAL_BEGIN_NAMESPACE

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

template < class _A >
class Chi2_flip_traits 
{
public:
    // -------------------------------------------------------------------------
    typedef _A Antichain;
    typedef typename Antichain::Vertex_handle Vertex_handle;
    typedef typename Antichain::Edge_handle   Edge_handle;
    typedef typename Antichain::Face_handle   Face_handle;
    // -------------------------------------------------------------------------
    // CONSTRUCTOR
    Chi2_flip_traits(/*const*/ Antichain* an) : a(an) { }
    // -------------------------------------------------------------------------
    template < class _Tr >
    Edge_handle   first_cusp_edge(Vertex_handle v, _Tr)    const;
    template < class _Tr >
    Edge_handle   second_cusp_edge(Vertex_handle v, _Tr)   const;
    template < class _Tr >
    Edge_handle   third_cusp_edge(Vertex_handle v, _Tr)    const;
    template < class _Tr >
    Edge_handle   start_edge(Vertex_handle v, _Tr)         const;
    template < class _Tr >
    Vertex_handle walk(Edge_handle r, Edge_handle l, _Tr)  const;
    // -------------------------------------------------------------------------
    template < class _Tr >
    Vertex_handle operator()(const Vertex_handle& v , _Tr tr) const
    {
	typename _Tr::Right_traits right;
	typename _Tr::Left_traits  left;
	return walk(start_edge(v,right),start_edge(v,left),tr);
    }
    // -------------------------------------------------------------------------
public:
    // -------------------------------------------------------------------------
    // Returns false if phiv is tangent to r
    template < class _Tr >
    bool advance(Edge_handle r, Vertex_handle phiv, const _Tr&) const;
    // -------------------------------------------------------------------------
    // Returns the vertex leaving the support disk of r and entering the support
    // disk of l. Check the equality to sup(r),inf(r),sup(l),inf(l) to avoid
    // duplicates
    template < class _Tr >
    Vertex_handle bit(Edge_handle r, Edge_handle l, const _Tr&) const;
  //template < class _Tr >
  //Vertex_handle bit(Edge_handle r, Edge_handle l, _Tr) const
  //{ return bit(r,l, _Tr()); }
    // -------------------------------------------------------------------------
protected:
    /*const*/ Antichain* a; 
};

// -----------------------------------------------------------------------------

template < class _A >
template < class _Tr >
inline bool
Chi2_flip_traits<_A>::advance(Edge_handle r, Vertex_handle phiv, const _Tr&) const
{
    typename _Tr::Sup sup; typename _Tr::Sign sign;
    typename _Tr::Less_vertex_handle chi1; 
    if (sup(r) == 0 || (sign(r) && !sup(r)->is_constraint())) return false;
    //typename _Tr::Target_object target_object; 
    //if (target_object(sup(r)) != r->object()) return false;
#ifdef FLIP_STATISTICS
    ++chi1_calls;
    if (chi1(phiv,sup(r))) ++chi1_naive_calls;
#endif
    return chi1(sup(r),phiv);
}

// -----------------------------------------------------------------------------

template < class _A >
template < class _Tr >
inline 
typename Chi2_flip_traits<_A>::Vertex_handle
Chi2_flip_traits<_A>::bit(Edge_handle r, Edge_handle l, const _Tr&) const
{
    typename _Tr::Sup sup; typename _Tr::Inf inf;
    typename _Tr::Vertex_creator vertex;
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

template < class _A >
template < class _Tr >
typename Chi2_flip_traits<_A>::Edge_handle
Chi2_flip_traits<_A>::first_cusp_edge(Vertex_handle v, _Tr ) const 
{
    // ------------------------------------------------------------------------- 
    // Operators used by this method
    typename _Tr::Dr dr; typename _Tr::Ul ul; 
    typename _Tr::Sup sup; 
    typename _Tr::CcR ccR; typename _Tr::CwR cwR; 
    typename _Tr::Ccw_source_edge ccw_source_edge;
    typename _Tr::Ccw_target_edge ccw_target_edge;
    typename _Tr::Cw_target_edge  cw_target_edge;
    typename _Tr::Cw_source_edge  cw_source_edge;
    typename _Tr::Is_left_xx is_left_xx;
    typename _Tr::Is_xx_left is_xx_left;
    // ------------------------------------------------------------------------- 
    // xx-right bitangents
    if (!is_xx_left(v)) return ccw_target_edge(v);
    // ------------------------------------------------------------------------- 
    Vertex_handle cusp1 = sup(dr(cw_source_edge(v)));
    CGAL_precondition(cusp1 != 0);
    // ------------------------------------------------------------------------- 
    // right-left bitangents
    if (!is_left_xx(v)) return (is_xx_left(cusp1)) ? ccw_source_edge(cusp1) :
						      cw_target_edge(cusp1);
    // ------------------------------------------------------------------------- 
    // left-left bitangents
    Vertex_handle cusp2 = sup(ul(cw_source_edge(v)));
    CGAL_precondition(cusp2 != 0);
    // Rtwo(v) is an arc ?
    if (ccR(cwR(cusp2)) == cusp2) return cw_source_edge(cusp2); 
    // Otherwise return edge on cusp1
    return (is_xx_left(cusp1)) ? ccw_source_edge(cusp1) :
				  cw_target_edge(cusp1);
    // ------------------------------------------------------------------------- 
}

// -----------------------------------------------------------------------------
// Returns the last arc of Rleave(b). 
// Constraint case not treated FIXME

template < class _A >
template < class _Tr >
typename Chi2_flip_traits<_A>::Edge_handle
Chi2_flip_traits<_A>::second_cusp_edge(Vertex_handle v, _Tr ) const 
{
    // ------------------------------------------------------------------------- 
    // Operators used by this method
    typename _Tr::Ul ul; 
    typename _Tr::Sup sup; typename _Tr::Inf inf; 
    typename _Tr::CwR cwR; 
    typename _Tr::Cw_edge         cw_edge;
    typename _Tr::Cw_source_edge  cw_source_edge;
    typename _Tr::Cw_target_edge  cw_target_edge;
    typename _Tr::Is_left_xx      is_left_xx;
    typename _Tr::Source_object   source_object;
    // ------------------------------------------------------------------------- 
    // left-xx bitangents
    if (is_left_xx(v)) {
	Vertex_handle cusp2 = sup(ul(cw_source_edge(v)));
	CGAL_precondition(cusp2 != 0);
	return cw_source_edge(cusp2);
    }
    // ------------------------------------------------------------------------- 
    // right-xx bitangents
    Edge_handle e = cw_edge(cwR(v)->pi(),source_object(v));
    return (e == 0 || inf(e) == 0) ? cw_target_edge(v->pi()) : e;
    // ------------------------------------------------------------------------- 
}

// -----------------------------------------------------------------------------
// Returns the last arc of Rleave(b). 
// Constraint case not treated FIXME

template < class _A >
template < class _Tr >
typename Chi2_flip_traits<_A>::Edge_handle
Chi2_flip_traits<_A>::third_cusp_edge(Vertex_handle v, _Tr ) const 
{
    // ------------------------------------------------------------------------- 
    // Operators used by this method
    typename _Tr::Inf inf; 
    typename _Tr::CwR cwR;
    typename _Tr::Cw_edge cw_edge;
    typename _Tr::Cw_target_edge  cw_target_edge;
    typename _Tr::Source_object source_object;
    // ------------------------------------------------------------------------- 
    Edge_handle e = cw_edge(cwR(v)->pi(),source_object(v));
    return (e == 0 || inf(e) == 0) ? cw_target_edge(v->pi()) : e;
    // ------------------------------------------------------------------------- 
}

// -----------------------------------------------------------------------------

template < class _A >
template < class _Tr >
typename Chi2_flip_traits<_A>::Edge_handle
Chi2_flip_traits<_A>::start_edge(Vertex_handle v, _Tr ) const 
{
    // ------------------------------------------------------------------------- 
    // Operators used by this method
    typename _Tr::Dr dr; typename _Tr::Ul ul; typename _Tr::Ur ur;
    typename _Tr::Inf inf; typename _Tr::Sup sup; 
    typename _Tr::CcR ccR; typename _Tr::CcL ccL;
    typename _Tr::CwR cwR; typename _Tr::CwL cwL;
    typename _Tr::Ccw_target_edge ccw_target_edge;
    typename _Tr::Ccw_source_edge ccw_source_edge;
    typename _Tr::Cw_target_edge cw_target_edge;
    typename _Tr::Cw_source_edge cw_source_edge;
    typename _Tr::Is_left_xx is_left_xx;
    typename _Tr::Is_xx_left is_xx_left;
    typename _Tr::Less_vertex_handle chi2;
    typename _Tr::Source_object source_object;
    typename _Tr::Target_object target_object;
    // ------------------------------------------------------------------------- 
    // If v is xx-right we have access to the first cusp point directly
    if (!is_xx_left(v)) return ccw_target_edge(v);
    // ------------------------------------------------------------------------- 
    // The minimal bitangent of the second side of R(b) is the sink of 
    // the face face dr(cw_source_edge(v)).
    Vertex_handle cusp1 = sup(dr(cw_source_edge(v)));
    CGAL_precondition(cusp1 != 0);
    // ------------------------------------------------------------------------- 
    // For a left-xx bitangent we have access to the first bitangent of the
    // third side of R(b). we use this information to see if we can jump
    // directly to the second cusp point.
    if (is_left_xx(v)) {
	Vertex_handle cusp2 = sup(ul(cw_source_edge(v)));
	if (!cusp2->is_constraint()) {
	    if (!is_left_xx(cusp2) && (a->is_on_convex_hull(cusp1) || 
				       chi2(cusp2,cusp1))) {
		while (!cwR(cusp2)->is_constraint() &&
		       (sup(cwR(cusp2)) != 0 || inf(cwR(cusp2)) == 0)) 
		    cusp2 = cwR(cusp2);
		// First case  : Rtwo(b) is an arc
		if (ccR(cwR(cusp2)) == cusp2) 
		    return cw_source_edge(cusp2);
		// Second case : Rtwo(b) contains only constraints 
		cusp1 = v;	
		while (is_xx_left(cusp1) && cwR(ccR(cusp1)) == cusp1)
		    cusp1 = ccR(cusp1);
		return (is_xx_left(cusp1)) ? ccw_source_edge(ccR(cusp1)):
					     ccw_target_edge(cusp1);
	    }
	    else if (is_left_xx(cusp2)) {
		// Compute the maximal non constraint bitangent of Rtwo(b)
		// to apply the chi3 predicate.
		while (cwR(cusp2)->is_constraint() && ccR(cwR(cusp2)) == cusp2)
		    cusp2 = cwR(cusp2); // constraints of Rthr(b)
		cusp2 = cwR(cusp2); // maximal bitangent of Rtwo(b)
		while (cusp2->is_constraint() && ccL(cwL(cusp2)) == cusp2)
		    cusp2 = cwL(cusp2); // constraints of Rtwo(b)
		typename _Tr::Chi3 chi3;
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
	while (is_xx_left(cusp1) && cwR(ccR(cusp1)) == cusp1)
	    cusp1 = ccR(cusp1);
	return (is_xx_left(cusp1)) ? ccw_source_edge(ccR(cusp1)):
				     ccw_target_edge(cusp1);
    }
    // ------------------------------------------------------------------------- 
    // if cusp1 is xx-left, jump directly to its tail
    if (is_xx_left(cusp1)) return ccw_source_edge(cusp1);
    // ------------------------------------------------------------------------- 
    // if the arc preceeding cusp1 is a cusp arc return this arc
    if (ccR(cwL(cusp1)) == cusp1)   return cw_target_edge(cusp1);
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
    if (sup(cusp1) == 0 && chi2(cusp1,phisvR->pi())) 
	return ccw_source_edge(cusp1);
    // ------------------------------------------------------------------------- 
    // Now we know that pi(phi_*(vR)) enters the first atom of Rtwo(v). 
    // The first atom of Rtwo(v) might be spliced by several bitangents in
    // pi(G_*).  Traverse all these bitangents in G_* (this will be done only
    // once).
    Vertex_handle w = cwL(cusp1);
    while (!w->is_constraint() && source_object(w) == target_object(cusp1)) 
	w = cwR(w);
    // Case 1: w is regular
    if (!w->is_constraint()) return ccw_target_edge(w);
    // Case 2: w is a constraint. Walk back on constraints till we find the cusp
    // point.
    while (!is_xx_left(w) && ccL(cwL(w)) == w) w = cwL(w);
    return (!is_xx_left(w)) ? cw_target_edge(w) : ccw_source_edge(w);
    // ------------------------------------------------------------------------- 
}

// -----------------------------------------------------------------------------

template < class _A >
template < class _Tr >
typename Chi2_flip_traits<_A>::Vertex_handle
Chi2_flip_traits<_A>::walk(Edge_handle r , Edge_handle l, _Tr) const 
{
    // ------------------------------------------------------------------------- 
    // Operators used by this method
    typename _Tr::Less_vertex_handle chi2; 
    typename _Tr::Ccw_target_edge ccw_target_edge;
    typename _Tr::Ccw_source_edge ccw_source_edge;
    typename _Tr::Target_object   target_object;
    typename _Tr::Sup sup; 
    typename _Tr::Splice splice; 
    // ------------------------------------------------------------------------- 
    typename _Tr::Right_traits right_traits;
    typename _Tr::Left_traits  left_traits;
    // ------------------------------------------------------------------------- 
    Vertex_handle phiv = bit(r,l, _Tr());
    bool advance_r = advance(r,phiv, right_traits);
    bool advance_l = advance(l,phiv, left_traits);
    while (advance_r || advance_l) {
	// If we have to advance on both sides we use the chi2 predicate to
	// choose.
	if (advance_r && advance_l) advance_r = chi2(sup(r),sup(l));
	if (advance_r) {
	    if (target_object(sup(r)) == r->object())
		 r = ccw_source_edge(sup(r));
	    else r = ccw_target_edge(sup(r));
	}
	else {
	    if (target_object(sup(l)) == l->object())
		 l = ccw_source_edge(sup(l));
	    else l = ccw_target_edge(sup(l));
	}
	phiv = bit(r,l, _Tr());
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

template < class _A >
class PhiR_flip_traits : public Chi2_flip_traits<_A>
{
public:
    // -------------------------------------------------------------------------
    typedef _A Antichain;
    typedef typename Antichain::Vertex_handle Vertex_handle;
    typedef typename Antichain::Edge_handle   Edge_handle;
    typedef typename Antichain::Edge          Edge;
    typedef typename Antichain::Face_handle   Face_handle;
    typedef Chi2_flip_traits<_A>          Base;
    // -------------------------------------------------------------------------
    // CONSTRUCTOR
    PhiR_flip_traits(/*const*/ Antichain* an) : Base(an) { }
    // -------------------------------------------------------------------------
    template < class _Tr >
    Vertex_handle operator()(const Vertex_handle& v , _Tr) const
    {
	CGAL_precondition(!a->is_on_convex_hull(v));
	typename _Tr::Right_traits right;
	typename _Tr::Left_traits  left;
	typename _Tr::Ccw_source_edge  ccw_source_edge;
	typename _Tr::Ccw_target_edge  ccw_target_edge;
	return walk(ccw_source_edge(phiR(v,right)),
		    ccw_target_edge(phiR(v,left)),_Tr());
    }
    // -------------------------------------------------------------------------
protected:
    // -------------------------------------------------------------------------
    template < class _Tr >
    Vertex_handle phiR(const Vertex_handle& v , _Tr tr) const;
    template < class _Tr >
    Vertex_handle phiR_walk(Edge_handle r , Edge_handle l , 
			    Vertex_handle phis,_Tr tr) const;
    // ------------------------------------------------------------------------
    template < class _Tr >
    Vertex_handle bit_with_pi(Edge_handle r, Edge_handle l,_Tr tr)  const
    {
	typename _Tr::Vertex_creator vertex;
	CGAL_precondition(l != 0);
	Edge_handle lpi = new Edge(*l); 
	lpi->set_sign(!l->sign());
	Vertex_handle v = vertex(r,lpi); delete lpi;
	return v;
    }
    // -------------------------------------------------------------------------
};

// -----------------------------------------------------------------------------

template < class _A >
template < class _Tr >
inline typename PhiR_flip_traits<_A>::Vertex_handle
PhiR_flip_traits<_A>::
phiR_walk(Edge_handle r , Edge_handle l , Vertex_handle phis,  _Tr tr) const
{
    // ------------------------------------------------------------------------- 
    // Operators used by this method
    typename _Tr::Sup sup;
    typename _Tr::Less_vertex_handle chi1;
    typename _Tr::Ccw_target_edge ccw_target_edge;
    typename _Tr::Ccw_source_edge ccw_source_edge;
    typename _Tr::Set_ccw_edge    set_ccw_edge;
    typename _Tr::Right_traits right_traits;
    // ------------------------------------------------------------------------- 
    Vertex_handle phiR = bit_with_pi(r,l,_Tr());
    bool advance_r = advance(r,phiR, right_traits);
    bool advance_l = (!advance_r && sup(l) != phis && chi1(sup(l)->pi(),phiR));
#ifdef FLIP_STATISTICS
    if (!advance_r && sup(l) != phis) ++chi1_calls;
#endif
    while (advance_r || advance_l) {
	if (advance_r) r = ccw_source_edge(sup(r));
	else           l = ccw_target_edge(sup(l));
	delete phiR; phiR = bit_with_pi(r,l,_Tr());
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

template < class _A >
template < class _Tr >
inline typename PhiR_flip_traits<_A>::Vertex_handle
PhiR_flip_traits<_A>::phiR(const Vertex_handle& v , _Tr tr) const
{
    // ------------------------------------------------------------------------- 
    // Operators used by this method
    typename _Tr::Inf inf;
    typename _Tr::CcR ccR; 
    typename _Tr::Ccw_source_edge ccw_source_edge;
    // ------------------------------------------------------------------------- 
    Edge_handle r = first_cusp_edge(v,_Tr());
    Edge_handle l = ccw_source_edge(v->pi());
    Vertex_handle phis = inf(inf(ccR(v)));
    // ------------------------------------------------------------------------- 
    return phiR_walk(r,l,phis,tr);
    // ------------------------------------------------------------------------- 
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

template < class _A >
class Chi1_flip_traits : public PhiR_flip_traits<_A>
{
public:
    // -------------------------------------------------------------------------
    typedef _A Antichain;
    typedef typename Antichain::Vertex_handle Vertex_handle;
    typedef typename Antichain::Edge_handle   Edge_handle;
    typedef typename Antichain::Edge          Edge;
    typedef typename Antichain::Face_handle   Face_handle;
    typedef typename Antichain::Right_ccw_traits Right_ccw_traits;
    typedef typename Antichain::Left_ccw_traits  Left_ccw_traits;
    typedef typename Antichain::Right_cw_traits  Right_cw_traits;
    typedef typename Antichain::Left_cw_traits   Left_cw_traits;
    typedef PhiR_flip_traits<_A>          Base;
    // -------------------------------------------------------------------------
    // CONSTRUCTOR
    Chi1_flip_traits(/*const*/ Antichain* an) : Base(an) { }
    // -------------------------------------------------------------------------
    template < class _Tr >
    Vertex_handle operator()(const Vertex_handle& v , _Tr) const;
    // -------------------------------------------------------------------------
private:
    // -------------------------------------------------------------------------
    template < class _Tr >
    Vertex_handle chi1_walk(Edge_handle l, Edge_handle r, 
			    bool end_known) const;
    // -------------------------------------------------------------------------
};

// -----------------------------------------------------------------------------

template < class _A >
template < class _Tr >
inline typename Chi1_flip_traits<_A>::Vertex_handle
Chi1_flip_traits<_A>::
chi1_walk(Edge_handle l , Edge_handle r, bool end_known = true) const
{
    // -------------------------------------------------------------------------
    // Operators used in this method
    typename _Tr::Sign            sign;
    typename _Tr::Splice          splice;
    typename _Tr::Sup sup; typename _Tr::Inf inf;
    typename _Tr::Target_object target_object;
    typename _Tr::Less_vertex_handle chi1;
    typename _Tr::Cw_target_edge  cw_target_edge;
    typename _Tr::Cw_source_edge  cw_source_edge;
    typename _Tr::Ccw_target_edge ccw_target_edge;
    typename _Tr::Ccw_source_edge ccw_source_edge;
    typename _Tr::Bottom_edge     bottom_edge;
    typename _Tr::Right_traits right_traits;
    typename _Tr::Left_cw_traits left_cw_traits;
    typename _Tr::Right_cw_traits right_cw_traits;
    // -------------------------------------------------------------------------
    // If l is negative then l is the first atom of Lthr(v)
    if ( !sign(l) || l == cw_target_edge(sup(l)) ) return walk(r,l,_Tr());
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
    Vertex_handle phiR = bit_with_pi(r,l,_Tr());
    bool advance_r = advance(r,phiR,right_traits);
    bool advance_l = (!advance_r && l != lmax && chi1(sup(l)->pi(),phiR));
#ifdef FLIP_STATISTICS
    if (!advance_r && l != lmax) ++chi1_calls;
#endif
    while (advance_r || advance_l) {
	if (advance_r) r = ccw_source_edge(sup(r));
	else           l = ccw_target_edge(sup(l));
	if (l == ccw_target_edge(phis)) break;
	delete phiR; phiR = bit_with_pi(r,l,_Tr());
	advance_r = advance(r,phiR,right_traits);
	advance_l = (!advance_r && l != lmax && chi1(sup(l)->pi(),phiR));
#ifdef FLIP_STATISTICS
	if (!advance_r && l != lmax) ++chi1_calls;
#endif
    }
    // -------------------------------------------------------------------------
    // Special case when phiR(sup(l)) leaves first atom of Lthr(v)
    if (l == ccw_target_edge(phis)) {
	CGAL_precondition(end_known == false);
	l = lmin;
	while (sign(l)) l = ccw_target_edge(sup(l));
	sup(lmin)->set_phiR(l);
	return walk(r,l,_Tr()); // phi(v) is xx-right enters Lthr(v)
    }
    // -------------------------------------------------------------------------
    if (!sign(l)) { 
	CGAL_precondition(l->object() == lmin->object());
	Vertex_handle phiv = bit(r,lmin,_Tr());
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

template < class _A >
template < class _Tr >
inline typename Chi1_flip_traits<_A>::Vertex_handle
Chi1_flip_traits<_A>::operator()(const Vertex_handle& v , _Tr tr) const
{
    // -------------------------------------------------------------------------
    // Operators used in this method
    typename _Tr::Is_xx_left      is_xx_left;
    typename _Tr::Is_left_xx      is_left_xx;
    typename _Tr::Sup sup; typename _Tr::Inf inf;
    typename _Tr::Less_vertex_handle chi1;
    typename _Tr::Vertex_creator  vertex;
    typename _Tr::Splice          splice;
    typename _Tr::Ccw_target_edge ccw_target_edge;
    typename _Tr::Ccw_source_edge ccw_source_edge;
    typename _Tr::Cw_target_edge  cw_target_edge;
    typename _Tr::Cw_source_edge  cw_source_edge;
    typename _Tr::CcR ccR; typename _Tr::CcL ccL;
    typename _Tr::CwR cwR; typename _Tr::CwL cwL;
    typename _Tr::Right_traits right_traits;
    typename _Tr::Left_traits  left_traits;
    typename _Tr::Right_cw_traits  right_cw_traits;
    typename _Tr::Left_cw_traits   left_cw_traits;
    // -------------------------------------------------------------------------
    // Starting point on Rtwo(v)
    Edge_handle r = 0;
    if (is_xx_left(v) && cwR(ccR(v)) == v) {
	if (ccR(v)->phiR() == 0) {
	    ccR(v)->set_phiR_vertex(phiR(v,right_traits));
	    ccR(v)->set_phiR(ccw_source_edge(ccR(v)->phiR_vertex()));
	}
	r = ccR(v)->phiR();
    }
    else r = first_cusp_edge(v,right_traits);
    CGAL_precondition(r != 0);
    // -------------------------------------------------------------------------
    // Starting point on Rtwo(v)
    Edge_handle l = 0;
    if (!is_left_xx(v) && cwL(ccL(v)) == v) { 
	if (ccL(v)->phiL() == 0) {
	    ccL(v)->set_phiL_vertex(phiR(v,left_traits));
	    ccL(v)->set_phiL(ccw_target_edge(ccL(v)->phiL_vertex()));
	}
	l = ccL(v)->phiL();
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
	    if (!is_left_xx(v) && ccR(cwR(v)) == v)
		 l = cw_source_edge(cwR(v));
	    else if (ccR(cwR(v)) == v)
		 l = third_cusp_edge(cwR(v),right_cw_traits);
	    else l = second_cusp_edge(cwR(v),left_cw_traits);
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
	    l = ccw_target_edge(sup(l));
	}
    }
    // -------------------------------------------------------------------------
    // If phiv was not found at the previuos step walk by computing successive
    // phiR.
    if (phiv == 0) {
	Edge_handle lmin = l;
	phiv = chi1_walk<_Tr>(l,r,false);
	while ( phiv == 0 ) {
	    l = ccw_target_edge(sup(l));
	    phiv = chi1_walk<_Tr>(l,r,false);
#ifdef FLIP_STATISTICS
	    ++psi_computed;
#endif
	}
	while (lmin != l) {
	    if (sup(lmin)->phiR() == cw_source_edge(phiv))
		sup(lmin)->set_phiR(ccw_source_edge(phiv));
	    lmin = ccw_target_edge(sup(lmin));
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

template < class _A >
class Backward_flip_traits : public Chi2_flip_traits<_A>
{
public:
    // -------------------------------------------------------------------------
    typedef _A Antichain;
    typedef typename Antichain::Vertex_handle Vertex_handle;
    typedef typename Antichain::Edge_handle   Edge_handle;
    typedef typename Antichain::Edge          Edge;
    typedef typename Antichain::Face_handle   Face_handle;
    typedef Chi2_flip_traits<_A>          Base;
    // -------------------------------------------------------------------------
    // CONSTRUCTOR
    Backward_flip_traits(/*const*/ Antichain* an) : Base(an) { }
    // -------------------------------------------------------------------------
    template < class _Tr >
    Vertex_handle operator()(const Vertex_handle& v , _Tr) const
    {
	CGAL_precondition(!a->is_on_convex_hull(v));
	typename _Tr::Right_traits right;
	typename _Tr::Left_traits  left;
	return walk_backward<_Tr>(second_cusp_edge(v,right),
				  second_cusp_edge(v,left));
			    
    }
    // -------------------------------------------------------------------------
private:
    template < class _Tr >
    Vertex_handle walk_backward(const Edge_handle& right,
			        const Edge_handle& left) const;
    template < class _Tr >
    bool advance(Edge_handle r, Vertex_handle phiv) const;
};

// -----------------------------------------------------------------------------
// Returns false if phiv is tangent to r

template < class _A >
template < class _Tr >
inline bool
Backward_flip_traits<_A>::advance(Edge_handle r, Vertex_handle phiv) const
{
    typename _Tr::Inf inf; 
    typename _Tr::Less_vertex_handle chi1; 
    typename _Tr::Target_object target_object; 
    typename _Tr::Is_xx_left is_xx_left;
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

template < class _A >
template < class _Tr >
typename Backward_flip_traits<_A>::Vertex_handle
Backward_flip_traits<_A>::walk_backward(const Edge_handle& right , 
					const Edge_handle& left) const 
{
    // ------------------------------------------------------------------------- 
    // Operators used by this method
    typename _Tr::Cw_target_edge cw_target_edge;
    typename _Tr::Cw_source_edge cw_source_edge;
    typename _Tr::Inf inf; 
    typename _Tr::Splice splice; 
    // ------------------------------------------------------------------------- 
    typename _Tr::Right_traits right_traits;
    typename _Tr::Left_traits  left_traits;
    // ------------------------------------------------------------------------- 
    Edge_handle r = right; Edge_handle l = left;
    Vertex_handle phiv = bit(r,l,_Tr());
    bool advance_r = advance(r,phiv,right_traits);
    bool advance_l = (!advance_r && advance(l,phiv,left_traits));
    while (advance_r || advance_l) {
	if (advance_r) r = cw_target_edge(inf(r));
	else           l = cw_source_edge(inf(l));
	delete phiv; phiv = bit(r,l,_Tr());
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

struct Visibility_complex_flip_traits {
    template < class _A >
    struct Flip_wrapper {
	typedef Chi2_flip_traits<_A> Flip_traits;
    };
};

// -----------------------------------------------------------------------------

struct Visibility_complex_phiR_flip_traits {
    template < class _A >
    struct Flip_wrapper {
	typedef PhiR_flip_traits<_A> Flip_traits;
    };
};

// -----------------------------------------------------------------------------

struct Visibility_complex_chi1_flip_traits {
    template < class _A >
    struct Flip_wrapper {
	typedef Chi1_flip_traits<_A> Flip_traits;
    };
};

// -----------------------------------------------------------------------------

struct Visibility_complex_backward_flip_traits {
    template < class _A >
    struct Flip_wrapper {
	typedef Backward_flip_traits<_A> Flip_traits;
    };
};

// -----------------------------------------------------------------------------

CGAL_END_NAMESPACE

#endif // VISIBILITY_COMPLEX_FLIP_TRAITS_H
