#ifndef VISIBILITY_COMPLEX_LINEAR_DECORATOR_H
#define VISIBILITY_COMPLEX_LINEAR_DECORATOR_H

#ifndef VISIBILITY_COMPLEX_DECORATOR_CCW_CW_TRAITS_H
#include<CGAL/Visibility_complex_ccw_cw_traits.h>
#endif

CGAL_BEGIN_NAMESPACE

// -----------------------------------------------------------------------------
// This decorator assumes linear space. In ohter words the only knonw vertices
// are the ones given from a pseudo-triangulation.

template <class _Vc>
class Visibility_complex_linear_decorator
{
public:
    // -------------------------------------------------------------------------
    typedef _Vc                                Antichain;
    typedef typename _Vc::Gt                   Gt;
    // -------------------------------------------------------------------------
    typedef typename Antichain::Vertex         Vertex;
    typedef typename Antichain::Vertex_handle  Vertex_handle;
    typedef typename Antichain::Edge           Edge;  
    typedef typename Antichain::Edge_handle    Edge_handle;
    typedef typename Antichain::Face           Face;
    typedef typename Antichain::Face_handle    Face_handle;
    // -------------------------------------------------------------------------
    typedef Visibility_complex_left_ccw_traits<Antichain>  Left_ccw_traits;
    typedef Visibility_complex_right_ccw_traits<Antichain> Right_ccw_traits;
    typedef Visibility_complex_left_cw_traits<Antichain>   Left_cw_traits;
    typedef Visibility_complex_right_cw_traits<Antichain>  Right_cw_traits;
    // -------------------------------------------------------------------------
private:
    // -------------------------------------------------------------------------
    // Internal reference to visibility complex
    Antichain*  a;
    // -------------------------------------------------------------------------
public:
    // -------------------------------------------------------------------------
    // Constructor
    Visibility_complex_linear_decorator( Antichain& v) : a(&v) {}
    // -------------------------------------------------------------------------
    // pi operator
    Vertex_handle pi(Vertex_handle b) { return b->pi(); }
    // -------------------------------------------------------------------------
    // phi operator
    template <class _Tr >
    Vertex_handle phi(Vertex_handle b,_Tr) {
	typename _Tr::Sup sup;
	if (b.is_constraint()) return b->pi();
	return sup(sup(b));
    }
    Vertex_handle phi(Vertex_handle b) { return b->sup()->sup(); }
    // -------------------------------------------------------------------------
    // phi_star operator
    template <class _Tr >
    Vertex_handle phi_star(Vertex_handle b,_Tr) {
	typename _Tr::Inf inf;
	return inf(inf(b));
    }
    Vertex_handle phi_star(Vertex_handle b) { return b->inf()->inf(); }
    // -------------------------------------------------------------------------
    // ccw and cw operators
    Vertex_handle ccR(Vertex_handle b) { return b->ccR(); }
    Vertex_handle ccL(Vertex_handle b) { return b->ccL(); }
    Vertex_handle cwR(Vertex_handle b) { return b->cwR(); }
    Vertex_handle cwL(Vertex_handle b) { return b->cwL(); }
    // -------------------------------------------------------------------------
    // Face Operators
    template <class _Tr > 
    Vertex_handle phirr(Vertex_handle b,_Tr);
    Vertex_handle phirr(Vertex_handle b) 
	{ return phirr(b,Right_ccw_traits()); }
    template <class _Tr > 
    Vertex_handle phill(Vertex_handle b,_Tr)
	{ return phirr(b,change_right_left(_Tr())); }
    Vertex_handle phill(Vertex_handle b) 
	{ return phill(b,Right_ccw_traits()); }
    // -------------------------------------------------------------------------
    // Face Operators
    template <class _Tr > 
    Vertex_handle phiback(Vertex_handle b,_Tr);
    Vertex_handle phiback(Vertex_handle b) 
	{ return phiback(b,Right_ccw_traits()); }
    // -------------------------------------------------------------------------
    template <class _Tr > 
    Vertex_handle phiforw(Vertex_handle b,_Tr) 
	{ return phiback(b,change_right_left(_Tr())); }
    Vertex_handle phiforw(Vertex_handle b) 
	{ return phiforw(b,Right_ccw_traits()); }
    // -------------------------------------------------------------------------
    template <class _Tr > 
    Vertex_handle phiback_star(Vertex_handle b,_Tr) 
	{ return phiback(b,change_ccw_cw(_Tr())); }
    Vertex_handle phiback_star(Vertex_handle b) 
	{ return phiback_star(b,Right_ccw_traits()); }
    // -------------------------------------------------------------------------
    template <class _Tr > 
    Vertex_handle phiforw_star(Vertex_handle b,_Tr) 
	{ return phiforw(b,change_ccw_cw(_Tr())); }
    Vertex_handle phiforw_star(Vertex_handle b) 
	{ return phiforw_star(b,Right_ccw_traits()); }
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // Face boundaries
    // FIXME
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // Greedy Pseudo-triangles and subsets of greedy pseudotriangles
    template < class _OutputIt , class _Tr >
    _OutputIt R1(Vertex_handle b, _OutputIt result , _Tr tr);
    template < class _OutputIt >
    _OutputIt R1(Vertex_handle b, _OutputIt result)
	{ return R1(b,result,Right_ccw_traits()); }
    template < class _OutputIt , class _Tr >
    _OutputIt R2_1(Vertex_handle b, _OutputIt result, _Tr tr);
    template < class _OutputIt >
    _OutputIt R2_1(Vertex_handle b, _OutputIt result)
	{ return R2_1(b,result,Right_ccw_traits()); }
    template < class _OutputIt , class _Tr >
    _OutputIt R2_2(Vertex_handle b, _OutputIt result, _Tr tr);
    template < class _OutputIt , class _Tr >
    _OutputIt R2(Vertex_handle b, _OutputIt result , _Tr tr);
    template < class _OutputIt >
    _OutputIt R2(Vertex_handle b, _OutputIt result)
	{ return R2(b,result,Right_ccw_traits()); }
    template < class _OutputIt , class _Tr >
    _OutputIt R3(Vertex_handle b, _OutputIt result, _Tr tr);
    template < class _OutputIt , class _Tr >
    _OutputIt RP(Vertex_handle b, _OutputIt result , _Tr tr);
    // -------------------------------------------------------------------------
    template < class _OutputIt , class _Tr >
    _OutputIt R(Vertex_handle b, _OutputIt result , _Tr tr);
    template < class _OutputIt >
    _OutputIt R(Vertex_handle b, _OutputIt result)
	{ return R(b,result,Right_ccw_traits()); }
    // -------------------------------------------------------------------------
    template < class _OutputIt , class _Tr >
    _OutputIt L(Vertex_handle b, _OutputIt result , _Tr) 
	{ return R(b,result,change_right_left(_Tr())); }
    template < class _OutputIt >
    _OutputIt L(Vertex_handle b, _OutputIt result) 
	{ return L(b,result,Right_ccw_traits()); }
    // -------------------------------------------------------------------------
    template < class _OutputIt , class _Tr >
    _OutputIt L_star(Vertex_handle b, _OutputIt result , _Tr) 
	{ return L(b,result,change_ccw_cw(_Tr())); }
    template < class _OutputIt >
    _OutputIt L_star(Vertex_handle b, _OutputIt result) 
	{ return L_star(b,result,Right_ccw_traits()); }
    // -------------------------------------------------------------------------
    template < class _OutputIt , class _Tr >
    _OutputIt R_star(Vertex_handle b, _OutputIt result , _Tr) 
	{ return R(b,result,change_ccw_cw(_Tr())); }
    template < class _OutputIt >
    _OutputIt R_star(Vertex_handle b, _OutputIt result) 
	{ return R_star(b,result,Right_ccw_traits()); }
    // -------------------------------------------------------------------------
    // phiR and phiL Operators
    // FIXME
    // -------------------------------------------------------------------------
};

// -----------------------------------------------------------------------------

template < class _Vc>
template < class _Tr > 
Visibility_complex_linear_decorator<_Vc>::Vertex_handle
Visibility_complex_linear_decorator<_Vc>::phirr(Vertex_handle b,_Tr) 
{
    typename _Tr::Sup              sup;
    typename _Tr::Dr               dr;
    typename _Tr::Cw_source_edge   cw_source_edge;

    return sup(dr(cw_source_edge(b)));
}

// -----------------------------------------------------------------------------

template < class _Vc>
template < class _Tr > 
Visibility_complex_linear_decorator<_Vc>::Vertex_handle
Visibility_complex_linear_decorator<_Vc>::phiback(Vertex_handle b,_Tr t) 
{
    typename _Tr::Is_left_xx       is_left_xx;
    typename _Tr::Sup              sup;
    typename _Tr::Ul               ul;
    typename _Tr::Cw_source_edge   cw_source_edge;
    typename _Tr::Source_object    source_object;

    if (is_left_xx(b)) return sup(ul(cw_source_edge(b)));

    Vertex_handle c = b;
    while (source_object(c) == source_object(b)) c = ccL(b,t);
    return c;
}

// -----------------------------------------------------------------------------

template < class _Vc>
template < class _OutputIt , class _Tr >
_OutputIt
Visibility_complex_linear_decorator<_Vc>::
R1(Vertex_handle b, _OutputIt result ,_Tr tr)
{
    typename _Tr::Is_xx_left       is_xx_left;
    typename _Tr::CcR ccR; typename _Tr::CwR cwR;

    CGAL_precondition(b != 0);
    *result++ = *b;
    if (is_xx_left(b) && a->is_on_convex_hull(b)) return result;

    Vertex_handle c = b;
    while (is_xx_left(c) && cwR(ccR(c)) == c) {
	c = ccR(c);
	*result++ = *c;
    }
    return result;
}

// -----------------------------------------------------------------------------

template < class _Vc>
template < class _OutputIt , class _Tr >
_OutputIt
Visibility_complex_linear_decorator<_Vc>::
R2_1(Vertex_handle b, _OutputIt result, _Tr tr)
{
    typename _Tr::Sup              sup;
    typename _Tr::Dr               dr;
    typename _Tr::Cw_source_edge   cw_source_edge;
    typename _Tr::Is_left_xx       is_left_xx;

    CGAL_precondition(b != 0);
    if (a->is_on_convex_hull(b) && is_left_xx(b)) return result;

    Vertex_handle phib = phi(b,tr);
    Vertex_handle phiback_phib = phiback(phib,tr);
    Vertex_handle c = sup(dr(cw_source_edge(b))); // phirr(b) 
    while (!is_left_xx(c) && c != phiback_phib) {
	*result++ = *c;
	c = phiback(c,tr);
    }
    if (is_left_xx(c) && c != phiback_phib) *result++ = *c;
    return result;
}

// -------------------------------------------------------------------------

template < class _Vc>
template < class _OutputIt , class _Tr >
_OutputIt
Visibility_complex_linear_decorator<_Vc>::
R2_2(Vertex_handle b, _OutputIt result, _Tr tr)
{
    typename _Tr::Sup              sup;
    typename _Tr::Is_left_xx      is_left_xx;

    CGAL_precondition(b != 0);
    if (a->is_on_convex_hull(b) && is_left_xx(b)) return result;

    Vertex_handle phib = sup(sup(b));
    if (is_left_xx(phib)) return result;

    Vertex_handle phibackb = phiback(b,tr);
    Vertex_handle phiback_phibackb = phiback(phibackb,tr);
    Vertex_handle c = phiback(phib,tr);
    while (!is_left_xx(c) && c != phiback_phibackb) {
	*result++ = *c;
	c = phiback(c,tr);
    }
    if (is_left_xx(c) && c != phiback_phibackb) *result++ = *c;
    return result;
}

// -------------------------------------------------------------------------

template < class _Vc>
template < class _OutputIt , class _Tr >
_OutputIt
Visibility_complex_linear_decorator<_Vc>::
R3(Vertex_handle b, _OutputIt result, _Tr tr)
{
    typename _Tr::Is_left_xx      is_left_xx;
    typename _Tr::Is_xx_left      is_xx_left;

    CGAL_precondition(b != 0);
    if (!is_left_xx(b) || a->is_on_convex_hull(b)) return result;

    Vertex_handle c = phiback(b,tr); // min R^3
    while (is_xx_left(c)) {
	*result++ = *c;
	c = phiforw(c,tr);
    }
    *result++ = *c;
    return result;
}

// -------------------------------------------------------------------------

template < class _Vc>
template < class _OutputIt , class _Tr >
_OutputIt
Visibility_complex_linear_decorator<_Vc>::
RP(Vertex_handle b, _OutputIt result, _Tr tr)
{
    typename _Tr::Sup              sup;
    typename _Tr::Dr dr; typename _Tr::Ul ul;
    typename _Tr::Cw_source_edge   cw_source_edge;
    typename _Tr::Cw_target_edge   cw_target_edge;
    typename _Tr::Is_left_xx       is_left_xx;

    CGAL_precondition(b != 0);
    if (a->is_on_convex_hull(b) && is_left_xx(b)) return result;

    Vertex_handle phiback_phib = phiback(sup(sup(b)),tr);
    Vertex_handle phirr = sup(dr(cw_source_edge(b))); // phirr(b) 
    Vertex_handle c = phirr;
    while (!is_left_xx(c) && c != phiback_phib) {
	if (c != phirr && ul(cw_target_edge(c)) != sup(b)) *result++ = *c;
	c = phiback(c,tr);
    }
    if (is_left_xx(c) && c != phiback_phib) *result++ = *c;
    return result;
}

// -------------------------------------------------------------------------

template < class _Vc>
template < class _OutputIt , class _Tr >
_OutputIt
Visibility_complex_linear_decorator<_Vc>::R(Vertex_handle b, 
				     _OutputIt result ,
				     _Tr tr) 
{
    R1  (b,result,tr);
    R2_1(b,result,tr);
    R2_2(b,result,tr);
    R3  (b,result,tr);
    return result;
}

// -----------------------------------------------------------------------------

CGAL_END_NAMESPACE

#endif // VISIBILITY_COMPLEX_DECORATOR_H
