// ============================================================================
//
// Copyright (c) 1999,2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Interval_arithmetic/IA_Lazy_exact_nt.h
// revision      : $Revision$
// revision_date : $Date$
// package       : Interval Arithmetic
// author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
// ============================================================================

#ifndef CGAL_IA_LAZY_EXACT_NT_H
#define CGAL_IA_LAZY_EXACT_NT_H

CGAL_BEGIN_NAMESPACE

template <class RT>
inline
Interval_nt_advanced
convert_from_to (const Interval_nt_advanced&, const Lazy_exact_nt<RT> & z)
{
	CGAL_expensive_assertion(FPU_empiric_test() == CGAL_FE_UPWARD);
	return z.approx_adv();
}

#ifndef CGAL_CFG_NO_PARTIAL_CLASS_TEMPLATE_SPECIALISATION
template <class RT>
struct converter<Interval_nt_advanced, Lazy_exact_nt<RT> >
{
    static inline Interval_nt_advanced do_it (const Lazy_exact_nt<RT> & z)
    {
	return convert_from_to(Interval_nt_advanced(), z);
    }
};
#endif // CGAL_CFG_NO_PARTIAL_CLASS_TEMPLATE_SPECIALISATION

CGAL_END_NAMESPACE

#endif // CGAL_IA_LAZY_EXACT_NT_H
