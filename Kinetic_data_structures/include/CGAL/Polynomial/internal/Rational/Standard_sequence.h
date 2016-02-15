// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_STANDARD_SEQUENCE_H
#define CGAL_STANDARD_SEQUENCE_H


namespace CGAL { namespace POLYNOMIAL { namespace internal {

template<class Sturm_sequence_t>
class Standard_sequence : public Sturm_sequence_t
{
    public:
        typedef Sturm_sequence_t                        Sturm_sequence;
        typedef typename Sturm_sequence::Kernel         Kernel;
        typedef typename Sturm_sequence::Polynomial     Polynomial;

        typedef typename Polynomial::NT                 NT;

    protected:
        typedef Sturm_sequence                          Base;

    public:
        Standard_sequence() : Base() {}
        Standard_sequence(const Polynomial& p, const Kernel& k=Kernel())
            : Base(p, k.differentiate_object()(p), k) {}

    protected:
        template<class NTRep>
            unsigned int number_of_real_roots_base(const NTRep& a,
            const NTRep& b) const
        {
            CGAL_precondition( b >= a );

            unsigned int Va = sign_variations(a);
            if ( Va == 0 ) { return 0; }

            unsigned int Vb = sign_variations(b);

//CGAL_assertion( Va > Vb );

            return Va - Vb;
        }

    public:
        template<class T>
            unsigned int
            number_of_real_roots(const T& a, const T& b) const
        {
            return number_of_real_roots_base(a, b);
        }

};

#if 0
template<class Kernel, class Ret, class M = CGAL::Integral_domain_without_division_tag>
class Standard_sequence_k
: public Sturm_sequence_k<Kernel, Ret,M>
{
//public:
    typedef Sturm_sequence_k<Kernel, Ret,M> Parent;
    typedef typename Parent::Polynomial                                       Polynomial;
//typedef M                                       Method_tag;

    typedef typename Parent::NT                 NT;

//protected:

    public:
        Standard_sequence_k() : Parent() {}
        Standard_sequence_k(const Polynomial& p)
            : Parent(p, p.derivative()) {}
};
#endif

} } } //namespace CGAL::POLYNOMIAL::internal
#endif                                            // CGAL_STANDARD_SEQUENCE_H
