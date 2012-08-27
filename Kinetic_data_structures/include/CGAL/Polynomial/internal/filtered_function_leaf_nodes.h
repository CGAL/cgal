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

#ifndef CGAL_POLYNOMIAL_INTERNAL_VIRTUAL_FUNCTION_REPS_H
#define CGAL_POLYNOMIAL_INTERNAL_VIRTUAL_FUNCTION_REPS_H
#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/internal/filtered_function_node_bases.h>
#include <CGAL/Polynomial/internal/interval_arithmetic.h>

namespace CGAL { namespace POLYNOMIAL { namespace internal {

template <class Traits>
class Filtered_function_node_constant: public Filtered_function_node<Traits>
{
    typedef Filtered_function_node_constant<Traits> This;
    typedef Filtered_function_node<Traits> P;
    public:
        Filtered_function_node_constant(const typename P::Exact_function::NT& val,
        const typename P::Interval_function_converter &ifc): P(ifc) {
            P::set_interval_function(typename P::Interval_function(ifc.nt_converter()(val)));
            P::set_exact_function(typename P::Exact_function(val/*this->num_*/));
        }

        virtual ~Filtered_function_node_constant(){}

        virtual void write(std::ostream &out) const
        {
            out << P::exact_function();
        }
    protected:

        virtual  void generate_exact_function() const
        {
            CGAL_Polynomial_assertion(0);
        }
};

template <class Traits>
class Filtered_function_node_double_constant: public Filtered_function_node<Traits>
{
    typedef Filtered_function_node_double_constant<Traits> This;
    typedef Filtered_function_node<Traits> P;
    public:
        Filtered_function_node_double_constant(double d, const typename P::Interval_function_converter &ifc): P(ifc) {
            P::set_interval_function(typename P::Interval_function(std::pair<double,double>(d,d)));

        }

        virtual ~Filtered_function_node_double_constant(){}

        virtual void write(std::ostream &out) const
        {
            out << P::interval_function()[0].sup();
        }
    protected:

        virtual  void generate_exact_function() const
        {
            P::set_exact_function(typename P::Exact_function(typename P::Exact_function::NT(P::interval_function()[0].sup())));
        }
};

template <class Traits>
class Filtered_function_node_explicit: public Filtered_function_node<Traits>
{
    typedef Filtered_function_node_explicit<Traits> This;
    typedef Filtered_function_node<Traits> P;
    public:
        Filtered_function_node_explicit(const typename P::Exact_function f, const typename P::Interval_function_converter &ifc): P(ifc) {
            P::set_interval_function(ifc(f));
            P::set_exact_function(f);
        }

        virtual ~Filtered_function_node_explicit(){}

        virtual void write(std::ostream &out) const
        {
            if (P::exact_function().degree() >0) {
                out << "("<< P::exact_function() << ")";
            }
            else {
                out << P::exact_function();
            }
        }
    protected:

        virtual  void generate_exact_function() const
        {
            CGAL_Polynomial_assertion(0);
        }
};

template <class Traits>
class Filtered_function_node_double_explicit: public Filtered_function_node<Traits>
{
    typedef Filtered_function_node_double_explicit<Traits> This;
    typedef Filtered_function_node<Traits> P;
    public:
        template <class DF>
        Filtered_function_node_double_explicit(const DF& f, const typename P::Interval_function_converter &ifc): P(ifc) {
            typename P::Interval_function fi;
            fi.set_nominal_degree(f.degree());
            for (int i=0; i<= f.degree(); ++i) {
                fi.set_coef(i, std::pair<double,double>(f[i], f[i]));
            }
            P::set_interval_function(fi);
        }
        virtual ~Filtered_function_node_double_explicit(){}

        virtual void write(std::ostream &out) const
        {
            if (P::interval_function().nominal_degree() >0) {
                out << "("<< P::exact_function() << ")";
            }
            else {
                out << P::exact_function();
            }
        }
    protected:

        virtual  void generate_exact_function() const
        {
            typename P::Exact_function fe;
            int deg= P::interval_function().nominal_degree();
            fe.set_nominal_degree(deg);
            for (int i=0; i<= this->f.degree(); ++i) {
                fe.set_coef(i, typename P::Exact_function::NT(this->f[i].sup()));
            }
            P::set_exact_function(fe);
        }
};

} } } //namespace CGAL::POLYNOMIAL::internal
#endif
