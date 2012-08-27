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

#ifndef CGAL_POLYNOMIAL_INTERNAL_FILTERED_RATIONAL_TRAITS_H
#define CGAL_POLYNOMIAL_INTERNAL_FILTERED_RATIONAL_TRAITS_H
#include <CGAL/Polynomial/basic.h>

#include <CGAL/Polynomial/internal/Filtered_function.h>
#include <CGAL/Polynomial/internal/Rational/Rational_traits_base.h>

#include <CGAL/Polynomial/internal/Filtered_rational/Construct_filtered_function.h>
#include <CGAL/Polynomial/internal/Filtered_rational/Filtered_are_negations.h>
#include <CGAL/Polynomial/internal/Filtered_rational/Filtered_sign_at_rational.h>
#include <CGAL/Polynomial/internal/Filtered_rational/Filtered_root_bound_evaluator.h>
#include <CGAL/Polynomial/internal/Filtered_rational/Filtered_rational_multiplicity.h>
#include <CGAL/Polynomial/internal/Rational/Sign_above_rational.h>
#include <CGAL/Polynomial/internal/Rational/Sign_below_rational.h>

#define CGAL_DF_UNARY_CONSTRUCTION(UCName, lcname) class UCName {\
    typedef typename Exact_traits::UCName ED;\
    typedef typename Interval_traits::UCName ID;\
    typedef Filtered_function_node_unary_transform< Filtering_traits, ED, ID> Node;\
    public:\
        UCName(const This &k): ed_(k.exact_traits_object().lcname##_object()),\
        id_(k.interval_traits_object().lcname##_object()){} \
        typedef Function result_type;\
        typedef result_type argument_type;\
        result_type operator()(const argument_type &f) const \
        { \
            return result_type(new Node(f.tree(), ed_, id_));\
        }\
        protected:\
            ED ed_;\
            ID id_;\
        };\
        UCName lcname##_object() const \
        { \
            return UCName(*this);\
        }

// std::cout << "UCName of " << f << std::endl;

#define CGAL_DF_UNARY_CONSTRUCTION_DATA(UCName, lcname, data_type) class UCName {\
    typedef typename Exact_traits::UCName ED;\
    typedef typename Interval_traits::UCName ID;\
    typedef Filtered_function_node_unary_transform< Filtering_traits, ED, ID> Node;\
    public:\
        UCName(const typename ED::NT& d, const This &k): ed_(k.exact_traits_object().lcname##_object(d)),\
        id_(k.interval_traits_object().lcname##_object(CGAL_POLYNOMIAL_NS::To_interval<typename ED::NT>()(d))){} \
        UCName(double d, const This &k): ed_(k.exact_traits_object().lcname##_object(typename ED::NT(d))),\
        id_(k.interval_traits_object().lcname##_object(CGAL_POLYNOMIAL_NS::To_interval<double>()(d))){} \
        UCName(int d, const This &k): ed_(k.exact_traits_object().lcname##_object(d)),\
        id_(k.interval_traits_object().lcname##_object(d)){} \
        typedef Function result_type;\
        typedef result_type argument_type;\
        result_type operator()(const argument_type &f) const \
        { \
            return result_type(new Node(f.tree(), ed_, id_));\
        }\
        protected:\
            ED ed_;\
            ID id_;\
        };\
        UCName lcname##_object(data_type d) const \
        { \
            return UCName(d, *this);\
        }

//std::cout << "UCName of " << f << std::endl;

#define CGAL_DF_UNARY_CONSTRUCTION_DATA_2(UCName, lcname, data_type) class UCName {\
    typedef typename Exact_traits::UCName ED;\
    typedef typename Interval_traits::UCName ID;\
    typedef Filtered_function_node_unary_transform< Filtering_traits, ED, ID> Node;\
    public:\
        UCName(const typename ED::NT& a, const typename ED::NT& b, const This &k): ed_(k.exact_traits_object().lcname##_object(a,b)),\
        id_(k.interval_traits_object().lcname##_object(CGAL_POLYNOMIAL_NS::To_interval<typename ED::NT>()(a,b))){} \
        typedef Function result_type;\
        typedef result_type argument_type;\
        result_type operator()(const argument_type &f) const \
        { \
            return result_type(new Node(f.tree(), ed_, id_));\
        }\
        protected:\
            ED ed_;\
            ID id_;\
        };\
        UCName lcname##_object(data_type a, data_type b) const \
        { \
            return UCName(a, b, *this);\
        }

#define CGAL_DF_BINARY_CONSTRUCTION(UCName, lcname) class UCName {\
    typedef typename Exact_traits::UCName ED;\
    typedef typename Interval_traits::UCName ID;\
    typedef Filtered_function_node_binary_transform< Filtering_traits, ED, ID> Node;\
    public:\
        UCName(const This &k): ed_(k.exact_traits_object().lcname##_object()),\
        id_(k.interval_traits_object().lcname##_object()){} \
        typedef Function result_type;\
        typedef result_type first_argument_type;\
        typedef result_type second_argument_type;\
        result_type operator()(const first_argument_type &f, const second_argument_type &fp) const \
        { \
            return result_type(new Node(f.tree(), fp.tree(), ed_, id_));\
        }\
        protected:\
            ED ed_;\
            ID id_;\
        };\
        UCName lcname##_object() const \
        { \
            return UCName(*this);\
        }

//std::cout << "UCName of " << f << " and " << fp << std::endl;

                                namespace CGAL { namespace POLYNOMIAL { namespace internal {

                                template <class Filter_traits_t, template<class Fn> class Rational_traits = internal::Rational_traits_base>
                                class Filtered_rational_traits
                                {
                                    typedef Filtered_rational_traits<Filter_traits_t> This;
                                    public:
                                        typedef Filter_traits_t Filtering_traits;
                                        typedef Filtered_function<Filtering_traits> Function;
                                        typedef typename Filtering_traits::Exact_to_interval_converter Exact_to_interval_converter;
                                        typedef typename Function::NT NT;

                                        typedef Rational_traits<typename Filtering_traits::Exact_function> Exact_traits;
                                        typedef Rational_traits<typename Filtering_traits::Interval_function> Interval_traits;

                                        Filtered_rational_traits(){}

                                        CGAL_DF_UNARY_CONSTRUCTION(Differentiate, differentiate);

                                        CGAL_DF_UNARY_CONSTRUCTION(Invert_variable, invert_variable);

                                        CGAL_DF_UNARY_CONSTRUCTION_DATA(Rational_translate_zero, rational_translate_zero, NT);

                                        CGAL_DF_UNARY_CONSTRUCTION_DATA(Shift_power, shift_power, int);

                                        CGAL_DF_UNARY_CONSTRUCTION(Negate_variable, negate_variable);

                                        CGAL_DF_BINARY_CONSTRUCTION(Quotient, quotient);

                                        CGAL_DF_BINARY_CONSTRUCTION(Remainder, remainder);

                                        CGAL_DF_BINARY_CONSTRUCTION(Pseudo_quotient, pseudo_quotient);

                                        CGAL_DF_BINARY_CONSTRUCTION(Pseudo_remainder, pseudo_remainder);

                                        CGAL_DF_BINARY_CONSTRUCTION(Quotient_remainder, quotient_remainder);

                                        CGAL_DF_UNARY_CONSTRUCTION_DATA_2(Map_rational_interval_to_positive_2, map_rational_interval_to_positive_2, NT);

                                        class Standard_sequence
                                        {
                                            public:
                                                Standard_sequence(const Function& d, const This &k):
                                                ss_(k.exact_traits_object().standard_sequence_object(d.exact_function())){}
                                                Standard_sequence(){}
                                                typedef unsigned int result_type;
                                                typedef NT argument_type;

                                                template<class T>
                                                    unsigned int
                                                    number_of_real_roots(const T& a, const T& b) const
                                                {
                                                    return ss_.number_of_real_roots_base(a, b);
                                                }
                                            protected:
                                                typename Exact_traits::Standard_sequence ss_;
                                        };
                                        Standard_sequence standard_sequence_object(Function d) const
                                        {
                                            return Standard_sequence(d, *this);
                                        }

                                        class Sturm_sequence: public Exact_traits::Sturm_sequence
                                        {
                                            public:
                                                Sturm_sequence(const Function& a, const Function& b, const This &k):
                                                Exact_traits::Sturm_sequence(k.exact_traits_object().Sturm_sequence_object(a.exact_function(),
                                                    b.exact_function())){}
                                                Sturm_sequence(){}
                                        };
                                        Sturm_sequence Sturm_sequence_object(const Function &a, const Function &b) const
                                        {
//std::cout << "Stur of " << a << " and " << b << std::endl;
                                            return Sturm_sequence(a,b, *this);
                                        }

                                        typedef Filtered_root_bound_evaluator<This> Root_bound;
                                        Root_bound root_bound_object(bool pow_of_2=false) const
                                        {
                                            return Root_bound(pow_of_2, *this);
                                        }

                                        typedef internal::Sign_above_rational<This> Sign_above;
                                        Sign_above sign_above_object(const Function &f) const
                                        {
//std::cout << "Sign above of " << f << std::endl;
                                            return Sign_above(f, *this);
                                        }

                                        typedef internal::Sign_below_rational<This> Sign_below;
                                        Sign_below sign_below_object(const Function &f) const
                                        {
//std::cout << "Sign below of " << f << std::endl;
                                            return Sign_below(f, *this);
                                        }

                                        typedef Filtered_sign_at_rational<Function> Sign_at;
                                        Sign_at sign_at_object(const Function &f) const
                                        {
//std::cout << "Sign at of " << f << std::endl;
                                            return Sign_at(f);
                                        }

                                        typedef Filtered_are_negations<Function> Are_negations;
                                        Are_negations are_negations_object() const
                                        {
                                            return Are_negations();
                                        }

                                        typedef Filtered_rational_multiplicity<This> Multiplicity;
                                        Multiplicity multiplicity_object(const Function &f) const
                                        {
//std::cout << "Multiplicity of " << f << std::endl;
                                            return Multiplicity(f, *this);
                                        }

                                        const Exact_traits &exact_traits_object() const
                                        {
                                            return et_;
                                        }

                                        const Interval_traits &interval_traits_object() const
                                        {
                                            return it_;
                                        }

                                        const Exact_to_interval_converter exact_to_interval_converter_object() const
                                        {
                                            return efc_;
                                        }

                                        typedef internal::Construct_filtered_function<Function> Construct_function;
                                        Construct_function construct_function_object() const
                                        {
                                            return Construct_function(exact_to_interval_converter_object());
                                        }

                                        struct Compare_isolated_roots_in_interval
                                        {
                                            Compare_isolated_roots_in_interval(Function a, Function b, const This &k):
                                            pred_(k.exact_traits_object().compare_isolated_roots_in_interval_object(a.exact_function(), b.exact_function())) {
//std::cout << "Compare ir of " << a << " and " << b << std::endl;
                                            }

                                            typedef Comparison_result result_type;
                                            typedef NT first_argument_type;
                                            typedef NT second_argument_type;

                                            Comparison_result operator()(const NT &a, const NT &b) const
                                            {
                                                return pred_(a,b);
                                            }

                                            typename Exact_traits::Compare_isolated_roots_in_interval pred_;
                                        };

                                        Compare_isolated_roots_in_interval compare_isolated_roots_in_interval_object(Function a, Function b) const
                                        {
                                            return Compare_isolated_roots_in_interval(a,b, *this);
                                        }

                                    protected:
                                        Exact_traits et_;
                                        Interval_traits it_;
                                        Exact_to_interval_converter efc_;
                                };

                                } } } //namespace CGAL::POLYNOMIAL::internal
//#undef DF_UNARY_CONSTUCTION
#endif
