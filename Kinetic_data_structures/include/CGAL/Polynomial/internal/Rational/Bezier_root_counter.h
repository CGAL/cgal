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

#ifndef CGAL_BEZIER_ROOT_COUNTER_H
#define CGAL_BEZIER_ROOT_COUNTER_H

#include <CGAL/Polynomial/basic.h>
#include <vector>
#include <CGAL/Polynomial/internal/Sign_variations_counter.h>

namespace CGAL { namespace POLYNOMIAL { namespace internal {

/*class Bezier_root_count {
public:
  Bezier_root_count(int i):i_(i){
  }
  Bezier_root_count():i_(-100000){}
  bool operator==(const Bezier_root_count &o) const {
    if (i_ <0 && o.i_ < 0) return i_==o.i_;
    else if (i_ <0 || o.i_ <0){
      return (i_+4)%2 == (o.i_+4)%2;
    }
    return i_== o.i_;
}
bool operator!=(const Bezier_root_count &o) const {
return !(*this == o);
}
static Bezier_root_count zero() {
return Bezier_root_count(0);
}
static Bezier_root_count one() {
return Bezier_root_count(1);
}
static Bezier_root_count two() {
return Bezier_root_count(2);
}
static Bezier_root_count even() {
return Bezier_root_count(-4);
}
static Bezier_root_count odd() {
return Bezier_root_count(-3);
}
static Bezier_root_count single_even() {
return Bezier_root_count(-2);
}
static Bezier_root_count single_odd() {
return Bezier_root_count(-1);
}
bool is_odd() const {
return (i_+4)%2==1;
}
protected:
int i_;
};*/

#if 0
template<class K>
class Bezier_root_counter
{
    public:
        typedef typename K::Function                        Polynomial;
        typedef typename Polynomial::NT  NT;

    protected:
        typedef POLYNOMIAL_NS::Sign                Sign;
        typedef POLYNOMIAL_NS::Comparison_result   Comparison_result;
        typedef typename K::Sign_at      Sign_at;

    public:
        typedef NT                        argument_type;
        typedef NT                        argument_type1;
        typedef NT                        argument_type2;

        typedef unsigned int              result_type;

    protected:
//=====================================================================
// METHODS FOR COMPUTING THE BEZIER REPRESENTATION OF THE POLYNOMIAL
//=====================================================================
        static std::vector<NT>
        update_binom_coefs(const std::vector<NT>& coefs) {
            std::vector<NT> new_coefs(coefs.size() + 1);
            new_coefs[0] = 1;
            for (unsigned int i = 1; i < coefs.size(); i++) {
                new_coefs[i] = coefs[i - 1] + coefs[i];
            }
            new_coefs[coefs.size()] = 1;

            return new_coefs;
        }

        static
        std::vector<NT> compute_bezier_rep(const Polynomial& q) {
            std::vector<NT>  control(q.degree() + 1);

            std::vector<NT> binom_coefs;
            binom_coefs.push_back(NT(1));
            binom_coefs.push_back(NT(1));

            control[0] = q[0];

            int counter = q.degree();

            NT qq = NT(1);

            int k = 1;
            while ( counter > 0 ) {
                qq *= NT(k) / NT(counter);

                control[k] = q[k] * qq;
                for (int j = k-1, m = 1; j >= 0; j--, m++) {
                    if ( m % 2 == 1 ) {
                        control[k] += binom_coefs[m] * control[j];
                    }
                    else {
                        control[k] -= binom_coefs[m] * control[j];
                    }
                }

                counter--;
                k++;
                binom_coefs = update_binom_coefs(binom_coefs);
            }

            return control;
        }

//=====================================================================
// METHODS FOR PERFORMING SUBDIVISION OF THE BEZIER CURVE
//=====================================================================
        std::vector<NT>
            subdivide_left(const std::vector<NT>& c, const NT& t) const
        {
            NT tt = NT(1) - t;

            std::vector<NT> sub(c);

            for (unsigned int i = 1; i < c.size(); i++) {
                for (unsigned int j = c.size() - 1; j >= i; j--) {
                    sub[j] = tt * sub[j-1] + t * sub[j];
                }
            }

            return sub;
        }

        std::vector<NT>
            subdivide_right(const std::vector<NT>& c, const NT& t) const
//  subdivide_right(std::vector<NT> c, const NT& t)
        {
            NT tt = NT(1) - t;

            std::vector<NT> sub(c);

            for (unsigned int i = 1; i < c.size(); i++) {
                for (unsigned int j = 0; j < c.size() - i; j++) {
                    sub[j] = tt * sub[j] + t * sub[j+1];
                }
            }

            return sub;
        }

//============================================================
// METHOD FOR TRANSFORMING THE VALUES TO THE (0,1) INTERVAL
//============================================================
        static NT transform_value(const NT& x, Sign s_x) {
            if ( s_x == POLYNOMIAL_NS::NEGATIVE ) {
                return  NT(1) / (NT(1) - x);
            }
            return  NT(1) / (NT(1) + x);
        }

    protected:
//===========================================
// METHODS FOR COMPUTING THE SIGN VARIATIONS
//===========================================
        template<class Iterator>
            static
            unsigned int
            sign_variations(const Iterator& first, const Iterator& beyond,
        bool stop_if_more_than_one = false) {
            return
                Sign_variations_counter::sign_variations(first, beyond,
                stop_if_more_than_one);
        }

        unsigned int sign_variations(const NT& x) const
        {
            Sign s_x = POLYNOMIAL_NS::sign(x);

            NT y = transform_value(x, s_x);

            std::vector<NT> sub;
            if ( s_x == POLYNOMIAL_NS::NEGATIVE ) {
                sub = subdivide(control_neg, y);
            }
            else {
                sub = subdivide(control_pos, y);
            }

            std::vector<Sign> signs;
            for (unsigned int i = 0; i < sub.size(); i++) {
                signs[i] = POLYNOMIAL_NS::sign( sub[i] );
            }

            return sign_variations(signs.begin(), signs.end());
        }

//==============================================================
// METHOD FOR COMPUTING THE NUMBER OF REAL ROOTS IN AN INTERVAL
//==============================================================

        unsigned int
            number_of_roots_base(const NT& a, const NT& b,
            bool is_positive_interval,
            unsigned int nr_ub) const
        {
            Sign_at sign_at(p);
            Sign s_a = sign_at(a);
            Sign s_b = sign_at(b);

            if ( s_a != POLYNOMIAL_NS::ZERO && s_b != POLYNOMIAL_NS::ZERO && nr_ub == 1 ) {
                return ( s_a == s_b ) ? 0 : 1;
            }

            return number_of_roots_base(a, b, is_positive_interval);
        }

        unsigned int
            number_of_roots_base(const NT& a, const NT& b,
            bool is_positive_interval) const
        {
            NT aa, bb;

            if ( is_positive_interval ) {
                aa = transform_value(a, POLYNOMIAL_NS::POSITIVE);
                bb = transform_value(b, POLYNOMIAL_NS::POSITIVE);
            }
            else {
                aa = transform_value(a, POLYNOMIAL_NS::NEGATIVE);
                bb = transform_value(b, POLYNOMIAL_NS::NEGATIVE);
            }

#if 1                                     // this may create problems when I do filtering
            if ( is_positive_interval ) {
                Polynomial_assertion( bb <= aa && bb >= NT(0) && aa <= NT(1) );
            }
            else {
                Polynomial_assertion( aa <= bb && bb >= NT(0) && aa <= NT(1) );
            }
#endif

            NT x, y;
            if ( is_positive_interval ) {
                x = aa;
                y = bb / aa;
            }
            else {
                x = bb;
                y = aa / bb;
            }

            std::vector<NT> sub;

            if ( is_positive_interval ) {
                sub = subdivide_left(control_pos, x);
            }
            else {
                sub = subdivide_left(control_neg, x);
            }

            std::vector<NT> sub2 = subdivide_right(sub, y);

            std::vector<Sign> signs(sub2.size());
            for (unsigned int i = 0; i < sub2.size(); i++) {
                signs[i] = POLYNOMIAL_NS::sign( sub2[i] );
            }

            return sign_variations(signs.begin(), signs.end(), false);
        }

//=====================================================
// TRANFORM THE POLYNOMIAL SO THAT ZERO IS NOT A ROOT;
// WE KEEP TRACK OF ZERO'S MULTIPLICITY
//=====================================================
        void remove_zero() {
            multiplicity_of_zero_ = 0;

            if ( p.degree() == -1 ) { return; }

            int k = 0;
            Sign s = POLYNOMIAL_NS::sign( p[k] );

            while ( k <= p.degree() && s == POLYNOMIAL_NS::ZERO ) {
                k++;
                s = POLYNOMIAL_NS::sign( p[k] );

                multiplicity_of_zero_++;
            }

            Polynomial q;

            int qdeg = p.degree() - multiplicity_of_zero_;

            q.hint_degree( qdeg );

            for (int i = 0; i < qdeg; i++) {
                q.set_coef(i, p[i + multiplicity_of_zero_]);
            }

            p = q;
        }

    public:
//================
// CONSTRUCTORS
//================
        Bezier_root_counter() : p() {bool please_switch_to_kernel_version;}

        Bezier_root_counter(const Polynomial& p)
        : p(p) {
//bool please_switch_to_kernel_version;
            remove_zero();

            if (p.degree() == -1) { return; }

//    CGAL_precondition( CGAL::sign(p[0]) != CGAL::ZERO );

            Polynomial q_pos(p);
            q_pos = translate_zero(q_pos, NT(-1));
            q_pos = invert_variable(q_pos);

            Polynomial q_neg = p.negated_variable();
            q_neg = translate_zero(q_neg, NT(-1));
            q_neg = invert_variable(q_neg);

            control_pos = compute_bezier_rep(q_pos);
            control_neg = compute_bezier_rep(q_neg);

#if 0
            std::cout << "positive control polygon:" << std::endl;
            for (int i = 0; i < control_pos.size(); i++) {
                std::cout << control_pos[i] << " ";
            }
            std::cout << std::endl;
            std::cout << std::endl;

            std::cout << "negative control polygon:" << std::endl;
            for (int i = 0; i < control_neg.size(); i++) {
                std::cout << control_neg[i] << " ";
            }
            std::cout << std::endl;
            std::cout << std::endl;
#endif
        }

//=================================================================
// OPERATOR FOR COMPUTING THE NUMBER OF REAL ROOTS IN AN INTERVAL
//=================================================================
#if 0
        unsigned int operator()(const NT& a, const NT& b,
            unsigned int nr_ub) const
        {
            POLYNOMIAL_NS::Sign s_a = POLYNOMIAL_NS::sign(a);
            POLYNOMIAL_NS::Sign s_b = POLYNOMIAL_NS::sign(b);

            if ( s_a == POLYNOMIAL_NS::ZERO ) {
                return number_of_roots_base(a, b, true, nr_ub);
            }

            if ( s_b == POLYNOMIAL_NS::ZERO ) {
                return number_of_roots_base(a, b, false, nr_ub);
            }

            if ( s_a == POLYNOMIAL_NS::NEGATIVE && s_b == POLYNOMIAL_NS::POSITIVE ) {
                unsigned int nneg = number_of_roots_base(a, NT(0), false, nr_ub);
                unsigned int npos = number_of_roots_base(NT(0), b, true, nr_ub);
                return nneg + npos + multiplicity_of_zero_;
            }

            return number_of_roots_base(a, b, s_a == POLYNOMIAL_NS::POSITIVE, nr_ub);
        }
#endif

        unsigned int operator()(const NT& a, const NT& b) const
        {
            Sign s_a = POLYNOMIAL_NS::sign(a);
            Sign s_b = POLYNOMIAL_NS::sign(b);

            if ( s_a == POLYNOMIAL_NS::ZERO ) {
                return result_type(number_of_roots_base(a, b, true));
            }

            if ( s_b == POLYNOMIAL_NS::ZERO ) {
                return result_type(number_of_roots_base(a, b, false));
            }

            if ( s_a == POLYNOMIAL_NS::NEGATIVE && s_b == POLYNOMIAL_NS::POSITIVE ) {
                unsigned int nneg = number_of_roots_base(a, NT(0), false);
                unsigned int npos = number_of_roots_base(NT(0), b, true);

                return result_type( nneg + npos + multiplicity_of_zero_);
            }

            return result_type(number_of_roots_base(a, b, s_a == POLYNOMIAL_NS::POSITIVE));
        }

#if 0
//=====================================================
// THE NUMBER OF REAL ROOTS IN THE INTERVAL (0, +oo)
//=====================================================
        unsigned int number_of_positive_roots() const
        {
            std::vector<POLYNOMIAL_NS::Sign> signs(control_pos.size());
            for (unsigned int i = 0; i < control_pos.size(); i++) {
                signs[i] = POLYNOMIAL_NS::sign( control_pos[i] );
            }

            int sv = sign_variations(signs.begin(), signs.end(), false);
            return static_cast<unsigned int>(sv);
        }
#endif

        unsigned int& multiplicity_of_zero() {
            return multiplicity_of_zero_;
        }
        unsigned int multiplicity_of_zero() const
        {
            return multiplicity_of_zero_;
        }

        Polynomial& polynomial() { return p; }
        const Polynomial& polynomial() const { return p; }

        std::vector<NT>& positive_control_polygon() { return control_pos; }
        const std::vector<NT>& positive_control_polygon() const
        {
            return control_pos;
        }

        std::vector<NT>& negative_control_polygon() { return control_neg; }
        const std::vector<NT>& negative_control_polygon() const
        {
            return control_neg;
        }

        bool is_valid() const { return is_valid_; }
        bool& is_valid() { return is_valid_; }

    protected:
        bool             is_valid_;
        unsigned int     multiplicity_of_zero_;
        Polynomial       p;
        std::vector<NT>  control_pos;
        std::vector<NT>  control_neg;
};
#endif

template<class Kernel>
class Bezier_root_counter
{
    public:
        typedef typename Kernel::Function Polynomial;
        typedef typename Polynomial::NT  NT;

    protected:
        typedef POLYNOMIAL_NS::Sign                Sign;
        typedef POLYNOMIAL_NS::Comparison_result   Comparison_result;
        typedef typename Kernel::Sign_at       Sign_at;

    public:
        typedef NT                        argument_type;
        typedef NT                        argument_type1;
        typedef NT                        argument_type2;

        typedef unsigned int              result_type;

    protected:
//=====================================================================
// METHODS FOR COMPUTING THE BEZIER REPRESENTATION OF THE POLYNOMIAL
//=====================================================================
        static std::vector<NT>
        update_binom_coefs(const std::vector<NT>& coefs) {
            std::vector<NT> new_coefs(coefs.size() + 1);
            new_coefs[0] = 1;
            for (unsigned int i = 1; i < coefs.size(); i++) {
                new_coefs[i] = coefs[i - 1] + coefs[i];
            }
            new_coefs[coefs.size()] = 1;

            return new_coefs;
        }

        static
        std::vector<NT> compute_bezier_rep(const Polynomial& q) {
            std::vector<NT>  control(q.degree() + 1);

            std::vector<NT> binom_coefs;
            binom_coefs.push_back(NT(1));
            binom_coefs.push_back(NT(1));

            control[0] = q[0];

            int counter = q.degree();

            NT qq = NT(1);

            int k = 1;
            while ( counter > 0 ) {
                qq *= NT(k) / NT(counter);

                control[k] = q[k] * qq;
                for (int j = k-1, m = 1; j >= 0; j--, m++) {
                    if ( m % 2 == 1 ) {
                        control[k] += binom_coefs[m] * control[j];
                    }
                    else {
                        control[k] -= binom_coefs[m] * control[j];
                    }
                }

                counter--;
                k++;
                binom_coefs = update_binom_coefs(binom_coefs);
            }

            return control;
        }

//=====================================================================
// METHODS FOR PERFORMING SUBDIVISION OF THE BEZIER CURVE
//=====================================================================
        std::vector<NT>
            subdivide_left(const std::vector<NT>& c, const NT& t) const
        {
            NT tt = NT(1) - t;

            std::vector<NT> sub(c);

            for (unsigned int i = 1; i < c.size(); i++) {
                for (unsigned int j = c.size() - 1; j >= i; j--) {
                    sub[j] = tt * sub[j-1] + t * sub[j];
                }
            }

            return sub;
        }

        std::vector<NT>
            subdivide_right(const std::vector<NT>& c, const NT& t) const
//  subdivide_right(std::vector<NT> c, const NT& t)
        {
            NT tt = NT(1) - t;

            std::vector<NT> sub(c);

            for (unsigned int i = 1; i < c.size(); i++) {
                for (unsigned int j = 0; j < c.size() - i; j++) {
                    sub[j] = tt * sub[j] + t * sub[j+1];
                }
            }

            return sub;
        }

//============================================================
// METHOD FOR TRANSFORMING THE VALUES TO THE (0,1) INTERVAL
//============================================================
        static NT transform_value(const NT& x, Sign s_x) {
            if ( s_x == POLYNOMIAL_NS::NEGATIVE ) {
                return  NT(1) / (NT(1) - x);
            }
            return  NT(1) / (NT(1) + x);
        }

    protected:
//===========================================
// METHODS FOR COMPUTING THE SIGN VARIATIONS
//===========================================
        template<class Iterator>
            static
            unsigned int
            sign_variations(const Iterator& first, const Iterator& beyond,
        bool stop_if_more_than_one = false) {
            return
                Sign_variations_counter::sign_variations(first, beyond,
                stop_if_more_than_one);
        }

        unsigned int sign_variations(const NT& x) const
        {
            Sign s_x = POLYNOMIAL_NS::sign(x);

            NT y = transform_value(x, s_x);

            std::vector<NT> sub;
            if ( s_x == POLYNOMIAL_NS::NEGATIVE ) {
                sub = subdivide(control_neg, y);
            }
            else {
                sub = subdivide(control_pos, y);
            }

            std::vector<Sign> signs;
            for (unsigned int i = 0; i < sub.size(); i++) {
                signs[i] = POLYNOMIAL_NS::sign( sub[i] );
            }

            return sign_variations(signs.begin(), signs.end());
        }

//==============================================================
// METHOD FOR COMPUTING THE NUMBER OF REAL ROOTS IN AN INTERVAL
//==============================================================

        unsigned int
            number_of_roots_base(const NT& a, const NT& b,
            bool is_positive_interval,
            unsigned int nr_ub) const
        {
            Sign_at sign_at= kernel_.sign_at_object(p);
            Sign s_a = sign_at(a);
            Sign s_b = sign_at(b);

            if ( s_a != POLYNOMIAL_NS::ZERO && s_b != POLYNOMIAL_NS::ZERO && nr_ub == 1 ) {
                return ( s_a == s_b ) ? 0 : 1;
            }

            return  number_of_roots_base(a, b, is_positive_interval);
        }

        unsigned int
            number_of_roots_base(const NT& a, const NT& b,
            bool is_positive_interval) const
        {
            NT aa, bb;

            if ( is_positive_interval ) {
                aa = transform_value(a, POLYNOMIAL_NS::POSITIVE);
                bb = transform_value(b, POLYNOMIAL_NS::POSITIVE);
            }
            else {
                aa = transform_value(a, POLYNOMIAL_NS::NEGATIVE);
                bb = transform_value(b, POLYNOMIAL_NS::NEGATIVE);
            }

#if 0                                     // this may create problems when I do filtering
            if ( is_positive_interval ) {
                Polynomial_assertion( bb <= aa && bb >= NT(0) && aa <= NT(1) );
            }
            else {
                Polynomial_assertion( aa <= bb && bb >= NT(0) && aa <= NT(1) );
            }
#endif

            NT x, y;
            if ( is_positive_interval ) {
                x = aa;
                y = bb / aa;
            }
            else {
                x = bb;
                y = aa / bb;
            }

            std::vector<NT> sub;

            if ( is_positive_interval ) {
                sub = subdivide_left(control_pos, x);
            }
            else {
                sub = subdivide_left(control_neg, x);
            }

            std::vector<NT> sub2 = subdivide_right(sub, y);

            std::vector<Sign> signs(sub2.size());
            for (unsigned int i = 0; i < sub2.size(); i++) {
                signs[i] = POLYNOMIAL_NS::sign( sub2[i] );
            }

            return sign_variations(signs.begin(), signs.end(), false);
        }

//=====================================================
// TRANFORM THE POLYNOMIAL SO THAT ZERO IS NOT A ROOT;
// WE KEEP TRACK OF ZERO'S MULTIPLICITY
//=====================================================
        void remove_zero() {
            multiplicity_of_zero_ = 0;

            if ( p.degree() == -1 ) { return; }

            int k = 0;
            Sign s = POLYNOMIAL_NS::sign( p[k] );

            while ( k <= p.degree() && s == POLYNOMIAL_NS::ZERO ) {
                k++;
                s = POLYNOMIAL_NS::sign( p[k] );

                multiplicity_of_zero_++;
            }

            Polynomial q;

            int qdeg = p.degree() - multiplicity_of_zero_;

            q.hint_degree( qdeg );

            for (int i = 0; i < qdeg; i++) {
                q.set_coef(i, p[i + multiplicity_of_zero_]);
            }

            p = q;
        }

    public:
//================
// CONSTRUCTORS
//================
        Bezier_root_counter() : p() {}

        Bezier_root_counter(const Polynomial& p, Kernel k= Kernel())
        : p(p), kernel_(k) {
            remove_zero();

            if (p.degree() == -1) { return; }

//    POLYNOMIAL_NS_precondition( POLYNOMIAL_NS::sign(p[0]) != POLYNOMIAL_NS::ZERO );

            Polynomial q_pos(p);
            q_pos = kernel_.translate_zero_object(q_pos)( NT(-1));
            q_pos = kernel_.invert_variable_object()(q_pos);

            Polynomial q_neg = p.negated_variable();
            q_neg = kernel_.translate_zero_object(q_neg)(NT(-1));
            q_neg = kernel_.invert_variable_object()(q_neg);

            control_pos = compute_bezier_rep(q_pos);
            control_neg = compute_bezier_rep(q_neg);

        }

//=================================================================
// OPERATOR FOR COMPUTING THE NUMBER OF REAL ROOTS IN AN INTERVAL
//=================================================================
        result_type operator()(const NT& a, const NT& b,
            POLYNOMIAL_NS::Sign = POLYNOMIAL_NS::POSITIVE,
            POLYNOMIAL_NS::Sign = POLYNOMIAL_NS::POSITIVE) const
        {
            unsigned int ret;

            POLYNOMIAL_NS::Sign s_a=POLYNOMIAL_NS::sign(a);
            POLYNOMIAL_NS::Sign s_b=POLYNOMIAL_NS::sign(b);

            if ( s_a == POLYNOMIAL_NS::ZERO ) {
                ret= result_type(number_of_roots_base(a, b, true));
            }
            else if ( s_b == POLYNOMIAL_NS::ZERO ) {
                ret= result_type(number_of_roots_base(a, b, false));
            }
            else if ( s_a == POLYNOMIAL_NS::NEGATIVE && s_b == POLYNOMIAL_NS::POSITIVE ) {
                unsigned int nneg = number_of_roots_base(a, NT(0), false);
                unsigned int npos = number_of_roots_base(NT(0), b, true);

                ret= result_type(nneg + npos + multiplicity_of_zero_);
            }
            else {
                ret= result_type(number_of_roots_base(a, b, s_a == POLYNOMIAL_NS::POSITIVE));
            }
            return ret;
        }
/*
unsigned int& multiplicity_of_zero() {
  return multiplicity_of_zero_;
}
unsigned int multiplicity_of_zero() const {
  return multiplicity_of_zero_;
  }*/

/*
Polynomial& polynomial() { return p; }
const Polynomial& polynomial() const { return p; }*/

        std::vector<NT>& positive_control_polygon() { return control_pos; }
        const std::vector<NT>& positive_control_polygon() const
        {
            return control_pos;
        }

        std::vector<NT>& negative_control_polygon() { return control_neg; }
        const std::vector<NT>& negative_control_polygon() const
        {
            return control_neg;
        }

/*bool is_valid() const { return is_valid_; }
  bool& is_valid() { return is_valid_; }*/

    protected:
        bool             is_valid_;
        unsigned int     multiplicity_of_zero_;
        Polynomial       p;
        std::vector<NT>  control_pos;
        std::vector<NT>  control_neg;
        Kernel kernel_;
};

} } } //namespace CGAL::POLYNOMIAL::internal
#endif                                            // CGAL_BEZIER_ROOT_COUNTER_H
