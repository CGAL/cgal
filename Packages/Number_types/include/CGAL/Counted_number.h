// ======================================================================
//
// Copyright (c) 2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Counted_number.h
// package       : Number_types
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Geert-Jan Giezeman
//
// coordinator   : MPI, Saarbruecken
// ======================================================================

#ifndef CGAL_COUNTED_NUMBER_H
#define CGAL_COUNTED_NUMBER_H

#include <CGAL/basic.h>
#include <CGAL/tags.h>
#include <iostream>

CGAL_BEGIN_NAMESPACE

template <class NT>
class Counted_number {
    static unsigned long s_neg_count, s_add_count, s_sub_count,
                         s_mul_count, s_div_count, s_sqrt_count,
			 s_eq_count, s_comp_count ;
    NT m_rep;
public:
    static void reset()
            { s_neg_count=0; s_add_count=0; s_sub_count=0;
              s_mul_count=0; s_div_count=0; s_sqrt_count=0;
	      s_eq_count=0; s_comp_count = 0;
            }
    static void inc_neg_count() {++s_neg_count;}
    static void inc_add_count() {++s_add_count;}
    static void inc_sub_count() {++s_sub_count;}
    static void inc_mul_count() {++s_mul_count;}
    static void inc_div_count() {++s_div_count;}
    static void inc_sqrt_count() {++s_sqrt_count;}
    static void inc_eq_count() {++s_eq_count;}
    static void inc_comp_count() {++s_comp_count;}

    static unsigned long neg_count() {return s_neg_count;}
    static unsigned long add_count() {return s_add_count;}
    static unsigned long sub_count() {return s_sub_count;}
    static unsigned long mul_count() {return s_mul_count;}
    static unsigned long div_count() {return s_div_count;}
    static unsigned long sqrt_count() {return s_sqrt_count;}
    static unsigned long eq_count() {return s_eq_count;}
    static unsigned long comp_count() {return s_comp_count;}
    static unsigned long count()
            { return s_neg_count + s_add_count + s_sub_count +
                     s_mul_count + s_div_count + s_sqrt_count +
	             s_eq_count + s_comp_count;
            }
    static void report(std::ostream &os);
    NT rep() const {return m_rep;}
    Counted_number() {}
    explicit Counted_number(int n) :m_rep(n){}
    explicit Counted_number(NT n) :m_rep(n){}
    Counted_number operator-() const
            {inc_neg_count();return Counted_number(-m_rep);}
    Counted_number const & operator+=(Counted_number const &n) 
            {
		inc_add_count();
		m_rep += n.m_rep;
		return *this;}
    Counted_number const & operator-=(Counted_number const &n) 
            {inc_sub_count(); m_rep -= n.m_rep; return *this;}
    Counted_number const & operator*=(Counted_number const &n) 
            {inc_mul_count(); m_rep *= n.m_rep; return *this;}
    Counted_number const & operator/=(Counted_number const &n) 
            {inc_div_count(); m_rep /= n.m_rep; return *this;}
};

template <class NT>
unsigned long Counted_number<NT>::s_neg_count=0;

template <class NT>
unsigned long Counted_number<NT>::s_add_count=0;

template <class NT>
unsigned long Counted_number<NT>::s_sub_count=0;

template <class NT>
unsigned long Counted_number<NT>::s_mul_count=0;

template <class NT>
unsigned long Counted_number<NT>::s_div_count=0;

template <class NT>
unsigned long Counted_number<NT>::s_sqrt_count=0;

template <class NT>
unsigned long Counted_number<NT>::s_eq_count=0;

template <class NT>
unsigned long Counted_number<NT>::s_comp_count=0;

template <class NT>
Counted_number<NT>
operator+(Counted_number<NT> const &n1, Counted_number<NT> const &n2)
{
    Counted_number<NT>::inc_add_count();
    return Counted_number<NT>(n1.rep() + n2.rep());
}

template <class NT>
Counted_number<NT>
operator-(Counted_number<NT> const &n1, Counted_number<NT> const &n2)
{
    Counted_number<NT>::inc_sub_count();
    return Counted_number<NT>(n1.rep() - n2.rep());
}

template <class NT>
Counted_number<NT>
operator*(Counted_number<NT> const &n1, Counted_number<NT> const &n2)
{
    Counted_number<NT>::inc_mul_count();
    return Counted_number<NT>(n1.rep() * n2.rep());
}

template <class NT>
Counted_number<NT>
operator/(Counted_number<NT> const &n1, Counted_number<NT> const &n2)
{
    Counted_number<NT>::inc_div_count();
    return Counted_number<NT>(n1.rep() / n2.rep());
}

template <class NT>
Counted_number<NT>
sqrt(Counted_number<NT> const &n1)
{
    Counted_number<NT>::inc_sqrt_count();
    return Counted_number<NT>(sqrt(n1.rep()));
}

template <class NT>
bool
operator==(Counted_number<NT> const &n1, Counted_number<NT> const &n2)
{
    Counted_number<NT>::inc_eq_count();
    return (n1.rep() == n2.rep());
}

template <class NT>
bool
operator!=(Counted_number<NT> const &n1, Counted_number<NT> const &n2)
{
    Counted_number<NT>::inc_eq_count();
    return (n1.rep() != n2.rep());
}

template <class NT>
bool
operator<(Counted_number<NT> const &n1, Counted_number<NT> const &n2)
{
    Counted_number<NT>::inc_comp_count();
    return (n1.rep() < n2.rep());
}

template <class NT>
bool
operator>(Counted_number<NT> const &n1, Counted_number<NT> const &n2)
{
    Counted_number<NT>::inc_comp_count();
    return (n1.rep() > n2.rep());
}

template <class NT>
bool
operator<=(Counted_number<NT> const &n1, Counted_number<NT> const &n2)
{
    Counted_number<NT>::inc_comp_count();
    return (n1.rep() <= n2.rep());
}

template <class NT>
bool
operator>=(Counted_number<NT> const &n1, Counted_number<NT> const &n2)
{
    Counted_number<NT>::inc_comp_count();
    return (n1.rep() >= n2.rep());
}

template <class NT>
inline double to_double(Counted_number<NT> n)
{
    return to_double(n.rep());
}

template <class NT>
inline bool is_finite(Counted_number<NT> n)
{
    return is_finite(n.rep());
}

template <class NT>
inline bool is_valid(Counted_number<NT> n)
{
    return is_valid(n.rep());
}

template <class NT>
void Counted_number<NT>::report(std::ostream &os)
{
    os << count() << " operations\n";
    if (neg_count() > 0)
        os << "  " << neg_count() << " negations\n";
    if (add_count() > 0)
        os << "  " << add_count() << " additions\n";
    if (sub_count() > 0)
        os << "  " << sub_count() << " subtractions\n";
    if (mul_count() > 0)
        os << "  " << mul_count() << " multiplications\n";
    if (div_count() > 0)
        os << "  " << div_count() << " divisions\n";
    if (sqrt_count() > 0)
        os << "  " << sqrt_count() << " square roots\n";
    if (eq_count() > 0)
        os << "  " << eq_count() << " equality tests\n";
    if (comp_count() > 0)
        os << "  " << comp_count() << " comparisons\n";
}

template <class NT>
inline io_Operator
io_tag(Counted_number<NT> const &n)
{ return io_tag(n.rep());}

template <class NT>
std::ostream& operator<<(std::ostream &os, Counted_number<NT> const &n)
{
    return os << n.rep();
}

template <class NT>
std::istream& operator>>(std::istream &is, Counted_number<NT> &n)
{
    NT num;
    is >> num;
    if (is) n = Counted_number<NT>(num);
    return is;
}

CGAL_END_NAMESPACE

#endif
