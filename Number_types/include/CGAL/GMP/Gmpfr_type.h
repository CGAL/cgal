// Copyright (c) 2007-2010 Inria Lorraine (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

#ifndef CGAL_GMPFR_TYPE_H
#define CGAL_GMPFR_TYPE_H

#include <CGAL/disable_warnings.h>

#include <CGAL/gmp.h>
#include <mpfr.h>
#include <boost/operators.hpp>
#include <CGAL/Handle_for.h>
#include <CGAL/GMP/Gmpz_type.h>
#include <CGAL/GMP/Gmpzf_type.h>
#include <limits>
#include <CGAL/Uncertain.h>
#include <CGAL/ipower.h>
#include <CGAL/IO/io.h>

#if MPFR_VERSION_MAJOR < 3
        typedef mp_rnd_t mpfr_rnd_t;
        typedef mp_prec_t mpfr_prec_t;
        typedef mp_exp_t mpfr_exp_t;
        #define MPFR_RNDN GMP_RNDN
        #define MPFR_RNDZ GMP_RNDZ
        #define MPFR_RNDU GMP_RNDU
        #define MPFR_RNDD GMP_RNDD
        #define CGAL_GMPFR_GET_Z_2EXP mpfr_get_z_exp
#else
        #define CGAL_GMPFR_GET_Z_2EXP mpfr_get_z_2exp
#endif

namespace CGAL{

class Gmpfr;

bool operator<(const Gmpfr&,const Gmpfr&);
bool operator==(const Gmpfr&,const Gmpfr&);

bool operator<(const Gmpfr&,long);
bool operator>(const Gmpfr&,long);
bool operator==(const Gmpfr&,long);

bool operator<(const Gmpfr&,unsigned long);
bool operator>(const Gmpfr&,unsigned long);
bool operator==(const Gmpfr&,unsigned long);

bool operator<(const Gmpfr&,int);
bool operator>(const Gmpfr&,int);
bool operator==(const Gmpfr&,int);

bool operator<(const Gmpfr&,double);
bool operator>(const Gmpfr&,double);
bool operator==(const Gmpfr&,double);

bool operator<(const Gmpfr&,long double);
bool operator>(const Gmpfr&,long double);
bool operator==(const Gmpfr&,long double);

bool operator<(const Gmpfr&,const Gmpz&);
bool operator>(const Gmpfr&,const Gmpz&);
bool operator==(const Gmpfr&,const Gmpz&);

struct Gmpfr_rep{
        mpfr_t floating_point_number;
        bool clear_on_destruction;
        Gmpfr_rep():clear_on_destruction(true){}
        ~Gmpfr_rep(){
                if(clear_on_destruction)
                        mpfr_clear(floating_point_number);
        }
};

namespace internal{
        template <>
        struct Minmax_traits<mpfr_rnd_t>{
                static const mpfr_rnd_t min=MPFR_RNDN;
        #if MPFR_VERSION_MAJOR < 3
                static const mpfr_rnd_t max=GMP_RND_MAX;
        #else
                static const mpfr_rnd_t max=MPFR_RNDF;
        #endif
        };
} // namespace internal

// The class Gmpfr is reference-counted using CGAL's handle mechanism. This
// behavior may be changed, setting the flag CGAL_GMPFR_NO_REFCOUNT. A
// non-reference-counted class is slightly more efficient in case the
// implementation does not need to copy numbers (this is not usually the
// case). Nevertheless, setting the flag may be useful for debugging
// purposes.

class Gmpfr:
#ifdef CGAL_GMPFR_NO_REFCOUNT
        Gmpfr_rep,
#else
        Handle_for<Gmpfr_rep>,
#endif
        boost::ordered_euclidian_ring_operators1<Gmpfr,
        boost::ordered_euclidian_ring_operators2<Gmpfr,long,
        boost::ordered_euclidian_ring_operators2<Gmpfr,unsigned long,
        boost::ordered_euclidian_ring_operators2<Gmpfr,int,
        boost::totally_ordered2<Gmpfr,double,
        boost::totally_ordered2<Gmpfr,long double,
        boost::ordered_euclidian_ring_operators2<Gmpfr,Gmpz
        > > > > > > >
{
        private:

#ifndef CGAL_GMPFR_NO_REFCOUNT
        typedef Handle_for<Gmpfr_rep>   Base;
#endif

        static Uncertain<mpfr_rnd_t> _gmp_rnd(std::float_round_style r){
                switch(r){
                        case std::round_toward_infinity: return MPFR_RNDU;
                        case std::round_toward_neg_infinity: return MPFR_RNDD;
                        case std::round_toward_zero: return MPFR_RNDZ;
                        case std::round_to_nearest: return MPFR_RNDN;
                        default: return Uncertain<mpfr_rnd_t>::indeterminate();
                }
        };

        static std::float_round_style _cgal_rnd(mpfr_rnd_t r){
                switch(r){
                        case MPFR_RNDU: return std::round_toward_infinity;
                        case MPFR_RNDD: return std::round_toward_neg_infinity;
                        case MPFR_RNDZ: return std::round_toward_zero;
                        case MPFR_RNDN: return std::round_to_nearest;
                        default: return std::round_indeterminate;
                }
        };

        public:

        typedef mpfr_prec_t             Precision_type;

        // access

        mpfr_srcptr fr()const{
#ifdef CGAL_GMPFR_NO_REFCOUNT
                return floating_point_number;
#else
                return Ptr()->floating_point_number;
#endif
        }

        mpfr_ptr fr(){
#ifdef CGAL_GMPFR_NO_REFCOUNT
                return floating_point_number;
#else
                return ptr()->floating_point_number;
#endif
        }

        // The function dont_clear_on_destruction() is used to tell the
        // object that the mpfr_t must not be cleared at object destruction
        // (the destructor contrasts with that of the GMP types, where the
        // structure is always cleared). The reason to do this is that,
        // sometimes, the object is constructed from a mpfr_t structure and
        // we know that the structure will be cleared somewhere else. The
        // mpfr_t will not need to be cleared in case the object A is
        // constructed from another Gmpfr object B, specifying the
        // precision, and it happens that this precision is the same as the
        // precision of A. In this case, a shallow copy is constructed.

        void dont_clear_on_destruction(){
#ifdef CGAL_GMPFR_NO_REFCOUNT
                clear_on_destruction=false;
#else
                ptr()->clear_on_destruction=false;
#endif
        }

        bool is_unique(){
#ifdef CGAL_GMPFR_NO_REFCOUNT
                return true;
#else
                return unique();
#endif
        }

        // construction

        Gmpfr(){
                mpfr_init(fr());
        }

        Gmpfr(mpfr_srcptr f){
                mpfr_custom_init_set(
                        fr(),
                        mpfr_custom_get_kind(f),
                        mpfr_custom_get_exp(f),
                        mpfr_get_prec(f),
                        mpfr_custom_get_mantissa(f));
                dont_clear_on_destruction();
                CGAL_assertion((mpfr_nan_p(f)!=0 && mpfr_nan_p(fr())!=0) ||
                               (mpfr_unordered_p(f,fr())==0 &&
                                mpfr_equal_p(f,fr())!=0));
        }

        Gmpfr(mpfr_srcptr f,
              std::float_round_style r,
              Gmpfr::Precision_type p=Gmpfr::get_default_precision()){
                CGAL_assertion(p>=MPFR_PREC_MIN&&p<=MPFR_PREC_MAX);
                if(p==mpfr_get_prec(f)){
                        mpfr_custom_init_set(
                                fr(),
                                mpfr_custom_get_kind(f),
                                mpfr_custom_get_exp(f),
                                mpfr_get_prec(f),
                                mpfr_custom_get_mantissa(f));
                        dont_clear_on_destruction();
                        CGAL_assertion((mpfr_nan_p(f)!=0&&mpfr_nan_p(fr())!=0)||
                                       (mpfr_unordered_p(f,fr())==0&&
                                        mpfr_equal_p(f,fr())!=0));
                }else{
                        mpfr_init2(fr(),p);
                        mpfr_set(fr(),f,_gmp_rnd(r));
                        CGAL_assertion(mpfr_get_prec(fr())<mpfr_get_prec(f)||
                                       mpfr_equal_p(fr(),f)!=0);
                }
        }

        Gmpfr(mpfr_srcptr f,Gmpfr::Precision_type p){
                CGAL_assertion(p>=MPFR_PREC_MIN&&p<=MPFR_PREC_MAX);
                if(p==mpfr_get_prec(f)){
                        mpfr_custom_init_set(
                                fr(),
                                mpfr_custom_get_kind(f),
                                mpfr_custom_get_exp(f),
                                mpfr_get_prec(f),
                                mpfr_custom_get_mantissa(f));
                        dont_clear_on_destruction();
                        CGAL_assertion((mpfr_nan_p(f)!=0&&mpfr_nan_p(fr())!=0)||
                                       (mpfr_unordered_p(f,fr())==0&&
                                        mpfr_equal_p(f,fr())!=0));
                }else{
                        mpfr_init2(fr(),p);
                        mpfr_set(fr(),f,mpfr_get_default_rounding_mode());
                        CGAL_assertion(p<mpfr_get_prec(f)||
                                       mpfr_equal_p(fr(),f)!=0);
                }
        }

        Gmpfr(const Gmpzf &f,
              std::float_round_style r,
              Gmpfr::Precision_type p=Gmpfr::get_default_precision()){
                CGAL_assertion(p>=MPFR_PREC_MIN&&p<=MPFR_PREC_MAX);
                mpfr_init2(fr(),p);
                mpfr_set_z(fr(),f.man(),_gmp_rnd(r));
                mpfr_mul_2si(fr(),fr(),f.exp(),_gmp_rnd(r));
        }

        Gmpfr(const Gmpzf &f,Gmpfr::Precision_type p){
                CGAL_assertion(p>=MPFR_PREC_MIN&&p<=MPFR_PREC_MAX);
                mpfr_init2(fr(),p);
                mpfr_set_z(fr(),f.man(),mpfr_get_default_rounding_mode());
                mpfr_mul_2si(fr(),
                             fr(),
                             f.exp(),
                             mpfr_get_default_rounding_mode());
        }

        Gmpfr(const Gmpzf &f){
                mpfr_init2(fr(),
                           static_cast<Gmpfr::Precision_type>(
                                   mpz_sizeinbase(f.man(),2)<MPFR_PREC_MIN?
                                   MPFR_PREC_MIN:
                                   mpz_sizeinbase(f.man(),2)));
                mpfr_set_z(fr(),f.man(),MPFR_RNDN);
                CGAL_assertion_msg(mpfr_cmp_z(fr(),f.man())==0,
                                   "inexact conversion of a Gmpzf mantissa");
                CGAL_assertion_code(int inexact=)
                mpfr_mul_2si(fr(),fr(),f.exp(),MPFR_RNDN);
                CGAL_assertion_msg(inexact==0,"inexact conversion from Gmpzf");
        }

        Gmpfr(const std::pair<Gmpz,long> &intexp,
              std::float_round_style r=Gmpfr::get_default_rndmode(),
              Gmpfr::Precision_type p=Gmpfr::get_default_precision()){
                CGAL_assertion(p>=MPFR_PREC_MIN&&p<=MPFR_PREC_MAX);
                mpfr_init2(fr(),p);
                mpfr_set_z(fr(),intexp.first.mpz(),_gmp_rnd(r));
                mpfr_mul_2si(fr(),fr(),intexp.second,_gmp_rnd(r));
        }

        Gmpfr(const std::pair<Gmpz,long> &intexp,Gmpfr::Precision_type p){
                CGAL_assertion(p>=MPFR_PREC_MIN&&p<=MPFR_PREC_MAX);
                mpfr_init2(fr(),p);
                mpfr_set_z(fr(),
                           intexp.first.mpz(),
                           mpfr_get_default_rounding_mode());
                mpfr_mul_2si(fr(),
                             fr(),
                             intexp.second,
                             mpfr_get_default_rounding_mode());
        }

#define CGAL_GMPFR_CONSTRUCTOR_FROM_TYPE(_type,_fun) \
        Gmpfr(_type x, \
              std::float_round_style r, \
              Gmpfr::Precision_type p=Gmpfr::get_default_precision()){ \
                CGAL_assertion(p>=MPFR_PREC_MIN&&p<=MPFR_PREC_MAX); \
                mpfr_init2(fr(),p); \
                _fun(fr(),x,_gmp_rnd(r)); \
        } \
        Gmpfr(_type x,Gmpfr::Precision_type p){ \
                CGAL_assertion(p>=MPFR_PREC_MIN&&p<=MPFR_PREC_MAX); \
                mpfr_init2(fr(),p); \
                _fun(fr(),x,mpfr_get_default_rounding_mode()); \
        } \
        Gmpfr(_type x){ \
                mpfr_init2(fr(),mp_bits_per_limb*sizeof(_type)); \
                _fun(fr(),x,mpfr_get_default_rounding_mode()); \
        }

        CGAL_GMPFR_CONSTRUCTOR_FROM_TYPE(int,mpfr_set_si);
        CGAL_GMPFR_CONSTRUCTOR_FROM_TYPE(long,mpfr_set_si);
        CGAL_GMPFR_CONSTRUCTOR_FROM_TYPE(unsigned,mpfr_set_ui);
        CGAL_GMPFR_CONSTRUCTOR_FROM_TYPE(unsigned long,mpfr_set_ui);
        CGAL_GMPFR_CONSTRUCTOR_FROM_TYPE(double,mpfr_set_d);

        // With the MSVC compiler, 'long double' and 'double' are two
        // different types, but with the same size: sizeof(long double)==8.
        // For that reason, the cast of a long double to a double is
        // exact.
        // What is more, if one compile the mpfr library with mingw(32|64),
        // on Windows, this compiler has sizeof(long double)==16, as
        // gcc/g++ on Linux, and the produces libmpfr-1.dll has a symbol
        // mpfr_set_ld which is binary incompatible with a call from MSVC.
        // For those two reason, the constructor from 'long
        // double' calls 'mpfr_set_l' on MSVC, instead of 'mpfr_set_ld'.
        // That should not modify the semantic of a CGAL program, but
        // only avoid the binary incompatibility of a CGAL program compiled
        // with MSVC with the libmpfr-1.dll compiled with mingw.
#ifdef _MSC_VER
#  pragma warning(push)
#  pragma warning(disable: 4244)
        CGAL_GMPFR_CONSTRUCTOR_FROM_TYPE(long double,mpfr_set_d);
#  pragma warning(pop)
#else
        CGAL_GMPFR_CONSTRUCTOR_FROM_TYPE(long double,mpfr_set_ld);
#endif
#undef CGAL_GMPFR_CONSTRUCTOR_FROM_TYPE

#define CGAL_GMPFR_CONSTRUCTOR_FROM_OBJECT(_class,_member,_fun,_preccode) \
        Gmpfr(const _class &x, \
              std::float_round_style r, \
              Gmpfr::Precision_type p=Gmpfr::get_default_precision()){ \
                CGAL_assertion(p>=MPFR_PREC_MIN&&p<=MPFR_PREC_MAX); \
                mpfr_init2(fr(),p); \
                _fun(fr(),x._member,_gmp_rnd(r)); \
        } \
        Gmpfr(const _class &x,Gmpfr::Precision_type p){ \
                CGAL_assertion(p<=MPFR_PREC_MAX); \
                mpfr_init2(fr(),MPFR_PREC_MIN<p?p:MPFR_PREC_MIN); \
                _fun(fr(),x._member,mpfr_get_default_rounding_mode()); \
        } \
        Gmpfr(const _class &x){ \
                Gmpfr::Precision_type p=(_preccode); \
                mpfr_init2(fr(),MPFR_PREC_MIN<p?p:MPFR_PREC_MIN); \
                _fun(fr(),x._member,MPFR_RNDN); \
        }

        CGAL_GMPFR_CONSTRUCTOR_FROM_OBJECT(Gmpz,
                                           mpz(),
                                           mpfr_set_z,
                                           static_cast<Gmpfr::Precision_type>(
                                                x.bit_size()));

#undef CGAL_GMPFR_CONSTRUCTOR_FROM_OBJECT

        // When Gmpfr is refence counted, we inherit the assignment
        // operator and the copy constructor from Handle_for.
#ifdef CGAL_GMPFR_NO_REFCOUNT
        Gmpfr& operator=(const Gmpfr &a){
                mpfr_set_prec(fr(),a.get_precision());
                mpfr_set(fr(),a.fr(),mpfr_get_default_rounding_mode());
                return *this;
        }

        Gmpfr(const Gmpfr &a){
                mpfr_init2(fr(),a.get_precision());
                mpfr_set(fr(),a.fr(),MPFR_RNDN);
        }
#endif

        Gmpfr(const Gmpfr &a,
              std::float_round_style r,
              Gmpfr::Precision_type p=Gmpfr::get_default_precision()){
                CGAL_assertion(p>=MPFR_PREC_MIN&&p<=MPFR_PREC_MAX);
#ifndef CGAL_GMPFR_NO_REFCOUNT
                if(p==a.get_precision()){
                        Gmpfr temp(a);
                        // we use dont_clear_on_destruction because the
                        // mpfr_t pointed to by fr() was never initialized
                        dont_clear_on_destruction();
                        swap(temp);
                }else
#endif
                {
                        mpfr_init2(fr(),p);
                        mpfr_set(fr(),a.fr(),_gmp_rnd(r));
                }
        }

        Gmpfr(const Gmpfr &a,Gmpfr::Precision_type p){
                CGAL_assertion(p>=MPFR_PREC_MIN&&p<=MPFR_PREC_MAX);
#ifndef CGAL_GMPFR_NO_REFCOUNT
                if(p==a.get_precision()){
                        Gmpfr temp(a);
                        // we use dont_clear_on_destruction because the
                        // mpfr_t pointed to by fr() was never initialized
                        dont_clear_on_destruction();
                        swap(temp);
                }else
#endif
                {
                        mpfr_init2(fr(),p);
                        mpfr_set(fr(),a.fr(),mpfr_get_default_rounding_mode());
                }
        }

        // default rounding mode

        static std::float_round_style get_default_rndmode();
        static std::float_round_style
                set_default_rndmode(std::float_round_style);

        // default precision

        static Gmpfr::Precision_type get_default_precision();
        static Gmpfr::Precision_type
                set_default_precision(Gmpfr::Precision_type);

        // precision of a single Gmpfr object

        Gmpfr::Precision_type get_precision()const;
        Gmpfr round(Gmpfr::Precision_type,std::float_round_style)const;

        // mpfr global inexact flags

        static void clear_flags();
        static bool underflow_flag();
        static bool overflow_flag();
        static bool nan_flag();
        static bool inex_flag();
        static bool erange_flag();

        // arithmetics

        Gmpfr operator+()const;
        Gmpfr operator-()const;

#define CGAL_GMPFR_DECLARE_OPERATORS(_type) \
        Gmpfr& operator+=(_type); \
        Gmpfr& operator-=(_type); \
        Gmpfr& operator*=(_type); \
        Gmpfr& operator/=(_type);

        CGAL_GMPFR_DECLARE_OPERATORS(const Gmpfr&)
        Gmpfr& operator%=(const Gmpfr&);
        CGAL_GMPFR_DECLARE_OPERATORS(long)
        CGAL_GMPFR_DECLARE_OPERATORS(unsigned long)
        CGAL_GMPFR_DECLARE_OPERATORS(int)
        CGAL_GMPFR_DECLARE_OPERATORS(const Gmpz&)

#undef CGAL_GMPFR_DECLARE_OPERATORS

#define CGAL_GMPFR_DECLARE_STATIC_FUNCTION(_f,_t1,_t2) \
        static Gmpfr _f (_t1, \
                         _t2, \
                         std::float_round_style=Gmpfr::get_default_rndmode()); \
        static Gmpfr _f (_t1, \
                         _t2, \
                         Gmpfr::Precision_type, \
                         std::float_round_style=Gmpfr::get_default_rndmode());

#define CGAL_GMPFR_DECLARE_STATIC_FUNCTIONS(_t2) \
        CGAL_GMPFR_DECLARE_STATIC_FUNCTION(add,const Gmpfr&,_t2) \
        CGAL_GMPFR_DECLARE_STATIC_FUNCTION(sub,const Gmpfr&,_t2) \
        CGAL_GMPFR_DECLARE_STATIC_FUNCTION(mul,const Gmpfr&,_t2) \
        CGAL_GMPFR_DECLARE_STATIC_FUNCTION(div,const Gmpfr&,_t2)

        CGAL_GMPFR_DECLARE_STATIC_FUNCTIONS(const Gmpfr&)
        CGAL_GMPFR_DECLARE_STATIC_FUNCTIONS(long)
        CGAL_GMPFR_DECLARE_STATIC_FUNCTIONS(unsigned long)
        CGAL_GMPFR_DECLARE_STATIC_FUNCTIONS(int)
        CGAL_GMPFR_DECLARE_STATIC_FUNCTIONS(const Gmpz&)

#undef CGAL_GMPFR_DECLARE_STATIC_FUNCTION
#undef CGAL_GMPFR_DECLARE_STATIC_FUNCTIONS

#define CGAL_GMPFR_DECLARE_NO_ARGUMENT_FUNCTION(_f) \
        Gmpfr _f (std::float_round_style=Gmpfr::get_default_rndmode()) const; \
        Gmpfr _f (Gmpfr::Precision_type,\
                  std::float_round_style=Gmpfr::get_default_rndmode()) const;

        CGAL_GMPFR_DECLARE_NO_ARGUMENT_FUNCTION(abs)
        CGAL_GMPFR_DECLARE_NO_ARGUMENT_FUNCTION(sqrt)
        CGAL_GMPFR_DECLARE_NO_ARGUMENT_FUNCTION(cbrt)
        Gmpfr kthroot
                (int,std::float_round_style=Gmpfr::get_default_rndmode()) const;
        Gmpfr kthroot
                (int,
                 Gmpfr::Precision_type,
                 std::float_round_style=Gmpfr::get_default_rndmode()) const;
        CGAL_GMPFR_DECLARE_NO_ARGUMENT_FUNCTION(square)

#undef CGAL_GMPFR_DECLARE_NO_ARGUMENT_FUNCTION

        // comparison and query functions

        bool is_zero()const;
        bool is_one()const;
        bool is_nan()const;
        bool is_inf()const;
        bool is_number()const;
        Sign sign()const;
        bool is_square()const;
        bool is_square(Gmpfr&)const;
        Comparison_result compare(const Gmpfr&)const;

        // conversion functions

        double to_double(std::float_round_style=Gmpfr::get_default_rndmode())
                const;
        std::pair<double,double> to_interval()const;
        std::pair<double,long> to_double_exp(std::float_round_style=
                                             Gmpfr::get_default_rndmode())const;
        std::pair<std::pair<double,double>,long> to_interval_exp()const;
        std::pair<Gmpz,long> to_integer_exp()const;
};




// --------------
// implementation
// --------------

// default rounding mode, handled by mpfr
inline
std::float_round_style Gmpfr::get_default_rndmode(){
        return _cgal_rnd(mpfr_get_default_rounding_mode());
}

inline
std::float_round_style
Gmpfr::set_default_rndmode(std::float_round_style rnd_mode){
        std::float_round_style old_rnd_mode=Gmpfr::get_default_rndmode();
        mpfr_set_default_rounding_mode(_gmp_rnd(rnd_mode));
        return old_rnd_mode;
}

// default precision, handled by mpfr
inline
Gmpfr::Precision_type Gmpfr::get_default_precision(){
        return static_cast<Gmpfr::Precision_type>(mpfr_get_default_prec());
}

inline
Gmpfr::Precision_type Gmpfr::set_default_precision(Gmpfr::Precision_type prec){
        Gmpfr::Precision_type old_prec=Gmpfr::get_default_precision();
        CGAL_assertion(prec>=MPFR_PREC_MIN&&prec<=MPFR_PREC_MAX);
        mpfr_set_default_prec(prec);
        return old_prec;
}

// precision of a single Gmpfr object

inline
Gmpfr::Precision_type Gmpfr::get_precision()const{
        return mpfr_get_prec(fr());
}

inline
Gmpfr Gmpfr::round(Gmpfr::Precision_type p,std::float_round_style r)const{
        CGAL_assertion(p>=MPFR_PREC_MIN&&p<=MPFR_PREC_MAX);
        return Gmpfr(*this,r,p);
}

// mpfr global inexact flags

inline
void Gmpfr::clear_flags(){
        mpfr_clear_flags();
}

inline
bool Gmpfr::underflow_flag(){
        return mpfr_underflow_p()!=0;
}

inline
bool Gmpfr::overflow_flag(){
        return mpfr_overflow_p()!=0;
}

inline
bool Gmpfr::nan_flag(){
        return mpfr_nanflag_p()!=0;
}

inline
bool Gmpfr::inex_flag(){
        return mpfr_inexflag_p()!=0;
}

inline
bool Gmpfr::erange_flag(){
        return mpfr_erangeflag_p()!=0;
}

// arithmetics

inline
Gmpfr Gmpfr::operator+()const{
        return(*this);
}

inline
Gmpfr Gmpfr::operator-()const{
        Gmpfr result(0,get_precision());
        mpfr_neg(result.fr(),fr(),MPFR_RNDN);
        return result;
}

// CGAL_GMPFR_MEMBER_PREC returns the precision to be used to operate between
// *this and a number of another type or class. Currently, the maximum of
// *this' precision and the default precision is returned.
#define CGAL_GMPFR_MEMBER_PREC() \
        (get_precision()>Gmpfr::get_default_precision()? \
         get_precision(): \
         Gmpfr::get_default_precision())

// CGAL_GMPFR_MEMBER_PREC_2 returns the precision for the operation between Gmpfr
// objects *this and _b. Currently, it is the maximum of the precisions of
// *this and _b and the default precision.
// TODO: maybe we can rewrite this define optimally, maybe with an inline
#define CGAL_GMPFR_MEMBER_PREC_2(_b) \
        ( get_precision() >= mpfr_get_prec(_b.fr()) ? \
                ( get_precision()>(Gmpfr::get_default_precision())? \
                        get_precision():(Gmpfr::get_default_precision())): \
                ( mpfr_get_prec(_b.fr())>(Gmpfr::get_default_precision())? \
                        mpfr_get_prec(_b.fr()): \
                        (Gmpfr::get_default_precision())) \
        )

// CGAL_GMPFR_OBJECT_BINARY_OPERATOR defines an overloaded binary operator of
// the Gmpfr class, where the second parameter of the operator is an
// object. It behaves differently when the Gmpfr class is reference-counted
// or not.
#ifdef CGAL_GMPFR_NO_REFCOUNT
#define CGAL_GMPFR_OBJECT_BINARY_OPERATOR(_op,_class,_member,_fun) \
        inline \
        Gmpfr& Gmpfr::_op(const _class &b){ \
                if(get_precision()>=Gmpfr::get_default_precision()) { \
                        _fun(fr(), \
                             fr(), \
                             b._member, \
                             mpfr_get_default_rounding_mode()); \
                }else{ \
                        Gmpfr _temp(0,Gmpfr::get_default_precision()); \
                        _fun(_temp.fr(), \
                             fr(), \
                             b._member, \
                             mpfr_get_default_rounding_mode()); \
                        mpfr_swap(_temp.fr(),fr()); \
                } \
                return *this; \
        }
#else
#define CGAL_GMPFR_OBJECT_BINARY_OPERATOR(_op,_class,_member,_fun) \
        inline \
        Gmpfr& Gmpfr::_op(const _class &b){ \
                if(unique()){ \
                        if(get_precision()>Gmpfr::get_default_precision()) { \
                                _fun(fr(), \
                                     fr(), \
                                     b._member, \
                                     mpfr_get_default_rounding_mode()); \
                        }else{ \
                                Gmpfr _temp(0,Gmpfr::get_default_precision()); \
                                _fun(_temp.fr(), \
                                     fr(), \
                                     b._member, \
                                     mpfr_get_default_rounding_mode()); \
                                swap(_temp); \
                        } \
                }else{ \
                        Gmpfr result(0,CGAL_GMPFR_MEMBER_PREC()); \
                        _fun(result.fr(), \
                             fr(), \
                             b._member, \
                             mpfr_get_default_rounding_mode()); \
                        swap(result); \
                } \
                return *this; \
        }
#endif

// CGAL_GMPFR_GMPFR_BINARY_OPERATOR is analogous to
// CGAL_GMPFR_OBJECT_BINARY_OPERATOR, and it is used when the second operand is
// another Gmpfr. The difference is the computation of the operation
// precision.
#ifdef CGAL_GMPFR_NO_REFCOUNT
#define CGAL_GMPFR_GMPFR_BINARY_OPERATOR(_op,_fun) \
        inline \
        Gmpfr& Gmpfr::_op(const Gmpfr &b){ \
                Gmpfr::Precision_type _p=CGAL_GMPFR_MEMBER_PREC_2(b); \
                if(_p==get_precision()) { \
                        _fun(fr(), \
                             fr(), \
                             b.fr(), \
                             mpfr_get_default_rounding_mode()); \
                }else{ \
                        Gmpfr _temp(0,_p); \
                        _fun(_temp.fr(), \
                             fr(), \
                             b.fr(), \
                             mpfr_get_default_rounding_mode()); \
                        mpfr_swap(_temp.fr(),fr()); \
                } \
                return *this; \
        }
#else
#define CGAL_GMPFR_GMPFR_BINARY_OPERATOR(_op,_fun) \
        inline \
        Gmpfr& Gmpfr::_op(const Gmpfr &b){ \
                Gmpfr::Precision_type _p=CGAL_GMPFR_MEMBER_PREC_2(b); \
                if(unique()&&(_p==get_precision())){ \
                        _fun(fr(), \
                             fr(), \
                             b.fr(), \
                             mpfr_get_default_rounding_mode()); \
                }else{ \
                        Gmpfr result(0,_p); \
                        _fun(result.fr(), \
                             fr(), \
                             b.fr(), \
                             mpfr_get_default_rounding_mode()); \
                        swap(result); \
                } \
                return *this; \
        }
#endif

// CGAL_GMPFR_TYPE_BINARY_OPERATOR is analogous to the
// CGAL_GMPFR_OBJECT_BINARY_OPERATOR, where the second parameter is a type
// instead of an object.
#ifdef CGAL_GMPFR_NO_REFCOUNT
#define CGAL_GMPFR_TYPE_BINARY_OPERATOR(_op,_type,_fun) \
        inline \
        Gmpfr& Gmpfr::_op(_type x){ \
                if(get_precision()>=Gmpfr::get_default_precision()) { \
                        _fun(fr(),fr(),x,mpfr_get_default_rounding_mode()); \
                }else{ \
                        Gmpfr _temp(0,Gmpfr::get_default_precision()); \
                        _fun(_temp.fr(), \
                             fr(), \
                             x, \
                             mpfr_get_default_rounding_mode()); \
                        mpfr_swap(_temp.fr(),fr()); \
                } \
                return *this; \
        }
#else
#define CGAL_GMPFR_TYPE_BINARY_OPERATOR(_op,_type,_fun) \
        inline \
        Gmpfr& Gmpfr::_op(_type x){ \
                if(unique()){ \
                        if(get_precision()>Gmpfr::get_default_precision()) { \
                                _fun(fr(), \
                                     fr(), \
                                     x, \
                                     mpfr_get_default_rounding_mode()); \
                        }else{ \
                                Gmpfr _temp(0,Gmpfr::get_default_precision()); \
                                _fun(_temp.fr(), \
                                     fr(), \
                                     x, \
                                     mpfr_get_default_rounding_mode()); \
                                swap(_temp); \
                        } \
                }else{ \
                        Gmpfr result(0,CGAL_GMPFR_MEMBER_PREC()); \
                        _fun(result.fr(), \
                             fr(), \
                             x, \
                             mpfr_get_default_rounding_mode()); \
                        swap(result); \
                } \
                return *this; \
        }
#endif

CGAL_GMPFR_GMPFR_BINARY_OPERATOR(operator+=,mpfr_add)
CGAL_GMPFR_GMPFR_BINARY_OPERATOR(operator-=,mpfr_sub)
CGAL_GMPFR_GMPFR_BINARY_OPERATOR(operator*=,mpfr_mul)
CGAL_GMPFR_GMPFR_BINARY_OPERATOR(operator/=,mpfr_div)
#if(defined(MPFR_VERSION)&&(MPFR_VERSION>=MPFR_VERSION_NUM(2,3,0)))
CGAL_GMPFR_GMPFR_BINARY_OPERATOR(operator%=,mpfr_remainder)
#else
//#  warning "Gmpfr::operator%= is optimized in MPFR 2.3.0."
inline
Gmpfr& Gmpfr::operator%=(const Gmpfr &b){
        Gmpfr::Precision_type _p=CGAL_GMPFR_MEMBER_PREC_2(b);
        Gmpfr result(*this,_p);
        result/=b;
        mpfr_trunc(result.fr(),result.fr());
        result*=b;
        result-=*this;
        mpfr_neg(result.fr(),result.fr(),MPFR_RNDN);
#  ifdef CGAL_GMPFR_NO_REFCOUNT
        mpfr_swap(result.fr(),fr());
#  else
        if(unique())
                mpfr_swap(result.fr(),fr());
        else
                swap(result);
#  endif
        return *this;
}
#endif

CGAL_GMPFR_TYPE_BINARY_OPERATOR(operator+=,long,mpfr_add_si)
CGAL_GMPFR_TYPE_BINARY_OPERATOR(operator-=,long,mpfr_sub_si)
CGAL_GMPFR_TYPE_BINARY_OPERATOR(operator*=,long,mpfr_mul_si)
CGAL_GMPFR_TYPE_BINARY_OPERATOR(operator/=,long,mpfr_div_si)

CGAL_GMPFR_TYPE_BINARY_OPERATOR(operator+=,unsigned long,mpfr_add_ui)
CGAL_GMPFR_TYPE_BINARY_OPERATOR(operator-=,unsigned long,mpfr_sub_ui)
CGAL_GMPFR_TYPE_BINARY_OPERATOR(operator*=,unsigned long,mpfr_mul_ui)
CGAL_GMPFR_TYPE_BINARY_OPERATOR(operator/=,unsigned long,mpfr_div_ui)

CGAL_GMPFR_TYPE_BINARY_OPERATOR(operator+=,int,mpfr_add_si)
CGAL_GMPFR_TYPE_BINARY_OPERATOR(operator-=,int,mpfr_sub_si)
CGAL_GMPFR_TYPE_BINARY_OPERATOR(operator*=,int,mpfr_mul_si)
CGAL_GMPFR_TYPE_BINARY_OPERATOR(operator/=,int,mpfr_div_si)

CGAL_GMPFR_OBJECT_BINARY_OPERATOR(operator+=,Gmpz,mpz(),mpfr_add_z)
CGAL_GMPFR_OBJECT_BINARY_OPERATOR(operator-=,Gmpz,mpz(),mpfr_sub_z)
CGAL_GMPFR_OBJECT_BINARY_OPERATOR(operator*=,Gmpz,mpz(),mpfr_mul_z)
CGAL_GMPFR_OBJECT_BINARY_OPERATOR(operator/=,Gmpz,mpz(),mpfr_div_z)

#undef CGAL_GMPFR_OBJECT_BINARY_OPERATOR
#undef CGAL_GMPFR_GMPFR_BINARY_OPERATOR
#undef CGAL_GMPFR_TYPE_BINARY_OPERATOR

// the static arithmetic functions are defined in a separate file
#include <CGAL/GMP/Gmpfr_type_static.h>

#define CGAL_GMPFR_ARITHMETIC_FUNCTION(_name,_fun) \
        inline \
        Gmpfr Gmpfr::_name (std::float_round_style r)const{ \
                Gmpfr result(0,CGAL_GMPFR_MEMBER_PREC()); \
                _fun(result.fr(),fr(),_gmp_rnd(r)); \
                return result; \
        } \
        inline \
        Gmpfr Gmpfr::_name (Gmpfr::Precision_type p, \
                            std::float_round_style r)const{ \
                CGAL_assertion(p>=MPFR_PREC_MIN&&p<=MPFR_PREC_MAX); \
                Gmpfr result(0,p); \
                _fun(result.fr(),fr(),_gmp_rnd(r)); \
                return result; \
        }

CGAL_GMPFR_ARITHMETIC_FUNCTION(abs,mpfr_abs)
CGAL_GMPFR_ARITHMETIC_FUNCTION(sqrt,mpfr_sqrt)
CGAL_GMPFR_ARITHMETIC_FUNCTION(cbrt,mpfr_cbrt)

inline
Gmpfr Gmpfr::kthroot(int k,std::float_round_style r)const{
        Gmpfr result(0,CGAL_GMPFR_MEMBER_PREC());
        #if(MPFR_VERSION_MAJOR < 4)
            mpfr_root(result.fr(),fr(),k,_gmp_rnd(r));
        #else
            mpfr_rootn_ui(result.fr(),fr(),k,_gmp_rnd(r));
        #endif
        return result;
}

inline
Gmpfr Gmpfr::kthroot(int k,
                     Gmpfr::Precision_type p,
                     std::float_round_style r)const{
        CGAL_assertion(p>=MPFR_PREC_MIN&&p<=MPFR_PREC_MAX);
        Gmpfr result(0,p);
        #if(MPFR_VERSION_MAJOR < 4)
            mpfr_root(result.fr(),fr(),k,_gmp_rnd(r));
        #else
            mpfr_rootn_ui(result.fr(),fr(),k,_gmp_rnd(r));
        #endif
        return result;
}

CGAL_GMPFR_ARITHMETIC_FUNCTION(square,mpfr_sqr)

#undef CGAL_GMPFR_ARITHMETIC_FUNCTION
#undef CGAL_GMPFR_MEMBER_PREC
#undef CGAL_GMPFR_MEMBER_PREC_2

// comparison and query functions

inline
bool Gmpfr::is_zero()const{
        return mpfr_zero_p(fr())!=0;
}

inline
bool Gmpfr::is_one()const{
        return mpfr_cmp_ui(fr(),1)==0;
}

inline
bool Gmpfr::is_nan()const{
        return mpfr_nan_p(fr())!=0;
}

inline
bool Gmpfr::is_inf()const{
        return mpfr_inf_p(fr())!=0;
}

inline
bool Gmpfr::is_number()const{
        return mpfr_number_p(fr())!=0;
}

inline
Sign Gmpfr::sign()const{
        int s=mpfr_sgn(fr());
        return(s==0?ZERO:(s>0?POSITIVE:NEGATIVE));
}

inline
bool Gmpfr::is_square()const{
        Sign s=sign();
        if(s==NEGATIVE)
                return false;
        if(s==ZERO)
                return true;
        std::pair<Gmpz,long> r=Gmpfr::to_integer_exp();
        if(r.second%2)
                r.first=r.first*2;
        return mpz_perfect_square_p(r.first.mpz())!=0;
}

inline
bool Gmpfr::is_square(Gmpfr &y)const{
        bool ret=is_square();
        if(ret)
                y=sqrt();
        return ret;
}

inline
Comparison_result Gmpfr::compare(const Gmpfr& b)const{
        int c=mpfr_cmp(fr(),b.fr());
        return(c?(c>0?LARGER:SMALLER):EQUAL);
}

// conversion functions

inline
double Gmpfr::to_double(std::float_round_style r)const{
        return mpfr_get_d(fr(),_gmp_rnd(r));
}

inline
std::pair<double,double>Gmpfr::to_interval()const{
        return std::make_pair(
                        mpfr_get_d(fr(),MPFR_RNDD),
                        mpfr_get_d(fr(),MPFR_RNDU));
}

inline
std::pair<double,long> Gmpfr::to_double_exp(std::float_round_style r)const{
        long e;
        double d=mpfr_get_d_2exp(&e,fr(),_gmp_rnd(r));
        return std::make_pair(d,e);
}

inline
std::pair<std::pair<double,double>,long> Gmpfr::to_interval_exp()const{
        long e1,e2;
        double d_low=mpfr_get_d_2exp(&e1,fr(),MPFR_RNDD);
        double d_upp=mpfr_get_d_2exp(&e2,fr(),MPFR_RNDU);
        CGAL_assertion(e1==e2);
        return std::make_pair(std::make_pair(d_low,d_upp),e1);
}

inline
std::pair<Gmpz,long> Gmpfr::to_integer_exp()const{
        if(this->is_zero())
                return std::make_pair(Gmpz(0),long(0));

        Gmpz z;
        long e=CGAL_GMPFR_GET_Z_2EXP(z.mpz(),this->fr());

        long zeros=mpz_scan1(z.mpz(),0);
        CGAL_assertion(z==(z>>zeros)<<zeros);
        z>>=zeros;
        CGAL_assertion(z%2!=0);
        e+=zeros;

        CGAL_postcondition_code(if(e>=0))
        CGAL_postcondition(
              (*this)==(Gmpfr(z,(mpfr_prec_t)z.bit_size())*CGAL::ipower(Gmpfr(2),e)));
        CGAL_postcondition_code(else)
        CGAL_postcondition(((*this)*(Gmpz(1)<<(-e)))==z);

        return std::make_pair(z,e);
}


// input/output

// This function was based on the Gmpq's one. It reads a number in the form
// MeE, where M and E are integers. The read number is M.2^E. The number
// may contain spaces between integers and the 'e', but not in the middle
// of the numbers.
inline
std::istream& operator>>(std::istream& is,Gmpfr &f){
        std::istream::int_type c;
        std::ios::fmtflags old_flags = is.flags();

        is.unsetf(std::ios::skipws);
        internal::eat_white_space(is);

        // 1. read the mantissa, it starts in +, - or a digit and ends in e
        Gmpz mant(0);           // the mantissa of the number
        Gmpz exp(0);            // the exponent of the number
        bool neg_mant=false;    // true iff the mantissa is negative
        bool neg_exp=false;     // true iff the exponent is negative
        c=is.peek();
        switch(c){
                case '-':
                        neg_mant=true;
                        is.get();
                        internal::eat_white_space(is);
                        break;
                case '+':
                        is.get();
                        internal::eat_white_space(is);
                        break;
                case 'n':       // this is NaN
                        is.get();
                        if(is.get()=='a'&&is.get()=='n'){
                                f=Gmpfr();
                                return is;
                        }
                        else
                                goto invalid_number;
                default:
                        if(c!='i'&&(c<'0'||c>'9')){     // invalid number
                                invalid_number:
                                is.setstate(std::ios_base::failbit);
                                is.flags(old_flags);
                                return is;
                        }
        }

        // at this point, we have the sign of the number and we are ready
        // to read the mantissa
        c=is.get();
        if(c=='i'){     // infinity comes
                if(is.get()=='n'&&is.get()=='f'){
                        f=Gmpfr();
                        mpfr_set_inf(f.fr(),neg_mant?-1:1);
                        return is;
                }
                else
                        goto invalid_number;
        }

        while(c>='0'&&c<='9'){
                mant=10*mant+(c-'0');
                c=is.get();
        }

        // set the correct sign of the mantissa
        if(neg_mant)
                mant=-mant;

        is.putback(static_cast<std::istream::char_type>(c));
        internal::eat_white_space(is);

        switch(c=is.get()){
                case 'e':
                        break;
                default:
                        is.setstate(std::ios_base::failbit);
                        is.flags(old_flags);
                        return is;
        }
        c=is.peek();
        switch(c){
                case '-':
                        neg_exp=true;
                        is.get();
                        internal::eat_white_space(is);
                        break;
                case '+':
                        is.get();
                        internal::eat_white_space(is);
                        break;
                default:
                        if(c<'0'||c>'9')
                                goto invalid_number;
        }
        internal::eat_white_space(is);
        while((c=is.get())>='0'&&c<='9')
                exp=10*exp+(c-'0');
        is.putback(static_cast<std::istream::char_type>(c));
        if(exp.bit_size()>8*sizeof(mpfr_exp_t))
                mpfr_set_erangeflag();

        // we have now both exponent and mantissa
        f=Gmpfr(mant,
                static_cast<Gmpfr::Precision_type>(
                        mant.bit_size()>MPFR_PREC_MIN?
                        mant.bit_size():
                        MPFR_PREC_MIN));
        if(neg_exp)
                mpfr_div_2ui(f.fr(),f.fr(),mpz_get_ui(exp.mpz()),MPFR_RNDN);
        else
                mpfr_mul_2ui(f.fr(),f.fr(),mpz_get_ui(exp.mpz()),MPFR_RNDN);

        // this expensive assertion checks that we didn't lose bits when
        // multiplying or dividing by 2^exp
        CGAL_expensive_assertion_code( \
                Gmpfr g(0,static_cast<Gmpfr::Precision_type>( \
                                MPFR_PREC_MIN<mant.bit_size()? \
                                mant.bit_size(): \
                                MPFR_PREC_MIN)); \
                if(neg_exp) \
                        mpfr_div_2ui(g.fr(), \
                                     f.fr(), \
                                     mpz_get_ui(exp.mpz()), \
                                     MPFR_RNDN); \
                else \
                        mpfr_mul_2ui(g.fr(), \
                                     f.fr(), \
                                     mpz_get_ui(exp.mpz()), \
                                     MPFR_RNDN); \
        )
        CGAL_expensive_assertion(g==mant);

        return is;
}

inline
std::ostream& operator<<(std::ostream& os,const Gmpfr &a){
        if(a.is_nan())
                return os<<"nan";
        if(a.is_inf())
                return os<<(a<0?"-inf":"+inf");
        // The rest of the function was written by George Tzoumas.
        if (!is_pretty(os)) {
                std::pair<Gmpz,long> ie=a.to_integer_exp();
                os << ie.first << 'e' << ie.second;
                return os;
        } else {
                // human-readable format
                mpfr_exp_t expptr;
                char *str = mpfr_get_str(nullptr, &expptr, 10, 0, a.fr(),
                                mpfr_get_default_rounding_mode());
                if (str == nullptr) return os << "@err@";
                std::string s(str);
                mpfr_free_str(str);
                int i = 0;
                size_t n = s.length();
                size_t k = 0;
                while (k < n && s[n-k-1] == '0') k++; // count trailing zeros
                if (k == n) return os << "0";
                else if (k) {
                        s.erase(n-k, k);  // remove trailing zeros
                        n = s.length();
                }
                bool exp = false;
                if(s[0] == '-') { os << "-"; i++; n--; } // sign
                if (expptr < -5) {              // .125e-99
                        s.insert(i, 1, '.'); exp = true;
                } else if (expptr < 0) {
                        s.insert(i, -expptr, '0');  // .00000125 -- .0125
                        s.insert(i, 1, '.');
                // The following cast of expptr is done for avoiding some
                // compiler warnings. The cast is exact, because we know
                // expptr is not negative here.
                } else if ((size_t)expptr < n) {        // .125 -- 12.5
                        s.insert(i+expptr, 1, '.');
                } else if (expptr - n <= 5) {   // 125 -- 12500000
                        s.append(expptr - n, '0');
                } else {                        // .125e99
                        s.insert(i, 1, '.'); exp = true;
                }
                os << s.substr(i);
                if (exp) os << "e" << expptr;
                return os;
        }
}

// comparisons

inline
bool operator<(const Gmpfr &a,const Gmpfr &b){
        return mpfr_less_p(a.fr(),b.fr())!=0;
}

inline
bool operator==(const Gmpfr &a,const Gmpfr &b){
        return mpfr_equal_p(a.fr(),b.fr())!=0;
}

inline
bool operator<(const Gmpfr &a,long b){
        return(mpfr_cmp_si(a.fr(),b)<0);
}

inline
bool operator>(const Gmpfr &a,long b){
        return(mpfr_cmp_si(a.fr(),b)>0);
}

inline
bool operator==(const Gmpfr &a,long b){
        return !mpfr_cmp_si(a.fr(),b);
}

inline
bool operator<(const Gmpfr &a,unsigned long b){
        return(mpfr_cmp_ui(a.fr(),b)<0);
}

inline
bool operator>(const Gmpfr &a,unsigned long b){
        return(mpfr_cmp_ui(a.fr(),b)>0);
}

inline
bool operator==(const Gmpfr &a,unsigned long b){
        return !mpfr_cmp_ui(a.fr(),b);
}

inline
bool operator<(const Gmpfr &a,int b){
        return(mpfr_cmp_si(a.fr(),b)<0);
}

inline
bool operator>(const Gmpfr &a,int b){
        return(mpfr_cmp_si(a.fr(),b)>0);
}

inline
bool operator==(const Gmpfr &a,int b){
        return !mpfr_cmp_si(a.fr(),b);
}

inline
bool operator<(const Gmpfr &a,double b){
        return(mpfr_cmp_d(a.fr(),b)<0);
}

inline
bool operator>(const Gmpfr &a,double b){
        return(mpfr_cmp_d(a.fr(),b)>0);
}

inline
bool operator==(const Gmpfr &a,double b){
        return !mpfr_cmp_d(a.fr(),b);
}

// See the comment about mpfr_set_ld and MSVC++, above.
#ifdef _MSC_VER
inline
bool operator<(const Gmpfr &a,long double b){
        return(mpfr_cmp_d(a.fr(),static_cast<double>(b))<0);
}

inline
bool operator>(const Gmpfr &a,long double b){
        return(mpfr_cmp_d(a.fr(),static_cast<double>(b))>0);
}

inline
bool operator==(const Gmpfr &a,long double b){
        return !mpfr_cmp_d(a.fr(),static_cast<double>(b));
}
#else
inline
bool operator<(const Gmpfr &a,long double b){
        return(mpfr_cmp_ld(a.fr(),b)<0);
}

inline
bool operator>(const Gmpfr &a,long double b){
        return(mpfr_cmp_ld(a.fr(),b)>0);
}

inline
bool operator==(const Gmpfr &a,long double b){
        return !mpfr_cmp_ld(a.fr(),b);
}
#endif

inline
bool operator<(const Gmpfr &a,const Gmpz &b){
        return(mpfr_cmp_z(a.fr(),b.mpz())<0);
}

inline
bool operator>(const Gmpfr &a,const Gmpz &b){
        return(mpfr_cmp_z(a.fr(),b.mpz())>0);
}

inline
bool operator==(const Gmpfr &a,const Gmpz &b){
        return !mpfr_cmp_z(a.fr(),b.mpz());
}

inline
Gmpfr min BOOST_PREVENT_MACRO_SUBSTITUTION(const Gmpfr& x,const Gmpfr& y){
        return (x<=y)?x:y;
}

inline
Gmpfr max BOOST_PREVENT_MACRO_SUBSTITUTION(const Gmpfr& x,const Gmpfr& y){
        return (x>=y)?x:y;
}

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif  // CGAL_GMPFR_TYPE_H
