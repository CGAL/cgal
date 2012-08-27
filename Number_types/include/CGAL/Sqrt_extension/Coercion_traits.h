// Copyright (c) 2006-2008 Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Michael Hemmer   <hemmer@mpi-inf.mpg.de>


#ifndef CGAL_SQRT_EXTENSION_COERCION_TRAITS_H
#define CGAL_SQRT_EXTENSION_COERCION_TRAITS_H

#include <CGAL/basic.h>

#define CGAL_int(T)    typename First_if_different<int,    T>::Type

namespace CGAL {

/////////// COERCION_TRAITS BEGIN

// <EXT,int> and vice versa 
template <class COEFF, class ROOT, class ACDE_TAG,class FP_TAG>
struct Coercion_traits_for_level<Sqrt_extension<COEFF, ROOT, ACDE_TAG,FP_TAG>,CGAL_int(COEFF),CTL_SQRT_EXT>{
public:
  typedef Sqrt_extension<COEFF,ROOT, ACDE_TAG,FP_TAG> Type;
  typedef CGAL::Tag_true Are_explicit_interoperable;
  typedef CGAL::Tag_true Are_implicit_interoperable;
  struct Cast{
    typedef Type result_type;
    Type operator()(const Type& x) const { return x;}
    Type operator()(int x) const { return Type(x);}
  };
};

template <class COEFF, class ROOT, class ACDE_TAG,class FP_TAG>
struct Coercion_traits_for_level<CGAL_int(COEFF), Sqrt_extension<COEFF, ROOT, ACDE_TAG,FP_TAG>,CTL_SQRT_EXT>
  : public Coercion_traits_for_level<Sqrt_extension<COEFF, ROOT, ACDE_TAG,FP_TAG>,CGAL_int(COEFF),CTL_SQRT_EXT>{};




template <class COEFF, class ROOT, class ACDE_TAG,class FP_TAG>
struct Coercion_traits_for_level<Sqrt_extension<COEFF, ROOT, ACDE_TAG,FP_TAG>,COEFF,CTL_SQRT_EXT>{
public:
  typedef Sqrt_extension<COEFF,ROOT, ACDE_TAG,FP_TAG> Type;
  typedef CGAL::Tag_true Are_explicit_interoperable;
  typedef CGAL::Tag_true Are_implicit_interoperable;
  struct Cast{
    typedef Type result_type;
    Type operator()(const Type& x) const { return x;}
    Type operator()(COEFF x) const { return Type(x);}
  };
};

template <class COEFF, class ROOT, class ACDE_TAG,class FP_TAG>
struct Coercion_traits_for_level<COEFF, Sqrt_extension<COEFF, ROOT, ACDE_TAG,FP_TAG>,CTL_SQRT_EXT>
  : public Coercion_traits_for_level<Sqrt_extension<COEFF, ROOT, ACDE_TAG,FP_TAG>,COEFF,CTL_SQRT_EXT>{};


// <EXT,EXT>
template <class A_coeff, class B_coeff, class ROOT, class ACDE_TAG, class FP_TAG>
struct Coercion_traits_for_level<Sqrt_extension<A_coeff, ROOT, ACDE_TAG,FP_TAG>,
                           Sqrt_extension<B_coeff, ROOT, ACDE_TAG,FP_TAG>,
                           CTL_SQRT_EXT>{
private:
    typedef Coercion_traits<A_coeff, B_coeff> CT;
    typedef Sqrt_extension<A_coeff,ROOT, ACDE_TAG,FP_TAG> A;
    typedef Sqrt_extension<B_coeff,ROOT, ACDE_TAG,FP_TAG> B;

public:
    typedef CGAL::Tag_true  Are_explicit_interoperable;
    typedef CGAL::Tag_false Are_implicit_interoperable;
    typedef Sqrt_extension<typename CT::Type, ROOT, ACDE_TAG,FP_TAG> Type;

    struct Cast{
    private:
        inline Type cast(const Type& x) const{ return x; }

        template <class T>
        inline Type cast(const T& x) const{
            typename CT::Cast cast;
            if (x.is_extended()) {
                return result_type(cast(x.a0()),cast(x.a1()),x.root());
            } else {
                return result_type(cast(x.a0()));
            }
        }
    public:
        typedef Type result_type;
        // this is in order to allow A and B only
        Type operator()(const A& x) const { return cast(x);}
        Type operator()(const B& x) const { return cast(x);}
    };
};


template <class COEFF, class ROOT_1, class ROOT_2, class ACDE_TAG, class FP_TAG>
struct Coercion_traits_for_level<Sqrt_extension<Sqrt_extension<COEFF,ROOT_1,ACDE_TAG,FP_TAG>,ROOT_2,ACDE_TAG,FP_TAG>,
                           Sqrt_extension<COEFF,ROOT_1,ACDE_TAG,FP_TAG>,
                           CTL_SQRT_EXT>{
private:
  typedef Sqrt_extension<Sqrt_extension<COEFF,ROOT_1,ACDE_TAG>, ROOT_2,ACDE_TAG,FP_TAG> A;
  typedef Sqrt_extension<COEFF,ROOT_1,ACDE_TAG,FP_TAG> B;
public:
    typedef CGAL::Tag_true Are_explicit_interoperable;
    typedef CGAL::Tag_true Are_implicit_interoperable;

    // Type = A
  typedef Sqrt_extension<Sqrt_extension<COEFF,ROOT_1,ACDE_TAG,FP_TAG>, ROOT_2,ACDE_TAG,FP_TAG> Type;
    struct Cast{
        typedef Type result_type;
        Type operator()(const A& x) const { return x;}
        Type operator()(const B& x) const { return Type(x);}
    };
};

template <class COEFF, class ROOT_1, class ROOT_2, class ACDE_TAG,class FP_TAG>
struct Coercion_traits_for_level
<
            Sqrt_extension<COEFF,ROOT_1,ACDE_TAG,FP_TAG>,
            Sqrt_extension<Sqrt_extension<COEFF,ROOT_1,ACDE_TAG,FP_TAG>, ROOT_2, ACDE_TAG,FP_TAG>
            ,CTL_SQRT_EXT>
    :public Coercion_traits_for_level
<
            Sqrt_extension<Sqrt_extension<COEFF,ROOT_1,ACDE_TAG,FP_TAG>, ROOT_2, ACDE_TAG,FP_TAG>,
            Sqrt_extension<COEFF,ROOT_1,ACDE_TAG,FP_TAG>
            ,CTL_SQRT_EXT>
{};


template <class COEFF, class ROOT_1,class ACDE_TAG,class FP_TAG>
struct Coercion_traits_for_level
<
            Sqrt_extension<Sqrt_extension<COEFF,ROOT_1,ACDE_TAG,FP_TAG>, ROOT_1, ACDE_TAG,FP_TAG>,
            Sqrt_extension<COEFF,ROOT_1,ACDE_TAG,FP_TAG>
            ,CTL_SQRT_EXT>{
private:
  typedef  Sqrt_extension<Sqrt_extension<COEFF,ROOT_1,ACDE_TAG,FP_TAG>, ROOT_1,ACDE_TAG,FP_TAG> A;
  typedef  Sqrt_extension<COEFF,ROOT_1,ACDE_TAG,FP_TAG> B;
public:
    typedef CGAL::Tag_true Are_explicit_interoperable;
    typedef CGAL::Tag_true Are_implicit_interoperable;

  typedef  Sqrt_extension<Sqrt_extension<COEFF,ROOT_1,ACDE_TAG,FP_TAG>, ROOT_1, ACDE_TAG,FP_TAG> Type;
    struct Cast{
        typedef Type result_type;
        Type operator()(const A& x) const { return x;}
        Type operator()(const B& x) const { return Type(x);}
    };
};

template <class COEFF, class ROOT_1, class ACDE_TAG, class FP_TAG>
struct Coercion_traits_for_level
<
            Sqrt_extension<COEFF,ROOT_1,ACDE_TAG,FP_TAG>,
            Sqrt_extension<Sqrt_extension<COEFF,ROOT_1,ACDE_TAG,FP_TAG>, ROOT_1,ACDE_TAG,FP_TAG>
            ,CTL_SQRT_EXT>
    :public Coercion_traits_for_level
<
            Sqrt_extension<Sqrt_extension<COEFF,ROOT_1,ACDE_TAG,FP_TAG>, ROOT_1,ACDE_TAG,FP_TAG>,
            Sqrt_extension<COEFF,ROOT_1,ACDE_TAG,FP_TAG>
            ,CTL_SQRT_EXT>
{};


namespace INTERN_CT{
// Coercion_traits for Sqrt_extenison to FieldWithSqrt
template <class A, class B> struct CT_ext_to_fwsqrt;
// Coercion_traits for Sqrt_extenison not with FieldWithSqrt
template <class A, class B> struct CT_ext_not_to_fwsqrt;
} // namespace INTERN_CT


//<EXT,ANY>
template <class COEFF, class ROOT, class ACDE_TAG,class FP_TAG, class B>
struct Coercion_traits_for_level<Sqrt_extension<COEFF, ROOT, ACDE_TAG,FP_TAG>, B , CTL_SQRT_EXT>
:public ::boost::mpl::if_c<
             // if B is fwsqrt
              ::boost::is_base_and_derived<
                  Field_with_sqrt_tag,
typename Algebraic_structure_traits<B>::Algebraic_category >::value ||
              ::boost::is_same<
                  Field_with_sqrt_tag,
typename Algebraic_structure_traits<B>::Algebraic_category >::value
            ,
            //then take Intern::Coercion_traits for fwsqrt
            INTERN_CT::CT_ext_to_fwsqrt<Sqrt_extension<COEFF,ROOT, ACDE_TAG,FP_TAG>, B>
            ,
            //else take Intern::Coercion_traits not for fwsqrt
            INTERN_CT::CT_ext_not_to_fwsqrt< Sqrt_extension<COEFF,ROOT, ACDE_TAG,FP_TAG> ,B>
              >::type
{};

// <ANY,EXT>
template <class COEFF, class ROOT, class ACDE_TAG,class FP_TAG, class B>
struct Coercion_traits_for_level
<B,Sqrt_extension<COEFF, ROOT, ACDE_TAG,FP_TAG>,CTL_SQRT_EXT >
    :public Coercion_traits_for_level<Sqrt_extension<COEFF,ROOT, ACDE_TAG,FP_TAG>,B,CTL_SQRT_EXT>
{};

namespace INTERN_CT{
// EXT coercion with FieldWithSqrt
template <class COEFF, class ROOT, class ACDE_TAG,class FP_TAG, class FieldWithSqrt>
struct CT_ext_to_fwsqrt<Sqrt_extension<COEFF,ROOT,ACDE_TAG,FP_TAG>,
                                         FieldWithSqrt>{
private:
    typedef Sqrt_extension<COEFF,ROOT,ACDE_TAG,FP_TAG> A;
    typedef FieldWithSqrt B;
public:
    typedef CGAL::Tag_true  Are_explicit_interoperable;
    typedef CGAL::Tag_false Are_implicit_interoperable;

    typedef FieldWithSqrt Type;
    struct Cast{
        typedef Type result_type;
        Type operator()(const A& x) const {
            typedef Coercion_traits<COEFF,FieldWithSqrt> CT_coeff;
            typedef Coercion_traits<ROOT ,FieldWithSqrt> CT_root;
            typename CT_coeff::Cast coeff_cast;
            typename CT_root::Cast root_cast;
            if (x.is_extended()) {
                typename CGAL::Algebraic_structure_traits<
                typename CT_root::Type>::Sqrt sqrt;
                return // a0+a1*sqrt(root)
                    coeff_cast(x.a0())+
                    coeff_cast(x.a1())*
                    sqrt(root_cast(x.root()));
            } else {
                return coeff_cast(x.a0());
            }
        }
        Type operator()(const B& x) const { return x;}
    };
};

// EXT coercion not with FieldWithSqrt
template <class COEFF, class ROOT, class ACDE_TAG,class FP_TAG, class B_>
struct CT_ext_not_to_fwsqrt<Sqrt_extension<COEFF,ROOT,ACDE_TAG,FP_TAG>, B_>{
private:
    typedef Sqrt_extension<COEFF,ROOT,ACDE_TAG,FP_TAG> A;
    typedef B_ B;
    typedef Coercion_traits<COEFF,B> CT;
public:
    typedef CGAL::Tag_true  Are_explicit_interoperable;
    typedef CGAL::Tag_false Are_implicit_interoperable;
    typedef Sqrt_extension<typename CT::Type,ROOT,ACDE_TAG,FP_TAG> Type;
    struct Cast{
        typedef Type result_type;
        Type operator()(const A& x) const {
            typename CT::Cast cast;
            if (x.is_extended()) {
                return Type(cast(x.a0()),cast(x.a1()),x.root());
            } else {
                return Type(cast(x.a0()));
            }
        }
        Type operator()(const B& x) const {
            typename CT::Cast cast;
            return Type(cast(x));
        }
    };
};
} // namespace INTERN_CT

/////////// COERCION_TRAITS END

} //namespace CGAL

#undef CGAL_int

#endif
