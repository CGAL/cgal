// Copyright (c) 2014
// INRIA Saclay-Ile de France (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Marc Glisse

#ifndef CGAL_KERNEL_D_CARTESIAN_CONVERTER_H
#define CGAL_KERNEL_D_CARTESIAN_CONVERTER_H

#include <CGAL/basic.h>
#include <CGAL/tuple.h>
#include <CGAL/typeset.h>
#include <CGAL/Object.h>
#include <CGAL/Origin.h>
#include <CGAL/NT_converter.h>
#include <CGAL/NewKernel_d/functor_tags.h>
#include <CGAL/Kernel/mpl.h>
#include <CGAL/is_iterator.h>
#include <CGAL/transforming_iterator.h>
#include <boost/utility/enable_if.hpp>
#include <boost/mpl/if.hpp>
#include <CGAL/NewKernel_d/store_kernel.h>
#include <CGAL/NewKernel_d/Kernel_object_converter.h>

namespace CGAL {
namespace internal {
// Reverses order, but that shouldn't matter.
template<class K,class T> struct Map_taglist_to_typelist :
  Map_taglist_to_typelist<K,typename T::tail>::type
  ::template add<typename Get_type<K, typename T::head>::type>
{};
template<class K> struct Map_taglist_to_typelist<K,typeset<> > : typeset<> {};
}

template<class List = typeset<> >
struct Object_converter {
        typedef Object result_type;
        template<class F>
        result_type operator()(Object const& o, F const& f) const {
          typedef typename List::head H;
                if (H const* ptr = object_cast<H>(&o))
                        return make_object(f(*ptr));
                else
                        return Object_converter<typename List::tail>()(o,f);
        }
};
template<>
struct Object_converter <typeset<> > {
        typedef Object result_type;
        template<class F>
        result_type operator()(Object const&,F const&)const {
                CGAL_error_msg("Cartesiand_converter is unable to determine what is wrapped in the Object");
                return Object();
        }
};


        //TODO: special case when K1==K2 (or they are very close?)
template<class Final_, class K1, class K2, class List>
class KernelD_converter_
: public KernelD_converter_<Final_,K1,K2,typename List::tail>
{
        typedef typename List::head Tag_;
        typedef typename List::tail Rest;
        typedef KernelD_converter_<Final_,K1,K2,Rest> Base;
        typedef typename Get_type<K1,Tag_>::type K1_Obj;
        typedef typename Get_type<K2,Tag_>::type K2_Obj;
        typedef typename Get_functor<K1, Convert_ttag<Tag_> >::type K1_Conv;
        typedef KO_converter<Tag_,K1,K2> KOC;
        typedef std::is_same<K1_Conv, Null_functor> no_converter;
        typedef typename internal::Map_taglist_to_typelist<K1,Rest>::type::template contains<K1_Obj> duplicate;

        // Disable the conversion in some cases:
        struct Do_not_use{};

        // Explicit calls to boost::mpl functions to avoid parenthesis
        // warning on some versions of GCC
        typedef typename boost::mpl::if_ <
                          // If Point==Vector, keep only one conversion
          boost::mpl::or_<boost::mpl::bool_<duplicate::value>,
                          // For iterator objects, the default is make_transforming_iterator
                          boost::mpl::bool_<(iterator_tag_traits<Tag_>::is_iterator && no_converter::value)> >,
          Do_not_use,K1_Obj>::type argument_type;
        //typedef typename KOC::argument_type K1_Obj;
        //typedef typename KOC::result_type K2_Obj;
        public:
  using Base::operator(); // don't use directly, just make it accessible to the next level
        K2_Obj helper(K1_Obj const& o, std::true_type)const{
                return KOC()(this->myself().kernel(),this->myself().kernel2(),this->myself(),o);
        }
        K2_Obj helper(K1_Obj const& o, std::false_type)const{
                return K1_Conv(this->myself().kernel())(this->myself().kernel2(),this->myself(),o);
        }
        K2_Obj operator()(argument_type const& o)const{
          return helper(o,no_converter());
        }
        template<class X,int=0> struct result:Base::template result<X>{};
        template<int i> struct result<Final_(argument_type),i> {typedef K2_Obj type;};
};

template<class Final_, class K1, class K2>
class KernelD_converter_<Final_,K1,K2,typeset<> > {
        public:
        struct Do_not_use2{};
        void operator()(Do_not_use2)const{}
        template<class T> struct result;
        Final_& myself(){return *static_cast<Final_*>(this);}
        Final_ const& myself()const{return *static_cast<Final_ const*>(this);}
};


// TODO: use the intersection of Kn::Object_list.
template<class K1, class K2, class List_=
typename typeset_intersection<typename K1::Object_list, typename K2::Object_list>::type
//typeset<Point_tag>::add<Vector_tag>::type/*::add<Segment_tag>::type*/
> class KernelD_converter
        : public Store_kernel<K1>, public Store_kernel2<K2>,
        public KernelD_converter_<KernelD_converter<K1,K2,List_>,K1,K2,List_>
{
        typedef KernelD_converter Self;
        typedef Self Final_;
        typedef KernelD_converter_<Self,K1,K2,List_> Base;
        typedef typename Get_type<K1, FT_tag>::type FT1;
        typedef typename Get_type<K2, FT_tag>::type FT2;
        typedef NT_converter<FT1, FT2> NTc;
        NTc c; // TODO: compressed storage as this is likely empty and the converter gets passed around (and stored in iterators)

        public:
        KernelD_converter(){}
        KernelD_converter(K1 const&a,K2 const&b):Store_kernel<K1>(a),Store_kernel2<K2>(b){}

        // For boost::result_of, used in transforming_iterator
        template<class T,int i=is_iterator<T>::value?42:0> struct result:Base::template result<T>{};
        template<class T> struct result<Final_(T),42> {
                typedef transforming_iterator<Final_,T> type;
        };
        template<int i> struct result<Final_(K1),i>{typedef K2 type;};
        template<int i> struct result<Final_(int),i>{typedef int type;};
        // Ideally the next 2 would come with Point_tag and Vector_tag, but that's hard...
        template<int i> struct result<Final_(Origin),i>{typedef Origin type;};
        template<int i> struct result<Final_(Null_vector),i>{typedef Null_vector type;};
        template<int i> struct result<Final_(Object),i>{typedef Object type;};
        template<int i> struct result<Final_(FT1),i>{typedef FT2 type;};

        using Base::operator();
        typename Store_kernel2<K2>::reference2_type operator()(K1 const&)const{return this->kernel2();}
        int operator()(int i)const{return i;}
        Origin operator()(Origin const&o)const{return o;}
        Null_vector operator()(Null_vector const&v)const{return v;}
        FT2 operator()(FT1 const&x)const{return c(x);}
        //RT2 operator()(typename First_if_different<RT1,FT1>::Type const&x)const{return cr(x);}

        typename Get_type<K2, Flat_orientation_tag>::type const&
        operator()(typename Get_type<K1, Flat_orientation_tag>::type const&o)const
        { return o; } // Both kernels should have the same, returning a reference should warn if not.

        template<class It>
        transforming_iterator<Final_,typename boost::enable_if<is_iterator<It>,It>::type>
        operator()(It const& it) const {
                return make_transforming_iterator(it,*this);
        }

        template<class T>
        //TODO: use decltype in C++11 instead of result
        std::vector<typename result<Final_(T)>::type>
        operator()(const std::vector<T>& v) const {
                return std::vector<typename result<Final_(T)>::type>(operator()(v.begin()),operator()(v.begin()));
        }

        //TODO: convert std::list and other containers?

        Object
        operator()(const Object &obj) const
        {
                typedef typename internal::Map_taglist_to_typelist<K1,List_>::type Possibilities;
                //TODO: add Empty, vector<Point>, etc to the list.
                return Object_converter<Possibilities>()(obj,*this);
        }

        //TODO: convert boost::variant

};

} //namespace CGAL

#endif // CGAL_KERNEL_D_CARTESIAN_CONVERTER_H
