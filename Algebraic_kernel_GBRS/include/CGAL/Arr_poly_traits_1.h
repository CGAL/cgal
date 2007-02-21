// Copyright (c) 2006 Inria Lorraine (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
// Author(s)     : Luis Peñaranda <penarand@loria.fr>

#ifndef CGAL_ARR_POLY_TRAITS_1_H
#define CGAL_ARR_POLY_TRAITS_1_H

#include <CGAL/Object.h>
#include <CGAL/Gbrs_algebraic_kernel.h>
#include <vector>
#include <iterator>

CGAL_BEGIN_NAMESPACE

template<class Kernel_>class Curve;
template<class Kernel_>class Point;

template<class Kernel_>
class Arr_poly_traits_1:public Kernel_{
  public:
    Arr_poly_traits_1(){};
    ~Arr_poly_traits_1(){};
    // XXX: remove later unused typedefs
    typedef Kernel_                       K;
    typedef typename K::Polynomial_1      Polynomial;
    typedef typename K::Algebraic_real_1  Algebraic;
    typedef Curve<K>                      X_monotone_curve_2;
    typedef X_monotone_curve_2            Curve_2;
    typedef Point<K>                      Point_2;

    typedef unsigned int                  Multiplicity;

    // these typedefs will both change
    typedef Tag_false                     Has_left_category;
    typedef Tag_false                     Has_merge_category;

    // functor types of ArrangementBasicTraits_2
    class Compare_x_2:public K{
      public:
        inline Comparison_result operator()(const Point_2 &p1,
            const Point_2 &p2)const{
          return(K::construct_compare_1_object()
              (const_cast<Algebraic&>(p1.get_x()),
               const_cast<Algebraic&>(p2.get_x())));
        };
    };

    class Compare_xy_2:public K{
      public:
        Comparison_result operator()(const Point_2 &p1,const Point_2 &p2)const{
          Comparison_result comp_x=K::construct_compare_1_object()
            (const_cast<Algebraic&>(p1.get_x()),
             const_cast<Algebraic&>(p2.get_x()));
          if(comp_x!=EQUAL)
            return comp_x;
          switch(sign_1((p1.pol()-p2.pol()),p2.get_x())){
            case ZERO:return EQUAL;break;
            case POSITIVE:return LARGER;break;
            default:return SMALLER;
          }
        };
    };

    class Construct_min_vertex_2{
      public:
        inline Point_2 operator()(const X_monotone_curve_2 &c)const{
          return c.get_left_point();
        };
    };

    class Construct_max_vertex_2{
      public:
        inline Point_2 operator()(const X_monotone_curve_2 &c)const{
          return c.get_right_point();
        };
    };

    class Is_vertical_2{
      public:
        inline bool operator()(const X_monotone_curve_2 &c)const{
          return false;
        };
    };

    class Compare_y_at_x_2{
      public:
        Comparison_result operator()(const Point_2 &p,
            const X_monotone_curve_2 &c)const{
          if(p.pol()==c.pol())
            return EQUAL;
          switch(sign_1((p.pol()-c.pol()),p.get_x())){
            case ZERO:return EQUAL;break;
            case POSITIVE:return LARGER;break;
            default:return SMALLER;
          }
        };
    };

    class Compare_y_at_x_left_2{
      public:
        Comparison_result operator()(const X_monotone_curve_2 &c1,
            const X_monotone_curve_2 &c2,const Point_2 &p)const{
          if(c1.pol()==c2.pol())
            return EQUAL;
          switch(sign_1(c1.pol().derive()-c2.pol().derive(),p.get_x())){
            case NEGATIVE:return LARGER;break;
            case POSITIVE:return SMALLER;break;
                          // TODO:
                          // -compare the signs of the second derivative, if they differ
                          // we know the convexities
                          // -if this fails, compare the curvatures, see
                          // http://mathworld.wolfram.com/RadiusofCurvature.html
            default:std::cout<<"c1="<<c1<<", c2="<<c2<<", p="<<p<<std::endl;
                    CGAL_assertion_msg(false,"not implemented");
          };
          return EQUAL;
        };
    };

    class Compare_y_at_x_right_2{
      public:
        Comparison_result operator()(const X_monotone_curve_2 &c1,
            const X_monotone_curve_2 &c2,const Point_2 &p)const{
          if(c1.pol()==c2.pol())
            return EQUAL;
          switch(sign_1(c1.pol().derive()-c2.pol().derive(),p.get_x())){
            case NEGATIVE:return SMALLER;break;
            case POSITIVE:return LARGER;break;
            // TODO: see the function above
            default:std::cout<<"c1="<<c1<<", c2="<<c2<<", p="<<p<<std::endl;
                    CGAL_assertion_msg(false,"not implemented");
          };
          return EQUAL;
        };
    };

    class Equal_2:public K{
      public:
        inline bool operator()(const Point_2 &p1,const Point_2 &p2)const{
          Comparison_result comp_x=K::construct_compare_1_object()
            (const_cast<Algebraic&>(p1.get_x()),
             const_cast<Algebraic&>(p2.get_x()));
          if(comp_x!=EQUAL){
            return false;}
          bool eq=(sign_1((p1.pol()-p2.pol()),p2.get_x())==ZERO);
          return eq;
        };
        inline bool operator()(const X_monotone_curve_2 &c1,
            const X_monotone_curve_2 &c2)const{
          return((c1.pol()==c2.pol())&&
              (K::construct_compare_1_object()
               (const_cast<Algebraic&>(c1.left()),
                const_cast<Algebraic&>(c2.left()))==EQUAL)&&
              (K::construct_compare_1_object()
               (const_cast<Algebraic&>(c1.right()),
                const_cast<Algebraic&>(c2.right()))==EQUAL));
        };
    };

    class Split_2{
      public:
        void operator()(const X_monotone_curve_2 &c,const Point_2 &p,
            X_monotone_curve_2 &c1,X_monotone_curve_2 &c2)const{
          c1=X_monotone_curve_2(c.pol(),c.left(),p.get_x());
          c2=X_monotone_curve_2(c.pol(),p.get_x(),c.right());
        };
    };

    // functor types of ArrangementXMonotoneTraits_2
    class Intersect_2:public K{
      public:
        template<class OutputIterator>
          OutputIterator operator()(const X_monotone_curve_2 &cv1,
              const X_monotone_curve_2 &cv2,OutputIterator oi)const{
            CGAL_assertion(cv1.is_consistent()&&cv2.is_consistent());
            // return the empty iterator if the curves don't overlap in x
            if(((K::construct_compare_1_object()
                    (const_cast<Algebraic&>(cv1.right()),
                     const_cast<Algebraic&>(cv2.left()))==SMALLER))||
                ((K::construct_compare_1_object()
                  (const_cast<Algebraic&>(cv2.right()),
                   const_cast<Algebraic&>(cv1.left())))==SMALLER)){
              return oi;
            }
            Algebraic left_endpoint,right_endpoint;
            // calculate the left endpoint of the intersection
            if((K::construct_compare_1_object()
                  (const_cast<Algebraic&>(cv1.left()),
                   const_cast<Algebraic&>(cv2.left())))==LARGER)
              left_endpoint=cv1.left();
            else
              left_endpoint=cv2.left();
            // calculate the right endpoint of the intersection
            if(K::construct_compare_1_object()
                (const_cast<Algebraic&>(cv1.right()),
                 const_cast<Algebraic&>(cv2.right()))==SMALLER)
              right_endpoint=cv1.right();
            else
              right_endpoint=cv2.right();
            // when polynomials are the same, it's easy
            if(cv1.pol()==cv2.pol()){
              X_monotone_curve_2 curve(cv1.pol(),left_endpoint,right_endpoint);
              Object obj=make_object(curve);
              return(*oi++=obj);
            }
            // sorry, we have to call RS
            Polynomial *temp=new Polynomial(cv1.pol()-cv2.pol());
            typedef typename std::vector<Algebraic>::iterator itv;
            std::vector<Algebraic>roots;
            std::insert_iterator<std::vector<Algebraic> >it_r(roots,roots.begin());
            itv rootsit;
            K::construct_solve_1_object()(*temp,it_r,false);
            for(rootsit=roots.begin();rootsit!=roots.end();++rootsit){
              if((K::construct_compare_1_object()
                    (left_endpoint,*rootsit)!=LARGER)&&
                  (K::construct_compare_1_object()
                   (right_endpoint,*rootsit)!=SMALLER)){
                std::pair<Point_2,Multiplicity>pm=
                  std::make_pair(Point_2(*rootsit,cv2),rootsit->mult());
                Object obj=make_object(pm);
                *oi++=obj;
              }
            }
            return oi;
          }
    };

    // functor types of ArrangementTraits_2
    class Make_x_monotone_2{
      public:
        template<class OutputIterator>
          OutputIterator operator()(const Curve_2 &c,
              OutputIterator oi)const{
            return(*oi++=make_object(c));
          };
    };

    // functions to acces ArrangementBasicTraits_2 functors
    inline Compare_x_2 compare_x_2_object()const{return Compare_x_2();};
    inline Compare_xy_2 compare_xy_2_object()const{return Compare_xy_2();};
    inline Construct_min_vertex_2 construct_min_vertex_2_object()const{
      return Construct_min_vertex_2();};
    inline Construct_max_vertex_2 construct_max_vertex_2_object()const{
      return Construct_max_vertex_2();};
    inline Is_vertical_2 is_vertical_2_object()const{return Is_vertical_2();};
    inline Compare_y_at_x_2 compare_y_at_x_2_object()const{
      return Compare_y_at_x_2();};
    inline Compare_y_at_x_left_2 compare_y_at_x_left_2_object()const{
      return Compare_y_at_x_left_2();};
    inline Compare_y_at_x_right_2 compare_y_at_x_right_2_object()const{
      return Compare_y_at_x_right_2();};
    inline Equal_2 equal_2_object()const{return Equal_2();};

    // functions to acces ArrangementXMonotoneTraits_2 functors
    inline Intersect_2 intersect_2_object()const{return Intersect_2();};
    inline Split_2 split_2_object()const{return Split_2();};

    // functions to acces ArrangementTraits_2 functors
    inline Make_x_monotone_2 make_x_monotone_2_object()const{
      return Make_x_monotone_2();};
};

template<class Kernel_>class Point{
  private:
    typedef Kernel_                       K;
    typedef Curve<K>                      Curve_2;
    typedef typename K::Algebraic_real_1  Algebraic;
    typedef typename K::Polynomial_1      Polynomial;
    Algebraic *x;
    double y;
    Polynomial *p;
  public:
    Point():x(NULL),y(0),p(NULL){};
    Point(const Algebraic& x_,const Curve_2& c_){
      set_pol(c_.pol());
      if(x_.is_consistent()){
        x=new Algebraic(x_);
        // we evaluate y in the center of the interval x
        y=pol().eval_d(mpfi_get_d(x->mpfi()));
      }else
        y=0;
    };
    Point(const Algebraic& x_,double y_):x(&x_),y(y_),p(NULL){};
    ~Point(){};
    inline const Algebraic& get_x()const{return *x;};
    inline double get_y()const{return y;};
    inline void set_pol(const Polynomial &poly){
      p=const_cast<Polynomial*>(&poly);};
    inline const Polynomial& pol()const{return *p;};
    inline std::ostream& show(std::ostream &o)const{
      return(o<<"("<<get_x()<<","<<get_y()<<")");}
};

template<class K>
inline std::ostream& operator<<(std::ostream &o,const Point<K> &p){return p.show(o);}

template<class Kernel_>class Curve{
  private:
    typedef Kernel_                       K;
    typedef typename K::Polynomial_1      Polynomial;
    typedef typename K::Coefficient       Coefficient;
    typedef typename K::Algebraic_real_1  Algebraic;
    typedef Point<K>                      Point;
    Polynomial *p;
    Algebraic *l,*r;
    bool has_l,has_r;
  public:
    Curve():p(NULL),l(NULL),r(NULL){};
    Curve(const Polynomial &poly):has_l(false),has_r(false){
      set_pol(poly);};
    Curve(const Polynomial &poly,const Algebraic &left,const Algebraic &right){
      set_pol(poly);
      l=new Algebraic(left);
      r=new Algebraic(right);
      has_l=has_r=true;
    };
    const Curve& operator=(const Curve &c){
      set_pol(c.pol());
      l=new Algebraic(c.left());
      r=new Algebraic(c.right());
      has_l=has_r=true;
      return *this;
    };
    template<typename T>
    Curve(const Polynomial &poly,T left,T right){
      set_pol(poly);
      l=new Algebraic(left);
      r=new Algebraic(right);
      has_l=has_r=true;
    };
    ~Curve(){};
    inline void set_pol(const Polynomial &poly){
      p=const_cast<Polynomial*>(&poly);};
    inline const Polynomial& pol()const{return *p;};
    inline const Algebraic& left()const{return *l;};
    inline const Algebraic& right()const{return *r;};
    inline const Point& get_left_point()const{
      return *(new Point(left(),Curve(pol())));};
    inline const Point& get_right_point()const{
      return *(new Point(right(),Curve(pol())));};
    inline bool is_consistent()const{return (p?true:false);};
    std::ostream& show(std::ostream &o)const{
      o<<"("<<pol()<<",";
      if(has_l)o<<left()<<",";
      else o<<"-inf,";
      if(has_r)o<<right()<<")";
      else o<<"+inf)";
      return o;
    };
};

template<class K>
inline std::ostream& operator<<(std::ostream &o,const Curve<K> &c){return c.show(o);}

CGAL_END_NAMESPACE

#endif  // CGAL_ARR_POLY_TRAITS_1_H

/* do you use emacs? sorry
vim:ts=2 shiftwidth=2 expandtab
*/
