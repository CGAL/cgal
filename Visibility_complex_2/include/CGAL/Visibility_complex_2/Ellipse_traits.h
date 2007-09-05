// Copyright (c) 2001-2004  ENS of Paris (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Luc Habert


#ifndef CGAL_VISIBILITY_COMPLEX_2_ELLIPSE_TRAITS_H
#define CGAL_VISIBILITY_COMPLEX_2_ELLIPSE_TRAITS_H

#include<CGAL/basic.h>
#include<CGAL/Point_2.h>
#include<CGAL/Conic_2.h>
#include<CGAL/Visibility_complex_2/Bitangent_2.h>
#include<CGAL/Visibility_complex_2/Arc_2.h>
#include<boost/shared_ptr.hpp>

#include<CGAL/Visibility_complex_2/CORE_Root_of_traits.h>

CGAL_BEGIN_NAMESPACE
namespace Visibility_complex_2_details {

template<class Kernel,class RootOfTraits> class Ellipse_traits;


template<class Kernel,class RootOfTraits>
class Bitangent_2<Ellipse_traits<Kernel,RootOfTraits> >
  : public Bitangent_base<
  typename Ellipse_traits<Kernel,RootOfTraits>::Disk> {
  typedef Ellipse_traits<Kernel,RootOfTraits> ET;
  typedef typename ET::Bitangent_2 Self;
public:
  typedef typename ET::Disk Disk;
  typedef typename ET::Arc_2 Arc_2;
  typedef typename Kernel::RT RT;
  typedef typename ET::Point_2 Point_2;
  typedef typename ET::Segment_2 Segment_2;
private:
  typedef Bitangent_base<Disk> Base;
  typedef typename RootOfTraits::Root_of_1 Root_of_1;
  typedef typename RootOfTraits::Root_of_2 Root_of_2;
  typedef typename RootOfTraits::Root_of_4 Root_of_4;
public:
  using Base::LL;
  using Base::LR;
  using Base::RL;
  using Base::RR;
  using Base::type;
  using Base::source_object;
  using Base::target_object;
  typedef typename Base::Disk_handle Disk_handle;
  typedef typename Base::Type Type;
  Bitangent_2(Type t,Disk_handle source,Disk_handle target) :
    Base(t,source,target),
    roots(new Cache(*source,*target)),
    reversed(false), type_in_cache(t) {
    roots->compute(type());
  };
  Bitangent_2(const Bitangent_2 &sibling,bool reverse,Type t) :
    Base(sibling,reverse,t),
    roots(sibling.roots),
    reversed(reverse?!sibling.reversed:sibling.reversed),
    type_in_cache(sibling.reversed?
                  Base::reverse(t):t) {
    roots->compute(type_in_cache);
  }
  Bitangent_2(Type t,const Arc_2&source,const Arc_2&target) {
    *this=Self(t,source.object(),target.object());
  }

  Bitangent_2() {}

  const Root_of_4 & slope() const {
    return roots->roots[roots->compute(type_in_cache)].x;
  }

  typename ET::Tangency_point_position source_position() const {
    if (reversed) {
      switch (type()) {
      case LL: case RR: return target_position();
      default:
        switch (target_position()) {
        case ET::Left: return ET::Right;
        case ET::Below: return ET::Above;
        case ET::Right: return ET::Left;
        case ET::Above: return ET::Below;
        default:
          CGAL_assertion(false);
          return ET::Above;
        }
      }
    } else return roots->roots[roots->compute(type_in_cache)].p;
  }

  typename ET::Tangency_point_position target_position() const {
    if (reversed)
      return roots->roots[roots->compute(type_in_cache)].p;
    else {
      switch (type()) {
      case LL: case RR: return source_position();
      default:
        switch (source_position()) {
        case ET::Left: return ET::Right;
        case ET::Below: return ET::Above;
        case ET::Right: return ET::Left;
        case ET::Above: return ET::Below;
        default:
          CGAL_assertion(false);
          return ET::Above;
        }
      }
    }
  }
  typename ET::Bitangent_orientation orientation() const {
    switch (type()) {
    case LL: case LR:
      switch (source_position()) {
      case ET::Left: return ET::Downwards;
      case ET::Below: return ET::Rightwards;
      case ET::Right: return ET::Upwards;
      case ET::Above: return ET::Leftwards;
      default:
        CGAL_assertion(false);
        return ET::Downwards;
      }
    case RL: case RR:
      switch (source_position()) {
      case ET::Left: return ET::Upwards;
      case ET::Below: return ET::Leftwards;
      case ET::Right: return ET::Downwards;
      case ET::Above: return ET::Rightwards;
      default:
        CGAL_assertion(false);
        return ET::Downwards;
      }
    default:
      CGAL_assertion(false);
      return ET::Downwards;
    }
  }


  
  Point_2 source() {
    double x,y;
    typename ET::Tangency_point tp;
    double slp=0;
    switch (source_position()) {
    case ET::Below: case ET::Above:
      slp=slope().to_double();
    default:;
    }
    if (reversed) {
      tp(*target_object(),&x,&y,slp,target_position());
    } else {
      tp(*source_object(),&x,&y,slp,source_position());
    }
    return Point_2(static_cast<int>(round(x)),
                   static_cast<int>(round(y)));
  }
  Point_2 target() {
    double x,y;
    typename ET::Tangency_point tp;
    double slp=0;
    switch (source_position()) {
    case ET::Below: case ET::Above:
      slp=slope().to_double();
    default:;
    }
    if (reversed) {
      tp(*source_object(),&x,&y,slp,source_position());
    } else {
      tp(*target_object(),&x,&y,slp,target_position());
    }
    return Point_2(static_cast<int>(round(x)),
                   static_cast<int>(round(y)));
  }

  operator Segment_2() {
    return Segment_2(source(),target());
  }


private:
  struct Cache {
  public:
    int root_nb;
    int computed_roots;
    typedef Ellipse_traits<Kernel,RootOfTraits> ET;
    typedef typename ET::Bitangent_2 Bitangent_2;
    typedef typename Bitangent_2::Type Type;
    typedef typename ET::Tangency_point_position Tangency_point_position;
    typename RootOfTraits::Polynom_1 S;
    typename RootOfTraits::Polynom_2 T;
    RT P4,P3,P2,P1,P0;
    struct plouf {
      typename RootOfTraits::Root_of_4 x;
      Tangency_point_position p;
      Type t;
    };
    plouf roots[4];
    bool 
    is_above(const Disk& D1,const Disk& D2,const Root_of_2& x) const {
      bool above;
      RT pleutch=D1.t()*D2.s()-
        D1.s()*D2.t();
      if (pleutch==0) {
        RT plitch=-D2.v()*D1.s()+D1.v()*D2.s();
        above=plitch>0;
      } else {
        typename RootOfTraits::Polynom_1 ploutch(pleutch,
                                                 D1.v()*D2.s()
                                                 -D1.s()*D2.v());
        CGAL::Sign platch=ploutch.sign_at(x);
        above=(platch==CGAL::POSITIVE);
      }
      if (D1.s()<0) above= !above;
      if (D2.s()<0) above= !above;
      return above;
    }
    struct bit_type {
      typename Bitangent_2::Type t;
      void incr(void) {
        switch (t) {
        case LL: t=LR; break;
        case LR: t=RL; break;
        case RL: t=RR; break;
        case RR: t=LL; break;
        }
      };
    };
    bit_type up,down;
    Cache() {};
    Cache(const Disk& D1,const Disk& D2) {
      typename RootOfTraits::compare_object compare;    

      RT t11=(D1.v()*D1.v()-4*D1.s()*D1.w());
      RT t12=(4*D1.u()*D1.s()-2*D1.t()*D1.v());
      RT t13=(D1.t()*D1.t()-4*D1.r()*D1.s());
      RT t14=(2*D1.u()*D1.v()-4*D1.t()*D1.w());
      RT t15=(2*D1.u()*D1.t()-4*D1.r()*D1.v());
      RT t16=D1.u()*D1.u()-4*D1.r()*D1.w();
      RT t21=(D2.v()*D2.v()-4*D2.s()*D2.w());
      RT t22=(4*D2.u()*D2.s()-2*D2.t()*D2.v());
      RT t23=(D2.t()*D2.t()-4*D2.r()*D2.s());
      RT t24=(2*D2.u()*D2.v()-4*D2.t()*D2.w());
      RT t25=(2*D2.u()*D2.t()-4*D2.r()*D2.v());
      RT t26=D2.u()*D2.u()-4*D2.r()*D2.w();

      RT Q2=t13;
      RT Q1=-t12;
      RT Q0=t11;
      RT R2=t23;
      RT R1=-t22;
      RT R0=t21;


      P4=t23*t23*t11*t11+t13*t13*t21*t21-2*t23*t11*t13*t21-t11*t23*t12*t22
        +t12*t12*t21*t23-t12*t13*t21*t22+t11*t13*t22*t22;
      P3=t14*t13*t22*t22-t11*t23*t12*t25-t12*t13*t24*t22-
        t12*t13*t21*t25+2*t12*t21*t23*t15-2*t23*t11*t13*t24+t12*t12*t24*t23+
        2*t11*t13*t22*t25+2*t13*t13*t21*t24-t11*t22*t23*t15+2*t23*t23*t11*t14-
        t12*t23*t14*t22-2*t13*t21*t23*t14-t15*t13*t21*t22;
      P2=2*t13*t13*t21*t26+t23*t23*t14*t14+2*t12*t24*t23*
        t15-2*t23*t11*t13*t26+t21*t23*t15*t15-t14*t22*t23*t15+t11*t13*t25*t25-
        t15*t13*t21*t25-2*t13*t21*t23*t16-t12*t23*t14*t25+2*t23*t23*t11*t16
        -t15*t13*t24*t22-t11*t23*t15*t25-t12*t13*t24*t25-t12*t23*t16*t22+
        2*t14*t13*t22*t25+t13*t13*t24*t24+t16*t13*t22*t22-t12*t13*t26*t22-
        2*t23*t14*t13*t24+t12*t12*t26*t23;
      P1=t24*t23*t15*t15-2*t23*t14*t13*t26+2*t13*t13*t24*t26-t15*t13*t24*t25+
        2*t12*t26*t23*t15+2*t23*t23*t14*t16-t12*t13*t26*t25+2*t16*t13*t22*t25-
        t14*t23*t15*t25+t14*t13*t25*t25-2*t13*t24*t23*t16-t12*t23*t16*t25-
        t15*t23*t16*t22-t15*t13*t26*t22;
      P0=t26*t23*t15*t15-t15*t13*t26*t25-2*t23*t16*t13*t26+t13*t13*t26*t26+
        t23*t23*t16*t16+t16*t13*t25*t25-t15*t23*t16*t25;


      S=typename RootOfTraits::Polynom_1 (t23*t12-t13*t22,t23*t15-t13*t25);
      T=typename RootOfTraits::Polynom_2 (
        (-2*t13*t13*t21-t12*t12*t23+t12*t13*t22+2*t13*t23*t11),
        (2*t13*t23*t14+t12*t13*t25+t15*t13*t22-2*t13*t13*t24-2*t15*t23*t12),
        2*t13*t23*t16-2*t13*t13*t26-t23*t15*t15+t15*t13*t25);


      Root_of_2 xmin1(Q2,Q1,Q0,0);
      Root_of_2 xmin2(R2,R1,R0,0);
      Root_of_2 xmax1(Q2,Q1,Q0,1);
      Root_of_2 xmax2(R2,R1,R0,1);

      bool TL=false; bool TR=false;

    
      switch (compare(xmin1,xmin2)) {
      case SMALLER:
        switch (compare(xmax1,xmin2)) {
        case SMALLER:
          down.t=LL;
          up.t=RL;
          break;
        case EQUAL:
          TR=true;
          if (is_above(D1,D2,xmax1)) {
            down.t=LL;
            up.t=LR;  
          } else {
            down.t=LL;
            up.t=RL;
          }
          break;
        case LARGER:
          switch (compare(xmax1,xmax2)) {
          case SMALLER:
            if (is_above(D1,D2,xmax1)) {
              down.t=LL;
              up.t=LR;          
            } else {
              down.t=LL;
              up.t=RR;
            }
            break;
          case EQUAL:
            TR=true;
            if (is_above(D1,D2,xmax1)) {
              down.t=LL;
              up.t=LL;          
            } else {
              down.t=LL;
              up.t=RL;
            }
            break;
          case LARGER:        
            up.t=LL;
            down.t=LL;
            break;
          }
        }
        break;
      case EQUAL:
        TL=true;
        switch (compare(xmax1,xmax2)) {
        case SMALLER:
          if (is_above(D1,D2,xmax1)) {
            down.t=RR;
            up.t=LR;          
          } else {
            down.t=LL;
            up.t=RR;
          }
          break;
        case EQUAL:
          TR=true;
          if (is_above(D1,D2,xmax1)) {
            down.t=RR;
            up.t=LL;          
          } else {
            down.t=LL;
            up.t=RR;
          }
          break;
        case LARGER:        
          if (is_above(D1,D2,xmin1)) {
            down.t=RR;
            up.t=LL;          
          } else {
            down.t=LL;
            up.t=LL;
          }
          break;
        }
        break;
      case LARGER:
        switch (compare(xmin1,xmax2)) {
        case SMALLER:
          switch (compare(xmax1,xmax2)) {
          case SMALLER:
            if (is_above(D1,D2,xmin1)) {
              down.t=RR;
              up.t=LR;
            } else {
              down.t=LR;
              up.t=RR;
            }
            break;
          case EQUAL:
            TR=true;
            if (is_above(D1,D2,xmin1)) {
              down.t=RR;
              up.t=LL;
            } else {
              down.t=LR;
              up.t=RR;
            }
            break;
          case LARGER:        
            if (is_above(D1,D2,xmin1)) {
              down.t=RR;
              up.t=LL;
            } else {
              down.t=LR;
              up.t=LL;
            }        
            break;
          }
          break;
        case EQUAL:
          TL=true;
          if (is_above(D1,D2,xmin1)) {
            down.t=RL;
            up.t=LL;
          } else {
            down.t=LR;
            up.t=LL;
          }      
          break;
        case LARGER:        
          down.t=RL;
          up.t=LL;
          break;
        }
        break;
      }


      root_nb=4;
      computed_roots=0;
      if (TL) {
        root_nb --;
        roots[computed_roots].p=ET::Left;
        roots[computed_roots].t=down.t;
        computed_roots++;
        down.incr();
      }
      if (TR) {
        root_nb --;
        roots[computed_roots].p=ET::Right;
        roots[computed_roots].t=up.t;
        computed_roots++;
        up.incr();
      }
    };
    int compute(Type ty) {
      for (int i=0;i<computed_roots;i++) {
        if (roots[i].t==ty) return i;
      }

      while (1) {
        CGAL_assertion(computed_roots<4);
        Sign s,t;
        roots[computed_roots].x=
          Root_of_4(P4,P3,P2,P1,P0,computed_roots-(4-root_nb));
        s=S.sign_at(roots[computed_roots].x);
        t=T.sign_at(roots[computed_roots].x);
        if (s==ZERO) {
          roots[computed_roots].t=down.t;
          roots[computed_roots].p=ET::Below;
          down.incr();
          computed_roots++;
          roots[computed_roots].x=roots[computed_roots-1].x;
          roots[computed_roots].t=up.t;
          roots[computed_roots].p=ET::Above;
          computed_roots++;
          up.incr();
          if (roots[computed_roots-2].t==ty) return computed_roots-2;
          if (roots[computed_roots-1].t==ty) return computed_roots-1;
        } else {
          if ((((s==POSITIVE)&&(t==POSITIVE))||
               ((s==NEGATIVE)&&(t==NEGATIVE)))) {
            roots[computed_roots].t=up.t;
            roots[computed_roots].p=ET::Above;
            computed_roots++;
            up.incr();
            if (roots[computed_roots-1].t==ty) return computed_roots-1;
          } else {
            roots[computed_roots].t=down.t;
            roots[computed_roots].p=ET::Below;
            down.incr();
            computed_roots++;
            if (roots[computed_roots-1].t==ty) return computed_roots-1;
          }
        }
      }
    }
  };

  boost::shared_ptr<Cache> roots;
  bool reversed;
  Type type_in_cache;
};


template<class Kernel,class RootOfTraits> class Ellipse_traits {
  typedef Ellipse_traits<Kernel,RootOfTraits> Self;
public:
  enum Tangency_point_position {
    Above, Below, Left, Right
  };
  enum Bitangent_orientation {
    Downwards, Rightwards, Upwards, Leftwards
  };
  typedef Kernel R;
  typedef typename R::RT                RT;
  typedef typename R::Point_2           Point_2;
  typedef typename R::Segment_2         Segment_2;
  typedef Conic_2<R>                    Disk;
  typedef Visibility_complex_2_details::Arc_2<Self>  Arc_2;
  typedef typename Visibility_complex_2_details::Bitangent_2<Self> Bitangent_2;
  struct Orientation_object {
    Orientation operator ()(const Bitangent_2& a,const Bitangent_2& b) const {
      typename RootOfTraits::compare_object compare;
      switch (a.orientation()) {
      case Downwards:
        switch (b.orientation()) {
        case Downwards: return COLLINEAR;
        case Rightwards: return LEFT_TURN;
        case Upwards: return COLLINEAR;
        case Leftwards: return RIGHT_TURN;
        }
      case Rightwards:
        switch (b.orientation()) {
        case Downwards: return RIGHT_TURN;
        case Rightwards: 
          switch (compare(a.slope(),b.slope())) {
          case SMALLER: return LEFT_TURN;
          case EQUAL: return COLLINEAR;
          case LARGER: return RIGHT_TURN;
        };
        case Upwards: return LEFT_TURN;
        case Leftwards: 
          switch (compare(a.slope(),b.slope())) {
          case SMALLER: return RIGHT_TURN;
          case EQUAL: return COLLINEAR;
          case LARGER: return LEFT_TURN;
          };
        }
      case Upwards:
        switch (b.orientation()) {
        case Downwards: return COLLINEAR;
        case Rightwards: return RIGHT_TURN;
        case Upwards: return COLLINEAR;
        case Leftwards: return LEFT_TURN;
        }
      case Leftwards:
        switch (b.orientation()) {
        case Downwards: return LEFT_TURN;
        case Rightwards: 
          switch (compare(a.slope(),b.slope())) {
          case SMALLER: return RIGHT_TURN;
          case EQUAL: return COLLINEAR;
          case LARGER: return LEFT_TURN;
        };
        case Upwards: return RIGHT_TURN;
        case Leftwards: 
          switch (compare(a.slope(),b.slope())) {
          case SMALLER: return LEFT_TURN;
          case EQUAL: return COLLINEAR;
          case LARGER: return RIGHT_TURN;
          };
        }
      default:
        CGAL_assertion(false);
        return COLLINEAR;
      }
    }
  };
  class Compare_extreme_yx {
  private:
    typedef Ellipse_traits<Kernel,RootOfTraits> ET;
  public:
    typedef typename ET::Disk Disk;
    typedef typename RootOfTraits::Root_of_2 Root_of_2;
    typedef typename Kernel::RT RT;
    Comparison_result operator() (bool sa , const Disk& a,
                                  bool sb , const Disk& b) const {
      if (a==b) return EQUAL;

      typename RootOfTraits::compare_object compare;
      

      RT Q2=(a.t()*a.t()-4*a.r()*a.s());
      RT Q1=(2*a.u()*a.t()-4*a.r()*a.v());
      RT Q0=a.u()*a.u()-4*a.r()*a.w();
      RT R2=(b.t()*b.t()-4*b.r()*b.s());
      RT R1=(2*b.u()*b.t()-4*b.r()*b.v());
      RT R0=b.u()*b.u()-4*b.r()*b.w();

      Root_of_2 ya(Q2,Q1,Q0,sa?0:1);
      Root_of_2 yb(R2,R1,R0,sb?0:1);
      switch (compare(ya,yb)) {
      case SMALLER: return SMALLER;
      case EQUAL: {
        bool above;
        RT pleutch=a.t()*b.r()-a.r()*b.t();
        if (pleutch==0) {
          RT plitch=-b.u()*a.r()+a.u()*b.r();
          above=plitch>0;
        } else {
          typename RootOfTraits::Polynom_1 ploutch(pleutch,
                                                   a.u()*b.r()
                                                   -a.r()*b.u());
          Sign platch=ploutch.sign_at(ya);
          above=(platch==CGAL::POSITIVE);
        }
        if (a.r()<0) above= !above;
        if (b.r()<0) above= !above;
        if (above) return SMALLER; else return LARGER; 
      };
      case LARGER: return LARGER;
      default:
        CGAL_assertion(false);
        return EQUAL;
      }
    }
  };

  struct Is_upward_directed {
    bool operator()(const Bitangent_2& b) const {
      typename RootOfTraits::compare_object compare;
      static typename RootOfTraits::Root_of_1 zero(1,0);
      switch (b.orientation()) {
      case Downwards: return false;
      case Rightwards:
        switch (compare(b.slope(),zero)) {
        case SMALLER: case EQUAL: return false;
        case LARGER: return true;
        default:
          CGAL_assertion(false);
          return false;
        }
      case Upwards: return true;
      case Leftwards:
        switch (compare(b.slope(),zero)) {
        case SMALLER: case EQUAL: return true;
        case LARGER: return false;
        default:
          CGAL_assertion(false);
          return false;
        }
      default:
        CGAL_assertion(false);
        return false;
      }
    }
  };

  typedef Tag_false Supports_chi3;

  struct Orientation_infinite {
    Orientation operator() (const Bitangent_2&, 
                            const Bitangent_2&) const { 
      CGAL_assertion(false);
      return COLLINEAR; }
  };
  struct Equal_as_segments {
    bool operator() (const Bitangent_2& a, const Bitangent_2& b) const {
      return (a==b);
    }
  };
  struct Is_point {
    bool operator() (const Disk& c) const 
    { return false; }
  };

  struct Tangency_point {
    template<class flt> void operator() (const Disk& D,flt *x,flt *y,flt U,
                                   Tangency_point_position p) const {
      flt a=to_double(D.r());
      flt b=to_double(D.t());
      flt c=to_double(D.s());
      flt d=to_double(D.u());
      flt e=to_double(D.v());
      flt f=to_double(D.w());
      switch (p) {
      case Left:
        {
          flt A=b*b-4*a*c;
          flt B=2*e*b-4*c*d;
          flt C=e*e-4*c*f;
          *x=(-B+CGAL_NTS sqrt(B*B-4*A*C))/2/A;
          *y=-(b**x+e)/c/2;
          break;
        }
      case Below:
        {
          flt t11=(e*e-4*c*f);
          flt t12=(4*d*c-2*b*e);
          flt t13=(b*b-4*a*c);
          flt t14=(2*d*e-4*b*f);
          flt t15=(2*d*b-4*a*e);
          flt t16=d*d-4*a*f;
          flt A=t13;
          flt B=t15+t12*U;
          flt C=t11*U*U+t14*U+t16;
          flt V=(-B+CGAL_NTS sqrt(B*B-4*A*C))/2/A;
          *x=-(2*c*V*U+e*U+b*V+d)/(a+b*U+c*U*U)/2;
          *y=U**x+V;
                
          break;
        }
      case Right:
        {
          flt A=b*b-4*a*c;
          flt B=2*e*b-4*c*d;
          flt C=e*e-4*c*f;
          *x=(-B-CGAL_NTS sqrt(B*B-4*A*C))/2/A;      
          *y=-(b**x+e)/c/2;
          break;
        }
      case Above:
        {
          flt t11=(e*e-4*c*f);
          flt t12=(4*d*c-2*b*e);
          flt t13=(b*b-4*a*c);
          flt t14=(2*d*e-4*b*f);
          flt t15=(2*d*b-4*a*e);
          flt t16=d*d-4*a*f;
          flt A=t13;
          flt B=t15+t12*U;
          flt C=t11*U*U+t14*U+t16;
          flt V=(-B-CGAL_NTS sqrt(B*B-4*A*C))/2/A;

          *x=-(2*c*V*U+e*U+b*V+d)/(a+b*U+c*U*U)/2;
          *y=U**x+V;
          break;
        }
      }
    }    
  };
};

// #endif
}

template<class Kernel,
         class RootOfTraits=Visibility_complex_2_details::
           CORE_Root_of_traits<typename Kernel::FT> >
class Visibility_complex_2_ellipse_traits
  :public Visibility_complex_2_details::Ellipse_traits<Kernel,RootOfTraits> {};
CGAL_END_NAMESPACE

#endif // CGAL_VISIBILITY_COMPLEX_2_ELLIPSE_TRAITS_H
