#ifndef CGAL_LINE_TO_EPOINT_H
#define CGAL_LINE_TO_EPOINT_H

CGAL_BEGIN_NAMESPACE

template <class Kernel_>
struct Line_to_epoint {
  typedef Kernel_ Kernel;
  typedef typename Kernel::RT RT;
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Line_2 Line_2;
  enum point_type { SWCORNER=1, LEFTFRAME, NWCORNER, 
                    LOWERFRAME, STANDARD, UPPERFRAME,
                    SECORNER, RIGHTFRAME, NECORNER };


  static RT dx(const Line_2& l) { return l.b(); }
  static RT dy(const Line_2& l) { return -l.a(); }

  static FT ordinate_distance(const Line_2& l)
  { return Kernel::make_FT(-l.c(),l.b()); }

  static point_type determine_type(const Line_2& l)
  {
    RT adx = CGAL_NTS abs(dx(l)), ady = CGAL_NTS abs(dy(l));
    int sdx = CGAL_NTS sign(dx(l)), sdy = CGAL_NTS sign(dy(l));
    int cmp_dx_dy = CGAL_NTS compare(adx,ady), s(1);
    if (sdx < 0 && ( cmp_dx_dy > 0 || cmp_dx_dy == 0 && 
        sdy != (s = CGAL_NTS sign(ordinate_distance(l))))) {
      if (0 == s) return ( sdy < 0 ? SWCORNER : NWCORNER );
      else        return LEFTFRAME;
    } else if (sdx > 0 && ( cmp_dx_dy > 0 || cmp_dx_dy == 0 && 
               sdy != (s = CGAL_NTS sign(ordinate_distance(l))))) { 
      if (0 == s) return ( sdy < 0 ? SECORNER : NECORNER );
      else        return RIGHTFRAME;
    } else if (sdy < 0 && ( cmp_dx_dy < 0 || cmp_dx_dy == 0 && 
               ordinate_distance(l) < FT(0))) {
      return LOWERFRAME;
    } else if (sdy > 0 && ( cmp_dx_dy < 0 || cmp_dx_dy == 0 && 
               ordinate_distance(l) > FT(0))) {
      return UPPERFRAME;
    }
    CGAL_assertion_msg(false," determine_type: degenerate line.");
    return (point_type)-1; // never come here
  }


};

CGAL_END_NAMESPACE
#endif //CGAL_LINE_TO_EPOINT_H

