#ifndef _COLOR_RAMP_H
#define _COLOR_RAMP_H

#include <QtCore/qglobal.h>
#ifdef scene_color_ramp_EXPORTS
#  define SCENE_COLOR_RAMP_EXPORT Q_DECL_EXPORT
#else
#  define SCENE_COLOR_RAMP_EXPORT Q_DECL_IMPORT
#endif

#include <list>

class SCENE_COLOR_RAMP_EXPORT Color_component
{
  typedef std::list<std::pair<double,double> > Values;

public:
  Color_component();
  Color_component(const double c0, const double c1);
  ~Color_component() {}

  double interpolate(const double v) const;
  void add(const double v, double color);
  void rebuild(const double c0, const double c1);
  void print() const;

private:
  inline Values::const_iterator next_it(const double v) const;
  inline Values::iterator next_it(const double v);

private:
  Values values_;
};


class SCENE_COLOR_RAMP_EXPORT Color_ramp
{
public :
        Color_ramp();
        Color_ramp(const double r0, const double r1, const double g0, const double g1,
                   const double b0, const double b1);

public :
  inline double r(double v) const;
  inline double g(double v) const;
  inline double b(double v) const;

        void build_red();
        void build_blue();
        void build_thermal();
        void build_rainbow();

  void print() const;

private :
  Color_component r_;
  Color_component g_;
  Color_component b_;
};


inline
Color_component::Values::const_iterator
Color_component::
next_it(const double v) const
{
  Values::const_iterator next = values_.begin();
  while ( next != values_.end() && v >= next->first ) { ++next; }
  return next;
}

inline
Color_component::Values::iterator
Color_component::
next_it(const double v)
{
  Values::iterator next = values_.begin();
  while ( next != values_.end() && v >= next->first ) { ++next; }
  return next;
}

inline
double
Color_ramp::r(double v) const
{ return r_.interpolate(v); }

inline
double
Color_ramp::g(double v) const
{ return g_.interpolate(v);  }

inline
double
Color_ramp::b(double v) const
{ return b_.interpolate(v);  }


#endif // _COLOR_RAMP_H
