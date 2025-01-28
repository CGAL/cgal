#include "Color_ramp.h"
#include <iostream>

// -----------------------------------
// Color_component
// -----------------------------------
Color_component::
Color_component()
{
  add(0,0);
  add(1,1);
}

Color_component::
Color_component(const double c0, const double c1)
{
  add(0,c0);
  add(1,c1);
}

double
Color_component::
interpolate(const double v) const
{
  Values::const_iterator next = next_it(v);

  // next is just after v
  Values::const_iterator prev = --next;
  ++next;

  if ( v>=1 && next != values_.end())
    std::cerr << ".";

  if ( next == values_.end() )
  {
    return prev->second;
  }

  const double& a = prev->first;
  const double& b = next->first;
  return (b-v)/(b-a) * prev->second + (v-a)/(b-a) * next->second;

}

void
Color_component::
add(const double v, double color)
{
  if ( color > 1 ) { color = 1; }
  if ( color < 0 ) { color = 0; }

  Values::iterator next = next_it(v);
  values_.insert(next, std::make_pair(v,color));
}


void
Color_component::
rebuild(const double c0, const double c1)
{
  values_.clear();
  add(1,c1);
  add(0,c0);
}

void
Color_component::
print() const
{
  for ( Values::const_iterator it = values_.begin(),
       end = values_.end() ; it != end ; ++it )
  {
    std::cout << "<" << it->first << "," << it->second << "> ";
  }

  std::cout << std::endl;
}


// -----------------------------------
// Color_ramp
// -----------------------------------
Color_ramp::Color_ramp()
  : r_()
  , g_()
  , b_()
{
}

Color_ramp::Color_ramp(const double r0, const double r1,
           const double g0, const double g1,
           const double b0, const double b1)
{
  r_.rebuild(r0,r1);
  g_.rebuild(g0,g1);
  b_.rebuild(b0,b1);

  r_.add(0.5,(std::max)(r0,r1));
  g_.add(0.5,(std::max)(g0,g1));
  b_.add(0.5,(std::max)(b0,b1));
  }

void
Color_ramp::build_red()
{
  r_.rebuild(1,0.5);
  g_.rebuild(1,0);
  b_.rebuild(1,0);

  r_.add(0.3,1);
  r_.add(0.05,1);
  g_.add(0.05,0.95);
  g_.add(0.3,0.5);
  b_.add(0.05,0);
}

void
Color_ramp::build_blue()
{
  r_.rebuild(1,0);
  g_.rebuild(1,0);
  b_.rebuild(1,0.5);

  b_.add(0.1,0.8);
  g_.add(0.1,0.4);
  r_.add(0.1,0.4);
}

void
Color_ramp::build_thermal()
{
  r_.rebuild(1,0.5);
  g_.rebuild(1,0);
  b_.rebuild(1,0);

  r_.add(0.0,0);
  r_.add(1,1);
  g_.add(0.05,0.8);
  g_.add(0.3,0.5);
  b_.add(0.0,1.0);
  b_.add(1.0,0.0);
}

void
Color_ramp::build_rainbow()
{
  r_.rebuild(1,0.5);
  g_.rebuild(1,0);
  b_.rebuild(1,0);

  r_.add(0,0.75);
  g_.add(0,0.75);
  b_.add(0,1);

  r_.add(0.2,0);
  g_.add(0.2,0);
  b_.add(0.2,1);

  r_.add(0.4,0);
  g_.add(0.4,1);
  b_.add(0.4,0);

  r_.add(0.6,1);
  g_.add(0.6,1);
  b_.add(0.6,0);

  r_.add(0.8,1);
  g_.add(0.8,0);
  b_.add(0.8,0);
}

void
Color_ramp::
print() const
{
  std::cout << "r: ";
  r_.print();
  std::cout << "g: ";
  g_.print();
  std::cout << "b: ";
  b_.print();
}

