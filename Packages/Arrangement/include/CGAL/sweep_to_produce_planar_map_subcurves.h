#ifndef SWEEP_TO_PRODUCE_PLANAR_MAP_SUBCURVES_H
#define SWEEP_TO_PRODUCE_PLANAR_MAP_SUBCURVES_H

#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif

#ifndef CGAL_SWEEP_CURVES_TO_SUBCURVES_H
#include <CGAL/Sweep_curves_to_subcurves.h>
#endif

CGAL_BEGIN_NAMESPACE

template <class Curve_iterator, class Traits, class Container>
void  sweep_to_produce_planar_map_subcurves(Curve_iterator curves_begin, 
                                            Curve_iterator curves_end,  
                                            Traits& traits, 
                                            Container &subcurves,
                                            bool overlapping = false)
{
  Sweep_curves_to_subcurves<Curve_iterator, Traits, Container>  sweep_line;
  
  sweep_line.sweep_curves_to_subcurves(curves_begin, curves_end, 
                                       subcurves, overlapping);
}

/*template <class Curve_iterator, class PM, class Notifier>
  void  sweep_curves_to_pm(Curve_iterator curves_begin, Curve_iterator curves_end, PM &result, Notifier* notifier)
  {
  Sweep_line<Curve_iterator, PM, Notifier>  sweep_line;
  
  sweep_line.sweep_curves_to_pm(curves_begin, curves_end, result, notifier);
  }*/

CGAL_END_NAMESPACE

#endif



