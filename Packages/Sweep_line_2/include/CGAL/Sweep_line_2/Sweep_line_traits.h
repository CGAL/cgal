
#ifndef SWEEP_LINE_TRAITS_H
#define SWEEP_LINE_TRAITS_H


CGAL_BEGIN_NAMESPACE
template<class Traits>
struct Sweep_line_traits
{
  static Traits* static_traits;

  Sweep_line_traits(Traits *traits)
  {
    static_traits = traits;
  }

 
  static Traits* get_traits()
  {
    return static_traits;
  }
   
  ~Sweep_line_traits()
  {
  }

};


//TODO : move the below two line to Sweep_line_traits.C file and write instead
 //#include<CGAL/Sweep_line_2/Sweep_line_traits.C>
template<class Traits>
Traits* Sweep_line_traits<Traits>::static_traits = NULL;

CGAL_END_NAMESPACE

#endif
