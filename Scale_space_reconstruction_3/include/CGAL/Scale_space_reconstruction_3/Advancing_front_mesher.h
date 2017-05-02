#ifndef CGAL_SCALE_SPACE_RECONSTRUCTION_3_ADVANCING_FRONT_MESHER_H
#define CGAL_SCALE_SPACE_RECONSTRUCTION_3_ADVANCING_FRONT_MESHER_H

#include <CGAL/Advancing_front_surface_reconstruction.h>

#include <CGAL/Union_find.h>

namespace CGAL
{

namespace Scale_space_reconstruction_3
{
  
template <typename Geom_traits>
class Advancing_front_mesher
{
public:
  typedef typename Geom_traits::FT FT;
  typedef typename Geom_traits::Point_3                        Point;          ///< defines the point type.
  
  typedef CGAL::cpp11::array< unsigned int, 3 >       Facet;
private:

  class Priority
  {
    double bound;
  public:
    Priority (double bound)
      : bound(bound)
    {}

    template <typename AdvancingFront, typename Cell_handle>
    double operator() (const AdvancingFront& adv, Cell_handle& c,
                       const int& index) const
    {
      if(bound == 0){
        return adv.smallest_radius_delaunay_sphere (c, index);
      }

      double d = 0.;
      d = sqrt(squared_distance(c->vertex((index+1)%4)->point(),
                                c->vertex((index+2)%4)->point()));
      if(d>bound) return adv.infinity();
      d = sqrt(squared_distance(c->vertex((index+2)%4)->point(),
                                c->vertex((index+3)%4)->point()));
      if(d>bound) return adv.infinity();
      d = sqrt(squared_distance(c->vertex((index+1)%4)->point(),
                                c->vertex((index+3)%4)->point()));
      if(d>bound) return adv.infinity();

      return adv.smallest_radius_delaunay_sphere (c, index);
    }
  };

  Priority m_priority;
  FT m_radius_ratio_bound;
  FT m_beta;
  
public:

  Advancing_front_mesher (FT maximum_facet_length = 0.,
                          FT radius_ratio_bound = 5,
                          FT beta = 0.52)
    : m_priority (maximum_facet_length), m_radius_ratio_bound (radius_ratio_bound), m_beta (beta)
  {

  }

  template <typename InputIterator, typename OutputIterator>
  void operator() (InputIterator begin, InputIterator end, OutputIterator output)
  {
    CGAL::advancing_front_surface_reconstruction (begin, end, output,
                                                  m_priority,
                                                  m_radius_ratio_bound,
                                                  m_beta);
  }

};

  
} // namespace Scale_space_reconstruction_3

} // namespace CGAL

#endif // CGAL_SCALE_SPACE_RECONSTRUCTION_3_ADVANCING_FRONT_MESHER_H
