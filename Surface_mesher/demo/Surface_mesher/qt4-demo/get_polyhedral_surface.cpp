#include "get_polyhedral_surface.h"
#include "polyhedral_surface.h"

Surface* get_polyhedral_surface(QObject* parent,
				double sharp_edges_angle_lower_bound,
				double sharp_edges_angle_upper_bound = 180.)
{
  return new Polyhedral_surface(parent,
				sharp_edges_angle_lower_bound,
				sharp_edges_angle_upper_bound);
}
