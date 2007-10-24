#include "get_polyhedral_surface.h"
#include "polyhedral_surface.h"

Surface* get_polyhedral_surface()
{
  return new Polyhedral_surface();
}
