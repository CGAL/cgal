#include <CEP/Vanilla/flavor.h>
#include <LEDA/string.h>

// is f one of the valid flavors
bool valid_flavor(Flavor f)
{
    return (VANILLA <= f && f <= PISTACHIO);
}

//
// implementation of the complicated flavor_enhancement algorithm of
// Irene, et al.
//
// Precondition: f is valid
//
Flavor flavor_enhance(Flavor f)
{
   return f = (PISTACHIO == f) ? VANILLA : Flavor(f+1);
}

std::ostream& operator<<(std::ostream& os, Flavor f)
{
   switch (f) 
   {
     case VANILLA:              os << "VANILLA"; break;
     case CHOCOLATE_CHIP:       os << "CHOCOLATE_CHIP"; break;
     case MINT_CHOCOLATE_CHIP:  os << "MINT_CHOCOLATE_CHIP"; break;
     case CHOCOLATE:            os << "CHOCOLATE"; break;
     case STRAWBERRY:           os << "STRAWBERRY"; break;
     case NEAPOLITAN:           os << "NEAPOLITAN"; break;
     case PEACH:                os << "PEACH"; break;
     case ROCKY_ROAD:           os << "ROCKY_ROAD"; break;
     case PISTACHIO:            os << "PISTACHIO"; break;
   }
   return os;
}

std::istream& operator>>(std::istream& is, Flavor& f)
{
   leda_string name;

   is >> name;
   if (name == "VANILLA")
      f = VANILLA;
   else if (name == "CHOCOLATE_CHIP")
      f = CHOCOLATE_CHIP;
   else if (name == "MINT_CHOCOLATE_CHIP")
      f = MINT_CHOCOLATE_CHIP;
   else if (name == "CHOCOLATE")
      f = CHOCOLATE;
   else if (name == "STRAWBERRY")
      f = STRAWBERRY;
   else if (name == "NEAPOLITAN")
      f = NEAPOLITAN;
   else if (name == "PEACH")
      f = PEACH;
   else if (name == "ROCKY_ROAD")
      f = ROCKY_ROAD;
   else if (name == "PISTACHIO")
      f = PISTACHIO;
   else
      is.clear(std::ios::badbit);
   return is;
}
