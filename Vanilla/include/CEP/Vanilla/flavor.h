#ifndef FLAVOR_H
#define FLAVOR_H

#include <CGAL/config.h>
#include <iostream>


enum Flavor {VANILLA, CHOCOLATE_CHIP, MINT_CHOCOLATE_CHIP, CHOCOLATE, 
             STRAWBERRY, NEAPOLITAN, PEACH,  ROCKY_ROAD, PISTACHIO}; 

bool valid_flavor(Flavor f);

Flavor flavor_enhance(Flavor f);

std::ostream& operator<<(std::ostream& os, Flavor f);

std::istream& operator>>(std::istream& is, Flavor& f);

#endif // FLAVOR_H
