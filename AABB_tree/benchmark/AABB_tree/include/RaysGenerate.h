#pragma once 

#include <vector>
#include <CGAL/Random.h>

struct directions{
    double _x, _y, _z;
};


class RaysGenerate{
private:
public:
    RaysGenerate(int n);
    ~RaysGenerate();

    int _n;
    std::vector<directions> rayDirections;
};

RaysGenerate::RaysGenerate(int n): _n(n){

    CGAL::Random rand;

    for(size_t i=0; i!=n; ++i){

        directions direction;
        direction._x = rand.get_double(-1.0, 1.0);
        direction._y = rand.get_double(-1.0, 1.0); 
        direction._z = rand.get_double(-1.0, 1.0); 
        rayDirections.push_back(direction);
    }
}

RaysGenerate::~RaysGenerate(){}
