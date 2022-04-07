// Author(s)    : Samuel Hornus

#ifndef PURE_COMPLEX_TEST_GRID_H
#define PURE_COMPLEX_TEST_GRID_H

#include <vector>

template< typename K >
class Tilted_grid
{
public:
    typedef typename K::Point_d Point;
    typedef typename K::FT  FT;
    typedef std::vector<FT> FTVec;
    typedef std::vector<FTVec> PVec; // array of point coordinates

    void operator()(const int D, const int N, PVec & g) const
    {   // create N^D grid
        if( 0 >= D )
        {
            g.push_back(FTVec());
            return;
        }
        PVec h;
        (*this)(D - 1, N, h);
        g.clear();
        typename PVec::iterator hit = h.begin();
        while( hit != h.end() )
        {
            hit->push_back(FT(0));
            ++hit;
        }
        g.insert(g.end(), h.begin(), h.end());
        for( int i = 1; i < N; ++i )
        {
            hit = h.begin();
            while( hit != h.end() )
            {
                (*hit)[D-1] = FT(i);
                ++hit;
            }
            g.insert(g.end(), h.begin(), h.end());
        }
    }
};

#endif // PURE_COMPLEX_TEST_GRID_H
