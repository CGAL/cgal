#include <CGAL/box_intersection_d.h>
#include <CGAL/Bbox_2.h>

#include <vector>
#include <cassert>
#include <iostream>

using ClassicBox = CGAL::Box_intersection_d::Box_d<double, 2>;
using HandleBox  = CGAL::Box_intersection_d::Box_with_handle_d<double, 2, std::vector<int>::iterator>;
using InfoBox  = CGAL::Box_intersection_d::Box_with_info_d<double, 2, int>;

int main()
{
    const int n = 5;

    std::vector<ClassicBox> classic;
    for (int i=0; i<n; ++i) {
        classic.emplace_back(CGAL::Bbox_2(i*10.0, i*10.0, i*10.0+5.0, i*10.0+5.0));
    }

    std::vector<InfoBox> infoBoxes;
    for (int i=0; i<n; ++i)
      infoBoxes.emplace_back(CGAL::Bbox_2(float(i)+0.5, 0, float(i)+1.5, 1), i);

    std::vector<int> values(n+1);
    std::iota(values.begin(), values.end(), 0);
    std::vector<HandleBox> handleBoxes;
    for (int i=0; i<=n; ++i)
      handleBoxes.emplace_back(CGAL::Bbox_2(i, 0, i+1, 1), values.begin()+i);

    size_t nb_intersections=0;
    auto callback = [&](const InfoBox &infoBox,
                        const HandleBox &handleBox){
        int i = infoBox.info();
        int idx = *(handleBox.handle());
        assert(idx == i || idx == i+1);
        ++nb_intersections;
    };
    CGAL::box_intersection_d(infoBoxes.begin(), infoBoxes.end(), handleBoxes.begin(), handleBoxes.end(), callback);
    assert(nb_intersections==2*n);

    nb_intersections=0;
    auto callback2 = [&](const HandleBox &handleBox,
                         const InfoBox &infoBox){
        int i = infoBox.info();
        int idx = *(handleBox.handle());
        assert(idx == i || idx == i+1);
        ++nb_intersections;
    };
    CGAL::box_intersection_d(handleBoxes.begin(), handleBoxes.end(), infoBoxes.begin(), infoBoxes.end(), callback2);
    assert(nb_intersections==2*n);

    return 0;
}
