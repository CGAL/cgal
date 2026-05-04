#include <CGAL/box_intersection_d.h>
#include <CGAL/Bbox_2.h>

#include <vector>
#include <cassert>
#include <iostream>

#include <CGAL/Random.h>

using ClassicBox = CGAL::Box_intersection_d::Box_d<double, 2>;
using HandleBox  = CGAL::Box_intersection_d::Box_with_handle_d<double, 2, std::vector<int>::iterator>;
using InfoBox  = CGAL::Box_intersection_d::Box_with_info_d<double, 2, int>;
using Bb = CGAL::Bbox_2;

double min = 0;
double max = 10;

Bb random_box(CGAL::Random& r){
    double x1 = r.get_double(min, max);
    double y1 = r.get_double(min, max);

    double x2 = r.get_double(min, max);
    double y2 = r.get_double(min, max);
    return Bb(std::min(x1, x2), std::min(y1, y2),
              std::max(x1, x2), std::max(y1, y2));
}

void random_test(CGAL::Random& r){
  std::vector<ClassicBox> vec_A;
  std::vector<ClassicBox> vec_B;
  std::vector<InfoBox> vec_C;

  const int n=40;

  for(int i=0; i<n; ++i)
    vec_A.emplace_back(random_box(r));
   std::vector<ClassicBox> vec_A_copy(vec_A.begin(), vec_A.end());

  for(int i=0; i<n; ++i){
    Bb bb = random_box(r);
    vec_B.emplace_back(bb);
    vec_C.emplace_back(bb, 0);
  }

  int nb_inter_AB = 0;
  int nb_inter_AC = 0;

  CGAL::box_intersection_d(vec_A.begin(), vec_A.end(), vec_B.begin(), vec_B.end(), [&](auto /*a*/, auto /*b*/){ nb_inter_AB++; });
  CGAL::box_intersection_d(vec_A_copy.begin(), vec_A_copy.end(), vec_C.begin(), vec_C.end(), [&](auto /*a*/, auto /*c*/){ nb_inter_AC++; });
  assert(nb_inter_AB == nb_inter_AC);
}

int main()
{
    const int n = 5;

    std::vector<InfoBox> infoBoxes;
    for (int i=0; i<n; ++i)
      infoBoxes.emplace_back(Bb(float(i)+0.5, 0, float(i)+1.5, 1), i);

    std::vector<int> values(n+1);
    std::iota(values.begin(), values.end(), 0);
    std::vector<HandleBox> handleBoxes;
    for (int i=0; i<=n; ++i)
      handleBoxes.emplace_back(Bb(i, 0, i+1, 1), values.begin()+i);

    std::size_t nb_intersections=0;
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


    CGAL::Random r;
    for(int j=0; j<2; ++j)
      random_test(r);

    return 0;
}
