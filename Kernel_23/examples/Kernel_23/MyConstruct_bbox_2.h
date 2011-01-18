#ifndef MYCONSTRUCT_BBOX_2_H
#define MYCONSTRUCT_BBOX_2_H


template <class ConstructBbox_2>
class MyConstruct_bbox_2 : public ConstructBbox_2 {
public:
  using ConstructBbox_2::operator();

  CGAL::Bbox_2 operator()(const MyPointC2& p) const {
    return CGAL::Bbox_2(p.x(), p.y(), p.x(), p.y());
  }
};

#endif //MYCONSTRUCT_BBOX_2_H
