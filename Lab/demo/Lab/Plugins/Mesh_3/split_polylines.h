#include "Scene_polylines_item.h"

typedef Scene_polylines_item::Polylines_container Polylines_container;

auto split_polylines(const Polylines_container& input)
  -> std::vector<Polylines_container::value_type>;
