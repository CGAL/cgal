#ifndef _COLOR_MAP_H
#define _COLOR_MAP_H

#include <QColor>

inline QColor generate_color(double h, double s_min = 0.35)
{
  std::size_t s_max=255;
  if(h >0.8 && h < 0.95) //span of ugly pink, desaturates make it less garish IMO
    s_max = 160;
  std::size_t s = std::rand() % (s_max-static_cast<std::size_t>(s_min*255)) + static_cast<int>(s_min*255);
  return QColor::fromHsvF(h,s/255.0,1.0);
}


template <typename Output_color_iterator>
Output_color_iterator
compute_color_map(QColor base_color,
                  std::size_t nb_of_colors,
                  Output_color_iterator out)
{
  qreal hue = base_color.hueF();
  const qreal step = ((qreal)1) / nb_of_colors;

  qreal h = hue==-1 ? 0
                    :hue;
  for(unsigned i = 0; i < nb_of_colors; ++i) {
    if (h!=-1) h += step;
    if ( h > 1 ) { h -= 1; }
    *out++ = generate_color(h);
  }
  return out;
}

inline QColor generate_random_color() {
  std::size_t h = static_cast<std::size_t>(std::rand() % 360);
  return generate_color(h/359.0);
}

#endif
