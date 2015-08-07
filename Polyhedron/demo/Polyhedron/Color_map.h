#ifndef _COLOR_MAP_H
#define _COLOR_MAP_H

#include <QColor>

template <typename Output_color_iterator>
Output_color_iterator
compute_color_map(QColor base_color,
                  unsigned nb_of_colors,
                  Output_color_iterator out)
{
  qreal hue = base_color.hueF();
  const qreal step = ((qreal)1) / nb_of_colors;

  qreal h = hue;
  for(unsigned i = 0; i < nb_of_colors; ++i) {
    if (h!=-1) h += step;
    if ( h > 1 ) { h -= 1; }

    *out++ = QColor::fromHsvF(h, 
                              base_color.saturationF(), 
                              base_color.valueF());
  }
  return out;
}

#endif
