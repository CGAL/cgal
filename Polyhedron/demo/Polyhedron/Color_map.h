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

  int s = base_color.saturation(),
      v = base_color.value();
  s = s<90 ? s+90 : s;
  v = v<90 ? v+90 : v;
  float sf=s/256.0f,
      vf=v/256.0f;
  
  qreal h = hue==-1 ? 0
                    :hue;
  for(unsigned i = 0; i < nb_of_colors; ++i) {
    if (h!=-1) h += step;
    if ( h > 1 ) { h -= 1; }
//avoid S and V near 0 to avoid having all the same colors
    *out++ = QColor::fromHsvF(h, 
                              sf,
                              vf);
  }
  return out;
}

#endif
