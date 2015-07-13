#ifndef CGAL_VOLUME_PLANE_THREAD_H
#define CGAL_VOLUME_PLANE_THREAD_H

#include <CGAL/Image_3.h>
#include "Volume_plane.h"

#include "Scene_segmented_image_item.h"

#include <QApplication>
#include <QThread>
#include <vector>

template<typename Word>
struct Clamp_to_one_zero_range {
  std::pair<Word, Word> min_max;
  Word operator()(const Word& inVal) {
    Word inValNorm = inVal - min_max.first;
    Word aUpperNorm = min_max.second - min_max.first;
    Word bValNorm = inValNorm / aUpperNorm;
    return bValNorm;
  }
};

class Volume_plane_thread : public QThread {
Q_OBJECT  
public:
  Volume_plane_thread(const CGAL::Image_3* img, const QString& name)
    : img(img), item(NULL), name(name) { }

  Volume_plane_interface* getItem() {
    return item;
  }

Q_SIGNALS:
  void finished(Volume_plane_thread*);

protected:
  const CGAL::Image_3* img;
  Volume_plane_interface* item;
  std::vector<float> buffer;
  QString name;
};

template<typename Word>
class X_plane_thread : public Volume_plane_thread {
  Clamp_to_one_zero_range<Word> clamper;
public:
  X_plane_thread(const CGAL::Image_3* img, const Clamp_to_one_zero_range<Word>& clamp, const QString& name)
    : Volume_plane_thread(img, name), clamper(clamp) { }
protected:
  void run();
};

template<typename Word>
class Y_plane_thread : public Volume_plane_thread {
  Clamp_to_one_zero_range<Word> clamper;
public:
  Y_plane_thread(const CGAL::Image_3* img, const Clamp_to_one_zero_range<Word>& clamp, const QString& name)
    : Volume_plane_thread(img, name), clamper(clamp) { }
protected:
  void run();
};

template<typename Word>
class Z_plane_thread : public Volume_plane_thread {
  Clamp_to_one_zero_range<Word> clamper;
public:
  Z_plane_thread(const CGAL::Image_3* img, const Clamp_to_one_zero_range<Word>& clamp, const QString& name)
    : Volume_plane_thread(img, name), clamper(clamp) { }
protected:
  void run();
};

template<typename Word>
void X_plane_thread<Word>::run() {
    buffer.reserve(img->size());
    for(unsigned int i = 0; i < img->xdim(); ++i) {
      for(unsigned int j = 0; j < img->ydim(); ++j) {
        for(unsigned int k = 0; k < img->zdim(); ++k) {
          float x = CGAL::IMAGEIO::static_evaluate<Word>(img->image(), i, j, k);
          x = clamper(x);
          buffer.push_back(x);
        }
      }
    }
    item = new Volume_plane<x_tag>(img->ydim(), img->zdim(), img->xdim(), 
                                   img->vx(), img->vy(), img->vz(), buffer);

    item->setName(name);
    item->moveToThread(QApplication::instance()->thread());
    Q_EMIT finished(this);
}

template<typename Word>
void Y_plane_thread<Word>::run() {
      buffer.reserve(img->size());
    for(unsigned int i = 0; i < img->ydim(); ++i) {
      for(unsigned int j = 0; j < img->xdim(); ++j) {
        for(unsigned int k = 0; k < img->zdim(); ++k) {
          float x = CGAL::IMAGEIO::static_evaluate<Word>(img->image(), j, i, k);
          x = clamper(x);
          buffer.push_back(x);
        }
      }
    }
    item = new Volume_plane<y_tag>(img->xdim(), img->zdim(), img->ydim(), 
                                   img->vx(), img->vy(), img->vz(), buffer);
    item->setName(name);
    item->moveToThread(QApplication::instance()->thread());
    Q_EMIT finished(this);
}

template<typename Word>
void Z_plane_thread<Word>::run() {
  for(unsigned int i = 0; i < img->zdim(); ++i) {
    for(unsigned int j = 0; j < img->xdim(); ++j) {
      for(unsigned int k = 0; k < img->ydim(); ++k) {
        float x = CGAL::IMAGEIO::static_evaluate<Word>(img->image(), j, k, i);
        x = clamper(x);
        buffer.push_back(x);
      }
    }
  }
  item = new Volume_plane<z_tag>(img->xdim(), img->ydim(), img->zdim(), 
                                 img->vx(), img->vy(), img->vz(), buffer);
  item->setName(name);
  item->moveToThread(QApplication::instance()->thread());
  Q_EMIT finished(this);
}

#endif /* CGAL_VOLUME_PLANE_THREAD_H */
