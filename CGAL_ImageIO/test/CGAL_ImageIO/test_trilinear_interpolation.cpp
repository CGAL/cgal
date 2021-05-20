#include <CGAL/config.h>
#include <CGAL/Image_3.h>
#include <CGAL/Random.h>
#include <CGAL/Timer.h>

#include <cassert>
#include <iostream>
#include <fstream>


#if defined(BOOST_MSVC)
#  pragma warning(disable:4244) // double float conversion loss of data
#endif

typedef unsigned char Word;

int main() {
  CGAL::Set_ieee_double_precision pfr;

  CGAL::Image_3 image(_createImage(2, 2, 2, 1,
                                   1., 1., 1.,
                                   1, WK_FIXED, SGN_UNSIGNED));

  Word data[8] = { 0, 0, 0, 0, 0, 0, 0, 0};

  ImageIO_free(image.data());
  image.set_data(&data[0]);

  std::cerr << std::setprecision(2) << std::fixed;

  const Word c_value = 100;

  for(int x = 0; x <= 1; ++x)
    for(int y = 0; y <= 1; ++y)
      for(int z = 0; z <= 1; ++z)
      {
        data[z * 4 + y * 2 + x] = c_value;

        std::cerr << "#### data"
                  << "[" << x << "]"
                  << "[" << y << "]"
                  << "[" << z << "]"
                  << " = " << c_value << "\n";
        for(float d_x = 0.f; d_x <= 1.f; d_x += 0.499999f)
          for(float d_y = 0.f; d_y <= 1.f; d_y += 0.499999f)
            for(float d_z = 0.f; d_z <= 1.f; d_z += 0.499999f)
            {
              assert((int)d_x == 0);
              assert((int)d_y == 0);
              assert((int)d_z == 0);
              const double value =
                image.trilinear_interpolation<Word, double, double>(d_x,
                                                                    d_y,
                                                                    d_z,
                                                                    255);

              const Word label =
                image.labellized_trilinear_interpolation<Word,
                  double>(d_x, d_y, d_z, '\xFF');

              std::cerr << "val(" << d_x << ", " << d_y << " , " << d_z << ") = "
                        << value << "   -- "
                        << "label = " << (int)label << std::endl;

              assert((label == c_value) == (value >= 0.5 * c_value));

              const double sq_dist =
                (d_x - x) * (d_x - x) +
                (d_y - y) * (d_y - y) +
                (d_z - z) * (d_z - z);

              if(sq_dist <= 0.001) {
                assert(value >= 0.9 * c_value);
                assert(label == c_value);
              }
              else if(sq_dist <= 0.51*0.51/2.) {
                assert(value >= 0.79 * c_value);
                assert(label == c_value);
              }
              else if(sq_dist <= 0.51*0.51) {
                assert(value >= 0.49 * c_value);
              }
              else if(sq_dist <= 2*0.51*0.51) {
                assert(value >= 0.24 * c_value);
                assert(label == 0);
              }
              else if(sq_dist <= 3*0.51*0.51) {
                assert(value >= 0.12 * c_value);
                assert(label == 0);
              }
              else {
                assert(value <= 0.001 * c_value);
                assert(label == 0);
              }

              const float value2 =
                image.trilinear_interpolation<Word, float, float>(d_x,
                                                                  d_y,
                                                                  d_z,
                                                                  Word(0));

              const float value3 = triLinInterp(image.image(),
                                                d_x,
                                                d_y,
                                                d_z,
                                                0.f);
              std::cerr << "tri(" << d_x << ", " << d_y << " , " << d_z << ") = "
                        << value3 << std::endl;
              if(value2 != value3)
                std::cerr << std::setprecision(30)
                          << "   " << value2 << "\n!= " << value3 << std::endl;
              assert(value2 == value3);
            }

        data[z * 4 + y * 2 + x] = 0;
      }

  // BENCH
  std::cerr.precision(10);

  float max_diff = 0.f;
  CGAL::Timer timer_new_implementation;
  CGAL::Timer timer_old_implementation;

  for(int image_nb = 0; image_nb < 1000; ++image_nb)
  {
    // fill the 2x2x2 image
    for(int i = 0; i < 8; ++i) {
      data[i] = CGAL::get_default_random().uniform_smallint('\x00','\xFF');
    }

    // test the difference between the two implementations
    for(float d_x = 0.f; d_x < 0.9f; d_x += 0.1f)
      for(float d_y = 0.f; d_y < 0.9f; d_y += 0.1f)
        for(float d_z = 0.f; d_z < 0.9f; d_z += 0.1f)
        {

          const float value1 =
            image.trilinear_interpolation<Word, float, float>(d_x,
                                                              d_y,
                                                              d_z,
                                                              0);
          const float value2 = triLinInterp(image.image(),
                                            d_x,
                                            d_y,
                                            d_z,
                                            0.f);
          float diff = value1 - value2;
          if(diff < 0) {
            diff =-diff;
          }
          assert(diff < 0.1f);
          if(diff > max_diff)
            max_diff = diff;
        }
    // bench new implementation
    float sum = 0.f;
    timer_new_implementation.start();
    for(float d_x = 0.f; d_x < 0.9f; d_x += 0.05f)
      for(float d_y = 0.f; d_y < 0.9f; d_y += 0.05f)
        for(float d_z = 0.f; d_z < 0.9f; d_z += 0.05f)
        {

          sum +=
            image.trilinear_interpolation<Word, float, float>(d_x,
                                                              d_y,
                                                              d_z,
                                                              0);
        }
    timer_new_implementation.stop();

    sum = 0.f;
    timer_old_implementation.start();
    // bench old implementation
    for(float d_x = 0.f; d_x < 0.9f; d_x += 0.05f)
      for(float d_y = 0.f; d_y < 0.9f; d_y += 0.05f)
        for(float d_z = 0.f; d_z < 0.9f; d_z += 0.05f)
        {
          sum += triLinInterp(image.image(),
                              d_x,
                              d_y,
                              d_z,
                              0.f);
        }
    timer_old_implementation.stop();
  }
  std::cerr << "max difference = " << max_diff << "\n"
            << "timer new implementation: " << timer_new_implementation.time()
            << "\ntimer old implementation: " << timer_old_implementation.time()
            << "\n";
  image.set_data(nullptr); // trick to avoid ~Image_3 segfault.


  const char* filenames[] = {
    "data/skull_2.9.inr",
    "../../examples/Surface_mesher/data/skull_2.9.inr",
    "../../../Surface_mesher/examples/Surface_mesher/data/skull_2.9.inr",
    "../Surface_mesher_Examples/data/skull_2.9.inr"
  };

  std::size_t file_index = 0;
  for(   ; file_index < sizeof(filenames); ++file_index)
  {
    std::ifstream image_file(filenames[file_index], std::ios_base::binary | std::ios_base::in );
    if(image_file) {
      break;
    }
  }

  assert(file_index < sizeof(filenames) );

  std::cerr << "Opening file " << filenames[file_index] << "...\n";

  CGAL::Image_3 image2;
  const bool result = image2.read(filenames[file_index]);
  assert(result);

  std::cerr << "Image info:"
            << "\n  dim = "
            << image2.xdim() << "x" << image2.ydim() << "x" << image2.zdim()
            << "\n  v = ( " << image2.vx() << ", "
            << image2.vy() << ", " << image2.vz()
            << " )\n";

  int counter = 0;
  for(float d_x = 0.; d_x < image2.xdim() * image.vx(); d_x += 1.f/3)
    for(float d_y = 0.; d_y < image2.ydim() * image.vy(); d_y += 3.f/5)
      for(float d_z = 0.; d_z < image2.zdim() * image.vz(); d_z += 1.f/3)
      {
        ++counter;
        const float value1 =
          ::triLinInterp(image2.image(),
                         d_x,
                         d_y,
                         d_z,
                         0);
        const float value2 =
          image2.trilinear_interpolation<float, float, float>(d_x,
                                                              d_y,
                                                              d_z,
                                                              0);

        float diff = value2 - value1;
        if(diff < 0) diff = -diff;
        if(diff > 0.0001 )
        {
          std::cerr.precision(20);
          const std::size_t i1 = static_cast<std::size_t>(d_z / image2.image()->vz);
          const std::size_t j1 = static_cast<std::size_t>(d_y / image2.image()->vy);
          const std::size_t k1 = static_cast<std::size_t>(d_x / image2.image()->vx);

          std::cerr << "pos = (" << d_x << ", " << d_y << ", " << d_z << ") => ("
                    << i1 << ", " << j1 << ", " << k1 << ")\n";
          std::cerr << "value1 = " << value1
                    << "\nvalue2 = " << value2 << std::endl;

          for(std::size_t di = 0; di < 2 ; ++di)
            for(std::size_t dj = 0; dj < 2 ; ++dj)
              for(std::size_t dk = 0; dk < 2 ; ++dk)
              {
                std::cerr << "value(" << i1 + di
                          << ", " << j1 + dj
                          << ", " << k1 + dk << ") = "
                          << image2.value(i1 + di, j1+ dj, k1 + dk) << "\n";
              }

          assert(false);
        }
      }
  std::cerr << counter << " tests. OK.";
}


