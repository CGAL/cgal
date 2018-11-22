#include <fstream>

#include <CGAL/Hyperbolic_random_points_in_disc_2.h>

// CGAL headers
#include <CGAL/Exact_circular_kernel_2.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_2.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_traits_2.h>
#include <CGAL/Timer.h>

typedef CGAL::Exact_circular_kernel_2 K;
typedef CGAL::Hyperbolic_Delaunay_triangulation_traits_2<K> Gt;

typedef K::Point_2 Point_2;
typedef K::FT FT;

typedef CGAL::Hyperbolic_Delaunay_triangulation_2<Gt> HDt;
typedef HDt::Vertex_handle Vertex_handle;
typedef HDt::Face_handle Face_handle;

int main(int argc, char *argv[])
{   
  FT eps = 0;
  for(int k = 0; k < 2; k++) {
    if(k == 0) {
      eps = 1e-3;
    }
    if(k == 1) {
      eps = 1e-7;
    }
    
    int trials_nb = 10;
    int start_nb = 10000;
    
    std::cout << std::endl << eps << std::endl << std::endl;
    
    for(int nb = start_nb, k = 0; k < 4; nb = nb*10, k++) {
      
      std::vector< std::vector<Point_2> > pts(trials_nb);
      for(int i = 0; i < trials_nb; i++) {
        if(argc > 1 && argv[1][0] == 'e') {
          Random_points_in_disc_2<K>(pts[i], nb, i, eps);
        } else {
          Hyperbolic_random_points_in_disc_2<K>(pts[i], nb, i, eps);
        }
      }
      
      double average_time = 0;
      double average_nb = 0;
      for(int trials = 0; trials < trials_nb; trials++) {
        
        HDt hdt;
        
        CGAL::Timer timer;
        
        // Hyperbolic
        
        timer.reset();
        timer.start();
        
        spatial_sort (pts[trials].begin(), pts[trials].end(), K());
        HDt::Face_handle f;
        
        for(int i = 0; i < nb; i++) {
          f = hdt.insert(pts[trials][i], f)->face();
        }
        
        timer.stop();
        
        //std::cout << "Time for H^2 " << timer.time() << std::endl;
        //std::cout << "points nb for H^2 " << hdt.number_of_vertices() << std::endl;
        //
        
        // Euclidean
        /*
        timer.reset();
        timer.start();
        
        Dt::Face_handle f2;
        for(int i = 0; i < nb; i++) {
          f2 = dt.insert(pts[trials][i], f2)->face();
        }
        
        timer.stop();
        std::cout << std::endl << "Time for E^2 " << timer.time() << std::endl;
        std::cout << "points nb for E^2 " << dt.number_of_vertices() << std::endl;
        */
        //
        
        average_time += timer.time();
        timer.reset();
        
        average_nb += hdt.number_of_vertices();
      }
      
      average_time = average_time/trials_nb;
      average_nb = average_nb/trials_nb;
    
      std::cout << "Eps: " << eps << std::endl;
      std::cout << "Number of points: " << average_nb << std::endl;
      std::cout << "Time: " << average_time << std::endl;
      std::cout << std::endl;
    }
  }
  
  return 0;
}
