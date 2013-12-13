//
#include <CGAL/internal/Surface_mesh_segmentation/K_means_clustering.h>

#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
/**
 * Generates sample points using a few gauissians.
 * Then applies k-means on these generated points.
 * Provides a heuristic score for each k-means clustering result.
 *
 * EXIT_FAILURE does not mean failure but if approximate matching is too low it is best to check
 */
int main(void)
{
    boost::mt19937 engine;
    engine.seed(1340818006);

    // generate random data using gauissians below
    std::vector< boost::normal_distribution<double> > distributions;
    distributions.push_back(boost::normal_distribution<double>(0.1, 0.05));
    distributions.push_back(boost::normal_distribution<double>(0.4, 0.1));
    distributions.push_back(boost::normal_distribution<double>(0.55, 0.05));
    distributions.push_back(boost::normal_distribution<double>(0.7, 0.1));
    distributions.push_back(boost::normal_distribution<double>(0.9, 0.05));
    distributions.push_back(boost::normal_distribution<double>(1.0, 0.05));
    
    std::vector<double> data;    
    for(std::vector< boost::normal_distribution<double> >::iterator it = distributions.begin();
      it != distributions.end(); ++it)
    {
        boost::variate_generator<boost::mt19937&, boost::normal_distribution<double> > var_nor(engine, *it);
        
        for(std::size_t i = 0; i < 300; ++i) { data.push_back(var_nor()); }
    }
    
    // calculate closest center (using above gauissians) for each generated points
    // we will compare it with k-means results
    std::vector<std::size_t> data_centers;
    for(std::vector<double>::iterator it = data.begin(); it != data.end(); ++it)
    {
        std::size_t center_id = -1, center_counter = 0;;
        double min_distance = (std::numeric_limits<double>::max)();        
        for(std::vector< boost::normal_distribution<double> >::iterator dis_it = distributions.begin();
          dis_it != distributions.end(); ++dis_it, ++center_counter)
        {
            double distance = std::abs(*it - dis_it->mean());            
            if(min_distance > distance) 
            {
                min_distance = distance;
                center_id = center_counter;
            }            
        }
        data_centers.push_back(center_id);
    }
    
    // apply k-means clustering
    typedef CGAL::internal::K_means_clustering K_means;
    std::vector<K_means> k_means;
    k_means.push_back(K_means(distributions.size(), data, K_means::PLUS_INITIALIZATION));
    k_means.push_back(K_means(distributions.size(), data, K_means::RANDOM_INITIALIZATION));

    std::vector< std::vector<std::size_t> > calculated_centers(k_means.size());
    std::vector< std::vector<std::size_t> >::iterator calc_centers_it = calculated_centers.begin();
    for(std::vector<K_means>::iterator it = k_means.begin(); it != k_means.end(); ++it, ++calc_centers_it)
    {
        it->fill_with_center_ids(*calc_centers_it);
    }
    
    std::cout << "Compare results of k-means with 'expected' (but be aware, it is not optimal result in terms of within-cluster error)" << std::endl;
    std::cout << "Another words a clustering which has smaller within-cluster error can result in worse score in here" << std::endl;
    for(std::vector< std::vector<std::size_t> >::iterator calc_centers_it = calculated_centers.begin();
        calc_centers_it != calculated_centers.end(); ++calc_centers_it)
    {
        std::size_t true_count = 0;
        std::vector<std::size_t>::iterator calculated_it = calc_centers_it->begin();
        for(std::vector<std::size_t>::iterator it = data_centers.begin(); it != data_centers.end(); ++it, ++calculated_it)
        {
            if( (*it) == (*calculated_it) ) { ++true_count; }
        }
        double app_fit = static_cast<double>(true_count) / data_centers.size();
        std::cout << "[0,1]: " << app_fit << std::endl;
        if(app_fit < 0.7) {
            std::cerr << "There might be a problem if above printed comparison is too low." << std::endl;
            return EXIT_FAILURE;
        }
    }
}
