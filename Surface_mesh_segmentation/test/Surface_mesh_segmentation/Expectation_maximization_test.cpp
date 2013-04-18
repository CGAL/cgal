
#include <CGAL/internal/Surface_mesh_segmentation/Expectation_maximization.h>

#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
/**
 * Generates sample points using a few gauissians.
 * Then applies gmm fitting on these generated points.
 * Provides a heuristic score for each gmm fitting result.
 *
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
        
        for(int i = 0; i < 300; ++i) { data.push_back(var_nor()); }
    }
    
    // calculate closest center (using above gauissians) for each generated points
    // we will compare it with gmm fitting results
    // also we might want to compute mixing coef for each center and select centers according to mixing_coef * prob(data)
    std::vector<int> data_centers;
    for(std::vector<double>::iterator it = data.begin(); it != data.end(); ++it)
    {
        int center_id = -1, center_counter = 0;;
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
    
    // apply gmm fitting clustering
    typedef CGAL::internal::Expectation_maximization E_M;
    std::vector<E_M> gmm_fitters;
    gmm_fitters.push_back(E_M(distributions.size(), data, E_M::PLUS_INITIALIZATION));
    gmm_fitters.push_back(E_M(distributions.size(), data, E_M::RANDOM_INITIALIZATION));
    gmm_fitters.push_back(E_M(distributions.size(), data, E_M::K_MEANS_INITIALIZATION));

    std::vector< std::vector<int> > calculated_centers(gmm_fitters.size());
    std::vector< std::vector<int> >::iterator calc_centers_it = calculated_centers.begin();
    for(std::vector<E_M>::iterator it = gmm_fitters.begin(); it != gmm_fitters.end(); ++it, ++calc_centers_it)
    {
        it->fill_with_center_ids(*calc_centers_it);
    }
    
    std::cout << "Compare results of EM with 'expected' (but be aware, it is not optimal result in terms of likelihood)" << std::endl;
    std::cout << "Another words a clustering which has higher likelihood can result in worse score in here" << std::endl;
    for(std::vector< std::vector<int> >::iterator calc_centers_it = calculated_centers.begin();
        calc_centers_it != calculated_centers.end(); ++calc_centers_it)
    {
        int true_count = 0;
        std::vector<int>::iterator calculated_it = calc_centers_it->begin();
        for(std::vector<int>::iterator it = data_centers.begin(); it != data_centers.end(); ++it, ++calculated_it)
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