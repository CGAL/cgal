// Copyright (c) 2014 Stefan Walk
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: LicenseRef-RFL
// License notice in Installation/LICENSE.RFL
//
// Author(s)     : Stefan Walk

// Modifications from original library:
//  * changed inclusion protection tag
//  * moved to namespace CGAL::internal::
//  * add parameter "reset_trees" to train() to be able to construct
//    forest with several iterations
//  * training algorithm has been parallelized with Intel TBB
//  * remove the unused feature "register_obb"
//  * add option to not count labels (if it's know before)
//  * fix the randomization of input (which was implicitly losing
//    samples)
//  * add method to get feature usage

#ifndef CGAL_INTERNAL_LIBLEARNING_RANDOMFOREST_FOREST_H
#define CGAL_INTERNAL_LIBLEARNING_RANDOMFOREST_FOREST_H
#include "common-libraries.hpp"
#include "tree.hpp"
#include <boost/ptr_container/serialize_ptr_vector.hpp>
#if VERBOSE_TREE_PROGRESS
#include <cstdio>
#endif

#include <CGAL/algorithm.h>
#include <CGAL/IO/binary_file_io.h>
#include <CGAL/tags.h>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/scalable_allocator.h>
#include <mutex>
#endif // CGAL_LINKED_WITH_TBB


namespace CGAL { namespace internal {



namespace liblearning {
namespace RandomForest {

template <typename NodeT, typename SplitGenerator>
class Tree_training_functor
{
  typedef typename NodeT::ParamType ParamType;
  typedef typename NodeT::FeatureType FeatureType;
  typedef Tree<NodeT> TreeType;

  std::size_t seed_start;
  const std::vector<int>& sample_idxes;
  boost::ptr_vector<Tree<NodeT> >& trees;
  DataView2D<FeatureType> samples;
  DataView2D<int> labels;
  std::size_t n_in_bag_samples;
  const SplitGenerator& split_generator;

public:

  Tree_training_functor(std::size_t seed_start,
                        const std::vector<int>& sample_idxes,
                        boost::ptr_vector<Tree<NodeT> >& trees,
                        DataView2D<FeatureType> samples,
                        DataView2D<int> labels,
                        std::size_t n_in_bag_samples,
                        const SplitGenerator& split_generator)
    : seed_start (seed_start)
    , sample_idxes (sample_idxes)
    , trees (trees)
    , samples (samples)
    , labels (labels)
    , n_in_bag_samples(n_in_bag_samples)
    , split_generator(split_generator)
  { }

#ifdef CGAL_LINKED_WITH_TBB
  void operator()(const tbb::blocked_range<std::size_t>& r) const
  {
    for (std::size_t s = r.begin(); s != r.end(); ++ s)
      apply(s);
  }
#endif // CGAL_LINKED_WITH_TBB

  inline void apply (std::size_t i_tree) const
  {
    // initialize random generator with sequential seeds (one for each
    // tree)
    RandomGen gen(seed_start + i_tree);
    std::vector<int> in_bag_samples = sample_idxes;

    // Bagging: draw random sample indexes used for this tree
    CGAL::cpp98::random_shuffle (in_bag_samples.begin(),in_bag_samples.end());

    // Train the tree
    trees[i_tree].train(samples, labels, &in_bag_samples[0], n_in_bag_samples, split_generator, gen);
  }

};

template <typename NodeT>
class RandomForest {
public:
    typedef typename NodeT::ParamType ParamType;
    typedef typename NodeT::FeatureType FeatureType;
    typedef Tree<NodeT> TreeType;
    ParamType params;

    boost::ptr_vector< Tree<NodeT> > trees;

    RandomForest() {}
    RandomForest(ParamType const& params) : params(params) {}

    template<typename ConcurrencyTag, typename SplitGenerator>
    void train(DataView2D<FeatureType> samples,
               DataView2D<int> labels,
               DataView2D<int> train_sample_idxes,
               SplitGenerator const& split_generator,
               size_t seed_start = 1,
               bool reset_trees = true,
               std::size_t n_classes = std::size_t(-1)
               )
    {
        if (reset_trees)
          trees.clear();

        if (n_classes == std::size_t(-1))
          params.n_classes = *std::max_element(&labels(0,0), &labels(0,0)+labels.num_elements()) + 1;
        else
          params.n_classes = n_classes;

        params.n_features = samples.cols;
        params.n_samples  = samples.rows;

        std::vector<int> sample_idxes;

        if (train_sample_idxes.empty()) {
            // no indexes were passed, generate vector with all indexes
            sample_idxes.resize(params.n_samples);
            for (size_t i_sample = 0; i_sample < params.n_samples; ++i_sample) {
                sample_idxes[i_sample] = i_sample;
            }
        } else {
            // copy indexes
            sample_idxes.assign(&train_sample_idxes(0,0), &train_sample_idxes(0,0)+train_sample_idxes.num_elements());
        }

        size_t n_idxes = sample_idxes.size();
        params.n_in_bag_samples = n_idxes * (1 - params.sample_reduction);

        std::size_t nb_trees = trees.size();
        for (std::size_t i_tree = nb_trees; i_tree < nb_trees + params.n_trees; ++ i_tree)
          trees.push_back (new TreeType(&params));

        Tree_training_functor<NodeT, SplitGenerator>
          f (seed_start, sample_idxes, trees, samples, labels, params.n_in_bag_samples, split_generator);

#ifndef CGAL_LINKED_WITH_TBB
        CGAL_static_assertion_msg (!(std::is_convertible<ConcurrencyTag, Parallel_tag>::value),
                                   "Parallel_tag is enabled but TBB is unavailable.");
#else
        if (std::is_convertible<ConcurrencyTag,Parallel_tag>::value)
        {
          tbb::parallel_for(tbb::blocked_range<size_t>(nb_trees, nb_trees + params.n_trees), f);
        }
        else
#endif
        {
          for (size_t i_tree = nb_trees; i_tree < nb_trees + params.n_trees; ++i_tree)
          {
#if VERBOSE_TREE_PROGRESS
            std::printf("Training tree %zu/%zu, max depth %zu\n", i_tree+1, nb_trees + params.n_trees, params.max_depth);
#endif
            f.apply(i_tree);
          }
        }
    }
    int evaluate(FeatureType const* sample, float* results) {
        // initialize output probabilities to 0
        std::fill_n(results, params.n_classes, 0.0f);
        // accumulate votes of the trees
        for (size_t i_tree = 0; i_tree < trees.size(); ++i_tree) {
            float const* tree_result = trees[i_tree].evaluate(sample);
            for (size_t i_cls = 0; i_cls < params.n_classes; ++i_cls) {
                results[i_cls] += tree_result[i_cls];
            }
        }
        float best_val   = 0.0;
        int   best_class = 0;
        float scale      = 1.0 / trees.size();
        for (size_t i_cls = 0; i_cls < params.n_classes; ++i_cls) {
            // divide by number of trees to normalize probability
            results[i_cls] *= scale;
            // determine best class
            if (results[i_cls] > best_val) {
                best_val = results[i_cls];
                best_class = i_cls;
            }
        }
        return best_class;
    }
#if 0
    float similarity_endnode(float const* sample_1, float const* sample_2) {
        double sum = 0.0;
        for (size_t i_tree = 0; i_tree < trees.size(); ++i_tree) {
            sum += trees[i_tree].similarity_endnode(sample_1, sample_2);
        }
        return sum/trees.size();
    }
    float similarity_path(float const* sample_1, float const* sample_2) {
        double sum = 0.0;
        for (size_t i_tree = 0; i_tree < trees.size(); ++i_tree) {
            sum += trees[i_tree].similarity_path(sample_1, sample_2);
        }
        return sum/trees.size();
    }
#endif
#if defined(CGAL_LINKED_WITH_BOOST_IOSTREAMS) && defined(CGAL_LINKED_WITH_BOOST_SERIALIZATION)
    template <typename Archive>
    void serialize(Archive& ar, unsigned /* version */)
    {
        ar & BOOST_SERIALIZATION_NVP(params);
        ar & BOOST_SERIALIZATION_NVP(trees);
    }
#endif

    void write (std::ostream& os)
    {
      params.write(os);

      I_Binary_write_size_t_into_uinteger32 (os, trees.size());
      for (std::size_t i_tree = 0; i_tree < trees.size(); ++i_tree)
        trees[i_tree].write(os);
    }

    void read (std::istream& is)
    {
      params.read(is);

      std::size_t nb_trees;
      I_Binary_read_size_t_from_uinteger32 (is, nb_trees);
      for (std::size_t i = 0; i < nb_trees; ++ i)
      {
        trees.push_back (new TreeType(&params));
        trees.back().read(is);
      }
    }

    void get_feature_usage (std::vector<std::size_t>& count) const
    {
      for (std::size_t i_tree = 0; i_tree < trees.size(); ++i_tree)
        trees[i_tree].get_feature_usage(count);
    }
};

}
}

}} // namespace CGAL::internal::

#endif
