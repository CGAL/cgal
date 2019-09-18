#!/bin/bash

ack -l --cpp boost::lookup_named_param_def | xargs sed -i 's/boost::lookup_named_param_def/internal_np::Lookup_named_param_def/g'
ack -l --cpp boost::param_not_found | xargs sed -i 's/boost::param_not_found/internal_np::Param_not_found/g'
ack -l --cpp boost::is_default_param | xargs sed -i 's/boost::is_default_param/parameters::is_default_parameter/g'
ack -l --cpp boost::get_param | xargs sed -i 's/boost::get_param/parameters::get_parameter/g'
ack -l --cpp boost::choose_param | xargs sed -i 's/boost::choose_param/parameters::choose_parameter/g'
ack -l --cpp cgal_bgl_named_params | xargs sed -i 's/cgal_bgl_named_params/Named_function_parameters/g'


ack -l --cpp choose_param | xargs sed -i -E 's/choose_param[ ]*\(/choose_parameter(/g'
ack -l --cpp get_param | xargs sed -i -E 's/get_param[ ]*\(/get_parameter(/g'
ack -l --cpp is_default_param | xargs sed -i -E 's/is_default_param[ ]*\(/is_default_parameter(/g'

ack vertex_index_t | grep get_param

git checkout BGL/include/CGAL/boost/graph/dijkstra_shortest_paths.hpp
git checkout BGL/include/CGAL/boost/graph/Named_function_parameters.h
