#ifndef CGAL_KSR_TERMINAL_PARSER_H
#define CGAL_KSR_TERMINAL_PARSER_H

#if defined(WIN32) || defined(_WIN32)
#define _SR_ "\\"
#else
#define _SR_ "/"
#endif

// STL includes.
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <unordered_map>

namespace CGAL {
namespace KSR {

  template<typename FT>
  class Terminal_parser {

  public:
    using Input_parameters = std::unordered_map<std::string, std::string>;

    Terminal_parser(
      const int num_parameters,
      const char** parameters,
      const std::string path_to_save = "") :
    m_path_to_save(path_to_save) {

      // Help.
      show_help(num_parameters, parameters);

      // Handle all input parameters.
      Input_parameters input_parameters;
      set_input_parameters(num_parameters, parameters, input_parameters);

      // Set here all required parameters.
      std::vector<std::string> required(1);
      required[0] = "-data";

      // Set parameters.
      set_parameters(input_parameters, required);

      // Set here all parameters that should not be saved.
      std::vector<std::string> exceptions(2);
      exceptions[0] = "-data";
      exceptions[1] = "-params";

      // Save parameters.
      save_parameters_to_file(input_parameters, exceptions);
    }

    inline Input_parameters& get_input_parameters() {
      return m_parameters;
    }

    inline const Input_parameters& get_input_parameters() const {
      return m_parameters;
    }

    template<typename Scalar>
    void add_val_parameter(
      const std::string parameter_name,
      Scalar& variable_value) {

      if (!does_parameter_exist(parameter_name, m_parameters))
        return;

      const std::string parameter_value = m_parameters.at(parameter_name);
      if (parameter_value != "default")
        variable_value = static_cast<Scalar>(std::stod(parameter_value.c_str()));
      std::cout << parameter_name << " : " << variable_value << std::endl;
    }

    void add_str_parameter(
      const std::string parameter_name,
      std::string& variable_value) {

      if (!does_parameter_exist(parameter_name, m_parameters))
        return;

      const std::string parameter_value = m_parameters.at(parameter_name);
      if (parameter_value != "default")
        variable_value = parameter_value;
      std::cout << parameter_name << " : " << variable_value << std::endl;
    }

    void add_bool_parameter(
      const std::string parameter_name,
      bool& variable_value) {

      if (!does_parameter_exist(parameter_name, m_parameters))
        return;

      variable_value = true;
      std::cout << parameter_name << " : " << (variable_value ? "true" : "false") << std::endl;
    }

  private:
    const std::string m_path_to_save;
    Input_parameters  m_parameters;

    // Help.
    void show_help(
      const int num_parameters,
      const char** parameters) {

      if (!is_asked_for_help(num_parameters, parameters))
        return;

      print_help();
      exit(EXIT_SUCCESS);
    }

    bool is_asked_for_help(
      const int num_parameters,
      const char** parameters) {

      for (int i = 0; i < num_parameters; ++i)
        if (std::strcmp(parameters[i], "-help") == 0)
          return true;
      return false;
    }

    void print_help() {

      std::cout << std::endl << "* HELP:" << std::endl;

      std::cout << std::endl << "* EXAMPLE:" << std::endl;
      std::cout <<
        "your terminal name $ ."
      << std::string(_SR_) <<
        "kinetic_reconstruction_example -data path_to_data"
      << std::string(_SR_) <<
        "data_name.ply -other_param_name -other_param_value"
      << std::endl << std::endl;

      std::cout << std::endl << "REQUIRED PARAMETERS:" << std::endl << std::endl;

      std::cout <<
        "parameter name: -data" << std::endl <<
        "parameter value: path_to_data" << std::string(_SR_) << "data_name.ply" << std::endl <<
        "description: path to the file with input data" << std::endl << std::endl;

      std::cout << std::endl << "OPTIONAL PARAMETERS:" << std::endl << std::endl;

      std::cout <<
        "parameter name: -silent" << std::endl <<
        "description: supress any intermediate output except for the final result" << std::endl << std::endl;

      std::cout <<
        "parameter name: -params" << std::endl <<
        "parameter value: path_to" << std::string(_SR_) << "parameters.ksr" << std::endl <<
        "description: load parameters from the file" << std::endl << std::endl;
    }

    // Setting parameters.
    void set_input_parameters(
      const int num_parameters,
      const char** parameters,
      Input_parameters& input_parameters) {

      assert(num_parameters > 0);
      for (int i = 1; i < num_parameters; ++i) {

        std::string str   = static_cast<std::string>(parameters[i]);
        auto first_letter = str[0];

        if (first_letter == '-') {
          if (i + 1 < num_parameters) {

            str = static_cast<std::string>(parameters[i + 1]);
            first_letter = str[0];

            if (first_letter != '-')
              input_parameters[parameters[i]] = parameters[i + 1];
            else
              input_parameters[parameters[i]] = "default";

          } else input_parameters[parameters[i]] = "default";
        }
      }
    }

    void set_parameters(
      Input_parameters& input_parameters,
      const std::vector<std::string>& /* required */) {
      if (parameters_should_be_loaded(input_parameters))
        load_parameters_from_file(input_parameters);
      m_parameters = input_parameters;
    }

    bool are_required_parameters_set(
      const Input_parameters& input_parameters,
      const std::vector<std::string>& required) {

      bool are_all_set = true;
      for (std::size_t i = 0; i < required.size(); ++i)
        if (!is_required_parameter_set(required[i], input_parameters))
          are_all_set = false;
      return are_all_set;
    }

    bool is_required_parameter_set(
      const std::string parameter_name,
      const Input_parameters& input_parameters) {

      const bool is_set = does_parameter_exist(parameter_name, input_parameters) &&
      !does_parameter_have_default_value(parameter_name, input_parameters);

      if (!is_set)
        std::cerr << std::endl <<
          parameter_name << " PARAMETER IS REQUIRED!"
        << std::endl;
      return is_set;
    }

    bool does_parameter_exist(
      const std::string parameter_name,
      const Input_parameters& input_parameters) {

      for (Input_parameters::const_iterator parameter = input_parameters.begin();
      parameter != input_parameters.end(); ++parameter)
        if ((*parameter).first == parameter_name)
          return true;
      return false;
    }

    bool does_parameter_have_default_value(
      const std::string parameter_name,
      const Input_parameters& input_parameters) {
      assert(does_parameter_exist(parameter_name, input_parameters));
      return input_parameters.at(parameter_name) == "default";
    }

    // Loading from a file.
    bool parameters_should_be_loaded(const Input_parameters& input_parameters) {
      if (does_parameter_exist("-params", input_parameters))
        return true;
      return false;
    }

    void load_parameters_from_file(Input_parameters& input_parameters) {
      const std::string filePath = input_parameters.at("-params");
      if (filePath == "default") {

        std::cerr << std::endl <<
          "ERROR: PATH TO THE FILE WITH PARAMETERS IS NOT DEFINED!"
        << std::endl << std::endl;
        exit(EXIT_FAILURE);
      }

      std::ifstream file(filePath.c_str(), std::ios_base::in);
      if (!file) {
        std::cerr << std::endl <<
          "ERROR: ERROR LOADING FILE WITH PARAMETERS!"
        << std::endl << std::endl;
        exit(EXIT_FAILURE);
      }

      Input_parameters tmp_parameters;
      while (!file.eof()) {
        std::string parameter_name, parameter_value;
        file >> parameter_name >> parameter_value;

        if (parameter_name == "" || parameter_value == "")
          continue;
        tmp_parameters[parameter_name] = parameter_value;
      }

      for (Input_parameters::const_iterator pit = tmp_parameters.begin();
      pit != tmp_parameters.end(); ++pit)
        input_parameters[(*pit).first] = (*pit).second;
      file.close();
    }

    // Saving to a file.
    void save_parameters_to_file(
      const Input_parameters& input_parameters,
      const std::vector<std::string>& exceptions) {

      save_input_parameters(m_path_to_save, input_parameters, exceptions);
      return;
    }

    void save_input_parameters(
      const std::string path_to_save,
      const Input_parameters& input_parameters,
      const std::vector<std::string>& exceptions) {

      const std::string file_path = path_to_save + "parameters.ksr";
      save_parameters(file_path, input_parameters, exceptions);
      std::cout << "* parameters are saved in: " << file_path << std::endl;
    }

    void save_parameters(
      const std::string file_path,
      const Input_parameters& input_parameters,
      const std::vector<std::string>& exceptions) {

      std::ofstream file(file_path.c_str(), std::ios_base::out);
      if (!file) {
        std::cerr << std::endl <<
          "ERROR: SAVING FILE WITH THE NAME " << file_path
        << std::endl << std::endl;
        exit(EXIT_FAILURE);
      }

      for (Input_parameters::const_iterator parameter = input_parameters.begin();
      parameter != input_parameters.end(); ++parameter)
        if (parameter_should_be_saved((*parameter).first, exceptions))
          file << (*parameter).first << " " << (*parameter).second << std::endl;
      file.close();
    }

    bool parameter_should_be_saved(
      const std::string parameter_name,
      const std::vector<std::string>& exceptions) {
      for (std::size_t i = 0; i < exceptions.size(); ++i)
        if (exceptions[i] == parameter_name)
          return false;
      return true;
    }
  };

} // KSR
} // CGAL

#endif // CGAL_KSR_TERMINAL_PARSER_H
