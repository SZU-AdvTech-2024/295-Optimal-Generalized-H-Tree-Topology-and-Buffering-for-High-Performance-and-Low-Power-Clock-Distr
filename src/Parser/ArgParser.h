#include <string>


class ArgParser {
public:
  enum Args {
    problemFile,
    constraintFile,
    outputFile
  };

  void parseCmd(int argc, char** argv);

  std::string get_problem_file()    { return _problem_file; }
  std::string get_constraint_file() { return _constraint_file; }
  std::string get_output_file()     { return _output_file; }


private:
  std::string _problem_file;
  std::string _constraint_file;
  std::string _output_file;
};
