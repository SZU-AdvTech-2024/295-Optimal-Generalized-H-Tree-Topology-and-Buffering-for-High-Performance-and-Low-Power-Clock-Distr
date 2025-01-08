#include "ArgParser.h"
#include <cassert>
#include <iostream>
#include <map>


void ArgParser::parseCmd(int argc, char** argv) {
  std::map<std::string, Args> arg_maps{
      {"-def_file", Args::problemFile},
      {"-constraint", Args::constraintFile},
      {"-output",Args::outputFile} };

  for (int i = 1; i < argc; ++i) {
    // if (!arg_maps.count(argv[i])) {
    //   assert(0);
    //   std::cout << "not support arg option " << argv[i];
    // }

    switch (arg_maps[argv[i]]) {
    case Args::problemFile: {
      _problem_file = argv[++i];  
      break;
    }
    case Args::constraintFile: {
      _constraint_file = argv[++i];  
      break;
    }
    case Args::outputFile: {
      _output_file = argv[++i];  
      break;
    }
    default:
      break;
    }
  }
}