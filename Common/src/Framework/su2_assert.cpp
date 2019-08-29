#include "../include/Framework/su2_assert.hpp"

void _My_Assert(const std::string& expr_str, bool expr, const std::string& file, int line, const std::string& msg) {
  if(!expr) {
    std::cerr << "Assert failed:\t" << msg << "\n"
              << "Expected:\t" << expr_str << "\n"
              << "Source:\t\t" << file << ", line " << line << "\n";
    std::exit(ERRORSTATUS);
  }
}
