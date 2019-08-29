#ifndef __cplusplus
  #error You must use C++
#endif
#include <iostream>
#include <string>


#ifdef SU2_Assert
  #undef SU2_Assert
#endif

#ifndef ERRORSTATUS
  #define ERRORSTATUS 2
#endif

#ifndef NDEBUG
  #define SU2_Assert(Expr, Msg) \
          _My_Assert(#Expr, Expr, __FILE__, __LINE__, Msg)
#else
  #define SU2_Assert(Expr, Msg);
#endif

extern void _My_Assert(const std::string& expr_str, bool expr, const std::string& file, int line, const std::string& msg);
