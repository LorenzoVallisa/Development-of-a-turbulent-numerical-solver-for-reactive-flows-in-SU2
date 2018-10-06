#ifndef SU2_NOT_IMPLEMENTED_EXCEPTION
#define SU2_NOT_IMPLEMENTED_EXCEPTION

#include <exception>
#include <string>

/*! \class NotImplemented
 *  \brief Class for throwing an exception for a non implemented function
 *  \author G. Orlando.
 */
namespace Common {
  class NotImplemented: public std::exception {
  public:
    explicit NotImplemented(const std::string& what_arg): arg(what_arg) {}
    explicit NotImplemented(const char* what_arg): arg(what_arg) {}
    virtual const char* what() const throw() {
      return arg.c_str();
    }
    virtual ~NotImplemented() {}

  protected:
    std::string arg;
  };
}

#endif
