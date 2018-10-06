#ifndef SU2_NOT_SETUP_EXCEPTION
#define SU2_NOT_SETUP_EXCEPTION

#include <exception>
#include <string>

/*! \class NotImplemented
 *  \brief Class for throwing an exception for a non implemented function
 *  \author G. Orlando.
 */
namespace Common {
  class NotSetup: public std::exception {
  public:
    explicit NotSetup(const std::string& what_arg): arg(what_arg) {}
    explicit NotSetup(const char* what_arg): arg(what_arg) {}
    virtual const char* what() const throw() {
      return arg.c_str();
    }
    virtual ~NotSetup() {}

  protected:
    std::string arg;
  };
}

#endif
