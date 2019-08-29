#ifndef SU2_NOT_SETUP_EXCEPTION
#define SU2_NOT_SETUP_EXCEPTION

#include <exception>
#include <string>

/*! \class NotImplemented
 *  \brief Class for throwing an exception for non setup library.
 *  \author G. Orlando.
 */
namespace Common {
  class NotSetup: public std::exception {
  public:
    /*!
      * \brief Exception constructor
      * \param[in] what_arg - standard library string with error message.
      */
    explicit NotSetup(const std::string& what_arg): arg(what_arg) {}

    /*!
      * \brief Exception constructor
      * \param[in] what_arg - C-string with error message.
      */
    explicit NotSetup(const char* what_arg): arg(what_arg) {}

    /*!
      * \brief Exception message output
      */
    virtual const char* what() const throw() {
      return arg.c_str();
    }

    /*!
      * \brief Exception destructor
      */
    virtual ~NotSetup() {}

  protected:
    std::string arg; /*!< \brief String to store error message. */
  };
}

#endif
