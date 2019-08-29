#ifndef SU2_FACTORY
#define SU2_FACTORY

#include "reacting_model_library.hpp"

#include <string>

#include <stdexcept>
#include <sstream>
#include <memory>

/*!
 * This namespace provides a factory class to load run-time the library that
 * will compute the physical and chemical properties of the considered mixture
*/

namespace Common {

  /*!
    * \class Factory
    * \brief Class for loading libraries at run-time.
    * \author G. Orlando
  */
  template<class Base>
  class Factory: public Common::NotCopyable<Factory<Base>> {
  public:

    /*
     * \brief Constructor of this simple factory
     * \param[in] lib_name - Name of the desired library
     * \param[in] config_name - Name of the file to read in order to configure the library
     * \param[in] lib_path - Path where the library is present
    */
    Factory(const std::string& lib_name, const std::string& config_name, const std::string& lib_path);

    /*
     * \brief Factory destrcutor
    */
    ~Factory() {}

    /*
     * \brief Get the library pointer and pass it to the solver shared pointer
    */
    std::shared_ptr<Base> GetLibraryPtr(void) const {
      return my_library;
    }

  private:
    std::shared_ptr<Base> my_library; /*!< \brief Pointer to Base in order to access concrete version. */

  }; /*--- End of class Factory ---*/

  //
  //
  /*--- Constrcutor implementation ---*/
  template<class Base>
  Factory<Base>::Factory(const std::string& lib_name, const std::string& config_name, const std::string& lib_path) {
    if(lib_name.compare("My_Library") == 0)
      my_library = std::make_shared<Framework::ReactingModelLibrary>(config_name, lib_path);
    else
      throw std::out_of_range("The library wanted is not present in the factory");
  }

} /*--- End of namespace Common ---*/

#endif
