#ifndef SU2_ABSTRACT_PROVIDER
#define SU2_ABSTRACT_PROVIDER

#include "factory.hpp"

#include <string>

namespace Common {

/*!
  * \brief Abstract class for provider types
  * \author G. Orlando
  */

  template<class Base>
  class AbstractProvider {

  public:

    /*!
      * \brief Explicit constructor
      */
    explicit AbstractProvider(const std::string& name):provider_name(name) {
      Common::Factory<Base>::GetInstance().Regist(this);
    }

    /*!
      * \brief Virtual destructor
      */
    virtual ~AbstractProvider() {}

    /*!
      * \brief Get the name of this provider
      */
    inline std::string GetProviderName(void) const {
      return provider_name;
    }

    /*!
      *\brief Free an instance created by the factory
      *\@param ptr pointer to be freed
      */
    virtual void FreeInstance(void* ptr) = 0;


  private:

    std::string provider_name;

}; /*-- End of class AbstractProvider ---*/

//////////////////////////////////////////////////////////////////////////////

} /*-- End of namespace Common ---*/

#endif
