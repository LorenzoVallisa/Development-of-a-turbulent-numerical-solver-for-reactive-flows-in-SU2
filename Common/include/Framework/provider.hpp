#ifndef SU2_PROVIDER
#define SU2_PROVIDER

#include "../not_copyable.hpp"

#include "factory.hpp"

namespace Common {

 /*!
  * \brief Class for provider types
  * \author G. Orlando
  */
  template<class Base>
  class Provider: public Common::NotCopyable<Provider<Base>> {

  public:
    /*!
     * \brief Explicit constructor
     */
    explicit Provider(const std::string& name): provider_name(name) {
      Common::Factory<Base>::GetInstance().Regist(this);
    }

    /*!
     * \brief Virtual destructor
     */
    virtual ~Provider() {}

    /*!
     * \brief Get the name of this provider
     */
    std::string GetProviderName(void) const {
      return provider_name;
    }

  private:

    std::string provider_name;

  }; /*-- End of class Provider ---*/

} /*-- End of namespace Common ---*/

#endif
