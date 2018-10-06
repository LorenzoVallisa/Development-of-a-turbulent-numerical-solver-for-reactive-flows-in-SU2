#ifndef SU2_BUILDER_PROVIDER
#define SU2_BUILDER_PROVIDER

#include <type_traits>
#include <memory>

#include "concrete_provider.hpp"

namespace Common {

/*!
  * \brief Abstract builder class. It forms the base for builders to be used with factories.
  * \author G. Orlando
  */
template<class Abstract,class Concrete>
class ProviderBuilder {
 public:
  static_assert(std::is_base_of<Abstract, Concrete>::value,
  		         "Builder requires that Abstract be a base of Concrete");

  /*!
   * \brief Default constructor
   */
  ProviderBuilder() = default;

  /*!
   * \brief Class destructor;
   */
  ~ProviderBuilder() {}

  /*!
   * \brief Create objects of type Concrete with zero arguments;
   */
  std::unique_ptr<Abstract>  Create(typename Abstract::Arg1 arg) {
    return std::unique_ptr<Abstract>(new Concrete(arg));
  }

}; /*-- End of class ProviderBuilder ---*/

} /*-- End of namespace Common ---*/
#endif
