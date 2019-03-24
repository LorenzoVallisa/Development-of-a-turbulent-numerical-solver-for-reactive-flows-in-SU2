#ifndef SU2_BUILDER_PROVIDER
#define SU2_BUILDER_PROVIDER

#include "concrete_provider.hpp"

#include <type_traits>
#include <memory>
#include "su2_assert.hpp"

namespace Common {

  /*!
    * \brief Concrete builder class. It forms the base for builders to be used with factories(one argument for constructor).
    * \author G. Orlando
  */
  template<class Abstract,class Concrete,int nargs = 1>
  class ProviderBuilder: public Abstract::Provider {
  public:
    static_assert(std::is_base_of<Abstract, Concrete>::value, "Builder requires Abstract to be a base of Concrete");

    /*!
     * \brief Default constructor
    */
    explicit ProviderBuilder(const std::string& name): Abstract::Provider(name) {}

    /*!
     * \brief Class destructor;
    */
    virtual ~ProviderBuilder() {}

    /*!
     * \brief Create objects of type Concrete
     * \param[in] arg - Argument to create the provider
    */
    inline std::unique_ptr<Abstract> Create(typename Abstract::Arg1 arg) override {
      return std::unique_ptr<Abstract>(new Concrete(arg));
    }

  }; /*-- End of class ProviderBuilder ---*/

  /*!
    * \brief Abstract builder class. It forms the base for builders to be used with factories(two arguments for constructor).
    * \author G. Orlando
  */
  template<class Abstract,class Concrete>
  class ProviderBuilder<Abstract,Concrete,2>: public Abstract::Provider {
  public:
    static_assert(std::is_base_of<Abstract, Concrete>::value, "Builder requires Abstract to be a base of Concrete");

    /*!
     * \brief Default constructor
    */
    explicit ProviderBuilder(const std::string& name): Abstract::Provider(name) {}

    /*!
     * \brief Class destructor;
    */
    virtual ~ProviderBuilder() {}

    /*!
     * \brief Create objects of type Concrete
     * \param[in] arg1 - First argument to create the provider
     * \param[in] arg2 - Second argument to create the provider
    */
    inline std::unique_ptr<Abstract> Create(typename Abstract::Arg1 arg1, typename Abstract::Arg2 arg2) override {
      return std::unique_ptr<Abstract>(new Concrete(arg1,arg2));
    }

  }; /*-- End of class ProviderBuilder ---*/

} /*-- End of namespace Common ---*/

#endif
