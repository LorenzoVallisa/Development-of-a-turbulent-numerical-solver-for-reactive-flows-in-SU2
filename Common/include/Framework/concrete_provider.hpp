#ifndef SU2_CONCRETE_PROVIDER
#define SU2_CONCRETE_PROVIDER

#include "provider.hpp"

namespace Common {

  /*!
    * \brief Concrete class for provider types (one argument for constructor)
    * \author G. Orlando
  */
  template<class Base,int nargs = 1>
  class ConcreteProvider: public Common::Provider<Base> {

  public:
    /*!
     * \brief Constructor of the class
    */
    explicit ConcreteProvider(const std::string& name): Common::Provider<Base>(name) {}

    /*!
     * \brief Virtual destructor
    */
    virtual ~ConcreteProvider() {}

    /*!
     *\brief Create an instance of provider: it must take exactly one argument
     *\param[in] arg - argument to construct the desired provider
    */
    virtual std::unique_ptr<Base> Create(typename Base::Arg1 arg) = 0;

  }; /*-- End of class ConcreteProvider ---*/

  /*!
    * \brief Concrete class for provider types (two arguments for constructor)
    * \author G. Orlando
  */
  template<class Base>
  class ConcreteProvider<Base,2>: public Common::Provider<Base> {

  public:
    /*!
     * \brief Constructor of the class
    */
    explicit ConcreteProvider(const std::string& name): Common::Provider<Base>(name) {}

    /*!
     * \brief Virtual destructor
    */
    virtual ~ConcreteProvider() {}

    /*!
     *\brief Create an instance of provider: it must take exactly two arguments
     *\param[in] arg1 - First argument to construct the desired provider
     *\param[in] arg2 - Second argument to construct the desired provider
    */
    virtual std::unique_ptr<Base> Create(typename Base::Arg1 arg1, typename Base::Arg2 arg2) = 0;

  }; /*-- End of class ConcreteProvider ---*/

} /*-- End of namespace Common ---*/

#endif
