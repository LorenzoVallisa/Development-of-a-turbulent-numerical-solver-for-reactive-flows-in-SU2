#ifndef SU2_FACTORY
#define SU2_FACTORY

#include "../not_copyable.hpp"

#include <string>
#include <vector>

#include <map>
#include <stdexcept>
#include <memory>
#include <sstream>

/*!
 * This namespace provides a factory class to load run-time the library that
 * will compute the physical and chemical properties of the considered mixture
*/

namespace Common {

  template<class Base> class Provider;

  /*!
    * \class Factory
    * \brief Class for loading libraries at run-time.
    * \author G. Orlando
  */
  template<class Base>
  class Factory: public Common::NotCopyable<Factory<Base>> {
  public:

    /*
     * \brief Destrcutor
    */
    virtual ~Factory() {}

    /*
     * \brief Get the instance of this singleton
    */
    static Common::Factory<Base>& GetInstance(void);

    /*
     * \brief Checks if a provider is registered
     * \param[in] name - Name of the provider
    */
    bool Exists(const std::string& name) const;

    /*
     * \brief Registers a provider
     * \param[in] provider - Pointer to the provider to be registered
    */
    void Regist(Provider<Base>* provider);

    /*
     * \brief Removes a registered provider (throw exception if it doesn't exist)
     * \param[in] provider_name - Name of the provider to be unregistered
    */
    void Unregist(const std::string& provider_name);

    /*
     * \brief Get the name of the Base class
    */
    inline const std::string GetBaseName() const {
      return Base::GetBaseName();
    }

    /*
     * \brief Get a desired provider (throw exception if it doesn't exist)
     * \param[in] provider_name - Name of the provider to be unregistered
    */
    typename Base::Provider* GetProvider(const std::string& provider_name) const;

    /*
     * \brief Returns all the providers in this factory
    */
    std::vector<std::string> GetAllProviders(void) const;

  private:

    /*
     * \brief Constructor is private because this is a singleton
    */
    Factory() {}

  protected:

    typedef std::map<std::string,std::unique_ptr<Provider<Base>>> Container_type;

    Container_type database; /*!< \brief Database to store providers */

  }; /*--- End of class Factory ---*/

  //
  //
  /*
   * \brief Mayer's trick implementation
  */
  template<class Base>
  Factory<Base>& Factory<Base>::GetInstance(void) {
    static Common::Factory<Base> obj;
    return obj;
  }

  //
  //
  /*
   * \brief Check if a provider already exists
  */
  template<class Base>
  bool Factory<Base>::Exists(const std::string& name) const {
    return (database.count(name) > 0);
  }

  //
  //
  /*
   * \brief Regist a provider
  */
  template<class Base>
  void Factory<Base>::Regist(Provider<Base>* provider) {
    if(Exists(provider->GetProviderName())) {
      std::ostringstream converter;
      converter<<"In the factory [" << GetBaseName() <<
      "] a provider with the name [" << provider->GetProviderName() <<
      "] was found when trying to regist it\n";
      throw std::invalid_argument(converter.str());
    }
    database.emplace(provider->GetProviderName(), std::unique_ptr<Provider<Base>>(provider));
  }

  //
  //
  /*
   * \brief Unregist a provider
  */
  template<class Base>
  void Factory<Base>::Unregist(const std::string& provider_name) {
    if (!Exists(provider_name)) {
      std::string out="Provider " + provider_name + " is not stored in the factory";
	    throw std::invalid_argument(out);
    }
    database.erase(provider_name);
  }

  //
  //
  /*
   * \brief Get a pointer to a deisred provider
  */
  template<class Base>
  typename Base::Provider* Factory<Base>::GetProvider(const std::string& provider_name) const {
    if (!Exists(provider_name)) {
      std::string out="Provider " + provider_name + " is not stored in the factory";
	    throw std::invalid_argument(out);
    }
    return dynamic_cast<typename Base::Provider*>(database.find(provider_name)->second.get());
  }

  //
  //
  /*
   * \brief Get a list with all providers name
  */
  template<class Base>
  std::vector<std::string> Factory<Base>::GetAllProviders(void) const {
    std::vector<std::string> res;
    res.reserve(database.size());
    for(auto i = database.begin(); i != database.end(); ++i)
      res.push_back(i->first);
    return res;
  }

} /*--- End of namespace Common ---*/

#endif
