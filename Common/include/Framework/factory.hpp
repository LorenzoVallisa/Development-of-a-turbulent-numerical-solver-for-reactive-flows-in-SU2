#ifndef SU2_FACTORY
#define SU2_FACTORY

#include "abstract_factory.hpp"

#include <map>
#include <stdexcept>
#include <memory>
#include <vector>
#include <sstream>


/*!
  * This namespace provides a factory class to load run-time the library that
  * will compute the physical and chemical properties of the considered mixture
*/

namespace Common {

  template<class Base> class AbstractProvider; /*!< \brief Forward declaration of the abstract provider */

  /*!
    * \class Factory
    * \brief Class for loading libraries at run-time.
    * \author G. Orlando
    */

  template<class Base>
  class Factory: public Common::AbstractFactory {
  public:

    /*
    * \brief Default destrcutor
    */
    ~Factory() = default;


    /*
    * \brief Checks if a provider is registered
    * \param[in] name - Name of the provider
    */
    bool Exists(const std::string& name);

    /*
    * \brief Registers a provider
    * \param[in] provider - Pointer to the provider to be registered
    */
    void Regist(std::unique_ptr<AbstractProvider<Base>> provider);

    /*
    * \brief Removes a registered provider (throw exception if it doesn't exist)
    * \param[in] provider_name - Name of the provider to be unregistered
    */
    void Unregist(const std::string& provider_name);

    /*
    * \brief Get the name of the Base class
    */
    inline const std::string GetBaseName() const {
      return Base::GetClassName();
    }

    /*
    * \brief Get a desired provider (throw exception if it doesn't exist)
    * \param[in] provider_name - Name of the provider to be unregistered
    */
    std::unique_ptr<typename Base::Provider> GetProvider(const std::string& provider_name) const;

    /*
    * \brief Returns all the providers in this factory
    */
    std::vector<std::unique_ptr<AbstractProvider<Base>>> GetAllProviders(void) const;


  private:

    typedef std::map<std::string,std::unique_ptr<AbstractProvider<Base>>> Container_type;

    Container_type database; /*!< \brief Database to store providers */

};


template<class Base>
bool Factory<Base>::Exists(const std::string& name) {
  return (database.count(name) > 0);
}

template<class Base>
void Factory<Base>::Regist(std::unique_ptr<AbstractProvider<Base>> provider) {
  if (Exists(provider->GetName())) {
    std::ostringstream converter;
    converter<<"In factory of [" << GetBaseName() <<
    "] a provider with the name [" << provider->GetProviderName() <<
    "] was found when trying to regist it\n";
    throw std::invalid_argument(converter.str());
  }
  database.insert(std::make_pair(provider->GetName(), provider));
}

template<class Base>
void Factory<Base>::Unregist(const std::string& provider_name) {
  if (!Exists(provider_name)) {
    std::string out="Provider " + provider_name + " is not stored in the factory";
	  throw std::invalid_argument(out);
  }
  database.erase(provider_name);
}


template<class Base>
std::unique_ptr<typename Base::Provider> Factory<Base>::GetProvider(const std::string& provider_name) const {
  if (!Exists(provider_name)) {
    std::string out="Provider " + provider_name + " is not stored in the factory";
	  throw std::invalid_argument(out);
  }
  return (database.find(provider_name)->second());
}

template<class Base>
std::vector<std::unique_ptr<AbstractProvider<Base>>> Factory<Base>::GetAllProviders(void) const {
  std::vector<std::unique_ptr<AbstractProvider<Base>>> res;
  res.reserve(database.size());
  for(auto i=database.begin(); i!=database.end();++i)
    res.push_back(i->second());
  return res;
}

}

#endif
