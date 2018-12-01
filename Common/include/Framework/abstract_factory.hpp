#ifndef SU2_ABSTRACT_FACTORY
#define SU2_ABSTRACT_FACTORY

#include "singleton.hpp"

#include <string>

/*!
  * This namespace provides a factory class to load run-time the library that
  * will compute the physical and chemical properties of the considered mixture
*/

namespace Common {

/*!
 * \class AbstractFactory
 * \brief Interface for factory to load libraries at run-time.
 * \author G. Orlando
 */

 class AbstractFactory: public Common::Singleton<AbstractFactory> {
 public:

   /*
    * \brief Default destrcutor
    */
    virtual ~AbstractFactory() = default;

    /*
     * \brief Get the name of the Base class
     */

    virtual const std::string GetBaseName() const = 0;

};

}

#endif
