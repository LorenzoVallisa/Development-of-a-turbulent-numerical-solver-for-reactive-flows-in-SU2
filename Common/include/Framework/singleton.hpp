#ifndef SU2_SINGLETON
#define SU2_SINGLETON

namespace Common {

  /*!
    * \brief This class provides a clear and clean way to specify, through inheritance, that a certain object can't be copied
    * \author G. Orlando
    */

  template<typename T>
  class Singleton {
  public:
    /*!
      * \brief Static initialization
      */
    static Singleton& GetInstance();

    /*!
      * \brief Virtual destructor
      */
    virtual ~Singleton() = default;

  private:

    /*!
      * \brief Private default constructor
      */
    Singleton() = default;

  public:
    /*!
      * \brief Deleted copy constructor to prevent copies
      */
    Singleton(const Singleton&) = delete;

    /*!
      * \brief Deleted assignment operator to prevent copy-assignment
      */
    Singleton& operator=(const Singleton&) = delete;

  }; /*-- End of class NonCopyable ---*/


template<typename T>
Singleton<T>& Singleton<T>::GetInstance() {
  static Singleton instance; /*! \brief Guaranteed to be destroyed and instantiated on first use */
  return instance;
}

} /*-- End of Namespace Common ---*/

#endif
