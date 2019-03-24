#ifndef SU2_NOTCOPYABLE
#define SU2_NOTCOPYABLE

namespace Common {

  /*!
   * \brief This class provides a clear and clean way to specify, through inheritance, that a certain object can't be copied
   * \author G. Orlando
  */

  template<class T>
  class NotCopyable {
  public:
    /*!
      * \brief Default constructor
      */
    NotCopyable() = default;

    /*!
      * \brief Virtual destructor
      */
    virtual ~NotCopyable() = default;

    /*!
      * \brief Deleted copy constructor to prevent copies
      */
    NotCopyable(const NotCopyable&) = delete;

    /*!
      * \brief Deleted assignment operator to prevent copy-assignment
      */
    NotCopyable& operator=(const NotCopyable&) = delete;

  }; /*-- End of class NotCopyable ---*/

} /*-- End of Namespace Common ---*/

#endif
