#ifndef INTEGER_SEQUENCE_HPP_
#define INTEGER_SEQUENCE_HPP_

#include <type_traits>

namespace Common {

  /*!
   * Definition of a class integer_sequence (implementation from https://en.cppreference.com/w/cpp/utility/integer_sequence)
   * \author G.Orlando
   */
  template <typename T, T... Ints>
  struct integer_sequence {
    static_assert(std::is_integral<T>::value,"integer_sequence struct is only for integral types");
    using value_type = T;
    std::size_t size() const {
      return sizeof...(Ints);
    }
  };

  template<std::size_t... Ints>
  using index_sequence = Common::integer_sequence<std::size_t,Ints...>;

  template<typename T, std::size_t N,T... Ints>
  struct make_integer_sequence : make_integer_sequence<T, N-1, N-1, Ints...> {};

  /*!
    * \brief Specialization for N = 0 of make_integer_sequence
  */
  template<typename T,T...Ints>
  struct make_integer_sequence<T,0,Ints...>: Common::integer_sequence<T,Ints...> {};

  template<std::size_t N>
  using make_index_sequence = Common::make_integer_sequence<std::size_t, N>;

  template<typename... T>
  using index_sequence_for = Common::make_index_sequence<sizeof...(T)>;

}

#endif
