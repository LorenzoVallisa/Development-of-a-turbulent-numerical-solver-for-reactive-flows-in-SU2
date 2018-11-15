#ifndef DEFAULT_INITIALIZATION_HPP
#define DEFAULT_INITIALIZATION_HPP

#include "index_sequence.hpp"

#include <tuple>
#include <iostream>

namespace Common {

  template<typename T, std::size_t... Is>
  auto repeat_impl(const T& x, Common::index_sequence<Is...>)->decltype(std::make_tuple((Is,x)...)) {
    return std::make_tuple((Is,x)...);
  }


  template <std::size_t N,typename T>
  auto repeat(const T& x)->decltype(repeat_impl(x, Common::make_index_sequence<N>())) {
    return repeat_impl(x, Common::make_index_sequence<N>());
  }

}

#endif
