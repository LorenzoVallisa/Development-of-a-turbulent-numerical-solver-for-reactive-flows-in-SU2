#ifndef SU2_MOVE_POINTER
#define SU2_MOVE_POINTER

namespace Common {
  template<typename T>
  std::unique_ptr<T[]> wrap_in_unique( T*&& t ) {
    auto* tmp = t;
    t = NULL;
    return std::unique_ptr<T[]>(tmp);
}

}

#endif
