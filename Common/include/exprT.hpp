#ifndef SU2_ExprT
#define SU2_ExprT

#include <iostream>

namespace Common {


template <typename T, int N, int M>
struct ETuple {
  typedef T Type;
  enum {Size1=N};
  enum {Size2=M};
};

/**
  * Definition of a wrapper base expression class ExprT: derived objects from ExprT use the CRTP technique
  * The expression template accepts two template parameters:
  * 1- the derived expression
  * 2- a tuple argument storing all arguments of the derived expression: these arguments
  *    arguments cannot be provided directly by 'Derived' (via typedefs or enum) because
  *    at the moment of instantiation, 'Derived' is still an incomplete type, being
  *    deriving from ExprT itself.
  * @author G. Orlando
 */

template <typename Derived, typename Tuple>
struct ExprT {

  typedef typename Tuple::Type Type;

  /*
   * Cast operator to derived class (const version)
   */
  operator const Derived&() const {
    return static_cast<const Derived&>(*this);
  }

  /*
   * Cast operator to derived class (const version)
   */
  operator Derived&() {
    return static_cast<Derived>(*this);
  }

  /*
   * Access operator (const version)
  */
  const Type& at(std::size_t i) const {
    return this->operator const Derived&().at(i);
  }

  /*
   * Access operator (non const version)
  */
  Type& at(std::size_t i) {
    return this->operator Derived&().at(i);
  }

  /*
   * Size of the derived object
  */
  std::size_t size() const {
    return this->operator const Derived&().size();
  }

}; /*--- End of class ExprT ----*/


/**
 * Definition of an expression template class for basic binary operations.
 * A macro is used to save code duplication.
 * The expression template accepts two parameters:
 * 1. first operand expression type
 * 2. second operand expression type.
 *
 * @author G. Orlando
 */

#define TUPLE_TYPE(a) typename a::Type
#define TUPLE_VEC(a,b) ETuple<TUPLE_TYPE(a),std::max(a::Size1, b::Size1),0>


#define ET_BINARY(OpName,operation) \
template <typename Left, typename Right>			\
class OpName : public ExprT<OpName<Left,Right>, TUPLE_VEC(Left,Right)> {	\
public:						\
  OpName(TUPLE_TYPE(Left) l, TUPLE_TYPE(Right) r) :	    \
    ExprT<OpName<Left,Right>, TUPLE_VEC(Left,Right)>(this), e1(l), e2(r) {}	\
    									\
  TUPLE_TYPE(Left) at(std::size_t i) const {
    return operation;\
  } \
									\
  std::size_t size() const {
    return Left.size();
  }		\
 private:								\
    TUPLE_TYPE(Left) e1;							\
    TUPLE_TYPE(Right) e2;							\
};

ET_BINARY(AddT, e1.at(i)+e2.at(i))
ET_BINARY(SubT, e1.at(i)-e2.at(i))
ET_BINARY(MulT, e1.at(i)*e2.at(i))
ET_BINARY(DivT, e1.at(i)/e2.at(i))
ET_BINARY(MaxT, std::max(e1.at(i),e2.at(i)))
ET_BINARY(MinT, std::min(e1.at(i),e2.at(i)))
#undef EET_BINARY

}   /*--- End of namespace Common ----*/

#endif
