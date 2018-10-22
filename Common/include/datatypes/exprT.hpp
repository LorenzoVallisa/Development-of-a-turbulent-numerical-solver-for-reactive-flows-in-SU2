#ifndef SU2_EXPRT
#define SU2_EXPRT

#include <iostream>
#include <type_traits>
#include <cmath>
#include <cstring>

#include "../su2_assert.hpp"

#ifdef __GNUG__
#define TF <>
#else
#define TF
#endif

namespace Common {

  /**
    * Definition of a wrapper base expression class ExprT: derived objects from ExprT use the CRTP technique
    * The expression template accepts two template parameters:
    * 1- the derived expression (the vector)
    * 2- the return type of the epxression
    * \author G. Orlando
  */

  template <class Derived, typename T>
  struct ExprT {

    using Type = T;

    /*!
      *\brief Cast operator to derived class (const version)
    */
    inline operator const Derived&() const {
      return static_cast<const Derived&>(*this);
    }

    /*!
     *\brief  Cast operator to derived class (non const version)
    */
    inline operator Derived&() {
      return static_cast<Derived&>(*this);
    }

   /*!
    * \brief Access operator (const version)
   */
   //inline const Type& operator[](std::size_t i) const
   inline Type operator[](std::size_t i) const {
     return this->operator const Derived&().operator[](i);
   }

   /*!
    * \brief Access operator (non const version)
   */
   inline Type& operator[](std::size_t i) {
     return this->operator Derived&().operator[](i);
   }

   /*!
    * \brief Access operator (const version)
   */
   inline Type& at(std::size_t i) {
     return this->operator Derived&().at(i);
   }

   /*!
    * \brief Access operator (non const version)
   */
   //inline const Type& at(std::size_t i) const
   inline Type at(std::size_t i) const {
     return this->operator const Derived&().at(i);
   }


   /*!
    * \brief Size of the derived object
   */
   std::size_t size() const {
     return this->operator const Derived&().size();
   }

}; /*--- End of class ExprT ----*/

  /**
   * Definition of an expression template class for basic binary operations between vectors
   * A macro is used to save code duplication.
   * The expression template accepts two parameters:
   * 1. left operand type
   * 2. right operand type.
   * \author G. Orlando
   */

   #define TYPE(a) typename a::Type
   #define VET_TYPE(a) ExprT<a,TYPE(a)>

   #define ETVEC_BINARY(OpName,__op__) \
   template <class Left, class Right>			\
   struct OpName : public ExprT<OpName<Left,Right>, TYPE(Left)> {	\
     static_assert(std::is_convertible<TYPE(Left),TYPE(Right)>::value,"The types of the two vectors are not convertible"); \
                              \
     OpName(const VET_TYPE(Left)& l, const VET_TYPE(Right)& r) :e1(l), e2(r) { \
       SU2_Assert(e1.size() == e2.size(),"The size of vectors is not the same"); \
     }	\
    	                         \
     inline TYPE(Left) operator[](std::size_t i) const { \
       return e1[i] __op__ e2[i];\
     } \
                                \
     inline TYPE(Left) at(std::size_t i) const { \
       SU2_Assert(i < size(),"Index is beyond the size of the vector");\
       return e1.at(i) __op__ e2.at(i);\
     } \
									\
    std::size_t size() const { \
      return e1.size();\
    }		\
                          \
   private:								\
      const VET_TYPE(Left)& e1;							\
      const VET_TYPE(Right)& e2;							\
  };

  ETVEC_BINARY(Add, +)
  ETVEC_BINARY(Sub, -)
  ETVEC_BINARY(Mul, *)
  //ETVEC_BINARY(Div, /)
  //ETVEC_BINARY(Max, std::max(e1.at(i),e2.at(i)))
  //ETVEC_BINARY(Min, std::min(e1.at(i),e2.at(i)))
  #undef ETVEC_BINARY

  template <class Left, class Right>
  inline Add<Left,Right> operator +(const VET_TYPE(Left)& l, const VET_TYPE(Right)& r){return Add<Left,Right>(l,r);}

  template <class Left, class Right>
  inline Sub<Left,Right> operator -(const VET_TYPE(Left)& l, const VET_TYPE(Right)& r){return Sub<Left,Right>(l,r);}

  template <class Left, class Right>
  inline Mul<Left,Right> operator *(const VET_TYPE(Left)& l, const VET_TYPE(Right)& r){return Mul<Left,Right>(l,r);}

  //template <class Left, class Right>
  //inline Div<Left,Right> operator /(const VET_TYPE(Left)& l, const VET_TYPE(Right)& r){return Div<Left,Right>(l,r);}

  /**
   * Definition of an expression template class for basic binary operations between vector and a constant
   * A macro is used to save code duplication.
   * \author G. Orlando
   */

  #define ETVEC_BINARY_CR(OpName,__op__) \
  template <class Left>			\
  struct OpName : public ExprT<OpName<Left>, TYPE(Left)> {	\
    static_assert(std::is_convertible<TYPE(Left),double>::value,"The type of the vector is not convertible to double"); \
                         \
    OpName(const VET_TYPE(Left)& l,const double& r): e(l),c(r) { \
      SU2_Assert(strcmp(#__op__ ,"/") && std::abs(c)>0,"You can't divide by zero"); \
    }	\
                              \
    inline TYPE(Left) operator[](std::size_t i) const { \
      return e[i] __op__ c;\
    } \
                              \
    inline TYPE(Left) at(std::size_t i) const { \
      SU2_Assert(i < size(),"Index is beyond the size of the vector");\
      return e.at(i) __op__ c;\
    } \
                 \
    std::size_t size() const { \
      return e.size();\
    }		\
                         \
   private:	\
     const VET_TYPE(Left)& e; \
     const double c; \
   };

   ETVEC_BINARY_CR(Add_CR, +)
   ETVEC_BINARY_CR(Sub_CR, -)
   ETVEC_BINARY_CR(Mul_CR, *)
   ETVEC_BINARY_CR(Div_CR, /)

   #undef ETVEC_BINARY_CR

   template <class Left>
   inline Add_CR<Left> operator +(const VET_TYPE(Left)& l,const double& r){return Add_CR<Left>(l,r);}

   template <class Left>
   inline Sub_CR<Left> operator -(const VET_TYPE(Left)& l,const double& r){return Sub_CR<Left>(l,r);}

   template <class Left>
   inline Mul_CR<Left> operator *(const VET_TYPE(Left)& l,const double& r){return Mul_CR<Left>(l,r);}

   template <class Left>
   inline Div_CR<Left> operator /(const VET_TYPE(Left)& l,const double& r){return Div_CR<Left>(l,r);}

   /**
    * Definition of an expression template class for basic binary operations between a constant and a vector
    * A macro is used to save code duplication.
    * \author G. Orlando
    */

   #define ETVEC_BINARY_CL(OpName,__op__) \
   template <class Right>			\
   struct OpName : public ExprT<OpName<Right>, TYPE(Right)> {	\
     static_assert(std::is_convertible<double,TYPE(Right)>::value,"The type of the vector is not convertible to double"); \
                          \
     OpName(const double& l,const VET_TYPE(Right)& r): c(l),e(r) { \
       /*SU2_Assert(strcmp(#__op__ ,"/") && std::all_of(reinterpret_cast<Right&>(*this).begin(),reinterpret_cast<Right&>(*this).end(),*/ \
      /*            [](const TYPE(Right)& x){return x!= TYPE(Right)();}),"The vector contains zero and so you can't compute its reciprocal");*/ \
     } \
     	                       \
     inline TYPE(Right) operator[](std::size_t i) const { \
       return c __op__ e[i];\
     } \
                               \
     inline TYPE(Right) at(std::size_t i) const { \
       SU2_Assert(i < size(),"Index is beyond the size of the vector");\
       return c __op__ e.at(i);\
     } \
                  \
     std::size_t size() const { \
       return e.size();\
     }		\
                          \
    private:	\
      const double c; \
      const VET_TYPE(Right)& e; \
    };

    ETVEC_BINARY_CL(Add_CL, +)
    ETVEC_BINARY_CL(Sub_CL, -)
    ETVEC_BINARY_CL(Mul_CL, *)
    //ETVEC_BINARY_CL(Div_CL, /)

    #undef ETVEC_BINARY_CL

    template <class Right>
    inline Add_CL<Right> operator +(const double& l,const VET_TYPE(Right)& r){return Add_CL<Right>(l,r);}

    template <class Right>
    inline Sub_CL<Right> operator -(const double& l,const VET_TYPE(Right)& r){return Sub_CL<Right>(l,r);}

    template <class Right>
    inline Mul_CL<Right> operator *(const double& l,const VET_TYPE(Right)& r){return Mul_CL<Right>(l,r);}

    //template <class Right>
    //inline Div_CL<Right> operator /(const double& l,const VET_TYPE(Right)& r){return Div_CL<Right>(l,r);}

  /**
   * Definition of an expression template class for basic unary operations for vectors.
   * A macro is used to save code duplication.
   * \author G. Orlando
   */
   #define ETVEC_UNARY(OpName,result) \
   template <class Operand>			\
   struct OpName : public ExprT<OpName<Operand>, TYPE(Operand)> { \
     OpName(const Operand& op) : e(op) {}	\
                                \
     inline const TYPE(Operand) operator[](std::size_t i) const { \
       return result(e[i]); \
     } \
     									              \
     inline const TYPE(Operand) at(std::size_t i) const {\
       SU2_Assert(i < size(), "Index is beyond the size of the vector"); \
       return result(e.at(i)); \
     } \
     									\
     std::size_t size() const { \
       return e.size();\
     }			\
  private:								\
     const Operand& e;							\
   };

   //ETVEC_UNARY(Cos,  std::cos)
   //ETVEC_UNARY(Sin,  std::sin)
   //ETVEC_UNARY(Tan, std::tan)
   ETVEC_UNARY(Log, std::log)
   ETVEC_UNARY(Log10, std::log10)
   ETVEC_UNARY(Exp,  std::exp)
   //ETVEC_UNARY(Sqrt, std::sqrt)
   //ETVEC_UNARY(Abs,  std::abs)
   ETVEC_UNARY(MinusOp, -)

   #undef ETVEC_UNARY

   //! Exponential
   template <class Operand>
   inline Exp<Operand> exp(const Operand& op){return Exp<Operand>(op);}

   //! Minus
   template <class Operand>
   inline MinusOp<Operand> operator-(const Operand& op){return MinusOp<Operand>(op);}

   //! Natural logarithm
   template <class Operand>
   inline Log<Operand> log(const Operand& op){return Log<Operand>(op);}

   //! 10 logarithm
   template <class Operand>
   inline Log10<Operand> log10(const Operand& op){return Log10<Operand>(op);}

}   /*--- End of namespace Common ----*/

#endif
