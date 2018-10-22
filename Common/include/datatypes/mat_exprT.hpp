#ifndef SU2_MATEXPRT_HPP
#define SU2_MATEXPRT_HPP

#include "exprT.hpp"

namespace Common {

  /**
    * Definition of a wrapper base expression class MatExprT which deruves from
    * ExprT. The expression template accepts two template parameters:
    * 1- the derived expression (the matrix)
    * 2- the return type of the epxression
    * \author G. Orlando
  */
  template<class Derived,typename T>
  struct MatExprT: public ExprT<Derived,T> {
    using Type = T;

    /*!
     * \brief Gets the number of rows
    */
    inline std::size_t nbRows() const {
      return this->operator const Derived&().nbRows();
    }

    /*!
     * \brief Gets the number of rows
    */
    inline std::size_t nbCols() const {
      return this->operator const Derived&().nbCols();
    }

    /*!
     * \brief Access to the individual entry (non const version)
    */
    inline T& operator()(std::size_t i, size_t j) {
      return this->operator Derived&().at(i,j);
    }

    /*!
     * \brief Access to the individual entry (const version)
    */
    //inline const T& operator()(std::size_t i, size_t j) const
    inline T operator()(std::size_t i, size_t j) const {
      return this->operator const Derived&().at(i,j);
    }


    /*!
     * \brief Access to the individual entry (non const version)
    */
    inline T& at(std::size_t i, size_t j) {
      return this->operator Derived&().at(i,j);
    }

    /*!
     * \brief Access to the individual entry (const version)
    */
    //inline const T& at(std::size_t i, size_t j) const
    inline T at(std::size_t i, size_t j) const {
      return this->operator const Derived&().at(i,j);
    }

  };

  /**
   * Definition of an expression template class for basic binary operations between matrixes
   * A macro is used to save code duplication.
   * The expression template accepts two parameters:
   * 1. left operand type
   * 2. right operand type.
   * \author G. Orlando
   */

  #define MET_TYPE(a) MatExprT<a,TYPE(a)>

  #define ETMAT_BINARY(OpName,__op__) \
  template <class Left, class Right>			\
  struct OpName: public MatExprT<OpName<Left,Right>, TYPE(Left)> { \
    static_assert(std::is_convertible<TYPE(Left),TYPE(Right)>::value,"The types of the two matrixes are not compatible"); \
                                                \
    OpName(const MET_TYPE(Left)& l, const MET_TYPE(Right)& r): e1(l), e2(r) { \
      SU2_Assert(e1.nbRows() == e2.nbRows(),"The number of rows of the two matrixes is not the same"); \
      SU2_Assert(e1.nbCols() == e2.nbCols(),"The number of columns of the two matrixes is not the same"); \
    }	\
                                          \
    inline TYPE(Left) operator()(std::size_t i, std::size_t j) const { \
      return e1(i,j) __op__ e2(i,j); \
    } \
                                      \
    inline TYPE(Left) at(std::size_t i, std::size_t j) const { \
      SU2_Assert(i < nbRows(),"Index is beyond the number of rows"); \
      SU2_Assert(j < nbCols(),"Index is beyond the number of columns"); \
      return e1.at(i,j) __op__ e2.at(i,j); \
    } \
    									               \
    std::size_t size() const { \
      return e1.size(); \
    }	\
                                        \
    std::size_t nbCols() const { \
      return e1.nbCols(); \
    } \
                                          \
    std::size_t nbRows() const { \
      return e1.nbRows(); \
    }  \
									                     \
  private: \
    const MET_TYPE(Left)& e1;	\
    const MET_TYPE(Right)& e2;	\
  };

  ETMAT_BINARY(MAdd, +)
  ETMAT_BINARY(MSub, -)

  #undef ETMAT_BINARY

  template <class Left, class Right>
  inline MAdd<Left,Right> operator +(const MET_TYPE(Left)& l, const MET_TYPE(Right)& r){return MAdd<Left,Right>(l,r);}

  template <class Left, class Right>
  inline MSub<Left,Right> operator -(const MET_TYPE(Left)& l, const MET_TYPE(Right)& r){return MSub<Left,Right>(l,r);}

  /**
   * Definition of an expression template class for basic binary operations between matrix and a double
   * A macro is used to save code duplication.
   * \author G. Orlando
   */
  #define ETMAT_BINARY_CR(OpName,__op__) \
  template <class Left>			\
  struct OpName: public MatExprT<OpName<Left>, TYPE(Left)> { \
    static_assert(std::is_convertible<TYPE(Left),double>::value,"The types of the matrix is not compatible with double"); \
                                                 \
    OpName (const MET_TYPE(Left)& l, const double& r): e(l), c(r) { \
      SU2_Assert(strcmp(#__op__ ,"/") && std::abs(c)>0,"You can't divide by zero"); \
    }	\
                                        \
    inline TYPE(Left) at(std::size_t i, std::size_t j) const { \
      SU2_Assert(i < nbRows(),"Index of row is beyond the number of rows"); \
      SU2_Assert(j < nbCols(),"Index of col is beyond the number of columns"); \
      return e.at(i,j) __op__ c; \
    } \
                                      \
    std::size_t size() const { \
      return e.size(); \
    }	\
                                         \
    std::size_t nbCols() const { \
      return e.nbCols(); \
    } \
                                           \
    std::size_t nbRows() const { \
       return e.nbRows(); \
    }  \
                                        \
  private: \
    const MET_TYPE(Left)& e;	\
    const double& c;	\
  };

  ETMAT_BINARY_CR(MAdd_CR, +)
  ETMAT_BINARY_CR(MSub_CR, -)
  ETMAT_BINARY_CR(MMul_CR, *)
  ETMAT_BINARY_CR(MDiv_CR, /)

  #undef ETMAT_BINARY_CR

  template <class Left>
  inline MAdd_CR<Left> operator +(const MET_TYPE(Left)& l, const double& r){return MAdd_CR<Left>(l,r);}

  template <class Left>
  inline MSub_CR<Left> operator -(const MET_TYPE(Left)& l, const double& r){return MSub_CR<Left>(l,r);}

  template <class Left>
  inline MMul_CR<Left> operator *(const MET_TYPE(Left)& l, const double& r){return MMul_CR<Left>(l,r);}

  template <class Left>
  inline MDiv_CR<Left> operator /(const MET_TYPE(Left)& l, const double& r){return MDiv_CR<Left>(l,r);}

  /**
   * Definition of an expression template class for basic binary operations between double and a matrix
   * A macro is used to save code duplication.
   * \author G. Orlando
   */
  #define ETMAT_BINARY_CL(OpName,__op__) \
  template <class Right>			\
  struct OpName: public MatExprT<OpName<Right>, TYPE(Right)> { \
    static_assert(std::is_convertible<double,TYPE(Right)>::value,"The types of the matrix is not compatible with double"); \
                                                 \
    OpName (const double& l,const MET_TYPE(Right)& r): c(l), e(r) {}	\
                                        \
    inline TYPE(Right) at(std::size_t i, std::size_t j) const { \
      SU2_Assert(i < nbRows(),"Index of row is beyond the number of rows"); \
      SU2_Assert(j < nbCols(),"Index of col is beyond the number of columns"); \
      return c __op__ e.at(i,j); \
    } \
                                      \
    std::size_t size() const { \
      return e.size(); \
    }	\
                                         \
    std::size_t nbCols() const { \
      return e.nbCols(); \
    } \
                                           \
    std::size_t nbRows() const { \
       return e.nbRows(); \
    }  \
                                        \
  private: \
    const double& c;	\
    const MET_TYPE(Right)& e;	\
  };

  ETMAT_BINARY_CL(MAdd_CL, +)
  ETMAT_BINARY_CL(MSub_CL, -)
  ETMAT_BINARY_CL(MMul_CL, *)
  //ETMAT_BINARY_CL(MDiv_CL, /)

  #undef ETMAT_BINARY_CL

  template <class Right>
  inline MAdd_CL<Right> operator +(const double& l,const MET_TYPE(Right)& r){return MAdd_CL<Right>(l,r);}

  template <class Right>
  inline MSub_CL<Right> operator -(const double& l,const MET_TYPE(Right)& r){return MSub_CL<Right>(l,r);}

  template <class Right>
  inline MMul_CL<Right> operator *(const double& l,const MET_TYPE(Right)& r){return MMul_CL<Right>(l,r);}

  //template <class Right>
  //inline MDiv_CL<Right> operator /(const double& l,const MET_TYPE(Right)& r){return MDiv_CL<Right>(l,r);}


  /**
    * Expression class for multiplication between two matrixes
    * \author G. Orlando
  */
  template <class Left,class Right>
  struct Mat_Mult :public MatExprT<Mat_Mult<Left,Right>, TYPE(Left)> {
    static_assert(std::is_convertible<TYPE(Left),TYPE(Right)>::value,"The types of the two matrixes are not compatible");

    /*!
     * \brief Class constructor
    */
    Mat_Mult(const MET_TYPE(Left)& l, const MET_TYPE(Right)& r): e1(l), e2(r) {
      SU2_Assert(e1.nbCols() == e2.nbRows(),"Matrix multiplication is not possible");
    }

    /*!
     * \brief Accessor to individual entry
    */
    TYPE(Left) operator()(std::size_t i, std::size_t j) const {
      TYPE(Left) res = TYPE(Left)();
      const std::size_t nbCols = e1.nbCols();
      for (std::size_t k = 0; k < nbCols; ++k)
        res += e1(i,k) * e2(k,j);
      return res;
     }

    /*!
     * \brief Accessor to individual entry
    */
    TYPE(Left) at(std::size_t i, std::size_t j) const {
      SU2_Assert(i < e1.nbRows(),"Index of row exceeds left matrix number of rows");
      SU2_Assert(j < e2.nbCols(),"Index of col exceeds right matrix number of columns");
      TYPE(Left) res = TYPE(Left)();
      const std::size_t nbCols = e1.nbCols();
      for (std::size_t k = 0; k < nbCols; ++k)
        res += e1.at(i,k) * e2.at(k,j);
      return res;
     }

    /*!
     * \brief Size of the resulting matrix
    */
    inline std::size_t size() const {
      return e1.nbRows()*e2.nbCols();
    }

    /*!
     * \brief Gets the number of columns of the resulting matrix
    */
    inline std::size_t nbCols() const {
      return e2.nbCols();
    }

    /*!
     * \brief Gets the number of rows of the resulting matrix
    */
    inline std::size_t nbRows() const {
      return e1.nbRows();
    }

  private:

    const MET_TYPE(Left)& e1;
    const MET_TYPE(Right)& e2;
  };

  template <class Left, class Right>
  inline Mat_Mult<Left,Right> operator*(const MET_TYPE(Left)& l, const MET_TYPE(Right)& r) {return Mat_Mult<Left,Right>(l,r);}

  /**
    * Expression class for matrix-vector multiplication
    * \author G. Orlando
  */
  template <class Mat, class Vec>
  struct MatVec_Mult :public MatExprT<MatVec_Mult<Mat,Vec>, TYPE(Vec)> {
    static_assert(std::is_convertible<TYPE(Mat),TYPE(Vec)>::value,"The types of matrix and vector are not compatible");

    /*!
     * \brief Class constructor
    */
    MatVec_Mult(const MET_TYPE(Mat)& l,const VET_TYPE(Vec)& r): e1(l), e2(r) {
      SU2_Assert(e1.nbCols() == e2.size(),"The dimension for matrix-vector multiplication is not compatible");
    }

    /*!
     * \brief Accessor to individual entry
    */
    TYPE(Vec) operator[](std::size_t i) const {
      auto res = TYPE(Vec)();
      for (std::size_t k = 0; k < e2.size(); ++k) {
        res += e1(i,k) * e2[k];
      }
      return res;
    }

    /*!
     * \brief Accessor to individual entry
    */
    TYPE(Vec) at(std::size_t i) const {
      SU2_Assert(i < size(),"The index exceed the dimension of the vector");
      auto res = TYPE(Vec)();
      for (std::size_t k = 0; k < e2.size(); ++k) {
        res += e1.at(i,k) * e2.at(k);
      }
      return res;
    }

    /*!
     * \brief Size of the resulting vector
    */
    inline std::size_t size() const {
      return e1.nbRows();
    }

  private:

    const MET_TYPE(Mat)& e1;
    const VET_TYPE(Vec)& e2;
  };

  template <class Mat, class Vec>
  inline MatVec_Mult<Mat,Vec> operator*(const MET_TYPE(Mat)& mat, const VET_TYPE(Vec)& vec) {return MatVec_Mult<Mat,Vec>(mat,vec);}

} /*--- End of namespace Common ---*/

#endif
