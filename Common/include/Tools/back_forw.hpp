#ifndef BACK_FORW_HPP
#define BACK_FORW_HPP

namespace MathTools {

  /*!
    * \brief Perform backward substition for an upper triangular matrix solving Ax = b
    * \param[in] A - Upper triangular matrix
    * \param[in] b - Known datum
  */
  template<class RealMatrix,class RealVec>
  RealVec Backward_Sub(const RealMatrix& A,const RealVec& b) {
    static_assert(std::is_arithmetic<typename RealMatrix::value_type>::value,
                  "The type of the matrix is not arithmetic. You can't use backward substition");
    static_assert(std::is_arithmetic<typename RealVec::value_type>::value,
                  "The type of the vector is not arithmetic. You can't use backward substition");

    std::size_t dim = A.nbRows();
    SU2_Assert(A.nbCols() == dim,"The matrix is not squared");
    SU2_Assert(b.size() == dim,"The dimension of known datum is not compatible with matrix");

    // We assume that the matrix is upper triangular without checking it
    RealVec x(dim);
    for(std::size_t j = 0; j < dim; ++j) {
      auto k = dim - j - 1;
      x[k] = b[k];
      for(std::size_t i = k + 1; i < dim; ++i)
        x[k] -= A(k,i)*x[i];
      x[k] /= A(k,k);
    }
    return x;
  }

  /*!
    * \brief Perform forward substition for a lower triangular matrix solving Ax = b
    * \param[in] A - Lower triangular matrix
    * \param[in] b - Known datum
  */
  template<class RealMatrix,class RealVec>
  RealVec Forward_Sub(const RealMatrix& A,const RealVec& b) {
    static_assert(std::is_arithmetic<typename RealMatrix::value_type>::value,
                  "The type of the matrix is not arithmetic. You can't use forward substition");
    static_assert(std::is_arithmetic<typename RealVec::value_type>::value,
                  "The type of the vector is not arithmetic. You can't use forward substition");

    std::size_t dim = A.nbRows();
    SU2_Assert(A.nbCols() == dim,"The matrix is not squared");
    SU2_Assert(b.size() == dim,"The dimension of known datum is not compatible with matrix");

    // We assume that the matrix is lower triangular without checking it
    RealVec x(dim);
    for(std::size_t j = dim; j > 0; --j) {
      auto k = dim - j;
      x[k] = b[k];
      for(int i = k - 1; i >= 0; --i)
        x[k] -= A(k,i)*x[i];
      x[k] /= A(k,k);
    }
    return x;
  }

} /*--- End of namespace Common ---*/

#endif
