#ifndef SU2_MATRIXT_HPP
#define SU2_MATRIXT_HPP

#include "mat_exprT.hpp"
#include "../Tools/back_forw.hpp"
#include "../option_structure.hpp"

#include <algorithm>
#include <iomanip>
#include <utility>

namespace Common {

  /*!
   * Declaration of output stream operator
  */
  template<class Derived,typename T>
  std::ostream& operator<<(std::ostream& out, const MatExprT<Derived,T>& A);

  /*!
   * Definition of a class SU2Mat that implements an expression template technique
   * \author G.Orlando
  */

  template <typename T>
  class SU2Mat : public MatExprT<SU2Mat<T>,T> {
  public:

    using value_type = T;
    using size_type  = std::size_t;

    /*!
      * \brief Default constructor
    */
    SU2Mat(): n_rows(), n_cols(), m_data(NULL) {}

    /*!
      * \brief Class constructor
      * \param[in] _n_rows - number of rows in the matrix
      * \param[in] _n_cols - number of columns in the matrix
      * \param[in] init - value to initialize (default provided)
    */
    SU2Mat(size_type _n_rows, size_type _n_cols, value_type init = value_type());

    /*!
      * \brief Class destructor
    */
    ~SU2Mat() {
      free_mem();
    }

    /*!
      * \brief Copy constructor
      * \param[in] SU2Mat - SU2Mat to copy from
    */
    SU2Mat(const SU2Mat& SU2Mat);

    /*!
      * \brief Move constructor
      * \param[in] SU2Mat - SU2Mat to move
    */
    SU2Mat(SU2Mat&& SU2Mat);

    /*!
      * \brief Copy assignment operator
      * \param[in] SU2Mat - SU2Mat to copy from
    */
    SU2Mat& operator=(const SU2Mat& SU2Mat);

    /*!
      * \brief Move assignment operator
      * \param[in] SU2Mat - SU2Mat to move
    */
    SU2Mat& operator=(SU2Mat&& SU2Mat);

    /*!
      * \brief Copy Constructor from already allocated memory
      * \param[in] data - pointer to pointer to value_type (stored elements)
      * \param[in] _n_rows - number of rows of the resulting matrix
      * \param[in] _n_cols - number of columns of the resulting matrix
    */
    SU2Mat(value_type** data, size_type _n_rows, size_type _n_cols);

    /*!
      * \brief Copy Constructor from standard library vector (matrix as vector of vectors)
      * \param[in] v - std::vector to copy from
    */
    explicit SU2Mat(const std::vector<std::vector<value_type>>& mat);

    /*!
      * \brief Copy constructor from MatExprT;
      * \param[in] e - MatExprT to copy from
     */
    template <class Derived>
    SU2Mat(const MatExprT<Derived,value_type>& e);

    /*!
      * \brief Assignment operator from MatExprT;
      * \param[in] e - MatExprT to copy from
     */
    template <class Derived>
    SU2Mat& operator=(const MatExprT<Derived,value_type>& e);

    /*!
      * \brief Overloaded assignment operator with constants (macro used);
    */
    #define SU2MAT_ASSIGN_OP_CONST(__op__) \
    SU2Mat& operator __op__ (value_type value) { \
      static_assert(std::is_arithmetic<value_type>::value,"The type of the matrix is not arithmetic"); \
                                                    \
      SU2_Assert((strcmp(#__op__ ,"/=") == 0 && std::abs(value) > EPS) || (strcmp(#__op__,"/=") != 0),"You can't divide by zero"); \
      for(size_type i = 0; i < size(); ++i) \
        m_data[i] __op__ value; \
      return *this; \
    }
    SU2MAT_ASSIGN_OP_CONST(+=)
    SU2MAT_ASSIGN_OP_CONST(-=)
    SU2MAT_ASSIGN_OP_CONST(*=)
    SU2MAT_ASSIGN_OP_CONST(/=)
    #undef SU2MAT_EASSIGN_OP_CONST

    /*!
      * \brief Overloaded assignment operator with MatExprT (macro used);
    */
    #define SU2MAT_ASSIGN_OP_EXPR(__op__) \
    template<class Derived> \
    SU2Mat& operator __op__ (const MatExprT<Derived,value_type>& e) { \
      static_assert(std::is_arithmetic<value_type>::value,"The type of the matrix is not arithmetic"); \
                                                    \
      SU2_Assert(this->nbRows() == e.nbRows() && this->nbCols() == e.nbCols(), \
                "The dimension of the matrixes is not compatible and so you can't operate with them"); \
      for(size_type i = 0; i < nbRows(); ++i) { \
        for(size_type j = 0; j < nbCols(); ++j) \
          (*this)(i,j) __op__ e(i,j); \
      } \
      return *this; \
    }
    SU2MAT_ASSIGN_OP_EXPR(+=)
    SU2MAT_ASSIGN_OP_EXPR(-=)
    #undef SU2MAT_ASSIGN_OP_EXPR

    /*!
      * \brief Overloaded *= operator with MatExprT (matrix-matrix multiplication);
    */
    template<class Derived>
    SU2Mat& operator*=(const MatExprT<Derived,value_type>& e) {
      static_assert(std::is_arithmetic<value_type>::value,"The type of the matrix is not arithmetic");

      SU2_Assert(this->nbCols() == e.nbRows(), "The dimension of the matrixes is not compatible and so you can't multiply them");
      std::vector<value_type> tmp(m_data, m_data + size());
      n_cols = e.nbCols();
      this->resize(n_rows,n_cols);
      for(size_type i = 0; i < n_rows; ++i) {
        for(size_type j = 0; j < n_cols; ++j) {
          for(size_type k = 0; k < e.nbRows(); ++k)
              (*this)(i,j) += tmp[sub2ind(i,k)] * e(k,j);
        }
      }
      return *this;
    }

    /*!
     * \brief Returns the (i,j)-th element (non const version)
     * \param[in] i - index of row
     * \param[in] j - index of columns
    */
    inline value_type& operator()(size_type i, size_type j) {
      return m_data[sub2ind(i,j)];
    }

    /*!
     * \brief Returns the (i,j)-th element (const version)
     * \param[in] i - index of row
     * \param[in] j - index of columns
    */
    inline const value_type& operator()(size_type i, size_type j) const {
      return m_data[sub2ind(i,j)];
    }

    /*!
     * \brief Returns the (i,j)-th element (non const version)
     * \param[in] i - index of row
     * \param[in] j - index of columns
    */
    inline value_type& at(size_type i, size_type j) {
      SU2_Assert(i < n_rows,"Index of row is beyond the number of rows in the matrix");
      SU2_Assert(j < n_cols,"Index of col is beyond the number of columns in the matrix");
      return m_data[sub2ind(i,j)];
    }

    /*!
     * \brief Returns the (i,j)-th element (const version)
     * \param[in] i - index of row
     * \param[in] j - index of columns
    */
    inline const value_type& at(size_type i, size_type j) const {
      SU2_Assert(i < n_rows,"Index of row is beyond the number of rows in the matrix");
      SU2_Assert(j < n_cols,"Index of col is beyond the number of columns in the matrix");
      return m_data[sub2ind(i,j)];
    }

    /*!
     * \brief Returns the number of rows in the matrix
    */
    inline size_type nbRows(void) const {
      return n_rows;
    }

    /*!
     * \brief Returns the number of rows in the matrix
    */
    inline size_type nbCols(void) const {
      return n_cols;
    }

    /*!
     * \brief Returns the size of the matrix
    */
    inline size_type size(void) const {
      return n_rows*n_cols;
    }

    /*!
     * \brief Clears the matrix
    */
    inline void clear(void) {
      free_mem();
      n_rows = 0;
      n_cols = 0;
    }

    /*!
     * \brief Swap the content of two matrixes
     * \param[in] su2mat - Matrix to swap with the current one
    */
    void swap(SU2Mat& su2mat);

    /*!
     * \brief Resize the matrix
     * \param[in] _n_rows - new number of rows in the matrix
     * \param[in] _n_cols - new number of columns in the matrix
     * \param[in] init - value to initialize (default provided)
    */
    void resize(size_type _n_rows, size_type _n_cols, value_type init = value_type());

    /*!
     * \brief Returns the underlined pointer
    */
    inline value_type* data(void) const {
      return m_data;
    }

    /*!
     * \brief Compute the determinant of a matrix
    */
    value_type determinant(void) const;

    /*!
     * \brief Compute the inverse of a matrix
    */
    SU2Mat inverse(void) const;

    /*!
     * \brief Compute the transpose of the matrix
    */
    SU2Mat transpose(void) const;

    /*!
     * \brief Returns the LU decomposition of the matrix
    */
    std::pair<SU2Mat,SU2Mat> LU_decomp(void) const;

    /*!
     * \brief Check if the matrix is diagonal
    */
    bool check_diag(void) const;

    /*!
     * \brief Compute the invers of a diagonal matrix
    */
    SU2Mat invert_diag(void) const;

    /*!
     * \brief Puts the copy of the desired row of the matrix in a supplied vector
     * \param[in] idx_row - index of desired row
     * \param[out] vec - the resulting vector
    */
    template <class Vec>
    Vec GetRow(size_type idx_row) const;

    /*!
     * \brief Puts the copy of the desired column of the matrix in a supplied vector
     * \param[in] idx_col - index of desired row
     * \param[out] vec - the resulting vector
    */
    template <class Vec>
    Vec GetColumn(size_type idx_col) const;

    /*!
     * \brief Sets the desired row of the matrix from a vector
     * \param[in] idx_row - index of desired row
     * \param[in] row - the vector where data comes from
    */
    template <class Vec>
    void SetRow(size_type idx_row, const Vec& row);

    /*!
     * \brief Sets the desired column of the matrix from a vector
     * \param[in] idx_col - index of desired column
     * \param[in] col - the vector where data comes from
    */
    template <class Vec>
    void SetColumn(size_type idx_col, const Vec& col);

    /*!
     * \brief Cast to standard library vector of vectors ("standard matrix")
    */
    operator std::vector<std::vector<value_type>>() const {
      auto result = std::vector<value_type>(n_rows,std::vector<value_type>(n_cols));
      for(size_type i = 0; i < n_rows; ++i) {
        for(size_type j = 0; j < n_cols; ++j)
          result[i][j] = m_data[sub2ind(i,j)];
      }
      return result;
    }

    /*!
     * \brief Begin iterator (non const version)
    */
    inline auto begin()->decltype(std::declval<std::vector<value_type>>().begin()) {
      return static_cast<decltype(std::declval<std::vector<value_type>>().begin())>(m_data);
    }

    /*!
     * \brief End iterator (non const version)
    */
    inline auto end()->decltype(std::declval<std::vector<value_type>>().end()) {
      return static_cast<decltype(std::declval<std::vector<value_type>>().end())>(m_data + n_rows*n_cols);
    }

    /*!
     * \brief Begin iterator (const version)
    */
    inline decltype(std::declval<std::vector<value_type>>().cbegin()) cbegin() const {
      return static_cast<decltype(std::declval<std::vector<value_type>>().cbegin())>(m_data);
    }

    /*!
     * \brief End iterator (const version)
    */
    inline decltype(std::declval<std::vector<value_type>>().cend()) cend() const {
      return static_cast<decltype(std::declval<std::vector<value_type>>().cend())>(m_data + n_rows*n_cols);
    }

    /*!
     * \brief Overloading of output stream operator
    */
    friend std::ostream& operator<< TF (std::ostream& out,const MatExprT<SU2Mat,value_type>& A);

    /*!
     * \brief Check if there are zeros in the matrix
    */
    inline bool check_notzero(void) const {
      return std::none_of(this->cbegin(), this->cend(), [&](const value_type& x){return x == value_type();});
    }

  private:

    size_type n_rows;  /*!< \brief Number of rows of the matrix. */
    size_type n_cols;  /*!< \brief Number of columns of the matrix. */
    value_type* m_data;  /*!< \brief Stored elements. */

    /*!
     * \brief Helper function to allocate memory
    */
    void allocate(void);

    /*!
     * \brief Helper function to initialize vector
     * \param[in] value - Value initialization
    */
    void initialize(value_type value);

    /*!
     * \brief Helper function to free memory
    */
    void free_mem(void);

    /*!
     * \brief sub2ind
     * \param[in] i - row index
     * \param[in] j - column index
    */
    inline size_type sub2ind(size_type i, size_type j) const {
      return i*n_cols + j;
    }

    /*!
     * \brief Compute the cofactor of a specific element of a matrix
     * \param[in] idx_row - Index of row
     * \param[in] idx_col - Index of column
    */
    SU2Mat<T> GetCoFactor(size_type idx_row, size_type idx_col) const;

    /*!
     * \brief Compute the adjoint of a matrix
    */
    SU2Mat<T> adjoint(void) const;


  }; /*--- End of class declaration SU2Mat ---*/

  //
  //
  /*--- Allocate memory---*/
  template<typename T>
  void SU2Mat<T>::allocate(void) {
    if(size() > 0)
      m_data = new value_type[size()];
  }

  //
  //
  /*--- Initialize matrix with a value ---*/
  template<typename T>
  void SU2Mat<T>::initialize(value_type value) {
    std::fill(m_data, m_data + size(), value);
  }

  //
  //
  /*--- Free memory ---*/
  template<typename T>
  void SU2Mat<T>::free_mem(void) {
    if(m_data != NULL) {
      delete [] m_data;
      m_data = NULL;
    }
  }

  //
  //
  /* --- Swap the content of two vectors ---*/
  template<typename T>
  void SU2Mat<T>::swap(SU2Mat& su2mat) {
    std::swap(n_rows, su2mat.n_rows);
    std::swap(n_cols, su2mat.n_cols);
    std::swap(m_data, su2mat.m_data);
  }

  //
  //
  /*--- SU2Mat value constructor---*/
  template<typename T>
  SU2Mat<T>::SU2Mat(size_type _n_rows, size_type _n_cols, T init): n_rows(_n_rows), n_cols(_n_cols) {
    allocate();
    initialize(init);
  }

  //
  //
  /*--- SU2Mat copy constructor---*/
  template<typename T>
  SU2Mat<T>::SU2Mat(const SU2Mat& su2mat): n_rows(su2mat.n_rows), n_cols(su2mat.n_cols) {
    allocate();
    std::copy(su2mat.m_data, su2mat.m_data + size() ,m_data);
  }

  //
  //
  /*--- SU2Mat move constructor---*/
  template<typename T>
  SU2Mat<T>::SU2Mat(SU2Mat&& su2mat): n_rows(0), n_cols(0), m_data(NULL) {
    su2mat.swap(*this);
  }

  //
  //
  /*--- SU2Mat copy assignment operator ---*/
  template<typename T>
  SU2Mat<T>& SU2Mat<T>::operator=(const SU2Mat& su2mat) {
    if(this!= &su2mat) {
      n_rows = su2mat.n_rows;
      n_cols = su2mat.n_cols;
      if(m_data != NULL)
        delete[] m_data;
      m_data = new value_type(size());
      std::copy(su2mat.m_data, su2mat.m_data + size() ,m_data);
    }
    return *this;
  }

  //
  //
  /*--- SU2Mat move assignment operator ---*/
  template<typename T>
  SU2Mat<T>& SU2Mat<T>::operator=(SU2Mat&& su2mat) {
    if(this!= &su2mat) {
      n_rows = 0;
      n_cols = 0;
      if(m_data != NULL) {
        delete[] m_data;
        m_data = NULL;
      }
      su2mat.swap(*this);
    }
    return *this;
  }

  //
  //
  /*--- SU2Mat copy constructor from MatExprT ---*/
  template<typename T>
  template<class Derived>
  SU2Mat<T>::SU2Mat(const MatExprT<Derived,T>& e): n_rows(e.nbRows()), n_cols(e.nbCols()) {
    //const Derived& et(e);
    allocate();
    for(size_type i = 0; i < n_rows; ++i) {
      for(size_type j = 0; j < n_cols; ++j)
        (*this)(i,j) = e(i,j);
    }
  }

  //
  //
  /*--- SU2Mat assignment operator from MatExprT ---*/
  template<typename T>
  template<class Derived>
  SU2Mat<T>& SU2Mat<T>::operator=(const MatExprT<Derived,T>& e) {
    //const Derived& et(e);
    n_rows = e.nbRows();
    n_cols = e.nbCols();
    if(m_data != NULL)
      delete[] m_data;
    allocate();
    for(size_type i = 0; i < n_rows; ++i) {
      for(size_type j = 0; j < n_cols; ++j)
        (*this)(i,j) = e(i,j);
    }
    return *this;
  }

  //
  //
  /*--- SU2Mat copy constructor from already allocated memory ---*/
  template<typename T>
  SU2Mat<T>::SU2Mat(T** data, size_type _n_rows, size_type _n_cols): n_rows(_n_cols), n_cols(_n_cols) {
    allocate();
    for(size_type i = 0; i < _n_rows; ++i)
      for(size_type j = 0; j < _n_cols; ++j)
        (*this)(i,j) = data[i][j];
  }

  //
  //
  /*--- SU2Mat copy constructor from standard library vector of vectors ---*/
  template<typename T>
  SU2Mat<T>::SU2Mat(const std::vector<std::vector<T>>& mat): n_rows(mat.size()) {
    SU2_Assert(std::all_of(mat.cbegin(), mat.cend(),[&](const std::vector<T>& x){return x.size() == mat.front().size();}),
               "The number of columns is not homogeneous");
    n_cols = mat.front().size();
    allocate();
    for(size_type i = 0; i < n_rows; ++i) {
      for(size_type j = 0; j < n_cols; ++j)
        (*this)(i,j) = mat[i][j];
    }
  }

  //
  //
  /*--- Resize container ---*/
  template<typename T>
  void SU2Mat<T>::resize(size_type _n_rows, size_type _n_cols, value_type init) {
    if(nbRows() == _n_rows && nbCols() == _n_cols)
      return;

    size_type new_size = _n_rows * _n_cols;
    auto new_data = new value_type[new_size];

    if(_n_rows > nbRows() || _n_cols > nbCols())
      std::fill(new_data, new_data + new_size, init);

    size_type rowsToCopy = _n_rows < nbRows() ? _n_rows : nbRows();
    size_type colsToCopy = _n_cols < nbCols() ? _n_cols : nbCols();

    for(size_type r=0; r < rowsToCopy; ++r) {
      auto srcBegin = m_data + (r*nbCols());
      auto srcEnd = srcBegin + colsToCopy;
      auto dstBegin = new_data + (r*_n_cols);
      std::copy(srcBegin, srcEnd, dstBegin);
    }

    delete[] m_data;
    m_data = new_data;
    n_rows = _n_rows;
    n_cols = _n_cols;
  }

  //
  //
  /*--- Puts the copy of the desired row of the matrix in a supplied vector ---*/
  template<typename T>
  template<class Vec>
  Vec SU2Mat<T>::GetRow(size_type idx_row) const {
    SU2_Assert(idx_row < n_rows, "The index of the desired row exceeds the number of rows in the matrix");
    Vec vec(n_cols);
    for(size_type i = 0; i < n_cols; ++i)
      vec[i] = m_data[sub2ind(idx_row,i)];
    return vec;
  }

  //
  //
  /*--- Puts the copy of the desired column of the matrix in a supplied vector ---*/
  template<typename T>
  template<class Vec>
  Vec SU2Mat<T>::GetColumn(size_type idx_col) const {
    SU2_Assert(idx_col < n_cols, "The index of the desired column exceeds the number of columns in the matrix");
    Vec vec(n_rows);
    for(size_type i = 0; i < n_rows; ++i)
      vec[i] = m_data[sub2ind(i,idx_col)];
    return vec;
  }

  //
  //
  /*---  Sets the desired row of the matrix from a vector ---*/
  template<typename T>
  template<class Vec>
  void SU2Mat<T>::SetRow(size_type idx_row, const Vec& row)  {
    SU2_Assert(idx_row < n_rows,"The index of row to set exceeds the number of rows in the matrix");
    SU2_Assert(row.size() == n_cols,"The vector of data to set the row is different from the number of columns in the matrix");
    for(size_type i = 0; i < n_cols; ++i) {
      m_data[sub2ind(idx_row,i)] = row[i];
    }
  }

  //
  //
  /*---  Sets the desired column of the matrix from a vector ---*/
  template<typename T>
  template<class Vec>
  void SU2Mat<T>::SetColumn(size_type idx_col, const Vec& col)  {
    SU2_Assert(idx_col < n_cols,"The index of col to set exceeds the number of columns in the matrix");
    SU2_Assert(col.size() == n_rows,"The vector of data to set the column is different from the number of rows in the matrix");
    for(size_type i = 0; i < n_rows; ++i) {
      m_data[sub2ind(i,idx_col)] = col[i];
    }
  }

  //
  //
  /*--- Determinant of the matrix ---*/
  template<typename T>
  T SU2Mat<T>::determinant(void) const {
    static_assert(std::is_arithmetic<T>::value,"The type of the matrix is not arithmetic. You can't compute determinant");

    size_type dim = nbRows();
    SU2_Assert(dim == nbCols(),"The matrix is not squared");
    SU2_Assert(dim > 0,"The matrix is empty");
    if(dim == 1)
      return (*this)(0,0);
    else if(dim == 2)
      return (*this)(0,0)*(*this)(1,1) - (*this)(0,1)*(*this)(1,0);
    else {
      auto LU = LU_decomp();
      value_type det = LU.second(0,0);
      for(size_type i = 1; i < dim; ++i)
        det *= LU.second(i,i);
      return det;
    }
  }

  //
  //
  /*--- Inverse of the matrix ---*/
  template<typename T>
  SU2Mat<T> SU2Mat<T>::inverse(void) const {
    static_assert(std::is_arithmetic<T>::value,"The type of the matrix is not arithmetic. You can't compute inverse");

    size_type dim = nbRows();
    SU2_Assert(dim == nbCols(),"The matrix is not squared");
    SU2_Assert(dim > 0,"The matrix is empty");

    auto LU = LU_decomp();
    value_type det = LU.second(0,0);
    for(size_type i = 1; i < dim; ++i)
      det *= LU.second(i,i);
    if(std::abs(det) < EPS)
      throw std::runtime_error("Matrix singular. Cannot invert");

    if(check_diag())
      return invert_diag();

    SU2Mat result(dim,dim);

    for(size_type i = 0; i < dim; ++i) {
      std::vector<value_type> b(dim);
      b[i] = 1.0;
      result.SetColumn(i,MathTools::Backward_Sub(LU.second,MathTools::Forward_Sub(LU.first,b)));
    }
    return result;
  }


  //
  //
  /*--- LU decomposition of the matrix ---*/
  template<typename T>
  std::pair<SU2Mat<T>,SU2Mat<T>> SU2Mat<T>::LU_decomp(void) const {
    static_assert(std::is_arithmetic<T>::value,"The type of the matrix is not arithmetic. You can't compute LU decomposition");

    size_type dim = nbRows();
    SU2_Assert(dim == nbCols(),"The matrix is not squared");
    SU2_Assert(dim > 0,"The matrix is empty");

    SU2Mat L(dim,dim),U(dim,dim);

    for(size_type i = 0; i < dim; ++i) {
      // Upper Triangular
      for(size_type k = i; k < dim; ++k) {
        // Summation of L(i, j) * U(j, k)
        value_type sum = value_type();
        for(size_type j = 0; j < i; ++j)
          sum += (L(i,j)*U(j,k));
          // Evaluating U(i, k)
          U(i,k) = (*this)(i,k) - sum;
      }

      // Lower Triangular
      L(i,i) = 1.0;
      for(size_type k = i + 1; k < dim; ++k) {
        // Summation of L(k, j) * U(j, i)
        value_type sum = value_type();
        for(size_type j = 0; j < i; ++j)
          sum += L(k,j) * U(j,i);

        // Evaluating L(k, i)
        L(k,i) = ((*this)(k,i) - sum) / U(i,i);
      }
    }
    return std::make_pair(L,U);
  }


  //
  //
  /*--- Adjoint of the matrix ---*/
  template<typename T>
  SU2Mat<T> SU2Mat<T>::adjoint(void) const {
    static_assert(std::is_arithmetic<T>::value,"The type of the matrix is not arithmetic. You can't compute adjoint");

    size_type dim = nbRows();
    SU2_Assert(dim == nbCols(),"The matrix is not squared");
    SU2_Assert(dim > 0,"The matrix is empty");

    SU2Mat result(dim,dim);
    if(dim == 1) {
      result(0,0) = 1;
      return result;
    }

    for(size_type i=0; i < dim; ++i) {
      for(size_type j=0; j < dim; ++j) {
        // Get cofactor of A(i,j)
        auto tmp = GetCoFactor(i,j);

        // Sign of adj(i,j) positive if sum of row and column indexes is even.
        int sign = ((i+j) % 2 == 0) ? 1: -1;

        // Interchanging rows and columns to get the transpose of the cofactor matrix
        result(j,i) = sign*tmp.determinant();
      }
    }
    return result;
  }

  //
  //
  /*--- Cofactor of element (i,j) of the matrix ---*/
  template<typename T>
  SU2Mat<T> SU2Mat<T>::GetCoFactor(size_type idx_row,size_type idx_col) const {
    SU2_Assert(idx_row < nbRows(),"The index of desired row for cofactor exceeds the number of rows in the matrix");
    SU2_Assert(idx_col < nbCols(),"The index of desired column for cofactor exceeds the number of columns in the matrix");

    SU2Mat result(nbRows() - 1 ,nbCols() - 1);
    size_type i = 0, j = 0;

    // Looping for each element of the matrix
    for(size_type curr_row = 0; curr_row < nbRows(); ++curr_row) {
      for(size_type curr_col = 0; curr_col < nbCols(); ++curr_col) {
        //  Copying into result matrix only those element
        //  which are not in given row and column
        if(curr_row != idx_row && curr_col != idx_col)  {
          result(i,j++) = (*this)(curr_row,curr_col);

          // Row is filled, so increase row index and
          // reset col index
          if(curr_col == nbCols() - 1) {
            j = 0;
            i++;
          }
        }
      }
    }
    return result;
  }


  //
  //
  /*--- Transpose of a matrix ---*/
  template<typename T>
  SU2Mat<T> SU2Mat<T>::transpose(void) const {
    SU2Mat result(nbCols(),nbRows());
    for(size_type i = 0; i < n_cols; ++i) {
      for(size_type j = 0; j < n_rows; ++j)
        result(i,j) = m_data[this->sub2ind(j,i)];
    }
    return result;
  }

  //
  //
  /*--- Check if a matrix is diagonal ---*/
  template<typename T>
  bool SU2Mat<T>::check_diag(void) const {
    if(nbRows() != nbCols())
      return false;
    else {
      for(size_type i = 0; i < n_rows; ++i) {
        for(size_type j = 0; j < n_cols; ++j) {
          if((*this)(i,j) != value_type() && i!=j)
            return false;
        }
      }
      return true;
    }
  }

  //
  //
  /*--- Invert a diagonal matrix ---*/
  template<typename T>
  SU2Mat<T> SU2Mat<T>::invert_diag(void) const {
    static_assert(std::is_arithmetic<T>::value,"The type of the matrix is not arithmetic. You can't compute inverse");

    SU2_Assert(check_diag(),"You can invert only diagonal matrixes");
    SU2Mat result(n_rows,n_rows);
    value_type temp = value_type();
    for(size_type i = 0; i < n_rows; ++i) {
      temp = (*this)(i,i);
      SU2_Assert(std::abs(temp) > EPS,std::string("The diagional entry " + std::to_string(i) + " is equal to 0"));
      result(i,i) = 1.0/temp;
    }
    return result;
  }

  //
  //
  /*--- Operator / between a constant and a matrix ---*/
  template<class Derived,typename T>
  SU2Mat<T> operator/(double l, const MatExprT<Derived,T>& r) {
    static_assert(std::is_arithmetic<T>::value,"The type of the matrix is not arithmetic");

    SU2Mat<T> result(r);
    SU2_Assert(result.check_notzero(),"The matrix on rhs contains zero and so you can't compute its reciprocal");
    for(std::size_t i = 0; i < result.nbRows(); ++i) {
      for(std::size_t j = 0; j < result.nbCols(); ++j)
        result(i,j) = l/result(i,j);
    }
    return result;
  }

  //
  //
  /*--- Output stream operator ---*/
  template<class Derived, typename T>
  std::ostream& operator<<(std::ostream& out, const MatExprT<Derived,T>& A) {
    out<<std::fixed<<std::setprecision(4);
    for(std::size_t i = 0; i < A.nbRows(); ++i) {
      for(std::size_t j = 0; j < A.nbCols(); ++j)
        out<<std::right<<std::setw(10)<<A(i,j);
      out<<'\n';
    }
    return out;
  }


  typedef SU2Mat<double> RealMatrix;

} /*--- End of namespace Common ---*/

#endif
