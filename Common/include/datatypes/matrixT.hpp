#ifndef SU2_MATRIXT_HPP
#define SU2_MATRIXT_HPP

#include "mat_exprT.hpp"
#include <vector>
#include <algorithm>
#include <iomanip>

namespace Common {

  /*!
   * Forward declaration of SU2Mat class
  */
  template<typename T> class SU2Mat;

  /*!
   * Declaration of output stream operator
  */
  template<typename T> std::ostream& operator<<(std::ostream& out,const SU2Mat<T>& A);

  /*!
   * Definition of a class SU2Mat that implements an expression template technique
   * \author G.Orlando
  */

  template <typename T>
  class SU2Mat : public MatExprT<SU2Mat<T>,T> {
  public:

    using Type = T;
    using value_type = T;
    typedef typename std::vector<std::vector<T>>::size_type size_type;

    /*!
      * \brief Default constructor
    */
    SU2Mat(): n_rows(),n_cols(),m_data(NULL) {}

    /*!
      * \brief Class constructor
      * \param[in] _n_rows - number of rows in the matrix
      * \param[in] _n_cols - number of columns in the matrix
      * \param[in] init - value to initialize (default provided)
    */
    SU2Mat(std::size_t _n_rows,std::size_t _n_cols,Type init = Type());

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
      * \param[in] data - pointer to pointer to Type (stored elements)
      * \param[in] _n_rows - number of rows of the resulting matrix
      * \param[in] _n_cols - number of columns of the resulting matrix
    */
    SU2Mat(Type** data,std::size_t _n_rows,std::size_t _n_cols):n_rows(_n_cols),n_cols(_n_cols) {
      allocate();
      for (std::size_t i = 0; i < _n_rows; ++i)
        for (std::size_t j = 0; j < _n_cols; ++j)
          (*this)(i,j) = data[i][j];
    }


    /*!
      * \brief Copy Constructor from standard library vector (matrix as vector of vectors)
      * \param[in] v - std::vector to copy from
    */
    explicit SU2Mat(const std::vector<std::vector<Type>>& mat);

    /*!
      * \brief Coyp constructor from MatExprT;
      * \param[in] e - MatExprT to copy from
     */
    template <class Derived>
    SU2Mat(const MatExprT<Derived,Type>& e):n_rows(e.nbRows()),n_cols(e.nbCols()) {
      //const Derived& et(e);
      allocate();
      for(std::size_t i = 0; i < n_rows; ++i) {
        for(std::size_t j = 0; j < n_cols; ++j)
          (*this)(i,j) = e(i,j);
      }
    }

    /*!
      * \brief Assignment operator from MatExprT;
      * \param[in] e - MatExprT to copy from
     */
    template <class Derived>
    SU2Mat& operator=(MatExprT<Derived,Type>& e) {
      //const Derived& et(e);
      n_rows = e.nbRows();
      n_cols = e.nbCols();
      if(m_data!=NULL)
        delete[] m_data;
      allocate();
      for(std::size_t i = 0; i < n_rows; ++i) {
        for(std::size_t j = 0; j < n_cols; ++j)
          (*this)(i,j) = e(i,j);
      }
      return *this;
    }

    /*!
      * \brief Overloaded assignment operator with constants (macro used);
    */

    #define SU2MAT_ASSIGN_OP_CONST(__op__) \
    SU2Mat& operator __op__ (Type value) { \
      SU2_Assert(strcmp(#__op__ ,"/=") && std::abs(value)>0,"You can't divide by zero"); \
      for(std::size_t i = 0; i < size(); ++i) \
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
    SU2Mat& operator __op__ (const MatExprT<Derived,Type>& e) { \
      SU2_Assert(this->nbRows() == e.nbRows() && this->nbCols() == e.nbCols(), \
                "The dimension of the matrixes is not compatible and so you can't operate with them"); \
      for(std::size_t i = 0; i < nbRows(); ++i) { \
        for(std::size_t j = 0; j < nbCols(); ++j) \
          (*this)(i,j) __op__ e(i,j); \
      } \
      return *this; \
    }
    //SU2MAT_ASSIGN_OP_EXPR(=)
    SU2MAT_ASSIGN_OP_EXPR(+=)
    SU2MAT_ASSIGN_OP_EXPR(-=)
    //SU2MAT_ASSIGN_OP_EXPR(*=)
    //SU2MAT_ASSIGN_OP_CONST(/=)
    #undef SU2MAT_ASSIGN_OP_EXPR

    /*!
      * \brief Overloaded *= operator with MatExprT (matrix-matrix multiplication);
    */

    template<class Derived>
    SU2Mat& operator*=(const MatExprT<Derived,Type>& e) {
      SU2_Assert(this->nbCols() == e.nbRows(), "The dimension of the matrixes is not compatible and so you can't multiply them"); \
      std::vector<Type> tmp(m_data, m_data + size());
      n_cols = e.nbCols();
      this->resize(n_rows,n_cols);
      for (std::size_t i = 0; i < n_rows; ++i) {
        for (std::size_t j = 0; j < n_cols; ++j) {
          for (std::size_t k = 0; k < e.nbRows(); ++k)
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
    inline Type& operator()(std::size_t i,std::size_t j) {
      return m_data[sub2ind(i,j)];
    }

    /*!
     * \brief Returns the (i,j)-th element (const version)
     * \param[in] i - index of row
     * \param[in] j - index of columns
    */
    //inline Type operator()(std::size_t i,std::size_t j) const
    inline const Type& operator()(std::size_t i,std::size_t j) const {
      return m_data[sub2ind(i,j)];
    }

    /*!
     * \brief Returns the (i,j)-th element (non const version)
     * \param[in] i - index of row
     * \param[in] j - index of columns
    */
    inline Type& at(std::size_t i,std::size_t j) {
      SU2_Assert(i < n_rows,"Index of row is beyond the number of rows in the matrix");
      SU2_Assert(j < n_cols,"Index of col is beyond the number of columns in the matrix");
      return m_data[sub2ind(i,j)];
    }

    /*!
     * \brief Returns the (i,j)-th element (const version)
     * \param[in] i - index of row
     * \param[in] j - index of columns
    */
    //inline Type at(std::size_t i,std::size_t j) const
    inline const Type& at(std::size_t i,std::size_t j) const {
      SU2_Assert(i < n_rows,"Index of row is beyond the number of rows in the matrix");
      SU2_Assert(j < n_cols,"Index of col is beyond the number of columns in the matrix");
      return m_data[sub2ind(i,j)];
    }

    /*!
     * \brief Returns the number of rows in the matrix
    */
    inline std::size_t nbRows(void) const {
      return n_rows;
    }

    /*!
     * \brief Returns the number of rows in the matrix
    */
    inline std::size_t nbCols(void) const {
      return n_cols;
    }

    /*!
     * \brief Returns the size of the matrix
    */
    inline std::size_t size(void) const {
      return n_rows*n_cols;
    }

    /*!
     * \brief Resize the matrix
     * \param[in] _n_rows - new number of rows in the matrix
     * \param[in] _n_cols - new number of columns in the matrix
     * \param[in] init - value to initialize (default provided)
    */
    void resize(std::size_t _n_rows, std::size_t _n_cols,Type init = Type());

    /*!
     * \brief Returns the underlined pointer
    */
    inline Type* data(void) const {
      return m_data;
    }

    /*!
     * \brief Compute the determinant of a 2 x 2 matrix
    */
    Type determ2(void) const;

    /*!
     * \brief Compute the determinant of a 3 x 3 matrix
    */
    Type determ3(void) const;

    /*!
     * \brief Returns the transpose of the matrix (the original matrix remains untouched).
     * \param[out] result - the transposed matrix
    */
    void transpose(SU2Mat<T>& result) const;

    /*!
     * \brief Check if the matrix is diagonal
    */
    bool check_diag(void) const;

    /*!
     * \brief Returns the invers of a diagonal matrix (the original matrix remains untouched).
     * \param[out] result - the transposed matrix
    */
    void invert_diag(SU2Mat<T>& result) const;

    /*!
     * \brief Puts the copy of the desired row of the matrix in a supplied vector
     * \param[in] idx_row - index of desired row
     * \param[out] vec - the resulting vector
    */
    template <class Vec>
    void GetRow(std::size_t idx_row, Vec& vec) const {
      SU2_Assert(idx_row < n_rows, "The index of the desired row exceeds the number of rows in the matrix");
      vec.resize(n_cols);
      for(std::size_t i = 0; i < n_cols; ++i)
        vec[i] = m_data[sub2ind(idx_row,i)];
    }

    /*!
     * \brief Puts the copy of the desired column of the matrix in a supplied vector
     * \param[in] idx_col - index of desired row
     * \param[out] vec - the resulting vector
    */
    template <class Vec>
    void GetColumn(std::size_t idx_col, Vec& vec) const  {
      SU2_Assert(idx_col < n_cols, "The index of the desired column exceeds the number of columns in the matrix");
      vec.resize(n_rows);
      for(std::size_t i = 0; i < n_cols; ++i)
        vec[i] = m_data[sub2ind(i,idx_col)];
    }

    /*!
     * \brief Sets the desired row of the matrix from a vector
     * \param[in] idx_row - index of desired row
     * \param[in] row - the vector where data comes from
    */
    template <class Vec>
    void SetRow(std::size_t idx_row, const Vec& row)  {
      SU2_Assert(idx_row < n_rows,"The index of row to set exceeds the number of rows in the matrix");
      SU2_Assert(row.size() == n_cols,"The vector of data to set the row is different from the number of columns in the matrix");
      for(std::size_t i = 0; i < n_cols; ++i) {
        m_data[sub2ind(idx_row,i)] = row[i];
      }
    }

    /*!
     * \brief Sets the desired column of the matrix from a vector
     * \param[in] idx_col - index of desired column
     * \param[in] col - the vector where data comes from
    */
    template <class Vec>
    void SetColumn(std::size_t idx_col, const Vec& col)  {
      SU2_Assert(idx_col < n_cols,"The index of col to set exceeds the number of columns in the matrix");
      SU2_Assert(col.size() == n_rows,"The vector of data to set the column is different from the number of rows in the matrix");
      for (std::size_t i = 0; i < n_rows; ++i) {
        m_data[sub2ind(i,idx_col)] = col[i];
      }
    }

    /*!
     * \brief Cast to standard library vector of vectors ("standard matrix")
    */
    operator std::vector<std::vector<Type>> () const {
      auto result = std::vector<Type>(n_rows,std::vector<Type>(n_cols));
      for(std::size_t i = 0; i < n_rows; ++i) {
        for(std::size_t j = 0; j < n_cols; ++j)
          result[i][j] = m_data[sub2ind(i,j)];
      }
      return result;
    }

    /*!
     * \brief Check if there are zeros in the matrix
    */
    inline bool check_notzero(void) {
      return std::none_of(this->cbegin(), this->cend(), [&](const Type& x){return x == Type();});
    }

    /*!
     * \brief Begin iterator (non const version)
    */
    inline auto begin()->decltype(std::declval<std::vector<Type>>().begin()) {
      return static_cast<decltype(std::declval<std::vector<Type>>().begin())>(m_data);
    }

    /*!
     * \brief End iterator (non const version)
    */
    inline auto end()->decltype(std::declval<std::vector<Type>>().end()) {
      return static_cast<decltype(std::declval<std::vector<Type>>().end())>(m_data + n_rows*n_cols);
    }

    /*!
     * \brief Begin iterator (const version)
    */
    inline auto cbegin()->decltype(std::declval<std::vector<Type>>().cbegin()) const {
      return static_cast<decltype(std::declval<std::vector<Type>>().cbegin())>(m_data);
    }

    /*!
     * \brief End iterator (const version)
    */
    inline auto cend()->decltype(std::declval<std::vector<Type>>().cend()) const {
      return static_cast<decltype(std::declval<std::vector<Type>>().cend())>(m_data + n_rows*n_cols);
    }

    /*!
     * \brief Overloading of output stream operator
    */
    friend std::ostream& operator<< TF (std::ostream& out,const SU2Mat<Type>& A);

  private:

    std::size_t n_rows;  /*!< \brief Number of rows of the matrix. */
    std::size_t n_cols;  /*!< \brief Number of columns of the matrix. */
    Type* m_data;  /*!< \brief Stored elements. */

    /*!
     * \brief Helper function to allocate memory
    */
    void allocate(void);

    /*!
     * \brief Helper function to initialize vector
     * \param[in] value - Value initialization
    */
    void initialize(Type value);

    /*!
     * \brief Helper function to free memory
    */
    void free_mem(void);

    /*!
     * \brief sub2ind
    */
    inline std::size_t sub2ind(std::size_t i, std::size_t j) const {
      return i*n_cols + j;
    }

  }; /*--- End of class declaration SU2Mat ---*/

  //
  //
  /*--- SU2Mat value constructor---*/
  template<typename T>
  SU2Mat<T>::SU2Mat(std::size_t _n_rows,std::size_t _n_cols,T init):n_rows(_n_rows),n_cols(_n_cols) {
    allocate();
    initialize(init);
  }

  //
  //
  /*--- SU2Mat copy constructor---*/
  template<typename T>
  SU2Mat<T>::SU2Mat(const SU2Mat& su2mat):n_rows(su2mat.n_rows),n_cols(su2mat.n_cols) {
    allocate();
    for(std::size_t i = 0; i<size(); ++i)
      m_data[i] = su2mat.m_data[i];
  }

  //
  //
  /*--- SU2Mat move constructor---*/
  template<typename T>
  SU2Mat<T>::SU2Mat(SU2Mat&& su2mat):n_rows(su2mat.n_rows),n_cols(su2mat.n_cols) {
    m_data = su2mat.m_data;
    su2mat.m_data = NULL;
    su2mat.n_rows = 0;
    su2mat.n_cols = 0;
  }

  //
  //
  /*--- SU2Mat copy assignment operator ---*/
  template<typename T>
  SU2Mat<T>& SU2Mat<T>::operator=(const SU2Mat& su2mat) {
      if(this!= &su2mat) {
          n_rows = su2mat.n_rows;
          n_cols = su2mat.n_cols;
          if(m_data!=NULL)
            delete[] m_data;
          m_data = new T(size());
          for(std::size_t i = 0; i<size(); ++i)
            m_data[i] = su2mat.m_data[i];
      }
      return *this;
  }

  //
  //
  /*--- SU2Mat move assignment operator ---*/
  template<typename T>
  SU2Mat<T>& SU2Mat<T>::operator=(SU2Mat&& su2mat) {
      if(this!= &su2mat) {
          n_rows = std::move(su2mat.n_rows);
          n_cols = std::move(su2mat.n_cols);
          if(m_data!=NULL)
            delete[] m_data;
          m_data = su2mat.m_data;
          su2mat.m_data = NULL;
          su2mat.n_rows = 0;
          su2mat.n_cols = 0;
      }
      return *this;
  }

  //
  //
  /*--- SU2Mat copy constructor from standard library vector of vectors ---*/
  template<typename T>
  SU2Mat<T>::SU2Mat(const std::vector<std::vector<T>>& mat):n_rows(mat.size()) {
    SU2_Assert(std::all_of(mat.cbegin(), mat.cend(),[&](const std::vector<T>& x){return x.size() == mat.front().size();}),
               "The number of columns is not homogeneous");
    n_cols = mat.front().size();
    allocate();
    for(std::size_t i = 0; i < n_rows; ++i) {
      for(std::size_t j = 0; j < n_cols; ++j)
        m_data[sub2ind(i,j)] = mat[i][j];
    }
  }

  //
  //
  /*--- Allocate memory---*/
  template<typename T>
  void SU2Mat<T>::allocate(void) {
    if (size() > 0)
      m_data = new T[size()];
  }

  //
  //
  /*--- Initialize matrix with a value ---*/
  template<typename T>
  void SU2Mat<T>::initialize(T value) {
    for(std::size_t i = 0; i < size(); ++i)
      m_data[i] = value;
  }

  //
  //
  /*--- free memory ---*/
  template<typename T>
  void SU2Mat<T>::free_mem(void) {
    if (m_data!= NULL) {
      delete [] m_data;
      m_data = NULL;
      n_rows = 0;
      n_cols = 0;
    }
  }

  //
  //
  /*--- Resize container ---*/
  template<typename T>
  void SU2Mat<T>::resize(std::size_t _n_rows,std::size_t _n_cols,T init) {
    free_mem();
    n_rows = _n_rows;
    n_cols = _n_cols;
    allocate();
    initialize(init);
  }

  //
  //
  /*--- Determinant of 2 x 2 matrix ---*/
  template<typename T>
  T SU2Mat<T>::determ2(void) const {
    SU2_Assert(nbRows() == 2,"The number of rows in the matrix is not equal to 2");
    SU2_Assert(nbCols() == 2,"The number of columns in the matrix is not equal to 2");
    return m_data[0]*m_data[3] - m_data[1]*m_data[2];
  }

  //
  //
  /*--- Determinant of 3 x 3 matrix ---*/
  template<typename T>
  T SU2Mat<T>::determ3(void) const {
    SU2_Assert(nbRows() == 3,"The number of rows in the matrix is not equal to 3");
    SU2_Assert(nbCols() == 3,"The number of columns in the matrix is not equal to 3");
    return m_data[0]*(m_data[4]*m_data[8] - m_data[5]*m_data[7]) -
           m_data[1]*(m_data[3]*m_data[8] - m_data[5]*m_data[6]) +
           m_data[2]*(m_data[3]*m_data[7] - m_data[4]*m_data[6]);
  }

  //
  //
  /*--- Transpose of a matrix ---*/
  template<typename T>
  void SU2Mat<T>::transpose(SU2Mat& result) const {
    result.resize(this->nbCols(),this->nbRows());
    for(std::size_t i = 0; i < n_cols; ++i) {
      for(std::size_t j = 0; j < n_rows; ++j)
        result(i,j) = m_data[this->sub2ind(j,i)];
    }
  }

  //
  //
  /*--- Check if a matrix is diagonal ---*/
  template<typename T>
  bool SU2Mat<T>::check_diag(void) const {
    if(nbRows() != nbCols())
      return false;
    else {
      for(std::size_t i = 0; i < n_rows; ++i) {
        for(std::size_t j = 0; j < n_cols; ++j) {
          if((*this)(i,j) != T() && i!=j)
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
  void SU2Mat<T>::invert_diag(SU2Mat& result) const {
    SU2_Assert(check_diag(),"You can invert only diagonal matrixes");
    result.resize(n_rows,n_rows);
    T temp = T();
    for(std::size_t i = 0; i < n_rows; ++i) {
      temp = (*this)(i,i);
      SU2_Assert(std::abs(temp) > 0.0,std::string("The diagional entry " + std::to_string(i) + " is equal to 0"));
      result(i,i) = 1.0/temp;
    }
  }

  //
  //
  /*--- Output stream operator ---*/
  template <typename T>
  std::ostream& operator<<(std::ostream& out, const SU2Mat<T>& A) {
    std::cout<<std::fixed<<std::setprecision(4);
    for (std::size_t i = 0; i < A.nbRows(); ++i) {
      for (std::size_t j = 0; j < A.nbCols(); ++j)
        out<<std::right<<A(i,j)<<std::setw(10);
    out<<'\n';
    }
    return out;
  }

  //
  //
  /*--- Operator / between a constant and a matrix ---*/
  template<class Derived,typename T>
  SU2Mat<T> operator/(double l, const MatExprT<Derived,T>& r) {
    static_assert(std::is_convertible<T,double>::value,"The type of the matrix is not convertible to double"); \
    SU2Mat<T> result(r);
    SU2_Assert(result.check_notzero(),"The matrix on rhs contains zero and so you can't compute its reciprocal");
    for (std::size_t i = 0; i < result.nbRows(); ++i) {
      for (std::size_t j = 0; j < result.nbCols(); ++j)
        result(i,j) = l/result(i,j);
    }
    return result;
  }

  using RealMatrix = SU2Mat<double>;


}

#endif
