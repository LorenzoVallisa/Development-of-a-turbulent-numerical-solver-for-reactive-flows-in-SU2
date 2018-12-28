#ifndef SU2_VECTORT_HPP
#define SU2_VECTORT_HPP

#include "exprT.hpp"
#include "../option_structure.hpp"

#include <vector>
#include <iomanip>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <numeric>

namespace Common {

  /*!
   * Declaration of output stream operator
   */
  template<class Derived,typename T>
  std::ostream& operator<<(std::ostream& out, const ExprT<Derived,T>& v);

  /*!
   * Definition of a class Su2Vec that implements an expression template technique
   * \author G.Orlando
  */

  template <typename T>
  class SU2Vec : public ExprT<SU2Vec<T>,T> {
  public:
    using value_type = T;
    using size_type = std::size_t;

    /*!
      * \brief Default constructor
    */
    SU2Vec(): N(), m_data(NULL) {}

    /*!
      * \brief Class constructor
      * \param[in] _N - size of the vector
      * \param[in] init - value to initialize (default provided)
    */
    SU2Vec(size_type _N,value_type init = value_type());

    /*!
      * \brief Class destructor
    */
    ~SU2Vec() {
      free_mem();
    }

    /*!
      * \brief Copy constructor
      * \param[in] su2vec - SU2Vec to copy from
    */
    SU2Vec(const SU2Vec& su2vec);

    /*!
      * \brief Move constructor
      * \param[in] su2vec - SU2Vec to copy from
    */
    SU2Vec(SU2Vec&& su2vec);

    /*!
      * \brief Copy assignment operator
      * \param[in] su2vec - SU2Vec to copy from
    */
    SU2Vec& operator=(const SU2Vec& su2vec);

    /*!
      * \brief Move assignment operator
      * \param[in] su2vec - SU2Vec to copy from
    */
    SU2Vec& operator=(SU2Vec&& su2vec);

    /*!
      * \brief Copy Constructor with range
      * \param[in] init - iterator to first datum
      * \param[in] end - iterator to one past last datum
    */
    template<typename Iterator>
    SU2Vec(Iterator b,Iterator e);

    /*!
      * \brief Copy Constructor from standard library vector
      * \param[in] v - std::vector to copy from
    */
    explicit SU2Vec(const std::vector<value_type>& v);

    /*!
      * \brief Copy constructor from ExprT;
      * \param[in] e - ExprT to copy from
     */
    template <class Derived>
    SU2Vec(const ExprT<Derived,value_type>& e);

    /*!
      * \brief Assignment operator from ExprT;
      * \param[in] e - ExprT to assign
     */
    template <class Derived>
    SU2Vec& operator=(const ExprT<Derived,value_type>& e);

    /*!
      * \brief Overloaded assignment operator with constants (macro used);
     */
    #define SU2VEC_ASSIGN_OP_CONST(__op__) \
    SU2Vec& operator __op__ (value_type value) { \
      static_assert(std::is_arithmetic<value_type>::value,"The type of the vector is not arithmetic"); \
                                                \
      SU2_Assert((strcmp(#__op__ ,"/=") == 0 && std::abs(value) > EPS) || (strcmp(#__op__,"/=") != 0),"You can't divide by zero"); \
      for(size_type i = 0; i < size(); ++i) \
        m_data[i] __op__ value; \
      return *this; \
    }
    SU2VEC_ASSIGN_OP_CONST(+=)
    SU2VEC_ASSIGN_OP_CONST(-=)
    SU2VEC_ASSIGN_OP_CONST(*=)
    SU2VEC_ASSIGN_OP_CONST(/=)
    #undef SU2VEC_ASSIGN_OP_CONST

    /*!
     * \brief Overloaded assignment operator with ExprT (macro used);
    */
    #define SU2VEC_ASSIGN_OP_EXPR(__op__) \
    template<class Derived> \
    SU2Vec& operator __op__ (const ExprT<Derived,value_type>& e) { \
      static_assert(std::is_arithmetic<value_type>::value,"The type of the vector is not arithmetic"); \
                                                \
      SU2_Assert(this->size() == e.size(),"The size of the vectors is not the same and so you can't operate with them"); \
      for(size_type i = 0; i < size(); ++i) \
        m_data[i] __op__ e[i]; \
      return *this; \
    }
    SU2VEC_ASSIGN_OP_EXPR(+=)
    SU2VEC_ASSIGN_OP_EXPR(-=)
    SU2VEC_ASSIGN_OP_EXPR(*=)
    #undef SU2VEC_ASSIGN_OP_EXPR

    /*!
     * \brief Overloaded /= operator with ExprT (checks if rhs contains zeros);
    */
    template<class Derived>
    SU2Vec& operator /= (const ExprT<Derived,value_type>& e) {
      static_assert(std::is_arithmetic<value_type>::value,"The type of the vector is not arithmetic");

      SU2_Assert(this->size() == e.size(),"The size of the vectors is not the same and so you can't operate with them");
      SU2Vec tmp(e);
      SU2_Assert(tmp.check_notzero(),"The vector on rhs contains zeros and so you can't divide by it");
      for(size_type i = 0; i < size(); ++i) \
        m_data[i] /= tmp[i]; \
      return *this;
    }

    /*!
     * \brief Returns i-th element (non const version)
     * \param[in] i - index of the element
    */
    inline value_type& operator[](size_type i) {
      return m_data[i];
    }

    /*!
     * \brief Returns i-th element (const version)
     * \param[in] i - index of the element
    */
    inline const value_type& operator[](size_type i) const {
      return m_data[i];
    }

    /*!
     * \brief Returns i-th element (non const version)
     * \param[in] i - index of the element
    */
    inline value_type& at(size_type i) {
      SU2_Assert(i < size(),"Index is beyond the size of the vector");
      return m_data[i];
    }

    /*!
     * \brief Returns i-th element (const version)
     * \param[in] i - index of the element
    */
    inline const value_type& at(size_type i) const {
      SU2_Assert(i < size(),"Index is beyond the size of the vector");
      return m_data[i];
    }

    /*!
     * \brief Returns the size of the vector
    */
    inline size_type size(void) const {
      return N;
    }

    /*!
     * \brief Clears the vector
    */
    inline void clear(void) {
      free_mem();
      N = 0;
    }

    /*!
     * \brief Swap the content of two vectors
     * \param[in] su2vec - Vector to swap with the current one
    */
    void swap(SU2Vec& su2vec);

    /*!
     * \brief Resize of the vector
     * \param[in] _N - new size of the vector
     * \param[in] init - value to initialize (default provided)
    */
    void resize(size_type _N, value_type init = value_type());

    /*!
     * \brief Add an element at the end of the vector
     * \param[in] value - value to add
    */
    void push_back(const value_type& value);

    /*!
     * \brief Returns the underlined pointer
    */
    inline value_type* data(void) const {
      return m_data;
    }

    inline value_type norm2(void) const {
      static_assert(std::is_arithmetic<T>::value,"The type of the vector is not arithmetic");
      return std::sqrt(std::inner_product(this->cbegin(),this->cend(),this->cbegin(),0.0));
    }

    /*!
     * \brief Cast to standard library vector
    */
    inline operator std::vector<value_type> () const {
      return std::vector<value_type>(m_data,m_data + N);
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
      return static_cast<decltype(std::declval<std::vector<value_type>>().end())>(m_data + N);
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
      return static_cast<decltype(std::declval<std::vector<value_type>>().cend())>(m_data + N);
    }

    /*!
     * \brief Check if there are zeros in the vector
    */
    inline bool check_notzero(void) const {
      return std::none_of(this->cbegin(), this->cend(), [&](const value_type& x){return x == value_type();});
    }

    /*!
     * \brief Overloading of output stream operator
    */
    friend std::ostream& operator<< TF (std::ostream& out, const ExprT<SU2Vec,value_type>& v);

  private:

    size_type N;  /*!< \brief Size of the vector. */
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

  }; /*--- End of class declaration SU2Vec ---*/

  //
  //
  /*--- Allocate memory---*/
  template<typename T>
  void SU2Vec<T>::allocate(void) {
    if(N > 0)
      m_data = new value_type[N];
  }

  //
  //
  /*--- Initialize vector ---*/
  template<typename T>
  void SU2Vec<T>::initialize(value_type value) {
    std::fill(m_data, m_data + N, value);
  }

  //
  //
  /*--- Free memory ---*/
  template<typename T>
  void SU2Vec<T>::free_mem(void) {
    if(m_data!= NULL) {
      delete[] m_data;
      m_data = NULL;
    }
  }

  //
  //
  /* --- Swap the content of two vectors ---*/
  template<typename T>
  void SU2Vec<T>::swap(SU2Vec& su2vec) {
    std::swap(N, su2vec.N);
    std::swap(m_data, su2vec.m_data);
  }

  //
  //
  /*--- SU2Vec value constructor---*/
  template<typename T>
  SU2Vec<T>::SU2Vec(size_type _N, value_type init): N(_N) {
    allocate();
    initialize(init);
  }

  //
  //
  /*--- SU2Vec copy constructor---*/
  template<typename T>
  SU2Vec<T>::SU2Vec(const SU2Vec& su2vec):N(su2vec.N) {
    allocate();
    std::copy(su2vec.m_data, su2vec.m_data + N, m_data);
  }

  //
  //
  /*--- SU2Vec move constructor---*/
  template<typename T>
  SU2Vec<T>::SU2Vec(SU2Vec&& su2vec): N(0), m_data(NULL) {
    su2vec.swap(*this);
  }

  //
  //
  /*--- SU2Vec copy assignment operator ---*/
  template<typename T>
  SU2Vec<T>& SU2Vec<T>::operator=(const SU2Vec& su2vec) {
    if(this!= &su2vec) {
      N = su2vec.N;
      if(m_data!= NULL)
        delete[] m_data;
      allocate();
      std::copy(su2vec.m_data, su2vec.m_data + N, m_data);
    }
    return *this;
  }

  //
  //
  /*--- SU2Vec range constructor---*/
  template<typename T>
  template<typename Iterator>
  SU2Vec<T>::SU2Vec(Iterator b,Iterator e): N(std::distance(b,e)) {
    allocate();
    for(size_type i = 0; i < N; ++i)
      m_data[i] = (*b++);
  }

  //
  //
  /*--- SU2Vec move assignment operator ---*/
  template<typename T>
  SU2Vec<T>& SU2Vec<T>::operator=(SU2Vec&& su2vec) {
    if(this!= &su2vec) {
      N = 0;
      if(m_data!= NULL)
        delete[] m_data;
      m_data = NULL;
      su2vec.swap(*this);
    }
    return *this;
  }


  //
  //
  /*--- SU2Vec copy constructor from ExprT ---*/
  template<typename T>
  template<class Derived>
  SU2Vec<T>::SU2Vec(const ExprT<Derived,T>& e) : N(e.size()) {
    //const Derived& et(e);
    allocate();
    for(size_type i = 0; i < N; ++i)
      m_data[i] = e[i];
  }

  //
  //
  /*--- SU2Vec assignment operator from ExprT ---*/
  template<typename T>
  template<class Derived>
  SU2Vec<T>& SU2Vec<T>::operator=(const ExprT<Derived,T>& e) {
    //const Derived& et(e);
    N = e.size();
    if(m_data != NULL)
      delete[] m_data;
    allocate();
    for(size_type i = 0; i < N; ++i)
      m_data[i] = e[i];
    return *this;
  }

  //
  //
  /*--- SU2Vec copy constructor from standard library vector ---*/
  template<typename T>
  SU2Vec<T>::SU2Vec(const std::vector<value_type>& v):N(v.size()) {
    allocate();
    std::copy(v.cbegin(), v.cend(), m_data);
  }

  //
  //
  /*--- Resize container ---*/
  template<typename T>
  void SU2Vec<T>::resize(size_type _N, value_type init) {
    if(_N == size())
      return;
    if(_N < size()) {
      N = _N;
      auto tmp = new T[size()];
      std::copy(m_data, m_data + N ,tmp);
      delete[] m_data;
      m_data = tmp;
    }
    else {
      auto tmp = new value_type[_N];
      size_type i;
      std::fill(tmp + size(),tmp + _N,init);
      std::copy(m_data, m_data + N ,tmp);
      N = _N;
      delete[] m_data;
      m_data = tmp;
    }
  }

  //
  //
  /*--- Add an element at the end of the container ---*/
  template<typename T>
  void SU2Vec<T>::push_back(const value_type& value) {
    value_type* tmp = new value_type[++N];
    std::copy(m_data, m_data + N - 1,tmp);
    tmp[N - 1] = value;
    if(m_data != NULL)
      delete[] m_data;
    m_data = tmp;
  }

  //
  //
  /*--- Operator / between two vectors ---*/
  template<class Derived_Left,class Derived_Right,typename T>
  SU2Vec<T> operator/(const ExprT<Derived_Left,T>& l, const ExprT<Derived_Right,T>& r) {
    static_assert(std::is_arithmetic<T>::value,"The type of the vector is not arithmetic"); \

    SU2_Assert(l.size() == r.size(),"The size of the vectors is not the same and so you can't operate with them");
    SU2Vec<T> result(r);
    SU2_Assert(result.check_notzero(),"The vector on rhs contains zero and so you can't divide by it");
    for(std::size_t i = 0; i < result.size(); ++i)
      result[i] = l[i]/result[i];
    return result;
  }

  //
  //
  /*--- Operator / between a constant and a vector ---*/
  template<class Derived,typename T>
  SU2Vec<T> operator/(double l, const ExprT<Derived,T>& r) {
    static_assert(std::is_arithmetic<T>::value,"The type of the vector is not arithmetic");

    SU2Vec<T> result(r);
    SU2_Assert(result.check_notzero(),"The vector on rhs contains zero and so you can't compute its reciprocal");
    for(std::size_t i = 0; i < result.size(); ++i)
      result[i] = l/result[i];
    return result;
  }

  //
  //
  /*--- Stream output operator ---*/
  template<class Derived,typename T>
  std::ostream& operator<<(std::ostream& out,const ExprT<Derived,T>& v) {
    out<<std::fixed<<std::setprecision(4);
    for(std::size_t i = 0; i < v.size(); ++i)
      out<<std::right<<std::setw(10)<<v[i]<<std::endl;
    out<<'\n';
    return out;
  }

  typedef SU2Vec<double> RealVec;

} /*--- End of namespace Common ---*/

#endif
