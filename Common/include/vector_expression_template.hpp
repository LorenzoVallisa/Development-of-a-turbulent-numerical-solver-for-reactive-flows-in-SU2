#ifndef SU2_VECTOR_ET
#define SU2_VECTOR_ET

//////////////////////////////////////////////////////////////////////////////

#include "exprT.hpp"

namespace Common {

  template  <typename T, int N> class su2vec;

  template <typename T, int N> std::ostream& operator<< (std::ostream& out, const su2vec<T,N>& v);
  template <typename T, int N> std::istream& operator>> (std::istream& in, su2vec<T,N>& v);
  template <typename T, int N> bool operator== (const su2vec<T,N>& v1, const su2vec<T,N>& v2);
  template <typename T, int N> bool operator== (const su2vec<T,N>& v1, T value);
  template <typename T, int N> bool operator!= (const su2vec<T,N>& v1, const su2vec<T,N>& v2);
  template <typename T, int N> bool operator!= (const su2vec<T,N>& v1, T value);

  template <typename T, int N> bool operator<  (const su2vec<T,N>& v1, const su2vec<T,N>& v2);
  template <typename T, int N> bool operator<  (const su2vec<T,N>& v1, T value);
  template <typename T, int N> bool operator>  (const su2vec<T,N>& v1, const su2vec<T,N>& v2);
  template <typename T, int N> bool operator>  (const su2vec<T,N>& v1, T value);
  template <typename T, int N> bool operator<= (const su2vec<T,N>& v1, const su2vec<T,N>& v2);
  template <typename T, int N> bool operator<= (const su2vec<T,N>& v1, T value);
  template <typename T, int N> bool operator>= (const su2vec<T,N>& v1, const su2vec<T,N>& v2);
  template <typename T, int N> bool operator>= (const su2vec<T,N>& v1, T value);
  template <typename T, int N> void copy       (const su2vec<T,N>& v1, su2vec<T,N>& v2);

//////////////////////////////////////////////////////////////////////////////

  /* Definition of a class su2vec that implements an expression template technique
  /* with statically defined size
  / @author G. Orlando
  */

template <typename T, int N = 0>
class su2vec : public ExprT<su2vec<T,N>,tuple_to_be_added> {
public:

  /// Default constructor
  HHOST_DEV su2vec() : EETVEC(ExprT,su2vec,T,N)(this) {this->initialize(T());}

  /// Copy constructor from expression
  template <typename EXPR>
  HHOST_DEV su2vec(EETYPEV(EXPR) expr) : EETVEC(ExprT,su2vec,T,N)(this) {EVECLOOP(i, 0, N, m_data[i] = expr.at(i));}

  /// Copy constructor from su2vec
  HHOST_DEV su2vec(const su2vec<T,N>& orig) : EETVEC(ExprT,su2vec,T,N)(this) {copy(orig,*this);}

  /// Default destructor
  HHOST_DEV ~su2vec() {}

  /// Overloading of assignment operator(s)
#define VEC_EASSIGN_OP(__op__) \
  template <typename EXPR>						\
  HHOST_DEV const su2vec<T,N>& operator __op__ (EETYPEV(EXPR) expr) {EVECLOOP(i, 0, N, m_data[i] __op__ expr.at(i)); return *this;}
VEC_EASSIGN_OP(=)
VEC_EASSIGN_OP(+=)
VEC_EASSIGN_OP(-=)
VEC_EASSIGN_OP(*=)
VEC_EASSIGN_OP(/=)
#undef VEC_EASSIGN_OP

  /// Overloading of assignment operator(s) with constants
#define VEC_EASSIGN_OP_CONST(__op__) \
  HHOST_DEV const su2vec<T,N>& operator __op__ (T expr) {EVECLOOP(i, 0, N, m_data[i] __op__ expr); return *this;}
VEC_EASSIGN_OP_CONST(=)
VEC_EASSIGN_OP_CONST(+=)
VEC_EASSIGN_OP_CONST(-=)
VEC_EASSIGN_OP_CONST(*=)
VEC_EASSIGN_OP_CONST(/=)
#undef VEC_EASSIGN_OP_CONST

  /// copy content of another array
  template <typename ARRAY>
  void copyFrom(const ARRAY& in) {EVECLOOP(i, 0, N, m_data[i] = in[i]);}

  /// copy content of another array
  template <typename ARRAY>
  void copyTo(ARRAY& in) {EVECLOOP(i, 0, N, in[i] = m_data[i]);}

  /// @return a vector slice with fixed size
  template <int NV>
  su2vecSlice<T,NV> slice(size_t start) {return su2vecSlice<T,NV>(&m_data[start]);}

  /// @return a vector slice
  su2vecSlice<T,0> slice(size_t start, size_t ns) {return su2vecSlice<T,0>(&m_data[start], ns);}

  /// @return the raw data
  T* ptr() {return &m_data[0];}

  /// return the array size
  size_t size() const {return N;}

  /// Overloading of the stream operator "<<" for the output.
  /// "\n"ine introduced at the end of every line of the matrix.
  friend std::ostream& operator<<
  LTGT (std::ostream& out, const su2vec<T,N>& v);

  /// Overloading of the stream operator ">>" for the input
  friend std::istream& operator>>
  LTGT (std::istream& in, su2vec<T,N>& v);

  /// Overloading of the "==" operator.
  /// @return true if all elements are equal elementwise
  friend bool operator== LTGT (const su2vec<T,N>& v1, const su2vec<T,N>& v2);

  /// Overloading of the "==" operator.
  /// @param v given array
  /// @param value value for the comparison
  /// @return true if all elements are equal to value
  friend bool operator== LTGT (const su2vec<T,N>& v1, T value);

  /// Overloading of the "!=" operator.
  /// @return true if all elements are not equal elementwise
  friend bool operator!= LTGT (const su2vec<T,N>& v1, const su2vec<T,N>& v2);

  /// Overloading of the "!=" operator.
  /// @param v given array
  /// @param value value for the comparison
  /// @return true if at least one element is not equal to value
  friend bool operator!= LTGT (const su2vec<T,N>& v1, T value);

  /// Overloading of the "<" operator.
  /// @return true if the norm of the first su2vec is < than the norm of the second su2vec.
  friend bool operator< LTGT (const su2vec<T,N>& v1, const su2vec<T,N>& v2);

  /// Overloading of the "<" operator.
  /// @return true if all entries of the first su2vec are < than the given value
  friend bool operator< LTGT (const su2vec<T,N>& v1, T value);

  /// Overloading of the ">" operator.
  /// @return true if the norm of the first su2vec is > than the norm of the second su2vec.
  friend bool operator> LTGT (const su2vec<T,N>& v1, const su2vec<T,N>& v2);

  /// Overloading of the ">" operator.
  /// @return true if all entries of the first su2vec are > than the given value
  friend bool operator> LTGT (const su2vec<T,N>& v1, T value);

  /// Overloading of the "<=" operator.
  /// @return true if the norm of the first su2vec is <= than the norm of the second su2vec.
  friend bool operator<= LTGT (const su2vec<T,N>& v1, const su2vec<T,N>& v2);

  /// Overloading of the "<=" operator.
  /// @return true if all entries of the first su2vec are <= than the given value
  friend bool operator<= LTGT (const su2vec<T,N>& v1, T value);

  /// Overloading of the ">=" operator.
  /// @return true if the norm of the first su2vec is >= than the norm of the second su2vec.
  friend bool operator>= LTGT (const su2vec<T,N>& v1, const su2vec<T,N>& v2);

  /// Overloading of the ">=" operator.
  /// @return true if all entries of the first su2vec are >= than the given value
  friend bool operator>= LTGT (const su2vec<T,N>& v1, T value);

  /// Copy one su2vec into another one
  /// @pre v1.size() == v2.size()
  /// @param v1 source vector
  /// @param v2 destination vector
  friend void copy LTGT (const su2vec<T,N>& orig, su2vec<T,N>& dest);

  /// Normalizes this vector so that the norm is equal to one.
  /// \f$ v = v/\|v\|\f$
  void normalize() {*this /= this->norm2();}

  /// Projection of one vector onto another.
  /// @param v1 1st su2vector
  /// @param v2 2nd su2vector
  /// @post this su2vec contains the projected vector of v2 onto v1.
  template <typename V1, typename V2>
  void proj(const V1& v1, const V2& v2)
  {
    assert(v1.size() = v2.size());
    T scale = (MathFunctions::innerProd(v1, v2) / v1.sqrNorm());
    *this = v1 * scale;
  }

  /// Projection of this su2vector onto a given one.
  /// @pre this and the given object must be of the same size.
  /// @param v1 the other su2vector
  /// @post this su2vector contains the projected vector of itself onto v1.
  template <typename V1>
  void proj(const V1& v1) { proj(*this,v1);}

private:

  /// array data
  T m_data[N];
};
