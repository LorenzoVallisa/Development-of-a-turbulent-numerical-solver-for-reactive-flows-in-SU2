#ifndef SU2_SPLINE_HPP
#define SU2_SPLINE_HPP

#include <vector>

namespace MathTools {

  typedef std::vector<double> RealVec;

  /*!
    * \brief Sets the second derivative coefficients for natural spline.
    * \param[in] x - The vector of values to be tabulated (input)
    * \param[in] y - The vector of values tabulated,i.e y_i = f(x_i) (input)
    * \param[in] yp1 - Value of first derivative in the first point (input)
    * \param[in] ypn - Value of first derivative in the last point (input)
    * \param[out] y2 - Coefficients of second derivative for interpolating function (output)
  */
  void SetSpline(const RealVec& x, const RealVec& y, const double yp1, const double ypn, RealVec& y2);

  /*!
    * \brief Get the value of a property for a certan temperature through natural spline.
    * \param[in] x - The vector of values to be tabulated
    * \param[in] y - The vector of values tabulated,i.e yi = f(xi)
    * \param[in] y2 - Coefficients of second derivative for interpolating function
    * \param[in] value - The value we want to tabulate through spline
  */
  double GetSpline(const RealVec& x, const RealVec& y, const RealVec& y2, const double value);
}

#endif
