#ifndef SU2_UTILITY_HPP
#define SU2_UTILITY_HPP

#include <set>
#include <map>
#include "su2_assert.hpp"
#include "../../externals/Eigen/Dense"

namespace Utility {

  typedef std::map<std::string,unsigned short> MyMap;
  typedef std::set<std::string> MySet;
  using RealMatrix = Eigen::MatrixXd;

  /*!
    * \brief Read a single chemical reaction
    * \param[in] line - line of the reaction
    * \param[in] n_reac - current number of reaction
    * \param[in] is_rev - flag whether the reaction is considered reversible
    * \param[in] is_reac - true is we are considering the reactants side
    * \param[in] map_names - map with all species name
    * \param[out] stoich_coeff - matrix of stechiometric coefficients
    * \param[out] stoich_coeff_exp_reac - matrix of stechiometric coefficients at the exponent of concetration for reactants
    * \param[out] stoich_coeff_exp_prod - matrix of stechiometric coefficients at the exponent of concetration for products
  */
  void Parse_Terms(std::string& line, unsigned n_reac, bool is_rev, bool is_reac, const MyMap& map_names, RealMatrix& stoich_coeff,
                   RealMatrix& stoich_coeff_exp_reac, RealMatrix& stoich_coeff_exp_prod);

} /*--- End of namespace Utility ---*/

#endif
