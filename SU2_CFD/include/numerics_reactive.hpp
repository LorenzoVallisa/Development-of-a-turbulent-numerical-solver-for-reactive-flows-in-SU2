#ifndef SU2_NUMERICS_REACTIVE
#define SU2_NUMERICS_REACTIVE

#include "numerics_structure.hpp"
#include "../../Common/include/physical_property_library.hpp"
#include "../../Common/include/reacting_model_library.hpp"
#include "variable_reactive.hpp"

#include <memory>

/*!
 * \class CUpwReactiveAUSM
 * \brief Class for solving using AUSM method.
 * \author G. Orlando
 */
class CUpwReactiveAUSM: public CNumerics {
public:

  using  RealVec = CReactiveEulerVariable::RealVec;
  using  RealMatrix = CReactiveEulerVariable::RealMatrix;
  using  LibraryPtr = CReactiveEulerVariable::LibraryPtr;

protected:

  LibraryPtr library; /*!< \brief Smart pointer to the library that computes physical-chemical properties. */

  bool implicit; /*!< \brief Flag for implicit scheme. */

  const unsigned short nSpecies; /*!< \brief Number of species in the mixture. */

  RealVec   Velocity_i,Velocity_j; /*!< \brief Velocity at node i and at node j. */
  RealVec   ProjFlux_i,ProjFlux_j; /*!< \brief Projected velocities. */

public:

  /*!
   * \brief Default constructor of the class.
   */
  CUpwReactiveAUSM():CNumerics(),implicit(),nSpecies() {}

  /*!
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
  CUpwReactiveAUSM(unsigned short val_nDim, unsigned short val_nVar, std::unique_ptr<CConfig> config);

  /*!
	 * \brief Destructor of the class.
	 */
	virtual ~CUpwReactiveAUSM() {};

  /*!
	 * \brief Compute the flux between two nodes i and j.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void ComputeResidual(su2double* val_residual, su2double** val_Jacobian_i,
                       su2double** val_Jacobian_j, CConfig* config) override;

};


/*!
 * \class CAvgGradReactive_Flow
 * \brief Class for computing viscous term using the average of gradients for a chemically reactive flow.
 * \author G. Orlando
 */

class CAvgGradReactive_Flow : public CNumerics {
public:

  typedef std::vector<su2double> RealVec;
  typedef std::vector<RealVec>  RealMatrix;
  typedef std::unique_ptr<Framework::PhysicalPropertyLibrary> LibraryPtr;

protected:

  LibraryPtr library; /*!< \brief Smart pointer to the library that computes physical-chemical properties. */

  bool implicit; /*!< \brief Flag for implicit computations. */
  bool limiter; /*!< \brief Flag for limiter computations. */

  const unsigned short nSpecies; /*!< \brief Total number of species. */

  RealVec Mean_PrimVar;           /*!< \brief Mean primitive variables. */
  RealVec PrimVar_i, PrimVar_j;   /*!< \brief Primitives variables at point i and j. */
  RealMatrix Mean_GradPrimVar;    /*!< \brief Mean value of the gradient. */
  RealVec Proj_Mean_GradPrimVar_Edge; /*!< \brief Mean gradient projection. */

public:

  /*!
   * \brief Default constructor of the class.
   */
  CAvgGradReactive_Flow():CNumerics(),implicit(),limiter(),nSpecies() {}

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGradReactive_Flow(unsigned short val_nDim, unsigned short val_nVar, std::unique_ptr<CConfig> config);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CAvgGradReactive_Flow() {};

 /*!
   * \brief Compute the viscous flow residual using an average of gradients.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double* val_residual, su2double** val_Jacobian_i,
                       su2double** val_Jacobian_j, CConfig* config) override;
};


class CSourceReactive: public CNumerics {

protected:

  bool implicit; /*!< \brief Flag for implicit computations. */

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceReactive(unsigned short val_nDim, unsigned short val_nVar, std::unique_ptr<CConfig> config);

  /*!
   * \brief Destructor of the class.
   */
  ~CSourceReactive() {};

  /*!
   * \brief Calculation of the chemistry source term
   * \param[in] config - Definition of the particular problem.
   * \param[out] val_residual - residual of the source terms
   * \param[out] val_Jacobian_i - Jacobian of the source terms
   */
  void ComputeChemistry(su2double* val_residual, su2double** val_Jacobian_i, CConfig* config) override;


  /*!
   * \brief Residual for source term integration.
   * \param[out] val_residual - Pointer to the source residual containing chemistry terms.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual_Axisymmetric(su2double* val_residual, CConfig* config) override;

};


#endif
