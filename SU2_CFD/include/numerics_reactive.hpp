#ifndef SU2_NUMERICS_REACTIVE
#define SU2_NUMERICS_REACTIVE

#include "numerics_structure.hpp"
#include "variable_reactive.hpp"

/*!
 * \class CUpwReactiveAUSM
 * \brief Class for solving using AUSM method.
 * \author G. Orlando
 */
class CUpwReactiveAUSM: public CNumerics {
public:

  using RealVec = CReactiveEulerVariable::RealVec;
  using RealMatrix = CReactiveEulerVariable::RealMatrix;
  using LibraryPtr = CReactiveEulerVariable::LibraryPtr;

protected:

  LibraryPtr library; /*!< \brief Smart pointer to the library that computes physical-chemical properties. */

  bool implicit; /*!< \brief Flag for implicit scheme. */

  const unsigned short nSpecies; /*!< \brief Number of species in the mixture. */

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
  CUpwReactiveAUSM(unsigned short val_nDim, unsigned short val_nVar, CConfig* config);

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

  /*!
   * \brief Set the simulation to explicit
   */
  inline void SetExplicit(void) {
    implicit = false;
  }

};


/*!
 * \class CAvgGradReactive_Flow
 * \brief Class for computing viscous term using the average of gradients for a chemically reactive flow.
 * \author G. Orlando
 */

class CAvgGradReactive_Flow : public CNumerics {
public:

  using RealVec = CReactiveEulerVariable::RealVec;
  using RealMatrix = CReactiveEulerVariable::RealMatrix;
  using LibraryPtr = CReactiveEulerVariable::LibraryPtr;
  using SmartArr = CReactiveEulerVariable::SmartArr;


protected:

  LibraryPtr library; /*!< \brief Smart pointer to the library that computes physical-chemical properties. */

  bool implicit; /*!< \brief Flag for implicit computations. */
  bool limiter; /*!< \brief Flag for limiter computations. */

  unsigned short nPrimVar; /*!< \brief Numbers of primitive variables. */
  unsigned short nPrimVarAvgGrad; /*!< \brief Numbers of primitive variables to compute gradient for average computation. */
  const unsigned short nSpecies; /*!< \brief Total number of species. */

  RealVec Mean_PrimVar;           /*!< \brief Mean primitive variables. */
  RealVec PrimVar_i, PrimVar_j;   /*!< \brief Primitives variables at point i and j. */
  RealMatrix GradPrimVar_i, GradPrimVar_j; /*!< \brief Gradient of primitives variables at point i and j for average gradient computation. */
  RealMatrix Mean_GradPrimVar;    /*!< \brief Mean value of the gradient. */

public:

  /*!
   * \brief Default constructor of the class.
   */
  CAvgGradReactive_Flow():CNumerics(),implicit(),limiter(),nPrimVar(),nPrimVarAvgGrad(),nSpecies() {}

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGradReactive_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig* config);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CAvgGradReactive_Flow() {};

  /*!
   * \brief Compute projection of the viscous fluxes into a direction
   * \param[in] val_primvar - Primitive variables.
   * \param[in] val_grad_primvar - Gradient of the primitive variables.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] val_diffusioncoeff - Diffusion coefficients for each species in the mixture.
   * \param[in] val_viscosity - Laminar viscosity.
   * \param[in] val_thermal_conductivity - Thermal Conductivity.
   * \param[in] config - Configuration file
   */
  virtual void GetViscousProjFlux(RealVec& val_primvar, RealMatrix& val_grad_primvar,
                                  SmartArr val_normal, RealVec& val_diffusioncoeff,
                                  su2double val_viscosity,su2double val_therm_conductivity);


  /*!
   * \brief Approximation of Viscous NS Jacobians in Thermochemical Non Equilibrium.
   * \param[in] val_Mean_PrimVar - Mean value of the primitive variables.
   * \param[in] val_Mean_GradPriVar - Mean value of the gradient of the primitive variables.
   * \param[in] val_diffusion_coeff - Value of diffusion coefficients for each species
   * \param[in] val_laminar_viscosity - Value of the laminar viscosity.
   * \param[in] val_thermal_conductivity - Value of the thermal conductivity.
   * \param[in] val_dist_ij - Distance between the points.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] val_dS - Area of the face between two nodes.
   * \param[in] val_Proj_Visc_Flux - Pointer to the projected viscous flux.
   * \param[out] val_Proj_Jac_Tensor_i - Pointer to the projected viscous Jacobian at point i.
   * \param[out] val_Proj_Jac_Tensor_j - Pointer to the projected viscous Jacobian at point j.
   * \param[in] config - Definition of the particular problem
  */
  virtual void GetViscousProjJacs(RealVec& val_Mean_PrimVar, RealMatrix& val_Mean_GradPrimVar,
                                  RealVec& val_diffusion_coeff, su2double val_laminar_viscosity,su2double val_thermal_conductivity,
                                  su2double val_dist_ij, SmartArr val_normal, su2double val_dS,
                                  su2double* val_Proj_Visc_Flux,su2double** val_Proj_Jac_Tensor_i,su2double** val_Proj_Jac_Tensor_j,
                                  CConfig* config);

 /*!
   * \brief Compute the viscous flow residual using an average of gradients.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double* val_residual, su2double** val_Jacobian_i,
                       su2double** val_Jacobian_j, CConfig* config) override;

  /*!
   * \brief Set the Gradient of primitives for computing residual
   */
  inline void SetGradient_AvgPrimitive(const RealMatrix& Grad_i,const RealMatrix& Grad_j) {
    GradPrimVar_i = Grad_i;
    GradPrimVar_j = Grad_j;
  }

  /*!
   * \brief Set the simulation to explicit
   */
  inline void SetExplicit(void) {
    implicit = false;
  }

};


class CSourceReactive: public CNumerics {
public:

  using RealVec = CReactiveEulerVariable::RealVec;
  using RealMatrix = CReactiveEulerVariable::RealMatrix;
  using LibraryPtr = CReactiveEulerVariable::LibraryPtr;

protected:

  LibraryPtr library; /*!< \brief Smart pointer to the library that computes physical-chemical properties. */

  bool implicit; /*!< \brief Flag for implicit scheme. */

  const unsigned short nSpecies; /*!< \brief Number of species in the mixture. */

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceReactive(unsigned short val_nDim, unsigned short val_nVar, CConfig* config);

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
  //void ComputeResidual_Axisymmetric(su2double* val_residual, CConfig* config) override;

  /*!
   * \brief Set the simulation to explicit
   */
  inline void SetExplicit(void) {
    implicit = false;
  }

};


#endif
