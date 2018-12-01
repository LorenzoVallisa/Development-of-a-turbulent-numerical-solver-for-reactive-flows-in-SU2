#ifndef SU2_NUMERICS_REACTIVE
#define SU2_NUMERICS_REACTIVE

#include "numerics_structure.hpp"
#include "variable_reactive.hpp"

/*!
 * \class CUpwReactiveAUSM
 * \brief Class for computing convective flux using AUSM method.
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

  unsigned short nSpecies; /*!< \brief Number of species in the mixture. */

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
	 * \brief Compute the residual associated to convective flux between two nodes i and j.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void ComputeResidual(su2double* val_residual, su2double** val_Jacobian_i, su2double** val_Jacobian_j, CConfig* config) override;

  /*!
   * \brief Set the simulation to explicit
   */
  inline void SetExplicit(void) {
    implicit = false;
  }

};


/*!
 * \class CAvgGradReactive_Flow
 * \brief Class for computing viscous flux using the average gradients for a chemically reactive flow.
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
  unsigned short nSpecies; /*!< \brief Total number of species. */

  RealVec Mean_PrimVar;           /*!< \brief Mean primitive variables. */
  RealVec PrimVar_i, PrimVar_j;   /*!< \brief Primitives variables at point i and j. */
  RealMatrix GradPrimVar_i, GradPrimVar_j; /*!< \brief Gradient of primitives variables at point i and j for average gradient computation. */
  RealMatrix Mean_GradPrimVar;    /*!< \brief Mean value of the gradient. */

public:

  /**
   * \brief Mapping between the primitive variable gradient name and its position in the physical data
   */
  static constexpr unsigned short T_INDEX_GRAD = 0;
  static constexpr unsigned short VX_INDEX_GRAD = 1;
  //static const unsigned short RHO_INDEX_GRAD;
  static const unsigned short RHOS_INDEX_GRAD;


  /**
   * Mapping between the primitive variable gradient name
   * for average computations and its position in the physical data
   */
  static constexpr unsigned short T_INDEX_AVGGRAD = 0;
  static constexpr unsigned short VX_INDEX_AVGGRAD = 1;
  static const unsigned short RHOS_INDEX_AVGGRAD;

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
   * \brief Compute the diffusive flux along a certain direction
   * \param[in] val_density - Density of the mixture
   * \param[in] val_xs - Molar fractions.
   * \param[in] val_grad_xs - Component along the desired dimension of gradient of molar fractions.
   * \param[in] val_diffusioncoeff - Corrected diffusion coefficients for each species ((1-Xs)/Ds).
   */
  RealVec Solve_SM(const su2double val_density, const RealVec& val_xs, const RealVec& val_grad_xs, const RealVec& val_diffusioncoeff);

  /*!
   * \brief Compute projection of the viscous fluxes using Ramshaw self-consisten modification
   * \param[in] val_primvar - Primitive variables.
   * \param[in] val_grad_primvar - Gradient of the primitive variables.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] val_viscosity - Laminar viscosity.
   * \param[in] val_thermal_conductivity - Thermal Conductivity.
   * \param[in] val_diffusioncoeff - Effective diffusion coefficients for each species in the mixture.
   */
  virtual void GetViscousProjFlux_Ramshaw(const RealVec& val_primvar, const RealMatrix& val_grad_primvar, SmartArr val_normal,
                                          const su2double val_viscosity, const su2double val_therm_conductivity, const RealVec& val_diffusioncoeff);

  /*!
   * \brief Compute projection of the viscous fluxes solving Stefan-Maxwell equations
   * \param[in] val_primvar - Primitive variables.
   * \param[in] val_grad_primvar - Gradient of the primitive variables.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] val_viscosity - Laminar viscosity.
   * \param[in] val_thermal_conductivity - Thermal Conductivity.
   * \param[in] val_diffusioncoeff - Effective diffusion coefficients for each species in the mixture.
   */
  virtual void GetViscousProjFlux_SM(const RealVec& val_primvar, const RealMatrix& val_grad_primvar, SmartArr val_normal,
                                     const su2double val_viscosity,const su2double val_therm_conductivity, const RealVec& val_diffusioncoeff);

  /*!
   * \brief Approximation of Viscous NS Jacobians in Thermochemical Non Equilibrium.
   * \param[in] val_Mean_PrimVar - Mean value of the primitive variables.
   * \param[in] val_Mean_GradPriVar - Mean value of the gradient of the primitive variables.
   * \param[in] val_diffusion_coeff - Value of diffusion coefficients for each species
   * \param[in] val_laminar_viscosity - Value of the laminar viscosity.
   * \param[in] val_thermal_conductivity - Value of the thermal conductivity.
   * \param[in] val_dist_ij - Distance between the points.
   * \param[in] val_normal - Normal vector
   * \param[in] val_dS - Area of the face between two nodes.
   * \param[in] val_Proj_Visc_Flux - Pointer to the projected viscous flux.
   * \param[out] val_Proj_Jac_Tensor_i - Pointer to the projected viscous Jacobian at point i.
   * \param[out] val_Proj_Jac_Tensor_j - Pointer to the projected viscous Jacobian at point j.
   * \param[in] config - Definition of the particular problem
  */
  virtual void GetViscousProjJacs(const RealVec& val_Mean_PrimVar, const RealMatrix& val_Mean_GradPrimVar,
                                  const RealVec& val_diffusion_coeff, const su2double val_laminar_viscosity, const su2double val_thermal_conductivity,
                                  const su2double val_dist_ij, SmartArr val_normal, const su2double val_dS, su2double* val_Proj_Visc_Flux,
                                  su2double** val_Proj_Jac_Tensor_i, su2double** val_Proj_Jac_Tensor_j, CConfig* config);

 /*!
   * \brief Compute the viscous flow residual using average gradient method.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double* val_residual, su2double** val_Jacobian_i, su2double** val_Jacobian_j, CConfig* config) override;


  /*!
   * \brief Set the Gradient of primitives for computing residual
   * \param[in] val_nvar - Index of desired variable
   * \param[in] val_nDim - Index of the desired dimension.
   * \param[in] val_grad_i - Value to assign at point i.
   * \param[in] val_grad_j - Value to assign at point j.
   */
  inline void SetGradient_AvgPrimitive(unsigned short val_nvar, unsigned short val_nDim, su2double val_grad_i, su2double val_grad_j) {
    GradPrimVar_i.at(val_nvar,val_nDim) = val_grad_i;
    GradPrimVar_j.at(val_nvar,val_nDim) = val_grad_j;
  }

  /*!
   * \brief Set the Gradient of primitives for computing residual
   * \param[in] Grad_i - Gradient of Primitive variables.
   * \param[in] Grad_j - Gradient of Primitive variables.
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
//const unsigned short CAvgGradReactive_Flow::RHO_INDEX_GRAD = CAvgGradReactive_Flow::VX_INDEX_GRAD + CReactiveNSVariable::GetnDim();
const unsigned short CAvgGradReactive_Flow::RHOS_INDEX_GRAD = CAvgGradReactive_Flow::VX_INDEX_GRAD + CReactiveNSVariable::GetnDim();

const unsigned short CAvgGradReactive_Flow::RHOS_INDEX_AVGGRAD = CAvgGradReactive_Flow::VX_INDEX_AVGGRAD + CReactiveNSVariable::GetnDim();

/*!
 * \class CSourceReactive
 * \brief Class for computing residual due to chemistry source term.
 * \author G. Orlando
 */

class CSourceReactive: public CNumerics {
public:

  using RealVec = CReactiveEulerVariable::RealVec;
  using RealMatrix = CReactiveEulerVariable::RealMatrix;
  using LibraryPtr = CReactiveEulerVariable::LibraryPtr;

protected:

  LibraryPtr library; /*!< \brief Smart pointer to the library that computes physical-chemical properties. */

  bool implicit; /*!< \brief Flag for implicit scheme. */

  unsigned short nSpecies; /*!< \brief Number of species in the mixture. */

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
  virtual ~CSourceReactive() {};

  /*!
   * \brief Calculation of the residual of chemistry source term
   * \param[out] val_residual - Residual of the source terms
   * \param[out] val_Jacobian_i - Jacobian of the source terms
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeChemistry(su2double* val_residual, su2double** val_Jacobian_i, CConfig* config) override;

  /*!
   * \brief Residual for source term integration in case of axisymmetric simulation.
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
