#ifndef SU2_NUMERICS_REACTIVE
#define SU2_NUMERICS_REACTIVE

#include "numerics_structure.hpp"
#include "variable_reactive.hpp"

/*!
 * \class CUpwReactiveAUSM
 * \brief Class for computing convective flux using AUSM+-up method for multispecies flows.
 * \author G. Orlando
 */
class CUpwReactiveAUSM: public CNumerics {
public:
  using RealVec = CReactiveEulerVariable::RealVec;
  using LibraryPtr = CReactiveEulerVariable::LibraryPtr;

protected:
  LibraryPtr library; /*!< \brief Smart pointer to the library that computes physical-chemical properties. */

  bool implicit; /*!< \brief Flag for implicit scheme. */

  unsigned short nSpecies; /*!< \brief Number of species in the mixture. */

  su2double mInfty;  /*!< \brief Freestream Mach number used to compute the reference Mach number. */

  RealVec Phi_i,
          Phi_j;  /*!< \brief Vectors to describe the variables of the problem used in the AUSM scheme. */

private:
  unsigned short T_INDEX_PRIM, VX_INDEX_PRIM,
                 P_INDEX_PRIM, RHO_INDEX_PRIM,
                 H_INDEX_PRIM, A_INDEX_PRIM,
                 RHOS_INDEX_PRIM;               /*!< \brief Mapping for position in primitives array. */

  unsigned short RHO_INDEX_SOL, RHOVX_INDEX_SOL,
                 RHOE_INDEX_SOL, RHOS_INDEX_SOL; /*!< \brief Mapping for position in conserved array. */

  unsigned short T_INDEX_GRAD, VX_INDEX_GRAD,
                 P_INDEX_GRAD;                 /*!< \brief Mapping for position in primitives gradient. */

  unsigned short T_INDEX_LIM, VX_INDEX_LIM,
                 P_INDEX_LIM;                 /*!< \brief Mapping for position for limited variables. */

public:
  /*!
   * \brief Default constructor of the class.
   */
  CUpwReactiveAUSM(): CNumerics(), implicit(), nSpecies(), mInfty(), T_INDEX_PRIM(), VX_INDEX_PRIM(), P_INDEX_PRIM(), RHO_INDEX_PRIM(),
                      H_INDEX_PRIM(), A_INDEX_PRIM(), RHOS_INDEX_PRIM(), RHO_INDEX_SOL(), RHOVX_INDEX_SOL(), RHOE_INDEX_SOL(),
                      RHOS_INDEX_SOL(), T_INDEX_GRAD(), VX_INDEX_GRAD(), P_INDEX_GRAD(), T_INDEX_LIM(), VX_INDEX_LIM(), P_INDEX_LIM() {}

  /*!
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
   * \param[in] lib_ptr - Pointer to the external library for physical-chemical properties
	 */
  CUpwReactiveAUSM(unsigned short val_nDim, unsigned short val_nVar, CConfig* config, LibraryPtr lib_ptr);

  /*!
	 * \brief Destructor of the class.
	 */
	virtual ~CUpwReactiveAUSM() {}

  /*!
	 * \brief Compute the residual associated to convective flux between two nodes i and j.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void ComputeResidual(su2double* val_residual, su2double** val_Jacobian_i, su2double** val_Jacobian_j, CConfig* config) override;

};

/*!
 * \class CAvgGradReactive_Boundary
 * \brief Class for computing viscous flux using the average of gradients for multispecies flows for boundary nodes.
 * \author G. Orlando
 */
class CAvgGradReactive_Boundary: public CNumerics {
public:
  using RealVec = CReactiveEulerVariable::RealVec;
  using RealMatrix = CReactiveNSVariable::RealMatrix;
  using LibraryPtr = CReactiveEulerVariable::LibraryPtr;
  using Vec = Eigen::VectorXd;

protected:
  LibraryPtr library; /*!< \brief Smart pointer to the library that computes physical-chemical properties. */

  bool implicit;    /*!< \brief Flag for implicit computations. */

  unsigned short nPrimVar; /*!< \brief Numbers of primitive variables. */
  unsigned short nPrimVarAvgGrad; /*!< \brief Numbers of primitive variables to compute gradient for average computation. */
  unsigned short nSpecies; /*!< \brief Total number of species. */

  su2double alpha;        /*!< \brief Artificial diffusion coefficient for Stefan-Maxwell equations. */

  Vec PrimVar_i,          /*!< \brief Primitive variables at node i. */
      PrimVar_j,          /*!< \brief Primitive variables at node j. */
      Mean_PrimVar;       /*!< \brief Mean primitive variables. */

  RealMatrix Mean_GradPrimVar;    /*!< \brief Mean value of the gradient. */

  RealVec Xs_i,               /*!< \brief Auxiliary vector for mole fractions at point i. */
          Xs_j;               /*!< \brief Auxiliary vector for mole fractions at point j. */

  RealVec Xs,                 /*!< \brief Auxiliary vector for mean mole fractions. */
          Ys;                 /*!< \brief Auxiliary vector for mean mass fractions. */

  unsigned short T_INDEX_PRIM, VX_INDEX_PRIM,
                 P_INDEX_PRIM, RHO_INDEX_PRIM,
                 H_INDEX_PRIM, A_INDEX_PRIM,
                 RHOS_INDEX_PRIM;               /*!< \brief Mapping for position in primitives array. */

  unsigned short RHO_INDEX_SOL, RHOVX_INDEX_SOL,
                 RHOE_INDEX_SOL, RHOS_INDEX_SOL; /*!< \brief Mapping for position in conserved array. */

  unsigned short T_INDEX_GRAD, VX_INDEX_GRAD,
                 P_INDEX_GRAD, RHOS_INDEX_GRAD;  /*!< \brief Mapping for position in primitives gradient. */

  unsigned short T_INDEX_LIM, VX_INDEX_LIM,
                 P_INDEX_LIM;                 /*!< \brief Mapping for position for limited variables. */

  RealMatrix Dij_i,             /*!< \brief Binary diffusion coefficients at point i. */
             Dij_j,             /*!< \brief Binary diffusion coefficients at point j. */
             Mean_Dij;          /*!< \brief Harmonic average of binary diffusion coefficients at point i and j. */

  RealMatrix Gamma,
             Gamma_tilde;  /*!< \brief Auxiliary matrices for solving Stefan-Maxwell equations. */

  RealVec hs,                   /*!< \brief Auxiliary vector to store partial enthalpy for species diffusion flux contribution. */
          Cps;                  /*!< \brief Auxiliary vector to store Cp for species diffusion flux Jacobian contribution. */

  Vec Jd;                       /*!< \brief Auxiliary vector to store S-M solution. */

  Vec Grad_Xs_norm;            /*!< \brief Auxiliary vector to store normal gradient of mole fractions. */

  Vec Ds_i,                    /*!< \brief Auxiliary vector to store Ramshaw diffusion coefficients at node i. */
      Ds_j,                    /*!< \brief Auxiliary vector to store Ramshaw diffusion coefficients at node j. */
      Ds;                      /*!< \brief Auxiliary vector to store average Ramshaw diffusion coefficients. */

  unsigned short T_INDEX_AVGGRAD,
                 VX_INDEX_AVGGRAD,
                 RHOS_INDEX_AVGGRAD;  /*!< \brief Mapping between the primitive variable gradient name
                                                  for average computations and its position in the physical data. */

public:
  /*!
   * \brief Default constructor of the class.
   */
  CAvgGradReactive_Boundary(): CNumerics(), implicit(), nPrimVar(), nPrimVarAvgGrad(), nSpecies(), T_INDEX_PRIM(), VX_INDEX_PRIM(),
                               P_INDEX_PRIM(), RHO_INDEX_PRIM(), H_INDEX_PRIM(), A_INDEX_PRIM(), RHOS_INDEX_PRIM(), RHO_INDEX_SOL(),
                               RHOVX_INDEX_SOL(), RHOE_INDEX_SOL(), RHOS_INDEX_SOL(), T_INDEX_GRAD(), VX_INDEX_GRAD(), P_INDEX_GRAD(),
                               RHOS_INDEX_GRAD(), T_INDEX_LIM(), VX_INDEX_LIM(), P_INDEX_LIM(), T_INDEX_AVGGRAD(), VX_INDEX_AVGGRAD(),
                               RHOS_INDEX_AVGGRAD() {}

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] lib_ptr - Pointer to the external library for physical-chemical properties
   */
  CAvgGradReactive_Boundary(unsigned short val_nDim, unsigned short val_nVar, CConfig* config, LibraryPtr lib_ptr);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CAvgGradReactive_Boundary() = default;

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
  inline void SetPrimVarGradient(unsigned short val_nvar, unsigned short val_nDim, su2double val_grad_i, su2double val_grad_j) {
    SU2_Assert(PrimVar_Grad_i[val_nvar] != NULL,
               std::string("The row " + std::to_string(val_nvar) + " of gradient primitive has not been allocated"));
    SU2_Assert(PrimVar_Grad_j[val_nvar] != NULL,
               std::string("The row " + std::to_string(val_nvar) + " of gradient primitive has not been allocated"));
    PrimVar_Grad_i[val_nvar][val_nDim] = val_grad_i;
    PrimVar_Grad_j[val_nvar][val_nDim] = val_grad_j;
  }

  /*!
   * \brief Set the binary diffusion coefficients at node i and j.
   * \param[in] Diff_i - Binary diffusion coefficients at node i.
   * \param[in] Diff_j - Binary diffusion coefficients at node j.
   */
  inline void SetBinaryDiffCoeff(const RealMatrix& Diff_i, const RealMatrix& Diff_j) {
    Dij_i = Diff_i;
    Dij_j = Diff_j;
  }

protected:
  /*!
   * \brief Compute projection of the viscous fluxes using Ramshaw self-consistent modification.
   * \param[in] val_primvar - Primitive variables.
   * \param[in] val_grad_primvar - Gradient of the primitive variables.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] val_viscosity - Laminar viscosity.
   * \param[in] val_thermal_conductivity - Thermal Conductivity.
   * \param[in] val_diffusioncoeff - Effective diffusion coefficients for each species in the mixture.
   * \param[in] config - Definition of the particular problem
   */
  void GetViscousProjFlux(const Vec& val_primvar, const RealMatrix& val_grad_primvar, su2double* val_normal,
                          const su2double val_viscosity, const su2double val_therm_conductivity,
                          const Vec& val_diffusioncoeff, CConfig* config);

  /*!
   * \brief Compute projection of the viscous fluxes solving Stefan-Maxwell equations.
   * \param[in] val_primvar - Primitive variables.
   * \param[in] val_grad_primvar - Gradient of the primitive variables.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] val_viscosity - Laminar viscosity.
   * \param[in] val_thermal_conductivity - Thermal Conductivity.
   * \param[in] val_Dij - Harmonic average of binary diffusion coefficients.
   * \param[in] config - Definition of the particular problem
   */
  void GetViscousProjFlux(const Vec& val_primvar, const RealMatrix& val_grad_primvar, su2double* val_normal,
                          const su2double val_viscosity, const su2double val_therm_conductivity,
                          const RealMatrix& val_Dij, CConfig* config);

  /*!
   * \brief Approximation of Viscous NS Jacobians in Thermochemical Non Equilibrium.
   * \param[in] val_Mean_PrimVar - Mean value of the primitive variables.
   * \param[in] val_laminar_viscosity - Value of the laminar viscosity.
   * \param[in] val_thermal_conductivity - Value of the thermal conductivity.
   * \param[in] val_alpha - Value of artificual diffusion coefficient to solve Stefan-Maxwell equations.
   * \param[in] val_grad_xs_norm - Value of normal gradient of mole fractions.
   * \param[in] val_diffusion_coeff - Value of binary diffusion coefficients at interface.
   * \param[in] val_dist_ij - Distance between the points.
   * \param[in] val_dS - Area of the current face.
   * \param[in] val_normal - Normal vector
   * \param[in] val_Proj_Visc_Flux - Pointer to the projected viscous flux.
   * \param[out] val_Proj_Jac_Tensor_i - Pointer to the projected viscous Jacobian at point i.
   * \param[out] val_Proj_Jac_Tensor_j - Pointer to the projected viscous Jacobian at point j.
   * \param[in] config - Definition of the particular problem
  */
  void GetViscousProjJacs(const Vec& val_Mean_PrimVar, const su2double val_laminar_viscosity, const su2double val_thermal_conductivity,
                          const su2double val_alpha, const Vec& val_grad_xs_norm, const Vec& val_diffusion_coeff,
                          const su2double val_dist_ij, const su2double val_dS, su2double* val_normal, su2double* val_Proj_Visc_Flux,
                          su2double** val_Proj_Jac_Tensor_i, su2double** val_Proj_Jac_Tensor_j, CConfig* config);

  /*!
   * \brief Compute the diffusive flux along a certain direction
   * \param[in] val_density - Density of the mixture.
   * \param[in] val_alpha - Parameter for artifical diffusion.
   * \param[in] val_Dij - Harmonic average of binary diffusion coefficients.
   * \param[in] val_xs - Molar fractions.
   * \param[in] val_grad_xs - Component along the desired dimension of gradient of molar fractions.
   * \param[in] val_ys - Mass fractions.
   */
  void Solve_SM(const su2double val_density, const su2double val_alpha, const RealMatrix& val_Dij,
                const RealVec& val_xs, const Vec& val_grad_xs, const RealVec& val_ys);
};

/*!
 * \class CAvgGradReactive_Flow
 * \brief Class for computing viscous flux using the corrected average gradients for multispecies flows for internal nodes.
 * \author G. Orlando
 */
class CAvgGradReactive_Flow: public CAvgGradReactive_Boundary {
protected:
  bool limiter;     /*!< \brief Flag for limiter computations. */

private:
  Vec Edge_Vector;        /*!< \brief Vector connecting point i to j. */

  Vec Diff_PrimVar;       /*!< \brief Difference of primitive varaibles involved in average computations. */

  Vec Proj_Mean_GradPrimVar_Edge; /*!< \brief Projected mean value of the gradient. */

public:
  /*!
   * \brief Default constructor of the class.
   */
  CAvgGradReactive_Flow(): CAvgGradReactive_Boundary(), limiter() {}

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] lib_ptr - Pointer to the external library for physical-chemical properties
   */
  CAvgGradReactive_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig* config, LibraryPtr lib_ptr);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CAvgGradReactive_Flow() = default;

protected:
  /*!
   * \brief Compute the viscous flow residual using corrected average of gradients method.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double* val_residual, su2double** val_Jacobian_i, su2double** val_Jacobian_j, CConfig* config) override;
};

/*!
 * \class CSourceReactive
 * \brief Class for computing residual term due to chemistry.
 * \author G. Orlando
 */
class CSourceReactive: public CNumerics {
public:
  using RealVec = CReactiveEulerVariable::RealVec;
  using LibraryPtr = CReactiveEulerVariable::LibraryPtr;
  using RealMatrix = CReactiveNSVariable::RealMatrix;

protected:
  LibraryPtr library; /*!< \brief Smart pointer to the library that computes physical-chemical properties. */

  bool implicit;    /*!< \brief Flag for implicit scheme. */

  unsigned short nSpecies;  /*!< \brief Number of species in the mixture. */

  RealVec Ys; /*!< \brief Auxiliary vector for mass fractions in the mixture. */

private:
  unsigned short T_INDEX_PRIM, VX_INDEX_PRIM,
                 P_INDEX_PRIM, RHO_INDEX_PRIM,
                 H_INDEX_PRIM, A_INDEX_PRIM,
                 RHOS_INDEX_PRIM;               /*!< \brief Mapping for position in primitives array. */

  unsigned short RHO_INDEX_SOL, RHOVX_INDEX_SOL,
                 RHOE_INDEX_SOL, RHOS_INDEX_SOL; /*!< \brief Mapping for position in conserved array. */

  unsigned short T_INDEX_GRAD, VX_INDEX_GRAD,
                 P_INDEX_GRAD;                 /*!< \brief Mapping for position in primitives gradient. */

  unsigned short T_INDEX_LIM, VX_INDEX_LIM,
                 P_INDEX_LIM;                 /*!< \brief Mapping for position for limited variables. */

private:
  RealVec omega; /*!< \brief Auxiliary vector for mass production term. */

  RealMatrix source_jac; /*!< \brief Auxiliary vector for mass production term. */

public:
  /*!
   * \brief Default constructor of the class.
   */
  CSourceReactive(): CNumerics(), implicit(), nSpecies(), T_INDEX_PRIM(), VX_INDEX_PRIM(), P_INDEX_PRIM(), RHO_INDEX_PRIM(),
                     H_INDEX_PRIM(), A_INDEX_PRIM(), RHOS_INDEX_PRIM(), RHO_INDEX_SOL(), RHOVX_INDEX_SOL(), RHOE_INDEX_SOL(),
                     RHOS_INDEX_SOL(), T_INDEX_GRAD(), VX_INDEX_GRAD(), P_INDEX_GRAD(), T_INDEX_LIM(), VX_INDEX_LIM(), P_INDEX_LIM() {}

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] lib_ptr - Pointer to the external library for physical-chemical properties
   */
  CSourceReactive(unsigned short val_nDim, unsigned short val_nVar, CConfig* config, LibraryPtr lib_ptr);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CSourceReactive() {}

  /*!
   * \brief Calculation of the residual of chemistry source term
   * \param[out] val_residual - Residual of the source terms
   * \param[out] val_Jacobian_i - Jacobian of the source terms
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeChemistry(su2double* val_residual, su2double** val_Jacobian_i, CConfig* config) override;

};

#endif
