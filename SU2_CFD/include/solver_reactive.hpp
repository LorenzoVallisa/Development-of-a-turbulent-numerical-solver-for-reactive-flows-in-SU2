#ifndef SU2_SOLVER_REACTIVE
#define SU2_SOLVER_REACTIVE

#include "solver_structure.hpp"
#include "variable_reactive.hpp"

#include <memory>

/*! \class CReactiveEulerSolver
 *  \brief Main class for defining a solver for chemically reacting inviscid flows.
 *  \author G. Orlando.
 */
class CReactiveEulerSolver:public CSolver {
public:
  typedef std::unique_ptr<su2double[]> SmartArr;

  using   RealVec = CReactiveEulerVariable::RealVec;
  using   RealMatrix = CReactiveEulerVariable::RealMatrix;
  using   LibraryPtr = CReactiveEulerVariable::LibraryPtr;

protected:
  LibraryPtr library; /*!< \brief Smart pointer to the library that computes physical-chemical properties. */

  bool  space_centered,  /*!< \brief True if space centered scheeme used. */
        implicit,      /*!< \brief True if euler implicit scheme used. */
        least_squares;        /*!< \brief True if computing gradients by least squares. */

  const unsigned short nSpecies; /*!< \brief Total number of species. */

  unsigned long nMarker;         /*!< \brief Total number of markers using the grid information. */
  std::unique_ptr<unsigned long[]> nVertex;         /*!< \brief Total number of vertices using the grid information. */

  su2double  Gamma,              /*!< \brief Mixture Cp/Cv. */
	           Gamma_Minus_One;	   /*!< \brief Mixture Cp/Cv - 1. */

  su2double  Mach_Inf,       	  /*!< \brief Free stream Mach number. */
             Density_Inf,       /*!< \brief Free stream density. */
             Pressure_Inf,		  /*!< \brief Free stream pressure. */
	           Temperature_Inf;   /*!< \brief Trans.-rot. free stream temperature. */

  SmartArr   Density,       /*!< \brief Free stream species density. */
             Velocity_Inf,  /*!< \brief Free stream flow velocity. */
             MassFrac_Inf;  /*!< \brief Free stream species mass fraction. */

  SmartArr   Sol_i,  /*!< \brief Auxiliary vector for storing the solution at point i. */
             Sol_j,      /*!< \brief Auxiliary vector for storing the solution at point j. */
             Primitive,    /*!< \brief Auxiliary nPrimVar vector. */
             Primitive_i,        /*!< \brief Auxiliary nPrimVar vector for storing the primitive at point i. */
             Primitive_j;        /*!< \brief Auxiliary nPrimVar vector for storing the primitive at point j. */

public:

  /*!
	 * \brief Default constructor of the class.
	 */
   CReactiveEulerSolver():CSolver(),nSpecies(),nMarker(),space_centered(),implicit(),least_squares(),Gamma(),Gamma_Minus_One(),
                          Mach_Inf(),Density_Inf(),Pressure_Inf(),Temperature_Inf() {}

	/*!
	 * \overloaded constructor of the class
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	 CReactiveEulerSolver(std::unique_ptr<CGeometry> geometry, std::unique_ptr<CConfig> config,unsigned short iMesh);

	/*!
	 * \brief Destructor of the class.
	 */
	 virtual ~CReactiveEulerSolver() {}

  /*!
   * \brief Set gradient primitive variables using Green Gauss.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
   void SetPrimitive_Gradient_GG(CGeometry* geometry, CConfig* config) override;

  /*!
   * \brief Set gradient primitive variables using weighted least squares.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
   void SetPrimitive_Gradient_LS(CGeometry* geometry, CConfig* config) override;

  /*!
   * \brief Impose the send-receive boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
   void Set_MPI_Solution(CGeometry* geometry, CConfig* config) override;

  /*!
   * \brief Impose the send-receive boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
   void Set_MPI_Solution_Old(CGeometry* geometry, CConfig* config) override;

  /*!
   * \brief Impose the send-receive boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
   void Set_MPI_Primitive(CGeometry* geometry, CConfig* config) override;

  /*!
   * \brief Impose the send-receive boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
   void Set_MPI_Solution_Gradient(CGeometry* geometry, CConfig* config) override;

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
   virtual void SetPrimitive_Gradient(std::unique_ptr<CConfig> config);


  /*!
   * \brief Impose the send-receive boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
   void Set_MPI_Primitive_Gradient(CGeometry* geometry, CConfig* config) override;

  /*!
   * \brief Compute the limiter of the primitive variables.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
   void SetPrimitive_Limiter(CGeometry* geometry, CConfig* config) override;

  /*!
   * \brief Impose the send-receive boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
   void Set_MPI_Primitive_Limiter(CGeometry* geometry, CConfig* config) override;

  /*!
   * \brief Impose the send-receive boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
   void Set_MPI_Solution_Limiter(CGeometry* geometry, CConfig* config) override;

  /*!
   * \brief Set the maximum value of the eigenvalue.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
   void SetMax_Eigenvalue(CGeometry* geometry, CConfig* config) override;

  /*!
   * \brief Impose the send-receive boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
   void Set_MPI_MaxEigenvalue(CGeometry* geometry, CConfig* config) override;

  /*!
   * \brief Impose the send-receive boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
   void Set_MPI_Undivided_Laplacian(CGeometry* geometry, CConfig* config) override;

  /*!
   * \brief Set the fluid solver nondimensionalization.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
   void SetNondimensionalization(CGeometry* geometry, CConfig* config, unsigned short iMesh) override;

  /*!
   * \brief Preprocessing.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] Output - boolean to determine whether to print output.
   */
   void Preprocessing(CGeometry* geometry, CSolver** solver_container, CConfig* config,
                     unsigned short iMesh, unsigned short iRKStep,
                     unsigned short RunTime_EqSystem, bool Output) override;

  /*!
   * \brief Compute the time step for solving the Euler equations.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] Iteration - Index of the current iteration.
   */
   void SetTime_Step(CGeometry* geometry, CSolver** solver_container, CConfig* config,
                    unsigned short iMesh, unsigned long Iteration) override;

  /*!
   * \brief Compute the spatial integration using a centered scheme.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   */
   void Centered_Residual(CGeometry* geometry, CSolver** solver_container, CNumerics* numerics,
                         CConfig* config, unsigned short iMesh, unsigned short iRKStep) override;

  /*!
   * \brief Compute the spatial integration using a upwind scheme.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
   void Upwind_Residual(CGeometry* geometry, CSolver** solver_container, CNumerics* numerics,
                       CConfig* config, unsigned short iMesh) override;

   /*!
    * \brief Source term integration.
    * \param[in] geometry - Geometrical definition of the problem.
    * \param[in] solver_container - Container vector with all the solutions.
    * \param[in] numerics - Description of the numerical method.
    * \param[in] second_numerics - Description of the second numerical method.
    * \param[in] config - Definition of the particular problem.
    * \param[in] iMesh - Index of the mesh in multigrid computations.
    */
    void Source_Residual(CGeometry* geometry, CSolver** solver_container, CNumerics* numerics, CNumerics* second_numerics,
                       CConfig* config, unsigned short iMesh) override;

   /*!
    * \brief Impose via the residual the Euler wall boundary condition.
    * \param[in] geometry - Geometrical definition of the problem.
    * \param[in] solver_container - Container vector with all the solutions.
    * \param[in] numerics - Description of the numerical method.
    * \param[in] config - Definition of the particular problem.
    * \param[in] val_marker - Surface marker where the boundary condition is applied.
    */
    void BC_Euler_Wall(CGeometry* geometry, CSolver** solver_container, CNumerics* numerics, CConfig* config,
                       unsigned short val_marker) override;

   /*!
    * \brief Impose the far-field boundary condition.
    * \param[in] geometry - Geometrical definition of the problem.
    * \param[in] solver_container - Container vector with all the solutions.
    * \param[in] conv_numerics - Description of the numerical method for convective term.
    * \param[in] visc_numerics - Description of the numerical method for viscous term.
    * \param[in] config - Definition of the particular problem.
    * \param[in] val_marker - Surface marker where the boundary condition is applied.
    */
    void BC_Far_Field(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics, CNumerics* visc_numerics,
                    CConfig* config, unsigned short val_marker) override;

   /*!
    * \brief Impose the inlet boundary condition.
    * \param[in] geometry - Geometrical definition of the problem.
    * \param[in] solver_container - Container vector with all the solutions.
    * \param[in] conv_numerics - Description of the numerical method.
    * \param[in] visc_numerics - Description of the numerical method.
    * \param[in] config - Definition of the particular problem.
    * \param[in] val_marker - Surface marker where the boundary condition is applied.
    */
    void BC_Inlet(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics, CNumerics* visc_numerics, CConfig* config,
                  unsigned short val_marker) override;

   /*!
    * \brief Impose the outlet boundary condition.
    * \param[in] geometry - Geometrical definition of the problem.
    * \param[in] solver_container - Container vector with all the solutions.
    * \param[in] conv_numerics - Description of the numerical method for convective term.
    * \param[in] visc_numerics - Description of the numerical method for viscous term.
    * \param[in] config - Definition of the particular problem.
    * \param[in] val_marker - Surface marker where the boundary condition is applied.
    */
    void BC_Outlet(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics, CNumerics* visc_numerics, CConfig* config,
                 unsigned short val_marker) override;

   /*!
    * \brief Impose the symmetry plane boundary condition.
    * \param[in] geometry - Geometrical definition of the problem.
    * \param[in] solver_container - Container vector with all the solutions.
    * \param[in] conv_numerics - Description of the numerical method for convective term.
    * \param[in] visc_numerics - Description of the numerical method for viscous term.
    * \param[in] config - Definition of the particular problem.
    * \param[in] val_marker - Surface marker where the boundary condition is applied.
    */
    void BC_Sym_Plane(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                    CNumerics* visc_numerics, CConfig* config, unsigned short val_marker) override;

   /*!
    * \brief Update the solution using a Runge-Kutta scheme.
    * \param[in] geometry - Geometrical definition of the problem.
    * \param[in] solver_container - Container vector with all the solutions.
    * \param[in] config - Definition of the particular problem.
    * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
    */
    void ExplicitRK_Iteration(CGeometry* geometry, CSolver** solver_container, CConfig* config,
                            unsigned short iRKStep) override;

   /*!
    * \brief Update the solution using the explicit Euler scheme.
    * \param[in] geometry - Geometrical definition of the problem.
    * \param[in] solver_container - Container vector with all the solutions.
    * \param[in] config - Definition of the particular problem.
    */
    void ExplicitEuler_Iteration(CGeometry* geometry, CSolver** solver_container, CConfig* config) override;

   /*!
    * \brief Update the solution using an implicit Euler scheme.
    * \param[in] geometry - Geometrical definition of the problem.
    * \param[in] solver_container - Container vector with all the solutions.
    * \param[in] config - Definition of the particular problem.
    */
    void ImplicitEuler_Iteration(CGeometry* geometry, CSolver** solver_container, CConfig* config) override;

};


/*! \class CReactiveNSSolver
 *  \brief Main class for defining a solver for chemically reacting viscous flows.
 *  \author G. Orlando.
 */
class CReactiveNSSolver:public CReactiveEulerSolver {
public:

  /*!
	 * \brief Default constructor of the class.
	 */
   CReactiveNSSolver(void):CReactiveEulerSolver() {}

	/*!
	 * \overloaded constructor of the class
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CReactiveNSSolver(CGeometry* geometry, CConfig* config);

	/*!
	 * \brief Destructor of the class.
	 */
	~CReactiveNSSolver(void) {}

  /*!
   * \brief Preprocessing.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] Output - boolean to determine whether to print output.
   */
   void Preprocessing(CGeometry* geometry, CSolver** solver_container, CConfig* config,
                        unsigned short iMesh, unsigned short iRKStep,
                        unsigned short RunTime_EqSystem, bool Output) override;

  /*!
  * \brief Set primitive variables.
  * \param[in] solver_container - Container vector with all the solutions.
  * \param[in] config - Definition of the particular problem.
  * \param[in] Output - boolean to determine whether to print output.
  * \return - The number of non-physical points.
  */
  unsigned long SetPrimitive_Variables(CSolver** solver_container, CConfig* config, bool Output) override;

  /*!
   * \brief Compute the time step for solving the Navier-Stokes equations.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] Iteration - Index of the current iteration.
   */
   void SetTime_Step(CGeometry* geometry, CSolver** solver_container, CConfig* config,
                     unsigned short iMesh, unsigned long Iteration) override;

  /*!
   * \brief Impose the Navier-Stokes boundary condition (strong).
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method for convective term.
   * \param[in] visc_numerics - Description of the numerical method for viscous term.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
  */
  void BC_HeatFlux_Wall(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                        CNumerics* visc_numerics, CConfig* config, unsigned short val_marker) override;

  /*!
   * \brief Impose an isothermal wall boundary condition (no-slip).
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Isothermal_Wall(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                          CNumerics* visc_numerics, CConfig* config,unsigned short val_marker) override;

  /*!
   * \brief Compute the viscous residuals.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   */
  void Viscous_Residual(CGeometry* geometry, CSolver** solver_container, CNumerics* numerics,
                        CConfig* config, unsigned short iMesh, unsigned short iRKStep) override;

};

#endif
