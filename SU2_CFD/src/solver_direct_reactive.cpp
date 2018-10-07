#include "../include/solver_reactive.hpp"
#include "../../Common/include/reacting_model_library.hpp"
#include "../../Common/include/not_implemented_exception.hpp"
#include "../../Common/include/move_pointer.hpp"

namespace {
  using SmartArr = CReactiveEulerSolver::SmartArr;
  /*!
   * \brief Compute area for the current normal
   */
  su2double ComputeArea(const SmartArr& Normal,const unsigned short nDim) {
    su2double Area = 0.0;
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      Area += Normal[iDim]*Normal[iDim];
    Area = sqrt(Area);
    return Area;
  }
} /*-- End of unnamed namespace ---*/

//
//
/*--- Class constructor ---*/
//
//
CReactiveEulerSolver::CReactiveEulerSolver(std::unique_ptr<CGeometry> geometry, std::unique_ptr<CConfig> config,unsigned short iMesh):
                      CSolver(),library(new Framework::ReactingModelLibrary(config->GetLibraryName())),nSpecies(library->GetNSpecies()) {

  unsigned long iPoint, iVertex;
  unsigned short iVar, iDim, iMarker;

  nMarker  = config->GetnMarker_All();

  nVar = geometry->GetnDim();
  nPrimVar = nSpecies + nDim + 2;
  nPrimVarGrad = nSpecies + nDim + 5;

  /*--- Perform the non-dimensionalization for the flow equations using the
   specified reference values. ---*/
  SetNondimensionalization(geometry.get(), config.get(), iMesh);

	/*--- Store the number of vertices on each marker for deallocation ---*/
	nVertex = std::unique_ptr<unsigned long[]>(new unsigned long[nMarker]);
	for (iMarker = 0; iMarker < nMarker; iMarker++)
		nVertex[iMarker] = geometry->nVertex[iMarker];

	/*--- Allocate a CVariable array for each node of the mesh ---*/
	node = new CVariable*[nPoint];

	/*--- Define some auxiliary vectors related to the residual ---*/
	Residual = new su2double[nVar];
  for (iVar = 0; iVar < nVar; ++iVar)
    Residual[iVar] = 0.0;

  Residual_RMS = new su2double[nVar];
  for (iVar = 0; iVar < nVar; ++iVar)
    Residual_RMS[iVar] = 0.0;

  Residual_Max = new su2double[nVar];
  for (iVar = 0; iVar < nVar; ++iVar)
    Residual_Max[iVar] = 0.0;

  Residual_i = new su2double[nVar];
  for (iVar = 0; iVar < nVar; ++iVar)
    Residual_i[iVar] = 0.0;

  Residual_j = new su2double[nVar];
  for (iVar = 0; iVar < nVar; ++iVar)
    Residual_j[iVar]   = 0.0;

  Res_Conv = new su2double[nVar];
  for (iVar = 0; iVar < nVar; ++iVar)
    Res_Conv[iVar] = 0.0;

  Res_Visc = new su2double[nVar];
  for (iVar = 0; iVar < nVar; ++iVar)
    Res_Visc[iVar] = 0.0;

  Res_Sour = new su2double[nVar];
  for (iVar = 0; iVar < nVar; ++iVar)
    Res_Sour[iVar] = 0.0;

	/*--- Define some structure for locating max residuals ---*/
	Point_Max = new unsigned long[nVar];
  for (iVar = 0; iVar < nVar; ++iVar)
    Point_Max[iVar] = 0;

  Point_Max_Coord = new su2double*[nVar];
	for (iVar = 0; iVar < nVar; ++iVar){
		Point_Max_Coord[iVar] = new su2double[nDim];
		for (iDim = 0; iDim < nDim; ++iDim)
      Point_Max_Coord[iVar][iDim] = 0.0;
	}

  /*--- Allocate vectors related to the solution ---*/
  Sol_i = SmartArr(new su2double[nVar]);
  Sol_j = SmartArr(new su2double[nVar]);
  Primitive = SmartArr(new su2double[nPrimVar]);
  Primitive_i = SmartArr(new su2double[nPrimVar]);
  Primitive_j = SmartArr(new su2double[nPrimVar]);

  /*--- Allocate arrays for conserved variable limits ---*/
  /*lowerlimit = new su2double[nVar];
  upperlimit = new su2double[nVar];
  using RHO_INDEX_SOL = CReactiveEulerVariable::RHO_INDEX_SOL
  using VX_INDEX_SOL = CReactiveEulerVariable::VX_INDEX_SOL;
  using RHOE_INDEX_SOL = CReactiveEulerVariable::RHOE_INDEX_SOL;
  lowerlimit[RHO_INDEX_SOL] = 0.0;
  upperlimit[RHO_INDEX_SOL] = 1E16;
  for (iVar = VX_INDEX_SOL; iSpecies < VX_INDEX_SOL + nDim; ++iSpecies) {
    lowerlimit[iSpecies] = -1E16;
    upperlimit[iSpecies] = 1E16;
  }
  for (iVar = RHOE_INDEX_SOL; iVar < nVar; ++iVar) {
    lowerlimit[iVar] = 0;
    upperlimit[iVar] = 1E16;
  }*/


  /*--- Initialize the solution & residual CVectors ---*/
 	LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);

  /*--- Create the structure for storing extra information ---*/
 	if (config->GetExtraOutput()) {
    nOutputVariables = nVar;
    OutputVariables.Initialize(nPoint, nPointDomain, nOutputVariables, 0.0);
  }

	/*--- Allocate Jacobians for implicit time-stepping ---*/
	if (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT) {
    implicit = true;
		Jacobian_i = new su2double* [nVar];
		Jacobian_j = new su2double* [nVar];
		for (iVar = 0; iVar < nVar; iVar++) {
			Jacobian_i[iVar] = new su2double [nVar];
			Jacobian_j[iVar] = new su2double [nVar];
		}
  }
  else
    implicit = false;

  if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED )
    space_centered = true;
  else
    space_centered = false;

  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES)
    least_squares = true;
  else
    least_squares = false;

  Gamma = library->Gamma();
  Gamma_Minus_One = Gamma - 1.0;

  Density_Inf        = config->GetDensity_FreeStreamND();
  Pressure_Inf       = config->GetPressure_FreeStreamND();
	Temperature_Inf    = config->GetTemperature_FreeStreamND();
  Mach_Inf           = config->GetMach();

  Velocity_Inf       = Common::wrap_in_unique(config->GetVelocity_FreeStreamND());
  //MassFrac_Inf       = library->GetMassFrac_FreeStream();

  //node_infty = new CReactiveEulerVariable(Pressure_Inf, MassFrac_Inf,Velocity_Inf, Temperature_Inf,
  //                                       nDim, nVar,nPrimVar, nPrimVarGrad, config);

  /*--- Initialize the solution to the far-field state everywhere. ---*/
  //for (iPoint = 0; iPoint < nPoint; iPoint++)
    //node[iPoint] = new CReactiveEulerVariable(Pressure_Inf, MassFrac_Inf, Velocity_Inf, Temperature_Inf,
    //                                          nDim, nVar, nPrimVar, nPrimVarGrad,config);

  /*--- USE ANOTHER FUNCTION TO CHECK THAT THE INITIAL SOLUTION IS PHYSICAL ---*/

  /*--- MPI solution ---*/
	Set_MPI_Solution(geometry.get(), config.get());

}

//
//
/*!
 * \brief Set the fluid solver nondimensionalization.
 */
 void CReactiveEulerSolver::SetNondimensionalization(CGeometry* geometry, CConfig* config, unsigned short iMesh) {

 }
//
//



//
//
/*--- Variables preprocessing ---*/
//
//
void CReactiveEulerSolver::Preprocessing(CGeometry* geometry, CSolver** solver_container, CConfig* config,
                                         unsigned short iMesh, unsigned short iRKStep,
                                         unsigned short RunTime_EqSystem, bool Output) {


unsigned long ErrorCounter = 0;

unsigned long ExtIter = config->GetExtIter();
bool implicit         = config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT;
//bool low_fidelity     = (config->GetLowFidelitySim() && (iMesh == MESH_1));
bool second_order     = (config->GetSpatialOrder_Flow() == SECOND_ORDER or config->GetSpatialOrder_Flow() == SECOND_ORDER_LIMITER);
bool limiter          = (config->GetSpatialOrder_Flow() == SECOND_ORDER_LIMITER and ExtIter <= config->GetLimiterIter());
                        //and !low_fidelity;
bool center           = config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED;

/*--- Set the primitive variables ---*/
ErrorCounter = SetPrimitive_Variables(solver_container, config, Output);

/*--- Upwind second order reconstruction ---*/

if ((second_order and !center) and iMesh == MESH_0 and !Output) {

  /*--- Gradient computation ---*/

  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) {
    SetPrimitive_Gradient_GG(geometry, config);
  }
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    SetPrimitive_Gradient_LS(geometry, config);
  }

  /*--- Limiter computation ---*/

  if (limiter and iMesh == MESH_0 and !Output)
    SetPrimitive_Limiter(geometry, config);
}

/*--- Artificial dissipation ---*/

if (center && !Output)
    throw Common::NotImplemented("Centered convective scheme not implemented\n");


/*--- Initialize the Jacobian matrices ---*/

if (implicit)
  Jacobian.SetValZero();

/*--- Error message ---*/
if (config->GetConsole_Output_Verb() == VERB_HIGH) {
  #ifdef HAVE_MPI
    unsigned long MyErrorCounter = ErrorCounter; ErrorCounter = 0;
    SU2_MPI::Allreduce(&MyErrorCounter, &ErrorCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  #endif
  if (iMesh == MESH_0)
    config->SetNonphysical_Points(ErrorCounter);
}

}

//
//
/*--- Set primitive variables ---*/
//
//

unsigned long CReactiveEulerSolver::SetPrimitive_Variables(CSolver** solver_container, CConfig* config, bool Output) {

  unsigned long iPoint, ErrorCounter = 0;
  bool RightSol = true;

  for (iPoint = 0; iPoint < nPoint; iPoint ++) {

    /*--- Initialize the non-physical points vector ---*/

    node[iPoint]->SetNon_Physical(false);

    /*--- Compressible flow, primitive variables nSpecies+nDim+5, (T,vx,vy,vz,P,rho,h,c,rho1,...rhoNs) ---*/

    RightSol = static_cast<CReactiveEulerVariable*>(node[iPoint])->SetPrimVar();

    if (!RightSol) {
      node[iPoint]->SetNon_Physical(true);
      ErrorCounter++;
    }

    /*--- Initialize the convective, source and viscous residual vector ---*/

    if (!Output)
      LinSysRes.SetBlock_Zero(iPoint);

  }

  return ErrorCounter;

}

//
//
/*--- Setting time step ---*/
//
//

void CReactiveEulerSolver::SetTime_Step(CGeometry* geometry, CSolver** solver_container, CConfig* config,
                                        unsigned short iMesh, unsigned long Iteration) {

su2double Area, Vol, Mean_SoundSpeed = 0.0, Mean_ProjVel = 0.0, Lambda;
su2double Local_Delta_Time, Global_Delta_Time = 1E6;
SmartArr Normal;
unsigned long iEdge, iVertex, iPoint, jPoint;
unsigned short iDim, iMarker;

bool time_steping = config->GetUnsteady_Simulation() == TIME_STEPPING;

Min_Delta_Time = 1.E6;
Max_Delta_Time = 0.0;

/*--- Set maximum inviscid eigenvalue to zero, and compute sound speed ---*/

for (iPoint = 0; iPoint < nPointDomain; iPoint++)
  node[iPoint]->SetMax_Lambda_Inv(0.0);


/*--- Loop interior edges ---*/
for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

  /*--- Point identification, Normal vector and area ---*/

  iPoint = geometry->edge[iEdge]->GetNode(0);
  jPoint = geometry->edge[iEdge]->GetNode(1);

  Normal = Common::wrap_in_unique(geometry->edge[iEdge]->GetNormal());
  Area = ::ComputeArea(Normal,nDim);

  /*--- Mean Values ---*/

  Mean_ProjVel = 0.5 * (node[iPoint]->GetProjVel(Normal.get()) + node[jPoint]->GetProjVel(Normal.get()));
  Mean_SoundSpeed = 0.5 * (node[iPoint]->GetSoundSpeed() + node[jPoint]->GetSoundSpeed()) * Area;

  /*--- Inviscid contribution ---*/
  Lambda = std::abs(Mean_ProjVel) + Mean_SoundSpeed;

  if (geometry->node[iPoint]->GetDomain())
    node[iPoint]->AddMax_Lambda_Inv(Lambda);
  if (geometry->node[jPoint]->GetDomain())
    node[jPoint]->AddMax_Lambda_Inv(Lambda);

}

/*--- Loop boundary edges ---*/

for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
  if (config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY) {
    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {

      /*--- Point identification, Normal vector and area ---*/

      iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
      Normal = Common::wrap_in_unique(geometry->vertex[iMarker][iVertex]->GetNormal());
      Area = ::ComputeArea(Normal,nDim);

      /*--- Mean Values ---*/

      Mean_ProjVel = node[iPoint]->GetProjVel(Normal.get());
      Mean_SoundSpeed = node[iPoint]->GetSoundSpeed() * Area;

      /*--- Inviscid contribution ---*/
      Lambda = std::abs(Mean_ProjVel) + Mean_SoundSpeed;
      if (geometry->node[iPoint]->GetDomain())
        node[iPoint]->AddMax_Lambda_Inv(Lambda);
    }
  }
}

/*--- Each element uses their own speed for a steady state simulation ---*/

for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

  Vol = geometry->node[iPoint]->GetVolume();

  if (Vol != 0.0) {
    Local_Delta_Time = config->GetCFL(iMesh)*Vol / node[iPoint]->GetMax_Lambda_Inv();
    Global_Delta_Time = std::min(Global_Delta_Time, Local_Delta_Time);
    Min_Delta_Time = std::min(Min_Delta_Time, Local_Delta_Time);
    Max_Delta_Time = std::max(Max_Delta_Time, Local_Delta_Time);
    if (Local_Delta_Time > config->GetMax_DeltaTime())
      Local_Delta_Time = config->GetMax_DeltaTime();
      node[iPoint]->SetDelta_Time(Local_Delta_Time);
  }
  else {
    node[iPoint]->SetDelta_Time(0.0);
  }
}

/*--- Compute the max and the min dt (in parallel) ---*/
if (config->GetConsole_Output_Verb() == VERB_HIGH) {
  #ifdef HAVE_MPI
    su2double rbuf_time, sbuf_time;
    sbuf_time = Min_Delta_Time;
    SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    Min_Delta_Time = rbuf_time;

    sbuf_time = Max_Delta_Time;
    SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MAX, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    Max_Delta_Time = rbuf_time;
  #endif
}

if(time_steping) {
  #ifdef HAVE_MPI
    su2double rbuf_time, sbuf_time;
    sbuf_time = Global_Delta_Time;
    SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    Global_Delta_Time = rbuf_time;
  #endif
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    /*--- Sets the regular CFL equal to the unsteady CFL ---*/
    config->SetCFL(iMesh,config->GetUnst_CFL());

    /*--- If the unsteady CFL is set to zero, it uses the defined unsteady time step, otherwise
    it computes the time step based on the unsteady CFL ---*/
    if (config->GetCFL(iMesh) == 0.0)
      node[iPoint]->SetDelta_Time(config->GetDelta_UnstTime());
    else
      node[iPoint]->SetDelta_Time(Global_Delta_Time);

  }

}

}

//
//
/*--- Iteration of implicit Euler method ---*/
//
//

void CReactiveEulerSolver::ImplicitEuler_Iteration(CGeometry* geometry, CSolver** solver_container, CConfig* config) {

  unsigned short iVar;
	unsigned long iPoint, total_index, IterLinSol = 0;
	SmartArr local_Res_TruncError;

  /*--- Set maximum residual to zero ---*/
    for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }

  /*--- Build implicit system ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    /*--- Read the residual ---*/
    local_Res_TruncError = Common::wrap_in_unique(node[iPoint]->GetResTruncError());
    /*--- Read the volume ---*/
    auto Vol = geometry->node[iPoint]->GetVolume();

    if (node[iPoint]->GetDelta_Time() != 0.0) {
      auto Delta = Vol / node[iPoint]->GetDelta_Time();
      Jacobian.AddVal2Diag(iPoint, Delta);
    }
    else {
      Jacobian.SetVal2Diag(iPoint, 1.0);
      for (iVar = 0; iVar < nVar; iVar++) {
        total_index = iPoint*nVar + iVar;
        LinSysRes[total_index] = 0.0;
        local_Res_TruncError[iVar] = 0.0;
      }
    }

    /*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/
    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar + iVar;
      LinSysRes[total_index] = - (LinSysRes[total_index] + local_Res_TruncError[iVar]);
      LinSysSol[total_index] = 0.0;
      AddRes_RMS(iVar, LinSysRes[total_index]*LinSysRes[total_index]);
      AddRes_Max(iVar, abs(LinSysRes[total_index]), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
    }
  }

  /*--- Initialize residual and solution at the ghost points ---*/
  for (iPoint = nPointDomain; iPoint < nPoint; iPoint++) {
    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar + iVar;
      LinSysRes[total_index] = 0.0;
      LinSysSol[total_index] = 0.0;
    }
  }

  /*--- Solve or smooth the linear system ---*/
  CSysSolve system;
  IterLinSol = system.Solve(Jacobian, LinSysRes, LinSysSol, geometry, config);

  /*--- The the number of iterations of the linear solver ---*/
  SetIterLinSolver(IterLinSol);

  /*--- Update solution (system written in terms of increments) ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++)
		for (iVar = 0; iVar < nVar; iVar++)
			node[iPoint]->AddSolution(iVar, LinSysSol[iPoint*nVar+iVar]);


  /*--- Compute the root mean square residual ---*/

  SetResidual_RMS(geometry, config);

}

//
//
/*--- Iteration of explicit Euler method ---*/
//
//

void CReactiveEulerSolver::ExplicitEuler_Iteration(CGeometry* geometry, CSolver** solver_container, CConfig* config) {

  su2double Res;
  SmartArr Residual, Res_TruncError;
  unsigned short iVar;
  unsigned long iPoint;

  for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }

  /*--- Update the solution ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    auto Vol = geometry->node[iPoint]->GetVolume();
    auto Delta = node[iPoint]->GetDelta_Time() / Vol;

    Res_TruncError = Common::wrap_in_unique(node[iPoint]->GetResTruncError());
    Residual = Common::wrap_in_unique(LinSysRes.GetBlock(iPoint));

    for (iVar = 0; iVar < nVar; iVar++) {
        Res = Residual[iVar] + Res_TruncError[iVar];
        node[iPoint]->AddSolution(iVar, -Res*Delta);
        AddRes_RMS(iVar, Res*Res);
        AddRes_Max(iVar, abs(Res), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
      }
    }

  /*--- Compute the root mean square residual ---*/

  SetResidual_RMS(geometry, config);

}

//
//
/*--- Iteration of explicit RK method ---*/
//
//

void CReactiveEulerSolver::ExplicitRK_Iteration(CGeometry* geometry, CSolver** solver_container,
                                        CConfig* config, unsigned short iRKStep) {
  su2double Res;
  SmartArr Residual, Res_TruncError;
  unsigned short iVar;
  unsigned long iPoint;

  auto RK_AlphaCoeff = config->Get_Alpha_RKStep(iRKStep);

  for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }

  /*--- Update the solution ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    auto Vol = geometry->node[iPoint]->GetVolume();
    auto Delta = node[iPoint]->GetDelta_Time() / Vol;

    Res_TruncError = Common::wrap_in_unique(node[iPoint]->GetResTruncError());
    Residual = Common::wrap_in_unique(LinSysRes.GetBlock(iPoint));

    for (iVar = 0; iVar < nVar; iVar++) {
      Res = Residual[iVar] + Res_TruncError[iVar];
      node[iPoint]->AddSolution(iVar, -Res*Delta*RK_AlphaCoeff);
      AddRes_RMS(iVar, Res*Res);
      AddRes_Max(iVar, abs(Res), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
    }
  }

  /*--- Compute the root mean square residual ---*/

  SetResidual_RMS(geometry, config);

}

//
//
/*--- Centered residual convective term ---*/
//
//

void CReactiveEulerSolver::Centered_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                             CConfig *config, unsigned short iMesh, unsigned short iRKStep) {

    throw Common::NotImplemented("Function not implemented: Centered Residual\n");
}

//
//
/*--- Upwind residual convective term ---*/
//
//

void CReactiveEulerSolver::Upwind_Residual(CGeometry* geometry, CSolver** solver_container, CNumerics* numerics,
                                           CConfig* config, unsigned short iMesh) {

    unsigned long iEdge, iPoint, jPoint;
    unsigned short iDim, iVar;

    SmartArr Limiter_i,Limiter_j;
    su2double Project_Grad_i, Project_Grad_j;

    bool implicit     = config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT;
    //bool low_fidelity = (config->GetLowFidelitySim() && (iMesh == MESH_1));
    bool second_order = ((config->GetSpatialOrder_Flow() == SECOND_ORDER or config->GetSpatialOrder_Flow() == SECOND_ORDER_LIMITER)
                        and iMesh == MESH_0);
    bool limiter      = config->GetSpatialOrder_Flow() == SECOND_ORDER_LIMITER;

    /*--- Loop over all the edges ---*/

    for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

      /*--- Points in edge and normal vectors ---*/

      iPoint = geometry->edge[iEdge]->GetNode(0);
      jPoint = geometry->edge[iEdge]->GetNode(1);
      numerics->SetNormal(geometry->edge[iEdge]->GetNormal());

      /*--- Get primitive variables ---*/

      Primitive_i = Common::wrap_in_unique(node[iPoint]->GetPrimitive());
      Primitive_j = Common::wrap_in_unique(node[jPoint]->GetPrimitive());

      if (second_order) {

        for (iDim = 0; iDim < nDim; iDim++) {
          Vector_i[iDim] = 0.5*(geometry->node[jPoint]->GetCoord(iDim) - geometry->node[iPoint]->GetCoord(iDim));
          Vector_j[iDim] = 0.5*(geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim));
        }

      //  node[iPoint]->GetGradient_Primitive();
      //  node[jPoint]->GetGradient_Primitive();

        if (limiter) {
          Limiter_i = Common::wrap_in_unique(node[iPoint]->GetLimiter_Primitive());
          Limiter_j = Common::wrap_in_unique(node[jPoint]->GetLimiter_Primitive());
        }

        for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
          Project_Grad_i = 0.0; Project_Grad_j = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim];
            Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim];
          }
          if (limiter) {
            Primitive_i[iVar] += Limiter_i[iVar]*Project_Grad_i;
            Primitive_j[iVar] += Limiter_j[iVar]*Project_Grad_j;
          }
          else {
            Primitive_i[iVar] += Project_Grad_i;
            Primitive_j[iVar] += Project_Grad_j;
          }
        }
        numerics->SetPrimitive(Primitive_i.get(), Primitive_j.get());
      }
      else {

        /*--- Set primitive variables without reconstruction ---*/

        numerics->SetPrimitive(Primitive_i.get(), Primitive_j.get());

      }
      /*--- Set conserative variables ---*/

      numerics->SetConservative(node[iPoint]->GetSolution(),node[jPoint]->GetSolution());

      /*--- Compute the residual ---*/

      numerics->ComputeResidual(Res_Conv, Jacobian_i, Jacobian_j, config);

      /*--- Update residual value ---*/

      LinSysRes.AddBlock(iPoint, Res_Conv);
      LinSysRes.SubtractBlock(jPoint, Res_Conv);

      /*--- Set implicit Jacobians ---*/

      if (implicit) {
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
        Jacobian.AddBlock(iPoint, jPoint, Jacobian_j);
        Jacobian.SubtractBlock(jPoint, iPoint, Jacobian_i);
        Jacobian.SubtractBlock(jPoint, jPoint, Jacobian_j);
      }
    }
}

//
//
/*--- Residual source term ---*/
//
//

void CReactiveEulerSolver::Source_Residual(CGeometry* geometry, CSolver** solver_container,
                                           CNumerics* numerics, CNumerics* second_numerics,
                                           CConfig* config, unsigned short iMesh) {

unsigned long iPoint;
unsigned short iVar;

bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

/*--- Initialize the source residual to zero ---*/
for (iVar = 0; iVar < nVar; iVar++)
  Residual[iVar] = 0.0;

/*--- Loop over all points ---*/
for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

  /*--- Load the conservative variables ---*/
  numerics->SetConservative(node[iPoint]->GetSolution(),
                            node[iPoint]->GetSolution());

  /*--- Load the volume of the dual mesh cell ---*/
  numerics->SetVolume(geometry->node[iPoint]->GetVolume());

  /*--- Compute the rotating frame source residual ---*/
  numerics->ComputeResidual(Residual, Jacobian_i, config);

  /*--- Add the source residual to the total ---*/
  LinSysRes.AddBlock(iPoint, Residual);

  /*--- Add the implicit Jacobian contribution ---*/
  if (implicit)
    Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

}

}

//
//
/*--- Euler_wall_BC ---*/
//
//

void CReactiveEulerSolver::BC_Euler_Wall(CGeometry* geometry, CSolver** solver_container, CNumerics* solver,
                                         CConfig* config, unsigned short val_marker) {

  unsigned long iPoint, iVertex;
	unsigned short iDim, iVar;
	su2double Pressure;
  SmartArr Normal;

	bool implicit = config->GetKind_TimeIntScheme() == EULER_IMPLICIT;

	/*--- Buckle over all the vertices ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {

			Normal = Common::wrap_in_unique(geometry->vertex[val_marker][iVertex]->GetNormal());
      su2double Area = ::ComputeArea(Normal,nDim);

      RealVec UnitaryNormal(3);
			for (iDim = 0; iDim < nDim; iDim++)
				UnitaryNormal[iDim] = -Normal[iDim]/Area;

				//Pressure = library->node[iPoint]->GetPressure(iSpecies);

				//Residual[0] = 0.0;
				for (iDim = 0; iDim < nDim; iDim++)
					Residual[iDim+1] =  Pressure*UnitaryNormal[iDim]*Area;

        //Residual[nDim+1] = 0.0;
				//Residual[nDim+2] = 0.0;


			/*--- Add value to the residual ---*/
			LinSysRes.AddBlock(iPoint, Residual);

			/*--- In case we are doing a implicit computation ---*/

			if (implicit) {
				su2double a2;
				su2double phi;
				su2double Energy_el=0.0;

				//a2 = library->config->GetSpecies_Gamma(iSpecies)-1.0;
				//phi = a2*(library->0.5*node[iPoint]->GetVelocity2(iSpecies) - config->GetEnthalpy_Formation(iSpecies) - Energy_el);

				for (iVar = 0; iVar < nVar; iVar++) {
					Jacobian_i[0][iVar] = 0.0;
					Jacobian_i[nDim+1][iVar] = 0.0;
					Jacobian_i[nDim+2][iVar] = 0.0;
        }

				for (iDim = 0; iDim < nDim; iDim++) {
					Jacobian_i[iDim+1][0] = -phi*Normal[iDim];

          for (unsigned short jDim = 0; jDim < nDim; jDim++)
					//	Jacobian_i[iDim+1][jDim+1] = a2*node[iPoint]->GetVelocity(jDim,iSpecies)*Normal[iDim];

            Jacobian_i[iDim+1][nDim+1] = -a2*Normal[iDim];
						Jacobian_i[iDim+1][nDim+2] = a2*Normal[iDim];
				}

			Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
			}
		}
	}
}

//
//
/*--- Far-field_BC ---*/
//
//

void CReactiveEulerSolver::BC_Far_Field(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                                        CNumerics* visc_numerics, CConfig* config, unsigned short val_marker) {

  unsigned long iVertex, iPoint;
	unsigned short iVar, iDim;
	RealVec U_domain(nVar), U_infty(nVar);

	bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

	/*--- Bucle over all the vertices ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {

			/*--- Interpolated solution at the wall ---*/
			for (iVar = 0; iVar < nVar; iVar++)
				U_domain[iVar] = node[iPoint]->GetSolution(iVar);

			/*--- Solution at the infinity ---*/
			/*U_infty[0] = GetDensity_Inf(iSpecies);

			for (iDim = 0; iDim < nDim; iDim++)
				U_infty[iDim+1] = GetDensity_Velocity_Inf(iDim, iSpecies);

        U_infty[nDim+1] = GetDensity_Energy_Inf(iSpecies);*/

      /*--- Set the normal vector ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
			for (iDim = 0; iDim < nDim; iDim++)
        Vector[iDim] = -Vector[iDim];

      conv_numerics->SetNormal(Vector);

			/*--- Compute the residual using an upwind scheme ---*/
			conv_numerics->SetConservative(U_domain.data(), U_infty.data());
			conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

      LinSysRes.AddBlock(iPoint, Residual);

			/*--- In case we are doing a implicit computation ---*/
			if (implicit)
				Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
		}
	}
}

//
//
/*--- Symmetry Plane BC ---*/
//
//

void CReactiveEulerSolver::BC_Sym_Plane(CGeometry* geometry, CSolver** solver_container,
                                        CNumerics* conv_numerics, CNumerics* visc_numerics,
                                        CConfig* config, unsigned short val_marker) {

  unsigned long iPoint, iVertex;
	unsigned short iDim, iVar;
	su2double Pressure, Energy_el;
  SmartArr Normal;

	bool implicit = config->GetKind_TimeIntScheme() == EULER_IMPLICIT;

	/*--- Buckle over all the vertices ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {

			/*--- Compute the projected residual ---*/
			Normal = Common::wrap_in_unique(geometry->vertex[val_marker][iVertex]->GetNormal());
      su2double Area = ::ComputeArea(Normal,nDim);

      RealVec UnitaryNormal(3);

			for (iDim = 0; iDim < nDim; iDim++)
				UnitaryNormal[iDim] = -Normal[iDim]/Area;

			//Residual[0] = 0.0;
			//Residual[nDim+1] = 0.0;
			//Residual[nDim+2] = 0.0;

      //Pressure = node[iPoint]->GetPressure(iSpecies);
			for (iDim = 0; iDim < nDim; iDim++)
				Residual[iDim+1] = Pressure*UnitaryNormal[iDim]*Area;


      /*--- Add value to the residual ---*/
			LinSysRes.AddBlock(iPoint, Residual);

			/*--- In case we are doing a implicit computation ---*/
			if (implicit) {
				su2double a2 = 0.0;
				su2double phi = 0.0;
				Energy_el = 0.0;

        //a2 = config->GetSpecies_Gamma(iSpecies)-1.0;
				//phi = a2*(0.5*node[iPoint]->GetVelocity2(iSpecies) - config->GetEnthalpy_Formation(iSpecies) - Energy_el);

				for (iVar =  0; iVar <  nVar; iVar++) {
					Jacobian_i[0][iVar] = 0.0;
					Jacobian_i[nDim+1][iVar] = 0.0;
					Jacobian_i[nDim+2][iVar] = 0.0;
				}
				for (iDim = 0; iDim < nDim; iDim++) {
					Jacobian_i[iDim+1][0] = -phi*Normal[iDim];
					//for (unsigned short jDim = 0; jDim < nDim; jDim++)
					//	Jacobian_i[iDim+1][loc + jDim+1] = a2*node[iPoint]->GetVelocity(jDim,iSpecies)*Normal[iDim];

          Jacobian_i[iDim+1][nDim+1] = -a2*Normal[iDim];
					Jacobian_i[iDim+1][nDim+2] = a2*Normal[iDim];
				}

				Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
			}

		}
	}
}


//
//
/*--- Preprocessing variables for NS ---*/
//
//

void CReactiveNSSolver::Preprocessing(CGeometry* geometry, CSolver** solver_container, CConfig* config,
                                      unsigned short iMesh, unsigned short iRKStep,
                                      unsigned short RunTime_EqSystem, bool Output) {

unsigned long ErrorCounter = 0;

unsigned long ExtIter = config->GetExtIter();
bool implicit         = config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT;
bool second_order     = ((config->GetSpatialOrder_Flow() == SECOND_ORDER or config->GetSpatialOrder_Flow() == SECOND_ORDER_LIMITER)
                         and iMesh == MESH_0);
bool limiter          = (config->GetSpatialOrder_Flow() == SECOND_ORDER_LIMITER and ExtIter <= config->GetLimiterIter());
bool center           = (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED);
bool limiter_visc     = config->GetViscous_Limiter_Flow();

/*--- Set the primitive variables ---*/

ErrorCounter = SetPrimitive_Variables(solver_container, config, Output);


/*--- Gradient computation ---*/

if (config->GetKind_Gradient_Method() == GREEN_GAUSS) {
    SetPrimitive_Gradient_GG(geometry, config);
}
if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    SetPrimitive_Gradient_LS(geometry, config);
}


/*--- Limiter computation ---*/

if ((limiter || limiter_visc) && (iMesh == MESH_0) && !Output)
  SetPrimitive_Limiter(geometry, config);


/*--- Artificial dissipation ---*/

if (center && !Output)
    throw Common::NotImplemented("Centered convective scheme not implemented\n");


/*--- Initialize the Jacobian matrices ---*/

if (implicit && !config->GetDiscrete_Adjoint())
  Jacobian.SetValZero();

}

//
//
/*--- Setting primitive variables NS ---*/
//
//

unsigned long CReactiveNSSolver::SetPrimitive_Variables(CSolver** solver_container, CConfig* config, bool Output) {

  unsigned long iPoint, ErrorCounter = 0;
  bool RightSol = true;

  for (iPoint = 0; iPoint < nPoint; iPoint ++) {

    /*--- Initialize the non-physical points vector ---*/

    node[iPoint]->SetNon_Physical(false);

    /*--- Compressible flow, primitive variables nSpecies+nDim+8, (T, vx, vy, vz, P, rho, h, c, lamMu, ThCond, Cp,rhos) ---*/

    RightSol = static_cast<CReactiveNSVariable*>(node[iPoint])->SetPrimVar();

    if (!RightSol) {
       node[iPoint]->SetNon_Physical(true);
       ErrorCounter++;
    }

    /*--- Initialize the convective, source and viscous residual vector ---*/

    if (!Output) LinSysRes.SetBlock_Zero(iPoint);

  }

  return ErrorCounter;

}

//
//
/*--- Setting time step NS ---*/
//
//

void CReactiveNSSolver::SetTime_Step(CGeometry* geometry, CSolver** solver_container, CConfig* config,
                                        unsigned short iMesh, unsigned long Iteration) {

su2double Area, Vol, Mean_SoundSpeed = 0.0, Mean_ProjVel = 0.0, Lambda, Local_Delta_Time,
          Lambda_1,Global_Delta_Time = 1E6,Mean_LaminarVisc = 0.0,Mean_Density = 0.0;
SmartArr Normal;
unsigned long iEdge, iVertex, iPoint, jPoint;
unsigned short iDim, iMarker;

bool implicit = config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT;
bool time_steping = config->GetUnsteady_Simulation() == TIME_STEPPING;

Min_Delta_Time = 1.E6; Max_Delta_Time = 0.0;

/*--- Set maximum inviscid eigenvalue to zero, and compute sound speed ---*/

for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
  node[iPoint]->SetMax_Lambda_Inv(0.0);
  node[iPoint]->SetMax_Lambda_Visc(0.0);
}

/*--- Loop interior edges ---*/
for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

  /*--- Point identification, Normal vector and area ---*/

  iPoint = geometry->edge[iEdge]->GetNode(0);
  jPoint = geometry->edge[iEdge]->GetNode(1);

  Normal = Common::wrap_in_unique(geometry->edge[iEdge]->GetNormal());
  Area = ::ComputeArea(Normal,nDim);

  /*--- Mean Values ---*/

  Mean_ProjVel = 0.5 * (node[iPoint]->GetProjVel(Normal.get()) + node[jPoint]->GetProjVel(Normal.get()));
  Mean_SoundSpeed = 0.5 * (node[iPoint]->GetSoundSpeed() + node[jPoint]->GetSoundSpeed()) * Area;


  /*--- Inviscid contribution ---*/
  Lambda = std::abs(Mean_ProjVel) + Mean_SoundSpeed;

  if (geometry->node[iPoint]->GetDomain())
    node[iPoint]->AddMax_Lambda_Inv(Lambda);
  if (geometry->node[jPoint]->GetDomain())
    node[jPoint]->AddMax_Lambda_Inv(Lambda);

  /*--- Viscous contribution ---*/

  Mean_LaminarVisc = 0.5*(node[iPoint]->GetLaminarViscosity() + node[jPoint]->GetLaminarViscosity());
  Mean_Density     = 0.5*(node[iPoint]->GetSolution(0) + node[jPoint]->GetSolution(0));

  Lambda_1 = (4.0/3.0)*Mean_LaminarVisc;
  Lambda = Lambda_1*Area*Area/Mean_Density;

  if (geometry->node[iPoint]->GetDomain())
    node[iPoint]->AddMax_Lambda_Visc(Lambda);
  if (geometry->node[jPoint]->GetDomain())
    node[jPoint]->AddMax_Lambda_Visc(Lambda);


}

/*--- Loop boundary edges ---*/

for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
  if (config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY)
  for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {

    /*--- Point identification, Normal vector and area ---*/

    iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
    Normal = Common::wrap_in_unique(geometry->vertex[iMarker][iVertex]->GetNormal());
    Area = ::ComputeArea(Normal,nDim);

    /*--- Mean Values ---*/

    Mean_ProjVel = node[iPoint]->GetProjVel(Normal.get());
    Mean_SoundSpeed = node[iPoint]->GetSoundSpeed() * Area;

    /*--- Inviscid contribution ---*/
    Lambda = std::abs(Mean_ProjVel) + Mean_SoundSpeed;
    if (geometry->node[iPoint]->GetDomain())
      node[iPoint]->AddMax_Lambda_Inv(Lambda);

    /*--- Viscous contribution ---*/

    Mean_LaminarVisc = 0.5*(node[iPoint]->GetLaminarViscosity() + node[jPoint]->GetLaminarViscosity());
    Mean_Density     = 0.5*(node[iPoint]->GetSolution(0) + node[jPoint]->GetSolution(0));

    Lambda_1 = (4.0/3.0)*Mean_LaminarVisc;
    Lambda = Lambda_1 *Area*Area/Mean_Density;

    if (geometry->node[iPoint]->GetDomain())
      node[iPoint]->AddMax_Lambda_Visc(Lambda);

  }
}

/*--- Each element uses their own speed, steady state simulation ---*/

for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

  Vol = geometry->node[iPoint]->GetVolume();

  if (Vol != 0.0) {
    Local_Delta_Time = config->GetCFL(iMesh)*Vol / node[iPoint]->GetMax_Lambda_Inv();
    Global_Delta_Time = std::min(Global_Delta_Time, Local_Delta_Time);
    Min_Delta_Time = std::min(Min_Delta_Time, Local_Delta_Time);
    Max_Delta_Time = std::max(Max_Delta_Time, Local_Delta_Time);
    if (Local_Delta_Time > config->GetMax_DeltaTime())
      Local_Delta_Time = config->GetMax_DeltaTime();
      node[iPoint]->SetDelta_Time(Local_Delta_Time);
  }
  else {
    node[iPoint]->SetDelta_Time(0.0);
  }
}

for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    if(time_steping) {
    /*--- Sets the regular CFL equal to the unsteady CFL ---*/
    config->SetCFL(iMesh,config->GetUnst_CFL());

    /*--- If the unsteady CFL is set to zero, it uses the defined unsteady time step, otherwise
    it computes the time step based on the unsteady CFL ---*/
    if (config->GetCFL(iMesh) == 0.0)
      node[iPoint]->SetDelta_Time(config->GetDelta_UnstTime());
    else
      node[iPoint]->SetDelta_Time(Global_Delta_Time);
    }

  }
}

//
//
/*--- Viscous Residual ---*/
//
//

void CReactiveNSSolver::Viscous_Residual(CGeometry* geometry, CSolver** solver_container, CNumerics* numerics,
                                         CConfig* config, unsigned short iMesh, unsigned short iRKStep) {

unsigned long iPoint, jPoint, iEdge;

bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

  /*--- Points, coordinates and normal vector in edge ---*/

  iPoint = geometry->edge[iEdge]->GetNode(0);
  jPoint = geometry->edge[iEdge]->GetNode(1);
  numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[jPoint]->GetCoord());
  numerics->SetNormal(geometry->edge[iEdge]->GetNormal());

  /*--- Primitive and secondary variables ---*/

  numerics->SetPrimitive(node[iPoint]->GetPrimitive(), node[jPoint]->GetPrimitive());
  numerics->SetSecondary(node[iPoint]->GetSecondary(), node[jPoint]->GetSecondary());

  /*--- Gradient and limiters ---*/

  numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[jPoint]->GetGradient_Primitive());
  numerics->SetPrimVarLimiter(node[iPoint]->GetLimiter_Primitive(), node[jPoint]->GetLimiter_Primitive());

  /*--- Compute and update residual ---*/

  numerics->ComputeResidual(Res_Visc, Jacobian_i, Jacobian_j, config);

  LinSysRes.SubtractBlock(iPoint, Res_Visc);
  LinSysRes.AddBlock(jPoint, Res_Visc);

  /*--- Implicit part ---*/

  if (implicit) {
    Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
    Jacobian.SubtractBlock(iPoint, jPoint, Jacobian_j);
    Jacobian.AddBlock(jPoint, iPoint, Jacobian_i);
    Jacobian.AddBlock(jPoint, jPoint, Jacobian_j);
  }

}

}
