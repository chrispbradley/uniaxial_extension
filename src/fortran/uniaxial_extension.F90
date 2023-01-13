PROGRAM UniaxialExtension

  USE OpenCMISS
  USE OpenCMISS_Iron
#ifndef NOMPIMOD
  USE MPI
#endif

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif

  INTEGER(CMISSIntg), PARAMETER :: CONTEXT_USER_NUMBER = 1

  TYPE(cmfe_ContextType) :: context
  
  INTEGER(CMISSIntg) :: err

  !Intialise OpenCMISS
  CALL cmfe_Initialise(err)
  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,err)
  !Create a context
  CALL cmfe_Context_Initialise(context,err)
  CALL cmfe_Context_Create(CONTEXT_USER_NUMBER,context,err)
  
  !Set all diganostic levels on for testing
  !CALL cmfe_DiagnosticsSetOn(CMFE_FROM_DIAG_TYPE,[1,2,3,4,5],"Diagnostics", &
  !  & ["FiniteElasticity_FiniteElementResidualEvaluateNew", &
  !  &  "FiniteElasticity_FiniteElementJacobianEvaluateNew"],Err)

  !Input arguments: context, compressible, useGeneratedMesh, zeroLoad, useSimplex, usePressureBasis

  CALL SolveModel(context, 3, .FALSE., .FALSE., .FALSE., .FALSE., .FALSE.)
  !CALL SolveModel(context, 3, .TRUE., .FALSE. , .FALSE., .FALSE., .FALSE.)
  CALL SolveModel(context, 3, .FALSE., .TRUE. , .FALSE., .FALSE., .FALSE.)
  CALL SolveModel(context, 3, .FALSE., .FALSE., .TRUE., .FALSE., .FALSE.)
  CALL SolveModel(context, 3, .FALSE., .FALSE., .FALSE., .TRUE., .FALSE.)
  CALL SolveModel(context, 3, .FALSE., .TRUE., .FALSE., .FALSE., .TRUE.)

  !Destroy the context
  CALL cmfe_Context_Destroy(context,err)
  !Finalise OpenCMISS
  CALL cmfe_Finalise(err)

  STOP

CONTAINS

  SUBROUTINE SolveModel(context, numberOfDimensions, compressible, useGeneratedMesh, zeroLoad, useSimplex, usePressureBasis)

    TYPE(cmfe_ContextType), INTENT(IN) :: context
    INTEGER(CMISSIntg), INTENT(IN)  :: numberOfDimensions
    LOGICAL, INTENT(IN) :: compressible
    LOGICAL, INTENT(IN) :: useGeneratedMesh
    LOGICAL, INTENT(IN) :: zeroLoad
    LOGICAL, INTENT(IN) :: useSimplex
    LOGICAL, INTENT(IN) :: usePressureBasis

    REAL(CMISSRP), PARAMETER :: HEIGHT = 1.0_CMISSRP
    REAL(CMISSRP), PARAMETER :: WIDTH = 1.0_CMISSRP
    REAL(CMISSRP), PARAMETER :: LENGTH = 1.0_CMISSRP

    INTEGER(CMISSIntg), PARAMETER :: COORDINATE_SYSTEM_USER_NUMBER = 1
    INTEGER(CMISSIntg), PARAMETER :: REGION_USER_NUMBER = 1
    INTEGER(CMISSIntg), PARAMETER :: BASIS_USER_NUMBER = 1
    INTEGER(CMISSIntg), PARAMETER :: PRESSURE_BASIS_USER_NUMBER = 2
    INTEGER(CMISSIntg), PARAMETER :: GENERATED_MESH_USER_NUMBER = 1
    INTEGER(CMISSIntg), PARAMETER :: MESH_USER_NUMBER = 1
    INTEGER(CMISSIntg), PARAMETER :: DECOMPOSITION_USER_NUMBER = 1
    INTEGER(CMISSIntg), PARAMETER :: DECOMPOSER_USER_NUMBER = 1
    INTEGER(CMISSIntg), PARAMETER :: GEOMETRIC_FIELD_USER_NUMBER = 1
    INTEGER(CMISSIntg), PARAMETER :: FIBRE_FIELD_USER_NUMBER = 2
    INTEGER(CMISSIntg), PARAMETER :: MATERIAL_FIELD_USER_NUMBER = 3
    INTEGER(CMISSIntg), PARAMETER :: DEPENDENT_FIELD_USER_NUMBER = 4
    INTEGER(CMISSIntg), PARAMETER :: EQUATIONS_SET_FIELD_USER_NUMBER = 5
    INTEGER(CMISSIntg), PARAMETER :: DEFORMED_FIELD_USER_NUMBER = 6
    INTEGER(CMISSIntg), PARAMETER :: EQUATIONS_SET_USER_NUMBER = 1
    INTEGER(CMISSIntg), PARAMETER :: PROBLEM_USER_NUMBER = 1
 
    INTEGER(CMISSIntg), PARAMETER :: NUMBER_OF_GLOBAL_X_ELEMENTS = 1
    INTEGER(CMISSIntg), PARAMETER :: NUMBER_OF_GLOBAL_Y_ELEMENTS = 1
    INTEGER(CMISSIntg), PARAMETER :: NUMBER_OF_GLOBAL_Z_ELEMENTS = 1
    
    INTEGER(CMISSIntg) :: totalNumberOfNodes = 8
    INTEGER(CMISSIntg) :: totalNumberOfElements = 1
    INTEGER(CMISSIntg) :: meshComponentNumber = 1

    INTEGER(CMISSIntg) :: numberOfComputationalNodes,computationalNodeNumber
    INTEGER(CMISSIntg) :: componentIdx,err,numberOfMaterialComponents
    INTEGER(CMISSIntg) :: numberOfXi,quadratureOrder,decompositionIndex,equationsSetIndex
    INTEGER(CMISSIntg) :: numberOfGaussXi
    INTEGER(CMISSIntg) :: interpolationType
    INTEGER(CMISSIntg) :: numberOfMeshComponents
    REAL(CMISSRP) :: load
    LOGICAL :: directoryExists
    CHARACTER(LEN=255) :: outputFile,suffix

    !CMISS variables

    TYPE(cmfe_BasisType) :: basis,pressureBasis
    TYPE(cmfe_BoundaryConditionsType) :: boundaryConditions
    TYPE(cmfe_ComputationEnvironmentType) :: computationEnvironment
    TYPE(cmfe_CoordinateSystemType)  :: coordinateSystem
    TYPE(cmfe_DecomposerType) :: decomposer
    TYPE(cmfe_DecompositionType) :: decomposition
    TYPE(cmfe_EquationsType) :: equations
    TYPE(cmfe_EquationsSetType) :: equationsSet
    TYPE(cmfe_FieldType) :: geometricField,equationsSetField,fibreField
    TYPE(cmfe_FieldType) :: dependentField,materialField,deformedField
    TYPE(cmfe_FieldsType) :: fields
    TYPE(cmfe_MeshType) :: mesh
    TYPE(cmfe_GeneratedMeshType) :: generatedMesh
    TYPE(cmfe_MeshElementsType) :: meshElements
    TYPE(cmfe_NodesType) :: nodes
    TYPE(cmfe_ProblemType) :: problem
    TYPE(cmfe_RegionType) :: region,worldRegion
    TYPE(cmfe_SolverType) :: solver,nonlinearSolver,linearSolver
    TYPE(cmfe_SolverEquationsType) :: solverEquations
    TYPE(cmfe_ControlLoopType) :: controlLoop
    TYPE(cmfe_WorkGroupType) :: worldWorkGroup

    WRITE(*,'(A)') "Program starting."

    IF(usePressureBasis) THEN
      numberOfMeshComponents = 2
    ELSE
      numberOfMeshComponents = 1
    ENDIF
    IF(NUMBER_OF_GLOBAL_Z_ELEMENTS==0) THEN
      numberOfXi = 2
    ELSE
      numberOfXi = 3
    ENDIF

    IF(useSimplex) THEN
      interpolationType = 7
      quadratureOrder = 3
      numberOfGaussXi = 0
    ELSE
      interpolationType = 1
      numberOfGaussXi = 2
      quadratureOrder = 0
   ENDIF

    CALL cmfe_Region_Initialise(worldRegion,err)
    CALL cmfe_Context_WorldRegionGet(context,worldRegion,err)
    
    !Get the number of computational nodes and this computational node number
    CALL cmfe_ComputationEnvironment_Initialise(computationEnvironment,err)
    CALL cmfe_Context_ComputationEnvironmentGet(context,computationEnvironment,err)
  
    CALL cmfe_WorkGroup_Initialise(worldWorkGroup,err)
    CALL cmfe_ComputationEnvironment_WorldWorkGroupGet(computationEnvironment,worldWorkGroup,err)
    CALL cmfe_WorkGroup_NumberOfGroupNodesGet(worldWorkGroup,numberOfComputationalNodes,err)
    CALL cmfe_WorkGroup_GroupNodeNumberGet(worldWorkGroup,computationalNodeNumber,err)

    CALL cmfe_CoordinateSystem_Initialise(coordinateSystem,err)
    CALL cmfe_CoordinateSystem_CreateStart(COORDINATE_SYSTEM_USER_NUMBER,context,coordinateSystem,err)
    CALL cmfe_CoordinateSystem_DimensionSet(coordinateSystem,3,err)
    CALL cmfe_CoordinateSystem_CreateFinish(coordinateSystem,err)

    !Create a region and assign the coordinate system to the region
    CALL cmfe_Region_Initialise(region,err)
    CALL cmfe_Region_CreateStart(REGION_USER_NUMBER,worldRegion,region,err)
    CALL cmfe_Region_LabelSet(region,"Region",err)
    CALL cmfe_Region_CoordinateSystemSet(region,coordinateSystem,err)
    CALL cmfe_Region_CreateFinish(region,err)

    !Define basis
    CALL cmfe_Basis_Initialise(basis,err)
    CALL cmfe_Basis_CreateStart(BASIS_USER_NUMBER,context,basis,err)
    SELECT CASE(interpolationType)
    CASE(CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
      & CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION, &
      & CMFE_BASIS_CUBIC_LAGRANGE_INTERPOLATION, &
      & CMFE_BASIS_CUBIC_HERMITE_INTERPOLATION)
      CALL cmfe_Basis_TypeSet(basis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,err)
      CALL cmfe_Basis_NumberOfXiSet(basis,numberOfXi,err)
      IF(NUMBER_OF_GLOBAL_Z_ELEMENTS==0) THEN
        CALL cmfe_Basis_InterpolationXiSet(basis, &
          & [CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
          &  CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION],err)
      ELSE
        CALL cmfe_Basis_InterpolationXiSet(basis, &
           & [CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
           &  CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
           &  CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION],err)
      ENDIF
      IF(numberOfGaussXi>0) THEN
        IF(NUMBER_OF_GLOBAL_Z_ELEMENTS==0) THEN
          CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(basis,[numberOfGaussXi,numberOfGaussXi],err)
        ELSE
          CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(basis,[numberOfGaussXi,numberOfGaussXi,numberOfGaussXi],err)
        ENDIF
      ENDIF
    CASE(CMFE_BASIS_LINEAR_SIMPLEX_INTERPOLATION, &
      & CMFE_BASIS_QUADRATIC_SIMPLEX_INTERPOLATION, &
      & CMFE_BASIS_CUBIC_SIMPLEX_INTERPOLATION)
      CALL cmfe_Basis_TypeSet(basis,CMFE_BASIS_SIMPLEX_TYPE,err)
      CALL cmfe_Basis_NumberOfXiSet(basis,numberOfXi,err)
      IF(NUMBER_OF_GLOBAL_Z_ELEMENTS==0) THEN
        CALL cmfe_Basis_InterpolationXiSet(basis, &
          & [CMFE_BASIS_LINEAR_SIMPLEX_INTERPOLATION, &
          &  CMFE_BASIS_LINEAR_SIMPLEX_INTERPOLATION],err)
      ELSE
        CALL cmfe_Basis_InterpolationXiSet(basis, &
           & [CMFE_BASIS_LINEAR_SIMPLEX_INTERPOLATION, &
           &  CMFE_BASIS_LINEAR_SIMPLEX_INTERPOLATION, &
           &  CMFE_BASIS_LINEAR_SIMPLEX_INTERPOLATION],err)
      ENDIF
      CALL cmfe_Basis_QuadratureOrderSet(basis,quadratureOrder,err)
    CASE DEFAULT
      CALL HandelError("Invalid interpolation type.")
    END SELECT
    CALL cmfe_Basis_CreateFinish(basis,err)

    IF(usePressureBasis) THEN
      !Define pressure basis
      CALL cmfe_Basis_Initialise(pressureBasis,err)
      CALL cmfe_Basis_CreateStart(PRESSURE_BASIS_USER_NUMBER,context,pressureBasis,err)
      SELECT CASE(interpolationType)
      CASE(CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
        & CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION, &
        & CMFE_BASIS_CUBIC_LAGRANGE_INTERPOLATION, &
        & CMFE_BASIS_CUBIC_HERMITE_INTERPOLATION)
        CALL cmfe_Basis_TypeSet(pressureBasis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,err)
        CALL cmfe_Basis_NumberOfXiSet(pressureBasis,numberOfXi,err)
        IF(NUMBER_OF_GLOBAL_Z_ELEMENTS==0) THEN
          CALL cmfe_Basis_InterpolationXiSet(pressureBasis, &
            & [CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
            &  CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION],err)
        ELSE
          CALL cmfe_Basis_InterpolationXiSet(pressureBasis, &
             & [CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
             &  CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
             &  CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION],err)
        ENDIF
        IF (numberOfGaussXi>0) THEN
          IF(NUMBER_OF_GLOBAL_Z_ELEMENTS==0) THEN
            CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(pressureBasis, &
              & [numberOfGaussXi,numberOfGaussXi],err)
          ELSE
            CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(pressureBasis, &
              & [numberOfGaussXi,numberOfGaussXi,numberOfGaussXi],err)
          END IF
        END IF
      CASE(CMFE_BASIS_LINEAR_SIMPLEX_INTERPOLATION, &
        & CMFE_BASIS_QUADRATIC_SIMPLEX_INTERPOLATION, &
        & CMFE_BASIS_CUBIC_SIMPLEX_INTERPOLATION)
        CALL cmfe_Basis_TypeSet(pressureBasis,CMFE_BASIS_SIMPLEX_TYPE,err)
        CALL cmfe_Basis_NumberOfXiSet(pressureBasis,numberOfXi,err)
        IF(NUMBER_OF_GLOBAL_Z_ELEMENTS==0) THEN
          CALL cmfe_Basis_InterpolationXiSet(pressureBasis, &
            & [CMFE_BASIS_LINEAR_SIMPLEX_INTERPOLATION, &
            &  CMFE_BASIS_LINEAR_SIMPLEX_INTERPOLATION],err)
        ELSE
          CALL cmfe_Basis_InterpolationXiSet(pressureBasis, &
             & [CMFE_BASIS_LINEAR_SIMPLEX_INTERPOLATION, &
             &  CMFE_BASIS_LINEAR_SIMPLEX_INTERPOLATION, &
             &  CMFE_BASIS_LINEAR_SIMPLEX_INTERPOLATION],err)
        ENDIF
      CALL cmfe_Basis_QuadratureOrderSet(pressureBasis,quadratureOrder,err)
      CASE DEFAULT
        CALL HandelError("Invalid interpolation type.")
      END SELECT
      CALL cmfe_Basis_CreateFinish(pressureBasis,err)
    ENDIF

    CALL cmfe_Mesh_Initialise(Mesh,err)
    IF(useGeneratedMesh) THEN
      !Start the creation of a generated mesh in the region
      CALL cmfe_GeneratedMesh_Initialise(generatedMesh,err)
      CALL cmfe_GeneratedMesh_CreateStart(GENERATED_MESH_USER_NUMBER,region,generatedMesh,err)
      CALL cmfe_GeneratedMesh_TypeSet(generatedMesh,CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,err)
      IF(usePressureBasis) THEN
        CALL cmfe_GeneratedMesh_BasisSet(generatedMesh,[basis,pressureBasis],err)
      ELSE
        CALL cmfe_GeneratedMesh_BasisSet(generatedMesh,basis,err)
      ENDIF
      CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,[WIDTH,LENGTH,HEIGHT],err)
      CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh, &
        & [NUMBER_OF_GLOBAL_X_ELEMENTS,NUMBER_OF_GLOBAL_Y_ELEMENTS, &
        & NUMBER_OF_GLOBAL_Z_ELEMENTS],err)
      CALL cmfe_GeneratedMesh_CreateFinish(GeneratedMesh,MESH_USER_NUMBER,mesh,err)
    ELSE
      !Start the creation of a manually generated mesh in the region
      CALL cmfe_Mesh_CreateStart(MESH_USER_NUMBER,region,numberOfXi,mesh,err)
      CALL cmfe_Mesh_NumberOfComponentsSet(mesh,numberOfMeshComponents,err)
      IF(useSimplex) THEN
        CALL cmfe_Mesh_NumberOfElementsSet(mesh,totalNumberOfElements*5,err)
      ELSE
        CALL cmfe_Mesh_NumberOfElementsSet(mesh,totalNumberOfElements,err)
      ENDIF

      !Define nodes for the mesh
      CALL cmfe_Nodes_Initialise(nodes,err)
      CALL cmfe_Nodes_CreateStart(region,totalNumberOfNodes,nodes,err)
      CALL cmfe_Nodes_CreateFinish(nodes,err)

      CALL cmfe_MeshElements_Initialise(meshElements,err)
      CALL cmfe_MeshElements_CreateStart(mesh,meshComponentNumber,basis,meshElements,err)
      IF(useSimplex) THEN
        CALL cmfe_MeshElements_NodesSet(meshElements,1,[1,2,4,6],err)
        CALL cmfe_MeshElements_NodesSet(meshElements,2,[1,4,3,7],err)
        CALL cmfe_MeshElements_NodesSet(meshElements,3,[1,6,7,5],err)
        CALL cmfe_MeshElements_NodesSet(meshElements,4,[6,4,7,8],err)
        CALL cmfe_MeshElements_NodesSet(meshElements,5,[1,6,4,7],err)
      ELSE
        CALL cmfe_MeshElements_NodesSet(meshElements,1,[1,2,3,4,5,6,7,8],err)
      ENDIF
      CALL cmfe_MeshElements_CreateFinish(meshElements,err)

      CALL cmfe_Mesh_CreateFinish(mesh,err)
    ENDIF

    !Create a decomposition for the mesh
    CALL cmfe_Decomposition_Initialise(decomposition,err)
    CALL cmfe_Decomposition_CreateStart(DECOMPOSITION_USER_NUMBER,mesh,decomposition,err)
    CALL cmfe_Decomposition_TypeSet(decomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,err)
    CALL cmfe_Decomposition_CreateFinish(decomposition,err)

    !Decompose
    CALL cmfe_Decomposer_Initialise(decomposer,err)
    CALL cmfe_Decomposer_CreateStart(DECOMPOSER_USER_NUMBER,region,worldWorkGroup,decomposer,err)
    !Add in the decomposition
    CALL cmfe_Decomposer_DecompositionAdd(decomposer,decomposition,decompositionIndex,err)
    !Finish the decomposer
    CALL cmfe_Decomposer_CreateFinish(decomposer,err)
    
    !Create a field for the geometry
    CALL cmfe_Field_Initialise(geometricField,err)
    CALL cmfe_Field_CreateStart(GEOMETRIC_FIELD_USER_NUMBER,region,geometricField,err)
    CALL cmfe_Field_DecompositionSet(geometricField,Decomposition,err)
    CALL cmfe_Field_TypeSet(geometricField,CMFE_FIELD_GEOMETRIC_TYPE,err)
    CALL cmfe_Field_VariableLabelSet(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,"Geometry",err)
    CALL cmfe_Field_ComponentMeshComponentSet(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,err)
    CALL cmfe_Field_ComponentMeshComponentSet(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,2,1,err)
    CALL cmfe_Field_ComponentMeshComponentSet(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,3,1,err)
    IF(interpolationType==CMFE_BASIS_CUBIC_HERMITE_INTERPOLATION) THEN
      CALL cmfe_Field_ScalingTypeSet(geometricField,CMFE_FIELD_ARITHMETIC_MEAN_SCALING,err)
    END IF
    CALL cmfe_Field_CreateFinish(geometricField,err)

    IF (useGeneratedMesh) THEN
      ! Update the geometric field parameters from generated mesh
      CALL cmfe_GeneratedMesh_GeometricParametersCalculate(generatedMesh,geometricField,err)
    ELSE
      ! Update the geometric field parameters manually
      CALL cmfe_Field_ParameterSetUpdateStart(geometricField, &
        & CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,err)
      ! node 1
      CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & CMFE_FIELD_VALUES_SET_TYPE,1,1,1,1,0.0_CMISSRP,err)
      CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & CMFE_FIELD_VALUES_SET_TYPE,1,1,1,2,0.0_CMISSRP,err)
      CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & CMFE_FIELD_VALUES_SET_TYPE,1,1,1,3,0.0_CMISSRP,err)
      ! node 2
      CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & CMFE_FIELD_VALUES_SET_TYPE,1,1,2,1,HEIGHT,err)
      CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & CMFE_FIELD_VALUES_SET_TYPE,1,1,2,2,0.0_CMISSRP,err)
      CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & CMFE_FIELD_VALUES_SET_TYPE,1,1,2,3,0.0_CMISSRP,err)
      ! node 3
      CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & CMFE_FIELD_VALUES_SET_TYPE,1,1,3,1,0.0_CMISSRP,err)
      CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & CMFE_FIELD_VALUES_SET_TYPE,1,1,3,2,WIDTH,err)
      CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & CMFE_FIELD_VALUES_SET_TYPE,1,1,3,3,0.0_CMISSRP,err)
      ! node 4
      CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & CMFE_FIELD_VALUES_SET_TYPE,1,1,4,1,HEIGHT,err)
      CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & CMFE_FIELD_VALUES_SET_TYPE,1,1,4,2,WIDTH,err)
      CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & CMFE_FIELD_VALUES_SET_TYPE,1,1,4,3,0.0_CMISSRP,err)
      ! node 5
      CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & CMFE_FIELD_VALUES_SET_TYPE,1,1,5,1,0.0_CMISSRP,err)
      CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & CMFE_FIELD_VALUES_SET_TYPE,1,1,5,2,0.0_CMISSRP,err)
      CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & CMFE_FIELD_VALUES_SET_TYPE,1,1,5,3,LENGTH,err)
      ! node 6
      CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & CMFE_FIELD_VALUES_SET_TYPE,1,1,6,1,HEIGHT,err)
      CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & CMFE_FIELD_VALUES_SET_TYPE,1,1,6,2,0.0_CMISSRP,err)
      CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & CMFE_FIELD_VALUES_SET_TYPE,1,1,6,3,LENGTH,err)
      ! node 7
      CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & CMFE_FIELD_VALUES_SET_TYPE,1,1,7,1,0.0_CMISSRP,err)
      CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & CMFE_FIELD_VALUES_SET_TYPE,1,1,7,2,WIDTH,err)
      CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & CMFE_FIELD_VALUES_SET_TYPE,1,1,7,3,LENGTH,err)
      ! node 8
      CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & CMFE_FIELD_VALUES_SET_TYPE,1,1,8,1,HEIGHT,err)
      CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & CMFE_FIELD_VALUES_SET_TYPE,1,1,8,2,WIDTH,err)
      CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & CMFE_FIELD_VALUES_SET_TYPE,1,1,8,3,LENGTH,err)
      CALL cmfe_Field_ParameterSetUpdateFinish(geometricField, &
        & CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,err)
    ENDIF

    !Create a fibre field and attach it to the geometric field
    CALL cmfe_Field_Initialise(fibreField,err)
    CALL cmfe_Field_CreateStart(FIBRE_FIELD_USER_NUMBER,region,fibreField,err)
    CALL cmfe_Field_TypeSet(fibreField,CMFE_FIELD_FIBRE_TYPE,err)
    CALL cmfe_Field_DecompositionSet(fibreField,decomposition,err)
    CALL cmfe_Field_GeometricFieldSet(fibreField,geometricField,err)
    CALL cmfe_Field_VariableLabelSet(fibreField,CMFE_FIELD_U_VARIABLE_TYPE,"Fibre",err)
    IF(interpolationType==CMFE_BASIS_CUBIC_HERMITE_INTERPOLATION) THEN
      CALL cmfe_Field_ScalingTypeSet(fibreField,CMFE_FIELD_ARITHMETIC_MEAN_SCALING,err)
    ENDIF
    CALL cmfe_Field_CreateFinish(fibreField,err)

    !Create the material field
    IF(compressible) THEN
      numberOfMaterialComponents = 3
    ELSE
      numberOfMaterialComponents = 2
    ENDIF
    CALL cmfe_Field_Initialise(materialField,err)
    CALL cmfe_Field_CreateStart(MATERIAL_FIELD_USER_NUMBER,Region,materialField,err)
    CALL cmfe_Field_TypeSet(materialField,CMFE_FIELD_MATERIAL_TYPE,err)
    CALL cmfe_Field_DecompositionSet(materialField,decomposition,err)
    CALL cmfe_Field_GeometricFieldSet(materialField,geometricField,err)
    CALL cmfe_Field_NumberOfVariablesSet(materialField,1,err)
    CALL cmfe_Field_NumberOfComponentsSet(materialField,CMFE_FIELD_U_VARIABLE_TYPE,numberOfMaterialComponents,err)
    CALL cmfe_Field_VariableLabelSet(materialField,CMFE_FIELD_U_VARIABLE_TYPE,"Material",err)
    CALL cmfe_Field_ComponentMeshComponentSet(materialField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,err)
    CALL cmfe_Field_ComponentMeshComponentSet(materialField,CMFE_FIELD_U_VARIABLE_TYPE,2,1,err)
    IF(compressible) THEN
      CALL cmfe_Field_ComponentMeshComponentSet(materialField,CMFE_FIELD_U_VARIABLE_TYPE,3,1,err)
    ENDIF
    IF(interpolationType==CMFE_BASIS_CUBIC_HERMITE_INTERPOLATION) THEN
      CALL cmfe_Field_ScalingTypeSet(materialField,CMFE_FIELD_ARITHMETIC_MEAN_SCALING,err)
    ENDIF
    CALL cmfe_Field_CreateFinish(materialField,err)

    !Set Mooney-Rivlin constants c10 and c01 respectively.
    CALL cmfe_Field_ComponentValuesInitialise(materialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
      & 1,2.0_CMISSRP,err)
    CALL cmfe_Field_ComponentValuesInitialise(materialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
      & 2,6.0_CMISSRP,err)
    IF(compressible) THEN
      CALL cmfe_Field_ComponentValuesInitialise(materialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
        & 3,1.0e9_CMISSRP,err)
    ENDIF

    !Create the dependent field
    IF(compressible) THEN
      numberOfMaterialComponents = 3
    ELSE
      numberOfMaterialComponents = 4
    ENDIF

    CALL cmfe_Field_Initialise(dependentField,err)
    CALL cmfe_Field_CreateStart(DEPENDENT_FIELD_USER_NUMBER,region,dependentField,err)
    CALL cmfe_Field_VariableLabelSet(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,"Dependent",err)
    CALL cmfe_Field_TypeSet(dependentField,CMFE_FIELD_GEOMETRIC_GENERAL_TYPE,err)
    CALL cmfe_Field_DecompositionSet(dependentField,decomposition,err)
    CALL cmfe_Field_GeometricFieldSet(dependentField,geometricField,err)
    CALL cmfe_Field_DependentTypeSet(dependentField,CMFE_FIELD_DEPENDENT_TYPE,err)
    CALL cmfe_Field_NumberOfVariablesSet(dependentField,2,err)
    CALL cmfe_Field_NumberOfComponentsSet(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,numberOfMaterialComponents,err)
    CALL cmfe_Field_NumberOfComponentsSet(dependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,numberOfMaterialComponents,err)
    CALL cmfe_Field_ComponentMeshComponentSet(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,err)
    CALL cmfe_Field_ComponentMeshComponentSet(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,2,1,err)
    CALL cmfe_Field_ComponentMeshComponentSet(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,3,1,err)
    CALL cmfe_Field_ComponentMeshComponentSet(dependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,1,err)
    CALL cmfe_Field_ComponentMeshComponentSet(dependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,2,1,err)
    CALL cmfe_Field_ComponentMeshComponentSet(dependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,3,1,err)
    IF(.NOT.compressible) THEN
      ! TODO we always had node based interpolation, right? --> check!
      CALL cmfe_Field_ComponentInterpolationSet(dependentField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & 4,CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,err)
      CALL cmfe_Field_ComponentInterpolationSet(dependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE, &
        & 4,CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,err)
      ! TODO end
      IF(usePressureBasis) THEN
        !Set the pressure to be nodally based and use the second mesh component
        IF(interpolationType==4) THEN
          CALL cmfe_Field_ComponentInterpolationSet(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,4, &
            & CMFE_FIELD_NODE_BASED_INTERPOLATION,err)
          CALL cmfe_Field_ComponentInterpolationSet(dependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,4, &
            & CMFE_FIELD_NODE_BASED_INTERPOLATION,err)
        ENDIF
        CALL cmfe_Field_ComponentInterpolationSet(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,4,2,err)
        CALL cmfe_Field_ComponentInterpolationSet(dependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,4,2,err)
      ENDIF
    ENDIF
    IF(interpolationType==CMFE_BASIS_CUBIC_HERMITE_INTERPOLATION) THEN
      CALL cmfe_Field_ScalingTypeSet(dependentField,CMFE_FIELD_ARITHMETIC_MEAN_SCALING,err)
    ENDIF
    CALL cmfe_Field_CreateFinish(dependentField,err)

    !Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
    CALL cmfe_Field_ParametersToFieldParametersComponentCopy( &
      & geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
      & dependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,err)
    CALL cmfe_Field_ParametersToFieldParametersComponentCopy( &
      & geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2, &
      & dependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,err)
    CALL cmfe_Field_ParametersToFieldParametersComponentCopy( &
      & geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3, &
      & dependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3,err)
    IF(.NOT.compressible) THEN
      CALL cmfe_Field_ComponentValuesInitialise(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
        & 4,0.0_CMISSRP,err)
    ENDIF

    !Create a deformed geometry field, as cmgui doesn't like displaying
    ! deformed fibres from the dependent field because it isn't a geometric field.
    CALL cmfe_Field_Initialise(deformedField,err)
    CALL cmfe_Field_CreateStart(DEFORMED_FIELD_USER_NUMBER,region,deformedField,err)
    CALL cmfe_Field_DecompositionSet(deformedField,decomposition,err)
    CALL cmfe_Field_TypeSet(deformedField,CMFE_FIELD_GEOMETRIC_TYPE,err)
    CALL cmfe_Field_VariableLabelSet(deformedField,CMFE_FIELD_U_VARIABLE_TYPE,"DeformedGeometry",err)
    DO componentIdx=1,3
      CALL cmfe_Field_ComponentMeshComponentSet(deformedField,CMFE_FIELD_U_VARIABLE_TYPE,componentIdx,1,err)
    ENDDO
    IF(interpolationType==CMFE_BASIS_CUBIC_HERMITE_INTERPOLATION) THEN
      CALL cmfe_Field_ScalingTypeSet(deformedField,CMFE_FIELD_ARITHMETIC_MEAN_SCALING,err)
    ENDIF
    CALL cmfe_Field_CreateFinish(deformedField,err)

    !Create the equations set
    CALL cmfe_Field_Initialise(equationsSetField,err)
    CALL cmfe_EquationsSet_Initialise(equationsSet,err)
    IF(compressible) THEN
      CALL cmfe_EquationsSet_CreateStart(EQUATIONS_SET_USER_NUMBER,region,fibreField, &
        & [CMFE_EQUATIONS_SET_ELASTICITY_CLASS, &
        &  CMFE_EQUATIONS_SET_FINITE_ELASTICITY_TYPE, &
        &  CMFE_EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE], &
        & EQUATIONS_SET_FIELD_USER_NUMBER,equationsSetField,equationsSet,err)
    ELSE
      CALL cmfe_EquationsSet_CreateStart(EQUATIONS_SET_USER_NUMBER,region,fibreField, &
        & [CMFE_EQUATIONS_SET_ELASTICITY_CLASS, &
        &  CMFE_EQUATIONS_SET_FINITE_ELASTICITY_TYPE, &
        &  CMFE_EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE], &
        & EQUATIONS_SET_FIELD_USER_NUMBER,equationsSetField,equationsSet,err)
    ENDIF
    CALL cmfe_EquationsSet_CreateFinish(equationsSet,err)
    CALL cmfe_EquationsSet_MaterialsCreateStart(equationsSet,MATERIAL_FIELD_USER_NUMBER,materialField,err)
    CALL cmfe_EquationsSet_MaterialsCreateFinish(equationsSet,err)
    CALL cmfe_EquationsSet_DependentCreateStart(equationsSet,DEPENDENT_FIELD_USER_NUMBER,dependentField,err)
    CALL cmfe_EquationsSet_DependentCreateFinish(equationsSet,err)

    !Create equations
    CALL cmfe_Equations_Initialise(equations,err)
    CALL cmfe_EquationsSet_EquationsCreateStart(equationsSet,equations,err)
    CALL cmfe_Equations_SparsityTypeSet(equations,CMFE_EQUATIONS_SPARSE_MATRICES,err)
    CALL cmfe_Equations_OutputTypeSet(equations,CMFE_EQUATIONS_NO_OUTPUT,err)
    CALL cmfe_EquationsSet_EquationsCreateFinish(equationsSet,err)

    CALL cmfe_Equations_JacobianCalculationTypeSet(equations,1,CMFE_FIELD_U_VARIABLE_TYPE, &
      & CMFE_EQUATIONS_JACOBIAN_ANALYTIC_CALCULATED,err)

    !Define the problem
    CALL cmfe_Problem_Initialise(problem,err)
    CALL cmfe_Problem_CreateStart(PROBLEM_USER_NUMBER,context, &
      & [CMFE_PROBLEM_ELASTICITY_CLASS, &
      &  CMFE_PROBLEM_FINITE_ELASTICITY_TYPE, &
      &  CMFE_PROBLEM_STATIC_FINITE_ELASTICITY_SUBTYPE],problem,err)
    CALL cmfe_Problem_CreateFinish(problem,err)

    !Create control loops
    CALL cmfe_Problem_ControlLoopCreateStart(problem,err)
    CALL cmfe_Problem_ControlLoopCreateFinish(problem,err)

    !Create problem solver
    CALL cmfe_Solver_Initialise(nonlinearSolver,err)
    CALL cmfe_Solver_Initialise(linearSolver,err)
    CALL cmfe_Problem_SolversCreateStart(problem,err)
    CALL cmfe_Problem_SolverGet(problem,CMFE_CONTROL_LOOP_NODE,1,nonLinearSolver,err)
    CALL cmfe_Solver_OutputTypeSet(nonlinearSolver,CMFE_SOLVER_PROGRESS_OUTPUT,err)
    !CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(nonlinearSolver, &
    !  & CMFE_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,err)
    CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(nonlinearSolver, &
      & CMFE_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED,err)
    CALL cmfe_Solver_NewtonLinearSolverGet(nonlinearSolver,linearSolver,err)
    CALL cmfe_Solver_NewtonAbsoluteToleranceSet(nonlinearSolver,1.0E-14_CMISSRP,err)
    CALL cmfe_Solver_NewtonSolutionToleranceSet(nonlinearSolver,1.0E-14_CMISSRP,err)
    CALL cmfe_Solver_NewtonRelativeToleranceSet(nonlinearSolver,1.0E-14_CMISSRP,err)
    CALL cmfe_Solver_LinearTypeSet(linearSolver,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,err)
    CALL cmfe_Problem_SolversCreateFinish(problem,err)

    !Create solver equations and add equations set to solver equations
    CALL cmfe_Solver_Initialise(solver,err)
    CALL cmfe_SolverEquations_Initialise(solverEquations,err)
    CALL cmfe_Problem_SolverEquationsCreateStart(problem,err)
    CALL cmfe_Problem_SolverGet(problem,CMFE_CONTROL_LOOP_NODE,1,solver,err)
    CALL cmfe_Solver_SolverEquationsGet(solver,solverEquations,err)
    CALL cmfe_SolverEquations_SparsityTypeSet(solverEquations,CMFE_SOLVER_SPARSE_MATRICES,err)
    CALL cmfe_SolverEquations_EquationsSetAdd(solverEquations,equationsSet,equationsSetIndex,err)
    CALL cmfe_Problem_SolverEquationsCreateFinish(problem,err)

    !Prescribe boundary conditions (absolute nodal parameters)
    CALL cmfe_BoundaryConditions_Initialise(boundaryConditions,err)
    CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(solverEquations,boundaryConditions,err)

    !Set x=0 nodes to no x displacment in x. Set x=WIDTH nodes to 10% x displacement
    CALL cmfe_BoundaryConditions_AddNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE, &
      & 1,1,1,1,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,err)
    CALL cmfe_BoundaryConditions_AddNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE, &
      & 1,1,3,1,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,err)
    CALL cmfe_BoundaryConditions_AddNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE, &
      & 1,1,5,1,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,err)
    CALL cmfe_BoundaryConditions_AddNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE, &
      & 1,1,7,1,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,err)
    IF (zeroLoad) THEN
      load = 0.0_CMISSRP
    ELSE
      load = 0.1_CMISSRP*WIDTH
    END IF
    CALL cmfe_BoundaryConditions_AddNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE, &
      & 1,1,2,1,CMFE_BOUNDARY_CONDITION_FIXED,load,err)
    CALL cmfe_BoundaryConditions_AddNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE, &
      & 1,1,4,1,CMFE_BOUNDARY_CONDITION_FIXED,load,err)
    CALL cmfe_BoundaryConditions_AddNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE, &
      & 1,1,6,1,CMFE_BOUNDARY_CONDITION_FIXED,load,err)
    CALL cmfe_BoundaryConditions_AddNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE, &
      & 1,1,8,1,CMFE_BOUNDARY_CONDITION_FIXED,load,err)

    ! Set y=0 nodes to no y displacement
    CALL cmfe_BoundaryConditions_AddNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE, &
      & 1,1,1,2,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,err)
    CALL cmfe_BoundaryConditions_AddNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE, &
      & 1,1,2,2,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,err)
    CALL cmfe_BoundaryConditions_AddNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE, &
      & 1,1,5,2,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,err)
    CALL cmfe_BoundaryConditions_AddNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE, &
      & 1,1,6,2,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,err)

    ! Set z=0 nodes to no y displacement
    CALL cmfe_BoundaryConditions_AddNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE, &
      & 1,1,1,3,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,err)
    CALL cmfe_BoundaryConditions_AddNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE, &
      & 1,1,2,3,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,err)
    CALL cmfe_BoundaryConditions_AddNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE, &
      & 1,1,3,3,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,err)
    CALL cmfe_BoundaryConditions_AddNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE, &
      & 1,1,4,3,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,err)

    CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(solverEquations,err)

    !Solve the problem
    CALL cmfe_Problem_Solve(problem,err)

    !Copy deformed geometry into deformed field
    DO componentIdx=1,3
      CALL cmfe_Field_ParametersToFieldParametersComponentCopy( &
        & dependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,componentIdx, &
        & deformedField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,componentIdx,err)
    ENDDO

    IF(useGeneratedMesh) THEN
      outputFile = "./results/unit_cube_generated_mesh"
    ELSE
      IF(useSimplex) THEN
        outputFile = "./results/unit_cube_manual_mesh_simplex"
      ELSE
        outputFile = "./results/unit_cube_manual_mesh"
      ENDIF
    ENDIF

    INQUIRE(FILE="./results",EXIST=directoryExists)
    IF(.NOT.directoryExists) CALL EXECUTE_COMMAND_LINE("mkdir ./results")

    suffix = ''
    IF(compressible) THEN
      IF(zeroLoad) THEN
        suffix = "_compressible_zero_load"
      ELSE
        suffix = "_compressible"
      ENDIF
    ELSE
      IF(zeroLoad) THEN
        suffix = "_zero_load"
      ENDIF
    ENDIF

    !Export results
    CALL cmfe_Fields_Initialise(fields,err)
    CALL cmfe_Fields_Create(region,fields,err)
    CALL cmfe_Fields_NodesExport(fields,trim(outputFile)//trim(suffix),"FORTRAN",err)
    CALL cmfe_Fields_ElementsExport(fields,trim(outputFile)//trim(suffix),"FORTRAN",err)
    CALL cmfe_Fields_Finalise(fields,err)

    CALL cmfe_Problem_Destroy(problem,err)
    IF(usePressureBasis) CALL cmfe_Basis_Destroy(pressureBasis,err)
    CALL cmfe_Basis_Destroy(basis,err)
    CALL cmfe_Region_Destroy(region,err)
    CALL cmfe_CoordinateSystem_Destroy(coordinateSystem,err)
   
    WRITE(*,'(A)') "Program successfully completed."

  END SUBROUTINE SolveModel

  SUBROUTINE HandelError(errorString)
    
    CHARACTER(LEN=*), INTENT(IN) :: errorString

    WRITE(*,'(">>ERROR: ",A)') errorString(1:LEN_TRIM(errorString))
    STOP
    
  END SUBROUTINE HandelError

END PROGRAM UniaxialExtension
