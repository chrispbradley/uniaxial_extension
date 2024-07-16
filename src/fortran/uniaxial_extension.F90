PROGRAM UniaxialExtension

  USE OpenCMISS

  IMPLICIT NONE

  INTEGER(OC_Intg), PARAMETER :: CONTEXT_USER_NUMBER = 1

  TYPE(OC_ContextType) :: context
  
  INTEGER(OC_Intg) :: err

  !Intialise OpenCMISS
  CALL OC_Initialise(err)
  CALL OC_ErrorHandlingModeSet(OC_ERRORS_TRAP_ERROR,err)
  !Create a context
  CALL OC_Context_Initialise(context,err)
  CALL OC_Context_Create(CONTEXT_USER_NUMBER,context,err)
  
  !Set all diganostic levels on for testing
  !CALL OC_DiagnosticsSetOn(OC_FROM_DIAG_TYPE,[1,2,3,4,5],"Diagnostics", &
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
  CALL OC_Context_Destroy(context,err)
  !Finalise OpenCMISS
  CALL OC_Finalise(err)

  STOP

CONTAINS

  SUBROUTINE SolveModel(context, numberOfDimensions, compressible, useGeneratedMesh, zeroLoad, useSimplex, usePressureBasis)

    TYPE(OC_ContextType), INTENT(IN) :: context
    INTEGER(OC_Intg), INTENT(IN)  :: numberOfDimensions
    LOGICAL, INTENT(IN) :: compressible
    LOGICAL, INTENT(IN) :: useGeneratedMesh
    LOGICAL, INTENT(IN) :: zeroLoad
    LOGICAL, INTENT(IN) :: useSimplex
    LOGICAL, INTENT(IN) :: usePressureBasis

    REAL(OC_RP), PARAMETER :: HEIGHT = 1.0_OC_RP
    REAL(OC_RP), PARAMETER :: WIDTH = 1.0_OC_RP
    REAL(OC_RP), PARAMETER :: LENGTH = 1.0_OC_RP

    INTEGER(OC_Intg), PARAMETER :: COORDINATE_SYSTEM_USER_NUMBER = 1
    INTEGER(OC_Intg), PARAMETER :: REGION_USER_NUMBER = 1
    INTEGER(OC_Intg), PARAMETER :: BASIS_USER_NUMBER = 1
    INTEGER(OC_Intg), PARAMETER :: PRESSURE_BASIS_USER_NUMBER = 2
    INTEGER(OC_Intg), PARAMETER :: GENERATED_MESH_USER_NUMBER = 1
    INTEGER(OC_Intg), PARAMETER :: MESH_USER_NUMBER = 1
    INTEGER(OC_Intg), PARAMETER :: DECOMPOSITION_USER_NUMBER = 1
    INTEGER(OC_Intg), PARAMETER :: DECOMPOSER_USER_NUMBER = 1
    INTEGER(OC_Intg), PARAMETER :: GEOMETRIC_FIELD_USER_NUMBER = 1
    INTEGER(OC_Intg), PARAMETER :: FIBRE_FIELD_USER_NUMBER = 2
    INTEGER(OC_Intg), PARAMETER :: MATERIAL_FIELD_USER_NUMBER = 3
    INTEGER(OC_Intg), PARAMETER :: DEPENDENT_FIELD_USER_NUMBER = 4
    INTEGER(OC_Intg), PARAMETER :: EQUATIONS_SET_FIELD_USER_NUMBER = 5
    INTEGER(OC_Intg), PARAMETER :: DEFORMED_FIELD_USER_NUMBER = 6
    INTEGER(OC_Intg), PARAMETER :: EQUATIONS_SET_USER_NUMBER = 1
    INTEGER(OC_Intg), PARAMETER :: PROBLEM_USER_NUMBER = 1
 
    INTEGER(OC_Intg), PARAMETER :: NUMBER_OF_GLOBAL_X_ELEMENTS = 1
    INTEGER(OC_Intg), PARAMETER :: NUMBER_OF_GLOBAL_Y_ELEMENTS = 1
    INTEGER(OC_Intg), PARAMETER :: NUMBER_OF_GLOBAL_Z_ELEMENTS = 1
    
    INTEGER(OC_Intg) :: totalNumberOfNodes = 8
    INTEGER(OC_Intg) :: totalNumberOfElements = 1
    INTEGER(OC_Intg) :: meshComponentNumber = 1

    INTEGER(OC_Intg) :: numberOfComputationalNodes,computationalNodeNumber
    INTEGER(OC_Intg) :: componentIdx,err,numberOfMaterialComponents
    INTEGER(OC_Intg) :: numberOfXi,quadratureOrder,decompositionIndex,equationsSetIndex
    INTEGER(OC_Intg) :: numberOfGaussXi
    INTEGER(OC_Intg) :: interpolationType
    INTEGER(OC_Intg) :: numberOfMeshComponents
    REAL(OC_RP) :: load
    LOGICAL :: directoryExists
    CHARACTER(LEN=255) :: outputFile,suffix

    !CMISS variables

    TYPE(OC_BasisType) :: basis,pressureBasis
    TYPE(OC_BoundaryConditionsType) :: boundaryConditions
    TYPE(OC_ComputationEnvironmentType) :: computationEnvironment
    TYPE(OC_CoordinateSystemType)  :: coordinateSystem
    TYPE(OC_DecomposerType) :: decomposer
    TYPE(OC_DecompositionType) :: decomposition
    TYPE(OC_EquationsType) :: equations
    TYPE(OC_EquationsSetType) :: equationsSet
    TYPE(OC_FieldType) :: geometricField,equationsSetField,fibreField
    TYPE(OC_FieldType) :: dependentField,materialField,deformedField
    TYPE(OC_FieldsType) :: fields
    TYPE(OC_MeshType) :: mesh
    TYPE(OC_GeneratedMeshType) :: generatedMesh
    TYPE(OC_MeshElementsType) :: meshElements
    TYPE(OC_NodesType) :: nodes
    TYPE(OC_ProblemType) :: problem
    TYPE(OC_RegionType) :: region,worldRegion
    TYPE(OC_SolverType) :: solver,nonlinearSolver,linearSolver
    TYPE(OC_SolverEquationsType) :: solverEquations
    TYPE(OC_ControlLoopType) :: controlLoop
    TYPE(OC_WorkGroupType) :: worldWorkGroup

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

    CALL OC_Region_Initialise(worldRegion,err)
    CALL OC_Context_WorldRegionGet(context,worldRegion,err)
    
    !Get the number of computational nodes and this computational node number
    CALL OC_ComputationEnvironment_Initialise(computationEnvironment,err)
    CALL OC_Context_ComputationEnvironmentGet(context,computationEnvironment,err)
  
    CALL OC_WorkGroup_Initialise(worldWorkGroup,err)
    CALL OC_ComputationEnvironment_WorldWorkGroupGet(computationEnvironment,worldWorkGroup,err)
    CALL OC_WorkGroup_NumberOfGroupNodesGet(worldWorkGroup,numberOfComputationalNodes,err)
    CALL OC_WorkGroup_GroupNodeNumberGet(worldWorkGroup,computationalNodeNumber,err)

    CALL OC_CoordinateSystem_Initialise(coordinateSystem,err)
    CALL OC_CoordinateSystem_CreateStart(COORDINATE_SYSTEM_USER_NUMBER,context,coordinateSystem,err)
    CALL OC_CoordinateSystem_DimensionSet(coordinateSystem,3,err)
    CALL OC_CoordinateSystem_CreateFinish(coordinateSystem,err)

    !Create a region and assign the coordinate system to the region
    CALL OC_Region_Initialise(region,err)
    CALL OC_Region_CreateStart(REGION_USER_NUMBER,worldRegion,region,err)
    CALL OC_Region_LabelSet(region,"Region",err)
    CALL OC_Region_CoordinateSystemSet(region,coordinateSystem,err)
    CALL OC_Region_CreateFinish(region,err)

    !Define basis
    CALL OC_Basis_Initialise(basis,err)
    CALL OC_Basis_CreateStart(BASIS_USER_NUMBER,context,basis,err)
    SELECT CASE(interpolationType)
    CASE(OC_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
      & OC_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION, &
      & OC_BASIS_CUBIC_LAGRANGE_INTERPOLATION, &
      & OC_BASIS_CUBIC_HERMITE_INTERPOLATION)
      CALL OC_Basis_TypeSet(basis,OC_BASIS_LAGRANGE_HERMITE_TP_TYPE,err)
      CALL OC_Basis_NumberOfXiSet(basis,numberOfXi,err)
      IF(NUMBER_OF_GLOBAL_Z_ELEMENTS==0) THEN
        CALL OC_Basis_InterpolationXiSet(basis, &
          & [OC_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
          &  OC_BASIS_LINEAR_LAGRANGE_INTERPOLATION],err)
      ELSE
        CALL OC_Basis_InterpolationXiSet(basis, &
           & [OC_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
           &  OC_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
           &  OC_BASIS_LINEAR_LAGRANGE_INTERPOLATION],err)
      ENDIF
      IF(numberOfGaussXi>0) THEN
        IF(NUMBER_OF_GLOBAL_Z_ELEMENTS==0) THEN
          CALL OC_Basis_QuadratureNumberOfGaussXiSet(basis,[numberOfGaussXi,numberOfGaussXi],err)
        ELSE
          CALL OC_Basis_QuadratureNumberOfGaussXiSet(basis,[numberOfGaussXi,numberOfGaussXi,numberOfGaussXi],err)
        ENDIF
      ENDIF
    CASE(OC_BASIS_LINEAR_SIMPLEX_INTERPOLATION, &
      & OC_BASIS_QUADRATIC_SIMPLEX_INTERPOLATION, &
      & OC_BASIS_CUBIC_SIMPLEX_INTERPOLATION)
      CALL OC_Basis_TypeSet(basis,OC_BASIS_SIMPLEX_TYPE,err)
      CALL OC_Basis_NumberOfXiSet(basis,numberOfXi,err)
      IF(NUMBER_OF_GLOBAL_Z_ELEMENTS==0) THEN
        CALL OC_Basis_InterpolationXiSet(basis, &
          & [OC_BASIS_LINEAR_SIMPLEX_INTERPOLATION, &
          &  OC_BASIS_LINEAR_SIMPLEX_INTERPOLATION],err)
      ELSE
        CALL OC_Basis_InterpolationXiSet(basis, &
           & [OC_BASIS_LINEAR_SIMPLEX_INTERPOLATION, &
           &  OC_BASIS_LINEAR_SIMPLEX_INTERPOLATION, &
           &  OC_BASIS_LINEAR_SIMPLEX_INTERPOLATION],err)
      ENDIF
      CALL OC_Basis_QuadratureOrderSet(basis,quadratureOrder,err)
    CASE DEFAULT
      CALL HandelError("Invalid interpolation type.")
    END SELECT
    CALL OC_Basis_CreateFinish(basis,err)

    IF(usePressureBasis) THEN
      !Define pressure basis
      CALL OC_Basis_Initialise(pressureBasis,err)
      CALL OC_Basis_CreateStart(PRESSURE_BASIS_USER_NUMBER,context,pressureBasis,err)
      SELECT CASE(interpolationType)
      CASE(OC_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
        & OC_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION, &
        & OC_BASIS_CUBIC_LAGRANGE_INTERPOLATION, &
        & OC_BASIS_CUBIC_HERMITE_INTERPOLATION)
        CALL OC_Basis_TypeSet(pressureBasis,OC_BASIS_LAGRANGE_HERMITE_TP_TYPE,err)
        CALL OC_Basis_NumberOfXiSet(pressureBasis,numberOfXi,err)
        IF(NUMBER_OF_GLOBAL_Z_ELEMENTS==0) THEN
          CALL OC_Basis_InterpolationXiSet(pressureBasis, &
            & [OC_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
            &  OC_BASIS_LINEAR_LAGRANGE_INTERPOLATION],err)
        ELSE
          CALL OC_Basis_InterpolationXiSet(pressureBasis, &
             & [OC_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
             &  OC_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
             &  OC_BASIS_LINEAR_LAGRANGE_INTERPOLATION],err)
        ENDIF
        IF (numberOfGaussXi>0) THEN
          IF(NUMBER_OF_GLOBAL_Z_ELEMENTS==0) THEN
            CALL OC_Basis_QuadratureNumberOfGaussXiSet(pressureBasis, &
              & [numberOfGaussXi,numberOfGaussXi],err)
          ELSE
            CALL OC_Basis_QuadratureNumberOfGaussXiSet(pressureBasis, &
              & [numberOfGaussXi,numberOfGaussXi,numberOfGaussXi],err)
          END IF
        END IF
      CASE(OC_BASIS_LINEAR_SIMPLEX_INTERPOLATION, &
        & OC_BASIS_QUADRATIC_SIMPLEX_INTERPOLATION, &
        & OC_BASIS_CUBIC_SIMPLEX_INTERPOLATION)
        CALL OC_Basis_TypeSet(pressureBasis,OC_BASIS_SIMPLEX_TYPE,err)
        CALL OC_Basis_NumberOfXiSet(pressureBasis,numberOfXi,err)
        IF(NUMBER_OF_GLOBAL_Z_ELEMENTS==0) THEN
          CALL OC_Basis_InterpolationXiSet(pressureBasis, &
            & [OC_BASIS_LINEAR_SIMPLEX_INTERPOLATION, &
            &  OC_BASIS_LINEAR_SIMPLEX_INTERPOLATION],err)
        ELSE
          CALL OC_Basis_InterpolationXiSet(pressureBasis, &
             & [OC_BASIS_LINEAR_SIMPLEX_INTERPOLATION, &
             &  OC_BASIS_LINEAR_SIMPLEX_INTERPOLATION, &
             &  OC_BASIS_LINEAR_SIMPLEX_INTERPOLATION],err)
        ENDIF
      CALL OC_Basis_QuadratureOrderSet(pressureBasis,quadratureOrder,err)
      CASE DEFAULT
        CALL HandelError("Invalid interpolation type.")
      END SELECT
      CALL OC_Basis_CreateFinish(pressureBasis,err)
    ENDIF

    CALL OC_Mesh_Initialise(Mesh,err)
    IF(useGeneratedMesh) THEN
      !Start the creation of a generated mesh in the region
      CALL OC_GeneratedMesh_Initialise(generatedMesh,err)
      CALL OC_GeneratedMesh_CreateStart(GENERATED_MESH_USER_NUMBER,region,generatedMesh,err)
      CALL OC_GeneratedMesh_TypeSet(generatedMesh,OC_GENERATED_MESH_REGULAR_MESH_TYPE,err)
      IF(usePressureBasis) THEN
        CALL OC_GeneratedMesh_BasisSet(generatedMesh,[basis,pressureBasis],err)
      ELSE
        CALL OC_GeneratedMesh_BasisSet(generatedMesh,basis,err)
      ENDIF
      CALL OC_GeneratedMesh_ExtentSet(GeneratedMesh,[WIDTH,LENGTH,HEIGHT],err)
      CALL OC_GeneratedMesh_NumberOfElementsSet(GeneratedMesh, &
        & [NUMBER_OF_GLOBAL_X_ELEMENTS,NUMBER_OF_GLOBAL_Y_ELEMENTS, &
        & NUMBER_OF_GLOBAL_Z_ELEMENTS],err)
      CALL OC_GeneratedMesh_CreateFinish(GeneratedMesh,MESH_USER_NUMBER,mesh,err)
    ELSE
      !Start the creation of a manually generated mesh in the region
      CALL OC_Mesh_CreateStart(MESH_USER_NUMBER,region,numberOfXi,mesh,err)
      CALL OC_Mesh_NumberOfComponentsSet(mesh,numberOfMeshComponents,err)
      IF(useSimplex) THEN
        CALL OC_Mesh_NumberOfElementsSet(mesh,totalNumberOfElements*5,err)
      ELSE
        CALL OC_Mesh_NumberOfElementsSet(mesh,totalNumberOfElements,err)
      ENDIF

      !Define nodes for the mesh
      CALL OC_Nodes_Initialise(nodes,err)
      CALL OC_Nodes_CreateStart(region,totalNumberOfNodes,nodes,err)
      CALL OC_Nodes_CreateFinish(nodes,err)

      CALL OC_MeshElements_Initialise(meshElements,err)
      CALL OC_MeshElements_CreateStart(mesh,meshComponentNumber,basis,meshElements,err)
      IF(useSimplex) THEN
        CALL OC_MeshElements_NodesSet(meshElements,1,[1,2,4,6],err)
        CALL OC_MeshElements_NodesSet(meshElements,2,[1,4,3,7],err)
        CALL OC_MeshElements_NodesSet(meshElements,3,[1,6,7,5],err)
        CALL OC_MeshElements_NodesSet(meshElements,4,[6,4,7,8],err)
        CALL OC_MeshElements_NodesSet(meshElements,5,[1,6,4,7],err)
      ELSE
        CALL OC_MeshElements_NodesSet(meshElements,1,[1,2,3,4,5,6,7,8],err)
      ENDIF
      CALL OC_MeshElements_CreateFinish(meshElements,err)

      CALL OC_Mesh_CreateFinish(mesh,err)
    ENDIF

    !Create a decomposition for the mesh
    CALL OC_Decomposition_Initialise(decomposition,err)
    CALL OC_Decomposition_CreateStart(DECOMPOSITION_USER_NUMBER,mesh,decomposition,err)
    CALL OC_Decomposition_TypeSet(decomposition,OC_DECOMPOSITION_CALCULATED_TYPE,err)
    CALL OC_Decomposition_CreateFinish(decomposition,err)

    !Decompose
    CALL OC_Decomposer_Initialise(decomposer,err)
    CALL OC_Decomposer_CreateStart(DECOMPOSER_USER_NUMBER,region,worldWorkGroup,decomposer,err)
    !Add in the decomposition
    CALL OC_Decomposer_DecompositionAdd(decomposer,decomposition,decompositionIndex,err)
    !Finish the decomposer
    CALL OC_Decomposer_CreateFinish(decomposer,err)
    
    !Create a field for the geometry
    CALL OC_Field_Initialise(geometricField,err)
    CALL OC_Field_CreateStart(GEOMETRIC_FIELD_USER_NUMBER,region,geometricField,err)
    CALL OC_Field_DecompositionSet(geometricField,Decomposition,err)
    CALL OC_Field_TypeSet(geometricField,OC_FIELD_GEOMETRIC_TYPE,err)
    CALL OC_Field_VariableLabelSet(geometricField,OC_FIELD_U_VARIABLE_TYPE,"Geometry",err)
    CALL OC_Field_ComponentMeshComponentSet(geometricField,OC_FIELD_U_VARIABLE_TYPE,1,1,err)
    CALL OC_Field_ComponentMeshComponentSet(geometricField,OC_FIELD_U_VARIABLE_TYPE,2,1,err)
    CALL OC_Field_ComponentMeshComponentSet(geometricField,OC_FIELD_U_VARIABLE_TYPE,3,1,err)
    IF(interpolationType==OC_BASIS_CUBIC_HERMITE_INTERPOLATION) THEN
      CALL OC_Field_ScalingTypeSet(geometricField,OC_FIELD_ARITHMETIC_MEAN_SCALING,err)
    END IF
    CALL OC_Field_CreateFinish(geometricField,err)

    IF (useGeneratedMesh) THEN
      ! Update the geometric field parameters from generated mesh
      CALL OC_GeneratedMesh_GeometricParametersCalculate(generatedMesh,geometricField,err)
    ELSE
      ! Update the geometric field parameters manually
      CALL OC_Field_ParameterSetUpdateStart(geometricField, &
        & OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,err)
      ! node 1
      CALL OC_Field_ParameterSetUpdateNode(geometricField,OC_FIELD_U_VARIABLE_TYPE, &
        & OC_FIELD_VALUES_SET_TYPE,1,1,1,1,0.0_OC_RP,err)
      CALL OC_Field_ParameterSetUpdateNode(geometricField,OC_FIELD_U_VARIABLE_TYPE, &
        & OC_FIELD_VALUES_SET_TYPE,1,1,1,2,0.0_OC_RP,err)
      CALL OC_Field_ParameterSetUpdateNode(geometricField,OC_FIELD_U_VARIABLE_TYPE, &
        & OC_FIELD_VALUES_SET_TYPE,1,1,1,3,0.0_OC_RP,err)
      ! node 2
      CALL OC_Field_ParameterSetUpdateNode(geometricField,OC_FIELD_U_VARIABLE_TYPE, &
        & OC_FIELD_VALUES_SET_TYPE,1,1,2,1,HEIGHT,err)
      CALL OC_Field_ParameterSetUpdateNode(geometricField,OC_FIELD_U_VARIABLE_TYPE, &
        & OC_FIELD_VALUES_SET_TYPE,1,1,2,2,0.0_OC_RP,err)
      CALL OC_Field_ParameterSetUpdateNode(geometricField,OC_FIELD_U_VARIABLE_TYPE, &
        & OC_FIELD_VALUES_SET_TYPE,1,1,2,3,0.0_OC_RP,err)
      ! node 3
      CALL OC_Field_ParameterSetUpdateNode(geometricField,OC_FIELD_U_VARIABLE_TYPE, &
        & OC_FIELD_VALUES_SET_TYPE,1,1,3,1,0.0_OC_RP,err)
      CALL OC_Field_ParameterSetUpdateNode(geometricField,OC_FIELD_U_VARIABLE_TYPE, &
        & OC_FIELD_VALUES_SET_TYPE,1,1,3,2,WIDTH,err)
      CALL OC_Field_ParameterSetUpdateNode(geometricField,OC_FIELD_U_VARIABLE_TYPE, &
        & OC_FIELD_VALUES_SET_TYPE,1,1,3,3,0.0_OC_RP,err)
      ! node 4
      CALL OC_Field_ParameterSetUpdateNode(geometricField,OC_FIELD_U_VARIABLE_TYPE, &
        & OC_FIELD_VALUES_SET_TYPE,1,1,4,1,HEIGHT,err)
      CALL OC_Field_ParameterSetUpdateNode(geometricField,OC_FIELD_U_VARIABLE_TYPE, &
        & OC_FIELD_VALUES_SET_TYPE,1,1,4,2,WIDTH,err)
      CALL OC_Field_ParameterSetUpdateNode(geometricField,OC_FIELD_U_VARIABLE_TYPE, &
        & OC_FIELD_VALUES_SET_TYPE,1,1,4,3,0.0_OC_RP,err)
      ! node 5
      CALL OC_Field_ParameterSetUpdateNode(geometricField,OC_FIELD_U_VARIABLE_TYPE, &
        & OC_FIELD_VALUES_SET_TYPE,1,1,5,1,0.0_OC_RP,err)
      CALL OC_Field_ParameterSetUpdateNode(geometricField,OC_FIELD_U_VARIABLE_TYPE, &
        & OC_FIELD_VALUES_SET_TYPE,1,1,5,2,0.0_OC_RP,err)
      CALL OC_Field_ParameterSetUpdateNode(geometricField,OC_FIELD_U_VARIABLE_TYPE, &
        & OC_FIELD_VALUES_SET_TYPE,1,1,5,3,LENGTH,err)
      ! node 6
      CALL OC_Field_ParameterSetUpdateNode(geometricField,OC_FIELD_U_VARIABLE_TYPE, &
        & OC_FIELD_VALUES_SET_TYPE,1,1,6,1,HEIGHT,err)
      CALL OC_Field_ParameterSetUpdateNode(geometricField,OC_FIELD_U_VARIABLE_TYPE, &
        & OC_FIELD_VALUES_SET_TYPE,1,1,6,2,0.0_OC_RP,err)
      CALL OC_Field_ParameterSetUpdateNode(geometricField,OC_FIELD_U_VARIABLE_TYPE, &
        & OC_FIELD_VALUES_SET_TYPE,1,1,6,3,LENGTH,err)
      ! node 7
      CALL OC_Field_ParameterSetUpdateNode(geometricField,OC_FIELD_U_VARIABLE_TYPE, &
        & OC_FIELD_VALUES_SET_TYPE,1,1,7,1,0.0_OC_RP,err)
      CALL OC_Field_ParameterSetUpdateNode(geometricField,OC_FIELD_U_VARIABLE_TYPE, &
        & OC_FIELD_VALUES_SET_TYPE,1,1,7,2,WIDTH,err)
      CALL OC_Field_ParameterSetUpdateNode(geometricField,OC_FIELD_U_VARIABLE_TYPE, &
        & OC_FIELD_VALUES_SET_TYPE,1,1,7,3,LENGTH,err)
      ! node 8
      CALL OC_Field_ParameterSetUpdateNode(geometricField,OC_FIELD_U_VARIABLE_TYPE, &
        & OC_FIELD_VALUES_SET_TYPE,1,1,8,1,HEIGHT,err)
      CALL OC_Field_ParameterSetUpdateNode(geometricField,OC_FIELD_U_VARIABLE_TYPE, &
        & OC_FIELD_VALUES_SET_TYPE,1,1,8,2,WIDTH,err)
      CALL OC_Field_ParameterSetUpdateNode(geometricField,OC_FIELD_U_VARIABLE_TYPE, &
        & OC_FIELD_VALUES_SET_TYPE,1,1,8,3,LENGTH,err)
      CALL OC_Field_ParameterSetUpdateFinish(geometricField, &
        & OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,err)
    ENDIF

    !Create a fibre field and attach it to the geometric field
    CALL OC_Field_Initialise(fibreField,err)
    CALL OC_Field_CreateStart(FIBRE_FIELD_USER_NUMBER,region,fibreField,err)
    CALL OC_Field_TypeSet(fibreField,OC_FIELD_FIBRE_TYPE,err)
    CALL OC_Field_DecompositionSet(fibreField,decomposition,err)
    CALL OC_Field_GeometricFieldSet(fibreField,geometricField,err)
    CALL OC_Field_VariableLabelSet(fibreField,OC_FIELD_U_VARIABLE_TYPE,"Fibre",err)
    IF(interpolationType==OC_BASIS_CUBIC_HERMITE_INTERPOLATION) THEN
      CALL OC_Field_ScalingTypeSet(fibreField,OC_FIELD_ARITHMETIC_MEAN_SCALING,err)
    ENDIF
    CALL OC_Field_CreateFinish(fibreField,err)

    !Create the material field
    IF(compressible) THEN
      numberOfMaterialComponents = 3
    ELSE
      numberOfMaterialComponents = 2
    ENDIF
    CALL OC_Field_Initialise(materialField,err)
    CALL OC_Field_CreateStart(MATERIAL_FIELD_USER_NUMBER,Region,materialField,err)
    CALL OC_Field_TypeSet(materialField,OC_FIELD_MATERIAL_TYPE,err)
    CALL OC_Field_DecompositionSet(materialField,decomposition,err)
    CALL OC_Field_GeometricFieldSet(materialField,geometricField,err)
    CALL OC_Field_NumberOfVariablesSet(materialField,1,err)
    CALL OC_Field_NumberOfComponentsSet(materialField,OC_FIELD_U_VARIABLE_TYPE,numberOfMaterialComponents,err)
    CALL OC_Field_VariableLabelSet(materialField,OC_FIELD_U_VARIABLE_TYPE,"Material",err)
    CALL OC_Field_ComponentMeshComponentSet(materialField,OC_FIELD_U_VARIABLE_TYPE,1,1,err)
    CALL OC_Field_ComponentMeshComponentSet(materialField,OC_FIELD_U_VARIABLE_TYPE,2,1,err)
    IF(compressible) THEN
      CALL OC_Field_ComponentMeshComponentSet(materialField,OC_FIELD_U_VARIABLE_TYPE,3,1,err)
    ENDIF
    IF(interpolationType==OC_BASIS_CUBIC_HERMITE_INTERPOLATION) THEN
      CALL OC_Field_ScalingTypeSet(materialField,OC_FIELD_ARITHMETIC_MEAN_SCALING,err)
    ENDIF
    CALL OC_Field_CreateFinish(materialField,err)

    !Set Mooney-Rivlin constants c10 and c01 respectively.
    CALL OC_Field_ComponentValuesInitialise(materialField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
      & 1,2.0_OC_RP,err)
    CALL OC_Field_ComponentValuesInitialise(materialField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
      & 2,6.0_OC_RP,err)
    IF(compressible) THEN
      CALL OC_Field_ComponentValuesInitialise(materialField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
        & 3,1.0e9_OC_RP,err)
    ENDIF

    !Create the dependent field
    IF(compressible) THEN
      numberOfMaterialComponents = 3
    ELSE
      numberOfMaterialComponents = 4
    ENDIF

    CALL OC_Field_Initialise(dependentField,err)
    CALL OC_Field_CreateStart(DEPENDENT_FIELD_USER_NUMBER,region,dependentField,err)
    CALL OC_Field_VariableLabelSet(dependentField,OC_FIELD_U_VARIABLE_TYPE,"Dependent",err)
    CALL OC_Field_TypeSet(dependentField,OC_FIELD_GEOMETRIC_GENERAL_TYPE,err)
    CALL OC_Field_DecompositionSet(dependentField,decomposition,err)
    CALL OC_Field_GeometricFieldSet(dependentField,geometricField,err)
    CALL OC_Field_DependentTypeSet(dependentField,OC_FIELD_DEPENDENT_TYPE,err)
    CALL OC_Field_NumberOfVariablesSet(dependentField,2,err)
    CALL OC_Field_NumberOfComponentsSet(dependentField,OC_FIELD_U_VARIABLE_TYPE,numberOfMaterialComponents,err)
    CALL OC_Field_NumberOfComponentsSet(dependentField,OC_FIELD_T_VARIABLE_TYPE,numberOfMaterialComponents,err)
    CALL OC_Field_ComponentMeshComponentSet(dependentField,OC_FIELD_U_VARIABLE_TYPE,1,1,err)
    CALL OC_Field_ComponentMeshComponentSet(dependentField,OC_FIELD_U_VARIABLE_TYPE,2,1,err)
    CALL OC_Field_ComponentMeshComponentSet(dependentField,OC_FIELD_U_VARIABLE_TYPE,3,1,err)
    CALL OC_Field_ComponentMeshComponentSet(dependentField,OC_FIELD_T_VARIABLE_TYPE,1,1,err)
    CALL OC_Field_ComponentMeshComponentSet(dependentField,OC_FIELD_T_VARIABLE_TYPE,2,1,err)
    CALL OC_Field_ComponentMeshComponentSet(dependentField,OC_FIELD_T_VARIABLE_TYPE,3,1,err)
    IF(.NOT.compressible) THEN
      ! TODO we always had node based interpolation, right? --> check!
      CALL OC_Field_ComponentInterpolationSet(dependentField,OC_FIELD_U_VARIABLE_TYPE, &
        & 4,OC_FIELD_ELEMENT_BASED_INTERPOLATION,err)
      CALL OC_Field_ComponentInterpolationSet(dependentField,OC_FIELD_T_VARIABLE_TYPE, &
        & 4,OC_FIELD_ELEMENT_BASED_INTERPOLATION,err)
      ! TODO end
      IF(usePressureBasis) THEN
        !Set the pressure to be nodally based and use the second mesh component
        IF(interpolationType==4) THEN
          CALL OC_Field_ComponentInterpolationSet(dependentField,OC_FIELD_U_VARIABLE_TYPE,4, &
            & OC_FIELD_NODE_BASED_INTERPOLATION,err)
          CALL OC_Field_ComponentInterpolationSet(dependentField,OC_FIELD_T_VARIABLE_TYPE,4, &
            & OC_FIELD_NODE_BASED_INTERPOLATION,err)
        ENDIF
        CALL OC_Field_ComponentInterpolationSet(dependentField,OC_FIELD_U_VARIABLE_TYPE,4,2,err)
        CALL OC_Field_ComponentInterpolationSet(dependentField,OC_FIELD_T_VARIABLE_TYPE,4,2,err)
      ENDIF
    ENDIF
    IF(interpolationType==OC_BASIS_CUBIC_HERMITE_INTERPOLATION) THEN
      CALL OC_Field_ScalingTypeSet(dependentField,OC_FIELD_ARITHMETIC_MEAN_SCALING,err)
    ENDIF
    CALL OC_Field_CreateFinish(dependentField,err)

    !Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
    CALL OC_Field_ParametersToFieldParametersComponentCopy( &
      & geometricField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1, &
      & dependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,err)
    CALL OC_Field_ParametersToFieldParametersComponentCopy( &
      & geometricField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,2, &
      & dependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,2,err)
    CALL OC_Field_ParametersToFieldParametersComponentCopy( &
      & geometricField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,3, &
      & dependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,3,err)
    IF(.NOT.compressible) THEN
      CALL OC_Field_ComponentValuesInitialise(dependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
        & 4,0.0_OC_RP,err)
    ENDIF

    !Create a deformed geometry field, as cmgui doesn't like displaying
    ! deformed fibres from the dependent field because it isn't a geometric field.
    CALL OC_Field_Initialise(deformedField,err)
    CALL OC_Field_CreateStart(DEFORMED_FIELD_USER_NUMBER,region,deformedField,err)
    CALL OC_Field_DecompositionSet(deformedField,decomposition,err)
    CALL OC_Field_TypeSet(deformedField,OC_FIELD_GEOMETRIC_TYPE,err)
    CALL OC_Field_VariableLabelSet(deformedField,OC_FIELD_U_VARIABLE_TYPE,"DeformedGeometry",err)
    DO componentIdx=1,3
      CALL OC_Field_ComponentMeshComponentSet(deformedField,OC_FIELD_U_VARIABLE_TYPE,componentIdx,1,err)
    ENDDO
    IF(interpolationType==OC_BASIS_CUBIC_HERMITE_INTERPOLATION) THEN
      CALL OC_Field_ScalingTypeSet(deformedField,OC_FIELD_ARITHMETIC_MEAN_SCALING,err)
    ENDIF
    CALL OC_Field_CreateFinish(deformedField,err)

    !Create the equations set
    CALL OC_Field_Initialise(equationsSetField,err)
    CALL OC_EquationsSet_Initialise(equationsSet,err)
    IF(compressible) THEN
      CALL OC_EquationsSet_CreateStart(EQUATIONS_SET_USER_NUMBER,region,fibreField, &
        & [OC_EQUATIONS_SET_ELASTICITY_CLASS, &
        &  OC_EQUATIONS_SET_FINITE_ELASTICITY_TYPE, &
        &  OC_EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE], &
        & EQUATIONS_SET_FIELD_USER_NUMBER,equationsSetField,equationsSet,err)
    ELSE
      CALL OC_EquationsSet_CreateStart(EQUATIONS_SET_USER_NUMBER,region,fibreField, &
        & [OC_EQUATIONS_SET_ELASTICITY_CLASS, &
        &  OC_EQUATIONS_SET_FINITE_ELASTICITY_TYPE, &
        &  OC_EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE], &
        & EQUATIONS_SET_FIELD_USER_NUMBER,equationsSetField,equationsSet,err)
    ENDIF
    CALL OC_EquationsSet_CreateFinish(equationsSet,err)
    CALL OC_EquationsSet_MaterialsCreateStart(equationsSet,MATERIAL_FIELD_USER_NUMBER,materialField,err)
    CALL OC_EquationsSet_MaterialsCreateFinish(equationsSet,err)
    CALL OC_EquationsSet_DependentCreateStart(equationsSet,DEPENDENT_FIELD_USER_NUMBER,dependentField,err)
    CALL OC_EquationsSet_DependentCreateFinish(equationsSet,err)

    !Create equations
    CALL OC_Equations_Initialise(equations,err)
    CALL OC_EquationsSet_EquationsCreateStart(equationsSet,equations,err)
    CALL OC_Equations_SparsityTypeSet(equations,OC_EQUATIONS_SPARSE_MATRICES,err)
    CALL OC_Equations_OutputTypeSet(equations,OC_EQUATIONS_NO_OUTPUT,err)
    CALL OC_EquationsSet_EquationsCreateFinish(equationsSet,err)

    CALL OC_Equations_JacobianCalculationTypeSet(equations,1,OC_FIELD_U_VARIABLE_TYPE, &
      & OC_EQUATIONS_JACOBIAN_ANALYTIC_CALCULATED,err)

    !Define the problem
    CALL OC_Problem_Initialise(problem,err)
    CALL OC_Problem_CreateStart(PROBLEM_USER_NUMBER,context, &
      & [OC_PROBLEM_ELASTICITY_CLASS, &
      &  OC_PROBLEM_FINITE_ELASTICITY_TYPE, &
      &  OC_PROBLEM_STATIC_FINITE_ELASTICITY_SUBTYPE],problem,err)
    CALL OC_Problem_CreateFinish(problem,err)

    !Create control loops
    CALL OC_Problem_ControlLoopCreateStart(problem,err)
    CALL OC_Problem_ControlLoopCreateFinish(problem,err)

    !Create problem solver
    CALL OC_Solver_Initialise(nonlinearSolver,err)
    CALL OC_Solver_Initialise(linearSolver,err)
    CALL OC_Problem_SolversCreateStart(problem,err)
    CALL OC_Problem_SolverGet(problem,OC_CONTROL_LOOP_NODE,1,nonLinearSolver,err)
    CALL OC_Solver_OutputTypeSet(nonlinearSolver,OC_SOLVER_PROGRESS_OUTPUT,err)
    !CALL OC_Solver_NewtonJacobianCalculationTypeSet(nonlinearSolver, &
    !  & OC_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,err)
    CALL OC_Solver_NewtonJacobianCalculationTypeSet(nonlinearSolver, &
      & OC_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED,err)
    CALL OC_Solver_NewtonLinearSolverGet(nonlinearSolver,linearSolver,err)
    CALL OC_Solver_NewtonAbsoluteToleranceSet(nonlinearSolver,1.0E-14_OC_RP,err)
    CALL OC_Solver_NewtonSolutionToleranceSet(nonlinearSolver,1.0E-14_OC_RP,err)
    CALL OC_Solver_NewtonRelativeToleranceSet(nonlinearSolver,1.0E-14_OC_RP,err)
    CALL OC_Solver_LinearTypeSet(linearSolver,OC_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,err)
    CALL OC_Problem_SolversCreateFinish(problem,err)

    !Create solver equations and add equations set to solver equations
    CALL OC_Solver_Initialise(solver,err)
    CALL OC_SolverEquations_Initialise(solverEquations,err)
    CALL OC_Problem_SolverEquationsCreateStart(problem,err)
    CALL OC_Problem_SolverGet(problem,OC_CONTROL_LOOP_NODE,1,solver,err)
    CALL OC_Solver_SolverEquationsGet(solver,solverEquations,err)
    CALL OC_SolverEquations_SparsityTypeSet(solverEquations,OC_SOLVER_SPARSE_MATRICES,err)
    CALL OC_SolverEquations_EquationsSetAdd(solverEquations,equationsSet,equationsSetIndex,err)
    CALL OC_Problem_SolverEquationsCreateFinish(problem,err)

    !Prescribe boundary conditions (absolute nodal parameters)
    CALL OC_BoundaryConditions_Initialise(boundaryConditions,err)
    CALL OC_SolverEquations_BoundaryConditionsCreateStart(solverEquations,boundaryConditions,err)

    !Set x=0 nodes to no x displacment in x. Set x=WIDTH nodes to 10% x displacement
    CALL OC_BoundaryConditions_AddNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE, &
      & 1,1,1,1,OC_BOUNDARY_CONDITION_FIXED,0.0_OC_RP,err)
    CALL OC_BoundaryConditions_AddNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE, &
      & 1,1,3,1,OC_BOUNDARY_CONDITION_FIXED,0.0_OC_RP,err)
    CALL OC_BoundaryConditions_AddNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE, &
      & 1,1,5,1,OC_BOUNDARY_CONDITION_FIXED,0.0_OC_RP,err)
    CALL OC_BoundaryConditions_AddNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE, &
      & 1,1,7,1,OC_BOUNDARY_CONDITION_FIXED,0.0_OC_RP,err)
    IF (zeroLoad) THEN
      load = 0.0_OC_RP
    ELSE
      load = 0.1_OC_RP*WIDTH
    END IF
    CALL OC_BoundaryConditions_AddNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE, &
      & 1,1,2,1,OC_BOUNDARY_CONDITION_FIXED,load,err)
    CALL OC_BoundaryConditions_AddNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE, &
      & 1,1,4,1,OC_BOUNDARY_CONDITION_FIXED,load,err)
    CALL OC_BoundaryConditions_AddNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE, &
      & 1,1,6,1,OC_BOUNDARY_CONDITION_FIXED,load,err)
    CALL OC_BoundaryConditions_AddNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE, &
      & 1,1,8,1,OC_BOUNDARY_CONDITION_FIXED,load,err)

    ! Set y=0 nodes to no y displacement
    CALL OC_BoundaryConditions_AddNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE, &
      & 1,1,1,2,OC_BOUNDARY_CONDITION_FIXED,0.0_OC_RP,err)
    CALL OC_BoundaryConditions_AddNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE, &
      & 1,1,2,2,OC_BOUNDARY_CONDITION_FIXED,0.0_OC_RP,err)
    CALL OC_BoundaryConditions_AddNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE, &
      & 1,1,5,2,OC_BOUNDARY_CONDITION_FIXED,0.0_OC_RP,err)
    CALL OC_BoundaryConditions_AddNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE, &
      & 1,1,6,2,OC_BOUNDARY_CONDITION_FIXED,0.0_OC_RP,err)

    ! Set z=0 nodes to no y displacement
    CALL OC_BoundaryConditions_AddNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE, &
      & 1,1,1,3,OC_BOUNDARY_CONDITION_FIXED,0.0_OC_RP,err)
    CALL OC_BoundaryConditions_AddNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE, &
      & 1,1,2,3,OC_BOUNDARY_CONDITION_FIXED,0.0_OC_RP,err)
    CALL OC_BoundaryConditions_AddNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE, &
      & 1,1,3,3,OC_BOUNDARY_CONDITION_FIXED,0.0_OC_RP,err)
    CALL OC_BoundaryConditions_AddNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE, &
      & 1,1,4,3,OC_BOUNDARY_CONDITION_FIXED,0.0_OC_RP,err)

    CALL OC_SolverEquations_BoundaryConditionsCreateFinish(solverEquations,err)

    !Solve the problem
    CALL OC_Problem_Solve(problem,err)

    !Copy deformed geometry into deformed field
    DO componentIdx=1,3
      CALL OC_Field_ParametersToFieldParametersComponentCopy( &
        & dependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,componentIdx, &
        & deformedField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,componentIdx,err)
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
    CALL OC_Fields_Initialise(fields,err)
    CALL OC_Fields_Create(region,fields,err)
    CALL OC_Fields_NodesExport(fields,trim(outputFile)//trim(suffix),"FORTRAN",err)
    CALL OC_Fields_ElementsExport(fields,trim(outputFile)//trim(suffix),"FORTRAN",err)
    CALL OC_Fields_Finalise(fields,err)

    CALL OC_Problem_Destroy(problem,err)
    IF(usePressureBasis) CALL OC_Basis_Destroy(pressureBasis,err)
    CALL OC_Basis_Destroy(basis,err)
    CALL OC_Region_Destroy(region,err)
    CALL OC_CoordinateSystem_Destroy(coordinateSystem,err)
   
    WRITE(*,'(A)') "Program successfully completed."

  END SUBROUTINE SolveModel

  SUBROUTINE HandelError(errorString)
    
    CHARACTER(LEN=*), INTENT(IN) :: errorString

    WRITE(*,'(">>ERROR: ",A)') errorString(1:LEN_TRIM(errorString))
    STOP
    
  END SUBROUTINE HandelError

END PROGRAM UniaxialExtension
