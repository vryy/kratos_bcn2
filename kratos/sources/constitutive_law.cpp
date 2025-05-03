//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Janosch Stascheit
//                   Nelson Maireni Lafontaine
//                   Josep Maria Carbonell
//

#include "includes/constitutive_law.h"


namespace Kratos
{

const unsigned int BaseConstitutiveLaw::msIndexVoigt3D6C [6][2] = { {0, 0}, {1, 1}, {2, 2}, {0, 1}, {1, 2}, {0, 2} };
const unsigned int BaseConstitutiveLaw::msIndexVoigt2D4C [4][2] = { {0, 0}, {1, 1}, {2, 2}, {0, 1} };
const unsigned int BaseConstitutiveLaw::msIndexVoigt2D3C [3][2] = { {0, 0}, {1, 1}, {0, 1} };

/**
 * Flags related to the Parameters of the Contitutive Law
 */
KRATOS_CREATE_LOCAL_FLAG( BaseConstitutiveLaw, COMPUTE_STRAIN,              0 );
KRATOS_CREATE_LOCAL_FLAG( BaseConstitutiveLaw, COMPUTE_STRESS,              1 );
KRATOS_CREATE_LOCAL_FLAG( BaseConstitutiveLaw, COMPUTE_CONSTITUTIVE_TENSOR, 2 );
KRATOS_CREATE_LOCAL_FLAG( BaseConstitutiveLaw, COMPUTE_STRAIN_ENERGY,       3 );

KRATOS_CREATE_LOCAL_FLAG( BaseConstitutiveLaw, ISOCHORIC_TENSOR_ONLY,       4 );
KRATOS_CREATE_LOCAL_FLAG( BaseConstitutiveLaw, VOLUMETRIC_TENSOR_ONLY,      5 );

KRATOS_CREATE_LOCAL_FLAG( BaseConstitutiveLaw, TOTAL_TENSOR,                6 );

KRATOS_CREATE_LOCAL_FLAG( BaseConstitutiveLaw, FINALIZE_MATERIAL_RESPONSE,  7 );

/**
 * Flags related to the Features of the Contitutive Law
 */
KRATOS_CREATE_LOCAL_FLAG( BaseConstitutiveLaw, FINITE_STRAINS,              8 );
KRATOS_CREATE_LOCAL_FLAG( BaseConstitutiveLaw, INFINITESIMAL_STRAINS,       9 );

KRATOS_CREATE_LOCAL_FLAG( BaseConstitutiveLaw, THREE_DIMENSIONAL_LAW,      10 );
KRATOS_CREATE_LOCAL_FLAG( BaseConstitutiveLaw, PLANE_STRAIN_LAW,           11 );
KRATOS_CREATE_LOCAL_FLAG( BaseConstitutiveLaw, PLANE_STRESS_LAW,           12 );
KRATOS_CREATE_LOCAL_FLAG( BaseConstitutiveLaw, AXISYMMETRIC_LAW,           13 );

KRATOS_CREATE_LOCAL_FLAG( BaseConstitutiveLaw, U_P_LAW,                    14 );
KRATOS_CREATE_LOCAL_FLAG( BaseConstitutiveLaw, ISOTROPIC,                  15 );
KRATOS_CREATE_LOCAL_FLAG( BaseConstitutiveLaw, ANISOTROPIC,                16 );

/**
 * Constructor.
 */
template<class TNodeType>
ConstitutiveLawImpl<TNodeType>::ConstitutiveLawImpl() : BaseConstitutiveLaw()
{
}

/**
 * Clone function (has to be implemented by any derived class)
 * @return a pointer to a new instance of this constitutive law
 * NOTE: implementation scheme:
 *      ConstitutiveLawImpl<TNodeType>::Pointer p_clone(new ConstitutiveLaw());
 *      return p_clone;
 */
template<class TNodeType>
typename ConstitutiveLawImpl<TNodeType>::Pointer ConstitutiveLawImpl<TNodeType>::Clone() const
{
    KRATOS_ERROR << "Called the virtual function for Clone";
}

/**
 * @brief It creates a new constitutive law pointer
 * @param NewParameters The configuration parameters of the new constitutive law
 * @return a Pointer to the new constitutive law
 */
template<class TNodeType>
typename ConstitutiveLawImpl<TNodeType>::Pointer ConstitutiveLawImpl<TNodeType>::Create(Kratos::Parameters NewParameters) const
{
    KRATOS_ERROR << "Called the virtual function for Clone";
}

/**
 * @brief It creates a new constitutive law pointer (version with properties)
 * @param NewParameters The configuration parameters of the new constitutive law
 * @param rProperties The properties of the material
 * @return a Pointer to the new constitutive law
 */
template<class TNodeType>
typename ConstitutiveLawImpl<TNodeType>::Pointer ConstitutiveLawImpl<TNodeType>::Create(
    Kratos::Parameters NewParameters,
    const Properties& rProperties
    ) const
{
    return this->Create(NewParameters);
}

/**
 * @return the working space dimension of the current constitutive law
 * NOTE: this function HAS TO BE IMPLEMENTED by any derived class
 */
template<class TNodeType>
typename ConstitutiveLawImpl<TNodeType>::SizeType ConstitutiveLawImpl<TNodeType>::WorkingSpaceDimension()
{
    KRATOS_ERROR << "Called the virtual function for WorkingSpaceDimension";
}

/**
 * returns the size of the strain vector of the current constitutive law
 * NOTE: this function HAS TO BE IMPLEMENTED by any derived class
 */
template<class TNodeType>
typename ConstitutiveLawImpl<TNodeType>::SizeType ConstitutiveLawImpl<TNodeType>::GetStrainSize() const
{
    KRATOS_ERROR << "Called the virtual function for GetStrainSize";
}

/**
 * returns whether this constitutive Law has specified variable
 * @param rThisVariable the variable to be checked for
 * @return true if the variable is defined in the constitutive law
 */
template<class TNodeType>
bool ConstitutiveLawImpl<TNodeType>::Has(const Variable<bool>& rThisVariable)
{
    return false;
}

/**
 * returns whether this constitutive Law has specified variable
 * @param rThisVariable the variable to be checked for
 * @return true if the variable is defined in the constitutive law
 */
template<class TNodeType>
bool ConstitutiveLawImpl<TNodeType>::Has(const Variable<int>& rThisVariable)
{
    return false;
}

/**
 * returns whether this constitutive Law has specified variable
 * @param rThisVariable the variable to be checked for
 * @return true if the variable is defined in the constitutive law
 */
template<class TNodeType>
bool ConstitutiveLawImpl<TNodeType>::Has(const Variable<DataType>& rThisVariable)
{
    return false;
}

/**
 * returns whether this constitutive Law has specified variable
 * @param rThisVariable the variable to be checked for
 * @return true if the variable is defined in the constitutive law
 */
template<class TNodeType>
bool ConstitutiveLawImpl<TNodeType>::Has(const Variable<VectorType>& rThisVariable)
{
    return false;
}

/**
 * returns whether this constitutive Law has specified variable
 * @param rThisVariable the variable to be checked for
 * @return true if the variable is defined in the constitutive law
 */
template<class TNodeType>
bool ConstitutiveLawImpl<TNodeType>::Has(const Variable<MatrixType>& rThisVariable)
{
    return false;
}

/**
 * returns whether this constitutive Law has specified variable
 * @param rThisVariable the variable to be checked for
 * @return true if the variable is defined in the constitutive law
 * NOTE: fixed size array of 3 DataTypes (e.g. for 2D stresses, plastic strains, ...)
 */
template<class TNodeType>
bool ConstitutiveLawImpl<TNodeType>::Has(const Variable<array_1d<DataType, 3> >& rThisVariable)
{
    return false;
}

/**
 * returns whether this constitutive Law has specified variable
 * @param rThisVariable the variable to be checked for
 * @return true if the variable is defined in the constitutive law
 * NOTE: fixed size array of 6 DataTypes (e.g. for stresses, plastic strains, ...)
 */
template<class TNodeType>
bool ConstitutiveLawImpl<TNodeType>::Has(const Variable<array_1d<DataType, 6> >& rThisVariable)
{
    return false;
}

/**
 * returns the value of a specified variable
 * @param rThisVariable the variable to be returned
 * @param rValue a reference to the returned value
 * @param rValue output: the value of the specified variable
 */
template<class TNodeType>
bool& ConstitutiveLawImpl<TNodeType>::GetValue(const Variable<bool>& rThisVariable, bool& rValue)
{
    return rValue;
}

/**
 * returns the value of a specified variable
 * @param rThisVariable the variable to be returned
 * @param rValue a reference to the returned value
 * @param rValue output: the value of the specified variable
 */
template<class TNodeType>
int& ConstitutiveLawImpl<TNodeType>::GetValue(const Variable<int>& rThisVariable, int& rValue)
{
    return rValue;
}

/**
 * returns the value of a specified variable
 * @param rThisVariable the variable to be returned
 * @param rValue a reference to the returned value
 * @param rValue output: the value of the specified variable
 */
template<class TNodeType>
typename ConstitutiveLawImpl<TNodeType>::DataType& ConstitutiveLawImpl<TNodeType>::GetValue(const Variable<DataType>& rThisVariable, DataType& rValue)
{
    return rValue;
}

/**
 * returns the value of a specified variable
 * @param rThisVariable the variable to be returned
 * @param rValue a reference to the returned value
 * @return the value of the specified variable
 */
template<class TNodeType>
typename ConstitutiveLawImpl<TNodeType>::VectorType& ConstitutiveLawImpl<TNodeType>::GetValue(const Variable<VectorType>& rThisVariable, VectorType& rValue)
{
    return rValue;
}

/**
 * returns the value of a specified variable
 * @param rThisVariable the variable to be returned
 * @return the value of the specified variable
 */
template<class TNodeType>
typename ConstitutiveLawImpl<TNodeType>::MatrixType& ConstitutiveLawImpl<TNodeType>::GetValue(const Variable<MatrixType>& rThisVariable, MatrixType& rValue)
{
    return rValue;
}

/**
 * returns the value of a specified variable
 * @param rThisVariable the variable to be returned
 * @return the value of the specified variable
 */
template<class TNodeType>
std::string& ConstitutiveLawImpl<TNodeType>::GetValue(const Variable<std::string>& rThisVariable, std::string& rValue)
{
    return rValue;
}

/**
 * returns the value of a specified variable
 * @param rThisVariable the variable to be returned
 * @param rValue a reference to the returned value
 * @return the value of the specified variable
 */
template<class TNodeType>
array_1d<typename ConstitutiveLawImpl<TNodeType>::DataType, 3 >& ConstitutiveLawImpl<TNodeType>::GetValue(const Variable<array_1d<DataType, 3> >& rThisVariable,
        array_1d<DataType, 3> & rValue)
{
    return rValue;
}

/**
 * returns the value of a specified variable
 * @param rThisVariable the variable to be returned
 * @param rValue a reference to the returned value
 * @return the value of the specified variable
 */
template<class TNodeType>
array_1d<typename ConstitutiveLawImpl<TNodeType>::DataType, 6> & ConstitutiveLawImpl<TNodeType>::GetValue(const Variable<array_1d<DataType, 6> >& rThisVariable,
        array_1d<DataType, 6> & rValue)
{
    return rValue;
}

/**
 * @brief Sets the value of a specified variable (bool)
 * @param rThisVariable the variable to be returned
 * @param Value new value of the specified variable
 * @param rCurrentProcessInfo the process info
 */
template<class TNodeType>
void ConstitutiveLawImpl<TNodeType>::SetValue(const Variable<bool>& rThisVariable,
                               const bool& Value,
                               const ProcessInfo& rCurrentProcessInfo)
{}

/**
 * @brief Sets the value of a specified variable (int)
 * @param rThisVariable the variable to be returned
 * @param Value new value of the specified variable
 * @param rCurrentProcessInfo the process info
 */
template<class TNodeType>
void ConstitutiveLawImpl<TNodeType>::SetValue(const Variable<int>& rThisVariable,
                               const int& Value,
                               const ProcessInfo& rCurrentProcessInfo)
{}

/**
 * @brief Sets the value of a specified variable (DataType)
 * @param rVariable the variable to be returned
 * @param rValue new value of the specified variable
 * @param rCurrentProcessInfo the process info
 */
template<class TNodeType>
void ConstitutiveLawImpl<TNodeType>::SetValue(const Variable<DataType>& rVariable,
                               const DataType& rValue,
                               const ProcessInfo& rCurrentProcessInfo)
{}

/**
 * @brief Sets the value of a specified variable (Vector)
 * @param rVariable the variable to be returned
 * @param rValue new value of the specified variable
 * @param rCurrentProcessInfo the process info
 */
template<class TNodeType>
void ConstitutiveLawImpl<TNodeType>::SetValue(const Variable<VectorType >& rVariable,
                               const VectorType& rValue, const ProcessInfo& rCurrentProcessInfo)
{}

/**
 * @brief Sets the value of a specified variable (Matrix)
 * @param rVariable the variable to be returned
 * @param rValue new value of the specified variable
 * @param rCurrentProcessInfo the process info
 */
template<class TNodeType>
void ConstitutiveLawImpl<TNodeType>::SetValue(const Variable<MatrixType >& rVariable,
                               const MatrixType& rValue, const ProcessInfo& rCurrentProcessInfo)
{}

/**
 * @brief Sets the value of a specified variable (array of 3 components)
 * @param rVariable the variable to be returned
 * @param rValue new value of the specified variable
 * @param rCurrentProcessInfo the process info
 */
template<class TNodeType>
void ConstitutiveLawImpl<TNodeType>::SetValue(const Variable<array_1d<DataType, 3> >& rVariable,
                               const array_1d<DataType, 3> & rValue,
                               const ProcessInfo& rCurrentProcessInfo)
{}

/**
 * @brief Sets the value of a specified variable (array of 6 components)
 * @param rVariable the variable to be returned
 * @param rValue new value of the specified variable
 * @param rCurrentProcessInfo the process info
 */
template<class TNodeType>
void ConstitutiveLawImpl<TNodeType>::SetValue(const Variable<array_1d<DataType, 6> >& rVariable,
                               const array_1d<DataType, 6> & rValue,
                               const ProcessInfo& rCurrentProcessInfo)
{}

/**
 * @brief Sets the value of a specified variable (std::string components)
 * @param rVariable the variable to be returned
 * @param rValue new value of the specified variable
 * @param rCurrentProcessInfo the process info
 */
template<class TNodeType>
void ConstitutiveLawImpl<TNodeType>::SetValue(const Variable<std::string >& rVariable,
                               const std::string& rValue, const ProcessInfo& rCurrentProcessInfo)
{}

/**
 * @brief Sets the value of a specified variable (ConstitutiveLaw component)
 * @param rVariable the variable to be returned
 * @param rValue new value of the specified variable
 * @param rCurrentProcessInfo the process info
 */
template<class TNodeType>
void ConstitutiveLawImpl<TNodeType>::SetValue(const Variable<ConstitutiveLawImpl<TNodeType>::Pointer>& rVariable,
                               ConstitutiveLawImpl<TNodeType>::Pointer rValue,
                               const ProcessInfo& rCurrentProcessInfo)
{}

/**
 * @brief Calculates the value of a specified variable (bool)
 * @param rParameterValues the needed parameters for the CL calculation
 * @param rThisVariable the variable to be returned
 * @param rValue a reference to the returned value
 * @param rValue output: the value of the specified variable
 */
template<class TNodeType>
bool& ConstitutiveLawImpl<TNodeType>::CalculateValue(Parameters& rParameterValues, const Variable<bool>& rThisVariable, bool& rValue)
{
    return rValue;
}

/**
 * @brief Calculates the value of a specified variable (int)
 * @param rParameterValues the needed parameters for the CL calculation
 * @param rThisVariable the variable to be returned
 * @param rValue a reference to the returned value
 * @param rValue output: the value of the specified variable
 */
template<class TNodeType>
int& ConstitutiveLawImpl<TNodeType>::CalculateValue(Parameters& rParameterValues, const Variable<int>& rThisVariable, int& rValue)
{
    return rValue;
}

/**
 * returns the value of a specified variable (DataType)
 * @param rParameterValues the needed parameters for the CL calculation
 * @param rThisVariable the variable to be returned
 * @param rValue a reference to the returned value
 * @param rValue output: the value of the specified variable
 */
template<class TNodeType>
typename ConstitutiveLawImpl<TNodeType>::DataType& ConstitutiveLawImpl<TNodeType>::CalculateValue(Parameters& rParameterValues, const Variable<DataType>& rThisVariable, DataType& rValue)
{
    return rValue;
}

/**
 * returns the value of a specified variable (Vector)
 * @param rParameterValues the needed parameters for the CL calculation
 * @param rThisVariable the variable to be returned
 * @param rValue a reference to the returned value
 * @param rValue output: the value of the specified variable
 */
template<class TNodeType>
typename ConstitutiveLawImpl<TNodeType>::VectorType& ConstitutiveLawImpl<TNodeType>::CalculateValue(Parameters& rParameterValues, const Variable<VectorType>& rThisVariable, VectorType& rValue)
{
    return rValue;
}

/**
 * returns the value of a specified variable (Matrix)
 * @param rParameterValues the needed parameters for the CL calculation
 * @param rThisVariable the variable to be returned
 * @param rValue a reference to the returned value
 * @param rValue output: the value of the specified variable
 */
template<class TNodeType>
typename ConstitutiveLawImpl<TNodeType>::MatrixType& ConstitutiveLawImpl<TNodeType>::CalculateValue(Parameters& rParameterValues, const Variable<MatrixType>& rThisVariable, MatrixType& rValue)
{
    return rValue;
}

/**
 * returns the value of a specified variable (array of 3 components)
 * @param rParameterValues the needed parameters for the CL calculation
 * @param rThisVariable the variable to be returned
 * @param rValue a reference to the returned value
 * @param rValue output: the value of the specified variable
 */
template<class TNodeType>
array_1d<typename ConstitutiveLawImpl<TNodeType>::DataType, 3> & ConstitutiveLawImpl<TNodeType>::CalculateValue(Parameters& rParameterValues, const Variable<array_1d<DataType, 3> >& rVariable,
        array_1d<DataType, 3> & rValue)
{
    return rValue;
}


  /**
 * returns the value of a specified variable (array of 6 components)
 * @param rParameterValues the needed parameters for the CL calculation
 * @param rThisVariable the variable to be returned
 * @param rValue a reference to the returned value
 * @param rValue output: the value of the specified variable
 */
template<class TNodeType>
array_1d<typename ConstitutiveLawImpl<TNodeType>::DataType, 6> & ConstitutiveLawImpl<TNodeType>::CalculateValue(Parameters& rParameterValues, const Variable<array_1d<DataType, 6> >& rVariable,
        array_1d<DataType, 6> & rValue)
{
    return rValue;
}

/**
 * Is called to check whether the provided material parameters in the Properties
 * match the requirements of current constitutive model.
 * @param rMaterialProperties the current Properties to be validated against.
 * @return true, if parameters are correct; false, if parameters are insufficient / faulty
 * NOTE: this has to implemented by each constitutive model. Returns false in base class since
 * no valid implementation is contained here.
 */
template<class TNodeType>
bool ConstitutiveLawImpl<TNodeType>::ValidateInput(const Properties& rMaterialProperties)
{
  return false;
}

/**
 * returns the expected strain measure of this constitutive law (by default linear strains)
 * @return the expected strain measure
 */
template<class TNodeType>
typename ConstitutiveLawImpl<TNodeType>::StrainMeasure ConstitutiveLawImpl<TNodeType>::GetStrainMeasure()
{
    return StrainMeasure_Infinitesimal;
}

/**
 * returns the stress measure of this constitutive law (by default 1st Piola-Kirchhoff stress in voigt notation)
 * @return the expected stress measure
 */
template<class TNodeType>
typename ConstitutiveLawImpl<TNodeType>::StressMeasure ConstitutiveLawImpl<TNodeType>::GetStressMeasure()
{
    return StressMeasure_PK1;
}

/**
 * returns whether this constitutive model is formulated in incremental strains/stresses
 * NOTE: by default, all constitutive models should be formulated in total strains
 * @return true, if formulated in incremental strains/stresses, false otherwise
 */
template<class TNodeType>
bool ConstitutiveLawImpl<TNodeType>::IsIncremental()
{
    return false;
}

/**
 * This is to be called at the very beginning of the calculation
 * (e.g. from InitializeElement) in order to initialize all relevant
 * attributes of the constitutive law
 * @param rMaterialProperties the Properties instance of the current element
 * @param rElementGeometry the geometry of the current element
 * @param rShapeFunctionsValues the shape functions values in the current integration point
 */
template<class TNodeType>
void ConstitutiveLawImpl<TNodeType>::InitializeMaterial(const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues)
{
}

/**
 * to be called at the beginning of each solution step
 * (e.g. from Element::InitializeSolutionStep)
 * @param rMaterialProperties the Properties instance of the current element
 * @param rElementGeometry the geometry of the current element
 * @param rShapeFunctionsValues the shape functions values in the current integration point
 * @param the current ProcessInfo instance
 */
template<class TNodeType>
void ConstitutiveLawImpl<TNodeType>::InitializeSolutionStep(const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry, //this is just to give the array of nodes
        const Vector& rShapeFunctionsValues,
        const ProcessInfo& rCurrentProcessInfo)
{
}

/**
 * to be called at the end of each solution step
 * (e.g. from Element::FinalizeSolutionStep)
 * @param rMaterialProperties the Properties instance of the current element
 * @param rElementGeometry the geometry of the current element
 * @param rShapeFunctionsValues the shape functions values in the current integration point
 * @param the current ProcessInfo instance
 */
template<class TNodeType>
void ConstitutiveLawImpl<TNodeType>::FinalizeSolutionStep(const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues,
        const ProcessInfo& rCurrentProcessInfo)
{
}

/**
 * to be called at the beginning of each step iteration
 * (e.g. from Element::InitializeNonLinearIteration)
 * @param rMaterialProperties the Properties instance of the current element
 * @param rElementGeometry the geometry of the current element
 * @param rShapeFunctionsValues the shape functions values in the current integration point
 * @param the current ProcessInfo instance
 */
template<class TNodeType>
void ConstitutiveLawImpl<TNodeType>::InitializeNonLinearIteration(const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues,
        const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR << "Calling virtual function for InitializeNonLinearIteration";
}

/**
 * to be called at the end of each step iteration
 * (e.g. from Element::FinalizeNonLinearIteration)
 * @param rMaterialProperties the Properties instance of the current element
 * @param rElementGeometry the geometry of the current element
 * @param rShapeFunctionsValues the shape functions values in the current integration point
 * @param the current ProcessInfo instance
 */
template<class TNodeType>
void ConstitutiveLawImpl<TNodeType>::FinalizeNonLinearIteration(const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues,
        const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR << "Calling virtual function for FinalizeNonLinearIteration";
}

/**
 * Computes the material response in terms of stresses and constitutive tensor
 * @see Parameters
 * @see StressMeasures
 */
template<class TNodeType>
void ConstitutiveLawImpl<TNodeType>::CalculateMaterialResponse(Parameters& rValues,const StressMeasure& rStressMeasure)
{
    switch(rStressMeasure)
    {
    case StressMeasure_PK1:
        CalculateMaterialResponsePK1(rValues);
        break;

    case StressMeasure_PK2:
        CalculateMaterialResponsePK2(rValues);
        break;

    case StressMeasure_Kirchhoff:
        CalculateMaterialResponseKirchhoff(rValues);
        break;

    case StressMeasure_Cauchy:
        CalculateMaterialResponseCauchy(rValues);
        break;

    default:
        KRATOS_ERROR << " Stress Measure not Defined ";
        break;
    }
}

/**
 * Computes the material response in terms of 1st Piola-Kirchhoff stresses and constitutive tensor
 * @see Parameters
 */
template<class TNodeType>
void ConstitutiveLawImpl<TNodeType>::CalculateMaterialResponsePK1(Parameters& rValues)
{
  KRATOS_ERROR << "Calling virtual function for CalculateMaterialResponsePK1";
}

/**
 * Computes the material response in terms of 2nd Piola-Kirchhoff stresses and constitutive tensor
 * @see Parameters
 */
template<class TNodeType>
void ConstitutiveLawImpl<TNodeType>::CalculateMaterialResponsePK2(Parameters& rValues)
{
  KRATOS_ERROR << "Calling virtual function for CalculateMaterialResponsePK2";
}

/**
 * Computes the material response in terms of Kirchhoff stresses and constitutive tensor
 * @see Parameters
 */
template<class TNodeType>
void ConstitutiveLawImpl<TNodeType>::CalculateMaterialResponseKirchhoff(Parameters& rValues)
{
    KRATOS_ERROR << "Calling virtual function for CalculateMaterialResponseKirchhoff";
}

/**
 * Computes the material response in terms of Cauchy stresses and constitutive tensor
 * @see Parameters
 */
template<class TNodeType>
void ConstitutiveLawImpl<TNodeType>::CalculateMaterialResponseCauchy(Parameters& rValues)
{
    KRATOS_ERROR << "Calling virtual function for CalculateMaterialResponseCauchy";
}

/**
 * Computes the material response in terms of Cauchy stresses and constitutive tensor
 * @see Parameters
 */
template<class TNodeType>
void ConstitutiveLawImpl<TNodeType>::CalculateStressResponse(Parameters& rValues, VectorType& rInternalVariables)
{
    KRATOS_ERROR << "Calling virtual function for CalculateMaterialResponseCauchy";
}

/**
 * @brief Initialize the material response,  called by the element in InitializeSolutionStep.
 * @see Parameters
 * @see StressMeasures
 */
template<class TNodeType>
void ConstitutiveLawImpl<TNodeType>::InitializeMaterialResponse(Parameters& rValues,const StressMeasure& rStressMeasure)
{
    switch(rStressMeasure)
    {
        case StressMeasure_PK1:         InitializeMaterialResponsePK1(rValues);
        break;

        case StressMeasure_PK2:         InitializeMaterialResponsePK2(rValues);
        break;

        case StressMeasure_Kirchhoff: 	InitializeMaterialResponseKirchhoff(rValues);
        break;

        case StressMeasure_Cauchy:	InitializeMaterialResponseCauchy(rValues);
        break;

        default:
        KRATOS_ERROR << " Stress Measure not Defined ";
        break;
    }
}

/**
 * @brief Initialize the material response in terms of 1st Piola-Kirchhoff stresses
 * @see Parameters
 */
template<class TNodeType>
void ConstitutiveLawImpl<TNodeType>::InitializeMaterialResponsePK1(Parameters& rValues)
{
    if (this->RequiresInitializeMaterialResponse())
        KRATOS_ERROR << "Calling virtual function for InitializeMaterialResponsePK1. "
            "Please implement InitializeMaterialResponsePK1 or RequiresInitializeMaterialResponse in case this CL does not require it";
}

/**
 * @brief Initialize the material response in terms of 2nd Piola-Kirchhoff stresses
 * @see Parameters
 */
template<class TNodeType>
void ConstitutiveLawImpl<TNodeType>::InitializeMaterialResponsePK2(Parameters& rValues)
{
    if (this->RequiresInitializeMaterialResponse())
        KRATOS_ERROR << "Calling virtual function for InitializeMaterialResponsePK2. "
            "Please implement InitializeMaterialResponsePK2 or RequiresInitializeMaterialResponse in case this CL does not require it";
}

/**
 * @brief Initialize the material response in terms of Kirchhoff stresses
 * @see Parameters
 */
template<class TNodeType>
void ConstitutiveLawImpl<TNodeType>::InitializeMaterialResponseKirchhoff(Parameters& rValues)
{
    if (this->RequiresInitializeMaterialResponse())
        KRATOS_ERROR << "Calling virtual function for InitializeMaterialResponseKirchhoff. "
            "Please implement InitializeMaterialResponseKirchhoff or RequiresInitializeMaterialResponse in case this CL does not require it";
}

/**
 * @brief Initialize the material response in terms of Cauchy stresses
 * @see Parameters
 */
template<class TNodeType>
void ConstitutiveLawImpl<TNodeType>::InitializeMaterialResponseCauchy(Parameters& rValues)
{
    if (this->RequiresInitializeMaterialResponse())
        KRATOS_ERROR << "Calling virtual function for InitializeMaterialResponseCauchy. "
            "Please implement InitializeMaterialResponseCauchy or RequiresInitializeMaterialResponse in case this CL does not require it";
}

/**
 * @brief Finalize the material response,  called by the element in FinalizeSolutionStep.
 * @see Parameters
 * @see StressMeasures
 */
template<class TNodeType>
void ConstitutiveLawImpl<TNodeType>::FinalizeMaterialResponse(Parameters& rValues,const StressMeasure& rStressMeasure)
{
    switch(rStressMeasure)
    {
    case StressMeasure_PK1:
        FinalizeMaterialResponsePK1(rValues);
        break;

    case StressMeasure_PK2:
        FinalizeMaterialResponsePK2(rValues);
        break;

    case StressMeasure_Kirchhoff:
        FinalizeMaterialResponseKirchhoff(rValues);
        break;

    case StressMeasure_Cauchy:
        FinalizeMaterialResponseCauchy(rValues);
        break;

    default:
        KRATOS_ERROR << " Stress Measure not Defined ";
        break;
    }
}

/**
 * Updates the material response in terms of 1st Piola-Kirchhoff stresses
 * @see Parameters
 */
template<class TNodeType>
void ConstitutiveLawImpl<TNodeType>::FinalizeMaterialResponsePK1(Parameters& rValues)
{
  KRATOS_ERROR << "Calling virtual function for FinalizeMaterialResponsePK1";
}

/**
 * Updates the material response in terms of 2nd Piola-Kirchhoff stresses
 * @see Parameters
 */
template<class TNodeType>
void ConstitutiveLawImpl<TNodeType>::FinalizeMaterialResponsePK2(Parameters& rValues)
{
  KRATOS_ERROR << "Calling virtual function for FinalizeMaterialResponsePK2";
}

/**
 * Updates the material response in terms of Kirchhoff stresses
 * @see Parameters
 */
template<class TNodeType>
void ConstitutiveLawImpl<TNodeType>::FinalizeMaterialResponseKirchhoff(Parameters& rValues)
{
  KRATOS_ERROR << "Calling virtual function for FinalizeMaterialResponseKirchhoff";
}

/**
 * Updates the material response in terms of Cauchy stresses
 * @see Parameters
 */
template<class TNodeType>
void ConstitutiveLawImpl<TNodeType>::FinalizeMaterialResponseCauchy(Parameters& rValues)
{
  KRATOS_ERROR << "Calling virtual function for FinalizeMaterialResponseCauchy";
}

/**
 * This can be used in order to reset all internal variables of the
 * constitutive law (e.g. if a model should be reset to its reference state)
 * @param rMaterialProperties the Properties instance of the current element
 * @param rElementGeometry the geometry of the current element
 * @param rShapeFunctionsValues the shape functions values in the current integration point
 * @param the current ProcessInfo instance
 */
template<class TNodeType>
void ConstitutiveLawImpl<TNodeType>::ResetMaterial(const Properties& rMaterialProperties,
                                    const GeometryType& rElementGeometry,
                                    const Vector& rShapeFunctionsValues)
{
    KRATOS_ERROR << "Calling virtual function for ResetMaterial";
}

/**
 * This can be used in order to rewind all internal variables of the
 * constitutive law (e.g. if a model should be rewinded to a previous state)
 * REMARK: EXPERIMENTAL METHOD
 * @param rMaterialProperties the Properties instance of the current element
 * @param rElementGeometry the geometry of the current element
 * @param rShapeFunctionsValues the shape functions values in the current integration point
 * @param the current ProcessInfo instance
 */
template<class TNodeType>
void ConstitutiveLawImpl<TNodeType>::RewindMaterial(const Properties& rMaterialProperties,
                                     const GeometryType& rElementGeometry,
                                     const Vector& rShapeFunctionsValues)
{
    KRATOS_ERROR << "Calling virtual function for RewindMaterial";
}

/**
 * Methods to transform strain Vectors:
 * @param rStrainVector the strain tensor in matrix which its stress measure will be changed
 * @param rF the DeformationGradientF matrix between the configurations
 * @param rdetF the determinant of the DeformationGradientF matrix between the configurations
 * @param rStrainInitial the measure of stress of the given  rStrainVector
 * @param rStrainFinal the measure of stress of the returned rStrainVector
 */
template<class TNodeType>
typename ConstitutiveLawImpl<TNodeType>::VectorType& ConstitutiveLawImpl<TNodeType>::TransformStrains(VectorType& rStrainVector,
                                              const MatrixType& rF,
                                              StrainMeasure rStrainInitial,
                                              StrainMeasure rStrainFinal)
{
    switch(rStrainInitial)
    {
    case StrainMeasure_GreenLagrange:

        switch(rStrainFinal)
        {
        case StrainMeasure_GreenLagrange:
            break;

        case StrainMeasure_Almansi:
        {
            MatrixType StrainMatrix = MathUtils<DataType>::StrainVectorToTensor( rStrainVector );

            CoVariantPushForward (StrainMatrix, rF);  //Almansi

            rStrainVector = MathUtils<DataType>::StrainTensorToVector( StrainMatrix, rStrainVector.size() );
        }
        break;

        case StrainMeasure_Hencky_Material:
            KRATOS_ERROR << "Hencky strain has no transformation coded";
            break;

        case StrainMeasure_Hencky_Spatial:
            KRATOS_ERROR << "Hencky strain has no transformation coded";
            break;

        default:
            KRATOS_ERROR << "FINAL STRAIN NOT DEFINED in StrainTransformation";
            break;
        }

        break;

    case StrainMeasure_Almansi:

        switch(rStrainFinal)
        {
        case StrainMeasure_GreenLagrange:
        {
            MatrixType StrainMatrix = MathUtils<DataType>::StrainVectorToTensor( rStrainVector );

            CoVariantPullBack (StrainMatrix,rF);  //GreenLagrange

            rStrainVector = MathUtils<DataType>::StrainTensorToVector( StrainMatrix, rStrainVector.size() );
        }
        break;

        case StrainMeasure_Almansi:
            break;

        case StrainMeasure_Hencky_Material:
            KRATOS_ERROR << "Hencky strain has no transformation coded";
            break;

        case StrainMeasure_Hencky_Spatial:
            KRATOS_ERROR << "Hencky strain has no transformation coded";
            break;

        default:
            KRATOS_ERROR << "FINAL STRAIN NOT DEFINED in StrainTransformation";
            break;
        }

        break;

    case StrainMeasure_Hencky_Material:
        KRATOS_ERROR << "Hencky strain has no transformation coded";
        break;

    case StrainMeasure_Hencky_Spatial:
        KRATOS_ERROR << "Hencky strain has no transformation coded";
        break;

    default:
        KRATOS_ERROR << "Measure of strain NOT DEFINED in Strains Transformation";
        break;
    }

    return rStrainVector;
}

/**
 * Methods to transform stress Matrices:
 * @param rStressMatrix the stress tensor in matrix which its stress measure will be changed
 * @param rF the DeformationGradientF matrix between the configurations
 * @param rdetF the determinant of the DeformationGradientF matrix between the configurations
 * @param rStressInitial the measure of stress of the given  rStressMatrix
 * @param rStressFinal the measure of stress of the returned rStressMatrix
 */
template<class TNodeType>
typename ConstitutiveLawImpl<TNodeType>::MatrixType& ConstitutiveLawImpl<TNodeType>::TransformStresses(MatrixType& rStressMatrix,
                                               const MatrixType& rF,
                                               const DataType& rdetF,
                                               StressMeasure rStressInitial,
                                               StressMeasure rStressFinal)
{
    VectorType StressVector;

    StressVector = MathUtils<DataType>::StressTensorToVector( rStressMatrix );

    StressVector=TransformStresses(StressVector,rF,rdetF,rStressInitial,rStressFinal);

    rStressMatrix = MathUtils<DataType>::StressVectorToTensor( StressVector );

    return rStressMatrix;
}

/**
 * Methods to transform stress Vectors:
 * @param rStressVector the stress tensor in matrix which its stress measure will be changed
 * @param rF the DeformationGradientF matrix between the configurations
 * @param rdetF the determinant of the DeformationGradientF matrix between the configurations
 * @param rStressInitial the measure of stress of the given  rStressVector
 * @param rStressFinal the measure of stress of the returned rStressVector
 */
template<class TNodeType>
typename ConstitutiveLawImpl<TNodeType>::VectorType& ConstitutiveLawImpl<TNodeType>::TransformStresses(VectorType& rStressVector,
                                               const MatrixType& rF,
                                               const DataType& rdetF,
                                               StressMeasure rStressInitial,
                                               StressMeasure rStressFinal)
{
    switch(rStressInitial)
    {
    case StressMeasure_PK1:

        TransformPK1Stresses(rStressVector,rF,rdetF,rStressFinal);

        break;

    case StressMeasure_PK2:

        TransformPK2Stresses(rStressVector,rF,rdetF,rStressFinal);

        break;

    case StressMeasure_Kirchhoff:

        TransformKirchhoffStresses(rStressVector,rF,rdetF,rStressFinal);

        break;

    case StressMeasure_Cauchy:

        TransformCauchyStresses(rStressVector,rF,rdetF,rStressFinal);

        break;

    default:
        KRATOS_ERROR <<"INITIAL STRESS NOT DEFINED in StressTransformation";
        break;
    }

    return rStressVector;
}

/**
 * Methods to transform stress Vectors specialized with the initial stress Measure PK1:
 * @param rStressVector the stress tensor in matrix which its stress measure will be changed
 * @param rF the DeformationGradientF matrix between the configurations
 * @param rdetF the determinant of the DeformationGradientF matrix between the configurations
 * @param rStressFinal the measure of stress of the returned rStressVector
 */
template<class TNodeType>
typename ConstitutiveLawImpl<TNodeType>::VectorType& ConstitutiveLawImpl<TNodeType>::TransformPK1Stresses(VectorType& rStressVector,
                                                  const MatrixType& rF,
                                                  const DataType& rdetF,
                                                  StressMeasure rStressFinal)
{
    unsigned int size = rF.size1(); //WorkingSpaceDimension();

    switch(rStressFinal)
    {
    case StressMeasure_PK1:
        break;

    case StressMeasure_PK2:
    {
        MatrixType StressMatrix = MathUtils<DataType>::StressVectorToTensor( rStressVector );
        MatrixType InvF ( size, size );
        DataType J;
        MathUtils<DataType>::InvertMatrix( rF, InvF, J );

        StressMatrix = prod( InvF, StressMatrix ); //PK2

        rStressVector = MathUtils<DataType>::StressTensorToVector( StressMatrix, rStressVector.size() );
    }
    break;

    case StressMeasure_Kirchhoff:
    {
        MatrixType StressMatrix = MathUtils<DataType>::StressVectorToTensor( rStressVector );
        MatrixType InvF ( size, size );
        DataType J;
        MathUtils<DataType>::InvertMatrix( rF, InvF, J );

        StressMatrix = prod( InvF, StressMatrix ); //PK2

        ContraVariantPushForward (StressMatrix,rF); //Kirchhoff

        rStressVector = MathUtils<DataType>::StressTensorToVector( StressMatrix, rStressVector.size() );
    }
    break;

    case StressMeasure_Cauchy:
    {
        MatrixType StressMatrix = MathUtils<DataType>::StressVectorToTensor( rStressVector );
        MatrixType InvF ( size, size );
        DataType J;
        MathUtils<DataType>::InvertMatrix( rF, InvF, J );

        StressMatrix = prod( InvF, StressMatrix ); //PK2

        ContraVariantPushForward (StressMatrix,rF); //Kirchhoff

        StressMatrix/=J; //Cauchy

        rStressVector = MathUtils<DataType>::StressTensorToVector( StressMatrix, rStressVector.size() );
    }
    break;

    default:
        KRATOS_ERROR <<"FINAL STRESS NOT DEFINED in StressTransformation";
        break;
    }

    return rStressVector;
}

/**
 * Methods to transform stress Vectors specialized with the initial stress Measure PK2:
 * @param rStressVector the stress tensor in matrix which its stress measure will be changed
 * @param rF the DeformationGradientF matrix between the configurations
 * @param rdetF the determinant of the DeformationGradientF matrix between the configurations
 * @param rStressFinal the measure of stress of the returned rStressVector
 */
template<class TNodeType>
typename ConstitutiveLawImpl<TNodeType>::VectorType& ConstitutiveLawImpl<TNodeType>::TransformPK2Stresses(VectorType& rStressVector,
                                                  const MatrixType& rF,
                                                  const DataType& rdetF,
                                                  StressMeasure rStressFinal)
{
    switch(rStressFinal)
    {
    case StressMeasure_PK1:
    {
        MatrixType StressMatrix = MathUtils<DataType>::StressVectorToTensor( rStressVector );

        StressMatrix = prod( rF, StressMatrix ); //PK1

        rStressVector = MathUtils<DataType>::StressTensorToVector( StressMatrix, rStressVector.size() );
    }
    break;

    case StressMeasure_PK2:
        break;

    case StressMeasure_Kirchhoff:
    {
        MatrixType StressMatrix = MathUtils<DataType>::StressVectorToTensor( rStressVector );

        ContraVariantPushForward (StressMatrix,rF); //Kirchhoff

        rStressVector = MathUtils<DataType>::StressTensorToVector( StressMatrix, rStressVector.size() );
    }
    break;

    case StressMeasure_Cauchy:
    {
        MatrixType StressMatrix = MathUtils<DataType>::StressVectorToTensor( rStressVector );

        ContraVariantPushForward (StressMatrix,rF); //Kirchhoff

        if(std::abs(rdetF) != 0)
            StressMatrix /= rdetF; //Cauchy

        rStressVector = MathUtils<DataType>::StressTensorToVector( StressMatrix, rStressVector.size() );

    }
    break;

    default:
        KRATOS_ERROR <<"FINAL STRESS NOT DEFINED in StressTransformation";
        break;
    }

    return rStressVector;
}

/**
 * Methods to transform stress Vectors specialized with the initial stress Measure Kirchooff:
 * @param rStressVector the stress tensor in matrix which its stress measure will be changed
 * @param rF the DeformationGradientF matrix between the configurations
 * @param rdetF the determinant of the DeformationGradientF matrix between the configurations
 * @param rStressFinal the measure of stress of the returned rStressVector
 */
template<class TNodeType>
typename ConstitutiveLawImpl<TNodeType>::VectorType& ConstitutiveLawImpl<TNodeType>::TransformKirchhoffStresses(VectorType& rStressVector,
                                                        const MatrixType& rF,
                                                        const DataType& rdetF,
                                                        StressMeasure rStressFinal)
{
    switch(rStressFinal)
    {
    case StressMeasure_PK1:
    {
        MatrixType StressMatrix = MathUtils<DataType>::StressVectorToTensor( rStressVector );

        ContraVariantPullBack (StressMatrix,rF);  //PK2

        StressMatrix = prod( rF, StressMatrix ); //PK1

        rStressVector = MathUtils<DataType>::StressTensorToVector( StressMatrix, rStressVector.size() );
    }
    break;

    case StressMeasure_PK2:
    {
        MatrixType StressMatrix = MathUtils<DataType>::StressVectorToTensor( rStressVector );

        ContraVariantPullBack (StressMatrix,rF);  //PK2

        rStressVector = MathUtils<DataType>::StressTensorToVector( StressMatrix, rStressVector.size() );
    }
    break;

    case StressMeasure_Kirchhoff:
        break;

    case StressMeasure_Cauchy:
    {
        if(std::abs(rdetF) != 0)
            rStressVector/=rdetF; //Cauchy
    }
    break;

    default:
        KRATOS_ERROR <<"FINAL STRESS NOT DEFINED in StressTransformation";
        break;
    }

    return rStressVector;
}

/**
 * Methods to transform stress Vectors specialized with the initial stress Measure Cauchy:
 * @param rStressVector the stress tensor in matrix which its stress measure will be changed
 * @param rF the DeformationGradientF matrix between the configurations
 * @param rdetF the determinant of the DeformationGradientF matrix between the configurations
 * @param rStressFinal the measure of stress of the returned rStressVector
 */
template<class TNodeType>
typename ConstitutiveLawImpl<TNodeType>::VectorType& ConstitutiveLawImpl<TNodeType>::TransformCauchyStresses(VectorType& rStressVector,
                                                     const MatrixType& rF,
                                                     const DataType& rdetF,
                                                     StressMeasure rStressFinal)
{

    switch(rStressFinal)
    {
    case StressMeasure_PK1:
    {
        rStressVector*=rdetF; //Kirchhoff

        MatrixType StressMatrix = MathUtils<DataType>::StressVectorToTensor( rStressVector );

        ContraVariantPullBack (StressMatrix,rF);  //PK2

        StressMatrix = prod( rF, StressMatrix ); //PK1

        rStressVector = MathUtils<DataType>::StressTensorToVector( StressMatrix, rStressVector.size() );
    }
    break;

    case StressMeasure_PK2:
    {
        rStressVector*=rdetF; //Kirchhoff

        MatrixType StressMatrix = MathUtils<DataType>::StressVectorToTensor( rStressVector );

        ContraVariantPullBack (StressMatrix,rF);  //PK2

        rStressVector = MathUtils<DataType>::StressTensorToVector( StressMatrix, rStressVector.size() );
    }
    break;

    case StressMeasure_Kirchhoff:

        rStressVector*=rdetF; //Kirchhoff

        break;

    case StressMeasure_Cauchy:
        break;

    default:
        KRATOS_ERROR <<"FINAL STRESS NOT DEFINED in StressTransformation";
        break;
    }

    return rStressVector;
}

/**
 * Methods to transform Constitutive Matrices:
 * @param rConstitutiveMatrix the constitutive matrix
 * @param rF the DeformationGradientF matrix between the configurations
 */

/**
 * This method performs a pull-back of the constitutive matrix
 */
template<class TNodeType>
void ConstitutiveLawImpl<TNodeType>::PullBackConstitutiveMatrix ( MatrixType& rConstitutiveMatrix, const MatrixType& rF )
{
    MatrixType OriginalConstitutiveMatrix = rConstitutiveMatrix;

    rConstitutiveMatrix.clear();

    MatrixType InverseF ( 3, 3 );
    DataType detF = 0;
    MathUtils<DataType>::InvertMatrix( rF, InverseF, detF);

    ConstitutiveMatrixTransformation( rConstitutiveMatrix, OriginalConstitutiveMatrix, InverseF );
}

/**
 * This method performs a push-forward of the constitutive matrix
 */
template<class TNodeType>
void ConstitutiveLawImpl<TNodeType>::PushForwardConstitutiveMatrix( MatrixType& rConstitutiveMatrix, const MatrixType& rF )
{
    MatrixType OriginalConstitutiveMatrix = rConstitutiveMatrix;

    rConstitutiveMatrix.clear();

    ConstitutiveMatrixTransformation( rConstitutiveMatrix, OriginalConstitutiveMatrix, rF );
}

/**
 * This function is designed to be called once to check compatibility with element
 * @param rFeatures
 */
template<class TNodeType>
void ConstitutiveLawImpl<TNodeType>::GetLawFeatures(Features& rFeatures)
{
    KRATOS_ERROR << "Calling virtual function for GetConstitutiveLawFeatures";
}

/**
 * This function is designed to be called once to perform all the checks needed
 * on the input provided. Checks can be "expensive" as the function is designed
 * to catch user's errors.
 * @param rMaterialProperties
 * @param rElementGeometry
 * @param rCurrentProcessInfo
 * @return
 */
template<class TNodeType>
int ConstitutiveLawImpl<TNodeType>::Check(const Properties& rMaterialProperties,
                           const GeometryType& rElementGeometry,
                           const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    return 0;
    KRATOS_CATCH("");
}


//*** PROTECTED METHODS: ***//

/**
 * This method performs a contra-variant push-forward between to tensors
 * i.e. 2nd PK stress to Kirchhoff stress
 */
template<class TNodeType>
void ConstitutiveLawImpl<TNodeType>::ContraVariantPushForward( MatrixType& rMatrix, const MatrixType& rF)  //i.e. 2nd PK stress to Kirchhoff stress
{
    unsigned int size = rF.size1(); //WorkingSpaceDimension();
    MatrixType temp ( size, size );

    noalias( temp )     = prod( rF, rMatrix );
    noalias( rMatrix )  = prod( temp, trans( rF ) );

}

/**
 * This method performs a contra-variant pull-back between to tensors
 * i.e. Kirchhoff stress to 2nd PK stress
 */
template<class TNodeType>
void ConstitutiveLawImpl<TNodeType>::ContraVariantPullBack( MatrixType& rMatrix, const MatrixType& rF)     //i.e. Kirchhoff stress to 2nd PK stress
{
    unsigned int size = rF.size1(); //WorkingSpaceDimension();
    MatrixType InvF ( size, size );
    DataType J;
    MathUtils<DataType>::InvertMatrix( rF, InvF, J );

    MatrixType temp ( size, size );

    noalias( temp )    = prod( InvF, rMatrix );
    noalias( rMatrix ) = prod( temp, trans( InvF ) );
}

/**
 * This method performs a co-variant push-forward between to tensors
 * i.e. Green-Lagrange strain to Almansi strain
 */
template<class TNodeType>
void ConstitutiveLawImpl<TNodeType>::CoVariantPushForward( MatrixType& rMatrix, const MatrixType& rF)      //i.e. Green-Lagrange strain to Almansi strain
{
    unsigned int size = rF.size1(); //WorkingSpaceDimension();
    MatrixType InvF ( size, size );
    DataType J;
    MathUtils<DataType>::InvertMatrix( rF, InvF, J );

    MatrixType temp ( size, size );

    noalias( temp )     = prod( trans( InvF ), rMatrix );
    noalias( rMatrix )  = prod( temp, InvF );
}

/**
 * This method performs a co-variant pull-back between to tensors
 * i.e. Almansi strain to Green-Lagrange strain
 */
template<class TNodeType>
void ConstitutiveLawImpl<TNodeType>::CoVariantPullBack( MatrixType& rMatrix, const MatrixType& rF)         //i.e. Almansi strain to Green-Lagrange strain
{
    unsigned int size = rF.size1(); //WorkingSpaceDimension();
    MatrixType temp ( size, size );

    noalias( temp )     = prod( trans( rF ), rMatrix );
    noalias( rMatrix )  = prod( temp, rF );
}

/**
 * This method performs a pull-back or a push-forward between two constitutive matrices
 */
template<class TNodeType>
void ConstitutiveLawImpl<TNodeType>::ConstitutiveMatrixTransformation ( MatrixType& rConstitutiveMatrix,
                                                         const MatrixType& rOriginalConstitutiveMatrix,
                                                         const MatrixType&  rF )
{
    unsigned int size = rOriginalConstitutiveMatrix.size1();
    if(  size == 6 )
    {

        for(unsigned int i=0; i<6; i++)
        {
            for(unsigned int j=0; j<6; j++)
            {
                rConstitutiveMatrix( i, j ) = TransformConstitutiveComponent(rConstitutiveMatrix( i, j ), rOriginalConstitutiveMatrix, rF,
                                              this->msIndexVoigt3D6C[i][0], this->msIndexVoigt3D6C[i][1], this->msIndexVoigt3D6C[j][0], this->msIndexVoigt3D6C[j][1]);
            }

        }
    }
    else if( size == 4 )
    {


        for(unsigned int i=0; i<4; i++)
        {
            for(unsigned int j=0; j<4; j++)
            {
                rConstitutiveMatrix( i, j ) = TransformConstitutiveComponent(rConstitutiveMatrix( i, j ), rOriginalConstitutiveMatrix, rF,
                                              this->msIndexVoigt2D4C[i][0], this->msIndexVoigt2D4C[i][1], this->msIndexVoigt2D4C[j][0], this->msIndexVoigt2D4C[j][1]);
            }

        }
    }
    else if( size == 3 )
    {


        for(unsigned int i=0; i<3; i++)
        {
            for(unsigned int j=0; j<3; j++)
            {
                rConstitutiveMatrix( i, j ) = TransformConstitutiveComponent(rConstitutiveMatrix( i, j ), rOriginalConstitutiveMatrix, rF,
                                              this->msIndexVoigt2D3C[i][0], this->msIndexVoigt2D3C[i][1], this->msIndexVoigt2D3C[j][0], this->msIndexVoigt2D3C[j][1]);
            }

        }
    }
}

/**
 * This method performs a pull-back or a push-forward between two constitutive tensor components
 */
template<class TNodeType>
typename ConstitutiveLawImpl<TNodeType>::DataType& ConstitutiveLawImpl<TNodeType>::TransformConstitutiveComponent(DataType& rCabcd,
                                                          const MatrixType& rConstitutiveMatrix,
                                                          const MatrixType& rF,
                                                          const unsigned int a, const unsigned int b,
                                                          const unsigned int c, const unsigned int d)
{

    rCabcd = 0;
    DataType Cijkl=0;

    unsigned int dimension = rF.size1();

    //Cabcd
    for(unsigned int j=0; j<dimension; j++)
    {
        for(unsigned int l=0; l<dimension; l++)
        {
            for(unsigned int k=0; k<dimension; k++)
            {
                for(unsigned int i=0; i<dimension; i++)
                {
                    //Cijkl
                    rCabcd +=rF(a,i)*rF(b,j)*rF(c,k)*rF(d,l)*GetConstitutiveComponent(Cijkl,rConstitutiveMatrix,i,j,k,l);
                }
            }
        }
    }

    return rCabcd;
}

/**
 * This method gets the constitutive tensor components
 * from a consitutive matrix supplied in voigt notation
 */
template<class TNodeType>
typename ConstitutiveLawImpl<TNodeType>::DataType& ConstitutiveLawImpl<TNodeType>::GetConstitutiveComponent(DataType& rCabcd,
                                                    const MatrixType& rConstitutiveMatrix,
                                                    const unsigned int a, const unsigned int b,
                                                    const unsigned int c, const unsigned int d)
{
    // matrix indices
    unsigned int k=0, l= 0;

    unsigned int size = rConstitutiveMatrix.size1();

    if( size == 3 )
    {
        //index k
        for(unsigned int i=0; i<3; i++)
        {
            if( a == b )
            {
                if( this->msIndexVoigt2D3C[i][0] == a && this->msIndexVoigt2D3C[i][1] == b )
                {
                    k = i;
                    break;
                }
            }
            else
            {
                if( (this->msIndexVoigt2D3C[i][0] == a && this->msIndexVoigt2D3C[i][1] == b) ||
                        (this->msIndexVoigt2D3C[i][1] == a && this->msIndexVoigt2D3C[i][0] == b) )
                {
                    k = i;
                    break;
                }
            }
        }

        //index l
        for(unsigned int i=0; i<3; i++)
        {
            if( c == d )
            {
                if( this->msIndexVoigt2D3C[i][0] == c && this->msIndexVoigt2D3C[i][1] == d )
                {
                    l = i;
                    break;
                }
            }
            else
            {
                if( (this->msIndexVoigt2D3C[i][0] == c && this->msIndexVoigt2D3C[i][1] == d) ||
                        (this->msIndexVoigt2D3C[i][1] == c && this->msIndexVoigt2D3C[i][0] == d) )
                {
                    l = i;
                    break;
                }
            }
        }
    }
    else if( size == 4 )
    {
        //index k
        for(unsigned int i=0; i<4; i++)
        {
            if( a == b )
            {
                if( this->msIndexVoigt2D4C[i][0] == a && this->msIndexVoigt2D4C[i][1] == b )
                {
                    k = i;
                    break;
                }
            }
            else
            {
                if( (this->msIndexVoigt2D4C[i][0] == a && this->msIndexVoigt2D4C[i][1] == b) ||
                        (this->msIndexVoigt2D4C[i][1] == a && this->msIndexVoigt2D4C[i][0] == b) )
                {
                    k = i;
                    break;
                }
            }
        }

        //index l
        for(unsigned int i=0; i<4; i++)
        {
            if( c == d )
            {
                if( this->msIndexVoigt2D4C[i][0] == c && this->msIndexVoigt2D4C[i][1] == d )
                {
                    l = i;
                    break;
                }
            }
            else
            {
                if( (this->msIndexVoigt2D4C[i][0] == c && this->msIndexVoigt2D4C[i][1] == d) ||
                        (this->msIndexVoigt2D4C[i][1] == c && this->msIndexVoigt2D4C[i][0] == d) )
                {
                    l = i;
                    break;
                }
            }
        }
    }
    else if( size == 6 )
    {
        //index k
        for(unsigned int i=0; i<6; i++)
        {
            if( a == b )
            {
                if( this->msIndexVoigt3D6C[i][0] == a && this->msIndexVoigt3D6C[i][1] == b )
                {
                    k = i;
                    break;
                }
            }
            else
            {
                if( (this->msIndexVoigt3D6C[i][0] == a && this->msIndexVoigt3D6C[i][1] == b) ||
                        (this->msIndexVoigt3D6C[i][1] == a && this->msIndexVoigt3D6C[i][0] == b) )
                {
                    k = i;
                    break;
                }
            }
        }

        //index l
        for(unsigned int i=0; i<6; i++)
        {
            if( c == d )
            {
                if( this->msIndexVoigt3D6C[i][0] == c && this->msIndexVoigt3D6C[i][1] == d )
                {
                    l = i;
                    break;
                }
            }
            else
            {
                if( (this->msIndexVoigt3D6C[i][0] == c && this->msIndexVoigt3D6C[i][1] == d) ||
                        (this->msIndexVoigt3D6C[i][1] == c && this->msIndexVoigt3D6C[i][0] == d) )
                {
                    l = i;
                    break;
                }
            }
        }
    }

    rCabcd = rConstitutiveMatrix(k,l);

    return rCabcd;
}

//*** OUTDATED METHODS: ***//

/**
 * Computes the material response in terms of stresses and algorithmic tangent
 * @param StrainVector the current strains (total strains, input)
 * @param DeformationGradient the current deformation gradient (can be an empty matrix if a linear strain measure is used)
 * @param StressVector the computed stresses (output)
 * @param algorithmicTangent the material tangent matrix (output)
 * @param rCurrentProcessInfo current ProcessInfo instance
 * @param rMaterialProperties the material's Properties object
 * @param rElementGeometry the element's geometry
 * @param rShapeFunctionsValues the shape functions values in the current integration pointer
 * @param CalculateStresses flag whether or not to compute the stress response
 * @param CalculateTangent flag to determine if to compute the material tangent
 * NOTE: the CalculateTangent flag is defined as int to allow for distinctive variants of the tangent
 * @param SaveInternalVariables flag whether or not to store internal (history) variables
 */
template<class TNodeType>
void ConstitutiveLawImpl<TNodeType>::CalculateMaterialResponse(const VectorType& StrainVector,
                                                const MatrixType& DeformationGradient,
                                                VectorType& StressVector,
                                                MatrixType& AlgorithmicTangent,
                                                const ProcessInfo& rCurrentProcessInfo,
                                                const Properties& rMaterialProperties,
                                                const GeometryType& rElementGeometry,
                                                const Vector& rShapeFunctionsValues,
                                                bool CalculateStresses,
                                                int CalculateTangent,
                                                bool SaveInternalVariables)
{
    KRATOS_ERROR << "Calling virtual function for CalculateMaterialResponse";
}

/**
 * Computes the volumetric part of the material response in terms of stresses and algorithmic tangent
 * @param StrainVector the current strains (total strains, input)
 * @param DeformationGradient the current deformation gradient (can be an empty matrix if a linear
 * @param StressVector the computed stresses (output)
 * @param algorithmicTangent the material tangent matrix (output)
 * @param rCurrentProcessInfo current ProcessInfo instance
 * @param rMaterialProperties the material's Properties object
 * @param rElementGeometry the element's geometry
 * @param rShapeFunctionsValues the shape functions values in the current integration pointer
 * @param CalculateStresses flag whether or not to compute the stress response
 * @param CalculateTangent flag to determine if to compute the material tangent
 * NOTE: the CalculateTangent flag is defined as int to allow for distinctive variants of the tangent
 * @param SaveInternalVariables flag whether or not to store internal (history) variables
 */
template<class TNodeType>
void ConstitutiveLawImpl<TNodeType>::CalculateVolumetricResponse(const DataType VolumetricStrain,
                                                  const MatrixType& DeformationGradient,
                                                  DataType& VolumetricStress,
                                                  DataType& AlgorithmicBulk,
                                                  const ProcessInfo& rCurrentProcessInfo,
                                                  const Properties& rMaterialProperties,
                                                  const GeometryType& rElementGeometry,
                                                  const Vector& rShapeFunctionsValues,
                                                  bool CalculateStresses,
                                                  int CalculateTangent,
                                                  bool SaveInternalVariables)
{
    KRATOS_ERROR << "Calling virtual function for CalculateVolumetricResponse";
}

/**
 * Computes the deviatoric part of the material response in terms of stresses and algorithmic tangent
 * @param StrainVector the current strains (total strains, input)
 * @param DeformationGradient the current deformation gradient (can be an empty matrix if a linear
 * @param StressVector the computed stresses (output)
 * @param algorithmicTangent the material tangent matrix (output)
 * @param rCurrentProcessInfo current ProcessInfo instance
 * @param rMaterialProperties the material's Properties object
 * TODO: add proper definition for algorithmic tangent
 */
template<class TNodeType>
void ConstitutiveLawImpl<TNodeType>::CalculateDeviatoricResponse(const VectorType& StrainVector,
                                                  const MatrixType& DeformationGradient,
                                                  VectorType& StressVector,
                                                  MatrixType& AlgorithmicTangent,
                                                  const ProcessInfo& rCurrentProcessInfo,
                                                  const Properties& rMaterialProperties,
                                                  const GeometryType& rElementGeometry,
                                                  const Vector& rShapeFunctionsValues,
                                                  bool CalculateStresses,
                                                  int CalculateTangent,
                                                  bool SaveInternalVariables)
{
    KRATOS_ERROR << "Calling virtual function for CalculateDeviatoricResponse";
}


// VM
template<class TNodeType>
void ConstitutiveLawImpl<TNodeType>::CalculateCauchyStresses(VectorType& Cauchy_StressVector,
                                              const MatrixType& F,
                                              const VectorType& PK2_StressVector,
                                              const VectorType& GreenLagrangeStrainVector)
{
}

template class KRATOS_API(KRATOS_CORE) ConstitutiveLawImpl<Node<3, KRATOS_DOUBLE_TYPE, Dof<KRATOS_DOUBLE_TYPE> > >;
template class KRATOS_API(KRATOS_CORE) ConstitutiveLawImpl<Node<3, KRATOS_DOUBLE_TYPE, Dof<KRATOS_COMPLEX_TYPE> > >;
template class KRATOS_API(KRATOS_CORE) ConstitutiveLawImpl<Node<3, KRATOS_COMPLEX_TYPE, Dof<KRATOS_COMPLEX_TYPE> > >;

template class KRATOS_API(KRATOS_CORE) KratosComponents<ConstitutiveLaw>;
template class KRATOS_API(KRATOS_CORE) KratosComponents<ComplexConstitutiveLaw>;
template class KRATOS_API(KRATOS_CORE) KratosComponents<GComplexConstitutiveLaw>;

} /* namespace Kratos.*/
