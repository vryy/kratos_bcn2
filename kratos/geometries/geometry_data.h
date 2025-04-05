//
//   Project Name:        Kratos
//   Last Modified by:    $Author:   JMCarbonell $
//   Date:                $Date:   December 2015 $
//   Revision:            $Revision:         1.7 $
//
//

#if !defined(KRATOS_GEOMETRY_DATA_H_INCLUDED )
#define  KRATOS_GEOMETRY_DATA_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/ublas_interface.h"
#include "integration/integration_point.h"


namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/** GeometryData base class.As a base class GeometryData has all the common
    interface of Kratos' geometries. Also it contains array of
    pointers to its points, reference to shape functions values in
    all integrations points and also local gradients of shape
    functions evaluated in all integrations points.

    @see Geometry
    @see Element
*/
class GeometryData
{
public:
    ///@name  Enum's
    ///@{

    /** Integration methods implemented in geometry. Each geometry
    supposed to have different integration method for
    integrating. These enums are used to refere to each
    different integration methods:

    - GI_GAUSS_1 gaussian integration with order 1.
    - GI_GAUSS_2 gaussian integration with order 2.
    - GI_GAUSS_3 gaussian integration with order 3.
    - GI_GAUSS_4 gaussian integration with order 4.
    - GI_GAUSS_5 gaussian integration with order 5.
    */
    enum class IntegrationMethod {
        GI_GAUSS_1,
        GI_GAUSS_2,
        GI_GAUSS_3,
        GI_GAUSS_4,
        GI_GAUSS_5,
        GI_EXTENDED_GAUSS_1,
        GI_EXTENDED_GAUSS_2,
        GI_EXTENDED_GAUSS_3,
        GI_EXTENDED_GAUSS_4,
        GI_EXTENDED_GAUSS_5,
        NumberOfIntegrationMethods
    };

    enum class KratosGeometryFamily
    {
        Kratos_NoElement,
        Kratos_Point,
        Kratos_Linear,
        Kratos_Triangle,
        Kratos_Quadrilateral,
        Kratos_Tetrahedra,
        Kratos_Hexahedra,
        Kratos_Prism,
        Kratos_NURBS,
        Kratos_Polygon,
        Kratos_Polyhedron,
        Kratos_generic_family
    };

    enum class KratosGeometryType
    {
        Kratos_Hexahedra3D20,
        Kratos_Hexahedra3D27,
        Kratos_Hexahedra3D8,
        Kratos_Prism3D15,
        Kratos_Prism3D6,
        Kratos_Quadrilateral2D4,
        Kratos_Quadrilateral2D8,
        Kratos_Quadrilateral2D9,
        Kratos_Quadrilateral3D4,
        Kratos_Quadrilateral3D8,
        Kratos_Quadrilateral3D9,
        Kratos_Tetrahedra3D10,
        Kratos_Tetrahedra3D4,
        Kratos_Triangle2D3,
        Kratos_Triangle2D6,
        Kratos_Triangle3D3,
        Kratos_Triangle3D6,
        Kratos_Line2D2,
        Kratos_Line2D3,
        Kratos_Line3D2,
        Kratos_Line3D3,
        Kratos_Point2D,
        Kratos_Point3D,
        Kratos_Sphere3D1,
        Kratos_Bezier1D,
        Kratos_Bezier1D2,
        Kratos_Bezier1D3,
        Kratos_Bezier2D,
        Kratos_Bezier2D3,
        Kratos_Bezier3D,
        Kratos_Tritagon, //3 // it's not a real name but need to avoid clashing with Kratos_Triangle
        Kratos_Tetragon, // 4
        Kratos_Pentagon, // 5
        Kratos_Hexagon, // 6
        Kratos_Heptagon, // 7
        Kratos_Octagon, // 8
        Kratos_Nonagon, // 9
        Kratos_Decagon, // 10
        Kratos_Hendecagon, // 11
        Kratos_Dodecagon, // 12
        Kratos_Triskaidecagon, // 13
        Kratos_Tetrakaidecagon, // 14
        Kratos_Pentadecagon, // 15
        Kratos_Hexakaidecagon, // 16
        Kratos_Heptadecagon, // 17
        Kratos_Octakaidecagon, // 18
        Kratos_Enneadecagon, // 19
        Kratos_Octahedron,
        Kratos_Dodecahedron,
        Kratos_Icosahedron,
        Kratos_generic_type
    };


    ///@}
    ///@name Type Definitions
    ///@{

    /// Pointer definition of GeometryData
    KRATOS_CLASS_POINTER_DEFINITION( GeometryData );

    /** Type used for indexing in geometry data class. Unsigned int used for indexing
    point or integration point access methods and also all other
    methods which need point or integration point index.
    */
    typedef KRATOS_INDEX_TYPE IndexType;

    /** This typed used to return size or dimension in
    geometry data. Dimension, WorkingSpaceDimension, PointsNumber and
    ... return this type as their results.
    */
    typedef KRATOS_SIZE_TYPE SizeType;

    /** This type used for representing an integration point in
    geometry data. This integration point is a point with an
    additional weight component.
    */
    typedef IntegrationPoint<3> IntegrationPointType;

    /** A Vector of IntegrationPointType which used to hold
    integration points related to an integration
    method. IntegrationPoints functions used this type to return
    their results.
    */
    typedef std::vector<IntegrationPointType> IntegrationPointsArrayType;

    /** A Vector of IntegrationPointsArrayType which used to hold
    integration points related to different integration method
    implemented in geometry.
    */
    typedef boost::array<IntegrationPointsArrayType, static_cast<int>(IntegrationMethod::NumberOfIntegrationMethods)> IntegrationPointsContainerType;

    /** A third order tensor used as shape functions' values
    continer.
    */
    typedef boost::array<Matrix, static_cast<int>(IntegrationMethod::NumberOfIntegrationMethods)> ShapeFunctionsValuesContainerType;

    /** A fourth order tensor used as shape functions' local
    gradients container in geometry data.
    */
    typedef boost::array<boost::numeric::ublas::vector<Matrix>, static_cast<int>(IntegrationMethod::NumberOfIntegrationMethods)> ShapeFunctionsLocalGradientsContainerType;

    /** A third order tensor to hold shape functions'
    gradients. ShapefunctionsLocalGradients function return this
    type as its result.
    */
    typedef boost::numeric::ublas::vector<Matrix> ShapeFunctionsGradientsType;

    typedef boost::numeric::ublas::vector<Matrix> ShapeFunctionsSecondDerivativesType;

    /**
     * fourth order tensor to hold the third order derivatives of the
     * shape functions
     */
    typedef boost::numeric::ublas::vector<boost::numeric::ublas::vector<Matrix> > ShapeFunctionsThirdDerivativesType;

    typedef boost::array<boost::numeric::ublas::vector<Matrix>, static_cast<int>(IntegrationMethod::NumberOfIntegrationMethods)> MassFactorsContainerType;

    typedef boost::array<boost::numeric::ublas::vector<Matrix>, static_cast<int>(IntegrationMethod::NumberOfIntegrationMethods)> LampedMassFactorsContainerType;
    ///@}
    ///@name Life Cycle
    ///@{

//    GeometryData(SizeType ThisDimension = 0,
//          SizeType ThisWorkingSpaceDimension = 0,
//          SizeType ThisLocalSpaceDimension = 0)
//  : mDimension(ThisDimension)
//  , mWorkingSpaceDimension(ThisWorkingSpaceDimension)
//  , mLocalSpaceDimension(ThisLocalSpaceDimension)
//  , mDefaultMethod(GI_GAUSS_1)
//  , mIntegrationPoints()
//  , mShapeFunctionsValues()
//  , mShapeFunctionsLocalGradients()
//  {
//  }

    /** Complete argument constructor. This constructor gives a
    complete set of arguments to pass all the initial value of
    all the member variables of geometry class. Also it has
    default value for integration variables to make it usefull
    in the case of constructing new geometry without mapping and
    integrating properties.

    @param ThisDimension Dimension of this geometry.

    @param ThisWorkingSpaceDimension Working space dimension. for
    example a triangle 3d is a 2 dimensional shape but can be used
    in 3 dimensional space.

    @param ThisLocalSpaceDimension Local space dimension.
    for example a triangle is a 2 dimensional shape but
    can have 3 dimensional area coordinates l1, l2, l3.

    @param ThisDefaultMethod Default integration method. Its
    default value is gaussian integration with orden one which
    make no deference while in this condition there is no shape
    function database exist and integrating is not possible
    including by default method.

    @param ThisIntegrationPoints All the integration points in
    all methods. This is a Vector of IntegrationPointsArrayType
    and It must have at least four component correspounding to
    four integration method defined now. If there is some
    geometry which don't have all this method implemented
    related points Vector must exist but with zero size. For
    example if a geometry don't have gaussian orden one
    ThisIntegrationPoints[GI_GAUSS_1] must be an empty
    IntegrationPointsArrayType.

    @param ThisShapeFunctionsValues Values of all the shape
    functions evaluated in all integrations points of all
    integration methods. It's a three dimensional array \f$
    F_{ijk} \f$ where i = GI_GAUSS_1,..., GI_GAUSS_4 and j is
    the integration point index and k is the shape function
    index. In the other word component \f$ f_{ijk} \f$ is the
    value of the shape function related to node k evaluated in
    integration point j of i integration method point set. Again
    if there is some integration method unsupported an empty
    Matrix must assigned to related place. For example if a
    geometry don't have gaussian orden four
    ThisShapeFunctionsValues[GI_GAUSS_4] must be an empty
    Matrix.

    @param ThisShapeFunctionsLocalGradients Values of local
    gradients respected to all local coordinates of all the
    shape functions evaluated in all integrations points of all
    integration methods. It's a four dimensional array \f$
    F_{ijkh} \f$ where i = GI_GAUSS_1,..., GI_GAUSS_4 and j is
    the integration point index and k is the shape function
    index and h is local coordinate index. In the other word
    component \f$ f_{ijkh} \f$ is the value of h'th component of
    local gradient of the shape function related to node k
    evaluated in integration point j of i integration method
    point set. Again if there is some integration method
    unsupported an empty ShapeFunctionsGradientsType must
    assigned to related place. For example if a geometry don't
    have gaussian order "2" ThisShapeFunctionsValues[GI_GAUSS_2]
    must be an empty ShapeFunctionsGradientsType.
    */
    GeometryData( SizeType ThisDimension,
                  SizeType ThisWorkingSpaceDimension,
                  SizeType ThisLocalSpaceDimension,
                  enum IntegrationMethod ThisDefaultMethod,
                  const IntegrationPointsContainerType& ThisIntegrationPoints,
                  const ShapeFunctionsValuesContainerType& ThisShapeFunctionsValues,
                  const ShapeFunctionsLocalGradientsContainerType& ThisShapeFunctionsLocalGradients )
        : mDimension( ThisDimension )
        , mWorkingSpaceDimension( ThisWorkingSpaceDimension )
        , mLocalSpaceDimension( ThisLocalSpaceDimension )
        , mDefaultMethod( ThisDefaultMethod )
        , mIntegrationPoints( ThisIntegrationPoints )
        , mShapeFunctionsValues( ThisShapeFunctionsValues )
        , mShapeFunctionsLocalGradients( ThisShapeFunctionsLocalGradients )
    {
        for ( unsigned int i = 0 ; i < static_cast<int>(IntegrationMethod::NumberOfIntegrationMethods) ; i++ )
        {
            boost::numeric::ublas::vector<Matrix> temp( mIntegrationPoints[i].size() );

            for ( unsigned int j = 0 ; j < mIntegrationPoints[i].size() ; j++ )
                //temp[j]=outer_prod(mShapeFunctionsValues[i][j], mShapeFunctionsValues[i][j]);
                temp[j] = outer_prod( row( mShapeFunctionsValues[i], j ), row( mShapeFunctionsValues[i], j ) );

            mMassFactors[i] = temp;
        }
    }

    /** Copy constructor.
    Construct this geometry data as a copy of given geometry data.
    */
    GeometryData( const GeometryData& rOther )
        : mDimension( rOther.mDimension )
        , mWorkingSpaceDimension( rOther.mWorkingSpaceDimension )
        , mLocalSpaceDimension( rOther.mLocalSpaceDimension )
        , mDefaultMethod( rOther.mDefaultMethod )
        , mIntegrationPoints( rOther.mIntegrationPoints )
        , mShapeFunctionsValues( rOther.mShapeFunctionsValues )
        , mShapeFunctionsLocalGradients( rOther.mShapeFunctionsLocalGradients )
        , mMassFactors( rOther.mMassFactors )
    {
    }



    /// Destructor. Do nothing!!!
    virtual ~GeometryData() {}


    ///@}
    ///@name Operators
    ///@{

    /** Assignment operator.

    @note This operator don't copy the points and this
    geometry shares points with given source geometry. It's
    obvious that any change to this geometry's point affect
    source geometry's points too.

    @see Clone
    @see ClonePoints
    */
    GeometryData& operator=( const GeometryData& rOther )
    {
        mDimension = rOther.mDimension;
        mWorkingSpaceDimension = rOther.mWorkingSpaceDimension;
        mLocalSpaceDimension = rOther.mLocalSpaceDimension;
        mDefaultMethod = rOther.mDefaultMethod;
        mIntegrationPoints = rOther.mIntegrationPoints;
        mShapeFunctionsValues = rOther.mShapeFunctionsValues;
        mShapeFunctionsLocalGradients = rOther.mShapeFunctionsLocalGradients;
        mMassFactors = rOther.mMassFactors;

        return *this;
    }

    ///@}
    ///@name Informations
    ///@{

    /** Dimension of the geometry for example a triangle2d is a 2
    dimensional shape

    @return SizeType, dimension of this geometry.
    @see WorkingSpaceDimension()
    @see LocalSpaceDimension()
    */
    SizeType Dimension() const
    {
        return mDimension;
    }

    /** Working space dimension. for example a triangle is a 2
    dimensional shape but can be used in 3 dimensional space.

    @return SizeType, working space dimension of this geometry.
    @see Dimension()
    @see LocalSpaceDimension()
    */
    SizeType WorkingSpaceDimension() const
    {
        return mWorkingSpaceDimension;
    }

    /** Local space dimension. for example a triangle is a 2
    dimensional shape but can have 3 dimensional area
    coordinates l1, l2, l3.

    @return SizeType, local space dimension of this geometry.
    @see Dimension()
    @see WorkingSpaceDimension()
    */
    SizeType LocalSpaceDimension() const
    {
        return mLocalSpaceDimension;
    }


    ///@}
    ///@name Inquiry
    ///@{

    /** This method confirm you if this geometry has a specific
    integration method or not. This method will be usefull to
    control the geometry before intagrating using a specific
    method. In GeometryData class this method controls if the
    integration points vector respecting to this method is empty
    or not.

    @return bool true if this integration method exist and false if this
    method is not imeplemented for this geometry.
    */
    bool HasIntegrationMethod( IntegrationMethod ThisMethod ) const
    {
        return ( !mIntegrationPoints[static_cast<int>(ThisMethod)].empty() );
    }

    ///@}
    ///@name Integration
    ///@{

    /** Number of integtation points for default integration
    method. This method just call IntegrationPointsNumber(enum
    IntegrationMethod ThisMethod) with default integration
    method.

    @return SizeType which is the number of integration points
    for default integrating method.
    */

    IntegrationMethod DefaultIntegrationMethod() const
    {
        return mDefaultMethod;
    }

    SizeType IntegrationPointsNumber() const
    {
        return mIntegrationPoints[static_cast<int>(mDefaultMethod)].size();
    }

    /** Number of integtation points for given integration
    method. This method use integration points data base to
    obtain size of the integration points Vector respected to
    given method.

    @return SizeType which is the number of integration points
    for given integrating method.
    */
    SizeType IntegrationPointsNumber( IntegrationMethod ThisMethod ) const
    {
        return mIntegrationPoints[static_cast<int>(ThisMethod)].size();
    }


    /** Integtation points for default integration
    method. This method just call IntegrationPoints(enum
    IntegrationMethod ThisMethod) with default integration
    method.

    @return const IntegrationPointsArrayType which is Vector of integration points
    for default integrating method.
    */
    const IntegrationPointsArrayType& IntegrationPoints() const
    {
        return mIntegrationPoints[static_cast<int>(mDefaultMethod)];
    }

    /** Integtation points for given integration
    method. This method use integration points data base to
    obtain integration points Vector respected to
    given method.

    @return const IntegrationPointsArrayType which is Vector of integration points
    for default integrating method.
    */
    const IntegrationPointsArrayType& IntegrationPoints( IntegrationMethod ThisMethod ) const
    {
        return mIntegrationPoints[static_cast<int>(ThisMethod)];
    }

    ///@}
    ///@name Shape Function
    ///@{

    /** This method gives all shape functions values evaluated in all
    integration points of default integration method. It just
    call ShapeFunctionsValues(enum IntegrationMethod ThisMethod)
    with default integration method.There is no calculation and
    it just give it from shape functions values container.

    \note There is no control if the return matrix is empty or not!

    @return Matrix of values of shape functions \f$ F_{ij} \f$
    where i is the integration point index and j is the shape
    function index. In other word component \f$ f_{ij} \f$ is value
    of the shape function corresponding to node j evaluated in
    integration point i of default integration method.

    @see ShapeFunctionValue
    @see ShapeFunctionsLocalGradients
    @see ShapeFunctionLocalGradient
    */
    const Matrix& ShapeFunctionsValues() const
    {
        return mShapeFunctionsValues[static_cast<int>(mDefaultMethod)];
    }

    /** This method gives all shape functions values evaluated in all
    integration points of given integration method. There is no
    calculation and it just give it from shape functions values
    container.

    \note There is no control if the return matrix is empty or not!

    @param ThisMethod integration method which shape functions
    evaluated in its integration points.

    @return Matrix of values of shape functions \f$ F_{ij} \f$
    where i is the integration point index and j is the shape
    function index. In other word component \f$ f_{ij} \f$ is value
    of the shape function corresponding to node j evaluated in
    integration point i of given integration method.

    @see ShapeFunctionValue
    @see ShapeFunctionsLocalGradients
    @see ShapeFunctionLocalGradient
    */
    const Matrix& ShapeFunctionsValues( IntegrationMethod ThisMethod ) const
    {
        return mShapeFunctionsValues[static_cast<int>(ThisMethod)];
    }

    /** This method gives value of given shape function evaluated in
    given integration point of default integration method. It just
    call ShapeFunctionValue(IndexType IntegrationPointIndex,
    IndexType ShapeFunctionIndex, enum IntegrationMethod
    ThisMethod) with default integration method. There is no
    calculation and it just give it from shape functions values
    container if they are existing. Otherwise it gives you error
    which this value is not exist.

    @param IntegrationPointIndex index of integration point
    which shape functions evaluated in it.

    @param ShapeFunctionIndex index of node which correspounding
    shape function evaluated in given integration point.

    @return Value of given shape function in given integration
    point of default integration method.

    @see ShapeFunctionsValues
    @see ShapeFunctionsLocalGradients
    @see ShapeFunctionLocalGradient
    */
    double ShapeFunctionValue( IndexType IntegrationPointIndex, IndexType ShapeFunctionIndex ) const
    {
        if ( mShapeFunctionsValues[static_cast<int>(mDefaultMethod)].size1() <= IntegrationPointIndex )
            KRATOS_THROW_ERROR( std::invalid_argument, "No existing integration point", "" );

        if ( mShapeFunctionsValues[static_cast<int>(mDefaultMethod)].size2() <= ShapeFunctionIndex )
            KRATOS_THROW_ERROR( std::invalid_argument, "No existing shape function value", "" );

        return mShapeFunctionsValues[static_cast<int>(mDefaultMethod)]( IntegrationPointIndex, ShapeFunctionIndex );
    }

    /** This method gives value of given shape function evaluated in given
    integration point of given integration method. There is no
    calculation and it just give it from shape functions values
    container if they are existing. Otherwise it gives you error
    which this value is not exist.

    @param IntegrationPointIndex index of integration point
    which shape functions evaluated in it.

    @param ShapeFunctionIndex index of node which correspounding
    shape function evaluated in given integration point.

    @param ThisMethod integration method which shape function
    evaluated in its integration point.

    @return Value of given shape function in given integration
    point of given integration method.

    @see ShapeFunctionsValues
    @see ShapeFunctionsLocalGradients
    @see ShapeFunctionLocalGradient
    */
    double ShapeFunctionValue( IndexType IntegrationPointIndex, IndexType ShapeFunctionIndex, IntegrationMethod ThisMethod ) const
    {
        if ( mShapeFunctionsValues[static_cast<int>(ThisMethod)].size1() <= IntegrationPointIndex )
            KRATOS_THROW_ERROR( std::invalid_argument, "No existing integration point", "" );

        if ( mShapeFunctionsValues[static_cast<int>(ThisMethod)].size2() <= ShapeFunctionIndex )
            KRATOS_THROW_ERROR( std::invalid_argument, "No existing shape function value", "" );

        return mShapeFunctionsValues[static_cast<int>(ThisMethod)]( IntegrationPointIndex, ShapeFunctionIndex );
    }

    /** This method gives all shape functions gradients evaluated in all
    integration points of default integration method. It just
    call ShapeFunctionsLocalGradients(enum IntegrationMethod ThisMethod)
    with default integration method. There is no calculation and
    it just give it from shape functions values container.

    \note There is no control if there is any gradient calculated or not!

    @return shape functions' gradients \f$ F_{ijk} \f$ where i
    is the integration point index and j is the shape function
    index and k is local coordinate index. In other word
    component \f$ f_{ijk} \f$ is k'th component of gradient of
    the shape function corresponding to node j evaluated in
    integration point i of default integration method.

    @see ShapeFunctionsValues
    @see ShapeFunctionValue
    @see ShapeFunctionLocalGradient
    */
    const ShapeFunctionsGradientsType& ShapeFunctionsLocalGradients() const
    {
        return mShapeFunctionsLocalGradients[static_cast<int>(mDefaultMethod)];
    }

    /** This method gives all shape functions gradients evaluated in
    all integration points of given integration method. There is
    no calculation and it just give it from shape functions
    values container.

    \note There is no control if there is any gradient calculated or not!

    @param ThisMethod integration method which shape functions
    gradients evaluated in its integration points.

    @return shape functions' gradients \f$ F_{ijk} \f$ where i
    is the integration point index and j is the shape function
    index and k is local coordinate index. In other word
    component \f$ f_{ijk} \f$ is k'th component of gradient of
    the shape function corresponding to node j evaluated in
    integration point i of given integration method.

    @see ShapeFunctionsValues
    @see ShapeFunctionValue
    @see ShapeFunctionLocalGradient
    */
    const ShapeFunctionsGradientsType& ShapeFunctionsLocalGradients( IntegrationMethod ThisMethod ) const
    {
        return mShapeFunctionsLocalGradients[static_cast<int>(ThisMethod)];
    }

    /** This method gives gradient of given shape function evaluated in
    given integration point of default integration method. It just
    call ShapeFunctionLocalGradient(IndexType IntegrationPointIndex,
    IndexType ShapeFunctionIndex, enum IntegrationMethod
    ThisMethod) with default integration method. There is no
    calculation and it just give it from shape functions values
    container if they are existing. Otherwise it gives you error
    which this value is not exist.

    @param IntegrationPointIndex index of integration point
    which shape function gradient evaluated in it.

    @param ShapeFunctionIndex index of node which correspounding
    shape function gradient evaluated in given integration point.

    @return Gradient of given shape function in given integration
    point of default integration method.

    @see ShapeFunctionsValues
    @see ShapeFunctionValue
    @see ShapeFunctionsLocalGradients
    */
    const Matrix& ShapeFunctionLocalGradient( IndexType IntegrationPointIndex ) const
    {
        if ( mShapeFunctionsLocalGradients[static_cast<int>(mDefaultMethod)].size() <= IntegrationPointIndex )
            KRATOS_THROW_ERROR( std::invalid_argument, "No existing integration point", "" );

        return mShapeFunctionsLocalGradients[static_cast<int>(mDefaultMethod)][IntegrationPointIndex];
    }

    /** This method gives gradient of given shape function evaluated
    in given integration point of given integration
    method. There is no calculation and it just give it from
    shape functions values container if they are
    existing. Otherwise it gives you error which this value is
    not exist.

    @param IntegrationPointIndex index of integration point
    which shape function gradient evaluated in it.

    @param ShapeFunctionIndex index of node which correspounding
    shape function gradient evaluated in given integration point.

    @param ThisMethod integration method which shape function gradient
    evaluated in its integration points.

    @return Gradient of given shape function in given integration
    point of given integration method.

    @see ShapeFunctionsValues
    @see ShapeFunctionValue
    @see ShapeFunctionsLocalGradients
    */
    const Matrix& ShapeFunctionLocalGradient( IndexType IntegrationPointIndex, IndexType ShapeFunctionIndex, IntegrationMethod ThisMethod ) const
    {
        if ( mShapeFunctionsLocalGradients[static_cast<int>(ThisMethod)].size() <= IntegrationPointIndex )
            KRATOS_THROW_ERROR( std::invalid_argument, "No existing integration point", "" );

        return mShapeFunctionsLocalGradients[static_cast<int>(ThisMethod)][IntegrationPointIndex];
    }

    boost::numeric::ublas::vector<Matrix> const& MassFactors() const
    {
        return mMassFactors[static_cast<int>(mDefaultMethod)];
    }

    boost::numeric::ublas::vector<Matrix> const& MassFactors( IntegrationMethod ThisMethod ) const
    {
        return mMassFactors[static_cast<int>(ThisMethod)];
    }

    Matrix const& MassFactors( IndexType IntegrationPointIndex ) const
    {
        return mMassFactors[static_cast<int>(mDefaultMethod)][IntegrationPointIndex];
    }

    Matrix const& MassFactors( IndexType IntegrationPointIndex, IntegrationMethod ThisMethod ) const
    {
        return mMassFactors[static_cast<int>(ThisMethod)][IntegrationPointIndex];
    }

    ///@}
    ///@name Input and output
    ///@{

    /** Turn back information as a string.

    @return String contains information about this geometry.
    @see PrintData()
    @see PrintInfo()
    */
    virtual std::string Info() const
    {
        return "geometry data";
    }

    /** Print information about this object.

    @param rOStream Stream to print into it.
    @see PrintData()
    @see Info()
    */
    virtual void PrintInfo( std::ostream& rOStream ) const
    {
        rOStream << "geometry data";
    }

    /** Print geometry's data into given stream. Prints it's points
    by the order they stored in the geometry and then center
    point of geometry.

    @param rOStream Stream to print into it.
    @see PrintInfo()
    @see Info()
    */
    virtual void PrintData( std::ostream& rOStream ) const
    {
        rOStream << "    Dimension               : " << mDimension << std::endl;
        rOStream << "    working space dimension : " << mWorkingSpaceDimension << std::endl;
        rOStream << "    Local space dimension   : " << mLocalSpaceDimension;
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    SizeType mDimension;

    SizeType mWorkingSpaceDimension;

    SizeType mLocalSpaceDimension;

    IntegrationMethod mDefaultMethod;

    IntegrationPointsContainerType mIntegrationPoints;

    ShapeFunctionsValuesContainerType mShapeFunctionsValues;

    ShapeFunctionsLocalGradientsContainerType mShapeFunctionsLocalGradients;

    MassFactorsContainerType mMassFactors;

    LampedMassFactorsContainerType mLampedMassFactors;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save( Serializer& rSerializer ) const
    {
        rSerializer.save( "Dimension", mDimension );
        rSerializer.save( "Working Space Dimension", mWorkingSpaceDimension );
        rSerializer.save( "Local Space Dimension", mLocalSpaceDimension );
//    rSerializer.save("Default Method", mDefaultMethod);
//    rSerializer.save("Integration Points", mIntegrationPoints);
//    rSerializer.save("Shape Functions Values", mShapeFunctionsValues);
//    rSerializer.save("Shape Functions Local Gradients", mShapeFunctionsLocalGradients);
//    rSerializer.save("Mass Factors", mMassFactors);
//    rSerializer.save("Lamped Mass Factors", mLampedMassFactors);
    }

    virtual void load( Serializer& rSerializer )
    {
        rSerializer.load( "Dimension", mDimension );
        rSerializer.load( "Working Space Dimension", mWorkingSpaceDimension );
        rSerializer.load( "Local Space Dimension", mLocalSpaceDimension );
//    rSerializer.load("Default Method", mDefaultMethod);
//    rSerializer.load("Integration Points", mIntegrationPoints);
//    rSerializer.load("Shape Functions Values", mShapeFunctionsValues);
//    rSerializer.load("Shape Functions Local Gradients", mShapeFunctionsLocalGradients);
//    rSerializer.load("Mass Factors", mMassFactors);
//    rSerializer.load("Lamped Mass Factors", mLampedMassFactors);
    }

    // Private default constructor for serialization
    GeometryData()
    {
    }


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Private Friends
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}

}; // Class GeometryData

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> ( std::istream& rIStream,
                                   GeometryData& rThis )
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << ( std::ostream& rOStream,
                                   const GeometryData& rThis )
{
    rThis.PrintInfo( rOStream );
    rOStream << std::endl;
    rThis.PrintData( rOStream );

    return rOStream;
}

/// output stream function for IntegrationMethod
inline std::ostream& operator << ( std::ostream& rOStream,
                                   const GeometryData::IntegrationMethod Value )
{
    switch(Value)
    {
        case GeometryData::IntegrationMethod::GI_GAUSS_1:           rOStream << "GI_GAUSS_1"; break;
        case GeometryData::IntegrationMethod::GI_GAUSS_2:           rOStream << "GI_GAUSS_2"; break;
        case GeometryData::IntegrationMethod::GI_GAUSS_3:           rOStream << "GI_GAUSS_3"; break;
        case GeometryData::IntegrationMethod::GI_GAUSS_4:           rOStream << "GI_GAUSS_4"; break;
        case GeometryData::IntegrationMethod::GI_GAUSS_5:           rOStream << "GI_GAUSS_5"; break;
        case GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_1:  rOStream << "GI_EXTENDED_GAUSS_1"; break;
        case GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_2:  rOStream << "GI_EXTENDED_GAUSS_2"; break;
        case GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_3:  rOStream << "GI_EXTENDED_GAUSS_3"; break;
        case GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_4:  rOStream << "GI_EXTENDED_GAUSS_4"; break;
        case GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_5:  rOStream << "GI_EXTENDED_GAUSS_5"; break;
        default:                                                    rOStream << "GI_Unknown"; break;
    }
    return rOStream;
}

/// output stream function for KratosGeometryFamily
inline std::ostream& operator << ( std::ostream& rOStream,
                                   const GeometryData::KratosGeometryFamily Value )
{
    switch(Value)
    {
        case GeometryData::KratosGeometryFamily::Kratos_NoElement:      rOStream << "NoElement";    break;
        case GeometryData::KratosGeometryFamily::Kratos_Point:          rOStream << "Point";        break;
        case GeometryData::KratosGeometryFamily::Kratos_Linear:         rOStream << "Linear";       break;
        case GeometryData::KratosGeometryFamily::Kratos_Triangle:       rOStream << "Triangle";     break;
        case GeometryData::KratosGeometryFamily::Kratos_Quadrilateral:  rOStream << "Quadrilateral";    break;
        case GeometryData::KratosGeometryFamily::Kratos_Tetrahedra:     rOStream << "Tetrahedra";       break;
        case GeometryData::KratosGeometryFamily::Kratos_Hexahedra:      rOStream << "Hexahedra";        break;
        case GeometryData::KratosGeometryFamily::Kratos_Prism:          rOStream << "Prism";            break;
        case GeometryData::KratosGeometryFamily::Kratos_NURBS:          rOStream << "NURBS";            break;
        case GeometryData::KratosGeometryFamily::Kratos_Polygon:        rOStream << "Polygon";          break;
        case GeometryData::KratosGeometryFamily::Kratos_Polyhedron:     rOStream << "Polyhedron";       break;
        default:                                                        rOStream << "generic";  break;
    }
    return rOStream;
}

/// output stream function for KratosGeometryType
inline std::ostream& operator << ( std::ostream& rOStream,
                                   const GeometryData::KratosGeometryType Value)
{
    switch (Value)
    {
        case GeometryData::KratosGeometryType::Kratos_Hexahedra3D20:      rOStream << "Hexahedra3D20";    break;
        case GeometryData::KratosGeometryType::Kratos_Hexahedra3D27:      rOStream << "Hexahedra3D27";    break;
        case GeometryData::KratosGeometryType::Kratos_Hexahedra3D8:       rOStream << "Hexahedra3D8";     break;
        case GeometryData::KratosGeometryType::Kratos_Prism3D15:          rOStream << "Prism3D15";        break;
        case GeometryData::KratosGeometryType::Kratos_Prism3D6:           rOStream << "Prism3D6";         break;
        case GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4:   rOStream << "Quadrilateral2D4"; break;
        case GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8:   rOStream << "Quadrilateral2D8"; break;
        case GeometryData::KratosGeometryType::Kratos_Quadrilateral2D9:   rOStream << "Quadrilateral2D9"; break;
        case GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4:   rOStream << "Quadrilateral3D4"; break;
        case GeometryData::KratosGeometryType::Kratos_Quadrilateral3D8:   rOStream << "Quadrilateral3D8"; break;
        case GeometryData::KratosGeometryType::Kratos_Quadrilateral3D9:   rOStream << "Quadrilateral3D9"; break;
        case GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10:     rOStream << "Tetrahedra3D10";   break;
        case GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4:      rOStream << "Tetrahedra3D4";    break;
        case GeometryData::KratosGeometryType::Kratos_Triangle2D3:        rOStream << "Triangle2D3";      break;
        case GeometryData::KratosGeometryType::Kratos_Triangle2D6:        rOStream << "Triangle2D6";      break;
        case GeometryData::KratosGeometryType::Kratos_Triangle3D3:        rOStream << "Triangle3D3";      break;
        case GeometryData::KratosGeometryType::Kratos_Triangle3D6:        rOStream << "Triangle3D6";      break;
        case GeometryData::KratosGeometryType::Kratos_Line2D2:            rOStream << "Line2D2";          break;
        case GeometryData::KratosGeometryType::Kratos_Line2D3:            rOStream << "Line2D3";          break;
        case GeometryData::KratosGeometryType::Kratos_Line3D2:            rOStream << "Line3D2";          break;
        case GeometryData::KratosGeometryType::Kratos_Line3D3:            rOStream << "Line3D3";          break;
        case GeometryData::KratosGeometryType::Kratos_Point2D:            rOStream << "Point2D";          break;
        case GeometryData::KratosGeometryType::Kratos_Point3D:            rOStream << "Point3D";          break;
        case GeometryData::KratosGeometryType::Kratos_Sphere3D1:          rOStream << "Sphere3D1";        break;
        case GeometryData::KratosGeometryType::Kratos_Bezier1D:           rOStream << "Bezier1D";         break;
        case GeometryData::KratosGeometryType::Kratos_Bezier1D2:          rOStream << "Bezier1D2";        break;
        case GeometryData::KratosGeometryType::Kratos_Bezier1D3:          rOStream << "Bezier1D3";        break;
        case GeometryData::KratosGeometryType::Kratos_Bezier2D:           rOStream << "Bezier2D";         break;
        case GeometryData::KratosGeometryType::Kratos_Bezier2D3:          rOStream << "Bezier2D3";        break;
        case GeometryData::KratosGeometryType::Kratos_Bezier3D:           rOStream << "Bezier3D";         break;
        case GeometryData::KratosGeometryType::Kratos_Tritagon:           rOStream << "Tritagon";         break;
        case GeometryData::KratosGeometryType::Kratos_Tetragon:           rOStream << "Tetragon";         break;
        case GeometryData::KratosGeometryType::Kratos_Pentagon:           rOStream << "Pentagon";         break;
        case GeometryData::KratosGeometryType::Kratos_Hexagon:            rOStream << "Hexagon";          break;
        case GeometryData::KratosGeometryType::Kratos_Heptagon:           rOStream << "Heptagon";         break;
        case GeometryData::KratosGeometryType::Kratos_Octagon:            rOStream << "Octagon";          break;
        case GeometryData::KratosGeometryType::Kratos_Nonagon:            rOStream << "Nonagon";          break;
        case GeometryData::KratosGeometryType::Kratos_Decagon:            rOStream << "Decagon";          break;
        case GeometryData::KratosGeometryType::Kratos_Hendecagon:         rOStream << "Hendecagon";       break;
        case GeometryData::KratosGeometryType::Kratos_Dodecagon:          rOStream << "Dodecagon";        break;
        case GeometryData::KratosGeometryType::Kratos_Triskaidecagon:     rOStream << "Triskaidecagon";   break;
        case GeometryData::KratosGeometryType::Kratos_Tetrakaidecagon:    rOStream << "Tetrakaidecagon";  break;
        case GeometryData::KratosGeometryType::Kratos_Pentadecagon:       rOStream << "Pentadecagon";     break;
        case GeometryData::KratosGeometryType::Kratos_Hexakaidecagon:     rOStream << "Hexakaidecagon";   break;
        case GeometryData::KratosGeometryType::Kratos_Heptadecagon:       rOStream << "Heptadecagon";     break;
        case GeometryData::KratosGeometryType::Kratos_Octakaidecagon:     rOStream << "Octakaidecagon";   break;
        case GeometryData::KratosGeometryType::Kratos_Enneadecagon:       rOStream << "Enneadecagon";     break;
        case GeometryData::KratosGeometryType::Kratos_Octahedron:         rOStream << "Octahedron";       break;
        case GeometryData::KratosGeometryType::Kratos_Dodecahedron:       rOStream << "Dodecahedron";     break;
        case GeometryData::KratosGeometryType::Kratos_Icosahedron:        rOStream << "Icosahedron";      break;
        default:                                                          rOStream << "generic";          break;
    }
    return rOStream;
}

///@}

}  // namespace Kratos.

#endif // KRATOS_GEOMETRY_DATA_H_INCLUDED  defined
