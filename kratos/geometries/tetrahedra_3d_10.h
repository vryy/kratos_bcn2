//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Janosch Stascheit
//                   Felix Nagel
//  contributors:    Hoang Giang Bui
//                   Josep Maria Carbonell
//



#if !defined(KRATOS_TETRAHEDRA_3D_10_H_INCLUDED )
#define  KRATOS_TETRAHEDRA_3D_10_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geometries/triangle_3d_6.h"
#include "integration/tetrahedron_gauss_legendre_integration_points.h"


namespace Kratos
{

/**
 * An eight node hexahedra geometry with linear shape functions
 */
template<class TPointType> class Tetrahedra3D10 : public Geometry<TPointType>
{
public:
    /**
     * Type Definitions
     */

    /**
     * Geometry as base class.
     */
    typedef Geometry<TPointType> BaseType;

    /**
     * Type of edge and face geometries
     */
    typedef Line3D3<TPointType> EdgeType;
    typedef Triangle3D6<TPointType> FaceType;

    /**
     * Pointer definition of Tetrahedra3D10
     */
    KRATOS_CLASS_POINTER_DEFINITION( Tetrahedra3D10 );

    /**
     * Integration methods implemented in geometry.
     */
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    /**
     * A Vector of counted pointers to Geometries. Used for
     * returning edges of the geometry.
     */
    typedef typename BaseType::GeometriesArrayType GeometriesArrayType;

    /**
     * Type used for coordinate value.
     */
    typedef typename BaseType::DataType DataType;

    /**
     * Type used for shape function value.
     */
    typedef typename BaseType::ValueType ValueType;

    /**
     * Type used for indexing in geometry class.std::size_t used for indexing
     * point or integration point access methods and also all other
     * methods which need point or integration point index.
     */
    typedef typename BaseType::IndexType IndexType;

    /**
     * This typed used to return size or dimension in
     * geometry. Dimension, WorkingDimension, PointsNumber and
     * ... return this type as their results.
     */
    typedef typename BaseType::SizeType SizeType;

    /**
     * Redefinition of template parameter TPointType.
     */
    typedef typename BaseType::PointType PointType;

    /**
     * Array of counted pointers to point. This type used to hold
     * geometry's points.
     */
    typedef typename BaseType::PointsArrayType PointsArrayType;

    /**
     * This type used for representing an integration point in
     * geometry. This integration point is a point with an
     * additional weight component.
     */
    typedef typename BaseType::IntegrationPointType IntegrationPointType;

    /**
     * A Vector of IntegrationPointType which used to hold
     * integration points related to an integration
     * method. IntegrationPoints functions used this type to return
     * their results.
     */
    typedef typename BaseType::IntegrationPointsArrayType IntegrationPointsArrayType;

    /**
     * A Vector of IntegrationPointsArrayType which used to hold
     * integration points related to different integration method
     * implemented in geometry.
     */
    typedef typename BaseType::IntegrationPointsContainerType IntegrationPointsContainerType;

    /**
     * A third order tensor used as shape functions' values
     * container.
     */
    typedef typename BaseType::ShapeFunctionsValuesContainerType ShapeFunctionsValuesContainerType;

    /**
     * A fourth order tensor used as shape functions' local
     * gradients container in geometry.
     */
    typedef typename BaseType::ShapeFunctionsLocalGradientsContainerType ShapeFunctionsLocalGradientsContainerType;

    /**
     * A third order tensor to hold jacobian matrices evaluated at
     * integration points. Jacobian and InverseOfJacobian functions
     * return this type as their result.
     */
    typedef typename BaseType::JacobiansType JacobiansType;

    /**
     * A third order tensor to hold shape functions' local gradients at all integration points.
     * ShapefunctionsLocalGradients function return this
     * type as its result.
     */
    typedef typename BaseType::ShapeFunctionsGradientsType ShapeFunctionsGradientsType;

    /**
     * A third order tensor to hold shape functions' local second derivatives at a point.
     * ShapeFunctionsSecondDerivatives function return this
     * type as its result.
     */
    typedef typename BaseType::ShapeFunctionsSecondDerivativesType ShapeFunctionsSecondDerivativesType;

    /**
    * A fourth order tensor to hold shape functions' local third derivatives at a point.
    * ShapeFunctionsThirdDerivatives function return this
    * type as its result.
    */
    typedef typename BaseType::ShapeFunctionsThirdDerivativesType ShapeFunctionsThirdDerivativesType;

    /**
     * Type of the normal vector used for normal to edges in geomety.
     */
    typedef typename BaseType::NormalType NormalType;

    /**
     * Type of coordinates array
     */
    typedef typename BaseType::CoordinatesArrayType CoordinatesArrayType;

    /** This type used for representing the local coordinates of
    an integration point
    */
    typedef typename BaseType::LocalCoordinatesArrayType LocalCoordinatesArrayType;

    /**
     * Type of Matrix
     */
    typedef typename BaseType::MatrixType MatrixType;
    typedef typename BaseType::ZeroMatrixType ZeroMatrixType;

    /**
     * Type of Vector
     */
    typedef typename BaseType::VectorType VectorType;
    typedef typename BaseType::ZeroVectorType ZeroVectorType;

    /**
     * Life Cycle
     */

    Tetrahedra3D10( const PointType& Point1, const PointType& Point2,
                    const PointType& Point3, const PointType& Point4 ,
                    const PointType& Point5, const PointType& Point6,
                    const PointType& Point7, const PointType& Point8 ,
                    const PointType& Point9, const PointType& Point10
                  )
        : BaseType( PointsArrayType(), &msGeometryData )
    {
        this->Points().reserve( 10 );
        this->Points().push_back( typename PointType::Pointer( new PointType( Point1 ) ) );
        this->Points().push_back( typename PointType::Pointer( new PointType( Point2 ) ) );
        this->Points().push_back( typename PointType::Pointer( new PointType( Point3 ) ) );
        this->Points().push_back( typename PointType::Pointer( new PointType( Point4 ) ) );
        this->Points().push_back( typename PointType::Pointer( new PointType( Point5 ) ) );
        this->Points().push_back( typename PointType::Pointer( new PointType( Point6 ) ) );
        this->Points().push_back( typename PointType::Pointer( new PointType( Point7 ) ) );
        this->Points().push_back( typename PointType::Pointer( new PointType( Point8 ) ) );
        this->Points().push_back( typename PointType::Pointer( new PointType( Point9 ) ) );
        this->Points().push_back( typename PointType::Pointer( new PointType( Point10 ) ) );
    }

    Tetrahedra3D10( typename PointType::Pointer pPoint1,
                    typename PointType::Pointer pPoint2,
                    typename PointType::Pointer pPoint3,
                    typename PointType::Pointer pPoint4,
                    typename PointType::Pointer pPoint5,
                    typename PointType::Pointer pPoint6,
                    typename PointType::Pointer pPoint7,
                    typename PointType::Pointer pPoint8,
                    typename PointType::Pointer pPoint9,
                    typename PointType::Pointer pPoint10
                  )
        : BaseType( PointsArrayType(), &msGeometryData )
    {
        this->Points().reserve( 10 );
        this->Points().push_back( pPoint1 );
        this->Points().push_back( pPoint2 );
        this->Points().push_back( pPoint3 );
        this->Points().push_back( pPoint4 );
        this->Points().push_back( pPoint5 );
        this->Points().push_back( pPoint6 );
        this->Points().push_back( pPoint7 );
        this->Points().push_back( pPoint8 );
        this->Points().push_back( pPoint9 );
        this->Points().push_back( pPoint10 );
    }

    Tetrahedra3D10( const PointsArrayType& ThisPoints )
        : BaseType( ThisPoints, &msGeometryData )
    {
        if ( this->PointsNumber() != 10 )
            KRATOS_THROW_ERROR( std::invalid_argument,
                                "Invalid points number. Expected 10, given " ,
                                this->PointsNumber() );
    }

    /**
     * Copy constructor.
     * Construct this geometry as a copy of given geometry.
     *
     * @note This copy constructor don't copy the points and new
     * geometry shares points with given source geometry. It's
     * obvious that any change to this new geometry's point affect
     * source geometry's points too.
     */
    Tetrahedra3D10( Tetrahedra3D10 const& rOther )
        : BaseType( rOther )
    {
    }

    /**
     * Copy constructor from a geometry with other point type.
     * Construct this geometry as a copy of given geometry which
     * has different type of points. The given goemetry's
     * TOtherPointType* must be implicity convertible to this
     * geometry PointType.
     *
     * @note This copy constructor don't copy the points and new
     * geometry shares points with given source geometry. It's
     * obvious that any change to this new geometry's point affect
     * source geometry's points too.
     */
    template<class TOtherPointType> Tetrahedra3D10( Tetrahedra3D10<TOtherPointType> const& rOther )
        : BaseType( rOther )
    {
    }

    /// Destructor. Does nothing!!!
    ~Tetrahedra3D10() override {}

    GeometryData::KratosGeometryFamily GetGeometryFamily() const final
    {
        return GeometryData::KratosGeometryFamily::Kratos_Tetrahedra;
    }

    GeometryData::KratosGeometryType GetGeometryType() const final
    {
        return GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10;
    }

    /**
     * Operators
     */

    /**
     * Assignment operator.
     *
     * @note This operator don't copy the points and this
     * geometry shares points with given source geometry. It's
     * obvious that any change to this geometry's point affect
     * source geometry's points too.
     *
     * @see Clone
     * @see ClonePoints
     */
    Tetrahedra3D10& operator=( const Tetrahedra3D10& rOther )
    {
        BaseType::operator=( rOther );
        return *this;
    }

    /**
     * Assignment operator for geometries with different point type.
     *
     * @note This operator don't copy the points and this
     * geometry shares points with given source geometry. It's
     * obvious that any change to this geometry's point affect
     * source geometry's points too.
     *
     * @see Clone
     * @see ClonePoints
     */
    template<class TOtherPointType>
    Tetrahedra3D10& operator=( Tetrahedra3D10<TOtherPointType> const & rOther )
    {
        BaseType::operator=( rOther );

        return *this;
    }


    /**
     * Operations
     */

    typename BaseType::Pointer Create( PointsArrayType const& ThisPoints ) const override
    {
        return typename BaseType::Pointer( new Tetrahedra3D10( ThisPoints ) );
    }

    typename Geometry< Point<3, DataType> >::Pointer Clone() const override
    {
        typename Geometry< Point<3, DataType> >::PointsArrayType NewPoints;

        //making a copy of the nodes TO POINTS (not Nodes!!!)
        for ( IndexType i = 0 ; i < this->Points().size() ; i++ )
            NewPoints.push_back( this->Points()[i] );

        //creating a geometry with the new points
        typename Geometry< Point<3, DataType> >::Pointer p_clone( new Tetrahedra3D10< Point<3, DataType> >( NewPoints ) );

        p_clone->ClonePoints();

        return p_clone;
    }

    //lumping factors for the calculation of the lumped mass matrix
    Vector& LumpingFactors( Vector& rResult ) const override
    {
        if(rResult.size() != 10)
            rResult.resize( 10, false );
        std::fill( rResult.begin(), rResult.end(), 1.00 / 10.00 );
        return rResult;
    }


    /**
     * Informations
     */

    /**
     * This method calculates and returns Length or charactereistic
     * length of this geometry depending on it's dimension. For one
     * dimensional geometry for example Line it returns length of it
     * and for the other geometries it gives Characteristic length
     * otherwise.
     *
     * @return DataType value contains length or Characteristic
     * length
     * @see Area()
     * @see Volume()
     * @see DomainSize()
     *
     * :TODO: might be necessary to reimplement
     */
    DataType Length() const override
    {
        return std::sqrt( std::abs( this->DeterminantOfJacobian( LocalCoordinatesArrayType() ) ) );
    }

    /**
     * This method calculates and returns area or surface area of
     * this geometry depending to it's dimension. For one dimensional
     * geometry it returns zero, for two dimensional it gives area
     * and for three dimensional geometries it gives surface area.
     *
     *
     * @return DataType value contains area or surface area.
     * @see Length()
     * @see Volume()
     * @see DomainSize()
     *
     * :TODO: might be necessary to reimplement
     */
    DataType Area() const override
    {
        return Volume();
    }

    DataType Volume() const override //Not a closed formula for a quadratic tetrahedra
    {
        VectorType temp;
        this->DeterminantOfJacobian( temp, msGeometryData.DefaultIntegrationMethod() );
        const IntegrationPointsArrayType& integration_points = this->IntegrationPoints( msGeometryData.DefaultIntegrationMethod() );

        DataType Volume = 0.00;
        for ( unsigned int i = 0; i < integration_points.size(); i++ )
        {
            Volume += temp[i] * integration_points[i].Weight();
        }

        return Volume;
    }

    /**
     * This method calculate and return length, area or volume of
     * this geometry depending to it's dimension. For one dimensional
     * geometry it returns its length, for two dimensional it gives area
     * and for three dimensional geometries it gives its volume.
     *
     * @return DataType value contains length, area or volume.
     * @see Length()
     * @see Area()
     * @see Volume()
     *
     * :TODO: might be necessary to reimplement
     */
    DataType DomainSize() const override
    {
        return  Volume();
    }

    /**
     * Returns whether local arbitrary point is inside the Geometry
     */
    bool IsInside( const LocalCoordinatesArrayType& rPoint ) const override
    {
        if( rPoint[0] >= 0.0-1.0e-8 && rPoint[0] <= 1.0+1.0e-8 )
            if( rPoint[1] >= 0.0-1.0e-8 && rPoint[1] <= 1.0 +1.0e-8)
                if( rPoint[2] >= 0.0-1.0e-8 && rPoint[2] <= 1.0+1.0e-8 )
                    if( ((1.0-(rPoint[0] + rPoint[1] + rPoint[2])) >= 0.0-1.0e-8)&&((1.0-(rPoint[0] + rPoint[1] + rPoint[2])) <= 1.0+1.0e-8) )
                        return true;
        return false;
    }


    /** This method gives you number of all edges of this
    geometry.
    @return SizeType containes number of this geometry edges.
    @see Edges()
    @see Edge()
     */
    // will be used by refinement algorithm, thus uncommented. janosch.
    SizeType EdgesNumber() const override
    {
        return 6;
    }

    SizeType FacesNumber() const override
    {
        return 4;
    }

    /** This method gives you all edges of this geometry.

    @return GeometriesArrayType containes this geometry edges.
    @see EdgesNumber()
    @see Edge()
    */
    GeometriesArrayType Edges( void ) const override
    {
        GeometriesArrayType edges = GeometriesArrayType();
        typedef typename Geometry<TPointType>::Pointer EdgePointerType;
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 0 ),
                                              this->pGetPoint( 4 ),
                                              this->pGetPoint( 1 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 1 ),
                                              this->pGetPoint( 5 ),
                                              this->pGetPoint( 2 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 2 ),
                                              this->pGetPoint( 6 ),
                                              this->pGetPoint( 0 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 0 ),
                                              this->pGetPoint( 7 ),
                                              this->pGetPoint( 3 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 1 ),
                                              this->pGetPoint( 8 ),
                                              this->pGetPoint( 3 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 2 ),
                                              this->pGetPoint( 9 ),
                                              this->pGetPoint( 3 ) ) ) );

        return edges;
    }

    GeometriesArrayType Faces( void ) const override
    {
        GeometriesArrayType faces = GeometriesArrayType();
        typedef typename Geometry<TPointType>::Pointer FacePointerType;
        faces.push_back( FacePointerType( new FaceType(
                                              this->pGetPoint( 0 ),
                                              this->pGetPoint( 2 ),
                                              this->pGetPoint( 1 ),
                                              this->pGetPoint( 6 ),
                                              this->pGetPoint( 5 ),
                                              this->pGetPoint( 4 ) ) ) );
        faces.push_back( FacePointerType( new FaceType(
                                              this->pGetPoint( 0 ),
                                              this->pGetPoint( 3 ),
                                              this->pGetPoint( 2 ),
                                              this->pGetPoint( 7 ),
                                              this->pGetPoint( 9 ),
                                              this->pGetPoint( 6 ) ) ) );
        faces.push_back( FacePointerType( new FaceType(
                                              this->pGetPoint( 0 ),
                                              this->pGetPoint( 1 ),
                                              this->pGetPoint( 3 ),
                                              this->pGetPoint( 4 ),
                                              this->pGetPoint( 8 ),
                                              this->pGetPoint( 7 ) ) ) );
        faces.push_back( FacePointerType( new FaceType(
                                              this->pGetPoint( 2 ),
                                              this->pGetPoint( 3 ),
                                              this->pGetPoint( 1 ),
                                              this->pGetPoint( 9 ),
                                              this->pGetPoint( 8 ),
                                              this->pGetPoint( 5 ) ) ) );
        return faces;
    }

    /**
     * Shape Function
     */

    /**
     * Calculates the value of a given shape function at a given point.
     *
     * @param ShapeFunctionIndex The number of the desired shape function
     * @param rPoint the given point in local coordinates at which the
     * value of the shape function is calculated
     *
     * @return the value of the shape function at the given point
     * TODO: TO BE VERIFIED
     * REF: Jacob Fish, A First Course in Finite Elements
     */
    ValueType ShapeFunctionValue( IndexType ShapeFunctionIndex, const LocalCoordinatesArrayType& rPoint ) const override
    {
        const ValueType zeta1 = rPoint[0];
        const ValueType zeta2 = rPoint[1];
        const ValueType zeta3 = rPoint[2];
        const ValueType zeta0 = 1. - zeta1 - zeta2 - zeta3;

        switch ( ShapeFunctionIndex )
        {
        case 0:
            return zeta0*(2.*zeta0 - 1.);
        case 1:
            return zeta1*(2.*zeta1 - 1.);
        case 2:
            return zeta2*(2.*zeta2 - 1.);
        case 3:
            return zeta3*(2.*zeta3 - 1.);
        case 4:
            return 4.*zeta0*zeta1;
        case 5:
            return 4.*zeta1*zeta2;
        case 6:
            return 4.*zeta2*zeta0;
        case 7:
            return 4.*zeta0*zeta3;
        case 8:
            return 4.*zeta1*zeta3;
        case 9:
            return 4.*zeta2*zeta3;
        default:
            KRATOS_THROW_ERROR( std::logic_error,
                                "Wrong index of shape function!" , *this );
        }

        return 0;
    }

    /** This method gives gradient of all shape functions evaluated
     * in given point.
     * There is no calculation and it just give it from
     * shape functions values container if they are
     * existing. Otherwise it gives you error which this value is
     * not exist.
     *
     * @param rResult the given Container that will be overwritten by the solution
     * @param rPoint the given local coordinates the gradients will be evaluated for
     * @return a matrix of gradients for each shape function
     */
    Matrix& ShapeFunctionsLocalGradients( Matrix& rResult, const LocalCoordinatesArrayType& rPoint ) const override
    {
        return ShapeFunctionsLocalGradientsImpl( rResult, rPoint );
    }

    /**
     * Input and output
     */

    /**
     * Turn back information as a string.
     *
     * @return String contains information about this geometry.
     * @see PrintData()
     * @see PrintInfo()
     */
    std::string Info() const override
    {
        return "3 dimensional tetrahedra with ten nodes in 3D space";
    }

    /**
     * Print information about this object.
     *
     * @param rOStream Stream to print into it.
     * @see PrintData()
     * @see Info()
     */
    void PrintInfo( std::ostream& rOStream ) const override
    {
        rOStream << "3 dimensional tetrahedra with ten nodes in 3D space";
    }

    /**
     * Print geometry's data into given stream.
     * Prints it's points by the order they stored in the geometry
     * and then center point of geometry.
     *
     * @param rOStream Stream to print into it.
     * @see PrintInfo()
     * @see Info()
     */
    void PrintData( std::ostream& rOStream ) const override
    {
        BaseType::PrintData( rOStream );
        std::cout << std::endl;
        MatrixType jacobian;
        this->Jacobian( jacobian, LocalCoordinatesArrayType() );
        rOStream << "    Jacobian in the origin\t : " << jacobian;
    }

protected:

    /**
     * there are no protected class members
     */

private:

    /**
     * Static Member Variables
     */
    static const GeometryData msGeometryData;

    ///@name Serialization
    ///@{

    friend class Serializer;

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, PointsArrayType );
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, PointsArrayType );
    }

    ///@}

    Tetrahedra3D10(): BaseType( PointsArrayType(), &msGeometryData ) {}

    /**
     * Private Operations
     */

    /**
     * TODO: TO BE VERIFIED
     */
    /**
     * Calculates the gradients in terms of local coordinateds
     * of all shape functions in a given point.
     *
     * @param rPoint the current point at which the gradients are calculated
     * @return the gradients of all shape functions
     * \f$ \frac{\partial N^i}{\partial \xi_j} \f$
     */
    static Matrix& ShapeFunctionsLocalGradientsImpl( Matrix& result,
            const LocalCoordinatesArrayType& rPoint )
    {
        // Area coordinates, pg. 205, Vol. I, Carey, Oden, Becker FEM
        const ValueType zeta1 = rPoint[0];
        const ValueType zeta2 = rPoint[1];
        const ValueType zeta3 = rPoint[2];
        const ValueType zeta0 = 1. - zeta1 - zeta2 - zeta3;

        const ValueType dzeta0dxi = -1.;
        const ValueType dzeta1dxi =  1.;
        const ValueType dzeta2dxi =  0.;
        const ValueType dzeta3dxi =  0.;

        const ValueType dzeta0deta = -1.;
        const ValueType dzeta1deta =  0.;
        const ValueType dzeta2deta =  1.;
        const ValueType dzeta3deta =  0.;

        const ValueType dzeta0dzeta = -1.;
        const ValueType dzeta1dzeta =  0.;
        const ValueType dzeta2dzeta =  0.;
        const ValueType dzeta3dzeta =  1.;

        if ( result.size1() != 10 || result.size2() != 3 )
            result.resize( 10, 3, false );

        result( 0, 0 ) = (4.*zeta0 - 1.)*dzeta0dxi;

        result( 1, 0 ) = (4.*zeta1 - 1.)*dzeta1dxi;

        result( 2, 0 ) = (4.*zeta2 - 1.)*dzeta2dxi;

        result( 3, 0 ) = (4.*zeta3 - 1.)*dzeta3dxi;

        result( 4, 0 ) = 4.*(zeta0*dzeta1dxi + dzeta0dxi*zeta1);

        result( 5, 0 ) = 4.*(zeta1*dzeta2dxi + dzeta1dxi*zeta2);

        result( 6, 0 ) = 4.*(zeta0*dzeta2dxi + dzeta0dxi*zeta2);

        result( 7, 0 ) = 4.*(zeta0*dzeta3dxi + dzeta0dxi*zeta3);

        result( 8, 0 ) = 4.*(zeta1*dzeta3dxi + dzeta1dxi*zeta3);

        result( 9, 0 ) = 4.*(zeta2*dzeta3dxi + dzeta2dxi*zeta3);

        result( 0, 1 ) = (4.*zeta0 - 1.)*dzeta0deta;

        result( 1, 1 ) = (4.*zeta1 - 1.)*dzeta1deta;

        result( 2, 1 ) = (4.*zeta2 - 1.)*dzeta2deta;

        result( 3, 1 ) = (4.*zeta3 - 1.)*dzeta3deta;

        result( 4, 1 ) = 4.*(zeta0*dzeta1deta + dzeta0deta*zeta1);

        result( 5, 1 ) = 4.*(zeta1*dzeta2deta + dzeta1deta*zeta2);

        result( 6, 1 ) = 4.*(zeta0*dzeta2deta + dzeta0deta*zeta2);

        result( 7, 1 ) = 4.*(zeta0*dzeta3deta + dzeta0deta*zeta3);

        result( 8, 1 ) = 4.*(zeta1*dzeta3deta + dzeta1deta*zeta3);

        result( 9, 1 ) = 4.*(zeta2*dzeta3deta + dzeta2deta*zeta3);

        result( 0, 2 ) = (4.*zeta0 - 1.)*dzeta0dzeta;

        result( 1, 2 ) = (4.*zeta1 - 1.)*dzeta1dzeta;

        result( 2, 2 ) = (4.*zeta2 - 1.)*dzeta2dzeta;

        result( 3, 2 ) = (4.*zeta3 - 1.)*dzeta3dzeta;

        result( 4, 2 ) = 4.*(zeta0*dzeta1dzeta + dzeta0dzeta*zeta1);

        result( 5, 2 ) = 4.*(zeta1*dzeta2dzeta + dzeta1dzeta*zeta2);

        result( 6, 2 ) = 4.*(zeta0*dzeta2dzeta + dzeta0dzeta*zeta2);

        result( 7, 2 ) = 4.*(zeta0*dzeta3dzeta + dzeta0dzeta*zeta3);

        result( 8, 2 ) = 4.*(zeta1*dzeta3dzeta + dzeta1dzeta*zeta3);

        result( 9, 2 ) = 4.*(zeta2*dzeta3dzeta + dzeta2dzeta*zeta3);

        return result;
    }

    /**
     * TODO: TO BE VERIFIED
     */
    /**
     * Calculates the values of all shape function in all integration points.
     * Integration points are expected to be given in local coordinates
     *
     * @param ThisMethod the current integration method
     * @return the matrix of values of every shape function in each integration point
     *
     */
    static Matrix CalculateShapeFunctionsIntegrationPointsValues(
        IntegrationMethod ThisMethod )
    {
        IntegrationPointsContainerType all_integration_points =
            AllIntegrationPoints();
        IntegrationPointsArrayType integration_points =
            all_integration_points[static_cast<int>(ThisMethod)];
        //number of integration points
        const int integration_points_number = integration_points.size();
        //number of nodes in current geometry
        const int points_number = 10;
        //setting up return matrix
        Matrix shape_function_values( integration_points_number, points_number );
        //loop over all integration points

        for ( int pnt = 0; pnt < integration_points_number; pnt++ )
        {
            const ValueType zeta1 = integration_points[pnt].X();
            const ValueType zeta2 = integration_points[pnt].Y();
            const ValueType zeta3 = integration_points[pnt].Z();
            const ValueType zeta0 = 1. - zeta1 - zeta2 - zeta3;

            shape_function_values( pnt, 0 ) = zeta0*(2.*zeta0 - 1.);

            shape_function_values( pnt, 1 ) = zeta1*(2.*zeta1 - 1.);

            shape_function_values( pnt, 2 ) = zeta2*(2.*zeta2 - 1.);

            shape_function_values( pnt, 3 ) = zeta3*(2.*zeta3 - 1.);

            shape_function_values( pnt, 4 ) = 4.*zeta0*zeta1;

            shape_function_values( pnt, 5 ) = 4.*zeta1*zeta2;

            shape_function_values( pnt, 6 ) = 4.*zeta2*zeta0;

            shape_function_values( pnt, 7 ) = 4.*zeta0*zeta3;

            shape_function_values( pnt, 8 ) = 4.*zeta1*zeta3;

            shape_function_values( pnt, 9 ) = 4.*zeta2*zeta3;
        }

        return shape_function_values;
    }

    /**
     * TODO: TO BE VERIFIED
     */
    /**
     * Calculates the local gradients of all shape functions in all integration points.
     * Integration points are expected to be given in local coordinates
     *
     * @param ThisMethod the current integration method
     * @return the vector of the gradients of all shape functions
     * in each integration point
     *
     */
    static ShapeFunctionsGradientsType CalculateShapeFunctionsIntegrationPointsLocalGradients(
        IntegrationMethod ThisMethod )
    {
        IntegrationPointsContainerType all_integration_points =
            AllIntegrationPoints();
        IntegrationPointsArrayType integration_points =
            all_integration_points[static_cast<int>(ThisMethod)];
        //number of integration points
        const int integration_points_number = integration_points.size();
        ShapeFunctionsGradientsType d_shape_f_values( integration_points_number );
        //initialising container
        //loop over all integration points

        for ( int pnt = 0; pnt < integration_points_number; pnt++ )
        {
            d_shape_f_values[pnt] = ShapeFunctionsLocalGradientsImpl(d_shape_f_values[pnt], integration_points[pnt]);
        }

        return d_shape_f_values;
    }

    static const IntegrationPointsContainerType AllIntegrationPoints()
    {
        IntegrationPointsContainerType integration_points =
        {
            {
                Quadrature < TetrahedronGaussLegendreIntegrationPoints1,
                3, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature < TetrahedronGaussLegendreIntegrationPoints2,
                3, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature < TetrahedronGaussLegendreIntegrationPoints3,
                3, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature < TetrahedronGaussLegendreIntegrationPoints4,
                3, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature < TetrahedronGaussLegendreIntegrationPoints5,
                3, IntegrationPoint<3> >::GenerateIntegrationPoints(),
            }
        };
        return integration_points;
    }

    static const ShapeFunctionsValuesContainerType AllShapeFunctionsValues()
    {
        ShapeFunctionsValuesContainerType shape_functions_values =
        {
            {
                Tetrahedra3D10<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::IntegrationMethod::GI_GAUSS_1 ),
                Tetrahedra3D10<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::IntegrationMethod::GI_GAUSS_2 ),
                Tetrahedra3D10<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::IntegrationMethod::GI_GAUSS_3 ),
                Tetrahedra3D10<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::IntegrationMethod::GI_GAUSS_4 ),
                Tetrahedra3D10<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::IntegrationMethod::GI_GAUSS_5 )
            }
        };
        return shape_functions_values;
    }

    static const ShapeFunctionsLocalGradientsContainerType AllShapeFunctionsLocalGradients()
    {
        ShapeFunctionsLocalGradientsContainerType shape_functions_local_gradients =
        {
            {
                Tetrahedra3D10<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(
                    GeometryData::IntegrationMethod::GI_GAUSS_1 ),
                Tetrahedra3D10<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(
                    GeometryData::IntegrationMethod::GI_GAUSS_2 ),
                Tetrahedra3D10<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(
                    GeometryData::IntegrationMethod::GI_GAUSS_3 ),
                Tetrahedra3D10<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(
                    GeometryData::IntegrationMethod::GI_GAUSS_4 ),
                Tetrahedra3D10<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(
                    GeometryData::IntegrationMethod::GI_GAUSS_5 )
            }
        };
        return shape_functions_local_gradients;
    }

    /**
     * Private Friends
     */

    template<class TOtherPointType> friend class Tetrahedra3D10;

    /**
     * Un accessible methods
     */

}; // Class Tetrahedra3D10

template<class TPointType> const
GeometryData Tetrahedra3D10<TPointType>::msGeometryData(
    3, 3, 3, GeometryData::IntegrationMethod::GI_GAUSS_2,
    Tetrahedra3D10<TPointType>::AllIntegrationPoints(),
    Tetrahedra3D10<TPointType>::AllShapeFunctionsValues(),
    AllShapeFunctionsLocalGradients()
);

} // namespace Kratos.

#endif // KRATOS_TETRAHEDRA_3D_4_H_INCLUDED  defined
