//
//   Project Name:        Kratos
//   Last Modified by:    $Author:   JMCarbonell $
//   Date:                $Date:   December 2015 $
//   Revision:            $Revision:         1.7 $
//
//

#if !defined(KRATOS_PRISM_GAUSS_LOBATTO_INTEGRATION_POINTS_H_INCLUDED )
#define  KRATOS_PRISM_GAUSS_LOBATTO_INTEGRATION_POINTS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "integration/quadrature.h"

//TO BE COMPLETED: Only the needed ones have been implemented

namespace Kratos
{

class PrismGaussLobattoIntegrationPoints1
{
public:
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 3;

    typedef IntegrationPoint<3> IntegrationPointType;

    typedef boost::array<IntegrationPointType, 3> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 3;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType( 0.0 , 0.0 , 0.5 , 1.00 / 6.00 ),
            IntegrationPointType( 1.0 , 0.0 , 0.5 , 1.00 / 6.00 ),
            IntegrationPointType( 0.0 , 1.0 , 0.5 , 1.00 / 6.00 )
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Prism Gauss-Lobatto quadrature 1 ";
        return buffer.str();
    }
}; // Class PrismGaussLobattoIntegrationPoints1

class PrismGaussLobattoIntegrationPoints2
{
public:
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 3;

    typedef IntegrationPoint<3> IntegrationPointType;

    typedef boost::array<IntegrationPointType, 6> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 6;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType( 0.0 , 0.0 , 0.0 , 1.00 / 12.00 ),
            IntegrationPointType( 1.0 , 0.0 , 0.0 , 1.00 / 12.00 ),
            IntegrationPointType( 0.0 , 1.0 , 0.0 , 1.00 / 12.00 ),
            IntegrationPointType( 0.0 , 0.0 , 1.0 , 1.00 / 12.00 ),
            IntegrationPointType( 1.0 , 0.0 , 1.0 , 1.00 / 12.00 ),
            IntegrationPointType( 0.0 , 1.0 , 1.0 , 1.00 / 12.00 )
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Prism Gauss-Lobatto quadrature 2 ";
        return buffer.str();
    }
}; // Class PrismGaussLobattoIntegrationPoints2

}

#endif // KRATOS_PRISM_GAUSS_LOBATTO_INTEGRATION_POINTS_H_INCLUDED defined
