//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//












// This define must be HERE
#define DKRATOS_EXPORT_INTERFACE_2 1

// System includes
#include <string>
#include <iostream>
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/cfd_variables.h"
#include "includes/kernel.h"

namespace Kratos
{
// Useful variables
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( ADVPROJ )
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( CONV_PROJ )
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( PRESS_PROJ )
KRATOS_CREATE_DOUBLE_VARIABLE( DIVPROJ )
KRATOS_CREATE_DOUBLE_VARIABLE( PRESSURE_OLD_IT )
KRATOS_CREATE_DOUBLE_VARIABLE( C_SMAGORINSKY )
KRATOS_CREATE_DOUBLE_VARIABLE( MOLECULAR_VISCOSITY )
KRATOS_CREATE_DOUBLE_VARIABLE( TURBULENT_VISCOSITY )
KRATOS_CREATE_DOUBLE_VARIABLE( Y_WALL)
KRATOS_CREATE_VARIABLE( int, OSS_SWITCH )

// Legacy variables
KRATOS_CREATE_DOUBLE_VARIABLE( DYNAMIC_TAU )
KRATOS_CREATE_DOUBLE_VARIABLE( DYNAMIC_VISCOSITY)
KRATOS_CREATE_DOUBLE_VARIABLE( EFFECTIVE_VISCOSITY )
KRATOS_CREATE_DOUBLE_VARIABLE( KINEMATIC_VISCOSITY)
KRATOS_CREATE_DOUBLE_VARIABLE( THAWONE )
KRATOS_CREATE_DOUBLE_VARIABLE( THAWTWO )
KRATOS_CREATE_DOUBLE_VARIABLE( M )

KRATOS_CREATE_DOUBLE_VARIABLE( CROSS_WIND_STABILIZATION_FACTOR )

void KratosApplication::RegisterCFDVariables()
{
    // Useful variables
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( ADVPROJ )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( CONV_PROJ )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( PRESS_PROJ )
    KRATOS_REGISTER_VARIABLE( DIVPROJ )
    KRATOS_REGISTER_VARIABLE( PRESSURE_OLD_IT )
    KRATOS_REGISTER_VARIABLE( C_SMAGORINSKY )
    KRATOS_REGISTER_VARIABLE( MOLECULAR_VISCOSITY )
    KRATOS_REGISTER_VARIABLE( TURBULENT_VISCOSITY )
    KRATOS_REGISTER_VARIABLE( Y_WALL)
    KRATOS_REGISTER_VARIABLE( OSS_SWITCH )

    // Legacy variables
    KRATOS_REGISTER_VARIABLE( DYNAMIC_TAU )
    KRATOS_REGISTER_VARIABLE( DYNAMIC_VISCOSITY)
    KRATOS_REGISTER_VARIABLE( EFFECTIVE_VISCOSITY )
    KRATOS_REGISTER_VARIABLE( KINEMATIC_VISCOSITY)
    KRATOS_REGISTER_VARIABLE( THAWONE )
    KRATOS_REGISTER_VARIABLE( THAWTWO )
    KRATOS_REGISTER_VARIABLE( M )

    KRATOS_REGISTER_VARIABLE( CROSS_WIND_STABILIZATION_FACTOR )
}

}  // namespace Kratos.

// This define must be HERE
#undef DKRATOS_EXPORT_INTERFACE_2
