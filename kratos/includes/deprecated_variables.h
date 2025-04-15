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











#if !defined(KRATOS_DEPRECATED_VARIABLES_H_INCLUDED )
#define  KRATOS_DEPRECATED_VARIABLES_H_INCLUDED



// System includes

// External includes


// Project includes
#include "includes/define.h"
#include "containers/variable_component.h"
#include "containers/vector_component_adaptor.h"
#include "includes/geometrical_object.h"
#include "containers/weak_pointer_vector.h"


#undef  KRATOS_EXPORT_MACRO
#define KRATOS_EXPORT_MACRO KRATOS_API

namespace Kratos
{

    //Define Variables by type:

    //bools

    //for Structural application:
    KRATOS_DEFINE_VARIABLE( bool, IS_INACTIVE )

    //for Level Set application:
    KRATOS_DEFINE_VARIABLE( bool, IS_DUPLICATED )
    KRATOS_DEFINE_VARIABLE( bool, SPLIT_ELEMENT )
    KRATOS_DEFINE_VARIABLE( bool, SPLIT_NODAL )


    //for PFEM fluids application:
    KRATOS_DEFINE_VARIABLE( int, IS_JACK_LINK )
    KRATOS_DEFINE_VARIABLE( int, IMPOSED_PRESSURE )
    KRATOS_DEFINE_VARIABLE( int, IMPOSED_VELOCITY_X )
    KRATOS_DEFINE_VARIABLE( int, IMPOSED_VELOCITY_Y )
    KRATOS_DEFINE_VARIABLE( int, IMPOSED_VELOCITY_Z )
    KRATOS_DEFINE_VARIABLE( int, IMPOSED_ANGULAR_VELOCITY_X )
    KRATOS_DEFINE_VARIABLE( int, IMPOSED_ANGULAR_VELOCITY_Y )
    KRATOS_DEFINE_VARIABLE( int, IMPOSED_ANGULAR_VELOCITY_Z )

    // For the DEM Application:
    KRATOS_DEFINE_DOUBLE_VARIABLE( IMPOSED_VELOCITY_X_VALUE )
    KRATOS_DEFINE_DOUBLE_VARIABLE( IMPOSED_VELOCITY_Y_VALUE )
    KRATOS_DEFINE_DOUBLE_VARIABLE( IMPOSED_VELOCITY_Z_VALUE )
    KRATOS_DEFINE_DOUBLE_VARIABLE( IMPOSED_ANGULAR_VELOCITY_X_VALUE )
    KRATOS_DEFINE_DOUBLE_VARIABLE( IMPOSED_ANGULAR_VELOCITY_Y_VALUE )
    KRATOS_DEFINE_DOUBLE_VARIABLE( IMPOSED_ANGULAR_VELOCITY_Z_VALUE )

    KRATOS_DEFINE_DOUBLE_VARIABLE( IS_INLET )
    KRATOS_DEFINE_DOUBLE_VARIABLE( IS_INTERFACE )
    KRATOS_DEFINE_DOUBLE_VARIABLE( IS_VISITED )
    KRATOS_DEFINE_DOUBLE_VARIABLE( IS_EROSIONABLE )

    KRATOS_DEFINE_DOUBLE_VARIABLE( IS_STRUCTURE )
    KRATOS_DEFINE_DOUBLE_VARIABLE( IS_POROUS )
    KRATOS_DEFINE_DOUBLE_VARIABLE( IS_WATER )
    KRATOS_DEFINE_DOUBLE_VARIABLE( IS_FLUID )
    KRATOS_DEFINE_DOUBLE_VARIABLE( IS_BOUNDARY )
    KRATOS_DEFINE_DOUBLE_VARIABLE( IS_FREE_SURFACE )
    KRATOS_DEFINE_DOUBLE_VARIABLE( IS_AIR_EXIT )
    KRATOS_DEFINE_DOUBLE_VARIABLE( IS_LAGRANGIAN_INLET )
    KRATOS_DEFINE_DOUBLE_VARIABLE( IS_WATER_ELEMENT )


    KRATOS_DEFINE_DOUBLE_VARIABLE( IS_BURN )
    KRATOS_DEFINE_DOUBLE_VARIABLE( IS_DRIPPING )
    KRATOS_DEFINE_DOUBLE_VARIABLE( IS_PERMANENT )
    KRATOS_DEFINE_DOUBLE_VARIABLE( IS_WALL )

    // KRATOS_DEFINE_DOUBLE_VARIABLE( Ypr ) //var name does not follow standard
    // KRATOS_DEFINE_DOUBLE_VARIABLE( Yox )
    // KRATOS_DEFINE_DOUBLE_VARIABLE( Yfuel )
    // KRATOS_DEFINE_DOUBLE_VARIABLE( Hfuel )
    // KRATOS_DEFINE_DOUBLE_VARIABLE( Hpr )
    // KRATOS_DEFINE_DOUBLE_VARIABLE( Hpr1 )
    // KRATOS_DEFINE_DOUBLE_VARIABLE( Hox )

    KRATOS_DEFINE_DOUBLE_VARIABLE( RADIATIVE_INTENSITY )
    KRATOS_DEFINE_DOUBLE_VARIABLE( RADIATIVE_INTENSITY_1 )
    KRATOS_DEFINE_DOUBLE_VARIABLE( RADIATIVE_INTENSITY_2 )
    KRATOS_DEFINE_DOUBLE_VARIABLE( RADIATIVE_INTENSITY_3 )
    KRATOS_DEFINE_DOUBLE_VARIABLE( RADIATIVE_INTENSITY_4 )
    KRATOS_DEFINE_DOUBLE_VARIABLE( RADIATIVE_INTENSITY_5 )
    KRATOS_DEFINE_DOUBLE_VARIABLE( RADIATIVE_INTENSITY_6 )
    KRATOS_DEFINE_DOUBLE_VARIABLE( RADIATIVE_INTENSITY_7 )
    KRATOS_DEFINE_DOUBLE_VARIABLE( RADIATIVE_INTENSITY_8 )
    KRATOS_DEFINE_DOUBLE_VARIABLE( RADIATIVE_INTENSITY_9 )
    KRATOS_DEFINE_DOUBLE_VARIABLE( RADIATIVE_INTENSITY_10 )
    KRATOS_DEFINE_DOUBLE_VARIABLE( RADIATIVE_INTENSITY_11 )
    KRATOS_DEFINE_DOUBLE_VARIABLE( RADIATIVE_INTENSITY_12 )
    KRATOS_DEFINE_DOUBLE_VARIABLE( RADIATIVE_INTENSITY_13 )
    KRATOS_DEFINE_DOUBLE_VARIABLE( RADIATIVE_INTENSITY_14 )
    KRATOS_DEFINE_DOUBLE_VARIABLE( RADIATIVE_INTENSITY_15 )
    KRATOS_DEFINE_DOUBLE_VARIABLE( RADIATIVE_INTENSITY_16 )
    KRATOS_DEFINE_DOUBLE_VARIABLE( RADIATIVE_INTENSITY_17 )
    KRATOS_DEFINE_DOUBLE_VARIABLE( RADIATIVE_INTENSITY_18 )
    KRATOS_DEFINE_DOUBLE_VARIABLE( RADIATIVE_INTENSITY_19 )
    KRATOS_DEFINE_DOUBLE_VARIABLE( RADIATIVE_INTENSITY_20 )
    KRATOS_DEFINE_DOUBLE_VARIABLE( RADIATIVE_INTENSITY_21 )
    KRATOS_DEFINE_DOUBLE_VARIABLE( RADIATIVE_INTENSITY_22 )
    KRATOS_DEFINE_DOUBLE_VARIABLE( RADIATIVE_INTENSITY_23 )
    KRATOS_DEFINE_DOUBLE_VARIABLE( RADIATIVE_INTENSITY_24 )

    // KRATOS_DEFINE_DOUBLE_VARIABLE( rhoD )
    // KRATOS_DEFINE_DOUBLE_VARIABLE( xi )
    // KRATOS_DEFINE_DOUBLE_VARIABLE( a )
    // KRATOS_DEFINE_DOUBLE_VARIABLE( b )


    KRATOS_DEFINE_DOUBLE_VARIABLE( IS_SLIP )

    //for Level Set application:
    KRATOS_DEFINE_DOUBLE_VARIABLE( IS_DIVIDED )

    // KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( xi_c )

    KRATOS_DEFINE_VARIABLE( WeakPointerVector< GeometricalObject<Node<3> > >, NEIGHBOUR_EMBEDDED_FACES )

}  // namespace Kratos.

#undef  KRATOS_EXPORT_MACRO
#define KRATOS_EXPORT_MACRO KRATOS_NO_EXPORT

#endif // KRATOS_DEPRECATED_VARIABLES_H_INCLUDED  defined
