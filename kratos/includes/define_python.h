//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

#if !defined(KRATOS_DEFINE_PYTHON_H_INCLUDED)
#define  KRATOS_DEFINE_PYTHON_H_INCLUDED

// System includes

// External includes
#include <boost/python.hpp>

// Project includes

// This function is used to print the ofstream-operator
// i.e. printing an object will give the same result in Python as in C++
// To be defined as the "__str__" function
// e.g. ".def("__str__", PrintObject<ProcessInfo>)"
// It replicates the function "self_ns::str(self))" of boost-python
template< class T>
std::string PrintObject(const T& rObject)
{
    std::stringstream ss;
    ss << rObject;
    return ss.str();
}

#ifdef KRATOS_REGISTER_IN_PYTHON_FLAG_IMPLEMENTATION
#undef KRATOS_REGISTER_IN_PYTHON_FLAG_IMPLEMENTATION
#endif
#define KRATOS_REGISTER_IN_PYTHON_FLAG_IMPLEMENTATION(flag) \
 scope().attr(#flag) = boost::ref(flag)      \

#ifdef KRATOS_REGISTER_IN_PYTHON_FLAG
#undef KRATOS_REGISTER_IN_PYTHON_FLAG
#endif
#define KRATOS_REGISTER_IN_PYTHON_FLAG(flag) \
    KRATOS_REGISTER_IN_PYTHON_FLAG_IMPLEMENTATION(flag);   \
    KRATOS_REGISTER_IN_PYTHON_FLAG_IMPLEMENTATION(NOT_##flag)

#ifdef KRATOS_REGISTER_IN_PYTHON_VARIABLE
#undef KRATOS_REGISTER_IN_PYTHON_VARIABLE
#endif
#define KRATOS_REGISTER_IN_PYTHON_VARIABLE(variable) \
    scope().attr(#variable) = boost::ref(variable);


#ifdef KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS
#undef KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS
#endif
#define KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(name) \
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(name) \
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(name##_X) \
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(name##_Y) \
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(name##_Z)

#endif // KRATOS_DEFINE_PYTHON_H_INCLUDED
