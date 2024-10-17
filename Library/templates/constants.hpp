///////////////////////////////////////////////////////////////////////////////
/// \file constants.hpp
/// \brief This file provides some preprocessor constants.
///
///  This file is part of the HILBERT program package for the numerical
///  solution of the Laplace equation with mixed boundary conditions by use of
///  BEM in 2D.
///
///  C++ adaptation of HILBERT V3.1 TUWien 2009-2013 for ANCSE17
///////////////////////////////////////////////////////////////////////////////
#ifndef _CONSTANTS_HPP_GUARD_mysolution
#define _CONSTANTS_HPP_GUARD_mysolution

#include <cmath>

#ifndef EPS
#define EPS 1e-12
#endif

#ifndef DEFAULT_ETA
#define DEFAULT_ETA 0.5
#endif

#ifndef GAUSS_ORDER
#define GAUSS_ORDER 16
#endif

#ifndef INFINITY
#define INFINITY (HUGE_VAL + HUGE_VAL)
#endif

#ifndef NAN
#define NAN (INFINITY - INFINITY)
#endif

#endif
