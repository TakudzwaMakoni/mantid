// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once
/**
 * @file test_bar.t.h
 * Test one for the joint test ehm, test.
 *
 * @author Gašper Ažman (GA), gasper.azman@gmail.com
 * @version 1.0
 * @since 2008-08-29 10:04:06 AM
 */

#include <cxxtest/TestSuite.h>
#include "requirement.hpp"

class TestBar : public CxxTest::TestSuite
{
public:
    void test_foo() {
        TS_ASSERT(call_a_requirement());
    }
};
