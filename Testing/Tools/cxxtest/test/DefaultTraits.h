// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#include <cxxtest/TestSuite.h>

//
// This test suite demonstrates the default ValueTraits
//

class DefaultTraits : public CxxTest::TestSuite
{
public:
    struct EightBytes
    {
        EightBytes() {}
        unsigned char data[8];
    };
    
    void testSmallDefaultTraits()
    {
        EightBytes x;
        for ( unsigned i = 0; i < sizeof(x.data); ++ i )
            x.data[i] = (unsigned char)i;
        TS_FAIL( x );
    }

    struct NineBytes
    {
        NineBytes() {}
        unsigned char data[9];
    };
    
    void testBigDefaultTraits()
    {
        NineBytes x;
        for ( unsigned i = 0; i < sizeof(x.data); ++ i )
            x.data[i] = (unsigned char)(0x98 + i);
        TS_FAIL( x );
    }
};

//
// Local Variables:
// compile-command: "perl test.pl"
// End:
//
