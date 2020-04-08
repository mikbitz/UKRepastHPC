#ifndef CONVERTOR
#define	CONVERTOR

#include <string>
#include <sstream>
#include <iostream>
//#include <stdlib.h>
#include <iomanip>
#include <vector>

class Convertor {
public:
    ~Convertor( );
    static Convertor* Get( );

    template< class T >
    const std::string ToString( const T& input ) const {
        std::stringstream stringStream;
        stringStream << input;

        return stringStream.str( );
    }
    
    double StringToNumber( const std::string& ) const;

    const std::vector<std::string> StringToWords( const std::string&, const char ) const;
    const std::string DoubleToPrecisionString( const double&, const unsigned& ) const;

    std::string ToLowercase( const std::string ) const;
    std::string RemoveWhiteSpace( const std::string ) const;

private:
    Convertor( );

    static Convertor* mThis;
};

#endif

