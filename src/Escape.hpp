#ifndef _GAMS_ESCAPE_HPP_
#define _GAMS_ESCAPE_HPP_

#include <string>

namespace gams {

/**
 * \brief Escape string so that it becomes suitable as an identifier in gams.
 * 
 * \param str the string to escape
 * \return the escaped string
 */
std::string escape_string( std::string str );

}

#endif


