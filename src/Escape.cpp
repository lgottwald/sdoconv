#include <algorithm>
#include "Escape.hpp"

namespace gams
{

struct isVocale
{
   bool operator()( char c ) const
   {
      return c == 'a' || c == 'e' || c == 'i' || c == 'o' || c == 'u';
   }
};

struct isIllegal
{
   bool operator()( char c ) const
   {
      return c == ' ' || c == '-';
   }
};

void replace_string_inplace( std::string& subject, const std::string& search,
                             const std::string& replace )
{
   std::size_t pos = 0;

   while( ( pos = subject.find( search, pos ) ) != std::string::npos )
   {
      subject.replace( pos, search.length(), replace );
      pos += replace.length();
   }
}

std::string escape_string( std::string str )
{
   str.erase( std::remove_if( str.begin(), str.end(), isIllegal() ), str.end() );

   //replace forbidden characters
   replace_string_inplace( str, "Ä", "AE" );
   replace_string_inplace( str, "Ö", "OE" );
   replace_string_inplace( str, "Ü", "UE" );
   replace_string_inplace( str, "ä", "ae" );
   replace_string_inplace( str, "ö", "oe" );
   replace_string_inplace( str, "ü", "ue" );
   replace_string_inplace( str, "ß", "ss" );

   //if too long remove vocales
   if( str.size() > 59 )
      str.erase( std::remove_if( str.begin(), str.end(), isVocale() ), str.end() );

   //if still too long cut off end
   if( str.size() > 59 )
      str.erase( str.begin() + 59, str.end() );

   return str;
}
}
