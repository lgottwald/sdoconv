#ifndef _GAMS_SET_INDEX_HPP_
#define _GAMS_SET_INDEX_HPP_

#include <string>
#include <limits>

namespace gams {

/**
 * Class representing a set index. the set is given by name and an optional
 * offset.
 */
class SetIndex 
{
public:
   /**
    * Contructor for implicit construction from string literal with offset 0.
    */
   SetIndex(const char* name) : name(name), offset(0) {}
   /**
    * Constructor for construction from string with offset 0
    */
   SetIndex(std::string name) : name(std::move(name)), offset(0) {}
   
   /**
    * Constructor for set index with an offset
    */
   SetIndex(std::string name, int offset) : name(std::move(name)), offset(offset) {}

   /**
    * Named "constructor" to get the SetIndex object that represents the first index
    * for the set with the given name.
    * 
    * \param name name of set
    * \return SetIndex object that represents first of set with given name
    */
   static SetIndex First(std::string name) {
      return SetIndex(std::move(name), FIRST);
   }

   /**
    * \brief Check if set index is the first of the set
    * \return true if the index represents the first element.
    */
   bool isFirst() const {
      return offset == FIRST;
   }

   /**
    * Get the name of the set
    * \return const reference to the name
    */
   const std::string& getName() const {
      return name;
   }

   /**
    * Get the offset of the set
    * \return the offset
    */
   int getOffset() const {
      return offset;
   }

private:
   /**
    * Constant that indicates the first value of a set if used as offset
    */
   constexpr static int FIRST = std::numeric_limits<int>::min();
   /**
    * Name of the set that is indexed
    */
   std::string name;

   /**
    * Offset for the index
    */
   int offset;
};

}

#endif