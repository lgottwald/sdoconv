#include <sdo/ButcherTableau.hpp>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <boost/lexical_cast.hpp>
#include <stack>
#include "SetIndex.hpp"
#include "GamsGenerator.hpp"
#include "Escape.hpp"

using namespace sdo;

using LkpTypePair = std::pair<LookupTable* const, gams::LookupData>;


static bool has_spline_type( const std::unordered_map<LookupTable*, gams::LookupData>& lkpData )
{
   auto cond = []( const LkpTypePair & val )
   {
      return val.second.type == gams::LookupFormulationType::SPLINE;
   };

   return std::any_of( lkpData.begin(), lkpData.end(), cond );
}


static bool has_sos2_type( const std::unordered_map<LookupTable*, gams::LookupData>& lkpData )
{
   return std::any_of(
             lkpData.begin(),
             lkpData.end(),
             []( const LkpTypePair & val )
               {
                  return val.second.type == gams::LookupFormulationType::SOS2;
               }
   );
}


namespace gams
{

void GamsGenerator::setLookupFormulationTypes( LookupFormulationType type )
{
   std::unordered_map<LookupTable*, LookupData> lkpData;

   for( auto & entry : exprGraph_.getSymbolTable() )
   {
      if( entry.second->op == ExpressionGraph::LOOKUP_TABLE )
      {
         lkpData[entry.second->lookup_table] = LookupData {entry.first, type};
      }
      else if( entry.second->op == ExpressionGraph::APPLY_LOOKUP )
      {
         lkpData[entry.second->child1->lookup_table] = LookupData {entry.first, type};
      }
   }

   lkpData_.swap( lkpData );
}

void GamsGenerator::setLookupFormulationType( LookupTable* lookup, LookupData data )
{
   lkpData_[lookup] = std::move( data );
}

void GamsGenerator::createStateSymbols()
{
   std::unordered_map<ExpressionGraph::Node*, Symbol> newSymbols;

   for( auto & entry : exprGraph_.getSymbolTable() )
   {
      std::stack<std::pair<int, ExpressionGraph::Node*>> stack;
      stack.emplace( 0, entry.second );
      Symbol symbol = entry.first;
      int level = 0;
      bool start = true;

      while( !stack.empty() )
      {
         std::pair<int, ExpressionGraph::Node*> top = stack.top();
         stack.pop();
         auto range = exprGraph_.getSymbol( top.second );

         if( !start && !range.empty() )
            continue;

         if( newSymbols.find( top.second ) != newSymbols.end() && top.first == 0 )
            continue;

         start = false;


         switch( top.second->op )
         {
         case ExpressionGraph::INTEG:
         {
            if( range.empty() )
            {
               if( top.first == 1 )
               {
                  std::string s( symbol.get() );
                  s += "_LV";
                  s += boost::lexical_cast<std::string>( ++level );
                  newSymbols[top.second] = Symbol( std::move( s ) );
                  continue;
               }
               else
               {
                  newSymbols.emplace( top.second, Symbol() );
                  stack.emplace( 1, top.second );
               }
            }

            stack.emplace( 0, top.second->child2 );
            stack.emplace( 0, top.second->child1 );
            continue;
         }

         case ExpressionGraph::IF:
         case ExpressionGraph::DELAY_FIXED:
         case ExpressionGraph::PULSE_TRAIN:
         case ExpressionGraph::RAMP:
            stack.emplace( 0, top.second->child3 );

         case ExpressionGraph::APPLY_LOOKUP:
         case ExpressionGraph::PULSE:
         case ExpressionGraph::ACTIVE_INITIAL:
         case ExpressionGraph::STEP:
         case ExpressionGraph::RANDOM_UNIFORM:
         case ExpressionGraph::PLUS:
         case ExpressionGraph::MINUS:
         case ExpressionGraph::MULT:
         case ExpressionGraph::DIV:
         case ExpressionGraph::G:
         case ExpressionGraph::GE:
         case ExpressionGraph::L:
         case ExpressionGraph::LE:
         case ExpressionGraph::EQ:
         case ExpressionGraph::NEQ:
         case ExpressionGraph::AND:
         case ExpressionGraph::OR:
         case ExpressionGraph::POWER:
         case ExpressionGraph::LOG:
         case ExpressionGraph::MIN:
         case ExpressionGraph::MAX:
         case ExpressionGraph::MODULO:
            stack.emplace( 0, top.second->child2 );

         case ExpressionGraph::INITIAL:
         case ExpressionGraph::UMINUS:
         case ExpressionGraph::SQRT:
         case ExpressionGraph::EXP:
         case ExpressionGraph::LN:
         case ExpressionGraph::ABS:
         case ExpressionGraph::INTEGER:
         case ExpressionGraph::NOT:
         case ExpressionGraph::SIN:
         case ExpressionGraph::COS:
         case ExpressionGraph::TAN:
         case ExpressionGraph::ARCSIN:
         case ExpressionGraph::ARCCOS:
         case ExpressionGraph::ARCTAN:
         case ExpressionGraph::SINH:
         case ExpressionGraph::COSH:
         case ExpressionGraph::TANH:
            stack.emplace( 0, top.second->child1 );

         case ExpressionGraph::TIME:
         case ExpressionGraph::CONSTANT:
         case ExpressionGraph::CONTROL:
         case ExpressionGraph::LOOKUP_TABLE:
         case ExpressionGraph::NIL:
            continue;
         };
      }
   }

   for( auto & entry : newSymbols )
   {
      exprGraph_.addSymbol( entry.second, entry.first );
   }
}

void GamsGenerator::controlSet( std::string set )
{
   auto& s = sets_[set];
   s.second = std::max( ++s.first, s.second );
}

void GamsGenerator::releaseSet( std::string set )
{
   auto& s = sets_[set];
   --s.first;
}

void GamsGenerator::createSet( std::string set )
{
   sets_.emplace( set, std::make_pair( 0, 0 ) );
}



std::string GamsGenerator::getSets( std::initializer_list< SetIndex > list ) const
{
   std::ostringstream stringstream;
   bool start = true;

   for( const SetIndex & idx : list )
   {
      if( !start )
      {
         stringstream << ", ";
      }

      start = false;

      if( idx.isFirst() )
      {
         stringstream << "'0'";
      }
      else
      {
         auto iter = sets_.find( idx.getName() );
         assert( iter != sets_.end() );
         stringstream << idx.getName();

         for( int i = 0; i < iter->second.first; ++i )
         {
            stringstream << idx.getName();
         }

         int offset = idx.getOffset();

         if( offset > 0 )
            stringstream << "+" << offset;
         else if( offset < 0 )
            stringstream << offset;
      }
   }

   return stringstream.str();
}

std::string GamsGenerator::getInitialSets() const
{
   if( tableau_.getName() == ButcherTableau::EULER )
   {
      return getSets( {SetIndex::First( "t" )} );
   }
   else
   {
      return getSets( {SetIndex::First( "t" ), SetIndex::First( "p" )} );
   }
}

std::string GamsGenerator::getVarSets() const
{
   if( tableau_.getName() == ButcherTableau::EULER )
   {
      return getSets( { "t" } );
   }
   else
   {
      return getSets( { "t", "p" } );
   }
}

void GamsGenerator::createDivisionGuards( std::vector<std::pair<int, std::string>>& varValues )
{
   std::stack<ExpressionGraph::Node*> stack;
   std::unordered_set<ExpressionGraph::Node*> nodes;

   for( auto & entry : exprGraph_.getSymbolTable() )
   {
      stack.push( entry.second );
   }

   int divisors = 0;
   std::ostringstream ss;

   while( !stack.empty() )
   {
      ExpressionGraph::Node* top = stack.top();
      stack.pop();

      if( nodes.find( top ) != nodes.end() || top->type != ExpressionGraph::DYNAMIC_NODE )
         continue;

      nodes.emplace( top );

      switch( top->op )
      {
      case ExpressionGraph::DIV:
      {
         if( top->child2->type == ExpressionGraph::DYNAMIC_NODE )
         {
            auto   range = exprGraph_.getSymbol( top->child2 );
            Symbol symb;

            if( range.empty() )
            {
               ss << "Divisor" << divisors++;
               symb = ss.str();
               ss.str( std::string() );
               exprGraph_.addSymbol( symb, top->child2 );
            }
            else
            {
               symb = range.begin()->second;
            }

            ss << escape_string(symb) << ".lo(" << getVarSets() << ") = EPSILON;\n";
            varValues.emplace_back( 0, ss.str() );
            ss.str( std::string() );
            stack.emplace( top->child2 );
         }

         stack.emplace( top->child1 );
         continue;
      }

      case ExpressionGraph::IF:
      case ExpressionGraph::DELAY_FIXED:
      case ExpressionGraph::PULSE_TRAIN:
      case ExpressionGraph::RAMP:
         stack.emplace( top->child3 );

      case ExpressionGraph::APPLY_LOOKUP:
      case ExpressionGraph::PULSE:
      case ExpressionGraph::ACTIVE_INITIAL:
      case ExpressionGraph::STEP:
      case ExpressionGraph::RANDOM_UNIFORM:
      case ExpressionGraph::PLUS:
      case ExpressionGraph::MINUS:
      case ExpressionGraph::MULT:
      case ExpressionGraph::G:
      case ExpressionGraph::GE:
      case ExpressionGraph::L:
      case ExpressionGraph::LE:
      case ExpressionGraph::EQ:
      case ExpressionGraph::NEQ:
      case ExpressionGraph::AND:
      case ExpressionGraph::OR:
      case ExpressionGraph::POWER:
      case ExpressionGraph::LOG:
      case ExpressionGraph::MIN:
      case ExpressionGraph::MAX:
      case ExpressionGraph::MODULO:
      case ExpressionGraph::INTEG:
         stack.emplace( top->child2 );

      case ExpressionGraph::INITIAL:
      case ExpressionGraph::UMINUS:
      case ExpressionGraph::SQRT:
      case ExpressionGraph::EXP:
      case ExpressionGraph::LN:
      case ExpressionGraph::ABS:
      case ExpressionGraph::INTEGER:
      case ExpressionGraph::NOT:
      case ExpressionGraph::SIN:
      case ExpressionGraph::COS:
      case ExpressionGraph::TAN:
      case ExpressionGraph::ARCSIN:
      case ExpressionGraph::ARCCOS:
      case ExpressionGraph::ARCTAN:
      case ExpressionGraph::SINH:
      case ExpressionGraph::COSH:
      case ExpressionGraph::TANH:
         stack.emplace( top->child1 );

      case ExpressionGraph::TIME:
      case ExpressionGraph::CONSTANT:
      case ExpressionGraph::CONTROL:
      case ExpressionGraph::LOOKUP_TABLE:
      case ExpressionGraph::NIL:
         continue;
      };
   }
}

void GamsGenerator::indexSos2Lookups()
{
   std::stack<ExpressionGraph::Node*> stack;
   std::unordered_set<ExpressionGraph::Node*> nodes;

   for( auto & entry : exprGraph_.getSymbolTable() )
   {
      stack.push( entry.second );
   }

   std::ostringstream ss;

   while( !stack.empty() )
   {
      ExpressionGraph::Node* top = stack.top();
      stack.pop();

      if( nodes.find( top ) != nodes.end() || top->type != ExpressionGraph::DYNAMIC_NODE )
         continue;

      nodes.emplace( top );

      switch( top->op )
      {
      case ExpressionGraph::APPLY_LOOKUP:
      {
         LookupData& lkpData = lkpData_[top->child1->lookup_table];
         if(lkpData.type == LookupFormulationType::SPLINE)
            continue;
         auto r = sos2LkpIds_.find( top );
         int id;

         if( r == sos2LkpIds_.end() )
         {
            id = ++lkpData.usages;
            sos2LkpIds_[top] = id;
         }
         else
         {
            id = r->second;
         }

         stack.emplace( top->child2 );
         continue;
      }

      case ExpressionGraph::IF:
      case ExpressionGraph::DELAY_FIXED:
      case ExpressionGraph::PULSE_TRAIN:
      case ExpressionGraph::RAMP:
         stack.emplace( top->child3 );

      case ExpressionGraph::PULSE:
      case ExpressionGraph::ACTIVE_INITIAL:
      case ExpressionGraph::STEP:
      case ExpressionGraph::RANDOM_UNIFORM:
      case ExpressionGraph::PLUS:
      case ExpressionGraph::MINUS:
      case ExpressionGraph::MULT:
      case ExpressionGraph::DIV:
      case ExpressionGraph::G:
      case ExpressionGraph::GE:
      case ExpressionGraph::L:
      case ExpressionGraph::LE:
      case ExpressionGraph::EQ:
      case ExpressionGraph::NEQ:
      case ExpressionGraph::AND:
      case ExpressionGraph::OR:
      case ExpressionGraph::POWER:
      case ExpressionGraph::LOG:
      case ExpressionGraph::MIN:
      case ExpressionGraph::MAX:
      case ExpressionGraph::MODULO:
      case ExpressionGraph::INTEG:
         stack.emplace( top->child2 );

      case ExpressionGraph::INITIAL:
      case ExpressionGraph::UMINUS:
      case ExpressionGraph::SQRT:
      case ExpressionGraph::EXP:
      case ExpressionGraph::LN:
      case ExpressionGraph::ABS:
      case ExpressionGraph::INTEGER:
      case ExpressionGraph::NOT:
      case ExpressionGraph::SIN:
      case ExpressionGraph::COS:
      case ExpressionGraph::TAN:
      case ExpressionGraph::ARCSIN:
      case ExpressionGraph::ARCCOS:
      case ExpressionGraph::ARCTAN:
      case ExpressionGraph::SINH:
      case ExpressionGraph::COSH:
      case ExpressionGraph::TANH:
         stack.emplace( top->child1 );

      case ExpressionGraph::TIME:
      case ExpressionGraph::CONSTANT:
      case ExpressionGraph::CONTROL:
      case ExpressionGraph::LOOKUP_TABLE:
      case ExpressionGraph::NIL:
         continue;
      };
   }
}

void GamsGenerator::initTableau( ButcherTableau::Name tableau )
{
   tableau_.setTableau( tableau );
   createSet( "t" );

   if( tableau != ButcherTableau::EULER )
      createSet( "p" );
}

void GamsGenerator::translateSymbol( std::ostream& stream, Symbol s, bool initial )
{
   auto node = exprGraph_.getNode( s );

   std::string varName = escape_string( s );

   switch( node->type )
   {
   case ExpressionGraph::CONSTANT_NODE:
      stream << varName;
      break;

   case ExpressionGraph::DYNAMIC_NODE:
      if( node->op == ExpressionGraph::CONTROL )
      {
         switch( node->control_size )
         {
         case 0:
            stream << varName;
            break;

         case 1:
            stream << varName;

            if( initial )
               stream << "(" << getSets( {SetIndex::First( "t" ) } ) << ")";
            else
               stream << "(" << getSets( {"t"} ) <<  ")";

            break;

         default:
            if( initial )
            {
               stream << "('0')";
            }
            else
            {
               std::string t = getSets( {"t"} );
               std::string csize = boost::lexical_cast<std::string>( node->control_size );
               stream << "sum(t" << csize << "$(ord(" << t << ") > (ord(t" << csize << ")-1)*" << csize
                      << " and ord(" << t << ") <= ord(t" << csize << ")*" << csize << "),"
                      << varName << "(t" << csize << "))";
            }
         }
      }
      else
      {
         if( initial )
         {
            stream << varName;

            if( node->init == ExpressionGraph::CONSTANT_INIT )
               stream << ".lo";

            stream << "(" << getInitialSets() << ")";
         }
         else
         {
            stream << varName << "(" << getVarSets() << ")";
         }
      }

      break;

   case ExpressionGraph::STATIC_NODE:
      stream << varName << "(" << getSets( {"t"} ) << ")";
      break;

   case ExpressionGraph::UNKNOWN:
      assert( false );
   }

}

void GamsGenerator::translate( std::ostream& stream, ExpressionGraph::Node* root, bool def, bool initial )
{
   std::stack<std::pair<int, ExpressionGraph::Node*>> stack;
   stack.emplace( 0, root );

   do
   {
      std::pair<int, ExpressionGraph::Node*>& top = stack.top();
      ExpressionGraph::Node* node = top.second;

      // if not translating a definition or not the root node of a definition
      // emit the symbol of a node if it exists
      if( !def || (def && node != root) )
      {
         auto range = exprGraph_.getSymbol( node );

         if( !range.empty() )
         {
            //symbol exists
            //for initial translation translate symbol only if it is a state, else use its initial value
            if ( !initial || ( initial && node->op == ExpressionGraph::INTEG ))
               translateSymbol(stream, range.begin()->second, initial );
            else
               stream << boost::lexical_cast<std::string>(node->value);
            stack.pop();
            continue;
         }
      }

      switch( node->op )
      {
      case ExpressionGraph::INTEG:
         stack.pop();
         if(initial)
            stack.emplace( 0, node->child2 );
         else
            stack.emplace( 0, node->child1 );
         continue;

      case ExpressionGraph::TIME:
         stream << "TIME(" << getSets( { "t" } ) << ")";
         stack.pop();
         continue;

      case ExpressionGraph::CONSTANT:
         stream << boost::lexical_cast<std::string>( node->value );
         stack.pop();
         continue;

      case ExpressionGraph::IF:
         switch( top.first )
         {
         case 0:
            stream << "(";
            ++top.first;
            stack.emplace( 0, node->child1 );
            continue;

         case 1:
            stream << ")*(";
            ++top.first;
            stack.emplace( 0, node->child2 );
            continue;

         case 2:
            stream << ")+(1-(";
            ++top.first;
            stack.emplace( 0, node->child1 );
            continue;

         case 3:
            stream << "))*(";
            ++top.first;
            stack.emplace( 0, node->child3 );
            continue;

         case 4:
            stream << ")";
            stack.pop();
            continue;
         }

      case ExpressionGraph::ACTIVE_INITIAL:
         stack.pop();

         if( initial )
            stack.emplace( 0, node->child2 );
         else
            stack.emplace( 0, node->child1 );

         continue;

      case ExpressionGraph::PULSE:
         switch( top.first )
         {
         case 0:
            stream << "( (TIME(" << getSets( { "t" } ) << ")+TIMESTEP/2) > ";
            ++top.first;
            stack.emplace( 0, node->child1 );
            continue;

         case 1:
            stream << " and (TIME(" << getSets( { "t" } ) << ")+TIMESTEP/2) < (";
            ++top.first;
            stack.emplace( 0, node->child1 );
            continue;

         case 2:
            stream << "+";
            ++top.first;
            stack.emplace( 0, node->child2 );
            continue;

         case 3:
            stream << ") )";
            stack.pop();
            continue;
         }

      case ExpressionGraph::PULSE_TRAIN:
         switch( top.first )
         {
         case 0:
            stream << "(mod(TIME(" << getSets( { "t" } ) << "), ";
            ++top.first;
            stack.emplace( 0, node->child2 );
            continue;

         case 1:
            stream << ")+TIMESTEP/2) > ";
            ++top.first;
            stack.emplace( 0, node->child1->child1 );
            continue;

         case 2:
            stream << " and (mod(TIME(" << getSets( { "t" } ) << "), ";
            ++top.first;
            stack.emplace( 0, node->child2 );
            continue;

         case 3:
            stream << ")+TIMESTEP/2) < (";
            ++top.first;
            stack.emplace( 0, node->child1->child1 );
            continue;

         case 4:
            stream << "+";
            ++top.first;
            stack.emplace( 0, node->child1->child2 );
            continue;

         case 5:
            stream << ") and ( TIME(" << getSets( { "t" } ) << ")+TIMESTEP/2 < ";
            ++top.first;
            stack.emplace( 0, node->child3 );
            continue;

         case 6:
            stream << ")";
            stack.pop();
         }

         continue;

      case ExpressionGraph::STEP:
         switch( top.first )
         {
         case 0:
            stream << "(TIME(" << getSets( {"t"} ) << ")+TIMESTEP/2 > ";
            ++top.first;
            stack.emplace( 0, node->child2 );
            continue;

         case 1:
            stream << ")*(";
            ++top.first;
            stack.emplace( 0, node->child1 );
            continue;

         case 2:
            stream << ")";
            stack.pop();
            continue;
         }

      case ExpressionGraph::RAMP:
         switch( top.first )
         {
         case 0:
            stream << "(";
            ++top.first;
            stack.emplace( 0, node->child1 );
            continue;

         case 1:
            stream << " * ( min(TIME(" << getSets( {"t"} ) << "),";
            ++top.first;
            stack.emplace( 0, node->child3 );
            continue;

         case 2:
            stream << ") - ";
            ++top.first;
            stack.emplace( 0, node->child2 );
            continue;

         case 3:
            stream << "))$(TIME(" << getSets( {"t"} ) << ") > ";
            ++top.first;
            stack.emplace( 0, node->child2 );
            continue;

         case 4:
            stream << ")";
            stack.pop();
         }

         continue;

      case ExpressionGraph::RANDOM_UNIFORM:
         stack.pop();
         stream << "TODO_RANDOM_UNIFORM";
         assert( false );
         continue;

      case ExpressionGraph::ABS:
         switch( top.first )
         {
         case 0:
            stream << "abs(";
            ++top.first;
            stack.emplace( 0, node->child1 );
            continue;

         case 1:
            stream << ")";
            stack.pop();
            continue;
         }

      case ExpressionGraph::SIN:
         switch( top.first )
         {
         case 0:
            stream << "sin(";
            ++top.first;
            stack.emplace( 0, node->child1 );
            continue;

         case 1:
            stream << ")";
            stack.pop();
            continue;
         }

      case ExpressionGraph::COS:
         switch( top.first )
         {
         case 0:
            stream << "cos(";
            ++top.first;
            stack.emplace( 0, node->child1 );
            continue;

         case 1:
            stream << ")";
            stack.pop();
            continue;
         }

      case ExpressionGraph::TAN:
         switch( top.first )
         {
         case 0:
            stream << "tan(";
            ++top.first;
            stack.emplace( 0, node->child1 );
            continue;

         case 1:
            stream << ")";
            stack.pop();
            continue;
         }

      case ExpressionGraph::ARCSIN:
         switch( top.first )
         {
         case 0:
            stream << "arcsin(";
            ++top.first;
            stack.emplace( 0, node->child1 );
            continue;

         case 1:
            stream << ")";
            stack.pop();
            continue;
         }

      case ExpressionGraph::ARCCOS:
         switch( top.first )
         {
         case 0:
            stream << "arccos(";
            ++top.first;
            stack.emplace( 0, node->child1 );
            continue;

         case 1:
            stream << ")";
            stack.pop();
            continue;
         }

      case ExpressionGraph::ARCTAN:
         switch( top.first )
         {
         case 0:
            stream << "arctan(";
            ++top.first;
            stack.emplace( 0, node->child1 );
            continue;

         case 1:
            stream << ")";
            stack.pop();
            continue;
         }

      case ExpressionGraph::SINH:
         switch( top.first )
         {
         case 0:
            stream << "sinh(";
            ++top.first;
            stack.emplace( 0, node->child1 );
            continue;

         case 1:
            stream << ")";
            stack.pop();
            continue;
         }

      case ExpressionGraph::COSH:
         switch( top.first )
         {
         case 0:
            stream << "cosh(";
            ++top.first;
            stack.emplace( 0, node->child1 );
            continue;

         case 1:
            stream << ")";
            stack.pop();
            continue;
         }

      case ExpressionGraph::TANH:
         switch( top.first )
         {
         case 0:
            stream << "tanh(";
            ++top.first;
            stack.emplace( 0, node->child1 );
            continue;

         case 1:
            stream << ")";
            stack.pop();
            continue;
         }

      case ExpressionGraph::EXP:
         switch( top.first )
         {
         case 0:
            stream << "exp(";
            ++top.first;
            stack.emplace( 0, node->child1 );
            continue;

         case 1:
            stream << ")";
            stack.pop();
            continue;
         }

      case ExpressionGraph::INTEGER:
         switch( top.first )
         {
         case 0:
            stream << "floor(";
            ++top.first;
            stack.emplace( 0, node->child1 );
            continue;

         case 1:
            stream << ")";
            stack.pop();
            continue;
         }

      case ExpressionGraph::LN:
         switch( top.first )
         {
         case 0:
            stream << "log(";
            ++top.first;
            stack.emplace( 0, node->child1 );
            continue;

         case 1:
            stream << ")";
            stack.pop();
            continue;
         }

      case ExpressionGraph::UMINUS:
         switch( top.first )
         {
         case 0:
            stream << "-(";
            ++top.first;
            stack.emplace( 0, node->child1 );
            continue;

         case 1:
            stream << ")";
            stack.pop();
            continue;
         }

      case ExpressionGraph::NOT:
         switch( top.first )
         {
         case 0:
            stream << "not (";
            ++top.first;
            stack.emplace( 0, node->child1 );
            continue;

         case 1:
            stream << ")";
            stack.pop();
            continue;
         }

      case ExpressionGraph::SQRT:
         switch( top.first )
         {
         case 0:
            stream << "sqrt(";
            ++top.first;
            stack.emplace( 0, node->child1 );
            continue;

         case 1:
            stream << ")";
            stack.pop();
            continue;
         }

      case ExpressionGraph::PLUS:
         switch( top.first )
         {
         case 0:
            ++top.first;
            stack.emplace( 0, node->child1 );
            continue;

         case 1:
            stream << "+";
            stack.pop();
            stack.emplace( 0, node->child2 );
            continue;
         }

      case ExpressionGraph::MINUS:
         switch( top.first )
         {
         case 0:
            ++top.first;
            stack.emplace( 0, node->child1 );
            continue;

         case 1:
            stream << "-(";
            ++top.first;
            stack.emplace( 0, node->child2 );
            continue;

         case 2:
            stream << ")";
            stack.pop();
            continue;
         }

      case ExpressionGraph::MULT:
         switch( top.first )
         {
         case 0:
            stream << "(";
            ++top.first;
            stack.emplace( 0, node->child1 );
            continue;

         case 1:
            stream << ")*(";
            ++top.first;
            stack.emplace( 0, node->child2 );
            continue;

         case 2:
            stream << ")";
            stack.pop();
            continue;
         }

      case ExpressionGraph::DIV:
         switch( top.first )
         {
         case 0:
            stream << "(";
            ++top.first;
            stack.emplace( 0, node->child1 );
            continue;

         case 1:
            stream << ")/(";
            ++top.first;
            stack.emplace( 0, node->child2 );
            continue;

         case 2:
            stream << ")";
            stack.pop();
            continue;
         }

      case ExpressionGraph::AND:
         switch( top.first )
         {
         case 0:
            stream << "(";
            ++top.first;
            stack.emplace( 0, node->child1 );
            continue;

         case 1:
            stream << " and ";
            ++top.first;
            stack.emplace( 0, node->child2 );
            continue;

         case 2:
            stream << ")";
            stack.pop();
            continue;
         }

      case ExpressionGraph::OR:
         switch( top.first )
         {
         case 0:
            stream << "(";
            ++top.first;
            stack.emplace( 0, node->child1 );
            continue;

         case 1:
            stream << " or ";
            ++top.first;
            stack.emplace( 0, node->child2 );
            continue;

         case 2:
            stream << ")";
            stack.pop();
            continue;
         }

      case ExpressionGraph::L:
         switch( top.first )
         {
         case 0:
            stream << "(";
            ++top.first;
            stack.emplace( 0, node->child1 );
            continue;

         case 1:
            stream << " < ";
            ++top.first;
            stack.emplace( 0, node->child2 );
            continue;

         case 2:
            stream << ")";
            stack.pop();
            continue;
         }

      case ExpressionGraph::LE:
         switch( top.first )
         {
         case 0:
            stream << "(";
            ++top.first;
            stack.emplace( 0, node->child1 );
            continue;

         case 1:
            stream << " <= ";
            ++top.first;
            stack.emplace( 0, node->child2 );
            continue;

         case 2:
            stream << ")";
            stack.pop();
            continue;
         }

      case ExpressionGraph::G:
         switch( top.first )
         {
         case 0:
            stream << "(";
            ++top.first;
            stack.emplace( 0, node->child1 );
            continue;

         case 1:
            stream << " > ";
            ++top.first;
            stack.emplace( 0, node->child2 );
            continue;

         case 2:
            stream << ")";
            stack.pop();
            continue;
         }

      case ExpressionGraph::GE:
         switch( top.first )
         {
         case 0:
            stream << "(";
            ++top.first;
            stack.emplace( 0, node->child1 );
            continue;

         case 1:
            stream << " >= ";
            ++top.first;
            stack.emplace( 0, node->child2 );
            continue;

         case 2:
            stream << ")";
            stack.pop();
            continue;
         }

      case ExpressionGraph::EQ:
         switch( top.first )
         {
         case 0:
            stream << "(";
            ++top.first;
            stack.emplace( 0, node->child1 );
            continue;

         case 1:
            stream << " eq ";
            ++top.first;
            stack.emplace( 0, node->child2 );
            continue;

         case 2:
            stream << ")";
            stack.pop();
            continue;
         }

      case ExpressionGraph::NEQ:
         switch( top.first )
         {
         case 0:
            stream << "(";
            ++top.first;
            stack.emplace( 0, node->child1 );
            continue;

         case 1:
            stream << " <> ";
            ++top.first;
            stack.emplace( 0, node->child2 );
            continue;

         case 2:
            stream << ")";
            stack.pop();
            continue;
         }

      case ExpressionGraph::LOG:
         switch( top.first )
         {
         case 0:
            stream << "log(";
            ++top.first;
            stack.emplace( 0, node->child1 );
            continue;

         case 1:
            stream << ")/log(";
            ++top.first;
            stack.emplace( 0, node->child2 );
            continue;

         case 2:
            stream << ")";
            stack.pop();
            continue;
         }

      case ExpressionGraph::POWER:
         switch( top.first )
         {
         case 0:
            stream << "(";
            ++top.first;
            stack.emplace( 0, node->child1 );
            continue;

         case 1:
            stream << ")**(";
            ++top.first;
            stack.emplace( 0, node->child2 );
            continue;

         case 2:
            stream << ")";
            stack.pop();
            continue;
         }

      case ExpressionGraph::MIN:
         switch( top.first )
         {
         case 0:
            stream << "min(";
            ++top.first;
            stack.emplace( 0, node->child1 );
            continue;

         case 1:
            stream << ", ";
            ++top.first;
            stack.emplace( 0, node->child2 );
            continue;

         case 2:
            stream << ")";
            stack.pop();
            continue;
         }

      case ExpressionGraph::MAX:
         switch( top.first )
         {
         case 0:
            stream << "max(";
            ++top.first;
            stack.emplace( 0, node->child1 );
            continue;

         case 1:
            stream << ", ";
            ++top.first;
            stack.emplace( 0, node->child2 );
            continue;

         case 2:
            stream << ")";
            stack.pop();
            continue;
         }

      case ExpressionGraph::MODULO:
         switch( top.first )
         {
            case 0:
               stream << "mod(";
               ++top.first;
               stack.emplace( 0, node->child1 );
               continue;
            case 1:
               stream << ", ";
               ++top.first;
               stack.emplace( 0, node->child2 );
               continue;
            case 2:
               stream << ")";
               stack.pop();
               continue;
         }

      case ExpressionGraph::INITIAL:
         switch( top.first )
         {
         case 0:
            if( initial )
            {
               stack.pop();
            }
            else
            {
               ++top.first;
               initial = true;
            }

            stack.emplace( 0, node->child1 );
            continue;

         case 1:
            initial = false;
            stack.pop();
            continue;
         }

      case ExpressionGraph::APPLY_LOOKUP:
      {
         //get lookup data
         LookupData& lkpData = lkpData_[node->child1->lookup_table];

         if( initial )
         {
            stream << boost::lexical_cast<std::string>( node->child1->lookup_table->operator()( node->child2->value ) );
            continue;
         }

         if( lkpData.type == LookupFormulationType::SPLINE )
         {
            //for spline lookups call lookuplib
            switch( top.first )
            {
            case 0:
               stream << "Lookup(";
               ++top.first;
               stack.emplace( 0, node->child2 );
               continue;

            case 1:
               stream << ", lkp_";
               stream << escape_string( lkpData.name );
               stream << ")";
               stack.pop();
               continue;
            }
         }
         else
         {
            std::string lkpName  = escape_string(lkpData.name);
            stream << "sum(lkp_" << lkpName << "_points, lkp_" << lkpName << sos2LkpIds_[node] << "_lambda(" << getVarSets()
                   << ", lkp_" << lkpName << "_points)*lkp_" << lkpName  << "_Y(lkp_" << lkpName << "_points) )";
            stack.pop();
            continue;
         }
      }

      case ExpressionGraph::DELAY_FIXED:
      {
         double timestep = exprGraph_.getNode( Symbol( "TIME STEP" ) )->value;
         double delaytime = std::min( node->child2->value, timestep );
         int dt = std::ceil( delaytime / timestep );

         switch( top.first )
         {
         case 0:
         {
            if( initial )
            {
               stack.pop();
               stack.emplace( 0, node->child3 );
               continue;
            }

            std::string t = getSets( {"t"} );
            controlSet( "t" );
            std::string tt = getSets( {"t"} );
            stream << "sum( " << tt << "$(ord(" << tt << ") eq ord(" << t << ") - " << dt << "), ";
            //input
            ++top.first;
            stack.emplace( 0, node->child1 );
            continue;
         }

         case 1:
            releaseSet( "t" );
            stream << ")+(";
            //initial
            initial = true;
            ++top.first;
            stack.emplace( 0, node->child3 );
            continue;

         case 2:
            initial = false;
            stream << ")$( ord(" << getSets( {"t"} ) << ") le " << dt << " )";
            stack.pop();
            continue;
         }
      }

      case ExpressionGraph::CONTROL:
         //err: control should have a symbol and thus be handled before this switch
         assert( false );

      case ExpressionGraph::LOOKUP_TABLE:
         //TODO
         assert( false );

      case ExpressionGraph::NIL:
         //err: nil node should not exists in expr graph and should be handled by front end
         assert( false );
      }
   }
   while( !stack.empty() );
}

void GamsGenerator::emitGams( std::ostream& stream )
{
   std::vector<std::pair<int, std::string>> parameters;
   std::vector<std::pair<int, std::string>> varValues;
   std::vector<std::pair<int, std::string>> equationDeclarations;
   std::vector<std::pair<int, std::string>> equations;

   auto parameter = [&parameters]( int level, std::ostringstream & ss )
   {
      parameters.emplace_back( level, ss.str() );
      ss.clear();
      ss.str( std::string() );
   };

   auto varValue = [&varValues]( int level, std::ostringstream & ss )
   {
      varValues.emplace_back( level, ss.str() );
      ss.clear();
      ss.str( std::string() );
   };

   auto equationDeclaration = [&equationDeclarations]( int level, std::ostringstream & ss )
   {
      equationDeclarations.emplace_back( level, ss.str() );
      ss.clear();
      ss.str( std::string() );
   };

   auto equation = [&equations]( int level, std::ostringstream & ss )
   {
      equations.emplace_back( level, ss.str() );
      ss.clear();
      ss.str( std::string() );
   };

   std::ostringstream ss;

   stream << "$offdigit\n";

   int lkp_line = 0;
   bool spline_type = has_spline_type( lkpData_ );

   if( spline_type )
   {
      stream << "$onecho > conopt.opt\n"
             << "lkdebg 0\n"
             << "$offecho\n"
             << "$onecho > lookups.dat\n"
             << "max_mixed_err     = 0.01\n"
             << "mixed_err_delta   = 1\n"
             << "min_knot_distance = 1e-6\n"
             << "obj_tolerance     = 1e-7\n";
   }

   //handle lookups
   for( LkpTypePair & pair : lkpData_ )
   {
      std::string lkp_name = escape_string( pair.second.name );

      //stream lookup data into lookups.dat and count lin
      if( pair.second.type == LookupFormulationType::SPLINE )
      {

         bool space = false;

         for( const auto & point : * ( pair.first ) )
         {
            if( space )
               stream << " ";

            space = true;
            double x = point.get<0>();
            double y = point.get<1>();
            stream << boost::lexical_cast<std::string>( x ) << " "
                   << boost::lexical_cast<std::string>( y );
         }

         stream << "\n";

         ss << "Parameter lkp_" << lkp_name << " / " << lkp_line++ << " /;\n";
         parameter( 0,  ss );
      }
      else
      {
         ss << "Parameter lkp_" << lkp_name << "_X(lkp_" << lkp_name << "_points) /";
         int i = 1;
         ss << "\n\t" << i++ << "\t-" << boost::lexical_cast<std::string>( lkp_infty_ );

         for( double xval : pair.first->getXvals() )
         {
            ss << "\n\t" << i++ << "\t" << boost::lexical_cast<std::string>( xval );
         }

         ss << "\n\t" << i++ << "\t" << boost::lexical_cast<std::string>( lkp_infty_ );
         ss << " /;\n";
         ss << "Parameter lkp_" << lkp_name << "_Y(lkp_" << lkp_name << "_points) /";
         i = 1;
         ss << "\n\t" << i++ << "\t" << boost::lexical_cast<std::string>( pair.first->getYvals().front() );

         for( double yval : pair.first->getYvals() )
         {
            ss << "\n\t" << i++ << "\t" << boost::lexical_cast<std::string>( yval );
         }

         ss << "\n\t" << i++ << "\t" << boost::lexical_cast<std::string>( pair.first->getYvals().back() );
         ss << " /;\n";
         parameter( 0, ss );
      }

   }

   if( spline_type )
   {
      stream << "$offecho\n";
      stream << "$funclibin liblookup %GAMS.workdir%liblookup.so\n"
             << "function Lookup / liblookup.Lookup /;\n\n";
   }

   double final_time = exprGraph_.getNode( Symbol( "FINAL TIME" ) )->value;
   double initial_time = exprGraph_.getNode( Symbol( "INITIAL TIME" ) )->value;
   double time_step = exprGraph_.getNode( Symbol( "TIME STEP" ) )->value;
   //stream sets
   stream << "Set t time periods / " << 0 << "*" << ( final_time - initial_time ) / time_step << " /;\n"
          << "Set tfirst(t) first period;\n"
          << "Set tlast(t) last period;\n\n";

   if( tableau_.getName() != ButcherTableau::EULER )
   {
      stream << "Set p discretization sampling points / 0*" << tableau_.columns() << " /;\n"
             << "Table coeff(p, p) discretization coefficients\n";

      for( int i = 1; i <= tableau_.columns(); ++i )
         stream << "\t" << i;

      stream << "\n";

      for( int i = 1; i <= tableau_.columns(); ++i )
      {
         stream << i;

         for( int j = 1; j <= tableau_.columns(); ++j )
         {
            stream << "\t" << boost::lexical_cast<std::string>( tableau_[i - 1][j - 1] );
         }

         stream << "\n";
      }

      stream << "Parameter weight(p) discretization point weights /\n";

      for( int i = 1; i <= tableau_.columns(); ++i )
         stream << "\t" << i << "\t" << boost::lexical_cast<std::string>( tableau_[tableau_.rows() - 1][i - 1] ) << std::endl;

      stream << "\t/;\n\n";
   }

   stream << "tfirst(t) = yes$(ord(t) eq 1);\n"
          << "tlast(t)  = yes$(ord(t) eq card(t));\n";
   std::set<std::size_t> control_step_sizes;

   for( auto & pair : exprGraph_.getSymbolTable() )
   {
      if( pair.second->op == ExpressionGraph::CONTROL )
      {
         if( pair.second->control_size )
         {
            control_step_sizes.insert( pair.second->control_size );
         }
      }
   }

   for( std::size_t control_step_size : control_step_sizes )
   {
      stream << "set t" << control_step_size << " time periods of " << control_step_size << " time steps / "
             << 0 << "*" << int( ( final_time - initial_time ) / time_step / control_step_size )
             << " /;\n";
   }

   for( const LkpTypePair & pair : lkpData_ )
   {
      if( pair.second.type == LookupFormulationType::SOS2 )
      {
         sdo::LookupTable& table = *( pair.first );
         stream << "set lkp_" << escape_string( pair.second.name ) << "_points / 1*" << table.size() + 2 << " /;\n";
      }
   }

   //create missing symbols
   createDivisionGuards( varValues );
   createStateSymbols();
   //fill map 'sos2LkpIds_'
   indexSos2Lookups();
   //create epsilon and time as parameter
   ss << "Parameter EPSILON / 1e-9 /;\n";
   parameter( 0, ss );

   ss << "Parameter TIME(t);\n"
      << "\tTIME(t) = INITIALTIME+(ord(t)-1)*TIMESTEP;\n";
   parameter( exprGraph_.getTimeNode()->level, ss );

   stream << "\n";

   //loop over all symbols. Emit variable declarations directly. Add parameters equations etc. to the corresponding vectors
   for( auto & entry : exprGraph_.getSymbolTable() )
   {
      if(entry.second->op == ExpressionGraph::LOOKUP_TABLE)
         continue;
      std::string var = escape_string( entry.first );
      std::string comment;
      {
         auto range = exprGraph_.getComments( entry.first );

         if( range.begin() != range.end() )
         {
            std::ostringstream stringstream;
            stringstream << " \"" << range.begin()->second << '"';
            comment = stringstream.str();
         }
      }

      switch( entry.second->type )
      {
      case ExpressionGraph::DYNAMIC_NODE: //node that depends on time and on states/controls
         if( entry.second->op == ExpressionGraph::CONTROL ) //control -> create variables and bounds
         {
            switch( entry.second->control_size )
            {
            case 0:
               stream << "Variable " << var  << comment << ";\n";

               if( entry.second->child1 )
               {
                  ss << var << ".lo = " << boost::lexical_cast<std::string>( entry.second->child1->value ) << ";\n";
                  varValue( 0, ss );
               }

               if( entry.second->child2 )
               {
                  ss << var << ".l = " << boost::lexical_cast<std::string>( entry.second->child2->value ) << ";\n";
                  varValue( 0, ss );
               }

               if( entry.second->child3 )
               {
                  ss << var << ".up = " << boost::lexical_cast<std::string>( entry.second->child3->value ) << ";\n";
                  varValue( 0, ss );
               }

               break;

            case 1:
               stream << "Variable " << var  << "(t)" << comment << ";\n";

               if( entry.second->child1 )
               {
                  ss << var << ".lo(t) = " << boost::lexical_cast<std::string>( entry.second->child1->value ) << ";\n";
                  varValue( 0, ss );
               }

               if( entry.second->child2 )
               {
                  ss << var << ".l(t) = " << boost::lexical_cast<std::string>( entry.second->child2->value ) << ";\n";
                  varValue( 0, ss );
               }

               if( entry.second->child3 )
               {
                  ss << var << ".up(t) = " << boost::lexical_cast<std::string>( entry.second->child3->value ) << ";\n";
                  varValue( 0, ss );
               }

               break;

            default:
               stream << "Variable " << var  << "(t"  << entry.second->control_size << ")" << comment << ";\n";

               if( entry.second->child1 )
               {
                  ss << var << ".lo(t" << entry.second->control_size << ") = " << boost::lexical_cast<std::string>( entry.second->child1->value ) << ";\n";
                  varValue( 0, ss );
               }

               if( entry.second->child2 )
               {
                  ss << var << ".l(t" << entry.second->control_size << ") = " << boost::lexical_cast<std::string>( entry.second->child2->value ) << ";\n";
                  varValue( 0, ss );
               }

               if( entry.second->child3 )
               {
                  ss << var << ".up(t" << entry.second->control_size << ") = " << boost::lexical_cast<std::string>( entry.second->child3->value ) << ";\n";
                  varValue( 0, ss );
               }
            }
         }
         else     //no control -> either state or algebraic
         {
            stream << "Variable " << var << "(" << getVarSets() << ")" << comment << ";\n";
            ss << "Equation eq_" << var << "(" << getVarSets() << ");\n";
            equationDeclaration( entry.second->level, ss );

            if( entry.second->op == ExpressionGraph::INTEG ) //for states create steps for discretization and initial values
            {
               if( tableau_.getName() == ButcherTableau::EULER )
               {
                  ss << "eq_" << var <<  "(t+1) ..\n\t" << var << "(t+1) =e= "
                     << var << "(t) + TIMESTEP * ( ";
                  translate( ss, entry.second->child1, false  );
                  ss << " );\n";
                  equation( entry.second->level, ss );
               }
               else
               {
                  //declare equation for Integration step which defines the value of state var as
                  //weighted sum of the intermediate time steps according to the butcher tableau
                  ss << "Equation eq_" << var << "IntegStep(" << getVarSets() << ");\n";
                  equationDeclaration( entry.second->level, ss );

                  //build definition of the integration step
                  ss << "eq_" << var << "IntegStep(" << getSets( {SetIndex( "t", 1 ), SetIndex::First( "p" )} ) << ") ..\n\t" << var << "(" << getSets( {SetIndex( "t", 1 ), SetIndex::First( "p" )} ) << ") =e= "
                     << var << "(" << getSets( {"t", SetIndex::First( "p" )} ) << ")+TIMESTEP*sum(p$( ord(p) > 1 ), weight(p)*(";
                  translate( ss, entry.second->child1, false );
                  ss << "));\n";
                  equation( entry.second->level, ss );
                  //define the intermediate steps according to the coefficients in the butcher tableau

                  ss << "eq_" << var << "(" << getVarSets() << ")$( ord(p) > 1 ) ..\n\t" << var << "(" << getVarSets() << ") =e= "
                     << var << "(" << getSets( {"t", SetIndex::First( "p" )} ) << ")+TIMESTEP*sum(pp$( ord(pp) > 1 ), coeff(p, pp)*(";
                  controlSet( "p" );
                  translate( ss, entry.second->child1, false );
                  ss << "));\n";
                  releaseSet( "p" );
                  equation( entry.second->level, ss );
               }

               if( entry.second->init == ExpressionGraph::CONSTANT_INIT )
               {
                  ss << var << ".fx(" << getInitialSets() << ") = ";
                  translate( ss, entry.second->child2, false, true );
                  ss << ";\n";
                  varValue( entry.second->level, ss );
               }
               else     //initial value is controled so enforce it by equation
               {
                  ss << "Equation eq_" << var << "Init;\n";
                  equationDeclaration( entry.second->level, ss );
                  ss << "eq_" << var << "Init ..\n\t" << var << "(" << getInitialSets() << ") =e= ";
                  translate( ss, entry.second->child2, false, true );
                  ss << ";\n";
                  equation( entry.second->level, ss );
               }
            }
            else      //no integ -> just add definition to equations
            {
               ss << "eq_" << var << "(" << getVarSets() << ") ..\n\t" << var << "(" << getVarSets() << ") =e= ";
               translate( ss, entry.second );
               ss << ";\n";
               equation( entry.second->level, ss );
            }
         } //end of case DYNAMIC_NODE

         break;

      case ExpressionGraph::STATIC_NODE:   // node that does depend on time but is constant at each time -> use a parameter
      {
         ss << "Parameter " << var << "(t)" << comment << ";\n";
         ss << "\t" << var << "(t) = ";
         translate( ss, entry.second );
         ss << ";\n";
         parameter( entry.second->level, ss );
         break;
      }

      case ExpressionGraph::CONSTANT_NODE: // constant node -> use a parameter
         ss << "Parameter " << var << comment << " / " << boost::lexical_cast<std::string>( entry.second->value ) << " /;\n";
         parameter( entry.second->level, ss );
         break;

      case ExpressionGraph::UNKNOWN:
         assert( false );
      }
   }

   //Add equations and variables for sos2 lookups that were found
   for( std::pair<ExpressionGraph::Node* const, int>& entry : sos2LkpIds_ )
   {
      LookupData& lkpData = lkpData_[entry.first->child1->lookup_table];
      std::string lkpName = escape_string(lkpData.name);
      stream << "sos2 Variable lkp_" << lkpName << entry.second << "_lambda(" << getVarSets() << ", lkp_" << lkpName << "_points);\n";

      ss << "Equation eq_lkp_" << lkpName << entry.second  << "_norm(" << getVarSets() << ");\n"
         << "Equation eq_lkp_" << lkpName << entry.second << "_arg(" << getVarSets() << ");\n";
      equationDeclaration( entry.first->level, ss );

      ss << "eq_lkp_" << lkpName << entry.second << "_norm(" << getVarSets() << ") ..\n\t"
         << "sum(lkp_" << lkpName << "_points, lkp_" << lkpName << entry.second << "_lambda(" << getVarSets() << ", lkp_" << lkpName << "_points)) =e= 1;\n";
      ss << "eq_lkp_" << lkpName << entry.second << "_arg(" << getVarSets() << ") ..\n\t";
      translate( ss, entry.first->child2 );
      ss << " =e= sum(lkp_" << lkpName << "_points, lkp_" << lkpName << entry.second << "_lambda(" << getVarSets() << ", lkp_" << lkpName << "_points)*lkp_" << lkpName  << "_X(lkp_" << lkpName << "_points) );\n";
      equation( entry.first->level, ss );
   }

   if(!objective_.empty()) {
      stream << "Variable objective;\n";
      ss << "Equation eq_objective;\n";
      equationDeclaration( std::numeric_limits<int>::max(), ss );
      ss << "eq_objective ..\n\tobjective =e= ";
      bool first = true;

      for( Objective::Summand & s : objective_.getSummands() )
      {
         if( !first )
            ss << "+";
         bool discrSet = sets_.find("p") != sets_.end();
         if( s.type == Objective::Summand::MAYER ) {
            if(discrSet)
               ss << "sum( (t, p)$(ord(p) eq 1 and ord(t) eq card(t)), ";
            else
               ss << "sum( t$(ord(t) eq card(t)), ";
         } else {
            if(discrSet)
               ss << "sum( (t, p)$(ord(p) eq 1), ";
            else
               ss << "sum( t, ";
         }

         translateSymbol(ss, s.variable);

         ss << ")";
      }

      ss << ";\n";
      equation( std::numeric_limits<int>::max(), ss );
   }
   //now variables have been declared

   // sort all data by their level to enforce that they are emitted
   // in the proper order, i.e.
   //    'Parameter A(t) = B(t)+5;'
   // is emitted after the definition of B is emitted. If that is not
   // the case gams will complain.
   std::sort( varValues.begin(), varValues.end() );
   std::sort( equations.begin(), equations.end() );
   std::sort( equationDeclarations.begin(), equationDeclarations.end() );
   std::sort( parameters.begin(), parameters.end() );

   stream << "\n";

   for( auto & s : sets_ )
   {
      const std::string& setName = s.first;
      int nAliases = s.second.second;

      for( int n = 1; n <= nAliases; ++n )
      {
         stream << "alias(" << setName << ", " << setName;

         for( int i = 0; i < n; ++i )
         {
            stream << setName;
         }

         stream << ");\n";
      }
   }

   //now emit everything that has been generated
   stream << "\n";

   for( auto & pair : parameters )
   {
      stream << pair.second;
   }

   stream << "\n";

   for( auto & pair : varValues )
   {
      stream << pair.second;
   }

   stream << "\n";

   for( auto & pair : equationDeclarations )
   {
      stream << pair.second;
   }

   stream << "\n";

   for( auto & pair : equations )
   {
      stream << pair.second;
   }

   stream << "\n";

   if(!objective_.empty()) {
      stream << "Model m / all /;\n";

      if(has_spline_type(lkpData_))
         stream << "m.optfile = 1;\n";
      if(objective_.isMinimized())
         stream << "Solve m min objective ";
      else
         stream << "Solve m max objective ";
      if(has_sos2_type(lkpData_))
         stream << "using minlp;\n";
      else
         stream << "using nlp;\n";

   }
}

void GamsGenerator::addObjective( Objective obj )
{
   objective_ = std::move( obj );
}

void GamsGenerator::addArbitraryObjective()
{
   for( auto & entry : exprGraph_.getSymbolTable() )
   {
      if( entry.second->op == ExpressionGraph::INTEG )
      {
         objective_.addSummand( Objective::Summand::MAYER, entry.first, 1. );
         return;
      }
   }
   objective_.setMinimized(true);
}


}
