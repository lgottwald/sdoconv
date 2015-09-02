#ifndef  _GAMS_HPP_
#define  _GAMS_HPP_

#include <sdo/ExpressionGraph.hpp>
#include <sdo/ButcherTableau.hpp>
#include <sdo/Objective.hpp>
#include <sdo/LookupTable.hpp>
#include <unordered_map>
#include <vector>
#include <ostream>
#include <string>



namespace gams {

/**
 * Forward declare class SetIndex.
 */
class SetIndex;

using namespace sdo;

/**
 * \brief Enum type to indicate how a lookup function should be formulated in the gams model.
 */
enum class LookupFormulationType {
   SPLINE, //< This type means, that the GamsGenerator will an extrinsic function that approximates the lookups by a spline function to evaluate the lookup in gams
   SOS2 //< This type means, that the GamsGenerator will formulate the lookup in gams using sos2 variables to model the piecewise linear function described by the lookup values.
};


/**
 * \brief Relevant data of a lookup.
 * 
 * Struct to represent relevant informations about a lookup table. I.e.
 * how it should be formulated and it's name.
 */
struct LookupData {
   LookupData() = default;
   LookupData(Symbol name, LookupFormulationType type) : name(std::move(name)), type(type) , usages(0) {}

   /**
    * \brief The name of the lookup.
    * 
    * The name of the lookup table. For lookups defined by A = WITH_LOOKUP( B, ((0,0)(1,1)) )
    * the name 'A' can be used although it is the same for the symbol defined by applying
    * the lookup on 'B'.
    * For lookups defined seperately it is the lookups name, e.g. 'C' if the lookup is defined as
    * C((0,0),(1,1)).
    */
   Symbol name;
   LookupFormulationType type = LookupFormulationType::SPLINE; //< formulation type of lookup, i.e. SPLINE or SOS2
   int usages; //< will be used during gams translation to count usages of sos2 lookups 
};


/**
 * \brief Generates gams from the intermediate represenation of a mdl file.
 * 
 * The GamsGenerator is used to emit gams for a mdl file to a output stream. The mdl file
 * is given as a sdo::ExpressionGraph together with an optional sdo::Objective generated
 * from a vpd file to also generate an objective in the gams output.
 */
class GamsGenerator
{
public:


   /**
    * \brief Construct a gams generator for a given expression graph.
    * 
    * \param exprGraph the sdo::ExpressionGraph to generate the gams output. It is stored by reference and not copied.
    * \param tableau an enum value of sdo::ButcherTableau::Name to identify the discretization method.
    * \param lkpType the value for the boundaries of the sos2 lookup. The lookup argument should stay in [-lkp_infty,lkp_infty] during optimization.
    */
   GamsGenerator(
      sdo::ExpressionGraph& exprGraph,
      sdo::ButcherTableau::Name tableau = sdo::ButcherTableau::RUNGE_KUTTA_2,
      LookupFormulationType lkpType = LookupFormulationType::SPLINE
   ) : exprGraph_( exprGraph ), lkp_infty_(1e4)
   {
      initTableau(tableau);
      setLookupFormulationTypes(lkpType);
      //create symbols for all states
      createStateSymbols();
      //fill map 'sos2LkpIds_'
      indexSos2Lookups();
   }

   /**
    * \brief Generate and emit gams output.
    * 
    * \param stream the output stream used to emit the generated gams output.
    */
   void emitGams( std::ostream& stream );

   /**
    * \brief Generates an arbritrary objective function containing some state.
    * 
    * Some arbitrary state variable is selected and then an objective is genererated
    * minimizing its value at the final time. The objective is not meant to make sense
    * for the model, but to be able to directly execute the generated output in gams.
    */
   void addArbitraryObjective();

   /**
    * \brief Add an objective function to gams output
    * 
    * \param obj the objective function to add.
    */
   void addObjective( sdo::Objective obj );

   /**
    * \brief Sets formulation types for all lookups to given value.
    * 
    * \param type the formulation type for all lookups.
    */
   void setLookupFormulationTypes(LookupFormulationType type);


   /**
    * \brief Set formulation type for specific lookup table.
    */
   void setLookupFormulationType(LookupTable* lookup, LookupData data );

   /**
    * \brief Set boundary value used for sos2 formulation of lookups.
    * 
    * \param val the boundary value.
    */
   void setSos2LookupBoundary( double val ) {
      lkp_infty_ = val;
   }

private:
   /**
    * Creates lower bounds slightly above zero for all expressions that are divisors
    * 
    * \param varValues  vector containing the bounds and fixed values of variables together with their levels.
    */
   void createDivisionGuards(std::vector<std::pair<int, std::string>> &varValues);
   /**
    * Creates symbols for all states since SMOOTH DELAY etc. may contain hidden states that do not have a symbol in the mdl file.
    */
   void createStateSymbols();

   /**
    * Creates an id for each call of an sos2 lookup. The three calls lookup(a+b) lookup(b+a) lookup(c)
    * will get two id's since the first two calls are identical.
    * This avoids the creation of unnecessary sos2 variables.
    */
   void indexSos2Lookups();
   /**
    * Initializes the butcher tableau with the given one.
    * 
    * \param tableau enum value of type sdo::ButcherTableau::Name identifying the butcher tableau.
    */
   void initTableau(sdo::ButcherTableau::Name tableau);

   /**
    * \brief Translates a node in the expression graph to gams.
    * 
    * Translates an expression represented by its node in the sdo::ExpressionGraph
    * to its gams expression. The gams will represent the value of the expression
    * for each time using a set that represents the time steps. If the initial flag is
    * set to true, the expression will represent the initial value.
    * 
    * 
    * 
    * \param stream The gams output will be emitted to this output stream.
    * \param root The expression graph node that will be translated.
    * \param def If true the given node is translated as its definition, i.e. its symbol is not used but expanded.
    *            This is required to avoid definitions like "some_var(t) =e= some_var(t)"
    * \param initial if true the gams output will represent the initial value of the expression
    */
   void translate(std::ostream& stream, ExpressionGraph::Node* root, bool def = true, bool initial = false);

   /**
    * \brief Translates a node in the expression graph to gams.
    * 
    * Translates a symbol to its gams call. Uses Information from the expressiongraph
    * and adds the corresponding sets but uses the currently controled alias of it.
    * E.g. if the symbol is a constant it will just be escaped to form a proper gams identifier
    * but if it is a state the time set and discretization set will be appended.
    * 
    * \param stream The gams identifier with its sets will be emitted to this output stream.
    * \param symbol The symbol that shouldbe translated
    * \param initial if true the gams output will represent the initial value of the symbol
    */
   void translateSymbol( std::ostream& stream, Symbol s, bool initial = false );

   /**
    * \brief Make the set with the given name controled, so that getSets will give an alias string.
    * 
    * If getSets({{"t", 0}) returns t it will return tt after one call to controlSet("t")
    * and ttt after two calls to controlSet("t"). Afterwards calling releaseSet("t") will recover
    * to the state before the last call of controlSet("t").
    * 
    * \param set the name of the set to control
    */
   void controlSet( std::string set );

   /**
    * \brief Revert the changes made by controlSet(set)
    * 
    * \param set the name of the set
    */
   void releaseSet( std::string set );

   /**
    * \brief Create a new set that can then be controled and released by {control, release}Set
    * 
    * \param set the name of the set
    */
   void createSet( std::string set );

   /**
    * \brief Get a string representing the sets described by the given initializer list.
    * 
    * The initializer list contains pairs of set names with offsets. So 
    * getSets({{"t",1}, "p"}) will return 't+1, p'. If a set is controled
    * before the call to this function an alias will be used instead of the set
    * name. So the string returned for the above example might also be 'tt+1, p'.
    * The first value of a set is assumed to be '0'.
    * 
    * \param list the sets with their offsets that should be in the string
    */
   std::string getSets(std::initializer_list<SetIndex> list) const;

   /**
    * \brief Get first set of time and discretization, i.e. "'0'" for euler method or else "'0', '0'".
    */
   std::string getInitialSets() const;

   /**
    * \brief Get current time and discretization, i.e. "t" for euler method or else "t, p".
    * 
    * Result can also be the currently controled alias of the sets t and p , e.g. tt or pp.
    */
   std::string getVarSets() const;
 

   sdo::ButcherTableau tableau_;
   std::unordered_map<LookupTable*, LookupData> lkpData_;
   std::unordered_map<ExpressionGraph::Node*, int> sos2LkpIds_;
   sdo::ExpressionGraph& exprGraph_;
   sdo::Objective objective_;
   double lkp_infty_;
   std::unordered_map<std::string, std::pair<int,int>> sets_;
   
};


}

#endif //  _GAMS_GAMSGENERATOR_HPP_
