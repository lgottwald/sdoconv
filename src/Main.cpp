#include <iostream>
#include <sstream>
#include <fstream>
#include <sdo/Parsers.hpp>
#include <sdo/ButcherTableau.hpp>
#include "GamsGenerator.hpp"
#include <boost/program_options.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <vector>
#include <string>

int main( int argc, char const* argv[] )
{
   namespace po = boost::program_options;
   po::options_description desc( "Allowed options" );
   desc.add_options()
   ( "help,h", "produce help message" )
   ( "discretization-method,d", po::value<std::string>()->default_value( "rk2" ), "Method used for discretization. Available: euler, rk2, rk3, rk4, imid2, igl4" )
   ( "input-files", po::value< std::vector<std::string> >(), "Input files" )
   ( "output-file,o", po::value<std::string>(), "File to write gams output. If not set gams is written to stdout." )
   ( "lookup-type,l", po::value<std::string>()->default_value( "interactive" ), "Formulation type of lookups. sos2, spline or interactive" )
   ( "lookup-infinity,f", po::value<double>()->default_value( 1e5 ), "Value for lookup boundaries. Too small values may yield an infeasible gams-model. Too big values may result in numerical instabilities." )
   ;
   po::positional_options_description p;
   p.add( "input-files", -1 );

   po::variables_map vm;

   try
   {
      po::store( po::command_line_parser( argc, argv ).
                 options( desc ).positional( p ).run(), vm );
      po::notify( vm );
   }
   catch( boost::program_options::required_option& e )
   {
      std::cerr << "Error: " << e.what() << "\n";
      exit( 0 );
   }
   catch( boost::program_options::error& e )
   {
      std::cerr << "Error: " << e.what() << "\n";
      exit( 0 );
   }

   if( vm.count( "help" ) )
   {
      std::cout << desc;
      exit( 0 );
   }

   if( !vm.count( "input-files" ) )
   {
      std::cerr << "Error: no input file specified\n" << desc;
      exit( 0 );
   }
   std::unique_ptr<std::ofstream> file;
   if(vm.count( "output-file" ))
      file = std::unique_ptr<std::ofstream>{ new std::ofstream(vm["output-file"].as<std::string >()) };

   std::ostream &out = file ? *file : std::cout;
   if(!out.good())
   {
      std::cerr << "Error: unable to write to file '" << vm["output-file"].as<std::string >() << "'\n";
      exit( 0 );
   }

   std::unordered_map<sdo::LookupTable*, gams::LookupData> lkpTypes;

   std::vector<std::string> input_files = vm["input-files"].as< std::vector<std::string> >();
   std::string discretization_method_name = vm["discretization-method"].as<std::string >();
   std::string lookup_type = vm["lookup-type"].as<std::string>();
   if(lookup_type != "sos2" && lookup_type != "spline" && lookup_type != "interactive")
   {
      std::cerr << "Error: unknown lookup-type '" << lookup_type << "'\n";
      exit( 0 );
   }

   sdo::ButcherTableau::Name discretization_method;

   if( discretization_method_name == "euler" )
      discretization_method = sdo::ButcherTableau::EULER;
   else if( discretization_method_name == "rk2" )
      discretization_method = sdo::ButcherTableau::RUNGE_KUTTA_2;
   else if( discretization_method_name == "rk3" )
      discretization_method = sdo::ButcherTableau::RUNGE_KUTTA_3;
   else if( discretization_method_name == "rk4" )
      discretization_method = sdo::ButcherTableau::RUNGE_KUTTA_4;
   else if( discretization_method_name == "imid2" )
      discretization_method = sdo::ButcherTableau::IMPLICIT_MIDPOINT_2;
   else if( discretization_method_name == "igl4" )
      discretization_method = sdo::ButcherTableau::GAUSS_LEGENDRE_4;
   else
   {
      std::cerr << "error: unknow discretization method '" << discretization_method_name << "'\n";
      exit( 0 );
   }

   std::vector<std::string> mdlFiles;
   std::vector<std::string> vocFiles;
   std::vector<std::string> vpdFiles;

   for( std::vector<std::string>::iterator iter = input_files.begin(); iter != input_files.end(); ++iter )
   {
      if( boost::algorithm::ends_with( *iter, ".mdl" ) )
         mdlFiles.push_back( *iter );
      else if( boost::algorithm::ends_with( *iter, ".voc" ) )
         vocFiles.push_back( *iter );
      else if( boost::algorithm::ends_with( *iter, ".vop" ) || boost::algorithm::ends_with( *iter, ".sdo" ) )
      {
         try
         {

            sdo::VopFile vopfile = sdo::parse_vop_file(*iter);

            if( !vopfile.getModelFile().empty() )
               mdlFiles.push_back( vopfile.getModelFile() );

            if( !vopfile.getControlFile().empty() )
               vocFiles.push_back( vopfile.getControlFile() );

            if( !vopfile.getObjectiveFile().empty() )
               vpdFiles.push_back( vopfile.getObjectiveFile() );
         }
         catch( const sdo::parse_error &err )
         {
            std::cerr << err.what();
         }
         catch( const std::ifstream::failure &err ) {
            std::cerr << "Error: cannot read file '" << *iter << "': " << err.what();
         }
      }
      else if( boost::algorithm::ends_with( *iter, ".vpd" ) )
      {
         vpdFiles.push_back( *iter );
      }
      else
      {
         std::cerr << "Error: unknown file type '" << *iter << "'\n";
         exit( 0 );
      }
   }

   std::string objectiveFile;

   if(vpdFiles.size() == 1)
   {
      objectiveFile = vpdFiles.front();
   }
   else if( vpdFiles.size() > 1 )
   {
      std::cerr << "Found multiple objective functions:\n";

      for( std::size_t i = 0; i < vpdFiles.size(); ++i )
         std::cerr << "[ " << i << " ]: " << vpdFiles[i] << "\n";

      while( true )
      {
         std::cerr << "Choose: ";
         unsigned choice;
         std::cin >> choice;

         if( std::cin.good() && choice < vpdFiles.size() )
         {
            objectiveFile = vpdFiles[choice];
            break;
         }
      }
   }

   try
   {
      sdo::ExpressionGraph exprGraph;
      exprGraph.useUniqueConstants(true);

      for(auto &vocFile : vocFiles)
      {
         try
         {
            sdo::parse_voc_file(vocFile, exprGraph);
         }
         catch( const std::ifstream::failure &err )
         {
            std::cerr << "Error: cannot read file '" << vocFile << "'\n";
         }
      }
      for(auto &mdlFile: mdlFiles)
      {
         try
         {
            sdo::parse_mdl_file(mdlFile, exprGraph);
         }
         catch( const std::ifstream::failure &err )
         {
            std::cerr << "Error: cannot read file '" << mdlFile << "'\n";
         }
      }

      exprGraph.analyze();

      gams::GamsGenerator gams(exprGraph, discretization_method);

      if(lookup_type == "sos2")
         gams.setLookupFormulationTypes(gams::LookupFormulationType::SOS2);
      else if(lookup_type == "interactive") {
         for(auto& entry : exprGraph.getSymbolTable()) {
            sdo::LookupTable* lkpTable;
            if(entry.second->op == sdo::ExpressionGraph::LOOKUP_TABLE) {
               lkpTable = entry.second->lookup_table;
            } else if (entry.second->op == sdo::ExpressionGraph::APPLY_LOOKUP) {
               auto range = exprGraph.getSymbol(entry.second->child1);
               if(!range.empty())
                  continue;
               lkpTable = entry.second->child1->lookup_table;
            } else {
               continue;
            }

            int type = -1;
            std::cout << "Found Lookup '" << entry.first << "' used at: \n";
            for( auto &usage : entry.second->usages )
               std::cout <<  "\t" << usage << "\n";
            do {
               std::cout << "Choose type [0=SPLINE, 1=SOS2]: ";
               std::cin >> type;
            } while(type != 0 && type != 1);

            if(type == 0) {
               gams.setLookupFormulationType(lkpTable, gams::LookupData{ entry.first, gams::LookupFormulationType::SPLINE});
            } else {
               gams.setLookupFormulationType(lkpTable, gams::LookupData{ entry.first, gams::LookupFormulationType::SOS2});
            }
         }
      }

      if( !objectiveFile.empty() )
      {
         sdo::Objective objective;
         sdo::parse_vpd_file(objectiveFile, objective);
         gams.addObjective(std::move(objective));
      }
      else
      {
         gams.addArbitraryObjective();
      }
      gams.emitGams( out );
   }
   catch( const sdo::parse_error &err )
   {
      std::cerr << err.what();
   }
   catch( const std::ifstream::failure &err )
   {
      std::cerr << "Error: cannot read file\n";
   }

   return 0;
}
