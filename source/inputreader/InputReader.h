/*-----------------------------------------------------------------------*\
|																		  |
|	_   _           _            _____ __  __  ____  _  ________  		  |
|	| \ | |         | |          / ____|  \/  |/ __ \| |/ /  ____| 		  |
|	|  \| | ___   __| |_   _ ___| (___ | \  / | |  | | ' /| |__    		  |
|	| . ` |/ _ \ / _` | | | / __|\___ \| |\/| | |  | |  < |  __|   		  |
|	| |\  | (_) | (_| | |_| \__ \____) | |  | | |__| | . \| |____  		  |
|	|_| \_|\___/ \__,_|\__,_|___/_____/|_|  |_|\____/|_|\_\______|		  |                                                              |
|                                                                         |
|   Author: Matteo Mensi <matteo.mensi@mail.polimi.it>                    |
|   CRECK Modeling Group <http://creckmodeling.chem.polimi.it>            |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
\*-----------------------------------------------------------------------*/



#ifndef NODUSSMOKE_INPUTREADER_H
#define NODUSSMOKE_INPUTREADER_H

// Boost filesystem handler library
#include "boost/filesystem.hpp"

// OpenSMOKE++ Definitions
#include "OpenSMOKEpp"

// CHEMKIN maps
#include "maps/Maps_CHEMKIN"

#include "NodusSMOKE"

// LEGACY INCLUDES - SILVIO TRESPI IS RESPONSIBLE
#include <sstream>

#include <fstream>

#include <iomanip>

#include <cctype>

#include <cmath>

#include <cstdlib>

#include "silviosource/Classes.h"

#include "silviosource/Functions.h"

namespace NodusSMOKE
{

  class InputReader
  {

    //!  A class to provide an easier interface for Silvio Trespi's input reader for reactor networks
    /*!
    The purpose of this class is to provide understandable functions to read a reactor network morphology and kinetic scheme.
    The author of this class is Matteo Mensi but the inner files are work of Silvio Trespi's.
    Maintenability of the back-end is not assured by Matteo.
    Please ask Silvio Trespi <silvio.trespi@mail.polimi.it> or look at his own documentation.
    */

    public:

        /**
        *@brief Default constructor
        *@param none
        */
		InputReader(boost::filesystem::path Working_Folder_Path,
					OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMapGasIN,
					OpenSMOKE::KineticsMap_CHEMKIN& kineticsMapGasIN,
					OpenSMOKE::ThermodynamicsMap_Solid_CHEMKIN& thermodynamicsMapSolidIN,
					OpenSMOKE::KineticsMap_Solid_CHEMKIN& kineticsMapSolidIN);
			

        /**
        *@brief Default destructor
        *@param none
        */
        ~InputReader();

        /**
        *@brief Creates the thermodynamics and kinetics maps
        */
        void CreateInfoFromMaps();

        /**
        *@brief Sets the type of simulation
        */
		void CheckSimulationType(bool &solid_case_, double &solid_density_main, bool &recycle_boolean_, bool &video_print_, bool &legacy_print_, bool &restart, bool &rtd, std::string &track);

        /**
        *@brief Read the input, fills the storage structures to be copied by import functions and draws the network - SILVIO TRESPI IS RESPONSIBLE FOR WHAT IS INSIDE
        *@param
        */
        int ReadAndDraw();

        /**
        *@brief Copies necessary data over the main file where they will be provided to the ReactorNetwork class object
        *@param
        */
        void ImportData(std::vector<NodusSMOKE::StreamInfo>& StreamStruct, std::vector<NodusSMOKE::UnitInfo> &UnitStruct, std::vector<int> &From, std::vector<int> &To, std::vector<double> &relative_split );

    protected:
        /**
        Functions used by ImportData to read info of streams and units
        */
		void CheckCompositions();
        void ImportStreams(std::vector<NodusSMOKE::StreamInfo>& StreamStruct);
        void ImportUnits(std::vector<NodusSMOKE::UnitInfo>& UnitStruct, std::vector<NodusSMOKE::StreamInfo> &StreamStruct);
        void ImportFlowVectors(std::vector<NodusSMOKE::StreamInfo>& StreamStruct, std::vector<NodusSMOKE::UnitInfo>& UnitStruct,std::vector<int> &MainFrom, std::vector<int> &MainTo, std::vector<double> &Main_relative_split );

    private:

        // Important files
        boost::filesystem::path kinetic_path;
        std::string input_path;
        std::string network_map_path;
        std::string track_this;
        
        // Simulation type boolean
        bool solid_case_input;
        bool NoRecycles;
		bool video_print_input;
		bool legacy_print_input;
        bool rtd_;
        bool restart_;

        // Simulation pressure
        double System_Pressure;
        std::string System_Pressure_unit_of_measurement;

		// Solid density
		double rho_solid;
		std::string rho_solid_unit_of_measurement;

        // Indexes
        int reactor_count;          //Unit and stream counters
        int mixer_count;
        int splitter_count;
        int phasesplitter_count;
        int stream_count;
        int n_species;
        int n_gas_species;
        int n_solid_species;

        // local gas maps
        OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMapGas;
  	    OpenSMOKE::KineticsMap_CHEMKIN& kineticsMapGas;   

        // local solid maps
        OpenSMOKE::ThermodynamicsMap_Solid_CHEMKIN& thermodynamicsMapSolid;
  	    OpenSMOKE::KineticsMap_Solid_CHEMKIN& kineticsMapSolid;            

        // Stream storage
        std::vector <Stream_class> Stream;

        // Mixer storage
        std::vector <Mixer_unit> Mixer;

        // Phase splitter storage
        std::vector <PS_unit> PhaseSplitter;
        
        // Splitter storage
        std::vector <S_unit> Splitter;

        // Reactor storage
        std::vector <Reactor_unit> Reactor;    
    
        // Incidence vectors storages
        std::vector <string> From;     //Unit FROM which the stream arrives
        std::vector <string> To;       //Unit TO which the stream arrives
        std::vector <double> Splitting_ratio; 
        std::vector <string> Stream_phase;
        std::vector <int> Stream_number;
  };


}

#include "InputReader.hpp"


#endif /*NODUSSMOKE_INPUTREADER_H*/
