
namespace NetSMOKE
{
    

  InputReader::InputReader(boost::filesystem::path Working_Folder_Path,
				OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMapGasIN,
				OpenSMOKE::KineticsMap_CHEMKIN& kineticsMapGasIN,
				OpenSMOKE::ThermodynamicsMap_Solid_CHEMKIN& thermodynamicsMapSolidIN,
				OpenSMOKE::KineticsMap_Solid_CHEMKIN& kineticsMapSolidIN)
				:
				thermodynamicsMapGas(thermodynamicsMapGasIN),
				kineticsMapGas(kineticsMapGasIN),
				thermodynamicsMapSolid(thermodynamicsMapSolidIN),
				kineticsMapSolid(kineticsMapSolidIN)
		  {
			// Folder containing the required xmls
			kinetic_path = Working_Folder_Path / "kinetics";
			// File containing the input
			boost::filesystem::path input_path_boost = Working_Folder_Path / "input.dic";
			input_path = input_path_boost.string();
			// Folder to print the network with graphviz
			boost::filesystem::path network_map_path_boost = Working_Folder_Path / "NetworkMap_graphvizDot.dot";
			network_map_path = network_map_path_boost.string();

			// Default options
			video_print_input = false; 
			legacy_print_input = false;
			NoRecycles = false;
			restart_ = false;
			rtd_ = false;

		 }

  
  InputReader::~InputReader()
  {
        
        //DESTRUCTOR since I have no idea how much memory will occupy this code
        //Please call this after having imported the data and the maps over your main code
        
		/*
        delete thermodynamicsMapGas;
        delete kineticsMapGas;
        
		if (solid_case_input == true) {
			delete thermodynamicsMapSolid;
			delete kineticsMapSolid;
		}
		*/
		

      cout<< "The input reader memory is being released. The associated object is now destroyed. " << endl; 
  }
  

  void InputReader::CheckSimulationType(bool &solid_case_, double &solid_density_main, bool &recycle_boolean_, bool &video_print_, bool &legacy_print_, bool &restart, bool &rtd, std::string &track)
  {

      solid_case_input = false;

      // Check if any input stream is solid or mix
      for (unsigned int i = 0; i < Stream.size(); ++i){
          for (unsigned int j = 0; j < Stream_number.size(); ++j){
            if (Stream_number[j] == Stream[i].Name){
                if (From[j] == "0"){ // if it is an inlet
                    if (Stream[i].Phase == "Solid" || Stream[i].Phase == "Mix"){
                        solid_case_input = true;
						solid_density_main = rho_solid;

                    }
                }
            }
          }
      }

	  solid_case_ = solid_case_input; // Transmit info to main
	  recycle_boolean_ = NoRecycles;
	  video_print_ = video_print_input;
	  legacy_print_ = legacy_print_input;
	  restart = restart_;
	  rtd=rtd_;
	  if(rtd_ == true){
		  track = track_this;
	  }

  }


 void InputReader::CreateInfoFromMaps()
  {

    if (solid_case_input == true){
        n_species = thermodynamicsMapSolid.NumberOfSpecies();
		n_gas_species = thermodynamicsMapSolid.number_of_gas_species();
		n_solid_species = thermodynamicsMapSolid.number_of_solid_species();
  	}
    else {
        n_species = thermodynamicsMapGas.NumberOfSpecies();
		n_gas_species = n_species;
		n_solid_species = 0;
	}

    std::cout << " DONE" << std::endl;

 };

  void InputReader::ImportData(vector<NetSMOKE::StreamInfo> &StreamStruct, vector<NetSMOKE::UnitInfo> &UnitStruct, std::vector<int> &From, std::vector<int> &To, std::vector<double> &relative_split )
  {


    if (System_Pressure_unit_of_measurement == "Pa")
        System_Pressure = System_Pressure;
    else if (System_Pressure_unit_of_measurement == "atm")
        System_Pressure = System_Pressure*101325.;
    else if (System_Pressure_unit_of_measurement == "bar")
        System_Pressure = System_Pressure*10.e5;
    else
    {
        std::cout << "ERROR: Pressure is absent, has wrong units or are not supported" << std::endl;
		OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
    }

    CheckCompositions(); // Need to fix that shit he did

    ImportStreams(StreamStruct);
    ImportUnits(UnitStruct, StreamStruct);
    ImportFlowVectors(StreamStruct, UnitStruct, From, To, relative_split);
    std::cout << " DATA SUCCESFULLY IMPORTED" << std::endl;

  }

  void InputReader::ImportStreams(vector<NetSMOKE::StreamInfo> &StreamStruct)
  {
        // Copy stream info
        /* It looks like shitty code but silvio just didn't use any kind of imposed order to save the data 
        so you have to constantly look for matching info and reorder your structures accordingly */

        std::cout << " Reading streams' info..." << std::endl;
        StreamStruct.resize(Stream_number.size());

        // Initialize all streams
        for (int j = 0; j < Stream_number.size(); ++j) {
			StreamStruct[j].name = Stream_number[j];
            StreamStruct[j].mass_flow_rate_gas = 0.;
			StreamStruct[j].mass_flow_rate_solid = 0.;
            StreamStruct[j].temperature = 1000.;
            StreamStruct[j].pressure = System_Pressure;
            OpenSMOKE::ChangeDimensions(n_gas_species, &StreamStruct[j].omega_gas, true);
			OpenSMOKE::ChangeDimensions(n_gas_species, &StreamStruct[j].x_gas, true);
			OpenSMOKE::ChangeDimensions(n_solid_species, &StreamStruct[j].omega_solid, true);
			OpenSMOKE::ChangeDimensions(n_solid_species, &StreamStruct[j].x_solid, true);
	    }
		

        // Copy info of declared streams
        for (unsigned int i = 0; i<Stream.size(); ++i){
            for (unsigned int j = 0; j <Stream_number.size(); ++j){
                if (Stream_number[j] == Stream[i].Name){
					// Phase import
					StreamStruct[j].phase = Stream[i].Phase;
					// Check temperature units
                    if (Stream[i].Temperature_unit_of_measurement == "K") { /* Do nothing */ StreamStruct[j].temperature = Stream[i].Temperature;}
                    else if (Stream[i].Temperature_unit_of_measurement == "°C") { /*convert to K*/ StreamStruct[j].temperature = Stream[i].Temperature + 273.15;}
                    else
                    {
                    std::cout << "ERROR: Stream " << Stream[i].Name << " T is absent, has wrong units or are not supported" << std::endl;
					OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
                    }

                    // Check mass flow units
					double tempvar_MassFlow;
                    if (Stream[i].MassFlowRate_unit_of_measurement == "kg/s")
                        tempvar_MassFlow = Stream[i].MassFlowRate;
                    else if (Stream[i].MassFlowRate_unit_of_measurement == "g/s")
                        tempvar_MassFlow = Stream[i].MassFlowRate/1000;
                    else if (Stream[i].MassFlowRate_unit_of_measurement == "kg/hr")
                        tempvar_MassFlow = Stream[i].MassFlowRate/3600;
                    else if (Stream[i].MassFlowRate_unit_of_measurement == "kg/min")
                        tempvar_MassFlow = Stream[i].MassFlowRate*60;
                    else
                    {
                    std::cout << "ERROR: Stream " << Stream[i].Name << " MassFlowRate is absent, has wrong units or are not supported" << std::endl;
					OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
                    }
                    
					// Import composition and mass flow, according to the declared phase
					if (Stream[i].Phase == "Gas"){
						for (unsigned int k = 1; k <= n_gas_species; ++k) {
							StreamStruct[j].omega_gas[k] = Stream[i].omega_teo[k - 1];
						}
						StreamStruct[j].omega_gas = StreamStruct[j].omega_gas.operator/=(StreamStruct[j].omega_gas.SumElements());
						StreamStruct[j].mass_flow_rate_gas = tempvar_MassFlow;
					}
					if (Stream[i].Phase == "Solid"){
						for (unsigned int k = 1; k <= n_solid_species; ++k) {
							StreamStruct[j].omega_solid[k] = Stream[i].omega_teo[n_gas_species + k - 1];
						}
						StreamStruct[j].omega_solid = StreamStruct[j].omega_solid.operator/=(StreamStruct[j].omega_solid.SumElements());
						StreamStruct[j].mass_flow_rate_solid = tempvar_MassFlow;
					}
					if (Stream[i].Phase == "Mix"){
						for (unsigned int k = 1; k <= n_gas_species; ++k) {
							StreamStruct[j].omega_gas[k] = Stream[i].omega_teo[k - 1];
						}
						for (unsigned int k = 1; k <= n_solid_species; ++k) {
							StreamStruct[j].omega_solid[k] = Stream[i].omega_teo[n_gas_species + k - 1];
						}
						StreamStruct[j].mass_flow_rate_gas = tempvar_MassFlow*StreamStruct[j].omega_gas.SumElements();
						if (StreamStruct[j].mass_flow_rate_gas >0.)				
							StreamStruct[j].omega_gas = StreamStruct[j].omega_gas.operator/=(StreamStruct[j].omega_gas.SumElements());
						StreamStruct[j].mass_flow_rate_solid = tempvar_MassFlow*StreamStruct[j].omega_solid.SumElements();
						if (StreamStruct[j].mass_flow_rate_solid >0.)	
							StreamStruct[j].omega_solid = StreamStruct[j].omega_solid.operator/=(StreamStruct[j].omega_solid.SumElements());
					} 
				}	//  End of importing stream[i]
            } // End of search for corresponding order 
        }	// End of for cycle over all streams

		
        // Make sure no stream is mix/solid if no inlet is solid
        if (solid_case_input == false){
            for (unsigned int i=0; i<StreamStruct.size(); ++i)
                StreamStruct[i].phase = "Gas";
		}
		

        std::cout << " DONE" << std::endl;
        
  }


  void InputReader::ImportUnits(vector<NetSMOKE::UnitInfo> &UnitStruct, vector<NetSMOKE::StreamInfo> &StreamStruct)
  {

        // Import Reactors
        std::cout << " Reading reactors..." << std::endl;
        for (int i = 0; i<Reactor.size(); ++i){
            NetSMOKE::UnitInfo TempUnit;
            TempUnit.tag = "Reactor";            
            TempUnit.name = Reactor[i].Name;
            TempUnit.phase = Reactor[i].Phase;
			TempUnit.type = Reactor[i].Type;
			TempUnit.UA = 0.;
			if (solid_case_input == false){
				TempUnit.Nequations = n_gas_species + 1;
			}
			else if (solid_case_input == true){
				if (Reactor[i].Phase == "Mix" || Reactor[i].Phase == "Solid"  )
					TempUnit.Nequations = n_species + 3;
				else if (Reactor[i].Phase == "Gas")
					TempUnit.Nequations = n_gas_species + 2;
			}


			if (Reactor[i].Energy == "Isothermal"){
				TempUnit.energy = Reactor[i].Energy;
			}
			else {
				TempUnit.energy = "NonIsothermal";
			}
            TempUnit.inlets_names.push_back(Reactor[i].In_stream);
            TempUnit.outlets_names.push_back(Reactor[i].Out_stream);

            // Check residence time
            if (Reactor[i].Residence_time_unit_of_measurement == "s")
            TempUnit.residence_time = Reactor[i].Residence_time;
            else if (Reactor[i].Residence_time_unit_of_measurement == "ms")
            TempUnit.residence_time = Reactor[i].Residence_time/1000.;
            else if (Reactor[i].Residence_time_unit_of_measurement == "min")
            TempUnit.residence_time = Reactor[i].Residence_time*60.;
            else if (Reactor[i].Residence_time_unit_of_measurement == "hr")
            TempUnit.residence_time = Reactor[i].Residence_time*3600.;
            else
            {
				std::cout << " the unit is " << Reactor[i].Residence_time_unit_of_measurement << std::endl;
                std::cout << "ERROR: Reactor " << Reactor[i].Name << " residence time is absent, has wrong units or are not supported" << std::endl;
				OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
            }

            // Check T if Isothermal or UA if not isothermal
            if (Reactor[i].Energy == "HeatExchanger"){
                if (Reactor[i].UA_unit_of_measurement == "W/K")
                    TempUnit.UA = Reactor[i].UA;
                else
                {
                    std::cout << "ERROR: Reactor " << Reactor[i].Name << " U*A is absent, has wrong units or are not supported" << std::endl;
					OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
                }
            }
            else if (Reactor[i].Energy == "Isothermal") {
                if (Reactor[i].Isothermal_temperature_unit_of_measurement == "K"){
                    TempUnit.temperature = Reactor[i].Isothermal_temperature;
                }
                else if (Reactor[i].Isothermal_temperature_unit_of_measurement == "°C"){
                    TempUnit.temperature = Reactor[i].Isothermal_temperature +273.15;
                }
                else
                {
                    std::cout << "ERROR: Reactor " << Reactor[i].Name << " T is absent, has wrong units or are not supported" << std::endl;
					OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
                }
            }           
            UnitStruct.push_back(TempUnit);
        }

        std::cout << " DONE" << std::endl;

        // Import Mixers
        std::cout << " Reading mixers..." << std::endl;
        for (int i = 0; i<Mixer.size(); ++i){
            NetSMOKE::UnitInfo TempUnit;
			TempUnit.tag = "Mixer";
            TempUnit.name = Mixer[i].Name;
            TempUnit.energy = Mixer[i].Energy;
			TempUnit.Nequations = 0;

            for (unsigned int j = 0; j<Mixer[i].In_stream.size(); ++j)
                TempUnit.inlets_names.push_back(Mixer[i].In_stream[j]);

            TempUnit.outlets_names.push_back(Mixer[i].Out_stream);
            
            // Check T if Isothermal
            if (Mixer[i].Energy == "Isothermal"){
                if (Mixer[i].Temperature_unit_of_measurement == "K"){
                    TempUnit.temperature = Mixer[i].Temperature;
                }
                else if (Mixer[i].Temperature_unit_of_measurement == "°C"){
                    TempUnit.temperature = Mixer[i].Temperature+273.15;
                }
                else
                {
                    std::cout << "ERROR: Mixer " << Mixer[i].Name << " T is absent, has wrong units or are not supported" << std::endl;
					OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
                }
            }           
            UnitStruct.push_back(TempUnit);
        }
        std::cout << " DONE" << std::endl;

        //Import Splitters
        std::cout << " Reading splitters..." << std::endl;
        for (int i = 0; i<Splitter.size(); ++i){
            NetSMOKE::UnitInfo TempUnit;
			TempUnit.tag = "Splitter";
            TempUnit.name = Splitter[i].Name;
			TempUnit.Nequations = 0;

            for (unsigned int j = 0; j<Splitter[i].Out_stream.size(); ++j){
                TempUnit.outlets_names.push_back(Splitter[i].Out_stream[j]);
                TempUnit.split_factor.push_back(Splitter[i].Splitting_ratio[j]);
            }
        
            TempUnit.inlets_names.push_back(Splitter[i].In_stream);
          
            UnitStruct.push_back(TempUnit);
        }
        std::cout << " DONE" << std::endl;

        //Import Phase Splitters
        std::cout << " Reading phase splitters..." << std::endl;
        for (int i = 0; i<PhaseSplitter.size(); ++i){
            NetSMOKE::UnitInfo TempUnit;
			TempUnit.tag = "PhaseSplitter";
            TempUnit.name = PhaseSplitter[i].Name;
            TempUnit.inlets_names.push_back(PhaseSplitter[i].In_stream);
			TempUnit.Nequations = 0;

            for (unsigned int j = 0; j<PhaseSplitter[i].Out_stream.size(); ++j){
                TempUnit.outlets_names.push_back(PhaseSplitter[i].Out_stream[j]);
				TempUnit.splitted_phase.push_back(PhaseSplitter[i].Splitted_phase[j]);
            }
          
            UnitStruct.push_back(TempUnit);
		}
	
        std::cout << " DONE" << std::endl;

        // Set pressures
        for (unsigned int i = 0; i<UnitStruct.size(); ++i)
            UnitStruct[i].pressure = System_Pressure;

        std::cout << "Parsed " << UnitStruct.size() << " units!" << std::endl;

        // Make sure the silvio method of using "default stuff" does not fuck up
        if (solid_case_input == false){
            for (unsigned int i=0; i<UnitStruct.size(); ++i)
                UnitStruct[i].phase == "Gas";
        }

		// Assign "unit index" parameter
		for (unsigned int i = 0; i < UnitStruct.size(); ++i)
			UnitStruct[i].flow_matrix_index = i;
		
		// Assigning correct inlet and outlets stream
		for (int j = 0; j < UnitStruct.size(); ++j){
			for (int i = 0; i < UnitStruct[j].outlets_names.size(); ++i){
				for (int k = 0; k < StreamStruct.size(); ++k){
					if (UnitStruct[j].outlets_names[i] == StreamStruct[k].name){
						UnitStruct[j].outlets.push_back(k);
					}
				}
			}
			for (int i = 0; i < UnitStruct[j].inlets_names.size(); ++i){
				for (int k = 0; k < StreamStruct.size(); ++k){
					if (UnitStruct[j].inlets_names[i] == StreamStruct[k].name){
						UnitStruct[j].inlets.push_back(k);
					}
				}
			}
		}
   
    }


  void InputReader::ImportFlowVectors(std::vector<NetSMOKE::StreamInfo>& StreamStruct, std::vector<NetSMOKE::UnitInfo>& UnitStruct, std::vector<int> &MainFrom, std::vector<int> &MainTo, std::vector<double> &Main_relative_split)
  {
    std::cout << " Reading Flow Vectors..." << std::endl;
    
    MainFrom.resize(Stream_number.size());
    MainTo.resize(Stream_number.size());
    Main_relative_split.resize(Stream_number.size());

	// Create temporary unit name. I fucking hate this method how can you come up with so much disorganized
	// data even if you are a bachelor only
	std::vector<string> UnitName(UnitStruct.size() + 1);
	UnitName[0] = "0";
	for (unsigned int i = 0; i < UnitStruct.size(); ++i) {
		UnitName[i + 1] = UnitStruct[i].name;
	}

    Main_relative_split = Splitting_ratio;

    for (unsigned int i = 0; i < From.size(); ++i){
        for (unsigned int j = 0; j < UnitStruct.size()+1; ++j){
            if (UnitName[j] == From[i])
                MainFrom[i]=j;
        }
    }

    for (unsigned int i = 0; i < To.size(); ++i){
        for (unsigned int j = 0; j < UnitStruct.size()+1; ++j){
            if (UnitName[j] == To[i])
                MainTo[i]=j;
        }
    }

    std::cout << " DONE" << std::endl;
    
  }

  void InputReader::CheckCompositions()
  {
      for (unsigned int i = 0; i<Stream.size(); ++i){
          Stream[i].omega_teo.resize(n_species);
      }

    
      for (unsigned int i = 0; i<Stream.size(); ++i){
          if (Stream[i].Phase == "Solid"){
              // NEED PAOLO FUNCTION
              std::cout << "ERROR: Currently no pure solid inlet stream is permitted (Stream "<<Stream[i].Name <<") " << std::endl;
			  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
          }
          else
		  {
            for (unsigned int j = 0; j < Stream[i].Compound.size(); ++j)
			{
			
				if (solid_case_input == false) {
					Stream[i].omega_teo[thermodynamicsMapGas.IndexOfSpecies(Stream[i].Compound[j])-1] = Stream[i].Composition[j];
				}
				else {
					Stream[i].omega_teo[thermodynamicsMapSolid.IndexOfSpecies(Stream[i].Compound[j])-1] = Stream[i].Composition[j];
				}
            }
          }
      }
  }

  int InputReader::ReadAndDraw() {

	  using namespace std;

	  vector <string> Reactor_type_dictionary;            //List of dictionaries used
	  Reactor_type_dictionary.push_back("PSR");
	  Reactor_type_dictionary.push_back("PFR");

	  vector <string> Reactor_phase_dictionary;
	  Reactor_phase_dictionary.push_back("Gas");
	  Reactor_phase_dictionary.push_back("Solid");
	  Reactor_phase_dictionary.push_back("Mix");

	  vector <string> Reactor_energy_dictionary;
	  Reactor_energy_dictionary.push_back("Adiabatic");
	  Reactor_energy_dictionary.push_back("Isothermal");
	  Reactor_energy_dictionary.push_back("HeatExchanger");

	  vector <string> Solid_compound_dictionary;
	  Solid_compound_dictionary.push_back("C");
	  Solid_compound_dictionary.push_back("H");
	  Solid_compound_dictionary.push_back("O");
	  Solid_compound_dictionary.push_back("N");
	  Solid_compound_dictionary.push_back("S");
	  Solid_compound_dictionary.push_back("H2O");
	  Solid_compound_dictionary.push_back("Ash");

	  vector <string> Stream_solid_type_dictionary;
	  Stream_solid_type_dictionary.push_back("Biomass");
	  Stream_solid_type_dictionary.push_back("Cellulose");
	  Stream_solid_type_dictionary.push_back("Hemicellulose");
	  Stream_solid_type_dictionary.push_back("Lignin");
	  Stream_solid_type_dictionary.push_back("Char");
	  Stream_solid_type_dictionary.push_back("Coal");
	  Stream_solid_type_dictionary.push_back("Waste");

	  vector <string> Stream_phase_dictionary;
	  Stream_phase_dictionary.push_back("Gas");
	  Stream_phase_dictionary.push_back("Solid");
	  Stream_phase_dictionary.push_back("Liquid");
	  Stream_phase_dictionary.push_back("Mix");

	  vector <string> Splitted_phase_dictionary;      //Phase_splitter Outlet_streams cannot be Mix!
	  Splitted_phase_dictionary.push_back("Solid");
	  Splitted_phase_dictionary.push_back("Liquid");
	  Splitted_phase_dictionary.push_back("Gas");

	  vector <string> Mixer_energy_dictionary;        //Mixer energy can only be Isothermal or Adiabatic
	  Mixer_energy_dictionary.push_back("Adiabatic");
	  Mixer_energy_dictionary.push_back("Isothermal");


	  vector <string> Main_keyword;       //Main keywords

	  Main_keyword.push_back("@Reactor");
	  Main_keyword.push_back("@Mixer");
	  Main_keyword.push_back("@Splitter");
	  Main_keyword.push_back("@Phase-splitter");
	  Main_keyword.push_back("@Stream");
	  Main_keyword.push_back("@GeneralOptions"); // Added by Mensi
	
	  vector <string> GeneralOptions_keyword; // Added by Mensi
	  GeneralOptions_keyword.push_back("SystemPressure");
	  GeneralOptions_keyword.push_back("NoRecycles");
	  GeneralOptions_keyword.push_back("SetOutletsVideoPrint");
	  GeneralOptions_keyword.push_back("SetLegacyOSppPrint");
	  GeneralOptions_keyword.push_back("SolidDensity");
	  GeneralOptions_keyword.push_back("Restart");
	  GeneralOptions_keyword.push_back("RTD");


	  vector <string> Reactor_keyword;        //List of units and stream keywords 

	  Reactor_keyword.push_back("Type");
	  Reactor_keyword.push_back("Energy");
	  Reactor_keyword.push_back("ResidenceTime");
	  Reactor_keyword.push_back("Phase");
	  Reactor_keyword.push_back("Temperature");

	  Reactor_keyword.push_back("UA");
	  Reactor_keyword.push_back("Inlet_stream");
	  Reactor_keyword.push_back("Outlet_stream");


	  vector <string> Stream_keyword;

	  Stream_keyword.push_back("Phase");
	  Stream_keyword.push_back("SolidType");
	  Stream_keyword.push_back("MassFlowRate");
	  Stream_keyword.push_back("Temperature");
	  Stream_keyword.push_back("MassFraction");

	  vector <string> Mixer_keyword;

	  Mixer_keyword.push_back("Energy");
	  Mixer_keyword.push_back("Inlet_stream");
	  Mixer_keyword.push_back("Outlet_stream");
	  Mixer_keyword.push_back("Temperature");


	  vector <string> Splitter_keyword;

	  Splitter_keyword.push_back("Inlet_stream");
	  Splitter_keyword.push_back("Outlet_stream");

	  vector <string> PhaseSplitter_keyword;

	  PhaseSplitter_keyword.push_back("Inlet_stream");
	  PhaseSplitter_keyword.push_back("Outlet_stream");



	  ifstream infile;            //Input_file
	  infile.open(input_path);



	  if (infile.fail()) //Checking file opening
	  {
		  cout << "Error while opening file\n";
		  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
	  }
	  cout << "Reading from the input file \n\n";


	  Reactor.resize(50);         //Allocating some memory

	  Mixer.resize(50);

	  Splitter.resize(50);

	  PhaseSplitter.resize(50);

	  Stream.resize(150);




	  string line;

	  int line_number = 0;    //number of the current line from the input file        

	  reactor_count = 0;          //Unit and stream counters
	  mixer_count = 0;
	  splitter_count = 0;
	  phasesplitter_count = 0;
	  stream_count = 0;
	  System_Pressure = 0;
	  System_Pressure_unit_of_measurement = "-";






	  vector <string> Splitted_line;

	  vector <string> Unit_name; //This vector will be filled with the name of every unit(Reactor,mixer,splitter ad phase-splitter)

	  Unit_name.push_back("0");  // Adding Input_output fictitious reactor

	  vector <int> Global_stream_number; //This vector will be filled with every input and output streams(there will be duplicates)

	  vector <string> Species;// This vector is filled with the species found


	  ofstream drawing;           //Reactor network digraph instructions 

	  drawing.open(network_map_path, ios::out);

	  if (drawing.fail()) //Checking drawing opening
	  {
		  cout << "Error while opening NetworkMap_graphvizDot!\n";
		  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
	  }

	  drawing << "digraph ReactorNetwork { " << endl << "size = \"8,11\";" << endl << "edge [fontsize = 400, style = \"setlinewidth(11)\",arrowsize = 20];" << endl;

	  drawing << "node [fontsize = 400];" << endl;
	  drawing << "ratio = fill;" << endl;
	  drawing << "ranksep = 12;" << endl;
	  drawing << "nodesep = 1;" << endl;

	  drawing << "INPUT[shape = polygon,sides = 5,peripheries = 3,color = lightblue,style = filled,height = 32,width = 32,rank = source];" << endl;
	  drawing << "OUTPUT[shape = polygon,sides = 5,peripheries = 3,color = lightblue,style = filled,height = 32,width = 32,rank = sink];" << endl;

	  NoRecycles = false;
	  while (getline(infile, line))  //Reading file

	  {

		  line_number++;



		  Splitted_line = String_comma_separator(line);  //Separating line elements using spaces and commas

		  if (Splitted_line.size() != 0)  //Checking that the line is not empty

		  {
				// This IF triggered by general options was added by Mensi. Quite the hack I know.
			  if (Splitted_line.at(0) == Main_keyword.at(5)) // GETTING SYSTEM PRESSURE AND METHOD
			  {


				  while (line.find("//") == string::npos) {

					  Splitted_line = String_comma_separator(line);

					  getline(infile, line);
					  line_number++;
					  if (Splitted_line.size() != 0) {
						  if (Splitted_line.at(0) == GeneralOptions_keyword[0])		// System pressure
						  {
							  double number;
							  stringstream(Splitted_line.at(1)) >> number;
							  System_Pressure = number;
							  System_Pressure_unit_of_measurement = Splitted_line.at(2);
						  }
						  else if (Splitted_line.at(0) == GeneralOptions_keyword[1]) // If NoRecyle options
						  {
							  NoRecycles = true;
						  }
						  else if (Splitted_line.at(0) == GeneralOptions_keyword[2]) // If VideoPrint options
						  {
							  video_print_input = true;
						  }
						  else if (Splitted_line.at(0) == GeneralOptions_keyword[3]) // If LegacyPrint options
						  {
							  legacy_print_input = true;
						  }
						  else if (Splitted_line.at(0) == GeneralOptions_keyword[4]) // Solid density
						  {
							  double number;
							  stringstream(Splitted_line.at(1)) >> number;
							  rho_solid_unit_of_measurement = Splitted_line.at(2);
							  if (rho_solid_unit_of_measurement == "kg/m3") {
								  rho_solid = number;
							  }
							  else {
								  OpenSMOKE::FatalErrorMessage("Please express the solid density in kg/m3!");
							  }
						  }
							else if (Splitted_line.at(0) == GeneralOptions_keyword[5]) // If restart 
						  {
							  restart_ = true;
						  }
						else if (Splitted_line.at(0) == GeneralOptions_keyword[6]) // If RTD 
						  {
							  rtd_ = true;
							  stringstream(Splitted_line.at(1)) >> track_this;

						  }
					  }
				  }
			  }


			  else if (Splitted_line.at(0) == Main_keyword.at(0))       //getting REACTOR information          

			  {
				  vector <bool> R_keyword_controlling(Reactor_keyword.size(), false); // if a reactor keyword is found, then the related boolean in this vector turns true


				  if (LengthControl(Splitted_line, 2, line_number) == 1) //CHECK: if reactor line does not have two elements (Reactor + Name), the program finishes
				  {
					  cout << "This line should have two elements! ";
					  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
				  }



				  if (FirstLetterThenNumberControl(Splitted_line.at(1), 'R') == "No") //CHECK: if reactor name does not start with 'R', the program finishes
				  {
					  cout << "Look at line " << line_number << ".";
					  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
				  }

				  (Reactor[reactor_count]).Name = Splitted_line.at(1);

				  Unit_name.push_back(Splitted_line.at(1));


				  for (int i = 0; i < reactor_count; i++) //CHECK: if two reactors have the same name, the program finishes
				  {
					  if (Reactor[reactor_count].Name == Reactor[i].Name)
					  {
						  cout << "ERROR at line " << line_number << "!\nTwo reactors are identified by the same name " << Reactor[reactor_count].Name << "!";
						  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
					  }
				  }

				  while (line.find("//") == string::npos)  //This while loop gets lines until the separation element( // ) is found

				  {
					  getline(infile, line);
					  line_number++;

					  vector <string> Splitted_line = String_comma_separator(line);


					  bool Unrecognized_keyword = true; // This boolean turns false if the first element of the line is one of the reactor keywords

					  if (Splitted_line.size() != 0) // Checking the line is not empty

					  {


						  for (int j = 0; j < Reactor_keyword.size(); j++)

						  {


							  if (Splitted_line.at(0) == Reactor_keyword.at(j))

							  {

								  Unrecognized_keyword = false;

								  if (R_keyword_controlling[j] == true)  // If the boolean related to eactor_keyword[j] has been already found, then the program finishes

								  {
									  cout << "ERROR at line " << line_number << "!\n";
									  cout << "Reactor " << Reactor[reactor_count].Name << " has more than one keyword " << Reactor_keyword[j] << "!";
									  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
								  }

								  R_keyword_controlling[j] = true;



								  switch (j)

								  {
								  case 0: //Type
								  case 1: // Energy
								  case 3: // Phase
									  if (LengthControl(Splitted_line, 2, line_number) == 1)
									  {
										  cout << "This line should be \"" << Reactor_keyword[j] << "\" followed by a word.";
										  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
									  }


									  break;

								  case 6: // Inlet_stream
								  case 7: // Outlet_stream
									  if (LengthControl(Splitted_line, 2, line_number) == 1)
									  {
										  cout << "This line should be \"" << Reactor_keyword[j] << "\" (stream number)";
										  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
									  }

									  if (IntControl(Splitted_line.at(1)) == "No")
									  {
										  cout << "ERROR at line " << line_number << "!\n";
										  cout << "Stream number must be a positive integer!";
										  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
									  }



									  break;

								  case 2: // Residence time
								  case 4: // Temperature
								  case 5: //UA

									  if (LengthControl(Splitted_line, 3, line_number) == 1)
									  {
										  cout << "This line should be " << Reactor_keyword[j] << " (number) (unit of measurement)";
										  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
									  }

									  if (DoubleControl(Splitted_line.at(1)) == "No")
									  {
										  cout << "Value of " << Reactor_keyword[j] << " is wrong!";
										  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
									  }
									  if (LetterControl(Splitted_line.at(2)) == "No")
									  {
										  cout << "Unit of measurement of " << Reactor_keyword[j] << " is wrong!";
										  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
									  }


								  }

								  if (j == 0)    //Type

								  {


									  if ((DictionaryControl(Splitted_line.at(1), Reactor_type_dictionary)) != 0)
									  {
										  cout << "\nLine " << line_number << ".";
										  return 1;
									  }

									  Reactor[reactor_count].Type = Splitted_line.at(1);



								  }

								  else if (j == 1)//Energy

								  {

									  if ((DictionaryControl(Splitted_line.at(1), Reactor_energy_dictionary)) != 0)
									  {
										  cout << "\nLine " << line_number << ".";
										  return 1;
									  }

									  Reactor[reactor_count].Energy = Splitted_line.at(1);


								  }

								  else if (j == 2)//ResidenceTime

								  {
									  double number;
									  stringstream(Splitted_line.at(1)) >> number;
									  Reactor[reactor_count].Residence_time = number;
									  Reactor[reactor_count].Residence_time_unit_of_measurement = Splitted_line.at(2);

									  if (Reactor[reactor_count].Residence_time < 0)
									  {
										  cout << "ERROR at line " << line_number << "!Residence time of reactor " << Reactor[reactor_count].Name << " is negative!\n";
										  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
									  }
								  }

								  else if (j == 3)//Phase

								  {

									  if ((DictionaryControl(Splitted_line.at(1), Reactor_phase_dictionary)) != 0)
									  {
										  cout << "\nline " << line_number << ".";
										  return 1;
									  }

									  Reactor[reactor_count].Phase = Splitted_line.at(1);


								  }

								  else if (j == 4)//Temperature

								  {
									  double number;
									  stringstream(Splitted_line.at(1)) >> number;
									  Reactor[reactor_count].Isothermal_temperature = number;
									  Reactor[reactor_count].Isothermal_temperature_unit_of_measurement = Splitted_line.at(2);

									  if (Reactor[reactor_count].Isothermal_temperature < 0)
									  {
										  if (Reactor[reactor_count].Isothermal_temperature_unit_of_measurement == "K")
										  {
											  cout << "ERROR at line " << line_number << "!Temperature of isothermal reactor " << Reactor[reactor_count].Name << " is under 0 K!";
											  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
										  }
									  }

								  }



								  else if (j == 5)//UA

								  {
									  double number;
									  stringstream(Splitted_line.at(1)) >> number;
									  Reactor[reactor_count].UA = number;
									  Reactor[reactor_count].UA_unit_of_measurement = Splitted_line.at(2);



								  }

								  else if (j == 6)//Inlet_stream

								  {
									  int number;
									  stringstream(Splitted_line.at(1)) >> number;

									  Reactor[reactor_count].In_stream = number;
									  Global_stream_number.push_back(number);

								  }

								  else if (j == 7)//Outlet_stream

								  {
									  int number;
									  stringstream(Splitted_line.at(1)) >> number;

									  Reactor[reactor_count].Out_stream = number;
									  Global_stream_number.push_back(number);

								  }



							  }

						  } // end of the for cycle


						  if (Unrecognized_keyword == true && Splitted_line.at(0) != "//") // CHECK: if Unrecognized_keyword == true then The first element of a line 
						  {                                                                   // is not a Reactor_keyword and therefore the program finishes

							  cout << "ERROR at line " << line_number << "!\n";

							  for (int i = 0; i < Main_keyword.size(); i++)
							  {
								  if (Splitted_line.at(0) == Main_keyword.at(i))
								  {
									  cout << "This is a new " << Main_keyword.at(i) << " declaration but the last one is not finished yet!\n";
									  cout << "Remember to end unit or stream declarations with \'//\'.";
									  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
								  }
							  }

							  cout << "\"" << Splitted_line.at(0) << "\" is not a Reactor keyword!";
							  return 1;
						  }




						  if (End_of_file((infile)) != 0)
						  {

							  return 1;
						  }


					  }


				  }// end of while loop



				  for (int j = 0; j < R_keyword_controlling.size(); j++)  // Checking the presence of all reactor keywords needed
				  {

					  if (R_keyword_controlling[j] == false)

					  {
						  if (j < 4 || j >= 6) // compulsory keywords
						  {
							  cout << "ERROR at reactor " << Reactor[reactor_count].Name << "! " << Reactor_keyword[j] << " is missing!";
							  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
						  }
						  else
						  {
							  if (j == 4) //Temperature, only needed for isothermal reactors
							  {
								  if (Reactor[reactor_count].Energy == "Isothermal")
								  {
									  cout << "ERROR at isothermal reactor " << Reactor[reactor_count].Name << "! " << Reactor_keyword[j] << " is missing!";
									  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
								  }


							  }



							  else if (j == 5) //UA, only needed for heatexchanger reactors
							  {
								  if (Reactor[reactor_count].Energy == "HeatExchanger")
								  {
									  cout << "ERROR at heat-exchanger reactor " << Reactor[reactor_count].Name << "! " << Reactor_keyword[j] << " is missing!";
									  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
								  }

							  }
						  }
					  }

					  else if (j == 4 && R_keyword_controlling[j] == true)  // if temperature is found but the reactor is not isothermal, then the program finishes

					  {
						  if (Reactor[reactor_count].Energy != "Isothermal")
						  {
							  cout << "ERROR!" << Reactor[reactor_count].Energy << " Reactor " << Reactor[reactor_count].Name << " does not need a Temperature!\n";
							  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
						  }
					  }

					  else if (j == 5 && R_keyword_controlling[j] == true) // if UA is found but the reactor is not heatexchanger, then the program finishes
					  {


						  if (Reactor[reactor_count].Energy != "HeatExchanger")
						  {
							  cout << "ERROR!" << Reactor[reactor_count].Energy << " Reactor " << Reactor[reactor_count].Name << " does not need UA!\n";
							  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
						  }

					  }


				  }

				  if (Reactor[reactor_count].Out_stream == Reactor[reactor_count].In_stream)  // if the input stream has the same name of the output stream, the the program finishes
				  {
					  cout << "ERROR at reactor " << Reactor[reactor_count].Name << "!\n";
					  cout << "Input and output stream have the same number \"" << Reactor[reactor_count].Out_stream << "\".";
					  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
				  }



				  drawing << Reactor[reactor_count].Name << " [shape=";  // Graphviz instructions


				  if (Reactor[reactor_count].Type == "PSR")
				  {
					  drawing << "circle,height =40,";
				  }
				  else if (Reactor[reactor_count].Type == "PFR")
				  {
					  drawing << "box,height = 48,width = 32,";
				  }

				  if (Reactor[reactor_count].Energy == "Isothermal")
				  {
					  drawing << "color = blue];" << endl;
				  }
				  else if (Reactor[reactor_count].Energy == "HeatExchanger")
				  {
					  drawing << "color = red, style = dashed];" << endl;
				  }
				  else if (Reactor[reactor_count].Energy == "Adiabatic")
				  {
					  drawing << "color = red];" << endl;
				  }
				  reactor_count++;

			  }//if(line.find("Reactor")




			  else if (Splitted_line.at(0) == Main_keyword.at(1))//Getting MIXER information

			  {



				  if (LengthControl(Splitted_line, 2, line_number) != 0)

				  {
					  return 1;
				  }

				  if (FirstLetterThenNumberControl(Splitted_line.at(1), 'M') == "No")

				  {

					  cout << "ERROR at line " << line_number << "!";
					  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
				  }

				  Mixer[mixer_count].Name = Splitted_line.at(1);
				  Unit_name.push_back(Splitted_line.at(1));

				  for (int i = 0; i < mixer_count; i++) // checking there are not two mixer with the same name
				  {
					  if (Mixer[mixer_count].Name == Mixer[i].Name)
					  {
						  cout << "ERROR at line " << line_number << "!\n";
						  cout << "Two mixers are identified by the same name " << Mixer[mixer_count].Name << ".";

						  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
					  }
				  }

				  vector <bool> M_keyword_controlling(Mixer_keyword.size(), false);//every elements of this vector turns true if the mixer keyword with the same index is found

				  while (line.find("//") == string::npos)
				  {

					  getline(infile, line);
					  line_number++;

					  Splitted_line = String_comma_separator(line);

					  bool Unrecognized_keyword = true;// This boolean turns false if the first element of the line is one of the mixer keywords

					  if (Splitted_line.size() != 0) // Checking the line is not empty
					  {

						  //                     if ( New_declaration(line,line_number) != 0)
						  //                    {
						  //                        return 1;
						  //                    }

						  //                    if (line.find(";")!= string::npos)
						  //                    {
						  //                        cout<<"ERROR at line "<<line_number<<"!You don't have to use \';\' .\nNo list is required for mixer "<<Mixer[mixer_count].Name<<".\n";
						  //                        cout<<"Remember that input streams must be separated by a comma.";
						  //                        return 1;
						  //                    }
						  for (int j = 0; j < Mixer_keyword.size(); j++) // for cycle: it is used to check that the first element of the line is a Mixer_keyword

						  {
							  if (Splitted_line.at(0) == Mixer_keyword[j])

							  {

								  Unrecognized_keyword = false;

								  if (M_keyword_controlling[j] == true)

								  {
									  cout << "ERROR at line " << line_number << "!\n";
									  cout << "Mixer " << Mixer[mixer_count].Name << " has more than one keyword " << Mixer_keyword[j] << "!";
									  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
								  }


								  M_keyword_controlling[j] = true;

								  //                        switch ( j )
								  //                                
								  //                        {
								  //                                case 0://Type
								  //                                    
								  //                                    if ( LengthControl(Splitted_line,2,line_number)!= 0)
								  //                                    {
								  //                                        cout<<"\nThis line should be :"<<Mixer_keyword[j]<<" (text string)";
								  //                                        return 1;
								  //                                    }
								  //                                    break;
								  //                                    
								  //                                case 1://Inlet_stream
								  //                                    if ( Splitted_line.size() < 2)
								  //                                    {
								  //                                        cout<<"ERROR at line "<<line_number<<endl;
								  //                                        cout<<"Input streams of mixer "<<Mixer[mixer_count].Name<<" must be more than one!";
								  //                                        return 1;
								  //                                    }
								  //                                    
								  //                                    break;
								  //                                    
								  //                                case 2://Outlet_stream
								  //                                    if ( LengthControl(Splitted_line,2,line_number)!= 0)
								  //                                    {
								  //                                        cout<<"Output stream of mixer "<<Mixer[mixer_count].Name<<" must be only one!";
								  //                                        return 1;
								  //                                    }
								  //                                    
								  //                        }

								  if (j == 0) //Energy        
								  {

									  if (LengthControl(Splitted_line, 2, line_number) != 0)

									  {
										  cout << "\nThis line should be :" << Mixer_keyword[j] << " (text string)";  //Energy Isothermal
										  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
									  }


									  if ((DictionaryControl(Splitted_line.at(1), Mixer_energy_dictionary)) != 0)
									  {
										  cout << "\nLine " << line_number << ".";
										  return 1;
									  }
									  Mixer[mixer_count].Energy = Splitted_line.at(1);


								  }

								  else if (j == 1) //Inlet_stream      
								  {

									  if (Splitted_line.size() < 2) // if the size is  < 2 , then  there is only one or no Inlet_stream and the program finishes

									  {
										  cout << "ERROR at line " << line_number << endl;
										  cout << "Input streams of mixer " << Mixer[mixer_count].Name << " must be more than one!";
										  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
									  }

									  for (int i = 1; i < Splitted_line.size(); i++)

									  {

										  if (IntControl(Splitted_line.at(i)) == "No")
										  {
											  cout << "ERROR at line " << line_number << "!\n";
											  cout << "Stream number \"" << Splitted_line.at(i) << "\" must be a positive integer!";
											  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
										  }

										  int number;
										  stringstream(Splitted_line.at(i)) >> number;

										  Mixer[mixer_count].In_stream.push_back(number);
										  Global_stream_number.push_back(number);


										  if (Mixer[mixer_count].InputStreamCheck() != 0)
										  {
											  cout << "\nLook at line " << line_number;
											  return 1;
										  }
									  }


								  }

								  else if (j == 2)//Outlet_stream
								  {

									  if (LengthControl(Splitted_line, 2, line_number) != 0)

									  {
										  cout << "Output stream of mixer " << Mixer[mixer_count].Name << " must be only one!";
										  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
									  }

									  if (IntControl(Splitted_line.at(1)) == "No")

									  {
										  cout << "ERROR at line " << line_number << "!\n";
										  cout << "Stream number \"" << Splitted_line.at(1) << "\" must be a positive integer!";
										  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
									  }

									  int number;
									  stringstream(Splitted_line.at(1)) >> number;

									  Mixer[mixer_count].Out_stream = number;
									  Global_stream_number.push_back(number);
								  }

								  else if (j == 3) // Temperature
								  {
									  if (LengthControl(Splitted_line, 3, line_number) != 0)

									  {
										  cout << "This line should be: Temperature (value) (unit of measurement)";
										  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
									  }

									  if (DoubleControl(Splitted_line.at(1)) == "No")
									  {
										  cout << "ERROR at line " << line_number << "!\n";
										  cout << "Temperature value is incorrect!";
										  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
									  }

									  if (LetterControl(Splitted_line.at(2)) == "No")
									  {
										  cout << "ERROR at line " << line_number << "!\n";
										  cout << "Temperature unit of measurement is incorrect!";
										  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
									  }


									  double Number;
									  stringstream(Splitted_line.at(1)) >> Number;

									  Mixer[mixer_count].Temperature = Number;

									  Mixer[mixer_count].Temperature_unit_of_measurement = Splitted_line.at(2);

									  if (Number < 0 && Splitted_line.at(2) == "K")

									  {
										  cout << "ERROR at line " << line_number << "!\n";
										  cout << "Temperature in K cannot be negative!";
										  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
									  }
								  }

							  }


						  } // end of for cycle

						  if (Unrecognized_keyword == true && Splitted_line.at(0) != "//") // CHECK: if Unrecognized_keyword == true then The first element of a line 
						  {                                                                   // is not a Mixer_keyword and therefore the program finishes

							  cout << "ERROR at line " << line_number << "!\n";

							  for (int i = 0; i < Main_keyword.size(); i++)
							  {
								  if (Splitted_line.at(0) == Main_keyword.at(i))
								  {
									  cout << "This is a new " << Main_keyword.at(i) << " declaration but the last one is not finished yet!\n";
									  cout << "Remember to end unit or stream declarations with \'//\'.";
									  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
								  }
							  }

							  cout << "\"" << Splitted_line.at(0) << "\" is not a mixer keyword!";
							  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
						  }

						  if (End_of_file((infile)) != 0)
						  {

							  return 1;
						  }


					  }

				  } // end of while loop

				  for (int i = 0; i < M_keyword_controlling.size(); i++) // Checking the presence of every compulsory mixer keyword
				  {
					  if (M_keyword_controlling[i] == false)
					  {

						  if (i != 3) // excluding Temperature
						  {
							  cout << "ERROR at mixer " << Mixer[mixer_count].Name << "!\n";
							  cout << Mixer_keyword[i] << " is missing!";
							  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
						  }
						  if (i == 3 && Mixer[mixer_count].Energy == "Isothermal") //if Temperature is absent and mixer is Isothermal, the program finishes
						  {
							  cout << "ERROR at isothermal mixer " << Mixer[mixer_count].Name << "!Temperature is required\n";
							  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
						  }

					  }

					  if (M_keyword_controlling[3] == true && Mixer[mixer_count].Energy != "Isothermal") //if Temperature is present but mixer is not Isothermal, the program finishes
					  {
						  cout << "ERROR! Mixer " << Mixer[mixer_count].Name << " is not isothermal and therefore does not need a Temperature!\n";

						  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
					  }


				  }

				  if (Mixer[mixer_count].InputOutputStreamCheck() != 0)
				  {

					  return 1;
				  }


				  drawing << Mixer[mixer_count].Name << " [shape=invtriangle,";  // Graphviz instructions
																				 //                drawing << Mixer[mixer_count].Name <<" [shape=point,";  // Graphviz

				  if (Mixer[mixer_count].Energy == "Isothermal")
				  {
					  drawing << "color = blue];" << endl;
				  }
				  else if (Mixer[mixer_count].Energy == "Adiabatic")
				  {
					  drawing << "color = red];" << endl;
				  }

				  mixer_count++;

			  } //End of the MIXER block


			  else if (Splitted_line.at(0) == Main_keyword.at(2))//Getting SPLITTER information

			  {

				  vector <bool> S_keyword_controlling(Splitter_keyword.size(), false);

				  if (LengthControl(Splitted_line, 2, line_number) != 0)

				  {
					  return 1;
				  }

				  if (FirstLetterThenNumberControl(Splitted_line.at(1), 'S') == "No")

				  {

					  cout << "ERROR at line " << line_number;
					  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
				  }

				  Splitter[splitter_count].Name = Splitted_line.at(1);
				  Unit_name.push_back(Splitted_line.at(1));



				  for (int i = 0; i < splitter_count; i++) // Checking that two splitters do not have the same name
				  {
					  if (Splitter[splitter_count].Name == Splitter[i].Name)
					  {
						  cout << "ERROR at line " << line_number << "!\n";
						  cout << "Two splitters are identified by the same name " << Splitter[splitter_count].Name << ".";

						  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
					  }
				  }




				  while (line.find("//") == string::npos)
				  {

					  getline(infile, line);
					  line_number++;

					  Splitted_line = String_comma_separator(line);



					  //                     if ( New_declaration(line,line_number) != 0)
					  //                    {
					  //                        return 1;
					  //                    }

					  if (Splitted_line.size() != 0)

					  {
						  bool Unrecognized_keyword = true;


						  for (int j = 0; j < Splitter_keyword.size(); j++) //for cycle

						  {
							  if (Splitted_line.at(0) == Splitter_keyword.at(j))

							  {

								  Unrecognized_keyword = false;

								  if (S_keyword_controlling.at(j) == true)

								  {
									  cout << "ERROR at line " << line_number << "!\n";
									  cout << "Splitter " << Splitter[splitter_count].Name << " has more than one keyword " << Splitter_keyword[j] << "!";
									  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
								  }

								  S_keyword_controlling.at(j) = true;

								  if (j == 0)//Inlet_stream     
								  {

									  if (LengthControl(Splitted_line, 2, line_number) != 0)

									  {
										  cout << "Input stream of splitter " << Splitter[splitter_count].Name << " must be only one!";
										  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
									  }

									  if (IntControl(Splitted_line.at(1)) == "No")

									  {
										  cout << "ERROR at line " << line_number << "!\n";
										  cout << "Stream number \"" << Splitted_line.at(1) << "\" must be a positive integer!";
										  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
									  }

									  int number;
									  stringstream(Splitted_line.at(1)) >> number;


									  Splitter[splitter_count].In_stream = number;
									  Global_stream_number.push_back(number);

								  }

								  else if (j == 1)//Outlet_stream
								  {

									  Splitted_line.erase(Splitted_line.begin());

									  double Summing_splitting_ratio = 0;

									  int loop_number = 0;



								  Splitter_Outlet_stream_loop: // this loop gets the stream number together with it splitting ratio

									  loop_number++;

									  if (LengthControl(Splitted_line, 2, line_number) != 0)
									  {
										  cout << "ERROR!Every line following Outlet_stream of splitter " << Splitter[splitter_count].Name << " should be (stream number),(splitting ratio)!\n";
										  cout << "End the list with a ;";
										  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
									  }

									  if (IntControl(Splitted_line.at(0)) == "No")
									  {
										  cout << "Stream number must be an integer!Remember to end the list of output streams of splitter " << Splitter[splitter_count].Name << " with a ';'";
										  cout << "\nERROR at line " << line_number;
										  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
									  }

									  if (DoubleControl(Splitted_line.at(1)) == "No")
									  {
										  cout << "Splitting ratio must be a double!Remember to end the list of output streams of splitter " << Splitter[splitter_count].Name << " with a ';'\n";
										  cout << "\nERROR at line " << line_number;
										  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
									  }

									  int int_number;
									  stringstream(Splitted_line.at(0)) >> int_number;

									  Splitter[splitter_count].Out_stream.push_back(int_number);
									  Global_stream_number.push_back(int_number);

									  if (Splitter[splitter_count].OutputStreamCheck() != 0)
									  {
										  cout << "\nLook at line " << line_number;
										  return 1;
									  }



									  double number;
									  stringstream(Splitted_line.at(1)) >> number;

									  if (number < 0)

									  {
										  cout << "ERROR at line " << line_number;
										  cout << "\nSplitting ratio cannot be negative!";
										  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
									  }

									  Summing_splitting_ratio += number;

									  Splitter[splitter_count].Splitting_ratio.push_back(number);

									  //                            if ( line.find("//") != string::npos)
									  //                            {
									  //                                cout<<"ERROR at splitter "<<Splitter[splitter_count].Name<<"!";
									  //                                cout<<"List of output streams does not end with a ;!";
									  //                                return 1;
									  //                            }

									  if (line.find(";") == string::npos) //Loop  condition: if  this line doesn't contain any ';', then the loop restarts

									  {
										  getline(infile, line);
										  line_number++;
										  Splitted_line = String_comma_separator(line);
										  goto Splitter_Outlet_stream_loop;

									  }

									  if (loop_number < 2) //Splitter requires at least two output streams, that means at least 2 loops
									  {
										  cout << "ERROR at Splitter " << Splitter[splitter_count].Name << "!\n";
										  cout << "Slitter requires at least two output streams!";
										  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
									  }


									  if (abs(Summing_splitting_ratio - 1) > 0.0000000001) //CHECK: splitting ratios must sum to 1
									  {
										  cout << "ERROR at splitter " << Splitter[splitter_count].Name << "!\n";
										  cout << "Splitting ratios do not sum to 1!";
										  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
									  }








								  }


							  }


						  } // end of for cycle



						  if (Unrecognized_keyword == true && Splitted_line.at(0) != "//") // CHECK: if Unrecognized_keyword == true then The first element of a line 
						  {                                                                   // is not a Splitter_keyword and therefore the program finishes

							  cout << "ERROR at line " << line_number << "!\n";

							  for (int i = 0; i < Main_keyword.size(); i++)
							  {
								  if (Splitted_line.at(0) == Main_keyword.at(i))
								  {
									  cout << "This is a new " << Main_keyword.at(i) << " declaration but the last one is not finished yet!\n";
									  cout << "Remember to end unit or stream declarations with \'//\'.";
									  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
								  }
							  }

							  cout << "\"" << Splitted_line.at(0) << "\" is not a splitter keyword!";
							  return 1;
						  }


						  if (End_of_file((infile)) != 0)
						  {

							  return 1;
						  }

					  }


				  }//end of while loop


				  for (int i = 0; i < S_keyword_controlling.size(); i++)
				  {
					  if (S_keyword_controlling[i] == false)
					  {
						  cout << "ERROR at splitter " << Splitter[splitter_count].Name << "!\n";
						  cout << Splitter_keyword[i] << " is missing!";
						  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");

					  }
				  }

				  if (Splitter[splitter_count].OutputInputStreamCheck() != 0)

				  {
					  return 1;
				  }

				  drawing << Splitter[splitter_count].Name << " [shape=triangle];" << endl; // Graphviz
																							//                 drawing << Splitter[splitter_count].Name <<" [shape=point];"<<endl; // Graphviz

				  splitter_count++;

			  }





			  else if (Splitted_line.at(0) == Main_keyword.at(3)) //Getting PHASE-SPLITTER information

			  {


				  vector <bool> PS_keyword_controlling(PhaseSplitter_keyword.size(), false);


				  if (LengthControl(Splitted_line, 2, line_number) != 0)

				  {
					  return 1;
				  }

				  if (FirstLetterThenNumberControl(Splitted_line.at(1), 'P') == "No")

				  {

					  cout << "ERROR at line " << line_number;
					  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
				  }

				  PhaseSplitter[phasesplitter_count].Name = Splitted_line.at(1);
				  Unit_name.push_back(Splitted_line.at(1));




				  for (int i = 0; i < phasesplitter_count; i++) // Checking two phase-splitter units do not have the same name

				  {
					  if (PhaseSplitter[phasesplitter_count].Name == PhaseSplitter[i].Name)
					  {
						  cout << "ERROR at line " << line_number << "!\n";
						  cout << "Two Phase-splitters are identified by the same name: " << PhaseSplitter[phasesplitter_count].Name << ".";

						  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
					  }
				  }




				  while (line.find("//") == string::npos)

				  {

					  getline(infile, line);
					  line_number++;

					  Splitted_line = String_comma_separator(line);



					  if (Splitted_line.size() != 0)

					  {


						  bool Unrecognized_keyword = true;


						  for (int j = 0; j < PhaseSplitter_keyword.size(); j++)

						  {
							  if (Splitted_line.at(0) == PhaseSplitter_keyword.at(j))

							  {

								  Unrecognized_keyword = false;

								  if (PS_keyword_controlling.at(j) == true)

								  {
									  cout << "ERROR at line " << line_number << "!\n";
									  cout << "Phase-splitter " << PhaseSplitter[phasesplitter_count].Name << " has more than one keyword " << PhaseSplitter_keyword[j] << "!";
									  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
								  }
								  PS_keyword_controlling[j] = true;



								  if (j == 0)//Inlet_stream     
								  {

									  if (LengthControl(Splitted_line, 2, line_number) != 0)

									  {
										  cout << "Input stream of Phase-Splitter " << PhaseSplitter[phasesplitter_count].Name << " must be only one!";
										  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
									  }

									  if (IntControl(Splitted_line.at(1)) == "No")

									  {
										  cout << "ERROR at line " << line_number << "!\n";
										  cout << "Stream number \"" << Splitted_line.at(1) << "\" must be a positive integer!";
										  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
									  }

									  int number;
									  stringstream(Splitted_line.at(1)) >> number;

									  PhaseSplitter[phasesplitter_count].In_stream = number;
									  Global_stream_number.push_back(number);

								  }

								  else if (j == 1)//Outlet_stream


								  {

									  Splitted_line.erase(Splitted_line.begin());


									  int loop_number = 0;

								  PhaseSplitter_Outlet_stream_loop: // This loop gets the stream name and the phase of this stream

									  loop_number++;

									  if (LengthControl(Splitted_line, 2, line_number) != 0)

									  {

										  cout << "ERROR!Every line following Outlet_stream of Phase-Splitter " << PhaseSplitter[phasesplitter_count].Name << " should be (stream number),(phase)!\n";
										  cout << "End the list with a ;";
										  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
									  }

									  if (IntControl(Splitted_line.at(0)) == "No")

									  {
										  cout << "Stream number must be an integer!Remember to end the list of output streams of Phase-Splitter " << PhaseSplitter[phasesplitter_count].Name << " with a ';'";
										  cout << "\nERROR at line " << line_number;
										  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
									  }

									  int number;
									  stringstream(Splitted_line.at(0)) >> number;

									  PhaseSplitter[phasesplitter_count].Out_stream.push_back(number);
									  Global_stream_number.push_back(number);



									  if ((DictionaryControl(Splitted_line.at(1), Splitted_phase_dictionary)) != 0)
									  {
										  cout << "\nLine " << line_number << ".";
										  return 1;
									  }

									  PhaseSplitter[phasesplitter_count].Splitted_phase.push_back(Splitted_line.at(1));

									  if (PhaseSplitter[phasesplitter_count].OutputStreamCheck() != 0)
									  {
										  cout << "\nLook at line " << line_number;
										  return 1;
									  }


									  if (line.find(";") == string::npos)

									  {
										  getline(infile, line);
										  line_number++;
										  Splitted_line = String_comma_separator(line);

										  goto PhaseSplitter_Outlet_stream_loop;

									  }


									  if (loop_number == 1)
									  {
										  cout << "ERROR at Phase-Splitter " << PhaseSplitter[phasesplitter_count].Name << "!\n";
										  cout << "Phase-Slitter requires at least two output streams!";
										  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
									  }


								  }


							  }


						  } // end of for cycle

						  if (Unrecognized_keyword == true && Splitted_line.at(0) != "//") // CHECK: if Unrecognized_keyword == true then The first element of a line 
						  {                                                                   // is not a Phase-splitter_keyword and therefore the program finishes

							  cout << "ERROR at line " << line_number << "!\n";

							  for (int i = 0; i < Main_keyword.size(); i++)
							  {
								  if (Splitted_line.at(0) == Main_keyword.at(i))
								  {
									  cout << "This is a new " << Main_keyword.at(i) << " declaration but the last one is not finished yet!\n";
									  cout << "Remember to end unit or stream declarations with \'//\'.";
									  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
								  }
							  }

							  cout << "\"" << Splitted_line.at(0) << "\" is not a Phase-splitter keyword!";
							  return 1;
						  }


						  if (End_of_file((infile)) != 0)
						  {

							  return 1;
						  }

					  }

				  }// end of while loop

				  for (int i = 0; i < PS_keyword_controlling.size(); i++)

				  {

					  if (PS_keyword_controlling[i] == false)

					  {

						  cout << "ERROR at Phase-Splitter " << PhaseSplitter[phasesplitter_count].Name << "!\n";
						  cout << PhaseSplitter_keyword[i] << " is missing!";
						  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");

					  }
				  }

				  if (PhaseSplitter[phasesplitter_count].OutputInputStreamCheck() != 0)
				  {
					  return 1;
				  }

				  drawing << PhaseSplitter[phasesplitter_count].Name << " [shape=triangle,color = midnightblue];" << endl;

				  phasesplitter_count++;

			  } //end of PHASE-SPLITTER block





			  else if (Splitted_line.at(0) == Main_keyword.at(4)) //Getting STREAM information

			  {

				  if (LengthControl(Splitted_line, 2, line_number) != 0)

				  {

					  return 1;
				  }

				  if (IntControl(Splitted_line.at(1)) == "No")

				  {
					  cout << "ERROR at line " << line_number << "!\n";
					  cout << "Stream number must be an integer!";
					  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");

				  }

				  int int_number;
				  stringstream(Splitted_line.at(1)) >> int_number;
				  Stream[stream_count].Name = int_number;
				  Stream[stream_count].GivingSizeToKeyword_presence(Stream_keyword.size());

				  for (int i = 0; i < stream_count; i++) // Checking two streams do not have the same number

				  {
					  if (Stream[stream_count].Name == Stream[i].Name)

					  {
						  cout << "ERROR at line " << line_number << "!\n";
						  cout << "Two Streams are identified by the same number: " << Stream[stream_count].Name << ".";

						  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
					  }
				  }


				  while (line.find("//") == string::npos)

				  {
					  getline(infile, line);
					  line_number++;




					  Splitted_line = String_comma_separator(line);

					  if (Splitted_line.size() != 0)

					  {

						  bool Unrecognized_keyword = true;



						  for (int j = 0; j < Stream_keyword.size(); j++)

						  {

							  if (Splitted_line.at(0) == Stream_keyword.at(j))

							  {

								  Unrecognized_keyword = false;



								  if (Stream[stream_count].TrueKeyword(j) != 0) //if a keyword is repeated twice, the program finishes
								  {
									  cout << "ERROR at line " << line_number << "!\n";
									  cout << "Stream " << Stream[stream_count].Name << " keyword \"" << Stream_keyword.at(j) << "\" is repeated twice!";
									  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");

								  }



								  if (j == 0)//Phase

								  {
									  if (LengthControl(Splitted_line, 2, line_number) != 0)

									  {
										  return 1;
									  }

									  if ((DictionaryControl(Splitted_line.at(1), Stream_phase_dictionary)) != 0)

									  {
										  cout << "\nLine " << line_number << ".";
										  return 1;
									  }

									  Stream[stream_count].Phase = Splitted_line.at(1);

								  }

								  if (j == 1)//SolidType

								  {


									  //                                    Splitted_line.erase(Splitted_line.begin());
									  //
									  //                                    
									  //                                    if ( LengthControl(Splitted_line,1,line_number)!= 0)
									  //                                    {
									  //
									  //                                        cout<<"Every line following SolidType of stream "<< Stream[stream_count].Name<<" should be (Cellulose,Lignin...)!\n";
									  //                                        cout<<"End the list with a ;";
									  //                                        return 1;
									  //                                    }

									  if (LengthControl(Splitted_line, 2, line_number) != 0)
									  {

										  cout << "This line should be " << Stream_keyword.at(j) << " (text string)!\n";

										  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
									  }


									  //                                    if((DictionaryControl(Splitted_line.at(0),Stream_solid_type_dictionary)) != 0)
									  //                                    {
									  //                                        cout<<"\nLine "<<line_number<<".";
									  //                                        return 1;
									  //                                    }

									  if ((DictionaryControl(Splitted_line.at(1), Stream_solid_type_dictionary)) != 0)
									  {
										  cout << "\nLine " << line_number << ".";
										  return 1;
									  }


									  Stream[stream_count].SolidType = Splitted_line.at(1);
									  //                                    Stream[stream_count].SolidType = Splitted_line.at(0);


									  //                                    Type_stream_loop:
									  //
									  //
									  //
									  //                                    if ( LengthControl(Splitted_line,1,line_number)!= 0)
									  //                                    {
									  //
									  //                                        cout<<"Every line following Type of stream "<< Stream[stream_count].Name<<" should be (Cellulose,Lignin...)!\n";
									  //                                        cout<<"End the list with a ;";
									  //                                        return 1;
									  //                                    }
									  //
									  //
									  //                                    if((DictionaryControl(Splitted_line.at(0),Stream_solid_type_dictionary)) != 0)
									  //                                    {
									  //                                        cout<<"\nLine "<<line_number<<".";
									  //                                        return 1;
									  //                                    }
									  //
									  //                                    Stream[stream_count].Type.push_back(Splitted_line.at(0));
									  //
									  //
									  //                                    if ( line.find(";") == string::npos)
									  //
									  //                                    {
									  //                                        getline(infile,line);
									  //                                        line_number++;
									  //                                        Splitted_line = String_comma_separator(line);
									  //
									  //                                        goto Type_stream_loop;
									  //
									  //                                    }





								  }

								  if (j == 2)//MassFlowRate

								  {

									  if (LengthControl(Splitted_line, 3, line_number) != 0)
									  {

										  cout << "\nThis line should contain (double) (text string)";

										  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
									  }

									  if (DoubleControl(Splitted_line.at(1)) == "No")
									  {
										  cout << "ERROR at line " << line_number << "!\n";
										  cout << "MassFlowRate value is incorrect!";
										  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
									  }

									  if (LetterControl(Splitted_line.at(2)) == "No")
									  {
										  cout << "ERROR at line " << line_number << "!\n";
										  cout << "MassFlowRate unit of measurement is incorrect!";
										  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
									  }

									  double Number;
									  stringstream(Splitted_line.at(1)) >> Number;
									  if (Number < 0)
									  {
										  cout << "ERROR at line " << line_number << "!\n";
										  cout << "MassFlowRate cannot be negative!";
										  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
									  }
									  Stream[stream_count].MassFlowRate = Number;
									  Stream[stream_count].MassFlowRate_unit_of_measurement = Splitted_line.at(2);

								  }

								  if (j == 3)//Temperature

								  {

									  if (LengthControl(Splitted_line, 3, line_number) != 0)
									  {

										  cout << "\nThis line should be Temperature (number) (unit of measurement)";

										  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
									  }

									  if (DoubleControl(Splitted_line.at(1)) == "No")
									  {
										  cout << "ERROR at line " << line_number << "!\n";
										  cout << "Temperature value is incorrect!";
										  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
									  }

									  if (LetterControl(Splitted_line.at(2)) == "No")
									  {
										  cout << "ERROR at line " << line_number << "!\n";
										  cout << "Temperature unit of measurement is incorrect!";
										  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
									  }


									  double Number;
									  stringstream(Splitted_line.at(1)) >> Number;

									  Stream[stream_count].Temperature = Number;

									  Stream[stream_count].Temperature_unit_of_measurement = Splitted_line.at(2);

									  if (Number < 0 && Splitted_line.at(2) == "K")

									  {
										  cout << "ERROR at line " << line_number << "!\n";
										  cout << "Temperature in K cannot be negative!";
										  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
									  }


								  }


								  if (j == 4)//MassFraction
								  {

									  Splitted_line.erase(Splitted_line.begin());



								  Mass_fraction_loop:



									  bool Controlling_species = false;

									  if (LengthControl(Splitted_line, 2, line_number) != 0)

									  {

										  cout << "Every line following MassFraction of stream " << Stream[stream_count].Name << " should be (Compound) (mass fraction)!\n";
										  cout << "End the list with a ;";
										  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
									  }




									  Stream[stream_count].Compound.push_back(Splitted_line.at(0));

									  if (Stream[stream_count].CompoundControl() != 0)
									  {
										  return 1;
									  }



									  for (int i = 0; i < Species.size(); i++) //Checking if the new compound has already been found

									  {
										  if ((Splitted_line.at(0)).compare(Species.at(i)) == 0)

										  {
											  Controlling_species = true;

											  break;
										  }

									  }

									  if (Controlling_species == false) // if not, then it is added to Species vector

									  {
										  Species.push_back(Splitted_line.at(0));
									  }


									  if (DoubleControl(Splitted_line.at(1)) == "No") //Checking Mass fraction value is a double
									  {
										  cout << "ERROR at line " << line_number << "!\n";
										  cout << "Mass Fraction value of compound \"" << Stream[stream_count].Compound.back() << "\" must be a double.";
										  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
									  }


									  double number;
									  stringstream(Splitted_line.at(1)) >> number;

									  if (number < 0)//Checking Mass fraction value is positive
									  {
										  cout << "ERROR at line " << line_number << "!\n";
										  cout << "Mass Fraction value of compound \"" << Stream[stream_count].Compound.back() << "\" must be a positive double.";
										  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
									  }
									  Stream[stream_count].Composition.push_back(number);



									  //                                    if (Stream[stream_count].SummingMassFractions() != 0)
									  //                                    {
									  //                                        return 1;
									  //                                    }


									  if (line.find(";") == string::npos)

									  {
										  getline(infile, line);
										  line_number++;
										  Splitted_line = String_comma_separator(line);

										  goto Mass_fraction_loop;

									  }
								  }





							  }

						  } //end of for cycle



						  if (Unrecognized_keyword == true && Splitted_line.at(0) != "//") // CHECK: if Unrecognized_keyword == true then The first element of a line 
						  {                                                                   // is not a Stream_keyword and therefore the program finishes

							  cout << "ERROR at line " << line_number << "!\n";

							  for (int i = 0; i < Main_keyword.size(); i++)
							  {
								  if (Splitted_line.at(0) == Main_keyword.at(i))
								  {
									  cout << "This is a new " << Main_keyword.at(i) << " declaration but the last one is not finished yet!\n";
									  cout << "Remember to end unit or stream declarations with \'//\'.";
									  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
								  }
							  }

							  cout << "\"" << Splitted_line.at(0) << "\" is not a Stream keyword!";
							  return 1;
						  }


						  if (End_of_file((infile)) != 0)
						  {

							  return 1;
						  }



					  }


				  }//end of while loop


				  if (abs(Stream[stream_count].MassFractionSum() - 1) > 0.00000000001)
				  {
					  cout << "ERROR! Mass fraction sum of stream " << Stream[stream_count].Name << " is more than 1!";

					  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
				  }

				  stream_count++;

			  } // End of STREAM block



			  else // if Splitted_line.at(0) is not Reactor,Mixer,Splitter,Phase-splitter,Stream

			  {


				  cout << "ERROR at line " << line_number << "!\n";
				  cout << "Main Keyword \"" << Splitted_line.at(0) << "\" not recognized!";
				  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");


			  }



		  }// if splitted_line.size() != 0

	  } //end of the INPUT FILE reading



	  infile.close(); //Closing input_file

	  Splitted_line.clear();



	  cout << endl << endl;


	  Reactor.resize(reactor_count);  //Resizing and shrinking vectors in order to occupy the right memory
	  Reactor.shrink_to_fit();

	  Mixer.resize(mixer_count);
	  Mixer.shrink_to_fit();

	  Splitter.resize(splitter_count);
	  Splitter.shrink_to_fit();

	  PhaseSplitter.resize(phasesplitter_count);
	  PhaseSplitter.shrink_to_fit();

	  Stream.resize(stream_count);
	  Stream.shrink_to_fit();


	  for (int i = 0; i < Global_stream_number.size(); i++) //Eliminating streams that are repeated twice from Global_stream_number
	  {

		  bool Check = false;


		  for (int j = 0; j < Stream_number.size(); j++)

		  {

			  if (Stream_number.at(j) == Global_stream_number.at(i))
			  {

				  Check = true;
				  break;
			  }

		  }

		  if (Check == false)

		  {

			  Stream_number.push_back(Global_stream_number.at(i));
		  }


	  }





	  for (int i = 0; i < stream_count; i++) //Checking the validity of every phase splitter
	  {
		  for (int j = 0; j < phasesplitter_count; j++)
		  {
			  if (Stream[i].Name == PhaseSplitter[j].In_stream)
			  {
				  if (Stream[i].Phase != "Mix")
				  {
					  cout << "ERROR! Stream " << PhaseSplitter[j].In_stream << " enters Phase-splitter " << PhaseSplitter[j].Name << " and therefore its Phase must be \"Mix\".";
					  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
				  }
			  }

			  for (int x = 0; x < PhaseSplitter[j].Out_stream.size(); x++)
			  {
				  if (Stream[i].Name == PhaseSplitter[j].Out_stream.at(x))
				  {
					  if (Stream[i].Phase != PhaseSplitter[j].Splitted_phase.at(x))
					  {
						  cout << "ERROR! Stream " << Stream[i].Name << " phase is " << Stream[i].Phase << " but Phase-splitter " << PhaseSplitter[j].Name << " declares it is " << PhaseSplitter[j].Splitted_phase.at(x) << "!";
						  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
					  }
				  }
			  }
		  }



	  }




	  //4 vectors related to stream

	  From.resize(Stream_number.size(), "0");     //Unit FROM which the stream arrives
	  To.resize(Stream_number.size(), "0");       //Unit TO which the stream arrives
	  Splitting_ratio.resize(Stream_number.size(), 1);
	  Stream_phase.resize(Stream_number.size(), "Mix"); //Phase of every system stream, by default it is Mix







	  for (int i = 0; i < Stream_number.size(); i++)   //Filling the 4 vectors taking information from the the vector of classes previously created 

	  {

		  for (int x = 0; x < splitter_count; x++)
		  {
			  for (int j = 0; j < Splitter[x].Out_stream.size(); j++)
			  {
				  if (Stream_number[i] == Splitter[x].Out_stream[j])
				  {
					  if (From[i] != "0")
					  {
						  cout << "ERROR! Units " << Splitter[x].Name << " and " << From[i] << " have the same outlet stream: " << Stream_number[i] << " !\n";
						  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
					  }

					  Splitting_ratio[i] = Splitter[x].Splitting_ratio[j];

					  From[i] = Splitter[x].Name;

				  }

			  }

			  if (Stream_number[i] == Splitter[x].In_stream)
			  {
				  if (To[i] != "0")
				  {
					  cout << "ERROR! Units " << Splitter[x].Name << " and " << To[i] << " have the same inlet stream: " << Stream_number[i] << " !\n";
					  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
				  }

				  To[i] = Splitter[x].Name;

			  }


		  }

		  for (int x = 0; x < reactor_count; x++)
		  {
			  if (Stream_number[i] == Reactor[x].Out_stream)
			  {
				  if (From[i] != "0")
				  {
					  cout << "ERROR! Units " << Reactor[x].Name << " and " << From[i] << " have the same output stream: " << Stream_number[i] << " !\n";
					  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
				  }
				  From[i] = Reactor[x].Name;

			  }


			  if (Stream_number[i] == Reactor[x].In_stream)
			  {
				  if (To[i] != "0")
				  {
					  cout << "ERROR! Units " << Reactor[x].Name << " and " << To[i] << " have the same input stream: " << Stream_number[i] << " !\n";
					  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
				  }
				  To[i] = Reactor[x].Name;

			  }

		  }



		  for (int x = 0; x < mixer_count; x++)

		  {
			  if (Stream_number[i] == Mixer[x].Out_stream)

			  {
				  if (From[i] != "0")

				  {
					  cout << "ERROR!Units " << Mixer[x].Name << " and " << From[i] << " have the same output stream: " << Stream_number[i] << " !\n";
					  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
				  }

				  From[i] = Mixer[x].Name;

			  }



			  for (int j = 0; j < Mixer[x].In_stream.size(); j++)

			  {
				  if (Stream_number[i] == Mixer[x].In_stream[j])

				  {
					  if (To[i] != "0")

					  {
						  cout << "ERROR! Units " << Mixer[x].Name << " and " << To[i] << " have the same input stream: " << Stream_number[i] << " !\n";
						  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
					  }

					  To[i] = Mixer[x].Name;

				  }

			  }
		  }

		  for (int x = 0; x < phasesplitter_count; x++)

		  {
			  if (Stream_number[i] == PhaseSplitter[x].In_stream)

			  {
				  if (To[i] != "0")
				  {
					  cout << "ERROR! Units " << PhaseSplitter[x].Name << " and " << To[i] << " have the same input stream: " << Stream_number[i] << " !\n";
					  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
				  }
				  To[i] = PhaseSplitter[x].Name;

			  }

			  for (int j = 0; j < PhaseSplitter[x].Out_stream.size(); j++)

			  {
				  if (Stream_number[i] == PhaseSplitter[x].Out_stream[j])
				  {

					  if (From[i] != "0")
					  {
						  cout << "ERROR! Units " << PhaseSplitter[x].Name << " and " << From[i] << " have the same output stream: " << Stream_number[i] << " !\n";
						  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
					  }
					  From[i] = PhaseSplitter[x].Name;

					  Stream_phase.at(i) = PhaseSplitter[x].Splitted_phase[j];


				  }

			  }
		  }

		  for (int x = 0; x < stream_count; x++)

		  {

			  if (Stream_number.at(i) == Stream[x].Name)

			  {

				  Stream_phase.at(i) = Stream[x].Phase;

			  }
		  }
	  }

	  vector <bool> Controlling_input_output(2, false); // The first element of this vector will turn true if there is a system input stream
														//The second element will turn true if there is a system output stream




	  for (int i = 0; i < Stream_number.size(); i++)

	  {

		  if (From.at(i) == "0") // if stream at index i is a system input, then it must have a declaration

		  {
			  Controlling_input_output.at(0) = true;

			  bool Input_stream_check = false;

			  for (int x = 0; x < stream_count; x++) // checking over every stream declaration
			  {
				  if (Stream_number.at(i) == Stream[x].Name)
				  {
					  Input_stream_check = true;
					  break;
				  }


			  }

			  if (Input_stream_check == false) // if this stream has not been found, then the program finishes
			  {
				  cout << "ERROR! System input stream " << Stream_number.at(i) << " does not have a declaration!\n";
				  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
			  }


		  }

		  if (To.at(i) == "0")

		  {
			  Controlling_input_output.at(1) = true;

		  }

	  }







	  for (int i = 0; i < 2; i++)
	  {
		  if (Controlling_input_output[i] == false)
		  {
			  if (i == 0)
			  {
				  cout << "ERROR! There is no system input!\n";
				  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
			  }

			  if (i == 1)
			  {
				  cout << "ERROR! There is no system output!\n";
				  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
			  }
		  }
	  }

	  for (int i = 0; i < Stream_number.size(); i++) //Checking the validity of the stream keywords found
	  {

		  for (int j = 0; j < stream_count; j++)

		  {


			  if (Stream_number.at(i) == Stream[j].Name)

			  {
				  if (From.at(i) == "0")  //Considering system input streams!

				  {
					  for (int x = 1; x < Stream_keyword.size(); x++)   //Starting from x = 1 since Phase is not considered. the default phase is Mix

					  {
						  if (Stream[j].Keyword_presence.at(x) == false)

						  {

							  if (x == 1) //if there is not type but the stream phase is Solid, then the program finishes

							  {

								  if (Stream[j].Phase == "Solid")
								  {
									  cout << "ERROR!Solid Stream " << Stream[j].Name << " is a system input and therefore needs: " << Stream_keyword.at(x) << ".";
									  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");

								  }
							  }

							  else //Every other stream keyword must be present in the declaration of a system input stream

							  {
								  cout << "ERROR!Stream " << Stream[j].Name << " is a system input and therefore needs: " << Stream_keyword.at(x) << ".";
								  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");

							  }



						  }



					  }


					  if (Stream[j].Keyword_presence.at(1) == true && Stream[j].Phase != "Solid") //if there is a type but the Stream is not solid, then the program finishes
					  {
						  cout << "ERROR!Stream " << Stream[j].Name << " is " << Stream[j].Phase << " and therefore does not need: " << Stream_keyword.at(1) << ".";
						  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
					  }

					  if (Stream[j].Phase == "Solid")//Checking that a solid system input stream mass fraction is expressed in terms of elements

					  {

						  for (int x = 0; x < Stream[j].Compound.size(); x++)

						  {


							  if ((DictionaryControl(Stream[j].Compound.at(x), Solid_compound_dictionary)) != 0)//This check is for system input streams only
							  {
								  cout << "Change the related MassFraction line of SOLID stream " << Stream[j].Name << ".";
								  return 1;
							  }
						  }
					  }

					  if (abs(Stream[j].MassFractionSum() - 1) > 0.00000000001)  //Checking the MassFractionSum is 1, only for system input streams

					  {

						  cout << "ERROR!Stream " << Stream[j].Name << " is a system input and therefore its mass-fraction-sum must be 1!\n";
						  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
					  }
				  }

				  else // Considering non-input streams!
				  {
					  if (Stream[j].Keyword_presence.at(1) == true)
					  {
						  cout << "ERROR!Stream " << Stream[j].Name << " is  not a system input stream, therefore does not need a SolidType!";
						  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");

					  }
				  }
			  }



		  }
	  }





	  for (int j = 0; j < stream_count; j++) // Checking that every declarated stream is part of the system

	  {
		  Stream[j].ShrinkingKeyword_presence();

		  bool Control = false;

		  for (int i = 0; i < Stream_number.size(); i++)
		  {

			  if (Stream[j].Name == Stream_number.at(i))
			  {
				  Control = true;
				  break;
			  }

		  }

		  if (Control == false)
		  {
			  cout << "ERROR!Stream " << Stream[j].Name << " has not been defined in the system!";
			  OpenSMOKE::FatalErrorMessage("Error detected! Program will now exit!");
		  }

	  }


	  for (int i = 0; i < splitter_count; i++) // Checking the validity of every splitter: it must be mono-phase
	  {
		  string In_stream_phase;

		  vector <string> Out_stream_phase;
		  vector <int> remembering_out_stream;

		  for (int j = 0; j < Stream_number.size(); j++)

		  {
			  if (Stream_number.at(j) == Splitter[i].In_stream)

			  {
				  In_stream_phase = Stream_phase.at(j);
			  }

			  else
			  {
				  for (int x = 0; x < Splitter[i].Out_stream.size(); x++)
				  {
					  if (Stream_number.at(j) == Splitter[i].Out_stream.at(x))
					  {
						  Out_stream_phase.push_back(Stream_phase.at(j));

						  remembering_out_stream.push_back(Stream_number.at(j));
					  }
				  }
			  }
		  }



		  for (int w = 0; w < Out_stream_phase.size(); w++)
		  {



			  if (In_stream_phase != Out_stream_phase.at(w))
			  {
				  cout << "ERROR!Splitters involve always only one phase! This might be due to default assignment and might not be an actual error \n";
				  cout << "Splitter " << Splitter[i].Name << " has: \"" << In_stream_phase << "\" stream " << Splitter[i].In_stream << " and \"" << Out_stream_phase.at(w) << "\" stream " << remembering_out_stream.at(w) << "!";
				  continue;
			  }
		  }


	  }







	  for (int i = 0; i < Stream_number.size(); i++)  //GraphViz instructions
	  {

		  if (From.at(i) == "0")
		  {


			  drawing << "INPUT -> " << To.at(i);

		  }
		  else if (To.at(i) == "0")
		  {
			  drawing << From.at(i) << " -> OUTPUT";
		  }
		  else
		  {
			  drawing << From.at(i) << " -> " << To.at(i);
		  }

		  drawing << "[ ";


		  if (Stream_phase.at(i) != "Mix")
		  {
			  if (Stream_phase.at(i) == "Liquid")
			  {
				  drawing << "color = lightblue,";
			  }

			  else if (Stream_phase.at(i) == "Gas")
			  {
				  drawing << "color = orange,";
			  }

			  else if (Stream_phase.at(i) == "Solid")
			  {
				  drawing << "color = green,";
			  }
		  }

		  drawing << " label =\" " << Stream_number.at(i);





		  if (*(From.at(i).begin()) == 'S')
		  {

			  drawing << ", " << Splitting_ratio.at(i) << " \" ];" << endl;
		  }
		  else
		  {
			  drawing << " \" ];" << endl;
		  }






	  }



	  drawing << "}" << endl;

	  drawing.close();

	  // Generate PDF
	  system("dot -Tpdf NetworkMap_graphvizDot.dot -o NetworkMap.pdf ");



	  return 0;

  }


} /* END OF NAMESPACE OPENSMOKE */
