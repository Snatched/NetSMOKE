/*
Author: Silvio Trespi
 */ 

#include "silviosource\Classes.h"

#include "silviosource\Functions.h"

int ReactorNetwork_InputReader::ReadAndDraw() {
    
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
        Main_keyword.push_back("@SystemPressure");
        
        
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
        
        
        
        if(infile.fail() ) //Checking file opening
        {   
            cout << "Error while opening file\n" ;
            return 1;
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
        
        drawing.open(network_map_path,ios::out);
        
        if(drawing.fail()) //Checking drawing opening
        {   
            cout << "Error while opening NetworkMap_graphvizDot!\n" ;
            return 1;
        }
        
        drawing << "digraph ReactorNetwork { "<<endl<<"size = \"8,11\";"<<endl<<"edge [fontsize = 400, style = \"setlinewidth(11)\",arrowsize = 20];"<<endl;
        
        drawing << "node [fontsize = 400];"<<endl;
        drawing << "ratio = fill;"<<endl;
        drawing <<"ranksep = 12;"<<endl;
        drawing <<"nodesep = 1;"<<endl;

        drawing <<"INPUT[shape = polygon,sides = 5,peripheries = 3,color = lightblue,style = filled,height = 32,width = 32,rank = source];"<<endl;
        drawing <<"OUTPUT[shape = polygon,sides = 5,peripheries = 3,color = lightblue,style = filled,height = 32,width = 32,rank = sink];"<<endl;
        
            
        while( getline(infile,line) )  //Reading file

        {
            
            line_number++;

            

            Splitted_line = String_comma_separator(line);  //Separating line elements using spaces and commas
            
            if (Splitted_line.size() != 0)  //Checking that the line is not empty
                
            {

                if (Splitted_line.at(0) == Main_keyword.at(5)) // GETTING SYSTEM PRESSURE
                {
                    while (line.find("//") == string::npos)
                    {
                        getline(infile, line);
                        line_number++;

                        double number;
                        stringstream(Splitted_line.at(1)) >> number;
                        System_Pressure = number;
                        System_Pressure_unit_of_measurement = Splitted_line.at(2);
                    }
                }

                else if (Splitted_line.at(0) == Main_keyword.at(0))       //getting REACTOR information          

                {
                    vector <bool> R_keyword_controlling (Reactor_keyword.size(),false); // if a reactor keyword is found, then the related boolean in this vector turns true


                    if (LengthControl(Splitted_line,2,line_number) == 1) //CHECK: if reactor line does not have two elements (Reactor + Name), the program finishes
                    {
                        cout<<"This line should have two elements! ";
                        return 1;
                    }

    

                    if (FirstLetterThenNumberControl(Splitted_line.at(1),'R') =="No") //CHECK: if reactor name does not start with 'R', the program finishes
                    {
                        cout<<"Look at line "<<line_number<<".";
                        return 1;
                    }

                    (Reactor[reactor_count]).Name = Splitted_line.at(1);

                    Unit_name.push_back(Splitted_line.at(1));


                    for ( int i = 0; i < reactor_count ; i++) //CHECK: if two reactors have the same name, the program finishes
                    {
                        if (Reactor[reactor_count].Name == Reactor[i].Name)
                        {
                            cout<<"ERROR at line "<<line_number<<"!\nTwo reactors are identified by the same name "<<Reactor[reactor_count].Name<<"!";
                            return 1;
                        }
                    }

                    while ( line.find("//") == string::npos)  //This while loop gets lines until the separation element( // ) is found

                    {
                        getline(infile,line);
                        line_number++;
                        
                        vector <string> Splitted_line = String_comma_separator (line);
                        
                        
                        bool Unrecognized_keyword = true; // This boolean turns false if the first element of the line is one of the reactor keywords

                        if (Splitted_line.size() != 0) // Checking the line is not empty

                        {


                            for ( int j = 0; j < Reactor_keyword.size() ; j++)  

                            {   


                                if ( Splitted_line.at(0) == Reactor_keyword.at(j) )

                                {   

                                    Unrecognized_keyword = false;

                                    if ( R_keyword_controlling[j] == true )  // If the boolean related to eactor_keyword[j] has been already found, then the program finishes

                                    {
                                        cout<<"ERROR at line "<<line_number<<"!\n";
                                        cout<<"Reactor "<<Reactor[reactor_count].Name<<" has more than one keyword "<<Reactor_keyword[j]<<"!";
                                        return 1;
                                    }

                                    R_keyword_controlling[j] = true;



                                    switch (j)   

                                    { 
                                            case 0: //Type
                                            case 1: // Energy
                                            case 3: // Phase
                                                if ( LengthControl(Splitted_line,2,line_number) == 1)
                                                {
                                                    cout<<"This line should be \""<<Reactor_keyword[j] <<"\" followed by a word.";
                                                    return 1;
                                                }


                                                break;

                                            case 6: // Inlet_stream
                                            case 7: // Outlet_stream
                                                if ( LengthControl(Splitted_line,2,line_number) == 1)
                                                {
                                                    cout<<"This line should be \""<<Reactor_keyword[j]<<"\" (stream number)";
                                                    return 1;
                                                }
                                                
                                                if ( IntControl(Splitted_line.at(1)) == "No" )
                                                {
                                                    cout<<"ERROR at line "<<line_number<<"!\n";
                                                    cout<<"Stream number must be a positive integer!";
                                                    return 1;
                                                }



                                                break;

                                            case 2: // Residence time
                                            case 4: // Temperature
                                            case 5: //UA

                                                if ( LengthControl(Splitted_line,3,line_number) == 1)
                                                {
                                                    cout<<"This line should be "<<Reactor_keyword[j] <<" (number) (unit of measurement)";
                                                    return 1;
                                                }

                                                if ( DoubleControl(Splitted_line.at(1)) == "No")
                                                {
                                                    cout<<"Value of "<<Reactor_keyword[j]<<" is wrong!";
                                                    return 1;
                                                }
                                                if ( LetterControl(Splitted_line.at(2)) == "No")
                                                {
                                                    cout<<"Unit of measurement of "<<Reactor_keyword[j]<<" is wrong!";
                                                    return 1;
                                                }


                                    }

                                    if ( j == 0)    //Type

                                    {


                                        if((DictionaryControl(Splitted_line.at(1),Reactor_type_dictionary)) != 0)
                                        {
                                            cout<<"\nLine "<<line_number<<".";
                                            return 1;
                                        }

                                        Reactor[reactor_count].Type = Splitted_line.at(1);



                                    }

                                    else if ( j == 1)//Energy

                                    {

                                        if((DictionaryControl(Splitted_line.at(1),Reactor_energy_dictionary)) != 0)
                                        {
                                            cout<<"\nLine "<<line_number<<".";
                                            return 1;
                                        }

                                        Reactor[reactor_count].Energy = Splitted_line.at(1);


                                    }

                                    else if ( j == 2)//ResidenceTime

                                    {
                                        double number;
                                        stringstream(Splitted_line.at(1))>> number;
                                        Reactor[reactor_count].Residence_time = number;
                                        Reactor[reactor_count].Residence_time_unit_of_measurement = Splitted_line.at(2);

                                        if ( Reactor[reactor_count].Residence_time < 0)
                                        {
                                            cout<<"ERROR at line "<<line_number<<"!Residence time of reactor "<<Reactor[reactor_count].Name<<" is negative!\n";
                                            return 1;
                                        }
                                    }

                                    else if ( j == 3)//Phase

                                    {

                                        if((DictionaryControl(Splitted_line.at(1),Reactor_phase_dictionary)) != 0)
                                        {
                                            cout<<"\nline "<<line_number<<".";
                                            return 1;
                                        }

                                        Reactor[reactor_count].Phase = Splitted_line.at(1);


                                    }

                                    else if ( j == 4)//Temperature

                                    {
                                        double number;
                                        stringstream(Splitted_line.at(1))>> number;
                                        Reactor[reactor_count].Isothermal_temperature = number;
                                        Reactor[reactor_count].Isothermal_temperature_unit_of_measurement = Splitted_line.at(2);

                                        if ( Reactor[reactor_count].Isothermal_temperature < 0)
                                        {
                                            if (Reactor[reactor_count].Isothermal_temperature_unit_of_measurement == "K")
                                            {
                                                cout<<"ERROR at line "<<line_number<<"!Temperature of isothermal reactor "<<Reactor[reactor_count].Name<<" is under 0 K!";        
                                                return 1;
                                            }
                                        }

                                    }



                                    else if ( j == 5)//UA

                                    {
                                        double number;
                                        stringstream(Splitted_line.at(1))>> number;
                                        Reactor[reactor_count].UA = number;
                                        Reactor[reactor_count].UA_unit_of_measurement = Splitted_line.at(2);



                                    }

                                    else if ( j == 6)//Inlet_stream

                                    {
                                        int number;
                                        stringstream(Splitted_line.at(1)) >> number;

                                        Reactor[reactor_count].In_stream = number;
                                        Global_stream_number.push_back(number);

                                    }

                                    else if ( j == 7)//Outlet_stream

                                    {
                                        int number;
                                        stringstream(Splitted_line.at(1)) >> number;
                                            
                                        Reactor[reactor_count].Out_stream = number;
                                        Global_stream_number.push_back(number);

                                    }



                                }
                                
                            } // end of the for cycle


                            if ( Unrecognized_keyword == true && Splitted_line.at(0) != "//" ) // CHECK: if Unrecognized_keyword == true then The first element of a line 
                            {                                                                   // is not a Reactor_keyword and therefore the program finishes
                                
                                cout<<"ERROR at line "<<line_number<<"!\n";
                                
                                for (int i = 0; i < Main_keyword.size() ; i++)
                                {
                                    if (Splitted_line.at(0) == Main_keyword.at(i) )
                                    {
                                        cout<<"This is a new "<<Main_keyword.at(i)<<" declaration but the last one is not finished yet!\n";
                                        cout<<"Remember to end unit or stream declarations with \'//\'.";
                                        return 1;
                                    }
                                }
                                
                                cout<<"\""<<Splitted_line.at(0)<<"\" is not a Reactor keyword!";
                                return 1;
                            }




                            if (End_of_file( (  infile) )!= 0)
                            {

                                return 1;
                            }


                        }


                    }// end of while loop



                    for (int j = 0; j < R_keyword_controlling.size() ; j++)  // Checking the presence of all reactor keywords needed
                    {

                        if ( R_keyword_controlling[j] == false)

                        {
                            if ( j < 4 || j >= 6) // compulsory keywords
                            {
                                cout<<"ERROR at reactor "<<Reactor[reactor_count].Name<<"! "<<Reactor_keyword[j]<<" is missing!";
                                return 1;
                            }
                            else
                            {
                                if ( j == 4) //Temperature, only needed for isothermal reactors
                                {
                                    if ( Reactor[reactor_count].Energy =="Isothermal")
                                    {
                                        cout<<"ERROR at isothermal reactor "<<Reactor[reactor_count].Name<<"! "<<Reactor_keyword[j]<<" is missing!";
                                        return 1;
                                    }


                                }

        

                                else if ( j == 5) //UA, only needed for heatexchanger reactors
                                {
                                    if ( Reactor[reactor_count].Energy =="HeatExchanger")
                                    {
                                        cout<<"ERROR at heat-exchanger reactor "<<Reactor[reactor_count].Name<<"! "<<Reactor_keyword[j]<<" is missing!";
                                        return 1;
                                    }
    
                                }
                            }
                        }

                        else if (j== 4 && R_keyword_controlling[j] == true )  // if temperature is found but the reactor is not isothermal, then the program finishes

                        {
                            if (Reactor[reactor_count].Energy != "Isothermal")
                            {
                                cout<<"ERROR!"<<Reactor[reactor_count].Energy<<" Reactor "<<Reactor[reactor_count].Name<<" does not need a Temperature!\n";
                                return 1;
                            }
                        }

                        else if ( j == 5 && R_keyword_controlling[j] == true) // if UA is found but the reactor is not heatexchanger, then the program finishes
                        {


                                if (Reactor[reactor_count].Energy != "HeatExchanger")
                                {
                                    cout<<"ERROR!"<<Reactor[reactor_count].Energy<<" Reactor "<<Reactor[reactor_count].Name<<" does not need UA!\n";
                                    return 1;
                                }

                        }


                    }

                    if ( Reactor[reactor_count].Out_stream ==   Reactor[reactor_count].In_stream)  // if the input stream has the same name of the output stream, the the program finishes
                    {
                        cout<<"ERROR at reactor "<<Reactor[reactor_count].Name<<"!\n";
                        cout<<"Input and output stream have the same number \""<<Reactor[reactor_count].Out_stream<<"\".";
                        return 1;
                    }



                    drawing << Reactor[reactor_count].Name <<" [shape=";  // Graphviz instructions


                    if (Reactor[reactor_count].Type == "PSR")
                    {
                        drawing <<"circle,height =40,";
                    }
                    else if ( Reactor[reactor_count].Type == "PFR")
                    {
                        drawing <<"box,height = 48,width = 32,";
                    }

                    if ( Reactor[reactor_count].Energy == "Isothermal")
                    {
                        drawing <<"color = blue];"<<endl;
                    }
                    else if ( Reactor[reactor_count].Energy == "HeatExchanger")
                    {
                        drawing <<"color = red, style = dashed];"<<endl;
                    }
                    else if ( Reactor[reactor_count].Energy == "Adiabatic")
                    {
                        drawing <<"color = red];"<<endl;
                    }
                    reactor_count++;

                }//if(line.find("Reactor")




                else if ( Splitted_line.at(0) == Main_keyword.at(1))//Getting MIXER information

                {



                    if (LengthControl(Splitted_line,2,line_number) != 0)

                    {
                        return 1;
                    }

                    if (FirstLetterThenNumberControl(Splitted_line.at(1),'M') == "No")

                    {

                        cout<<"ERROR at line "<<line_number<<"!";
                        return 1;
                    }

                    Mixer[mixer_count].Name = Splitted_line.at(1);
                    Unit_name.push_back(Splitted_line.at(1));

                    for ( int i = 0; i < mixer_count; i++) // checking there are not two mixer with the same name
                    {
                        if ( Mixer[mixer_count].Name == Mixer[i].Name)
                        {
                            cout<<"ERROR at line "<<line_number<<"!\n";
                            cout<<"Two mixers are identified by the same name "<< Mixer[mixer_count].Name <<".";

                            return 1;
                        }
                    }

                    vector <bool> M_keyword_controlling (Mixer_keyword.size(),false);//every elements of this vector turns true if the mixer keyword with the same index is found

                    while ( line.find("//") == string::npos)
                    {

                        getline(infile,line);
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
                        for ( int j = 0; j < Mixer_keyword.size() ; j++) // for cycle: it is used to check that the first element of the line is a Mixer_keyword

                        {
                            if ( Splitted_line.at(0) == Mixer_keyword[j] )

                            {
                                
                                Unrecognized_keyword = false;
                                
                                if ( M_keyword_controlling[j] == true ) 

                                {
                                    cout<<"ERROR at line "<<line_number<<"!\n";
                                    cout<<"Mixer "<<Mixer[mixer_count].Name<<" has more than one keyword "<<Mixer_keyword[j]<<"!";
                                    return 1;
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

                                if ( j == 0 ) //Energy        
                                {

                                    if ( LengthControl(Splitted_line,2,line_number)!= 0)

                                        {
                                            cout<<"\nThis line should be :"<<Mixer_keyword[j]<<" (text string)";  //Energy Isothermal
                                            return 1;
                                        }


                                    if((DictionaryControl(Splitted_line.at(1),Mixer_energy_dictionary)) != 0) 
                                    {
                                        cout<<"\nLine "<<line_number<<".";
                                        return 1;
                                    }
                                    Mixer[mixer_count].Energy = Splitted_line.at(1);


                                }

                                else if ( j == 1 ) //Inlet_stream      
                                {

                                    if ( Splitted_line.size() < 2) // if the size is  < 2 , then  there is only one or no Inlet_stream and the program finishes

                                        {
                                            cout<<"ERROR at line "<<line_number<<endl;
                                            cout<<"Input streams of mixer "<<Mixer[mixer_count].Name<<" must be more than one!";
                                            return 1;
                                        }

                                    for ( int i = 1; i < Splitted_line.size() ; i++)

                                    { 

                                        if ( IntControl(Splitted_line.at(i)) == "No" )
                                        {
                                            cout<<"ERROR at line "<<line_number<<"!\n";
                                            cout<<"Stream number \""<<Splitted_line.at(i)<<"\" must be a positive integer!";
                                            return 1;
                                        }

                                        int number;
                                        stringstream(Splitted_line.at(i))>>number;

                                        Mixer[mixer_count].In_stream.push_back(number);
                                        Global_stream_number.push_back(number);
                                        
                                        
                                        if ( Mixer[mixer_count].InputStreamCheck() != 0)
                                        {
                                            cout<<"\nLook at line "<<line_number;
                                            return 1;
                                        }
                                    }


                                }

                                else if ( j == 2 )//Outlet_stream
                                {

                                        if ( LengthControl(Splitted_line,2,line_number)!= 0)

                                    {
                                        cout<<"Output stream of mixer "<<Mixer[mixer_count].Name<<" must be only one!";
                                        return 1;
                                    }

                                        if ( IntControl(Splitted_line.at(1)) == "No" )

                                    {
                                        cout<<"ERROR at line "<<line_number<<"!\n";
                                        cout<<"Stream number \""<<Splitted_line.at(1)<<"\" must be a positive integer!";
                                        return 1;
                                    }

                                    int number;
                                    stringstream(Splitted_line.at(1))>>number;

                                    Mixer[mixer_count].Out_stream = number;
                                    Global_stream_number.push_back(number);
                                }
                                
                                else if ( j == 3) // Temperature
                                {
                                    if ( LengthControl(Splitted_line,3,line_number)!= 0)

                                    {
                                        cout<<"This line should be: Temperature (value) (unit of measurement)";
                                        return 1;
                                    }
                                    
                                    if (DoubleControl(Splitted_line.at(1)) =="No" )
                                        {
                                            cout<<"ERROR at line "<<line_number<<"!\n";
                                            cout<<"Temperature value is incorrect!";
                                            return 1;
                                        }

                                        if ( LetterControl(Splitted_line.at(2) )== "No" )
                                        {
                                            cout<<"ERROR at line "<<line_number<<"!\n";
                                            cout<<"Temperature unit of measurement is incorrect!";
                                            return 1;
                                        }


                                        double Number;
                                        stringstream(Splitted_line.at(1))>>Number;

                                        Mixer[mixer_count].Temperature = Number;

                                        Mixer[mixer_count].Temperature_unit_of_measurement = Splitted_line.at(2);

                                        if ( Number < 0 && Splitted_line.at(2) =="K")

                                        {
                                        cout<<"ERROR at line "<<line_number<<"!\n";
                                            cout<<"Temperature in K cannot be negative!";
                                            return 1;
                                        }
                                }

                            }


                        } // end of for cycle
                        
                            if ( Unrecognized_keyword == true && Splitted_line.at(0) != "//" ) // CHECK: if Unrecognized_keyword == true then The first element of a line 
                            {                                                                   // is not a Mixer_keyword and therefore the program finishes
                                
                                cout<<"ERROR at line "<<line_number<<"!\n";
                                
                                for (int i = 0; i < Main_keyword.size() ; i++)
                                {
                                    if (Splitted_line.at(0) == Main_keyword.at(i) )
                                    {
                                        cout<<"This is a new "<<Main_keyword.at(i)<<" declaration but the last one is not finished yet!\n";
                                        cout<<"Remember to end unit or stream declarations with \'//\'.";
                                        return 1;
                                    }
                                }
                                
                                cout<<"\""<<Splitted_line.at(0)<<"\" is not a mixer keyword!";
                                return 1;
                            }

                            if (End_of_file( (  infile) )!= 0)
                            {

                                return 1;
                            }


                    }
                        
                    } // end of while loop

                    for (int i = 0; i < M_keyword_controlling.size() ; i++) // Checking the presence of every compulsory mixer keyword
                    {
                        if ( M_keyword_controlling[i] == false )
                        {
                            
                            if ( i != 3) // excluding Temperature
                            {
                            cout<<"ERROR at mixer "<<Mixer[mixer_count].Name<<"!\n";
                            cout<<Mixer_keyword[i] <<" is missing!";
                            return 1;
                            }
                            if ( i == 3 && Mixer[mixer_count].Energy == "Isothermal") //if Temperature is absent and mixer is Isothermal, the program finishes
                            {
                                cout<<"ERROR at isothermal mixer "<<Mixer[mixer_count].Name<<"!Temperature is required\n";
                                return 1;
                            }

                        }
                        
                        if ( M_keyword_controlling[3] == true && Mixer[mixer_count].Energy != "Isothermal" ) //if Temperature is present but mixer is not Isothermal, the program finishes
                        {
                            cout<<"ERROR! Mixer "<<Mixer[mixer_count].Name<<" is not isothermal and therefore does not need a Temperature!\n";
                            
                            return 1;
                        }
                        
                        
                    }

                    if ( Mixer[mixer_count].InputOutputStreamCheck()  != 0)
                    {
                        
                        return 1;
                    }
                    

                    drawing << Mixer[mixer_count].Name <<" [shape=invtriangle,";  // Graphviz instructions
    //                drawing << Mixer[mixer_count].Name <<" [shape=point,";  // Graphviz

                    if ( Mixer[mixer_count].Energy == "Isothermal")
                    {
                        drawing <<"color = blue];"<<endl;
                    }
                    else if (Mixer[mixer_count].Energy == "Adiabatic")
                    {
                        drawing <<"color = red];"<<endl;
                    }

                    mixer_count++;

                } //End of the MIXER block


                else if ( Splitted_line.at(0) == Main_keyword.at(2) )//Getting SPLITTER information

                {

                    vector <bool> S_keyword_controlling (Splitter_keyword.size(),false);

                    if (LengthControl(Splitted_line,2,line_number) != 0)

                    {
                        return 1;
                    }

                    if (FirstLetterThenNumberControl(Splitted_line.at(1),'S') == "No")

                    {

                        cout<<"ERROR at line "<<line_number;
                        return 1;
                    }

                    Splitter[splitter_count].Name = Splitted_line.at(1);
                    Unit_name.push_back(Splitted_line.at(1));

                    
                    
                    for ( int i = 0; i < splitter_count; i++) // Checking that two splitters do not have the same name
                    {
                        if ( Splitter[splitter_count].Name == Splitter[i].Name)
                        {
                            cout<<"ERROR at line "<<line_number<<"!\n";
                            cout<<"Two splitters are identified by the same name "<< Splitter[splitter_count].Name <<".";

                            return 1;
                        }
                    }

                    


                    while ( line.find("//") == string::npos)
                    {

                        getline(infile,line);
                        line_number++;
                        
                        Splitted_line = String_comma_separator(line);



    //                     if ( New_declaration(line,line_number) != 0)
    //                    {
    //                        return 1;
    //                    }
                        
                        if ( Splitted_line.size() != 0)
                            
                        {
                            bool Unrecognized_keyword = true;

                            
                            for ( int j = 0; j < Splitter_keyword.size() ; j++) //for cycle

                            {
                                if ( Splitted_line.at(0) == Splitter_keyword.at(j) )

                                {
                                    
                                    Unrecognized_keyword = false;
                                    
                                    if ( S_keyword_controlling.at(j) == true ) 

                                    {
                                        cout<<"ERROR at line "<<line_number<<"!\n";
                                        cout<<"Splitter "<<Splitter[splitter_count].Name<<" has more than one keyword "<<Splitter_keyword[j]<<"!";
                                        return 1;
                                    }

                                    S_keyword_controlling.at(j) = true;

                                    if ( j == 0 )//Inlet_stream     
                                    {

                                        if ( LengthControl(Splitted_line,2,line_number)!= 0)

                                        {
                                            cout<<"Input stream of splitter "<<Splitter[splitter_count].Name<<" must be only one!";
                                            return 1;
                                        }

                                        if ( IntControl(Splitted_line.at(1)) == "No" )

                                        {
                                            cout<<"ERROR at line "<<line_number<<"!\n";
                                            cout<<"Stream number \""<<Splitted_line.at(1)<<"\" must be a positive integer!";
                                            return 1;
                                        }

                                        int number;
                                        stringstream(Splitted_line.at(1)) >> number;


                                        Splitter[splitter_count].In_stream = number;
                                        Global_stream_number.push_back(number);

                                    }

                                    else if ( j == 1 )//Outlet_stream
                                    {

                                        Splitted_line.erase(Splitted_line.begin());

                                        double Summing_splitting_ratio = 0;

                                        int loop_number = 0;



                                        Splitter_Outlet_stream_loop: // this loop gets the stream number together with it splitting ratio

                                        loop_number++;

                                        if ( LengthControl(Splitted_line,2,line_number)!= 0)
                                        {
                                            cout<<"ERROR!Every line following Outlet_stream of splitter "<<Splitter[splitter_count].Name<<" should be (stream number),(splitting ratio)!\n";
                                            cout<<"End the list with a ;";
                                            return 1;
                                        }

                                        if ( IntControl( Splitted_line.at(0) ) == "No")
                                        {
                                            cout<<"Stream number must be an integer!Remember to end the list of output streams of splitter "<<Splitter[splitter_count].Name<<" with a ';'";
                                            cout<<"\nERROR at line "<<line_number;
                                            return 1;
                                        }

                                        if ( DoubleControl( Splitted_line.at(1) ) == "No")
                                        {
                                            cout<<"Splitting ratio must be a double!Remember to end the list of output streams of splitter "<<Splitter[splitter_count].Name<<" with a ';'\n";
                                            cout<<"\nERROR at line "<<line_number;
                                            return 1;
                                        }

                                        int int_number;
                                        stringstream(Splitted_line.at(0) )>>int_number;

                                        Splitter[splitter_count].Out_stream.push_back(int_number );
                                        Global_stream_number.push_back(int_number);
                                        
                                        if ( Splitter[splitter_count].OutputStreamCheck() != 0)
                                        {
                                            cout<<"\nLook at line "<<line_number;
                                            return 1;
                                        }



                                        double number;
                                        stringstream( Splitted_line.at(1) ) >> number;

                                        if ( number < 0 )

                                        {
                                            cout<<"ERROR at line "<<line_number;
                                            cout<<"\nSplitting ratio cannot be negative!";
                                            return 1;
                                        }

                                        Summing_splitting_ratio+=number;

                                        Splitter[splitter_count].Splitting_ratio.push_back(number);

            //                            if ( line.find("//") != string::npos)
            //                            {
            //                                cout<<"ERROR at splitter "<<Splitter[splitter_count].Name<<"!";
            //                                cout<<"List of output streams does not end with a ;!";
            //                                return 1;
            //                            }

                                        if ( line.find(";") == string::npos) //Loop  condition: if  this line doesn't contain any ';', then the loop restarts

                                        {
                                            getline(infile,line);
                                            line_number++;
                                            Splitted_line = String_comma_separator(line);
                                            goto Splitter_Outlet_stream_loop;

                                        }

                                        if (loop_number < 2) //Splitter requires at least two output streams, that means at least 2 loops
                                        {
                                            cout<<"ERROR at Splitter "<<Splitter[splitter_count].Name<<"!\n";
                                            cout<<"Slitter requires at least two output streams!";
                                            return 1;
                                        }


                                        if (abs(Summing_splitting_ratio-1) > 0.0000000001) //CHECK: splitting ratios must sum to 1
                                        {
                                            cout<<"ERROR at splitter "<<Splitter[splitter_count].Name<<"!\n";
                                            cout<<"Splitting ratios do not sum to 1!";
                                            return 1;
                                        }

                                    
                                            

                                            
                                        


                                    }


                                }


                            } // end of for cycle
                            
                            
                            
                            if ( Unrecognized_keyword == true && Splitted_line.at(0) != "//" ) // CHECK: if Unrecognized_keyword == true then The first element of a line 
                            {                                                                   // is not a Splitter_keyword and therefore the program finishes
                                
                                cout<<"ERROR at line "<<line_number<<"!\n";
                                
                                for (int i = 0; i < Main_keyword.size() ; i++)
                                {
                                    if (Splitted_line.at(0) == Main_keyword.at(i) )
                                    {
                                        cout<<"This is a new "<<Main_keyword.at(i)<<" declaration but the last one is not finished yet!\n";
                                        cout<<"Remember to end unit or stream declarations with \'//\'.";
                                        return 1;
                                    }
                                }
                                
                                cout<<"\""<<Splitted_line.at(0)<<"\" is not a splitter keyword!";
                                return 1;
                            }


                            if (End_of_file( (  infile) )!= 0)
                            {

                                return 1;
                            }

                        }
                        
                        
                    }//end of while loop


                    for (int i = 0; i < S_keyword_controlling.size() ; i++)
                    {
                        if ( S_keyword_controlling[i] == false)
                        {
                            cout<<"ERROR at splitter "<<Splitter[splitter_count].Name<<"!\n";
                            cout<<Splitter_keyword[i] <<" is missing!";
                            return 1;

                        }
                    }
                    
                    if ( Splitter[splitter_count].OutputInputStreamCheck() != 0)
                        
                    {
                        return 1;
                    }

                    drawing << Splitter[splitter_count].Name <<" [shape=triangle];"<<endl; // Graphviz
    //                 drawing << Splitter[splitter_count].Name <<" [shape=point];"<<endl; // Graphviz

                    splitter_count++;
                    
                }





                else if ( Splitted_line.at(0) == Main_keyword.at(3) ) //Getting PHASE-SPLITTER information

                {
                    
                    
                    vector <bool> PS_keyword_controlling (PhaseSplitter_keyword.size(),false);

                    
                    if (LengthControl(Splitted_line,2,line_number) != 0)

                    {
                        return 1;
                    }

                    if (FirstLetterThenNumberControl(Splitted_line.at(1),'P') == "No")

                    {

                        cout<<"ERROR at line "<<line_number;
                        return 1;
                    }

                    PhaseSplitter[phasesplitter_count].Name = Splitted_line.at(1);
                    Unit_name.push_back(Splitted_line.at(1));

                    
                

                    for ( int i = 0; i < phasesplitter_count; i++) // Checking two phase-splitter units do not have the same name

                    {
                        if ( PhaseSplitter[phasesplitter_count].Name == PhaseSplitter[i].Name)
                        {
                            cout<<"ERROR at line "<<line_number<<"!\n";
                            cout<<"Two Phase-splitters are identified by the same name: "<< PhaseSplitter[phasesplitter_count].Name<<".";

                            return 1;
                        }
                    }

                    


                    while ( line.find("//") == string::npos)
                        
                    {

                        getline(infile,line);
                        line_number++;
                        
                        Splitted_line = String_comma_separator(line);

                        
                        
                        if ( Splitted_line.size() != 0)
                            
                        {
                            
                        
                            bool Unrecognized_keyword = true;


                            for ( int j = 0; j < PhaseSplitter_keyword.size() ; j++)

                            {
                                if ( Splitted_line.at(0) == PhaseSplitter_keyword.at(j)  )

                                {
                                    
                                    Unrecognized_keyword = false;

                                    if ( PS_keyword_controlling.at(j) == true ) 

                                    {
                                        cout<<"ERROR at line "<<line_number<<"!\n";
                                        cout<<"Phase-splitter "<<PhaseSplitter[phasesplitter_count].Name<<" has more than one keyword "<<PhaseSplitter_keyword[j]<<"!";
                                        return 1;
                                    }
                                    PS_keyword_controlling[j] = true;



                                    if ( j == 0 )//Inlet_stream     
                                    {

                                        if ( LengthControl(Splitted_line,2,line_number)!= 0)

                                        {
                                            cout<<"Input stream of Phase-Splitter "<<PhaseSplitter[phasesplitter_count].Name<<" must be only one!";
                                            return 1;
                                        }

                                        if ( IntControl(Splitted_line.at(1)) == "No" )

                                        {
                                            cout<<"ERROR at line "<<line_number<<"!\n";
                                            cout<<"Stream number \""<<Splitted_line.at(1)<<"\" must be a positive integer!";
                                            return 1;
                                        }

                                        int number;
                                        stringstream(Splitted_line.at(1)) >> number;

                                        PhaseSplitter[phasesplitter_count].In_stream = number;
                                        Global_stream_number.push_back(number);

                                    }

                                    else if ( j == 1 )//Outlet_stream
                                        
                                        
                                    {

                                        Splitted_line.erase(Splitted_line.begin());


                                        int loop_number = 0;

                                        PhaseSplitter_Outlet_stream_loop: // This loop gets the stream name and the phase of this stream

                                        loop_number++;

                                        if ( LengthControl(Splitted_line,2,line_number)!= 0)

                                        {

                                            cout<<"ERROR!Every line following Outlet_stream of Phase-Splitter "<< PhaseSplitter[phasesplitter_count].Name<<" should be (stream number),(phase)!\n";
                                            cout<<"End the list with a ;";
                                            return 1;
                                        }

                                        if ( IntControl( Splitted_line.at(0) ) == "No")

                                        {
                                            cout<<"Stream number must be an integer!Remember to end the list of output streams of Phase-Splitter "<< PhaseSplitter[phasesplitter_count].Name<<" with a ';'";
                                            cout<<"\nERROR at line "<<line_number;
                                            return 1;
                                        }

                                        int number;
                                        stringstream(Splitted_line.at(0))>>number;

                                        PhaseSplitter[phasesplitter_count].Out_stream.push_back( number );
                                        Global_stream_number.push_back(number);



                                        if((DictionaryControl(Splitted_line.at(1),Splitted_phase_dictionary)) != 0)
                                        {
                                            cout<<"\nLine "<<line_number<<".";
                                            return 1;
                                        }

                                        PhaseSplitter[phasesplitter_count].Splitted_phase.push_back(Splitted_line.at(1));

                                        if ( PhaseSplitter[phasesplitter_count].OutputStreamCheck() != 0)
                                        {
                                            cout<<"\nLook at line "<<line_number;
                                            return 1;
                                        }


                                        if ( line.find(";") == string::npos)

                                        {
                                            getline(infile,line);
                                            line_number++;
                                            Splitted_line = String_comma_separator(line);

                                            goto PhaseSplitter_Outlet_stream_loop;

                                        }


                                        if (loop_number == 1)
                                        {
                                            cout<<"ERROR at Phase-Splitter "<<PhaseSplitter[phasesplitter_count].Name<<"!\n";
                                            cout<<"Phase-Slitter requires at least two output streams!";
                                            return 1;
                                        }


                                    }


                                }


                            } // end of for cycle
                            
                            if ( Unrecognized_keyword == true && Splitted_line.at(0) != "//" ) // CHECK: if Unrecognized_keyword == true then The first element of a line 
                            {                                                                   // is not a Phase-splitter_keyword and therefore the program finishes
                                
                                cout<<"ERROR at line "<<line_number<<"!\n";
                                
                                for (int i = 0; i < Main_keyword.size() ; i++)
                                {
                                    if (Splitted_line.at(0) == Main_keyword.at(i) )
                                    {
                                        cout<<"This is a new "<<Main_keyword.at(i)<<" declaration but the last one is not finished yet!\n";
                                        cout<<"Remember to end unit or stream declarations with \'//\'.";
                                        return 1;
                                    }
                                }
                                
                                cout<<"\""<<Splitted_line.at(0)<<"\" is not a Phase-splitter keyword!";
                                return 1;
                            }


                            if (End_of_file( (  infile) )!= 0)
                            {

                                return 1;
                            }

                        }
                        
                    }// end of while loop

                    for (int i = 0; i < PS_keyword_controlling.size() ; i++)

                    {

                        if ( PS_keyword_controlling[i] == false)

                        {

                            cout<<"ERROR at Phase-Splitter "<<PhaseSplitter[phasesplitter_count].Name<<"!\n";
                            cout<<PhaseSplitter_keyword[i] <<" is missing!";
                            return 1;

                        }
                    }
                    
                    if ( PhaseSplitter[phasesplitter_count].OutputInputStreamCheck() != 0)
                    {
                        return 1;
                    }

                    drawing << PhaseSplitter[phasesplitter_count].Name <<" [shape=triangle,color = midnightblue];"<<endl;

                    phasesplitter_count++;

                } //end of PHASE-SPLITTER block

                
                
                
                
                else if ( Splitted_line.at(0) == Main_keyword.at(4) ) //Getting STREAM information

                {
                    
                    if ( LengthControl(Splitted_line,2,line_number) != 0)

                    {

                        return 1;
                    }

                    if ( IntControl(Splitted_line.at(1))== "No")

                    {
                        cout<<"ERROR at line "<<line_number<<"!\n";
                        cout<<"Stream number must be an integer!";
                        return 1;

                    }
                    
                    int int_number;
                    stringstream(Splitted_line.at(1))>> int_number;
                    Stream[stream_count].Name = int_number;
                    Stream[stream_count].GivingSizeToKeyword_presence( Stream_keyword.size() );
                    
                    for ( int i = 0; i < stream_count; i++) // Checking two streams do not have the same number

                    {
                        if ( Stream[stream_count].Name == Stream[i].Name)
                            
                        {
                            cout<<"ERROR at line "<<line_number<<"!\n";
                            cout<<"Two Streams are identified by the same number: "<< Stream[stream_count].Name<<".";

                            return 1;
                        }
                    }


                    while (line.find("//")== string::npos)

                    {
                        getline(infile,line);
                        line_number++;

                        

                        
                        Splitted_line = String_comma_separator(line);

                        if ( Splitted_line.size() != 0)
                            
                        {
                            
                            bool Unrecognized_keyword = true;
                        


                            for ( int j = 0; j < Stream_keyword.size(); j++)

                            {

                                if ( Splitted_line.at(0) == Stream_keyword.at(j) )
                                    
                                {
                                    
                                    Unrecognized_keyword = false;
                                    

                                    
                                    if(Stream[stream_count].TrueKeyword(j) != 0) //if a keyword is repeated twice, the program finishes
                                    {
                                        cout<<"ERROR at line "<<line_number<<"!\n";
                                        cout<<"Stream "<<Stream[stream_count].Name<<" keyword \""<<Stream_keyword.at(j)<<"\" is repeated twice!";
                                        return 1;
                                        
                                    }



                                    if ( j == 0)//Phase
                                        
                                    {
                                        if ( LengthControl(Splitted_line,2,line_number) != 0)

                                        {
                                            return 1;
                                        }

                                        if((DictionaryControl(Splitted_line.at(1),Stream_phase_dictionary)) != 0)
                                            
                                        {
                                            cout<<"\nLine "<<line_number<<".";
                                            return 1;
                                        }

                                        Stream[stream_count].Phase = Splitted_line.at(1);

                                    }

                                    if ( j == 1 )//SolidType
                                        
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
                                        
                                        if ( LengthControl(Splitted_line,2,line_number)!= 0)
                                        {

                                            cout<<"This line should be "<<Stream_keyword.at(j)<<" (text string)!\n";
                                            
                                            return 1;
                                        }


    //                                    if((DictionaryControl(Splitted_line.at(0),Stream_solid_type_dictionary)) != 0)
    //                                    {
    //                                        cout<<"\nLine "<<line_number<<".";
    //                                        return 1;
    //                                    }
                                        
                                        if((DictionaryControl(Splitted_line.at(1),Stream_solid_type_dictionary)) != 0)
                                        {
                                            cout<<"\nLine "<<line_number<<".";
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

                                    if ( j == 2)//MassFlowRate

                                    {

                                        if ( LengthControl(Splitted_line,3,line_number)!= 0)
                                        {

                                            cout<<"\nThis line should contain (double) (text string)";

                                            return 1;
                                        }

                                        if (DoubleControl(Splitted_line.at(1)) =="No" )
                                        {
                                            cout<<"ERROR at line "<<line_number<<"!\n";
                                            cout<<"MassFlowRate value is incorrect!";
                                            return 1;
                                        }

                                        if ( LetterControl(Splitted_line.at(2) )== "No" )
                                        {
                                            cout<<"ERROR at line "<<line_number<<"!\n";
                                            cout<<"MassFlowRate unit of measurement is incorrect!";
                                            return 1;
                                        }

                                        double Number;
                                        stringstream(Splitted_line.at(1))>>Number;
                                        if ( Number < 0 )
                                        {
                                            cout<<"ERROR at line "<<line_number<<"!\n";
                                            cout<<"MassFlowRate cannot be negative!";
                                            return 1;
                                        }
                                        Stream[stream_count].MassFlowRate = Number;
                                        Stream[stream_count].MassFlowRate_unit_of_measurement = Splitted_line.at(2);

                                    }

                                    if ( j == 3)//Temperature

                                    {

                                        if ( LengthControl(Splitted_line,3,line_number)!= 0)
                                        {

                                            cout<<"\nThis line should be Temperature (number) (unit of measurement)";

                                            return 1;
                                        }

                                        if (DoubleControl(Splitted_line.at(1)) =="No" )
                                        {
                                            cout<<"ERROR at line "<<line_number<<"!\n";
                                            cout<<"Temperature value is incorrect!";
                                            return 1;
                                        }

                                        if ( LetterControl(Splitted_line.at(2) )== "No" )
                                        {
                                            cout<<"ERROR at line "<<line_number<<"!\n";
                                            cout<<"Temperature unit of measurement is incorrect!";
                                            return 1;
                                        }


                                        double Number;
                                        stringstream(Splitted_line.at(1))>>Number;

                                        Stream[stream_count].Temperature = Number;

                                        Stream[stream_count].Temperature_unit_of_measurement = Splitted_line.at(2);

                                        if ( Number < 0 && Splitted_line.at(2) =="K")

                                        {
                                        cout<<"ERROR at line "<<line_number<<"!\n";
                                            cout<<"Temperature in K cannot be negative!";
                                            return 1;
                                        }


                                    }


                                    if ( j == 4)//MassFraction
                                    {

                                        Splitted_line.erase(Splitted_line.begin());
                                        
                                        

                                        Mass_fraction_loop:
                                                
                                                

                                        bool Controlling_species = false;

                                        if ( LengthControl(Splitted_line,2,line_number)!= 0)

                                        {

                                            cout<<"Every line following MassFraction of stream "<< Stream[stream_count].Name<<" should be (Compound) (mass fraction)!\n";
                                            cout<<"End the list with a ;";
                                            return 1;
                                        }
                                        
                                        


                                        Stream[stream_count].Compound.push_back(Splitted_line.at(0));

                                        if (Stream[stream_count].CompoundControl() != 0)
                                        {
                                            return 1;
                                        }



                                        for ( int i = 0; i < Species.size() ; i++) //Checking if the new compound has already been found

                                        {
                                            if ( ( Splitted_line.at(0) ).compare( Species.at(i) ) == 0 )

                                            {
                                                Controlling_species = true; 

                                                break;
                                            }

                                        }

                                        if ( Controlling_species == false) // if not, then it is added to Species vector

                                        {
                                        Species.push_back( Splitted_line.at(0) );
                                        }

                                        
                                        if (DoubleControl(Splitted_line.at(1) ) == "No" ) //Checking Mass fraction value is a double
                                        {
                                            cout<<"ERROR at line "<<line_number<<"!\n";
                                            cout<<"Mass Fraction value of compound \""<<Stream[stream_count].Compound.back()<<"\" must be a double.";
                                            return 1;
                                        }


                                        double number;
                                        stringstream(Splitted_line.at(1))>> number;
                                        
                                        if( number < 0 )//Checking Mass fraction value is positive
                                        {
                                            cout<<"ERROR at line "<<line_number<<"!\n";
                                            cout<<"Mass Fraction value of compound \""<<Stream[stream_count].Compound.back()<<"\" must be a positive double.";
                                            return 1;
                                        }
                                        Stream[stream_count].Composition.push_back(number);



    //                                    if (Stream[stream_count].SummingMassFractions() != 0)
    //                                    {
    //                                        return 1;
    //                                    }


                                        if ( line.find(";") == string::npos)

                                        {
                                            getline(infile,line);
                                            line_number++;
                                            Splitted_line = String_comma_separator(line);

                                            goto Mass_fraction_loop;

                                        }
                                    }





                                }

                            } //end of for cycle
                            
                            
                            
                            if ( Unrecognized_keyword == true && Splitted_line.at(0) != "//" ) // CHECK: if Unrecognized_keyword == true then The first element of a line 
                            {                                                                   // is not a Stream_keyword and therefore the program finishes
                                
                                cout<<"ERROR at line "<<line_number<<"!\n";
                                
                                for (int i = 0; i < Main_keyword.size() ; i++)
                                {
                                    if (Splitted_line.at(0) == Main_keyword.at(i) )
                                    {
                                        cout<<"This is a new "<<Main_keyword.at(i)<<" declaration but the last one is not finished yet!\n";
                                        cout<<"Remember to end unit or stream declarations with \'//\'.";
                                        return 1;
                                    }
                                }
                                
                                cout<<"\""<<Splitted_line.at(0)<<"\" is not a Stream keyword!";
                                return 1;
                            }


                            if (End_of_file( (  infile) )!= 0)
                            {

                                return 1;
                            }



                        }


                    }//end of while loop

                    
                    if ( Stream[stream_count].MassFractionSum() > 1)
                    {
                        cout<<"ERROR! Mass fraction sum of stream "<<Stream[stream_count].Name<<" is more than 1!";

                        return 1;
                    }
                    
                    stream_count++;  

                } // End of STREAM block

                
                
                else // if Splitted_line.at(0) is not Reactor,Mixer,Splitter,Phase-splitter,Stream
                    
                {


                cout<<"ERROR at line "<<line_number<<"!\n";
                cout<<"Main Keyword \""<<Splitted_line.at(0)<<"\" not recognized!";
                return 1;


                }
                
                

            }// if splitted_line.size() != 0
                
        } //end of the INPUT FILE reading
        
        
        
        infile.close(); //Closing input_file
        
        Splitted_line.clear(); 
        

        
        cout<<endl<<endl;
        
        
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
        
        
        vector <int> Stream_number;
    

        for (int i = 0; i < Global_stream_number.size() ; i++) //Eliminating streams that are repeated twice from Global_stream_number
        {
            
            bool Check = false;
            
            
                for ( int j = 0; j <  Stream_number.size();j++)
                    
                {
                    
                    if ( Stream_number.at(j) ==  Global_stream_number.at(i) )
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
        
        
        
        
        
        for ( int i = 0; i < stream_count; i++) //Checking the validity of every phase splitter
        {
            for ( int j = 0; j < phasesplitter_count; j++)
            {
                if ( Stream[i].Name == PhaseSplitter[j].In_stream)
                {
                    if (Stream[i].Phase != "Mix")
                    {
                        cout<<"ERROR! Stream "<<PhaseSplitter[j].In_stream<<" enters Phase-splitter "<<PhaseSplitter[j].Name<<" and therefore its Phase must be \"Mix\".";
                        return 1;
                    }
                }
                
                for ( int x = 0; x < PhaseSplitter[j].Out_stream.size() ; x++)
                {
                    if (Stream[i].Name == PhaseSplitter[j].Out_stream.at(x) )
                    {
                        if (Stream[i].Phase != PhaseSplitter[j].Splitted_phase.at(x) )
                        {
                            cout<<"ERROR! Stream "<<Stream[i].Name<<" phase is "<<Stream[i].Phase<<" but Phase-splitter "<<PhaseSplitter[j].Name<<" declares it is "<<PhaseSplitter[j].Splitted_phase.at(x)<<"!";
                            return 1;
                        }
                    }
                }
            }
            
        
                
        }
    
        
        
        
        //4 vectors related to stream
        
        vector <string> From(Stream_number.size(),"0");     //Unit FROM which the stream arrives
        vector <string> To(Stream_number.size(),"0");       //Unit TO which the stream arrives
        vector <double> Splitting_ratio(Stream_number.size(),1); 
        vector <string> Stream_phase(Stream_number.size(),"Mix"); //Phase of every system stream, by default it is Mix
        
        
        
        
        
        
        
        for ( int i = 0; i < Stream_number.size() ; i++)   //Filling the 4 vectors taking information from the the vector of classes previously created 
            
        {
            
            for ( int x = 0; x < splitter_count; x++)
            {
                for (int j = 0; j < Splitter[x].Out_stream.size();j++)
                {
                    if (Stream_number[i] == Splitter[x].Out_stream[j])
                    {
                        if (From[i] != "0")
                        {
                            cout<<"ERROR! Units "<<Splitter[x].Name<<" and "<<From[i]<<" have the same outlet stream: "<<Stream_number[i]<<" !\n";
                            return 1;
                        }
                        
                        Splitting_ratio[i] = Splitter[x].Splitting_ratio[j];
                        
                        From[i] = Splitter[x].Name;
                    
                    }
                
                }
                
                if (Stream_number[i] == Splitter[x].In_stream)
                {
                    if (To[i] != "0")
                    {
                        cout<<"ERROR! Units "<<Splitter[x].Name<<" and "<<To[i]<<" have the same inlet stream: "<<Stream_number[i]<<" !\n";
                        return 1;
                    }
                    
                    To[i] = Splitter[x].Name;
                
                }
                
                
            }
            
            for ( int x = 0; x < reactor_count ; x++)
            {
                if (Stream_number[i] == Reactor[x].Out_stream)
                {
                    if (From[i] != "0")
                    {
                        cout<<"ERROR! Units "<<Reactor[x].Name<<" and "<<From[i]<<" have the same output stream: "<<Stream_number[i]<<" !\n";
                        return 1;
                    }
                    From[i] = Reactor[x].Name;
                    
                }
                
                
                if (Stream_number[i] == Reactor[x].In_stream)
                {
                    if (To[i] != "0")
                    {
                        cout<<"ERROR! Units "<<Reactor[x].Name<<" and "<<To[i]<<" have the same input stream: "<<Stream_number[i]<<" !\n";
                        return 1;
                    }
                    To[i] = Reactor[x].Name;
                    
                }
                
            }
            

            
            for ( int x = 0; x < mixer_count ; x++)
                
            {
                if (Stream_number[i] == Mixer[x].Out_stream)
                    
                {
                    if (From[i] != "0")
                        
                    {
                        cout<<"ERROR!Units "<<Mixer[x].Name<<" and "<<From[i]<<" have the same output stream: "<<Stream_number[i]<<" !\n";
                        return 1;
                    }
                    
                    From[i] = Mixer[x].Name;
                
                }
                
                
                
                for ( int j = 0; j < Mixer[x].In_stream.size() ; j++)
                    
                {
                    if (Stream_number[i] == Mixer[x].In_stream[j])
                        
                    {
                        if (To[i] != "0")
                            
                        {
                            cout<<"ERROR! Units "<<Mixer[x].Name<<" and "<<To[i]<<" have the same input stream: "<<Stream_number[i]<<" !\n";
                            return 1;
                        }
                        
                        To[i] = Mixer[x].Name;
                    
                    }
                    
                }
            }
            
            for ( int x = 0; x < phasesplitter_count ; x++)
                
            {
                if (Stream_number[i] == PhaseSplitter[x].In_stream)
                    
                {
                    if (To[i] != "0")
                    {
                        cout<<"ERROR! Units "<<PhaseSplitter[x].Name<<" and "<<To[i]<<" have the same input stream: "<<Stream_number[i]<<" !\n";
                        return 1;
                    }
                    To[i] = PhaseSplitter[x].Name;
                    
                }
                
                for ( int j = 0; j < PhaseSplitter[x].Out_stream.size() ; j++)
                    
                {
                    if (Stream_number[i] == PhaseSplitter[x].Out_stream[j])
                    {
                        
                        if (From[i] != "0")
                        {
                            cout<<"ERROR! Units "<<PhaseSplitter[x].Name<<" and "<<From[i]<<" have the same output stream: "<<Stream_number[i]<<" !\n";
                            return 1;
                        }
                        From[i] = PhaseSplitter[x].Name;
                        
                        Stream_phase.at(i) = PhaseSplitter[x].Splitted_phase[j];
                        
                        
                    }
                    
                }
            }
            
            for ( int x = 0; x < stream_count ; x++)
                
            {
                
                if ( Stream_number.at(i) == Stream[x].Name)
                    
                {
                    
                    Stream_phase.at(i) = Stream[x].Phase;
                
                }
            }
        }
        
        vector <bool> Controlling_input_output(2,false); // The first element of this vector will turn true if there is a system input stream
                                                        //The second element will turn true if there is a system output stream
        
        
        
        
        for (int i = 0;  i < Stream_number.size() ; i++)
            
        {
            
            if (From.at(i) == "0") // if stream at index i is a system input, then it must have a declaration
                
            {
                Controlling_input_output.at(0) = true; 
                
                bool Input_stream_check = false;
                
                for (int x = 0; x < stream_count; x++) // checking over every stream declaration
                {
                    if (Stream_number.at(i) == Stream[x].Name )
                    {
                        Input_stream_check = true;
                        break;
                    }
                    
                    
                }
                
                if ( Input_stream_check == false ) // if this stream has not been found, then the program finishes
                {
                    cout<<"ERROR! System input stream "<<Stream_number.at(i)<<" does not have a declaration!\n";
                    return 1;
                }
                
                
            }
            
            if (To.at(i) == "0")
                
            {
                Controlling_input_output.at(1) = true;

            }
            
        }
        
        
        
        
        
        
        
        for ( int i = 0; i < 2 ; i++)
        {
            if ( Controlling_input_output[i] == false)
            {
                if ( i == 0)
                {
                    cout<<"ERROR! There is no system input!\n";
                    return 1;
                }
                
                if ( i == 1)
                {
                    cout<<"ERROR! There is no system output!\n";
                    return 1;
                }
            }
        }
    
        for ( int i = 0; i < Stream_number.size() ; i++) //Checking the validity of the stream keywords found
        {
            
            for ( int j = 0; j < stream_count; j++)
                
            {
                
                
                if ( Stream_number.at(i) == Stream[j].Name)
                    
                {
                    if (  From.at(i) == "0")  //Considering system input streams!
                        
                    {
                        for ( int x = 1; x < Stream_keyword.size() ; x++)   //Starting from x = 1 since Phase is not considered. the default phase is Mix

                        {
                            if (Stream[j].Keyword_presence.at(x) == false)

                            {

                                if ( x == 1) //if there is not type but the stream phase is Solid, then the program finishes
                                    
                                {

                                    if ( Stream[j].Phase == "Solid")
                                    {
                                        cout<<"ERROR!Solid Stream "<<Stream[j].Name<<" is a system input and therefore needs: "<< Stream_keyword.at(x)<<".";
                                        return 1;

                                    }
                                }

                                else //Every other stream keyword must be present in the declaration of a system input stream

                                {            
                                        cout<<"ERROR!Stream "<<Stream[j].Name<<" is a system input and therefore needs: "<<Stream_keyword.at(x)<<".";
                                        return 1;

                                }    

                                

                            }



                        }
                        
                        
                        if ( Stream[j].Keyword_presence.at(1) == true && Stream[j].Phase != "Solid") //if there is a type but the Stream is not solid, then the program finishes
                        {
                            cout<<"ERROR!Stream "<<Stream[j].Name<<" is "<<Stream[j].Phase<<" and therefore does not need: "<<Stream_keyword.at(1)<<".";
                            return 1;
                        }

                        if (Stream[j].Phase == "Solid")//Checking that a solid system input stream mass fraction is expressed in terms of elements

                        {

                            for ( int x = 0; x < Stream[j].Compound.size() ; x++) 

                            {


                                if((DictionaryControl(Stream[j].Compound.at(x),Solid_compound_dictionary)) != 0)//This check is for system input streams only
                                {
                                    cout<<"Change the related MassFraction line of SOLID stream "<<Stream[j].Name<<".";
                                    return 1;
                                }
                            }
                        }

                        if (abs(Stream[j].MassFractionSum()-1) > 0.00000000001)  //Checking the MassFractionSum is 1, only for system input streams

                        {

                            cout<<"ERROR!Stream "<<Stream[j].Name<<" is a system input and therefore its mass-fraction-sum must be 1!\n";
                            return 1;
                        }
                    }
                    
                    else // Considering non-input streams!
                    {
                        if (Stream[j].Keyword_presence.at(1) == true)
                        {
                            cout<<"ERROR!Stream "<<Stream[j].Name<<" is  not a system input stream, therefore does not need a SolidType!";
                            return 1;
                            
                        }
                    }
                }
                
                
                
            }
        }
        
        
        
        
        
        for ( int j = 0; j < stream_count; j++) // Checking that every declarated stream is part of the system
            
        {
            Stream[j].ShrinkingKeyword_presence();
            
            bool Control = false;
            
            for ( int i = 0; i < Stream_number.size() ; i++)
            {
            
                if ( Stream[j].Name == Stream_number.at(i) )
                {
                    Control = true;
                    break;
                }
            
            }
            
            if ( Control == false)
            {
                cout<<"ERROR!Stream "<<Stream[j].Name<<" has not been defined in the system!";
                return 1;
            }
            
        }
        
        
        for ( int i = 0;i < splitter_count; i++) // Checking the validity of every splitter: it must be mono-phase
        {
            string In_stream_phase;
            
            vector <string> Out_stream_phase;
            vector <int> remembering_out_stream;
            
            for( int j = 0; j < Stream_number.size(); j++)
                
            {
                if(Stream_number.at(j) == Splitter[i].In_stream)
                    
                {
                    In_stream_phase = Stream_phase.at(j);
                }
                
                else
                {
                    for( int x = 0; x < Splitter[i].Out_stream.size(); x++)
                    {
                        if ( Stream_number.at(j) == Splitter[i].Out_stream.at(x) )
                        {
                            Out_stream_phase.push_back( Stream_phase.at(j) );
                            
                            remembering_out_stream.push_back(Stream_number.at(j));
                        }
                    }
                }
            }
            
            
            
            for ( int w = 0; w < Out_stream_phase.size() ; w++)
            {
                
                
                
                if (In_stream_phase != Out_stream_phase.at(w))
                {
                    cout<<"ERROR!Splitters involve always only one phase!\n";
                    cout<<"Splitter "<<Splitter[i].Name<<" has: \""<<In_stream_phase<<"\" stream "<<Splitter[i].In_stream<<" and \""<<Out_stream_phase.at(w)<<"\" stream "<<remembering_out_stream.at(w)<<"!";
                    return 1;
                }
            }
            
            
        }
        
        

        
        
        
        
        for ( int i = 0; i < Stream_number.size(); i++)  //GraphViz instructions
        {
            
            if (From.at(i) == "0")
            {
                
                
                drawing <<"INPUT -> "<<To.at(i);
                
            }
            else if (To.at(i) == "0")
            {
                drawing <<From.at(i)<<" -> OUTPUT";
            }
            else
            {
                drawing <<From.at(i)<<" -> "<<To.at(i);
            }
            
            drawing <<"[ ";

            
            if ( Stream_phase.at(i) != "Mix")
            {
                if ( Stream_phase.at(i) == "Liquid")
                {
                    drawing <<"color = lightblue,";
                }
                
                else if ( Stream_phase.at(i) == "Gas")
                {
                    drawing <<"color = orange,";
                }
                
                else if ( Stream_phase.at(i) == "Solid")
                {
                    drawing <<"color = green,";
                }
            }
            
            drawing <<" label =\" "<<Stream_number.at(i);
            
            
        
                
                
            if ( *( From.at(i).begin() ) == 'S')
            {
                
                drawing <<", "<<Splitting_ratio.at(i)<<" \" ];"<<endl;
            }
            else
            {
                drawing <<" \" ];"<<endl;
            }
                
            
            
            
            
            
        }
        
        
        
        drawing << "}"<<endl;
        
        drawing.close();
        
        // Generate PDF
        system("dot -Tpdf NetworkMap_graphvizDot.dot -o NetworkMap.pdf ");
    
    

    return 0;
    
}







