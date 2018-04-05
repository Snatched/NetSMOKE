


/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Classes.h
 * Author: silvio14
 *
 * Created on June 27, 2017, 4:03 AM
 */

#ifndef CLASSES_H
#define CLASSES_H


#include <vector>

#include <iostream>

#include <string>

using namespace std;


class Reactor_unit

{
public:
    
    Reactor_unit();
    
    string Name;
    string Type;
    string Phase;
    string Energy;
    double Residence_time;
    string Residence_time_unit_of_measurement;
    double Volume;
    string Volume_unit_of_measurement;
    double Diameter;
    string Diameter_unit_of_measurement;
    double Length;
    string Length_unit_of_measurement;
    double Isothermal_temperature;
    string Isothermal_temperature_unit_of_measurement;

    
    double UA;
    string UA_unit_of_measurement;
    
    int In_stream;
    int Out_stream;
    
    
    
};




class Mixer_unit

{
  
public:
       
    
    string Name;
    int Out_stream;
    vector <int> In_stream;
    string Energy;
    double Temperature;
    string Temperature_unit_of_measurement;
   
    int InputStreamCheck(); // This function checks the last input stream taken has a different number then the other  previous input streams
    int InputOutputStreamCheck(); // This function checks that the output stream has a different number then every input stream
   
    
};  






class S_unit

{
    
public:
    
    string Name;
    int In_stream;
    vector <int> Out_stream;
    
    
    vector <double> Splitting_ratio;
    
    int OutputStreamCheck(); // This function checks the last output stream taken has a different number then the other  previous output streams
    int OutputInputStreamCheck(); // This function checks that the output stream has a different number then every input stream
    
};





class PS_unit

{
    
public:
    
    string Name;
    int In_stream;
    vector <int> Out_stream;
    
    vector <string> Splitted_phase;
    
    int OutputStreamCheck();
    int OutputInputStreamCheck();
    

};


class Stream_class

{
    
public:
    
    
    Stream_class();

   
    int CompoundControl(); // Checks that a compound is not declarated twice in the same stream
   
    double MassFractionSum();// returns the mass fraction sum
    
    int TrueKeyword(int); //if Keyword_presence at index j is not already true, Keyword_presence at index j becomes true and returns 0
    
    void ShrinkingKeyword_presence();// shrinking(minimizing occupied memory) the Keyword_presence vector
    
    void GivingSizeToKeyword_presence(int);//Keyword_presence is filled with the exact number of "false", equal to the size of stream_keyword
    

    
    int Name;
    string Phase;

    string SolidType;
    
    double MassFlowRate;
    string MassFlowRate_unit_of_measurement;
    double Temperature;
    string Temperature_unit_of_measurement;
    
    vector <string> Compound;
    vector <double> Composition;
	vector<double> omega_teo;
    
    
    vector <bool> Keyword_presence;
    
    
    

    
    
    
};

#include "Classes.hpp"

#endif /* CLASSES_H */




