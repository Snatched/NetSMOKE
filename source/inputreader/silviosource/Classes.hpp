/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */



Reactor_unit::Reactor_unit()
{
    
     
     
     Isothermal_temperature = 0;
     Isothermal_temperature_unit_of_measurement = "--";
     

     
     UA = 0;
     UA_unit_of_measurement = "--";
    
}






int Mixer_unit::InputStreamCheck()

{
     for ( int j = 0; j < In_stream.size()-1; j++)
    {
        if ( In_stream.back() == In_stream.at(j) )
        {
            
                
                cout<<"Two Mixer input streams have the same number: "<< In_stream.back()<<".";
                return 1;
            
        }
    }
     
     return 0;
}

int Mixer_unit::InputOutputStreamCheck()

{
   for ( int i = 0; i < In_stream.size() ; i++ )
    {


        if ( In_stream.at(i) == Out_stream)
        {
            cout<<"ERROR at mixer "<<Name<<"!\n";
            cout<<"Output stream "<<Out_stream<<" has the same number of one of the input streams.";
            return 1;
        }
    }
   
   return 0;
}








int S_unit::OutputStreamCheck()
{
     for ( int j = 0; j < Out_stream.size()-1; j++)
    {
        if ( Out_stream.back() == Out_stream.at(j) )
        {
            
                
                cout<<"Two Splitter output streams have the same number: "<< Out_stream.back()<<".";
                return 1;
            
        }
    }
     
     return 0;
}

int S_unit::OutputInputStreamCheck()
{
    for ( int i = 0; i < Out_stream.size() ; i++ )
    {


        if ( Out_stream.at(i) == In_stream)
        {
            cout<<"ERROR at splitter "<<Name<<"!\n";
            cout<<"Input stream "<<In_stream<<" has the same number of one of the output streams.";
            return 1;
        }
    }
   
   return 0;
}










int PS_unit::OutputStreamCheck()
{
    for ( int j = 0; j < Out_stream.size()-1; j++)
    {
        if ( Out_stream.back() == Out_stream.at(j) )
        {
            
                
                cout<<"Two Phase-splitter output streams have the same number: "<< Out_stream.back()<<".";
                return 1;
            
        }

        if ( Splitted_phase.back() == Splitted_phase.at(j) )
        {
            
               
                cout<<"Two output streams "<<Out_stream.back()<<" and "<<Out_stream.at(j)<<" have the same phase: "<< Splitted_phase.back()<<".";
                return 1;
            
        }
    }
    
    return 0;
}

int PS_unit::OutputInputStreamCheck()
{
    for ( int i = 0; i < Out_stream.size() ; i++ )
    {


        if ( Out_stream.at(i) == In_stream)
        {
            cout<<"ERROR at Phase-splitter "<<Name<<"!\n";
            cout<<"Input stream "<<In_stream<<" has the same number of one of the output streams.";
            return 1;
        }
    }
   
   return 0;
}










Stream_class::Stream_class()
{
    
    
    Phase = "Mix";//default stream phase

    MassFlowRate = 0;
    MassFlowRate_unit_of_measurement = "--";
    
    Temperature = 0;
    Temperature_unit_of_measurement = "--";
    
    
      
    

    
    
}



void Stream_class::ShrinkingKeyword_presence()
{
    Keyword_presence.clear();
    Keyword_presence.shrink_to_fit();

}

void Stream_class::GivingSizeToKeyword_presence(int size_of_Stream_keyword)
{
    for( int i = 0; i < size_of_Stream_keyword; i++)
    {
        Keyword_presence.push_back(false);

    }
    

    
}

int Stream_class::CompoundControl()
{
    for ( int i = 0; i < Compound.size()-1 ; i++)
    {
        if ( Compound.back() == Compound.at(i))
        {
            cout<<"ERROR at stream "<<Name<<"!Compound "<<Compound.at(i)<<" is repeated twice!";
            return 1;
        }
    }
    
    return 0;
}



double Stream_class::MassFractionSum()
{
    
    double Sum_of_mass_fractions = 0;
    
    for ( int i = 0; i < Composition.size() ; i++)
        
    {
        
        Sum_of_mass_fractions += Composition.at(i);
        
    }
    
    return Sum_of_mass_fractions;
}



int Stream_class::TrueKeyword(int index)
{
    if (Keyword_presence.at(index) == true)
    {
        return 1;
    }
    
    else
    {
        Keyword_presence.at(index) = true;
        
        return 0;
    }
    
}



