/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Function.h
 * Author: silvio14
 *
 * Created on June 25, 2017, 9:20 PM
 */

#ifndef FUNCTIONS_H
#define FUNCTIONS_H



#include <vector>

#include <string>

#include <iostream>

using namespace std;


vector <string> String_comma_separator(const string string1)
{
    
    string::const_iterator iterator1 = string1.begin(); //iterator that points to the first letter of string1
    
    string string2; //auxiliary string
    
    vector <string> Splitted_string; //this vector will be filled with the words of string1
    
    while ( iterator1 != string1.end() ) //iterating till the last letter of string1
        
    {
         
            ( string2 ).push_back( *iterator1 ); //Adding the  letter to string2

            if ( isspace( *iterator1 ) != 0 || *iterator1 == ',' )  //if this "letter" is a space or a comma, then remove it

            {

                ( string2 ).erase ( ( string2 ).end()-1 );   

                if (( string2 ).empty() == false) //if string2 is not empty(that means it is filled with some letters which are not separated by a space or a comma)

                {


                Splitted_string.push_back( string2 ); //adding the first word to the vector 

                ( string2 ).clear(); // Cleaning the auxiliary string


                }

            }
            
            if ( *iterator1 == ';')  
            
            {
                
                ( string2 ).erase ( ( string2 ).end()-1 );  //if a ';' is found, then it is removed because the previous letters build up a word

                if (( string2 ).empty() == false) 

                {

                Splitted_string.push_back( string2 );
                ( string2 ).clear();


                }  

            }
            
            else if ( iterator1 == string1.end()-1 ) // if iterator has arrived to the last letter of string1 and if the auxiliary string2 is not empty, then the string2 builds up another word

            {

                if (( string2 ).empty() == false) 

                {

                Splitted_string.push_back( string2 );
                 ( string2 ).clear();


                }
            }
            
            

            iterator1++;
            
        
    }
    
    return Splitted_string; 
    
}


int LengthControl(const vector <string> string1,const int length, const int number_of_line) // This function checks if vector <string> size is equal to length

{
    if ( string1.size() != length )
        
    {
        cout<<"ERROR! Information at line "<<number_of_line<<" is not correct!\n";     
        return 1;
        
    }
    
return 0;
    
}

string IntControl (const string string1)  // This function checks if the string is made of numbers only
    
    {
       string::const_iterator iterator1 = string1.begin();
       
      
       
       while ( iterator1 != string1.end() )
           
       {
           if ( isdigit(*iterator1) == false )
               
           {
               cout<<"Conversion of "<<string1<<" to integer failed!\n";
               
               return "No";
           }
           
        iterator1++;
        
       }
       
       
     
       return string1;
    }

string DoubleControl (const string string1)    //This function checks if the string is a double
    
    {
       string::const_iterator iterator1 = string1.begin();
       
      
       
       
       
       int point_count;
       
       point_count = 0;
       
       while ( iterator1 != string1.end() )
           
       {
           if ( isdigit(*iterator1) == false )
               
           {
                if ( iterator1 == string1.begin() && *iterator1 == '-' )
                   
                {
                    iterator1++;
                    break;
                }
               
                if ( *iterator1 == '.')

                {
                    if ( point_count == 0)

                    {
                        point_count = 1;
                    }

                    else

                    {
                        cout<<"Conversion of \""<<string1<<"\" to double failed!\n";
                        return "No";
                    }

                }

                else

                {
                cout<<"Conversion of \""<<string1<<"\" to double failed!\n";
                return "No";
                }
            
               
           }
           
        
        
        iterator1++;
        
       }
       
       
       
     
       return string1;
    }


 string LetterControl(const string string1) //This function checks if the string has letters only
  
    {
       string::const_iterator iterator1 = string1.begin();
      
       int letter_count = 0;
            
        while ( iterator1 != string1.end() )
           
       {
            
            if ( isalpha(*iterator1) != 0)
            {
                
                letter_count++;
            }  
             
           
        iterator1++;
        
       }
            
       if ( letter_count == 0 )     
       {
           
           return "No";
       }
             
            
       else
       {
         return string1;  
       }
            
                
    }
 
  string FirstLetterThenNumberControl(const string string1,char first_letter) //This function checks if the string is made of an initial letter followed by numbers
  
  {
       string::const_iterator iterator1 = string1.begin();
       
       if (*iterator1 != first_letter)
       {
           cout<<"Name \""<<string1<<"\" is wrong! It must start with letter: "<<first_letter<<"\n";
               
           return "No";
       }
       
       iterator1++;
       
       if ( first_letter == 'P' )  // if the first letter is a 'P' then the second must be an 'S' ( Phase-splitter declaration)
       {
           
           
           if (*iterator1 != 'S')
           {
               cout<<"ERROR!Unit \""<<string1<<"\" not recognized! \n";
                      
               return "No";        
           }
           
           else
           {
               iterator1++;
           }
       }
        while ( iterator1 != string1.end() )
           
       {
           if ( isdigit(*iterator1) == false )
               
           {
               cout<<"\""<<*iterator1<<"\" is not a number!\n";
               cout<<"\n"<<string1<<"\" must be made of one letter(or two, for Phase-splitters) followed by numbers!\n";
               
               return "No";
           }
           
        iterator1++;
        
       }
       
       return string1;
       
  }
  
  
  int End_of_file(ifstream& file) // this function checks if file is finished
  {
      if (file.eof())
      {
          cout<<"ERROR! \"//\" is missing at the end of file!";
          return 1;
      }
      
      return 0;
  }
  
  int New_declaration(string line,int line_number)
  {
      
    vector <string> Main_keyword;
    
    Main_keyword.push_back("Reactor");
    Main_keyword.push_back("Mixer");
    Main_keyword.push_back("Splitter");
    Main_keyword.push_back("Phase-Splitter");
    Main_keyword.push_back("Stream");
       
    
    for ( int i = 0; i <Main_keyword.size() ; i++)
    {
        if (line.find(Main_keyword.at(i)) != string::npos)
        {
            cout<<"ERROR!\n";
            cout<<"There is a new unit or stream declaration at line "<<line_number<<"!\n";
            cout<<"Remember to end unit or stream declarations with \'//\'.";
            return 1;
        }
    }
    
    return 0;
  }
  
  int DictionaryControl(string Statement, vector <string> Dictionary) // this function checks if Statement is one of the elements of Dictionary
{
    for ( int i = 0 ; i < Dictionary.size() ; i++)
    {
        if (Statement == Dictionary.at(i) )
        {
            
            return 0;
        }

    }

    cout <<"ERROR! \""<<Statement<<"\" not recognized!\n";
    cout<<"Maybe you have to update the related keyword dictionary (";

    for ( int i = 0 ; i < Dictionary.size() ; i++)
    {
        cout<<Dictionary.at(i)<<" ";
    }
    cout<<")\n";
    return 1;
}







#endif /* FUNCTIONS_H */



