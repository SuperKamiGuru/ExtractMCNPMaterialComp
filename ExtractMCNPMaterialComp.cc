#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include "include/ElementNames.hh"
#include <iomanip>
#include "include/IsotopeMass.hh"

using namespace std;

// this program uses the isotope mass data from ENDFVII contained in the IsotopeMass class file by default
// make sure that the material consistently use wt% or abundance to define their composition otherwise errors may occur in the output

// X create a static container class that contains the isotope mass used by the MCNP input file using another program
// add ability for this code to determine the ENDF version that is being used by the isotope and have it run the previously mentioned program and create a temporary
// static container with the isotope mass list for that version of ENDF and use it
void GetDataStream( string geoFileName, std::stringstream& ss);
void FormatData(std::stringstream& stream, bool wtPer);
void Swap(std::vector<int> &matNumVec, std::vector<string> &matDensVec,int pos1, int pos2);
void GetIsotopeData(stringstream &stream, stringstream &streamOut, char &letter, int matNum, int count, std::vector<string> &matDensVec, int matDegen, bool wtPer);
string CreateMacroName(string geoFileName, string outDirName);
void SetDataStream( string macroFileName, std::stringstream& ss);


int main(int argc, char** argv)
{

    string outDirName, mcnpSource, wtPercent;
    bool wtPer=true;

    ElementNames elementNames;
    elementNames.SetElementNames();

    IsotopeMass isotopeMass;
    isotopeMass.SetIsotopeMass();

    std::stringstream streamS;
    string macroFileName;

    if(argc==4)
    {
        outDirName = argv[1];
        mcnpSource = argv[2];
        wtPercent = argv[3];

        if((wtPercent=="false")||(wtPercent=="False")||(wtPercent=="abundance")||(wtPercent=="Abundance")||(wtPercent=="atom%")||(wtPercent=="Atom%")||(wtPercent=="atomic%")||(wtPercent=="Atomic%"))
        {
            wtPer=false;
        }
    }
    else if(argc==3)
    {
        outDirName = argv[1];
        mcnpSource = argv[2];
    }
    else
    {
        cout << "\nGive the the output directory, then specify the MCNP input file that you want to convert that you want to convert\n" <<  endl;
        elementNames.ClearStore();
        isotopeMass.ClearStore();
        return 1;
    }

    GetDataStream(mcnpSource, streamS);

    FormatData(streamS, wtPer);

    macroFileName = CreateMacroName(mcnpSource, outDirName);

    SetDataStream(macroFileName, streamS);

    elementNames.ClearStore();
    isotopeMass.ClearStore();
    return 0;
}

void GetDataStream( string geoFileName, std::stringstream& ss)
{
    string* data=NULL;

    // Use regular text file
    std::ifstream thefData( geoFileName.c_str() , std::ios::in | std::ios::ate );
    if ( thefData.good() )
    {
        int file_size = thefData.tellg();
        thefData.seekg( 0 , std::ios::beg );
        char* filedata = new char[ file_size ];
        while ( thefData )
        {
            thefData.read( filedata , file_size );
        }
        thefData.close();
        data = new string ( filedata , file_size );
        delete [] filedata;
    }
    else
    {
    // found no data file
    //                 set error bit to the stream
        ss.setstate( std::ios::badbit );
    }
    if (data != NULL)
    {
        ss.str(*data);
        if(data->back()!='\n')
            ss << "\n";
        ss.seekg( 0 , std::ios::beg );
    }

    delete data;
}

void FormatData(std::stringstream& stream, bool wtPer)
{
    char line[256];
    char letter;
    std::vector<int> matNumVec;
    std::vector<string> matDensVec;
    std::vector<int> matDegenVec;
    stringstream numConv, streamOut;
    bool found=false;
    string densUnits[2] = {"10^24 atoms/cm3", "g/cm3"};
    int num, type, degen=1, count, matNum;

    while(stream)
    {
        letter = stream.get();
        if((letter!='c')&&(letter!='$'))
        {
            stream.getline(line, 256);
            break;
        }
        else
        {
            stream.getline(line, 256);
        }
    }

    while(stream)
    {
        letter = stream.get();
        if((letter>='0')&&(letter<='9'))
        {
            while((letter>='0')&&(letter<='9'))
            {
                letter = stream.get();
            }
            while(!(((letter>='0')&&(letter<='9'))||((letter>='a')&&(letter<='z'))||((letter>='A')&&(letter<='Z'))))
            {
                letter = stream.get();
            }
            if((letter>='0')&&(letter<='9'))
            {
                while((letter>='0')&&(letter<='9'))
                {
                    numConv << letter;
                    letter = stream.get();
                }
                numConv >> num;
                matNumVec.push_back(num);
                numConv.clear();
                numConv.str("");
            }
            else
            {
                stream.getline(line, 256);
                break;
            }
            while(!((letter=='+')||(letter=='-')||((letter>='0')&&(letter<='9'))||(letter=='.')))
            {
                letter = stream.get();
            }
            if(((letter>='0')&&(letter<='9'))||(letter=='.'))
            {
                type=0;
            }
            else if(letter == '+')
            {
                type=0;
                letter = stream.get();
            }
            else
            {
                type=1;
                letter = stream.get();
            }

            while(((letter>='0')&&(letter<='9'))||(letter=='.'))
            {
                numConv << letter;
                letter = stream.get();
            }
            matDensVec.push_back(numConv.str()+densUnits[type]);
            numConv.clear();
            numConv.str("");
            stream.getline(line, 256);
        }
        else
        {
            stream.getline(line, 256);
        }
    }

    if(matNumVec.size()>0)
    {
        for(int i=0; i<int(matNumVec.size()); i++)
        {
            for(int j=i+1; j<int(matNumVec.size()); j++)
            {
                if(matNumVec[j]<matNumVec[i])
                {
                    Swap(matNumVec, matDensVec, i, j);
                }
                else if(matNumVec[j]==matNumVec[i])
                {
                    if(matDensVec[j]==matDensVec[i])
                    {
                        matNumVec.erase(matNumVec.begin()+j);
                        matDensVec.erase(matDensVec.begin()+j);
                        j--;
                    }
                    else if(matDensVec[j]<matDensVec[i])
                    {
                        Swap(matNumVec, matDensVec, i, j);
                    }
                }
            }
        }

        for(int i=1; i<int(matNumVec.size()); i++)
        {
            if(matNumVec[i]==matNumVec[i-1])
            {
                degen++;
            }
            else
            {
                matDegenVec.insert(matDegenVec.end(), degen, degen);
                degen=1;
            }
        }
        matDegenVec.insert(matDegenVec.end(), degen, degen);
    }
    else
    {
        cout << "\nError: Could not find cell card" << endl;
        return;
    }

    while(stream)
    {
        letter = stream.get();
        if(letter=='m')
        {
            letter = stream.get();
            if((letter>='0')&&(letter<='9'))
            {
                while((letter>='0')&&(letter<='9'))
                {
                    numConv << letter;
                    letter = stream.get();
                }
                numConv >> matNum;
                numConv.clear();
                numConv.str("");
                count=0;

                while((matNumVec[count]!=matNum)&&(count!=int(matNumVec.size())))
                {
                    count++;
                }
                if(count!=int(matNumVec.size()))
                {
                    GetIsotopeData(stream, streamOut, letter, matNum, count, matDensVec, matDegenVec[count], wtPer);
                    found=true;
                }
            }
            else
            {
                stream.getline(line, 256);
            }
        }
        else
        {
            stream.getline(line, 256);
        }
    }

    if(found)
    {
        stream.clear();
        stream.str(streamOut.str());
    }
    else
    {
        cout << "\nError: Could not find material card" << endl;
    }
}

void Swap(std::vector<int> &matNumVec, std::vector<string> &matDensVec,int pos1, int pos2)
{
    int swapNum;
    string swapDens;

    swapNum = matNumVec[pos1];
    swapDens = matDensVec[pos1];

    matNumVec[pos1] = matNumVec[pos2];
    matDensVec[pos1] = matDensVec[pos2];

    matNumVec[pos2] = swapNum;
    matDensVec[pos2] = swapDens;
}

void GetIsotopeData(stringstream &stream, stringstream &streamOut, char &letter, int matNum, int count, std::vector<string> &matDensVec, int matDegen, bool wtPer)
{
    ElementNames* elementNames;
    IsotopeMass* isotopeMass;
    std::vector<string> isoNameVec;
    std::vector<double> isoTempVec;
    std::vector<double> isoAmountVec;
    std::vector<double> isoAbunVec;
    std::vector<int> isoAmountTypeVec;
    std::vector<double> isoMassVec;
    double libTemp[5] = {293.6, 600, 900, 1200, 2500};
    string compUnits[2] = {"atomic%", "weight%"};
    string name;
    stringstream numConv;
    char line[256];
    int A, Z, num, index;
    double amount;

    while(true)
    {
        while(!((letter>='0')&&(letter<='9')))
        {
            letter=stream.get();
        }
        while((letter>='0')&&(letter<='9'))
        {
            numConv << letter;
            letter = stream.get();
        }
        numConv >> num;
        numConv.clear();
        numConv.str("");

        Z = floor(num/1000);
        A = num-Z*1000;

        numConv << Z << '_' << A << '_' << elementNames->GetName(Z);
        isoNameVec.push_back(numConv.str());
        isoMassVec.push_back(isotopeMass->GetIsotopeMass(Z, A));
        numConv.clear();
        numConv.str("");

        letter=stream.get();
        letter=stream.get();
        numConv << letter;
        numConv >> index;
        isoTempVec.push_back(libTemp[index]);
        numConv.clear();
        numConv.str("");

        letter=stream.get();
        letter=stream.get();
        while(!((letter=='+')||(letter=='-')||((letter>='0')&&(letter<='9'))||(letter=='.')))
        {
            letter = stream.get();
        }
        if(((letter>='0')&&(letter<='9'))||(letter=='.'))
        {
            isoAmountTypeVec.push_back(0);
        }
        else if(letter == '+')
        {
            isoAmountTypeVec.push_back(0);
            letter = stream.get();
        }
        else
        {
            isoAmountTypeVec.push_back(1);
            letter = stream.get();
        }

        while(((letter>='0')&&(letter<='9'))||(letter=='.'))
        {
            numConv << letter;
            letter = stream.get();
        }

        numConv >> amount;
        isoAmountVec.push_back(amount);
        numConv.clear();
        numConv.str("");

        while((letter!='&')&&(letter!='\n')&&(stream))
        {
            letter=stream.get();
        }
        stream.getline(line, 256);
        if(letter!='&')
        {
            break;
        }
    }

    double sum=0.;
    int state=0;

    for(int i=0; i<int(isoAmountVec.size()); i++)
    {
        if(isoAmountTypeVec[i]==1)
        {
            state=1;
            sum += isoAmountVec[i]/isoMassVec[i];
        }
        else
        {
            state=-1;
            sum += isoAmountVec[i]*isoMassVec[i];
            isoAbunVec.push_back(isoAmountVec[i]);
        }
    }

    if(state==1)
    {
        for(int i=0; i<int(isoAmountVec.size()); i++)
        {
            isoAbunVec.push_back(isoAmountVec[i]/(isoMassVec[i]*sum));
        }
    }

    if(wtPer)
    {
        if(state==-1)
        {
            for(int i=0; i<int(isoAmountVec.size()); i++)
            {
                isoAmountVec[i]=(isoAmountVec[i]*isoMassVec[i]/sum);
                isoAmountTypeVec[i]=1;
            }
        }
    }

    else
    {
        for(int i=0; i<int(isoAmountVec.size()); i++)
        {
            isoAmountVec[i]=isoAbunVec[i];
            isoAmountTypeVec[i]=0;
        }
    }

    for(int i=0; i<int(isoNameVec.size()); i++)
    {
        for(int j=i+1; j<int(isoNameVec.size()); j++)
        {
            if(isoNameVec[i]==isoNameVec[j])
            {
                isoTempVec[i] = (isoTempVec[i]*isoAbunVec[i]+isoTempVec[j]*isoAbunVec[j])/(isoAbunVec[i]+isoAbunVec[j]);
                isoAmountVec[i] += isoAmountVec[j];
                isoAbunVec[i] += isoAbunVec[j];
                isoNameVec.erase(isoNameVec.begin()+j);
                isoAmountVec.erase(isoAmountVec.begin()+j);
                isoAbunVec.erase(isoAbunVec.begin()+j);
                isoAmountTypeVec.erase(isoAmountTypeVec.begin()+j);
                isoMassVec.erase(isoMassVec.begin()+j);
                isoTempVec.erase(isoTempVec.begin()+j);
                j--;
            }
        }
    }
    sum=0.;
    for(int i=0; i<int(isoAmountVec.size()); i++)
    {
        sum += isoAmountVec[i];
    }
    for(int i=0; i<int(isoAmountVec.size()); i++)
    {
        isoAmountVec[i]/=sum;
    }

    numConv.clear();
    numConv.str("");
    numConv << "Material #: " << matNum;

    streamOut.fill('-');
    streamOut << std::setw(84) << std::left << numConv.str() << '\n' << endl;
    streamOut << "Densities Used: ";
    for(int i=0; i<matDegen; i++)
    {
        streamOut << matDensVec[count+i] << " ";
    }
    streamOut << '\n' << endl;
    streamOut.fill(' ');
    streamOut << std::setw(20) << std::left << "Isotope Name:" << std::setw(20) << std::left << "Amount:" << std::setw(20) << std::left << "Mass:"
              << std::setw(20) << std::left << "Temperature:" << '\n' << endl;

    for(int i=0; i<int(isoNameVec.size()); i++)
    {
        numConv.clear();
        numConv.str("");
        numConv << isoAmountVec[i] << compUnits[isoAmountTypeVec[i]];

        streamOut << std::setw(20) << std::left << isoNameVec[i] << std::setw(20) << std::left << numConv.str();

        numConv.clear();
        numConv.str("");
        numConv << isoMassVec[i] << "amu";

        streamOut << std::setw(20) << std::left << numConv.str();

        numConv.clear();
        numConv.str("");
        numConv << isoTempVec[i] << 'k';

        streamOut << std::setw(20) << std::left << numConv.str() << endl;
    }

    streamOut << '\n';

}

string CreateMacroName(string geoFileName, string outDirName)
{
    size_t pos = geoFileName.find_last_of('/');
    size_t pos2 = std::string::npos;
    if(pos == std::string::npos)
        pos=0;
    else
        pos++;

    return (outDirName+"MatComp"+geoFileName.substr(pos, pos2-pos));
}

void SetDataStream( string macroFileName, std::stringstream& ss)
{
  std::ofstream out( macroFileName.c_str() , std::ios::out | std::ios::trunc );
  if ( ss.good() )
  {
     ss.seekg( 0 , std::ios::end );
     int file_size = ss.tellg();
     ss.seekg( 0 , std::ios::beg );
     char* filedata = new char[ file_size ];

     while ( ss )
     {
        ss.read( filedata , file_size );
        if(!file_size)
        {
            cout << "\n #### Error the size of the stringstream is invalid ###" << endl;
            break;
        }
     }

     out.write(filedata, file_size);
     if (out.fail())
    {
        cout << endl << "writing the ascii data to the output file " << macroFileName << " failed" << endl
             << " may not have permission to delete an older version of the file" << endl;
    }
     out.close();
     delete [] filedata;
  }
  else
  {
// found no data file
//                 set error bit to the stream
     ss.setstate( std::ios::badbit );

     cout << endl << "### failed to write to ascii file " << macroFileName << " ###" << endl;
  }
   ss.str("");
}
