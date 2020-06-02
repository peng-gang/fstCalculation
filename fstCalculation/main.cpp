//
//  main.cpp
//  fstCalculation
//
//  Created by Gang Peng on 6/1/20.
//  Copyright © 2020 Gang Peng. All rights reserved.
//

#include <iostream>
#include <fstream>

#include "normal.h"

using namespace std;

int main(int argc, char ** argv) {
/*
 * -i input.txt // input file containing MAF informaiton
 * -o output.txt // output file
 * -r race1 race2 race3 ... // race groups to compare. these groups should be in the first row (header) of the input file. -r all: comparing all pairs  -r other: comparing every race to other
 */
    
    vector<string> mustOptions = {"-i", "-o", "-r"};
    vector<string> allOptions = {"-i", "-o", "-r"};
    
    map<string, vector<string> > cmLine = parseCMLine(argc, argv, allOptions, mustOptions);
    if(cmLine.size()==0){
        return -1;
    }
    
    string fnameIn = cmLine["-i"][0];
    string fnameOut = cmLine["-o"][0];
    
    ifstream fin(fnameIn.c_str());
    if(!fin.is_open()){
        cout << "Cannot open file " << fnameIn << ". Please check the input file name." << endl;
        return -1;
    }
    
    ofstream fout(fnameOut.c_str());
    if(!fout.is_open()){
        cout << "Cannot open file " << fnameOut << ". Please check the output file name." << endl;
        return -1;
    }
    
    string header;
    getline(fin, header);
    
    //check header and race in the header
    vector<string> vsHeader = split(header, "\t");
    
    vector<size_t> idxRace;
    if(cmLine["-r"].size()==1){
        if(cmLine["-r"][0] == "other"){
            idxRace.clear();
        } else if(cmLine["-r"][0] == "all"){
            idxRace.clear();
            // There is a bug here
            // I assume race info from 7 to 11 // 6-10
            for(size_t i=6; i<11; i++){
                idxRace.push_back(i);
            }
        } else {
            vector<string>::iterator it = find(vsHeader.begin(), vsHeader.end(), cmLine["-r"][0]);
            if(it!=vsHeader.end()){
                idxRace.push_back(distance(vsHeader.begin(), it));
            } else {
                cout << "Cannot find race \"" << cmLine["-r"][0] << "\" in the header of the input file. Plese take a look at the header of the input file." << endl;
                return -1;
            }
        }
    } else {
        for(size_t i=0; i<cmLine["-r"].size(); i++){
            vector<string>::iterator it = find(vsHeader.begin(), vsHeader.end(), cmLine["-r"][i]);
            if(it!=vsHeader.end()){
                idxRace.push_back(distance(vsHeader.begin(), it));
            } else {
                cout << "Cannot find race \"" << cmLine["-r"][i] << "\" in the header of the input file. Plese take a look at the header of the input file." << endl;
                return -1;
            }
        }
    }
    
    
    size_t idxSymbol=0;
    vector<string>::iterator itSymbol = find(vsHeader.begin(), vsHeader.end(), "SYMBOL");
    if(itSymbol == vsHeader.end()){
        cout << "Cannot find SYMBOL column in input file. Please make sure \"" << fnameIn <<"\" is the right input file." <<endl;
        return -1;
    } else {
        idxSymbol = distance(vsHeader.begin(), itSymbol);
    }
    
    cout << "Start" <<endl;
    long long numVar = 0;
    
    if(idxRace.size() == 0){
        // There is a bug here
        // I assume race info from 7 to 11 // 6-10
        fout<<"Chr\tPos\tSymbol";
        for(size_t i=6; i<11; i++){
            fout << "\t" << vsHeader[i] << "-Other";
        }
        fout << endl;
        
        
        while(!fin.eof()){
            string fline;
            getline(fin, fline);
            if(fline.size() < 2){
                break;
            }
            
            numVar++;
            if(numVar % 100000 == 0){
                cout << numVar << endl;
            }
            
            vector<string> vsLine = split(fline, "\t");
            
            fout << vsLine[0] << "\t" << vsLine[1] << "\t" << vsLine[idxSymbol];
            for(size_t i=6; i<11; i++){
                if(vsLine[i]=="NA"){
                    fout << "\tNA";
                    continue;
                } else {
                    double maf1 = stof(vsLine[i]);
                    double maf2 = 0;
                    int numOther = 0;
                    for(size_t j=6; j<11; j++){
                        if(j!=i){
                            if(vsLine[j] != "NA"){
                                maf2 += stof(vsLine[j]);
                                numOther++;
                            }
                        }
                    }
                    if(numOther==0){
                        fout << "\tNA";
                        continue;
                    }
                    
                    maf2 = maf2/numOther;
                    //calculate fst
                    double mafAll = (maf1 + maf2)/2.0;
                    double hAll = 1.0 - mafAll*mafAll - (1.0-mafAll)*(1.0-mafAll);
                    double h1 = 1.0 - maf1*maf1 - (1.0-maf1)*(1.0-maf1);
                    double h2 = 1.0 - maf2*maf2 - (1.0-maf2)*(1.0-maf2);
                    if((hAll - (h1+h2)/2.0)==0){
                        fout <<"\t" << 0.0;
                    } else {
                        double fst = (hAll - (h1+h2)/2.0)/hAll;
                        fout << "\t" << fst;
                    }
                }
            }
            fout << endl;
        }
    } else if(idxRace.size() == 1) {
        fout<<"Chr\tPos\tSymbol";
        fout<<"\t"<<cmLine["-r"][0]<<"-Other"<<endl;
        
        while(!fin.eof()){
            string fline;
            getline(fin, fline);
            if(fline.size() < 2){
                break;
            }
            
            numVar++;
            if(numVar % 100000 == 0){
                cout << numVar << endl;
            }
                   
            vector<string> vsLine = split(fline, "\t");
            fout << vsLine[0] << "\t" << vsLine[1] << "\t" << vsLine[idxSymbol];
            
            // There is a bug here
            // I assume race info from 7 to 11 // 6-10
            if(vsLine[idxRace[0]] == "NA"){
                fout << "\tNA";
                continue;
            }
            
            double maf1 = stof(vsLine[idxRace[0]]);
            
            double maf2 = 0;
            int numOther = 0;
            for(size_t i=6; i<11; i++){
                if(i != idxRace[0]){
                    if(vsLine[i] != "NA"){
                        maf2 += stof(vsLine[i]);
                        numOther++;
                    }
                }
            }
            
            if(numOther==0){
                fout << "\tNA";
                continue;
            }
            double mafAll = (maf1 + maf2)/2.0;
            double hAll = 1.0 - mafAll*mafAll - (1.0-mafAll)*(1.0-mafAll);
            double h1 = 1.0 - maf1*maf1 - (1.0-maf1)*(1.0-maf1);
            double h2 = 1.0 - maf2*maf2 - (1.0-maf2)*(1.0-maf2);
            if((hAll - (h1+h2)/2.0)==0){
                fout <<"\t" << 0.0;
            } else {
                double fst = (hAll - (h1+h2)/2.0)/hAll;
                fout << "\t" << fst;
            }
            fout << endl;
        }
    } else {
        fout<<"Chr\tPos\tSymbol";
        for(size_t i=0; i<(idxRace.size()-1); i++){
            for(size_t j=(i+1); j<idxRace.size(); j++){
                fout<<"\t"<<vsHeader[idxRace[i]]<<"-"<<vsHeader[idxRace[j]];
            }
        }
        fout<<endl;
        
        while(!fin.eof()){
            string fline;
            getline(fin, fline);
            
            if(fline.size() < 2){
                break;
            }
            
            numVar++;
            if(numVar % 100000 == 0){
                cout << numVar << endl;
            }
            
            vector<string> vsLine = split(fline, "\t");
            fout << vsLine[0] << "\t" << vsLine[1] << "\t" << vsLine[idxSymbol];
            
            for(size_t i=0; i<(idxRace.size()-1); i++){
                for(size_t j=(i+1); j<idxRace.size(); j++){
                    if(vsLine[idxRace[i]] != "NA" && vsLine[idxRace[j]] != "NA"){
                        double maf1 = stof(vsLine[idxRace[i]]);
                        double maf2 = stof(vsLine[idxRace[j]]);
                        double mafAll = (maf1 + maf2)/2.0;
                        double hAll = 1.0 - mafAll*mafAll - (1.0-mafAll)*(1.0-mafAll);
                        double h1 = 1.0 - maf1*maf1 - (1.0-maf1)*(1.0-maf1);
                        double h2 = 1.0 - maf2*maf2 - (1.0-maf2)*(1.0-maf2);
                        if((hAll - (h1+h2)/2.0)==0){
                            fout <<"\t" << 0.0;
                        } else {
                            double fst = (hAll - (h1+h2)/2.0)/hAll;
                            fout << "\t" << fst;
                        }
                       
                    } else {
                        fout <<"\tNA";
                    }
                }
            }
            fout<<endl;
        }
    }
    
    fin.close();
    fout.close();
    
    cout << "Total Variants: " << numVar << endl;
    
    cout << "Finished!" << endl;
    
    return 0;
}
