#ifndef CSVIO_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define CSVIO_H

# include <iostream>
using std::cout;
using std::endl;

#include <string>
using std::string;


// For use with the HAGRID experiment

void writeOut(ofstream &fileOUT, double peak1Low[],double peak1High[],double peak2Low[],double peak2High[]){

    for(int ii = 0;ii<13;ii++){
        fileOUT<<peak1Low[ii]<<",";
    }fileOUT<<endl;
    for(int ii = 0;ii<13;ii++){
        fileOUT<<peak1High[ii]<<",";
    }fileOUT<<endl;
    for(int ii = 0;ii<13;ii++){
        fileOUT<<peak2Low[ii]<<",";
    }fileOUT<<endl;
    for(int ii = 0;ii<13;ii++){
        fileOUT<<peak2High[ii]<<",";
    }fileOUT<<endl;
}


void readIn(ifstream &fileIN,double peaks[4][13]){

    string num = "";
    int row = 0;
    int column = 0;
    while (!fileIN.eof()) {
        char curr = fileIN.get();

        if (curr == ',') {
            peaks[row][column] = std::stod(num);
            // Reset
            num = "";
            column++;
            if (column > 12) {
                row++;
                column = 0;
            }
        }
        else {
            num += curr;
        }
    }
}


#endif
