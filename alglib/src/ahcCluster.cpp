#include "stdafx.h"
#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <vector>
#include <fstream>
#include "dataanalysis.h"

using namespace alglib;
using namespace std;

class RowVector {
    public:
        unsigned int row_id;
        vector<int> row_elements;
        RowVector(unsigned int cur_id, vector<int> elements) {
            row_id = cur_id;
            row_elements = elements;
        }
};

double euclidean_distance(RowVector vec1, RowVector vec2, int num_points) {
    double distance = 0;
    for (int i=0; i<num_points; i++) {
        double diff = vec1.row_elements[i] - vec2.row_elements[i];
        distance += diff;
    }
    return sqrt(distance);
}

// command-line: [exec] [input-file]
int main(int argc, char **argv) {

    

    ifstream ifile;
    ifile.open(argv[1]);
    if (!ifile) {
        cerr << "INPUT ERROR: Enter valid CSV heatmap file" << endl;
        return 1;
    }

    unsigned int row_num, col_num = 0;
    vector<RowVector> row_vectors;

    string infile = argv[1];
    ifstream data(infile);
    string line;

    // extract values from rows in heatmap csv data
    while (getline(data,line)) {
        stringstream lineStream(line);
        string cell;

        (row_num == 0) ? row_num++ : continue;

        col_num = 0;
        RowVector current_row;
        vector<int> row_values;

        while (getline(lineStream,cell,',')) {
            if (col_num = 0) {
                current_row = RowVector(row_num,row_values);
                continue;
            }

            current_row.row_elements.push_back(stof(cell));
            col_num++;
        }

        row_vectors.push_back(current_row);
        row_num++;
    }

    unsigned int num_vectors = row_num - 1;
    unsigned int num_points = col_num - 1;

    // compute condensed distance matrix 
    double* matrix = new double[(num_vectors * (num_vectors-1))/2];
    k = 0;
    for (int i=0; i<num_vectors; i++) {
        for (int j=i+1; j<num_vectors; j++) {
            matrix[k] = euclidean_distance(row_vectors[i], row_vectors[j], num_points);
            k++;
        }
    }

    
    /*
        ALGLIB HIERARCHICAL CLUSTERING 
    */

    clusterizerstate s;
    ahcreport report;
    real_2d_array xy = "[[1,0],[1,0],[4,0],[2,0],[4,0]]";

    clusterizercreate(s);
    clusterizersetahcalgo(s,1); // single linkage
    clusterizersetpoints(s, xy, 2);
    clusterizerrunahc(s, report);

    printf("%s\n", report.z.tostring().c_str());
    printf("%s\n", report.p.tostring().c_str());
    printf("%s\n", report.pm.tostring().c_str());
    return 0;

}