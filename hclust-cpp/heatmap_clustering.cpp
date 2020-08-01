#include <math.h>
#include <string.h>

#include <string>
#include <vector>
#include <stdio.h>
#include <algorithm>

#include "fastcluster.h"


class RowVector {
public:
  std::string row_name;
  unsigned int row_id; // Starts at 1
  std::vector<int> row_elements;
  RowVector(std::string name, unsigned int cur_id, std::vector<int> elements) { row_name=name; row_id=cur_id; row_elements=elements; }
};


double euclidean_distance(RowVector vec1, RowVector vec2, int num_points) {
  double distance = 0;
  for (int i = 0; i < num_points; i++) {
    double diff = vec1.row_elements[i]-vec2.row_elements[i];
    distance += diff*diff;
  }
  distance = sqrt(distance);
  return distance;
}


// main program -- takes in CSV, including both column and row labels
int main(int argc, char** argv) 
{

  int i,j,k,npoints;

  // parse command line
  std::string opt_infile;
  int opt_method = HCLUST_METHOD_SINGLE;
  const char* usagemsg = "Usage: hclust-heatmap  <infile> [-m (single|complete|average|median)]\n";
  for (i=1; i<argc; i++) {
    if (0 == strcmp(argv[i], "-m")) {
      i++;
      if (i<argc) {
        if (0 == strcmp(argv[i], "single"))
            opt_method = HCLUST_METHOD_SINGLE;
        else if (0 == strcmp(argv[i], "complete"))
            opt_method = HCLUST_METHOD_COMPLETE;
        else if (0 == strcmp(argv[i], "average"))
            opt_method = HCLUST_METHOD_AVERAGE;
        else if (0 == strcmp(argv[i], "median"))
            opt_method = HCLUST_METHOD_MEDIAN;
        else {
          fputs(usagemsg, stderr);
          return 1;
        }
      } else {
        fputs(usagemsg, stderr);
        return 1;
      }
    }
    else if (argv[i][0] == '-') {
      fputs(usagemsg, stderr);
      return 1;
    }
    else {
      opt_infile = argv[i];
    }
  }
  if (opt_infile == "") {
    fputs(usagemsg, stderr);
    return 1;
  }
  
  // read row_vectors from input file
  unsigned int row_num, col_num = 0;
  std::vector<RowVector> row_vectors;

  std::ifstream ifile;
  ifile.open(opt_infile);
  if(!ifile) {
      std::cerr << "CSV file does not exist at given filepath" << std::endl;
      return 1;
  }

  std::ifstream  data(opt_infile);
  std::string line;
  while(std::getline(data,line))
  {
      std::stringstream  lineStream(line);
      std::string        cell;
      
      if (row_num == 0){
          row_num += 1;
          continue;
      }

      col_num = 0;
      RowVector cur_row_vector;
      std::vector<int> cur_values;
      while(std::getline(lineStream,cell,','))
      {
          if (col_num == 0){
              std::string cur_row_name = std::stof(cell);
              cur_row_vector = RowVector(cur_row_name, row_num, cur_values);
              continue;
          }

          cur_row_vector.row_elements.push_back(std::stof(cell));
          col_num += 1;
      }
      row_vectors.push_back(cur_row_vector);
      row_num += 1;
  }

  unsigned int num_vectors = row_num - 1; // Number of vectors. Exclude column labels in count
  unsigned int num_points = col_num - 1; // Number of points in each vector. Exclude row labels in count.

  // computation of condensed distance matrix
  double* distmat = new double[(num_vectors*(num_vectors-1))/2];
  k = 0;
  for (i=0; i<num_vectors; i++) {
    for (j=i+1; j<num_vectors; j++) {
      distmat[k] = euclidean_distance(row_vectors[i], row_vectors[j], num_points);
      k++;
    }
  }

  // clustering call
  int* merge = new int[2*(num_vectors-1)];
  double* height = new double[num_vectors-1];
  hclust_fast(num_vectors, distmat, opt_method, merge, height); // merge & height represents a stepwise dendrogram

  // TODO: Use 'merge' and 'height' to rearrange CSV. Also use it to visualize dendrogram. 
  


  // clean up
  delete[] distmat;
  delete[] merge;
  delete[] height;
  delete[] labels;

  return 0;
}