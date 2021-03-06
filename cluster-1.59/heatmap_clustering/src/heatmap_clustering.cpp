/* heatmap - High performance heatmap creation in C.
 *
 * The MIT License (MIT)
 *
 * Copyright (c) 2013 Lucas Beyer
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the "Software"), to deal in
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
 * the Software, and to permit persons to whom the Software is furnished to do so,
 * subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
 * FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
 * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
 * IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>

#include "lodepng.h"
#include "heatmap.h"
extern "C" {
    #include "cluster.h"
}


#include "colorschemes/gray.h"
#include "colorschemes/Blues.h"
#include "colorschemes/BrBG.h"
#include "colorschemes/BuGn.h"
#include "colorschemes/BuPu.h"
#include "colorschemes/GnBu.h"
#include "colorschemes/Greens.h"
#include "colorschemes/Greys.h"
#include "colorschemes/Oranges.h"
#include "colorschemes/OrRd.h"
#include "colorschemes/PiYG.h"

#include "colorschemes/PRGn.h"
#include "colorschemes/PuBuGn.h"
#include "colorschemes/PuBu.h"
#include "colorschemes/PuOr.h"
#include "colorschemes/PuRd.h"
#include "colorschemes/Purples.h"
#include "colorschemes/RdBu.h"
#include "colorschemes/RdGy.h"
#include "colorschemes/RdPu.h"
#include "colorschemes/RdYlBu.h"
#include "colorschemes/RdYlGn.h"
#include "colorschemes/Reds.h"
#include "colorschemes/Spectral.h"
#include "colorschemes/YlGnBu.h"
#include "colorschemes/YlGn.h"
#include "colorschemes/YlOrBr.h"
#include "colorschemes/YlOrRd.h"


std::map<std::string, const heatmap_colorscheme_t*> g_schemes = {
    {"b2w", heatmap_cs_b2w},
    {"b2w_opaque", heatmap_cs_b2w_opaque},
    {"w2b", heatmap_cs_w2b},
    {"w2b_opaque", heatmap_cs_w2b_opaque},
    {"Blues_discrete", heatmap_cs_Blues_discrete},
    {"Blues_soft", heatmap_cs_Blues_soft},
    {"Blues_mixed", heatmap_cs_Blues_mixed},
    {"Blues_mixed_exp", heatmap_cs_Blues_mixed_exp},
    {"BrBG_discrete", heatmap_cs_BrBG_discrete},
    {"BrBG_soft", heatmap_cs_BrBG_soft},
    {"BrBG_mixed", heatmap_cs_BrBG_mixed},
    {"BrBG_mixed_exp", heatmap_cs_BrBG_mixed_exp},
    {"BuGn_discrete", heatmap_cs_BuGn_discrete},
    {"BuGn_soft", heatmap_cs_BuGn_soft},
    {"BuGn_mixed", heatmap_cs_BuGn_mixed},
    {"BuGn_mixed_exp", heatmap_cs_BuGn_mixed_exp},
    {"BuPu_discrete", heatmap_cs_BuPu_discrete},
    {"BuPu_soft", heatmap_cs_BuPu_soft},
    {"BuPu_mixed", heatmap_cs_BuPu_mixed},
    {"BuPu_mixed_exp", heatmap_cs_BuPu_mixed_exp},
    {"GnBu_discrete", heatmap_cs_GnBu_discrete},
    {"GnBu_soft", heatmap_cs_GnBu_soft},
    {"GnBu_mixed", heatmap_cs_GnBu_mixed},
    {"GnBu_mixed_exp", heatmap_cs_GnBu_mixed_exp},
    {"Greens_discrete", heatmap_cs_Greens_discrete},
    {"Greens_soft", heatmap_cs_Greens_soft},
    {"Greens_mixed", heatmap_cs_Greens_mixed},
    {"Greens_mixed_exp", heatmap_cs_Greens_mixed_exp},
    {"Greys_discrete", heatmap_cs_Greys_discrete},
    {"Greys_soft", heatmap_cs_Greys_soft},
    {"Greys_mixed", heatmap_cs_Greys_mixed},
    {"Greys_mixed_exp", heatmap_cs_Greys_mixed_exp},
    {"Oranges_discrete", heatmap_cs_Oranges_discrete},
    {"Oranges_soft", heatmap_cs_Oranges_soft},
    {"Oranges_mixed", heatmap_cs_Oranges_mixed},
    {"Oranges_mixed_exp", heatmap_cs_Oranges_mixed_exp},
    {"OrRd_discrete", heatmap_cs_OrRd_discrete},
    {"OrRd_soft", heatmap_cs_OrRd_soft},
    {"OrRd_mixed", heatmap_cs_OrRd_mixed},
    {"OrRd_mixed_exp", heatmap_cs_OrRd_mixed_exp},
    {"PiYG_discrete", heatmap_cs_PiYG_discrete},
    {"PiYG_soft", heatmap_cs_PiYG_soft},
    {"PiYG_mixed", heatmap_cs_PiYG_mixed},
    {"PiYG_mixed_exp", heatmap_cs_PiYG_mixed_exp},
    {"PRGn_discrete", heatmap_cs_PRGn_discrete},
    {"PRGn_soft", heatmap_cs_PRGn_soft},
    {"PRGn_mixed", heatmap_cs_PRGn_mixed},
    {"PRGn_mixed_exp", heatmap_cs_PRGn_mixed_exp},
    {"PuBuGn_discrete", heatmap_cs_PuBuGn_discrete},
    {"PuBuGn_soft", heatmap_cs_PuBuGn_soft},
    {"PuBuGn_mixed", heatmap_cs_PuBuGn_mixed},
    {"PuBuGn_mixed_exp", heatmap_cs_PuBuGn_mixed_exp},
    {"PuBu_discrete", heatmap_cs_PuBu_discrete},
    {"PuBu_soft", heatmap_cs_PuBu_soft},
    {"PuBu_mixed", heatmap_cs_PuBu_mixed},
    {"PuBu_mixed_exp", heatmap_cs_PuBu_mixed_exp},
    {"PuOr_discrete", heatmap_cs_PuOr_discrete},
    {"PuOr_soft", heatmap_cs_PuOr_soft},
    {"PuOr_mixed", heatmap_cs_PuOr_mixed},
    {"PuOr_mixed_exp", heatmap_cs_PuOr_mixed_exp},
    {"PuRd_discrete", heatmap_cs_PuRd_discrete},
    {"PuRd_soft", heatmap_cs_PuRd_soft},
    {"PuRd_mixed", heatmap_cs_PuRd_mixed},
    {"PuRd_mixed_exp", heatmap_cs_PuRd_mixed_exp},
    {"Purples_discrete", heatmap_cs_Purples_discrete},
    {"Purples_soft", heatmap_cs_Purples_soft},
    {"Purples_mixed", heatmap_cs_Purples_mixed},
    {"Purples_mixed_exp", heatmap_cs_Purples_mixed_exp},
    {"RdBu_discrete", heatmap_cs_RdBu_discrete},
    {"RdBu_soft", heatmap_cs_RdBu_soft},
    {"RdBu_mixed", heatmap_cs_RdBu_mixed},
    {"RdBu_mixed_exp", heatmap_cs_RdBu_mixed_exp},
    {"RdGy_discrete", heatmap_cs_RdGy_discrete},
    {"RdGy_soft", heatmap_cs_RdGy_soft},
    {"RdGy_mixed", heatmap_cs_RdGy_mixed},
    {"RdGy_mixed_exp", heatmap_cs_RdGy_mixed_exp},
    {"RdPu_discrete", heatmap_cs_RdPu_discrete},
    {"RdPu_soft", heatmap_cs_RdPu_soft},
    {"RdPu_mixed", heatmap_cs_RdPu_mixed},
    {"RdPu_mixed_exp", heatmap_cs_RdPu_mixed_exp},
    {"RdYlBu_discrete", heatmap_cs_RdYlBu_discrete},
    {"RdYlBu_soft", heatmap_cs_RdYlBu_soft},
    {"RdYlBu_mixed", heatmap_cs_RdYlBu_mixed},
    {"RdYlBu_mixed_exp", heatmap_cs_RdYlBu_mixed_exp},
    {"RdYlGn_discrete", heatmap_cs_RdYlGn_discrete},
    {"RdYlGn_soft", heatmap_cs_RdYlGn_soft},
    {"RdYlGn_mixed", heatmap_cs_RdYlGn_mixed},
    {"RdYlGn_mixed_exp", heatmap_cs_RdYlGn_mixed_exp},
    {"Reds_discrete", heatmap_cs_Reds_discrete},
    {"Reds_soft", heatmap_cs_Reds_soft},
    {"Reds_mixed", heatmap_cs_Reds_mixed},
    {"Reds_mixed_exp", heatmap_cs_Reds_mixed_exp},
    {"Spectral_discrete", heatmap_cs_Spectral_discrete},
    {"Spectral_soft", heatmap_cs_Spectral_soft},
    {"Spectral_mixed", heatmap_cs_Spectral_mixed},
    {"Spectral_mixed_exp", heatmap_cs_Spectral_mixed_exp},
    {"YlGnBu_discrete", heatmap_cs_YlGnBu_discrete},
    {"YlGnBu_soft", heatmap_cs_YlGnBu_soft},
    {"YlGnBu_mixed", heatmap_cs_YlGnBu_mixed},
    {"YlGnBu_mixed_exp", heatmap_cs_YlGnBu_mixed_exp},
    {"YlGn_discrete", heatmap_cs_YlGn_discrete},
    {"YlGn_soft", heatmap_cs_YlGn_soft},
    {"YlGn_mixed", heatmap_cs_YlGn_mixed},
    {"YlGn_mixed_exp", heatmap_cs_YlGn_mixed_exp},
    {"YlOrBr_discrete", heatmap_cs_YlOrBr_discrete},
    {"YlOrBr_soft", heatmap_cs_YlOrBr_soft},
    {"YlOrBr_mixed", heatmap_cs_YlOrBr_mixed},
    {"YlOrBr_mixed_exp", heatmap_cs_YlOrBr_mixed_exp},
    {"YlOrRd_discrete", heatmap_cs_YlOrRd_discrete},
    {"YlOrRd_soft", heatmap_cs_YlOrRd_soft},
    {"YlOrRd_mixed", heatmap_cs_YlOrRd_mixed},
    {"YlOrRd_mixed_exp", heatmap_cs_YlOrRd_mixed_exp},
};

void reorder_strings(std::vector<std::string> &label_names, int* indices, int n) 
{ 
    std::vector<std::string> temp; 
    // indices's elements are the labels/indices of label_names
    // indices's positions (element index) are the new positions (element index)
    for (int i=0; i<n; i++) {
        temp.push_back(label_names.at(indices[i])); 
    }
  
    // Copy temp[] to col_names[] 
    for (int i=0; i<n; i++) 
    {  
       label_names.at(i) = temp.at(i); 
    } 
} 

template <typename T>
void reorder_matrix(T** &matrix, int* index, int num_data_rows, int num_data_cols, char axis) 
{ 
    std::string possible_axes = "rc";
    if (possible_axes.find(axis) == std::string::npos) {
        std::cerr << std::endl << "Axis value must be 'r' (row) or 'c' (column)" << std::endl;
        return;
    }

    double **temp = new double*[num_data_rows];
    for (int i = 0; i < num_data_rows; i++ )
    {
        temp[i] = new double[num_data_cols];
        for (int j=0; j< num_data_cols; j++) 
        {
            if (axis == 'r') {
                temp[i][j] = matrix[index[i]][j]; 
            }
            else {
                temp[i][j] = matrix[i][index[j]]; 
            }
        }
    }
  
    // Copy temp[] to col_names[] 
    for (int i=0; i<num_data_rows; i++) {
        for (int j=0; j< num_data_cols; j++) {
            matrix[i][j] = temp[i][j]; 
            // index[i] = i;
        }
    }

    // Deallocate memory
    for(int i = 0; i < num_data_rows; ++i) {
        free(temp[i]);
    }
    free(temp);

} 

int main(int argc, char* argv[])
{  
    clock_t start, mid, mid2, end;
    start = clock();
    
    if(argc == 2 && std::string(argv[1]) == "-l") {
        for(auto& scheme : g_schemes) {
            std::cout << "  " << scheme.first << std::endl;
        }
        return 0;
    }
#ifdef FIT_IMAGE
    if(argc < 8 || 9 < argc) {
#else      
    if(argc < 10 || 11 < argc) {
#endif 
        std::cerr << std::endl << "Invalid number of arguments!" << std::endl;
        std::cerr << "Usage:" << std::endl;
#ifdef FIT_IMAGE
        std::cerr << "  " << argv[0] << " image_width image_height num_data_cols num_data_rows distance_function linkage_function csv_data_filename [colorscheme]" << std::endl;
#else      
        std::cerr << "  " << argv[0] << " image_width image_height tile_ratio_x tile_ratio_y num_data_cols num_data_rows distance_funcition linkage_function csv_data_filename [colorscheme]" << std::endl;
#endif 
        std::cerr << std::endl;
        std::cerr << "  To get a list of available colorschemes, run" << std::endl;
        std::cerr << "  " << argv[0] << " -l" << std::endl;
        std::cerr << "  The default colorscheme is Spectral_mixed." << std::endl << std::endl;

        return 1;
    }
    
    const unsigned image_width = atoi(argv[1]), image_height = atoi(argv[2]);

#ifdef FIT_IMAGE
    const unsigned num_data_cols = atoi(argv[3]); 
    const unsigned num_data_rows = atoi(argv[4]);   
    const char distance_function = argv[5][0]; 
    const char linkage_function = argv[6][0];     
    const char* csv_path = argv[7];    

    if (image_width < num_data_cols || image_height < num_data_rows) {
        std::cerr << std::endl << "Image dimensions must be at least the dimensions of the data." << std::endl;
        std::cerr << "Specifically, the following must be true: " << std::endl;
        std::cerr << " image_width > num_data_cols" << std::endl;
        std::cerr << " image_height > num_data_rows" << std::endl << std::endl;

        return 1;
    }        

    // Calculate appropriate sizing for tile
    unsigned tile_width =  (int) (image_width / (num_data_cols));
    // std::cerr << "tile_width: " << tile_width << std::endl;
    unsigned tile_height = (int) (image_height / (num_data_rows));
    // std::cerr << "tile_height: " << tile_height << std::endl;

    if(argc >= 9 && g_schemes.find(argv[8]) == g_schemes.end()) {
        std::cerr << "Unknown colorscheme. Run " << argv[0] << " -l for a list of valid ones." << std::endl;
        return 1;
    }
    const heatmap_colorscheme_t* colorscheme = argc == 9 ? g_schemes[argv[8]] : heatmap_cs_default;

#else
    const unsigned tile_ratio_x = argc >= 10 ? atoi(argv[3]) : 1; 
    const unsigned tile_ratio_y = argc >= 10 ? atoi(argv[4]) : 1; 
    const unsigned num_data_cols = argc >= 10 ? atoi(argv[5]) : atoi(argv[3]); 
    const unsigned num_data_rows = argc >= 10 ? atoi(argv[6]) : atoi(argv[4]);     
    const char distance_function = argc >= 10 ? argv[7][0] : argv[5][0]; 
    const char linkage_function = argc >= 10 ? argv[8][0] : argv[6][0];     
    const char* csv_path = argv[9];

    if (image_width < (tile_ratio_x * num_data_cols) || image_height < (tile_ratio_y * num_data_rows)) {
        std::cerr << std::endl << "Image dimensions are not enough to accomodate tile dimensions and amount of data." << std::endl;
        std::cerr << "Specifically, the following must be true: " << std::endl;
        std::cerr << " image_width >= (tile_ratio_x * num_data_cols)" << std::endl;
        std::cerr << " image_height >= (tile_ratio_y * num_data_rows)" << std::endl << std::endl;

        return 1;
    }

    // Calculate appropiate sizing for tile
    unsigned max_x_scaling_factor =  (int) (image_width / (num_data_cols * tile_ratio_x));
    // std::cerr << "max_x_scaling_factor: " << max_x_scaling_factor << std::endl;
    unsigned max_y_scaling_factor = (int) (image_height / (num_data_rows * tile_ratio_y));
    // std::cerr << "max_y_scaling_factor: " << max_y_scaling_factor << std::endl;
    unsigned scaling_factor = min(max_x_scaling_factor,max_y_scaling_factor);
    // std::cerr << "scaling_factor: " << scaling_factor << std::endl;
    unsigned tile_width = scaling_factor * tile_ratio_x;
    // std::cerr << "tile_width: " << tile_width << std::endl;
    unsigned tile_height = scaling_factor * tile_ratio_y;
    // std::cerr << "tile_height: " << tile_height << std::endl;

    if(argc >= 11 && g_schemes.find(argv[10]) == g_schemes.end()) {
        std::cerr << "Unknown colorscheme. Run " << argv[0] << " -l for a list of valid ones." << std::endl;
        return 1;
    }
    const heatmap_colorscheme_t* colorscheme = argc == 11 ? g_schemes[argv[10]] : heatmap_cs_default;

#endif

    std::string possible_distance_functions = "cauxskeb";
    if (possible_distance_functions.find(distance_function) == std::string::npos) {
        std::cerr << std::endl << "Distance function given is not an option." << std::endl;
        std::cerr << "See readme for more info" << std::endl;

        return 1;
    }

    std::string possible_linkgae_functions = "smac";
    if (possible_linkgae_functions.find(linkage_function) == std::string::npos) {
        std::cerr << std::endl << "Linkage function given is not an option." << std::endl;
        std::cerr << "See readme for more info" << std::endl;

        return 1;
    }

    unsigned updated_image_width = tile_width * num_data_cols;
    // std::cerr << "updated_image_width: " << updated_image_width << std::endl;
    unsigned updated_image_height = tile_height * num_data_rows;
    // std::cerr << "updated_image_height: " << updated_image_height << std::endl;

    // Processing heatmap CSV data

    double **heatmap_data = new double*[num_data_rows];
    for(unsigned i = 0; i < num_data_rows; i++) {
        heatmap_data[i] = new double[num_data_cols];
    }

    // Allocate array of data's column names
    std::vector<std::string> col_names;
    // Allocate array of data's row names
    std::vector<std::string> row_names;

    // Mask for missing data, needed by the hierarchical clustering algorithm
    int **mask = new int*[num_data_rows];
    for(unsigned i = 0; i < num_data_rows; i++) {
        mask[i] = new int[num_data_cols];
    }

    std::ifstream ifile;
    ifile.open(csv_path);
    if(!ifile) {
        std::cerr << "CSV file does not exist at given filepath" << std::endl;
        return 1;
    }

    float weight;
    unsigned row_num = 0;
    unsigned col_num = 0;

    std::ifstream  data(csv_path);
    std::string line;
    while(std::getline(data,line))
    {
        std::stringstream  lineStream(line);
        std::string        cell;
        
        col_num = 0;

        while(std::getline(lineStream,cell,','))
        {
            if (row_num == 0 && col_num == 0) {
                col_num += 1;
                continue;
            }
            else if (row_num == 0){
                col_names.push_back(std::string(cell));
            }
            else if (col_num == 0){
                row_names.push_back(std::string(cell));

            }
            else {
                weight = std::stof(cell);
                heatmap_data[row_num-1][col_num-1] = weight;
                if (!weight) {
                    mask[row_num-1][col_num-1] = 0;
                }
                else {
                    mask[row_num-1][col_num-1] = 1;
                }
            }

            col_num += 1;
        }
        row_num += 1;
    }

    /* =========================== Hierarchical clustering =========================== */
    
    unsigned col_nnodes = num_data_cols-1;

    // Get dendrogram tree for column data
    double *col_weight = new double[num_data_cols];
    for(unsigned i = 0; i < num_data_cols; ++i) {
        col_weight[i] = 1.0;
    }
    Node* col_tree = treecluster(num_data_rows, num_data_cols, heatmap_data, mask, col_weight, 1, distance_function, linkage_function, 0); ///// Hardcoding Euclidean distance and Average linkage for now
    if (!col_tree)
    {
        std::cerr << ("treecluster routine failed due to insufficient memory") << std::endl;
        free(col_weight);
        return 1;
    }

    
    /*
    // Print tree data
    std::cerr << "Node     Item 1   Item 2    Distance\n" << std::endl;
    for(unsigned i=0; i<col_nnodes; i++){
        std::cerr << -i-1 << "     " << col_tree[i].left << "     " << col_tree[i].right << "     " << col_tree[i].distance << std::endl;
    }
    */
    

    // Sort column tree nodes
    int *col_sorted_indices = new int[num_data_cols];

    unsigned col_sorted_index = 0;
    for (unsigned i = 0; i < col_nnodes; i++) {
        if (col_tree[i].left >= 0) {
            col_sorted_indices[col_sorted_index] = col_tree[i].left;
            col_sorted_index++;
        }
        if (col_tree[i].right >= 0) {
            col_sorted_indices[col_sorted_index] = col_tree[i].right;
            col_sorted_index++;
        }
    }

    /*
    std::cerr << "BEFORE" << std::endl;
    for (unsigned i = 0; i < num_data_rows; i++)
    {
        for (unsigned j = 0; j < num_data_cols; j++)
        {
            std::cerr << heatmap_data[i][j] << ' ';
        }
        std::cerr << std::endl;
    }
    */

    // Reorder column labels
    reorder_strings(col_names, col_sorted_indices, num_data_cols);
    // Reorder heatmap columns
    reorder_matrix(heatmap_data, col_sorted_indices, num_data_rows, num_data_cols, 'c');
    // Reorder mask columns
    reorder_matrix(mask, col_sorted_indices, num_data_rows, num_data_cols, 'c');

    free(col_tree);
    free(col_weight);

    unsigned row_nnodes = num_data_rows-1;

    ///Get dendrogram for row data
    double *row_weight = new double[num_data_rows];
    for(unsigned i = 0; i < num_data_rows; ++i) {
        row_weight[i] = 1.0;
    }
    Node* row_tree = treecluster(num_data_rows, num_data_cols, heatmap_data, mask, row_weight, 0, distance_function, linkage_function, 0); ///// Hardcoding Euclidean distance and Average linkage for now

    if (!row_tree)
    {
        std::cerr << ("treecluster routine failed due to insufficient memory\n");
        free(row_weight);
        return 1;
    }

    // Print tree data
    std::cerr << "Node     Item 1   Item 2    Distance\n" << std::endl;
    for(unsigned i=0; i<row_nnodes; i++){
        std::cerr << i << "     " << row_tree[i].left << "     " << row_tree[i].right << "     " << row_tree[i].distance << std::endl;
    }

    // Sort row tree nodes
    int *row_sorted_indices = new int[num_data_rows];

    unsigned row_sorted_index = 0;
    for (unsigned i = 0; i < row_nnodes; i++) {
        if (row_tree[i].left >= 0) {
            row_sorted_indices[row_sorted_index] = row_tree[i].left;
            row_sorted_index++;
        }
        if (row_tree[i].right >= 0) {
            row_sorted_indices[row_sorted_index] = row_tree[i].right;
            row_sorted_index++;
        }
    }


    // Reorder row labels
    reorder_strings(row_names, row_sorted_indices, num_data_rows);
    // Reorder heatmap rows
    reorder_matrix(heatmap_data, row_sorted_indices, num_data_rows, num_data_cols, 'r');
    // Reorder mask
    reorder_matrix(mask, row_sorted_indices, num_data_rows, num_data_cols, 'r');

    for (unsigned i=0; i<num_data_rows; i++) {
        std::cerr << i << row_names[i] << std::endl;
    }

    // Free memory used during hierarchical clustering
    free(row_tree);
    free(row_weight);

    /* ============================ Generating heatmap pixels ============================ */
    
    heatmap_t* hm = heatmap_new(updated_image_width, updated_image_height);
    heatmap_stamp_t* tile = heatmap_stamp_gen(tile_width, tile_height);

    unsigned x, y;

    for (unsigned cur_row = 0; cur_row < num_data_rows; cur_row++) {

        y = ((cur_row+1) - 0.5) * tile_height;

        for (unsigned cur_col = 0; cur_col < num_data_cols; cur_col++) {
            x = ((cur_col+1) - 0.5) * tile_width;
            weight = heatmap_data[cur_row][cur_col];

            if(x < updated_image_width && y < updated_image_height) {
                heatmap_add_weighted_point_with_stamp(hm, x, y, weight, tile);
            } else {
                std::cerr << "Warning: Skipping out-of-bound input coordinate: (" << x << "," << y << ")." << std::endl;
            }
        }
    }

    heatmap_stamp_free(tile);

    std::vector<unsigned char> image(updated_image_width*updated_image_height*4);
    heatmap_render_to(hm, colorscheme, &image[0]);
    heatmap_free(hm);

    std::vector<unsigned char> png;
    mid = clock();
    double time_taken_1 = double(mid - start) / double(CLOCKS_PER_SEC);
    mid2 = clock();

    if(unsigned error = lodepng::encode(png, image, updated_image_width, updated_image_height)) {
        std::cerr << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
        return 1;
    }

    // lodepng::save_file(png, output_png_name);
    std::cout.write((char*)&png[0], png.size());


    // De-allocating data arrays
    for(unsigned
     i = 0; i < num_data_rows; ++i) {
        delete [] heatmap_data[i];
    }
    delete [] heatmap_data;

    end = clock();

    double time_taken_2 = double(end - mid2) / double(CLOCKS_PER_SEC);

    double time_taken_overall = time_taken_1 + time_taken_2;

    std::cerr << "Overall time taken by program is : " << std::fixed 
         << time_taken_overall << std::setprecision(5); 
    std::cerr << " sec " << std::endl; 

    std::cerr << "Without png generation (first half) : " << std::fixed 
         << time_taken_1 << std::setprecision(5); 
    std::cerr << " sec " << std::endl; 

    std::cerr << "PNG generation (second half) : " << std::fixed 
         << time_taken_2 << std::setprecision(5); 
    std::cerr << " sec " << std::endl; 
    
    return 0;
}
