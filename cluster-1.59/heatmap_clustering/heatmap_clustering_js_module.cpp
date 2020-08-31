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
#include <vector>

extern "C" {
    #include "cluster.h"
}

class TreeNode {
    public:
        int NodeId;
        float Height;
        std::vector<int> Indices;
        std::vector<TreeNode> Children;
};

class Leaf : public TreeNode {
    public:
        std::string Label;
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
    (void)argc;
    clock_t start, mid, mid2, end;
    start = clock();

    std::string input = argv[1];

    /* =========================== Input Parsing (Destringifying) =========================== */

    // Args: heatmap_input, distance_function, linkage_function, axes (which dendrograms -- rows/cols/both)) -- parse these args from incoming string

    // Intaking data
    std::string heatmap_input;
    std::string _distance_function;
    std::string _linkage_function;
    std::string axes;

    std::stringstream  input_stream(input);
    std::string tmp_string = "";
    int arg_num = 0;

    for(std::string line; std::getline(input_stream, line); ) {
        tmp_string = std::string(line);
        if (arg_num == 0) {
            std::string heatmap_input(tmp_string.begin()+15, tmp_string.end()-2);  // also getting rid of the '[' and ']' at the beginning and the end -- makes parsing later easier
        }
        if (arg_num == 1) {
            std::string distance_function(tmp_string.begin()+18, tmp_string.end()-1);
        }
        if (arg_num == 2) {
            std::string linkage_function(tmp_string.begin()+17, tmp_string.end()-1);
        }
        if (arg_num == 3) {
            std::string axes(tmp_string.begin()+5, tmp_string.end()-1);
        }
        arg_num += 1;
    }    

    // Checking for valid inputs and setting up additional parameters
    char distance_function = _distance_function[0];
    char linkage_function = _linkage_function[0];
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

    bool col_dendro_flag = false;
    bool row_dendro_flag = false;
    if (axes == "c" || axes == "b") {
        col_dendro_flag = true;
    }
    if (axes == "r" || axes == "b") {
        row_dendro_flag = true;
    }

    int _second_bracket_idx = heatmap_input.find('[');
    int second_bracket_idx = heatmap_input.find('[', _second_bracket_idx+1);
    int num_data_cols = std::count(heatmap_input.begin(), heatmap_input.begin()+second_bracket_idx, ',');
    int num_data_rows = std::count(heatmap_input.begin(), heatmap_input.end(), '[') - 2; // subtract 2 to account for the extra '[' for a 2D array, and for the label row

    // Processing heatmap CSV data
    double **heatmap_data = new double*[num_data_rows];
    for(int i = 0; i < num_data_rows; i++) {
        heatmap_data[i] = new double[num_data_cols];
    }

    // Allocate array of data's column names
    std::vector<std::string> col_names;
    // Allocate array of data's row names
    std::vector<std::string> row_names;

    // Mask for missing data, needed by the hierarchical clustering algorithm
    int **mask = new int*[num_data_rows];
    for(int i = 0; i < num_data_rows; i++) {
        mask[i] = new int[num_data_cols];
    }

    float weight;
    int row_num = 0;
    int col_num = 0;

    std::stringstream  data(heatmap_input);
    std::string line;
    while(std::getline(data,line,']'))
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
                std::string _row_label = std::string(cell);
                std::string row_label(_row_label.begin()+2, _row_label.end()); // to get rid of the extra ',[' at the front
                row_names.push_back(row_label);

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

    mid = clock();

    /* =========================== Hierarchical clustering =========================== */
    std::map<int, TreeNode> col_node_dict;
    int cur_col_node_id = -1;
    std::map<int, TreeNode> row_node_dict;
    int cur_row_node_id = -1;

    if (col_dendro_flag) {
        int col_nnodes = num_data_cols-1;

        // Get dendrogram tree for column data
        double *col_weight = new double[num_data_cols];
        for(int i = 0; i < num_data_cols; ++i) {
            col_weight[i] = 1.0;
        }
        Node* col_tree = treecluster(num_data_rows, num_data_cols, heatmap_data, mask, col_weight, 1, distance_function, linkage_function, 0);
        if (!col_tree)
        {
            std::cerr << ("treecluster routine failed due to insufficient memory") << std::endl;
            free(col_weight);
            return 1;
        }

        // Print tree data
        std::cerr << "Node     Item 1   Item 2    Distance\n" << std::endl;
        for(int i=0; i<col_nnodes; i++){
            std::cerr << -i-1 << "     " << col_tree[i].left << "     " << col_tree[i].right << "     " << col_tree[i].distance << std::endl;
        }

        // Sort column tree nodes
        int *col_sorted_indices = new int[num_data_cols];

        int col_sorted_index = 0;
        for (int i = 0; i < col_nnodes; i++) {
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
        for (int i = 0; i < num_data_rows; i++)
        {
            for (int j = 0; j < num_data_cols; j++)
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

        // Creating TreeNode dict

        //// TODO: May need to remap leaf indices due to sorting heatmap data

        // Add leaves
        for(int i = 0; i < num_data_cols; ++i) {
            Leaf  new_leaf;
            new_leaf.NodeId = i;
            new_leaf.Height = 0;
            new_leaf.Indices = {(int)i};
            new_leaf.Children = {};
            new_leaf.Label = col_names[i];
            col_node_dict[i] = new_leaf;
        }

        // Add other nodes
        int cur_height = 1;
        
        for(int i=0; i<col_nnodes; i++){
            TreeNode new_tree_node;
            new_tree_node.NodeId = cur_col_node_id;
            new_tree_node.Height = cur_height;
            
            int left_child_id = col_tree[i].left;
            int right_child_id = col_tree[i].right;

            // Add all descendents to Children
            new_tree_node.Children.push_back(col_node_dict[left_child_id]);
            new_tree_node.Children.push_back(col_node_dict[right_child_id]);
            std::vector<TreeNode> left_child_children = col_node_dict[left_child_id].Children;
            std::vector<TreeNode> right_child_children = col_node_dict[right_child_id].Children;
            std::copy (left_child_children.begin(), left_child_children.end(), std::back_inserter(new_tree_node.Children));
            std::copy (right_child_children.begin(), right_child_children.end(), std::back_inserter(new_tree_node.Children));

            // Add itself and all descendents to Indices
            new_tree_node.Indices.push_back(cur_col_node_id);
            std::vector<int> left_child_indices = col_node_dict[left_child_id].Indices;
            std::vector<int> right_child_indices = col_node_dict[right_child_id].Indices;
            std::copy (left_child_indices.begin(), left_child_indices.end(), std::back_inserter(new_tree_node.Indices));
            std::copy (right_child_indices.begin(), right_child_indices.end(), std::back_inserter(new_tree_node.Indices));

            /*
            std::cerr << "New node: " << cur_col_node_id << std::endl;
            std::cerr << "children" << std::endl;
            for(int j=0; j < new_tree_node.Indices.size(); j++) {
                std::cerr << new_tree_node.Indices[j]
                << std::endl;
            }
            std::cerr << std::endl;
            */

            col_node_dict[cur_col_node_id] = new_tree_node;
            cur_col_node_id--;
            cur_height++;
        }

        free(col_tree);
        free(col_weight);

    }

    if (row_dendro_flag) {
        int row_nnodes = num_data_rows-1;

        double *row_weight = new double[num_data_rows];
        for(int i = 0; i < num_data_rows; ++i) {
            row_weight[i] = 1.0;
        }
        Node* row_tree = treecluster(num_data_rows, num_data_cols, heatmap_data, mask, row_weight, 0, distance_function, linkage_function, 0);

        if (!row_tree)
        {
            std::cerr << ("treecluster routine failed due to insufficient memory\n");
            free(row_weight);
            return 1;
        }

        // Sort row tree nodes
        int *row_sorted_indices = new int[num_data_rows];

        int row_sorted_index = 0;
        for (int i = 0; i < row_nnodes; i++) {
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

        // Creating TreeNode dict

        //// TODO: May need to remap leaf indices due to sorting heatmap data

        // Add leaves
        for(int i = 0; i < num_data_rows; ++i) {
            Leaf  new_leaf;
            new_leaf.NodeId = i;
            new_leaf.Height = 0;
            new_leaf.Indices = {(int)i};
            new_leaf.Children = {};
            new_leaf.Label = row_names[i];
            row_node_dict[i] = new_leaf;
        }

        // Add other nodes
        int cur_height = 1;
        
        for(int i=0; i<row_nnodes; i++){
            TreeNode new_tree_node;
            new_tree_node.NodeId = cur_row_node_id;
            new_tree_node.Height = cur_height;
            
            int left_child_id = row_tree[i].left;
            
            int right_child_id = row_tree[i].right;

            // Add all descendents to Children
            new_tree_node.Children.push_back(row_node_dict[left_child_id]);
            new_tree_node.Children.push_back(row_node_dict[right_child_id]);
            std::vector<TreeNode> left_child_children = row_node_dict[left_child_id].Children;
            std::vector<TreeNode> right_child_children = row_node_dict[right_child_id].Children;
            std::copy (left_child_children.begin(), left_child_children.end(), std::back_inserter(new_tree_node.Children));
            std::copy (right_child_children.begin(), right_child_children.end(), std::back_inserter(new_tree_node.Children));

            // Add itself and all descendents to Indices
            new_tree_node.Indices.push_back(cur_row_node_id);
            std::vector<int> left_child_indices = row_node_dict[left_child_id].Indices;
            std::vector<int> right_child_indices = row_node_dict[right_child_id].Indices;
            std::copy (left_child_indices.begin(), left_child_indices.end(), std::back_inserter(new_tree_node.Indices));
            std::copy (right_child_indices.begin(), right_child_indices.end(), std::back_inserter(new_tree_node.Indices));
            
            row_node_dict[cur_row_node_id] = new_tree_node;
            cur_row_node_id--;
            cur_height++;
        }

        // Free memory used during hierarchical clustering
        free(row_tree);
        free(row_weight);
    }

    mid2 = clock();

    /* =========================== Output Generation (Stringifying) =========================== */

    std::string output = "{heatmap:[";
    for (int i = 0; i < num_data_rows; i++ )
    {
        output.append("[");
        for (int j=0; j< num_data_cols; j++) 
        {
            output.append(std::to_string(heatmap_data[i][j]) + ",");
        }
        output.pop_back();
        output.append("]");
    }

    output.append("],\ncol_labels:[");
    for (int i = 0; i < num_data_cols; i++) {
        output.append(col_names.at(i) + ",");
    }
    output.pop_back();

    output.append("],\nrow_labels:[");
    for (int i = 0; i < num_data_rows; i++) {
        output.append(row_names.at(i) + ",");
    }
    output.pop_back();

    output.append("],\ncol_tree:");
    if (col_dendro_flag) {
        /// Stringify col tree and append to output
        // output.append(std::string(col_node_dict[cur_col_node_id+1]))
        std::cerr << "";
    }
    else {
        output.append("None");
    }

    output.append(",\nrow_tree:");
    if (row_dendro_flag) {
        /// Stringify row tree and append to output
        // output.append(std::string(row_node_dict[cur_row_node_id+1]))
        std::cerr << "";
    }
    else {
        output.append("None");
    }

    output.append("}");

    std::cerr << output << std::endl;

    // De-allocating data arrays
    for(int i = 0; i < num_data_rows; ++i) {
        delete [] heatmap_data[i];
    }
    delete [] heatmap_data;

    end = clock();
    double time_input = double(mid-start)/double(CLOCKS_PER_SEC);
    double time_clustering = double(mid2-mid)/double(CLOCKS_PER_SEC);
    double time_output = double(end-mid2)/double(CLOCKS_PER_SEC);
    double time_taken_overall = double(end-start)/double(CLOCKS_PER_SEC);

    std::cerr << "Input destringifying time : " << std::fixed 
         << time_input << std::setprecision(5); 
    std::cerr << " sec " << std::endl; 
    std::cerr << "Clustering time : " << std::fixed 
         << time_clustering << std::setprecision(5); 
    std::cerr << " sec " << std::endl; 
    std::cerr << "Output stringifying time : " << std::fixed 
         << time_output << std::setprecision(5); 
    std::cerr << " sec " << std::endl; 
    std::cerr << "Overall time taken by program is : " << std::fixed 
         << time_taken_overall << std::setprecision(5); 
    std::cerr << " sec " << std::endl; 

    return 0;
}
