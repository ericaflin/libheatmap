const addon = require('./build/Release/cclust');

// Input
var matrix_csvpath = "../../sample_data/smallGenesFile.csv";
var distance_function = "e";
var linkage_function = "a";
var axes = "b";

// Output after clustering
var output_object =  addon.ccluster(matrix_csvpath, distance_function, linkage_function, axes);
console.log(`\n The output is `, JSON.stringify(output_object));


