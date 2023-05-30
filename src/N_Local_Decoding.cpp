/*

    This script contains all necessary methods used along the project. A description of each method can be found above them.
    Briefly, with the developed methods, we were able to rewrite a set of sequences unsing N-local decoding, and build a dissimilarity
    matrix along with a dendogram built from hierarchical clustering. 
    To assess the reliability of the built trees, there is also some code to generate bootstrap samples from the original dataset, 
    which can be used to build trees that will be merged together using consensus.exe from the PHYLIP package.

*/

//Libraries used in the work
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <tuple>
#include <vector>
#include <string.h>
#include <map>
#include <algorithm>
#include <sstream>
#include <bits/stdc++.h>
#include <algorithm>

using namespace std;

//Defining some types that will be used later on
//This type will store the rewritten sequences
typedef vector<tuple<char,int>> tuple_list;
//This type will store a vector of rewritten sequences
typedef vector<vector<tuple<char, int>>> tuple_matrix;
//Dissimilarity matrix
typedef vector<vector<float>> diss_matrix;


//Global variables
vector<string> sequences;
tuple_matrix rewritten_sequences;
vector<string> sequence_names;
int n_sequences = 0;
map <int, vector<tuple<int, int>>> classes_map;

//The only parameter needed is N, for the word size used for the N-local decoding algorithm
int N = 15;

//This method can be used to read the original sequences. The output will be stored in the golbal variable "sequences"
void readFile(char* path){
    ifstream file;
    file.open(path);
    if(!file.is_open()){
        exit(EXIT_FAILURE);
    }
    
    string line;
    while(getline(file, line)){
        if(line[0] != '>'){
            if(!line.empty()) sequences.push_back("");
            while(!line.empty()){
                sequences.at(n_sequences) += line;
                getline(file, line);
            }
            n_sequences++;
        }else{
            sequence_names.push_back(line);
        }
    }
    file.close();
}

//Cluster structure for hierarchichal clustering algorithm
struct Cluster {
    int id;
    vector<int> members;
    float dist;
    Cluster *left;
    Cluster *right;
};

//Merging 2 clusters together
Cluster *merge_clusters(Cluster *c1, Cluster *c2, float dist) {
    //Merged clusters will have an id = -1. We take all members from both clusters, and keep the distance between those 2 clusters
    Cluster *new_cluster = new Cluster;
    new_cluster->id = -1;
    new_cluster->members.reserve(c1->members.size() + c2->members.size());
    new_cluster->members.insert(new_cluster->members.end(), c1->members.begin(), c1->members.end());
    new_cluster->members.insert(new_cluster->members.end(), c2->members.begin(), c2->members.end());
    new_cluster->dist = dist;
    new_cluster->left = c1;
    new_cluster->right = c2;
    return new_cluster;
}


//Calculating distance between 2 sequences. The way this is done is just as described in the original paper.
float distance(tuple_list s1, tuple_list s2){
    vector<int> compared_classes = {};
    float len, sum = 0;
    tuple_list shortest;

    //We will iterate over the shortest sequence of the 2
    if(s1.size() > s2.size()) shortest = s2;
    else shortest = s1;

    len = shortest.size();

    for(int i  = 0; i < len; i++){
        tuple<char, int> x = shortest.at(i);
        //If the class at position i has already been taken into account for the distance, go on to the next one
        if(count(compared_classes.begin(), compared_classes.end(), get<1>(x))) continue;
        //Otherwise, we take the minimun number of instances of class i in both sequences and add it to the sum
        compared_classes.push_back(get<1>(x));
        int u = count(s1.begin(), s1.end(), x);
        int v = count(s2.begin(), s2.end(), x);
        int min = std::min(u,v);
        sum += min;
    }

    //The sum will give an idea of how close both sequences are from each other. To normalize this value, we must divide it between the 
    //length of the shortest chain. To transform it into a distance measure, we just substract it from 1
    return (1-sum/len);
}

//Getting a sequence representative of a cluster, using the majority rule
tuple_list majority_chain(Cluster cluster){
    tuple_list majority;
    int min_size = std::numeric_limits<int>::max();

    //First we get the shortest sequence from the cluster
    for(int i = 0; i < cluster.members.size(); i++){
        int idx = cluster.members.at(i);
        if(rewritten_sequences.at(idx).size() < min_size) {min_size = rewritten_sequences.at(idx).size();}
    }


    //We iterate over the minimun chain size
    for(int i = 0; i < min_size; i++){
        map<tuple<char, int>, int> frequency_map;
        //Get the frequency of the different classes for position i. The frequency of each class is stored in the frequency_map
        for(int j = 0; j < cluster.members.size(); j++){
            int idx = cluster.members.at(j);
            tuple<char, int> class_i = rewritten_sequences.at(idx).at(i);
            if(frequency_map.count(class_i) > 0){ frequency_map.at(class_i) = frequency_map.at(class_i) + 1;}
            else{frequency_map[class_i] = 1;}
        }
        //Get the most frequent class, and add it to the ith position of the majority sequence
        int max_frequency = 0;
        tuple<char, int> most_frequent;
        for(auto k : frequency_map){
            if(k.second > max_frequency){
                max_frequency = k.second;
                most_frequent = k.first;
            }
        }
        majority.push_back(most_frequent);
    }

    return majority;

}

//Hierarchical clustering algorithm
Cluster *hierarchical_clustering(const vector<vector<float>> &dissimilarity_matrix) {
    size_t n = dissimilarity_matrix.size();
    vector<Cluster*> clusters(n);

    // First, we will start with a cluster for every single sequence. Then, according to the distance between the sequences, we will 
    // start merging them together
    for (size_t i = 0; i < n; ++i) {
        clusters[i] = new Cluster;
        clusters[i]->id = i;
        clusters[i]->members.push_back(i);
        clusters[i]->dist = 0;
        clusters[i]->left = nullptr;
        clusters[i]->right = nullptr;
    }
 
    //We stop when there is only one big cluster left with all sequences
    while (clusters.size() > 1) {
        cout << "Looping" << endl;

        //For the next merge, we will loop through every possible pair of clusters, and take the pair that is the closest to each other
        float min_dist = numeric_limits<float>::max();
        size_t min_i = 0, min_j = 0;

        for (size_t i = 0; i < clusters.size(); ++i) {
            for (size_t j = i + 1; j < clusters.size(); ++j) {
                float dist;
                //If the clusters have a size of 1, then we compute their distance by simply looking at the dissimilarity matrix
                if (clusters.at(i)->members.size() == 1 && clusters.at(j)->members.size() == 1){
                    int idx_i = clusters.at(i)->members.at(0);
                    int idx_j = clusters.at(j)->members.at(0);
                    dist = dissimilarity_matrix.at(idx_i).at(idx_j);
                }else{
                    //Otherwise, we calculate the majority chain for each of them and compute the distance
                    tuple_list new_chain_i = majority_chain(*clusters.at(i));
                    tuple_list new_chain_j = majority_chain(*clusters.at(j));
                    dist = distance(new_chain_i, new_chain_j);
                }

                // Keep a reference to the 2 closest clusters so far
                if (dist < min_dist) {
                    min_dist = dist;
                    min_i = i;
                    min_j = j;
                }
                
            }
        }

        //We merge both clusters the position of min_i, and erase the other cluster.
        Cluster *merged_cluster = merge_clusters(clusters[min_i], clusters[min_j], min_dist);
        clusters[min_i] = merged_cluster;
        clusters.erase(clusters.begin() + min_j);
    }
 
    return clusters[0];
}

//Visualizing hierarchical clustering tree. This is the first method we used to do so, although it's not the easiest format to 
//interprete. However, we decided to keep it in the script. Better ways of visualizing the results will be shown later.
void print_cluster(const Cluster *cluster, ofstream &file, int level = 0) {

    if (!cluster) return;

    for (int i = 0; i < level; ++i) {
        file << "  ";
    }
    
    file << "Cluster " << (cluster->id != -1 ? to_string(cluster->id) : "(merged)");
    file << " - Members: ";
    for (int member : cluster->members) {
        file << sequence_names.at(member) << " ";
    }
    file << " - Distance: " << cluster->dist << endl;
    print_cluster(cluster->left, file, level + 1);
    print_cluster(cluster->right, file, level + 1);
}

//Rewritting sequences using N-local decoding. No longer used. rewriteSequenceList is the more developed version, adapted to work 
//on several sequences at the same time.
tuple_list rewriteSequence(string sequence){

    //Clasify all sites to class 0 (unclassified)
    tuple_list res;
    int n_N_words = sequence.length() - (N-1);
    string N_words[n_N_words];
    for(int i = 0; i < sequence.length(); i++){
        res.push_back(tuple <char, int>(sequence[i], 0)); //Initialising the res <A,0>, <C,0>...
        //Storing every single N_word in an array
        if(i < n_N_words)
            N_words[i] = sequence.substr(i, N);

    }


    //Classifying sites
    int n = 1;
    //Bigger loop goes through every single N word
    for(int i = 0; i < n_N_words; i++){
        string N_word = N_words[i];
        //Second loop goes through all N words starting from i+1
        for(int j = i+1; j < n_N_words; j++){
            if(N_word == N_words[j]){
                //This loop goes through all letters in the matching N-words
                for(int k = 0; k < N; k++){
                    //Get the current class of each site
                    int class_i = get<1>(res.at(i+k));
                    int class_j = get<1>(res.at(j+k));
                    cout << "Letras: " << get<0>(res.at(i+k)) << "  " << get<0>(res.at(j+k)) << endl;

                    //In case one of them has been classifed and the other hasn't, we match the class of the unclassified one to the first one
                    if(class_i != 0 && class_j == 0)
                        class_j = class_i;
                    else if(class_i == 0 && class_j != 0)
                        class_i = class_j;
                    //If none have been classified, we create a new class for both of them
                    else if(class_i == 0 && class_j == 0){
                        class_j = n;
                        class_i = n;
                        n++;
                    }else{ //When both have been classified
                        //We go through every single site in the sequence, and we will preserve the class of the site i
                        //Whenever we find a site classified as class_j, we make its class equal to class_i
                        for(int l = 0; l < sequence.length(); l++){
                            if(get<1>(res.at(l)) == class_j){
                                get<1>(res.at(l)) = class_i;
                            }
                        }
                        class_j = class_i;
                    }
                    //Update the values of the classes in res
                    get<1>(res.at(i+k)) = class_i;
                    get<1>(res.at(j+k)) = class_j;
                }
            }
        }
    }

    return res;
}

//Rewritting all sequences using N-local decoding
tuple_matrix rewriteSequenceList(vector<string> sequences){

    //We have a matrix for storing the results, each row corresponding to the sites of a sequence
    tuple_matrix res;
    int n_N_words;
    //This vector stores al N_words, along with the sequence they are in, and their index within this sequence
    vector<tuple<string, tuple<int, int>>> N_words;

    //The bigger loop goes through all sequences
    for(int i = 0; i < sequences.size(); i++){
        n_N_words = sequences[i].length() - (N-1); 
        res.push_back({});
        //The smaller loop goes through every character in each sequence
        for(int j = 0; j < sequences[i].length(); j++){
            //Initialising res, with each site classified as class 0 (unclassified)
            res.at(i).push_back(tuple <char, int>(sequences[i][j], 0));
            //Storing N_words
            if(j < n_N_words)
                N_words.push_back({sequences[i].substr(j, N), {i,j}});
        }
    }

    //Total number of N_words in all sequences
    n_N_words = N_words.size();

    //Classifying sites, iterating over the set of all N_words
    int n = 1;
    for(int i = 0; i < n_N_words; i++){
        tuple<string, tuple<int, int>> N_word_i = N_words[i];

        for(int j = i+1; j < n_N_words; j++){
            tuple<string, tuple<int,int>> N_word_j = N_words[j];

            //For every pair of N_words, we check wether they are the same
            if(get<0>(N_word_i) == get<0>(N_word_j)){
                //If so get the index of their sequence and their index within those sequences
                int sequence_i = get<0>(get<1>(N_word_i));
                int sequence_j = get<0>(get<1>(N_word_j)); 
                int index_i = get<1>(get<1>(N_word_i));
                int index_j = get<1>(get<1>(N_word_j)); 
                
                //We loop over every site in the N_word
                for(int k = 0; k < N; k++){
                    //Here I get the class of the matching sites in both sequences
                    int class_i = get<1>(res.at(sequence_i).at(index_i+k));
                    int class_j = get<1>(res.at(sequence_j).at(index_j+k));

                    //If only one of them has already been classified, then classify the unclassified one 
                    if(class_i != 0 && class_j == 0){
                        class_j = class_i;
                        classes_map.at(class_i).push_back({sequence_j, index_j+k});
                    }
                    else if(class_i == 0 && class_j != 0){
                        class_i = class_j;
                        classes_map.at(class_j).push_back({sequence_i, index_i+k});
                    }
                    //If none of them have been classified, then create a new class
                    else if(class_i == 0 && class_j == 0){
                        class_j = n;
                        class_i = n;
                        classes_map.insert({n, {{sequence_i, index_i+k}, {sequence_j, index_j+k}}});
                        n++;
                    }
                    //If both have already been classified, we keep one class and discard the other
                    //In this case we keep class_i. Therefore, we look for all sites classified as class_j, and change them to be 
                    //part of class_i
                    else if(class_i != class_j){
                        vector<tuple<int, int>> to_change;
                        to_change = classes_map.at(class_j);
                        for(int l = 0; l < to_change.size(); l++){
                            int sequence = get<0>(to_change.at(l));
                            int index = get<1>(to_change.at(l));
                            get<1>(res.at(sequence).at(index)) = class_i;
                            classes_map.at(class_i).push_back(to_change.at(l));
                        }
                        classes_map.erase(class_j);
                        class_j = class_i;
                    }
                    get<1>(res.at(sequence_i).at(index_i+k)) = class_i;
                    get<1>(res.at(sequence_j).at(index_j+k)) = class_j;
                }
            }
        }
    }

    return res;
}

//Computing de dissimilarity matrix over all sequences
diss_matrix ComputeDissimilarityMatrix(tuple_matrix sequences){
    diss_matrix res;
    int len = sequences.size();
    //We just do a double for loop to iterate over all pairs of sequences and calculate their distance
    for(int i = 0; i < len; i++){
        vector<float> v(len, 0);
        res.push_back(v);
        for(int j = i+1; j < len; j++){
             float d = distance(sequences.at(i), sequences.at(j));
            res.at(i).at(j) = d;
        }
    }
    return res;
}

//Reading rewritten sequences from a file
tuple_matrix read_rewritten_sequences(string path){
    ifstream file;
    tuple_matrix res = {};

    file.open(path);
    if(!file.is_open()) exit(EXIT_FAILURE);
    
    string line;
    n_sequences = 0;
    while(getline(file, line)){
        res.push_back({});
        while(!line.empty()){
            string tmp;
            stringstream ss(line);
            while(getline(ss, tmp, ' ')){
                char letter = tmp[0];
                string numbers = tmp.substr(1, tmp.size()-1);
                int n =  stoi(numbers);
                res.at(n_sequences).push_back({letter,n});
            }
            getline(file, line);
        }
        n_sequences++;
    }

    file.close();
    return res;
}

//Reading dissimilarity matrix from a file
diss_matrix read_diss_matrix(string path){
    ifstream file;

    diss_matrix matrix;
    file.open(path);

    if(!file.is_open()) exit(EXIT_FAILURE);
    int i = 0;
    string line;
    while(getline(file, line)){
        matrix.push_back({});
        string tmp;
        stringstream ss(line);
        //int j = 0;
        while(getline(ss, tmp, ' ')){
            float n = stof(tmp);
            matrix.at(i).push_back(n);
        }
        i++;
    }
    file.close();
    return matrix;
}


// function to generate a single replicate sequence
tuple_list generateReplicate(const tuple_list& original) {
    tuple_list replicate;

    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_int_distribution<> distrib(0, original.size() - 1);

    for (int i = 0; i < original.size(); ++i) {
        int random_index = distrib(gen);
        replicate.push_back(original[random_index]);
    }

    return replicate;
}

// function to generate replicates for all sequences in the original set
tuple_matrix generateReplicates(tuple_matrix& originalSet) {
    tuple_matrix replicates;

    for (const auto& sequence : originalSet) {
        replicates.push_back(generateReplicate(sequence));
    }

    return replicates;
}


// Generate a representation of the hierarchical clustering trees in the Newick format. This will later be useful for representing the 
//trees and using the Consensus function of the PHYLIP package
string generateNewick(const Cluster* cluster) {
    if (cluster->left == nullptr && cluster->right == nullptr) {
        // Leaf node
        return sequence_names.at(cluster->id);
    } else {
        // Internal node
        string leftStr = generateNewick(cluster->left);
        string rightStr = generateNewick(cluster->right);
        return "(" + leftStr + "," + rightStr + ")";
    }
}

// Function to write trees in Newick format to a text file
void writeTreesToNewickFile(vector<Cluster*>& trees,const std::string& filename) {

    ofstream outputFile(filename);

    if (!outputFile.is_open()) {
        cout << "Error opening file";
        return;
    }

    for (const auto& tree : trees) {
        string newickStr = generateNewick(tree);
        outputFile << newickStr << ";" << endl;
    }

    outputFile.close();
    cout << "Trees written to file: " << filename << endl;
}

void WriteDissMatrix(diss_matrix matrix, string path){
    ofstream file1;
    file1.open(path);
    for(int i = 0; i < matrix.size(); i++){
        for(int j = 0; j < matrix.size(); j++){
            file1 << fixed << setprecision(5) << matrix.at(i).at(j) << " ";
        }
        file1 << endl;
    }
    file1.close();
}

void WriteRewrittenSequences(tuple_matrix rewritten_sequences, string path){
    ofstream file;
    file.open(path);
    for(int i = 0; i < rewritten_sequences.size(); i++){
        int j;
        if (rewritten_sequences.at(i).empty()){
            break;
        }
        for (j = 0; j < rewritten_sequences.at(i).size(); j++)
        {
            if(j % 20 == 0 && j != 0){
                file << get<0>(rewritten_sequences.at(i).at(j)) << get<1>(rewritten_sequences.at(i).at(j));
                file << endl;
                cout << get<0>(rewritten_sequences.at(i).at(j)) << get<1>(rewritten_sequences.at(i).at(j));
                cout << endl;
            }else{
                file << get<0>(rewritten_sequences.at(i).at(j)) << get<1>(rewritten_sequences.at(i).at(j)) << " ";
                cout << get<0>(rewritten_sequences.at(i).at(j)) << get<1>(rewritten_sequences.at(i).at(j)) << " ";

            }

        }
        if((j-1) % 20 == 0){
            file << endl;
        }else{
            file << endl << endl;
        }
    }
    file.close();
}


int main(){
    //Reading the sequences from the original dataset, storing the names in sequence_names
    readFile("nef.fsa");

    /*
    for(int i = 0; i < sequences.size(); i++){
        cout << sequences.at(i) << endl;
    }
    */
    
    //The following commented code will rewrite all sequences read from the readFile() method, and write them into a txt file.
    //The path of the output file can be changed as one lines, as well as the value for the parameter N.
    //This has already been done some sets of sequences using a value of N=15. The computation takes time, but for testing the
    //code you can uncomment the lines and run this
    
    /*
    rewritten_sequences = rewriteSequenceList(sequences);
    WriteRewrittenSequences(rewritten_sequences, "rewritten_66_sequences.txt");
    */
    

    //We read the rewritten sequences
    rewritten_sequences = read_rewritten_sequences("../data/Rewritten sequences/rewritten_66_sequences.txt");

    //This code will compute the dissimilarity matrix over a set of rewritten sequences (always use the original set of sequences,
    //not the bootstrap replicates) and save them to a txt file.
    
    /*
    diss_matrix matrix =  ComputeDissimilarityMatrix(rewritten_sequences);
    WriteDissMatrix(matrix, "dissmatrix_66_sequences.txt")
    */
    
    //This part of the code will generate a certain number of bootstrap replicates. The process of hierarchical clustering takes 
    //over 300 seconds for each set of 66 sequences, so a big number of replicates is not recommended unless you have a lot of 
    //time and/or computational power. Once again, several results have already been computed, so it's not necessary to run everything
    //again. Setting N_replicates = 0, this could all be run but just using the original dataset, which can be done in feasible time.
    
    /*
    vector<tuple_matrix> bootstrap_replicates;
    int N_replicates = 10;
    for(int i = 0; i < N_replicates; i++){
        bootstrap_replicates.push_back(generateReplicates(rewritten_sequences));
    }
    */
    

    //Computing dissimilarity matrices for the bootstrap samples
    
    /*
    vector<diss_matrix> bootstrap_diss;
    for(int i = 0; i < N_replicates; i++){
        diss_matrix m = ComputeDissimilarityMatrix(bootstrap_replicates.at(i));
        bootstrap_diss.push_back(m);
    }
    */

    //Reading the already computed dissimilarity matrix
    diss_matrix matrix = read_diss_matrix("../data/Dissimilarity Matrices/dissmatrix_66_sequences.txt");

    //Performing hierarchical clustering, over the original set of sequences as well as over the bootstrap replicates.
    //All the results are stored in the vector trees.
    
    
    Cluster *root = hierarchical_clustering(matrix);
    vector<Cluster*> trees;
    cout << "Finished clustering" << endl;
    trees.push_back(root);

    //If computing bootstrap trees, uncomment this part to add the trees 
    
    /*
    for(int i = 0; i < N_replicates; i++){
        Cluster * clust = hierarchical_clustering(bootstrap_diss.at(i));
        trees.push_back(clust);
    }
    */
    
    //Final part of the code, once we have the trees we want to run Consensus on, we save them to a txt file in the Newick format
    writeTreesToNewickFile(trees, "newick_trees_66_sequences.txt");
    

    return 0;
    
}
