/*

    INCLUDE DESCRIPTION

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



#define MAX_SEQUENCES 100
using namespace std;

//Defining some types that will be used later on
typedef vector<tuple<char,int>> tuple_list;
typedef vector<vector<tuple<char, int> >> tuple_matrix;
typedef vector<vector<float>> diss_matrix;

//Global variables
vector<string> sequences;
tuple_matrix rewritten_sequences;
vector<string> sequence_names;
int N;
int n_sequences = 0;
map <int, vector<tuple<int, int> >> classes_map;

//Reading original sequences
void readFile(char* path){
    ifstream file;
    file.open(path);
    if(!file.is_open()){
        exit(EXIT_FAILURE);
    }
    
    string line;
    while(getline(file, line)){
        sequences.push_back("");
        if(line[0] != '>'){
            while(!line.empty()){
                sequences[n_sequences] += line;
                getline(file, line);
            }
            n_sequences++;
        }else{
            sequence_names.push_back(line);
        }
    }
}

//Cluster structure for hierarchichal clustering
struct Cluster {
    int id;
    vector<int> members;
    float dist;
    Cluster *left;
    Cluster *right;
};

//Merging 2 clusters together
Cluster *merge_clusters(Cluster *c1, Cluster *c2, float dist) {
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


//Calculating distance between 2 sequences
float distance(tuple_list s1, tuple_list s2){
    vector<int> compared_classes = {};
    float len, sum = 0;
    tuple_list shortest;

    if(s1.size() > s2.size()) shortest = s2;
    else shortest = s1;

    len = shortest.size();

    for(int i  = 0; i < len; i++){
        tuple<char, int> x = shortest.at(i);
        if(count(compared_classes.begin(), compared_classes.end(), get<1>(x))) continue;
        compared_classes.push_back(get<1>(x));
        int u = count(s1.begin(), s1.end(), x);
        int v = count(s2.begin(), s2.end(), x);
        int min = std::min(u,v);
        sum += min;
    }
    return (1-sum/len);
}

//Getting a sequence representative of a cluster.
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
        //Get the frequency of the different classes for position i
        for(int j = 0; j < cluster.members.size(); j++){
            int idx = cluster.members.at(j);
            tuple<char, int> class_i = rewritten_sequences.at(idx).at(i);
            if(frequency_map.count(class_i) > 0){ frequency_map.at(class_i) = frequency_map.at(class_i) + 1;}
            else{frequency_map[class_i] = 1;}
        }
        //Get the most frequent class
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
    for (size_t i = 0; i < n; ++i) {
        clusters[i] = new Cluster;
        clusters[i]->id = i;
        clusters[i]->members.push_back(i);
        clusters[i]->dist = 0;
        clusters[i]->left = nullptr;
        clusters[i]->right = nullptr;
    }
 
    while (clusters.size() > 1) {
        cout << "Looping" << endl;
        float min_dist = numeric_limits<float>::max();
        size_t min_i = 0, min_j = 0;

        for (size_t i = 0; i < clusters.size(); ++i) {
            for (size_t j = i + 1; j < clusters.size(); ++j) {
                float dist;
                if (clusters.at(i)->members.size() == 1 && clusters.at(j)->members.size()){
                    int idx_i = clusters.at(i)->members.at(0);
                    int idx_j = clusters.at(j)->members.at(0);
                    dist = dissimilarity_matrix.at(idx_i).at(idx_j);
                }else{
                    tuple_list new_chain_i = majority_chain(*clusters.at(i));
                    tuple_list new_chain_j = majority_chain(*clusters.at(j));
                    dist = distance(new_chain_i, new_chain_j);
                }
                /*
                int i_id = clusters[i]->members[0];
                int j_id = clusters[j]->members[0];
                float dist = dissimilarity_matrix[i_id][j_id];
                */
                if (dist < min_dist) {
                    min_dist = dist;
                    min_i = i;
                    min_j = j;
                }
                
            }
        }

        Cluster *merged_cluster = merge_clusters(clusters[min_i], clusters[min_j], min_dist);
        clusters[min_i] = merged_cluster;
        clusters.erase(clusters.begin() + min_j);
    }
 
    return clusters[0];
}


//Visualizing hierarchical clustering tree
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

//Rewritting sequences using N-local decoding
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

    //Classifying sites
    int n = 1;
    for(int i = 0; i < n_N_words; i++){
        tuple<string, tuple<int, int>> N_word_i = N_words[i];

        for(int j = i+1; j < n_N_words; j++){
            tuple<string, tuple<int,int>> N_word_j = N_words[j];

            if(get<0>(N_word_i) == get<0>(N_word_j)){
                int sequence_i = get<0>(get<1>(N_word_i));
                int sequence_j = get<0>(get<1>(N_word_j)); 
                int index_i = get<1>(get<1>(N_word_i));
                int index_j = get<1>(get<1>(N_word_j)); 
                
                for(int k = 0; k < N; k++){
                    //Here I get the sequence and index of each site to be classified
                    int class_i = get<1>(res.at(sequence_i).at(index_i+k));
                    int class_j = get<1>(res.at(sequence_j).at(index_j+k));

                    if(class_i != 0 && class_j == 0){
                        class_j = class_i;
                        classes_map.at(class_i).push_back({sequence_j, index_j+k});
                    }
                    else if(class_i == 0 && class_j != 0){
                        class_i = class_j;
                        classes_map.at(class_j).push_back({sequence_i, index_i+k});
                    }
                    else if(class_i == 0 && class_j == 0){
                        class_j = n;
                        class_i = n;
                        classes_map.insert({n, {{sequence_i, index_i+k}, {sequence_j, index_j+k}}});
                        n++;
                    }else if(class_i != class_j){
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
    for(int i = 0; i < len; i++){
        vector<float> v(len, 0);
        res.push_back(v);
        for(int j = i+1; j < len; j++){
            res.at(i).at(j) = distance(sequences.at(i), sequences.at(j));
        }
    }
    return res;
}

//Reading rewritten sequences
tuple_matrix read_rewritten_sequences(){
    ifstream file;

    tuple_matrix res = {};

    file.open("rewritten1.txt");
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

    return res;
}

//Reading dissimilarity matrix
diss_matrix read_diss_matrix(){
    ifstream file;

    diss_matrix matrix;
    file.open("dissmatrix.txt");

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

    return matrix;
}

int main(){
    //Reading the sequences to save the names
    readFile("nef.fsa");
    /*
    N = 15;
    
    tuple_matrix res = rewriteSequenceList(sequences);
    ofstream file;
    file.open("rewritten1.txt");
    for(int i = 0; i < res.size(); i++){
        int j;
        for (j = 0; j < res.at(i).size(); j++)
        {
            if(j % 20 == 0 && j != 0){
                file << get<0>(res.at(i).at(j)) << get<1>(res.at(i).at(j));
                file << endl;
            }else{
                file << get<0>(res.at(i).at(j)) << get<1>(res.at(i).at(j)) << " ";
            }

        }
        if((j-1) % 20 == 0){
            file << endl;
        }else{
            file << endl << endl;
        }
    }
    file.close();*/
    
    //We read the rewritten sequences
    rewritten_sequences = read_rewritten_sequences();

    
    
    //float d = distance(rewritten_sequences.at(0), rewritten_sequences.at(3));

    /*
    ofstream file;
    file.open("dissmatrix.txt");
    diss_matrix diss =  ComputeDissimilarityMatrix(rewritten_sequences);
    for(int i = 0; i < diss.size(); i++){
        for(int j = 0; j < diss.size(); j++){
            file << fixed << setprecision(5) << diss.at(i).at(j) << " ";
        }
        file << endl;
    }
    file.close();
    */

    //Reading the already computed dissimilarity matrix
    diss_matrix matrix = read_diss_matrix();

    /*
    for(int i = 0; i < matrix.size(); i++){
        for (int j = 0; j < matrix.size(); j++){
            cout << matrix.at(i).at(j) << " ";
        }
        cout << endl;
    }
    */

    //Performing hierarchical clustering
    Cluster *root = hierarchical_clustering(matrix);

    ofstream file;
    file.open("dendogram.txt");

    print_cluster(root, file);

    file.close();


    return 0;
    
}
