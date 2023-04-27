#include <iostream>
#include <fstream>
#include <cstdlib>
#include <tuple>
#include <vector>
#include <string.h>
#include <map>
#include <algorithm>
#include <sstream>



#define MAX_SEQUENCES 100

using namespace std;

typedef vector<tuple<char,int> > tuple_list;
typedef vector<vector<tuple<char, int> >> tuple_matrix;

ifstream file;
vector<string> sequences;

int N;
int n_sequences = 0;
map <int, vector<tuple<int, int> >> classes_map;


void readFile(char* path){
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
            //cout << "sequence  " << endl << sequences[n_sequences] << endl << endl; 
            n_sequences++;
        }
    }
}



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

tuple_matrix read_rewritten_sequences(){
    tuple_matrix res = {};

    file.open("rewritten.txt");
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

int main(){
    //readFile("nef.fsa");
    //N = 5;
    
    /*
    tuple_matrix res = rewriteSequenceList(sequences);
    ofstream file;
    file.open("rewritten.txt");
    for(int i = 0; i < res.size(); i++){
        for (int j = 0; j < res.at(i).size(); j++)
        {
            file << get<0>(res.at(i).at(j)) << get<1>(res.at(i).at(j)) << " ";
            if(j % 20 == 0 && j != 0) file << endl;
        }
        file << endl << endl;
    }
    file.close();
    */
    tuple_matrix rewritten_sequences;
    rewritten_sequences = read_rewritten_sequences();

    float d = distance(rewritten_sequences.at(0), rewritten_sequences.at(1));
    cout << "DISTANCE BETWEEN S0 AND S1: " << d << endl;
    return 0;
}
