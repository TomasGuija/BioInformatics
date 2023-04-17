#include <iostream>
#include <fstream>
#include <cstdlib>
#include <tuple>
#include <vector>
#include <string.h>
#define MAX_SEQUENCES 100
#define MAX
using namespace std;

typedef vector<tuple<char,int> > tuple_list;

ifstream file;
string sequences[MAX_SEQUENCES];
int N;

void readFile(char* path){
    file.open(path);
    if(!file.is_open()){
        exit(EXIT_FAILURE);
    }
    
    string line;
    int i = 0;
    while(getline(file, line)){
        if(line[0] != '>'){
            while(!line.empty()){
                sequences[i] += line;
                getline(file, line);
            }
            i++;
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
                    cout << N_word[k] << " " << class_i << endl;
                }
            }
        }
    }

    return res;
}

int main(){
    readFile("nef.fsa");
    N = 5;
    tuple_list res = rewriteSequence(sequences[0]);
    cout << "ORIGINAL SEQUENCE: " << endl << sequences[0] << endl;
    cout << endl << "REWRITTEN SEQUENCE:" << endl;
    for(int i = 0; i < res.size(); i++){
        cout << get<0>(res.at(i)) << get<1>(res.at(i)) << " ";
    }
    return 0;
}
