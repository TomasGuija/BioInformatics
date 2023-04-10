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

/*
tuple_list classifier(tuple_list res, string sequence, int index, int n, string N_words[]){
    //Base case
    if(index >= sequence.length())
        return res;
    //If site has already been classified, go on to next site
    if(get<1>(res.at(index)) != 0)
        res = classifier(res, sequence, index++, n, N_words);
    else{
        //Create a new class and get all N_words containing this site
        get<1>(res.at(index)) = n;

        string words[N];
        int indexes[N];
        int cont = 0;
        for(int j = 0; j < N-1; j++){
            if(index-j >= 0 && index-j+N < sequence.length()){
                words[cont] = sequence.substr(index-j, N);
                indexes[cont] = j;
            }
            cont++;
        }

        //Look for instances of these same N_words in the sequence
        for(int j = 0; j < cont; j++){
            for(int k = 0; k < sequence.length() - (N-1); k++){
                if(words[j] == N_words[k]){
                    res = classifier(res, sequence, k+indexes[j], n, N_words);
                }
            }
        }
        
        res = classifier(res, sequence, index++, n++, N_words);
        return res;

    }
        
}
*/



tuple_list rewriteSequence(string sequence){

    //Clasify all sites to class 0 (unclassified)
    tuple_list res;
    for(int i = 0; i < sequence.length(); i++)
        res.push_back(tuple <char, int>(sequence[i], 0));

    //Get all N_words of the sequence
    int n_N_words = sequence.length() - (N-1);
    string N_words[n_N_words];

    for(int i = 0; i < n_N_words; i++){
        N_words[i] = sequence.substr(i, N);
    }

    //Classifying sites
    int n = 1;
    for(int i = 0; i < n_N_words; i++){
        string N_word = N_words[i];
        for(int j = i+1; j < n_N_words; j++){
            if(N_word == N_words[j]){
                //cout << N_word << "==" << N_words[j] << endl;

                for(int k = 0; k < N; k++){
                    int class_i = get<1>(res.at(i+k));
                    int class_j = get<1>(res.at(j+k));
                    if(class_i != 0 && class_j == 0)
                        class_j = class_i;
                    else if(class_i == 0 && class_j != 0)
                        class_i = class_j;
                    else if(class_i == 0 && class_j == 0){
                        class_j = n;
                        class_i = n;
                        n++;
                    }else{
                        for(int l = 0; l < sequence.length(); l++){
                            if(get<1>(res.at(l)) == class_j){
                                get<1>(res.at(l)) = class_i;
                            }
                        }
                        class_j = class_i;
                    }
                    get<1>(res.at(i+k)) = class_i;
                    get<1>(res.at(j+k)) = class_j;
                    //cout << N_word[k] << " " << class_i << endl;
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