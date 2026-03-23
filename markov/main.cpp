
#include<iostream>
#include<vector>
#include<algorithm>
#include<fstream>
#include<sstream>
#include<cstdlib>
#include<cstring>
#include<string>
#include<math.h>
#include"Markov.h"

using namespace std;
int main(int argc,char*argv[])
{
    char* file=argv[1];
    sequence seqs;
    load_data(file,seqs);
    cout<<"seq_size= "<<seqs.name.size()<<endl;
    vector<sequence> seqs_bg;
    sequence sq_bg;
    markov(seqs,sq_bg,3);
    char* outfile=file;
    write_seqs(outfile,sq_bg);
    return 0;

}
