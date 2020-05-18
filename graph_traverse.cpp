#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <map>
#include <fstream>
#include <sstream>
#include <vector>
#include <list>
#include <algorithm>
#include <set>
#include <ctime>
#include <vector>
#include <stdlib.h>
#include <typeinfo>
#include <string>
#include <bits/stdc++.h>
#include <climits>



#define NDEBUG
#include <cassert>

#define Pi 3.14159265
#define pi 3.14156

using namespace std;

vector< vector <int> > graphf,graphr;
int max_lvl_reached;
vector< vector<int> > nodefrag;

vector<int> nodevisted;

vector<int> thispath,thiscontig,longestcontig,longestpath;

vector<int> nodestatus;

int kmersize;
int nothit;
int maxlvl;

//void DFS(int source, int snode){
//
//    connodes.push_back(source);
//
//    if(nodestatus[source]==snode || graphf[source].size()==0 || graphf[source].size()>1){
//
//        for(int i=5;i<nodefrag[source].size();i++){
//            thiscontig.push_back(nodefrag[source][i]);
//        }
//
//        nodestatus[source]=snode;
//
//
//        if(nodestatus[source]==snode){
//            cout<<"Cycle detected."<<endl;
//
//        }
//
//
//        return;
//
//    }
//    else{
//        thiscontig.push_back(nodefrag[source][5]);
//        nodestatus[source]=snode;
//        DFS(graphf[source][0],snode);
//    }
//
//}

void DFS(int source, int thisnode, int thislvl, int toplvl){


    if(nodestatus[thisnode]==source){
        if(longestpath.size()<thispath.size()){
            longestcontig=thiscontig;
            longestpath=thispath;
        }
        return;
    }

    thispath.push_back(thisnode);
    thiscontig.push_back(nodefrag[thisnode][kmersize-1]);
    nodestatus[thisnode]=source;

    int iterationcontigsize=thiscontig.size();

    if(graphf[thisnode].size()==1){
        DFS(source,graphf[thisnode][0],thislvl+1,toplvl);
    }
    else{
            for(int ii=kmersize;ii<nodefrag[thisnode].size();ii++){
                thiscontig.push_back(nodefrag[thisnode][ii]);
            }

            longestcontig=thiscontig;
            longestpath=thispath;
    }

//    for(int nextnode: graphf[thisnode]){
//
//
//
//        thispath.resize(thislvl);
//        thiscontig.resize(iterationcontigsize);
//
//    }

//    if(longestpath.size()<thispath.size()){
//
//
//    }

}

void DFS_trim(int source, int thisnode, int thislvl, int toplvl){


//    if(nodestatus[thisnode]==source){
//        if(longestpath.size()<thispath.size()){
//            longestcontig=thiscontig;
//            longestpath=thispath;
//        }
//        return;
//    }

    if(nodestatus[thisnode]==source){
        nothit++;
    }

    if(nodestatus[thisnode]==source || thislvl>=toplvl){
        if(thislvl>maxlvl)
            maxlvl=thislvl;
        return;
    }

    nodestatus[thisnode]=source;

    for(int nextnode: graphf[thisnode]){

        DFS_trim(source,nextnode,thislvl+1,toplvl);

    }

    if(thislvl>maxlvl)
            maxlvl=thislvl;

}



int main(int argc, char *argv[]){



    ifstream infile1("graph.txt");
    ifstream infile2("graph_rev.txt");

    vector<vector<int> > fwd,rev;
    map<string,int> nodemap;
    int lastnode=0;


    vector< vector <string> > stringfwd;

    int dist=10000;

    string str;

    int totalin=0,totalout=0;


    vector<int> fwd_deeg(100,0);

    ofstream notenodes("nodenotes.txt");

    vector<string> storeallstr;

    while(getline(infile1,str)){

        stringstream sss(str);
        storeallstr.push_back(str);

        string instr;

        sss>>instr;

        notenodes<<"nn_"<<lastnode<<" "<<str<<endl;

        nodemap[instr]=lastnode++;

        vector<int> temvec;

        while(sss>>instr){

            if(instr.compare("fwd_conns:")==0)
                break;
            temvec.push_back(stoi(instr));
        }

        nodefrag.push_back(temvec);

        vector<string> vecstr;

        while(sss>>instr){
            vecstr.push_back(instr);
            totalout++;
        }

        if(vecstr.size()<100){
            fwd_deeg[vecstr.size()]++;
        }

        stringfwd.push_back(vecstr);

    }

    vector<int> sources,terminals,multiedges;

    int nodecnt=0;

    int problem=0;

    while(getline(infile2,str)){

        stringstream sss(str);

        string instr;

        sss>>instr;

        while(sss>>instr){

            if(instr.compare("rev_conns:")==0)
                break;
        }


        vector<string> vecstr;
        vector<int> temvec;

        while(sss>>instr){


            if(nodemap.find(instr)==nodemap.end()){
                cout<<"inproblem"<<nodecnt<<" "<<instr<<endl;
                problem++;
            }
            else{
                temvec.push_back(nodemap.find(instr)->second);
            }

            totalin++;
        }

        if(temvec.size()==0)
            sources.push_back(nodecnt);

        graphr.push_back(temvec);
        nodecnt++;

    }

    cout<<"Graph read. Total in: "<<totalin<<" Total out: "<<totalout<<endl;

//    cout<<"Inproblem: "<<problem;

    problem=0;

    int countextra=0,oneoutmorein=0;
    vector<int> singleton;

    ofstream graphcut("graphcut.txt");


    for(int i=0;i<stringfwd.size();i++){

        vector<int> temvec;
        for(int ii=0;ii<stringfwd[i].size();ii++){
            if(nodemap.find(stringfwd[i][ii])==nodemap.end()){
                cout<<"outprobelm: "<<i<<" "<<stringfwd[i][ii]<<endl;
                problem++;
            }
            else{
                temvec.push_back(nodemap.find(stringfwd[i][ii])->second);
            }

        }
        if(temvec.size()>15){
            graphcut<<storeallstr[i]<<endl;
            temvec.clear();
        }


        graphf.push_back(temvec);

        if(temvec.size()==0)
            terminals.push_back(i);

        if(temvec.size()>1)
            multiedges.push_back(i);


        if(temvec.size()==1 && graphr[i].size()==1)
            singleton.push_back(i);

        if(temvec.size()==1 && graphr[i].size()>1)
            oneoutmorein++;
    }
    cout<<endl;

//    cout<<"extra: "<<countextra<<endl;
    cout<<"Sources: "<<sources.size()<<" Terminals: "<<terminals.size()<< " Multiedges: "<<multiedges.size()<<" Singletons: "<< singleton.size()<<" Oneoutmorein: "<<oneoutmorein<<endl;
//    cout<<"Outprobelm: "<<problem<<endl;

//    for(int i=0;i<fwd_deeg.size();i++){
//
//        cout<<"deegree_"<<i<<" "<<fwd_deeg[i]<<endl;
//    }

    if(problem>0)
        return 1;


    int maxlevel=dist;

    ofstream maxreached("maxreached.txt");

    int resolved=0;

    int goodlength=0;

    vector< vector<int> > all_contigs;
    vector< vector<int> > contignodes;

    vector< vector <int> > source_contigs;

    vector< pair<int,int> > node_to_con(graphf.size(),make_pair(-1,-1));

    nodestatus.resize(graphf.size(),-1);

    kmersize = (nodefrag[0].size()-1)/2;

    ofstream finalcons("assembly.txt");
    ofstream outnodes("connodes.txt");


    int longestcon=0,totlen=0;

    vector< pair<int,int> > nodereachedat(graphf.size(),make_pair(-1,-1));

    int resol=0;

    for(int i=0;i<multiedges.size();i++){

        nodestatus[multiedges[i]]=multiedges[i];

        int maxlvlreached=0,maxedge=-1;


        vector<int> updatededges;
        maxreached<<multiedges[i]<<endl;



        for(int outedge: graphf[multiedges[i]]){

            nothit=0;
            maxlvl=0;
            DFS_trim(multiedges[i],outedge,1,50);
            maxreached<<"\t"<<outedge<<" "<<maxlvl;

            if(maxlvl==50){
                if(!updatededges.empty()){
                    maxreached<<" Keep";
                }
                updatededges.push_back(outedge);
                maxlvlreached=maxlvl;

            }
            else{
                if(maxlvlreached<maxlvl){
                    maxlvlreached=maxlvl;
                    maxedge=outedge;
                }

            }

            maxreached<<endl;

        }

        if(updatededges.empty()){
            updatededges.push_back(maxedge);
        }

        if(updatededges.size()==1)
            resol++;
        else{
            for(int keepedge: updatededges)
                sources.push_back(keepedge);
        }

        graphf[multiedges[i]]=updatededges;
    }

    cout<<"Multiedges unresolved: "<<multiedges.size()-resol<<endl;

//    return 1;

    vector<int> skipped (sources.size(),0);

    for(int i=0;i<sources.size();i++){

        thiscontig.clear();
        thispath.clear();

        longestcontig.clear();
        longestpath.clear();

//        tobeincluded=0;

        for(int ii=0;ii<kmersize-1;ii++){

            thiscontig.push_back(nodefrag[sources[i]][ii]);
        }
        thispath.push_back(sources[i]);

        nodestatus[sources[i]]=sources[i];

        for(int nextnode : graphf[sources[i]]){

            DFS(sources[i],nextnode,2,10000);

        }



        for(int ii=0;ii<longestpath.size();ii++){

            int nodethis=longestpath[ii];

            if(nodereachedat[nodethis].first==-1){
                nodereachedat[nodethis]=make_pair(i,ii);
            }
            else{

                if(nodereachedat[nodethis].second<ii){
                    skipped[nodereachedat[nodethis].first]++;
                    nodereachedat[nodethis]=make_pair(i,ii);
                }
                else{
                    skipped[i]++;
                }

                break;


            }

        }



        all_contigs.push_back(longestcontig);
        contignodes.push_back(longestpath);


    }


    vector<int> nodesnotvisited;

    vector<int> nodesincluded(graphf.size(),0);


    for(int ii=0;ii<nodestatus.size();ii++){

        if(nodestatus[ii]==-1)
            nodesnotvisited.push_back(ii);
    }

    cout<<"Nodes not visited: "<<nodesnotvisited.size()<<endl;

    int countedcontigs=0;
    double longestgenome=0;
    float totgenome=0;
    int filsize=30;

    cout<<"numcons: "<<all_contigs.size()<<" "<<endl;


//    for(int i=0;i<all_contigs.size();i++){
//        cout<<i<<" "<<all_contigs[i].size()<<endl;
//    }

//    return 1;



    for(int i=0;i<all_contigs.size();i++){

        if(skipped[i]>0 || all_contigs[i].size()<filsize)
            continue;

        countedcontigs++;

        finalcons<<"contig_"<<i<<endl<<"\t"<<all_contigs[i].size()<<" enzyme ";
        outnodes<<"contig_"<<i<<" ";

        for(int ii:contignodes[i]){
            outnodes<<" "<<ii;
            nodesincluded[ii]++;
        }


        outnodes<<endl;


        if(all_contigs[i].size()>longestcon)
            longestcon=all_contigs[i].size();

        totlen+=all_contigs[i].size();

        int thisgenlen=0;

        for(int ii:all_contigs[i]){

            finalcons<<(float)ii/1000<<" ";
            thisgenlen+=ii;
        }

        if(thisgenlen>longestgenome)
            longestgenome=thisgenlen;


        totgenome+=((float)thisgenlen/1000000);

        finalcons<<endl<<endl;

    }

    int notinc=0;
    for(int ii:nodesincluded){
        if(ii==0)
            notinc++;
    }


    cout<<"Contigs below size "<< filsize<<" filtered."<<endl;
    cout<<"Number of contigs: "<<countedcontigs<<endl;
    cout<<"Longest contig size: "<<longestcon<<endl;
    cout<<"Mean contig size: "<<totlen/countedcontigs<<endl;
    cout<<"Longest contig length in Mbp: "<<longestgenome/1000000<<endl;
    cout<<"Mean contig length in Mbp: "<<totgenome/(countedcontigs)<<endl;
    cout<<"Nodes not included: "<<notinc<< " Num nodes:"<<graphf.size()<<endl;


}
