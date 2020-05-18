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
 #include <climits>

#include "KDTree.hpp"


#define NDEBUG
#include <cassert>

#define Pi 3.14159265
#define pi 3.14156

// d as genomic distance

using namespace std;


class a_node{
public:
    int edge_id;
    point_t frags;
    int isprefix;
};

class an_edge{
public:
    vector<int> bilabels;
    vector<int> loc_freq;
    vector<pair<int,int > > locs;
    int rep_bil;

    a_node suffix;

    vector<a_node> prefixes;
    vector<int> prefreq;

    point_t frags;
};


int partition(int arr[], int l, int r);

// This function returns k'th smallest element in arr[l..r] using
// QuickSort based method.  ASSUMPTION: ALL ELEMENTS IN ARR[] ARE DISTINCT
int kthSmallest(int arr[], int l, int r, int k)
{
    // If k is smaller than number of elements in array
    if (k > 0 && k <= r - l + 1) {
        // Partition the array around last element and get
        // position of pivot element in sorted array
        int pos = partition(arr, l, r);

        // If position is same as k
        if (pos - l == k - 1)
            return arr[pos];
        if (pos - l > k - 1) // If position is more, recur for left subarray
            return kthSmallest(arr, l, pos - 1, k);

        // Else recur for right subarray
        return kthSmallest(arr, pos + 1, r, k - pos + l - 1);
    }

    // If k is more than number of elements in array
    return INT_MAX;
}

void swap(int* a, int* b)
{
    int temp = *a;
    *a = *b;
    *b = temp;
}

// Standard partition process of QuickSort().  It considers the last
// element as pivot and moves all smaller element to left of it
// and greater elements to right
int partition(int arr[], int l, int r)
{
    int x = arr[r], i = l;
    for (int j = l; j <= r - 1; j++) {
        if (arr[j] <= x) {
            swap(&arr[i], &arr[j]);
            i++;
        }
    }
    swap(&arr[i], &arr[r]);
    return i;
}


class bilabel{
    public:
    int kmer_size;
    int locus;
    int orientation;
    int chromosome;
    string rmap_from;
    int rmap_id;
    int i_rmap,q_rmap;
    int sums[3];
    vector<int> fragments;
};


point_t uncover_frags(int kmer_size, int i_rmap, int q_rmap,vector<int>& rmap){

    point_t retvar((kmer_size)*2);
    int i;
    for(i=0;i<kmer_size;i++){

        retvar[i]=rmap[i_rmap++];
    }

    for(int j=q_rmap+kmer_size;j>q_rmap+kmer_size-(kmer_size);j--){
        retvar[i++]=rmap[j];
    }

    return(retvar);
}



bilabel extract_one_bilabel(int index, int kmer_size, int dsize, int start_locus, vector<int>& temp, vector<int>& temp_sum, int orientation, int chromosome, string rmapname,int rmap_id){

    bilabel tembi;
    int locus=start_locus;

    int i=index;
    tembi.kmer_size=kmer_size;
    tembi.rmap_from=rmapname;
    tembi.rmap_id=rmap_id;
    tembi.i_rmap=i;
    tembi.chromosome=chromosome;
    tembi.locus=locus;
    tembi.orientation=orientation;

    int k1_end=i+kmer_size-1;
    int p_rmap=k1_end+1;
    int q_rmap=p_rmap;


    if(i>0)
        tembi.sums[0]=temp_sum[k1_end]-temp_sum[i-1];
    else
        tembi.sums[0]=temp_sum[k1_end];

    while(temp_sum[q_rmap]-temp_sum[k1_end]<dsize){
        q_rmap++;
        if(q_rmap+kmer_size>=temp.size()){
            tembi.kmer_size=-1;
            return(tembi);
        }

    }

    tembi.sums[1]=temp_sum[q_rmap]-temp_sum[k1_end];
    tembi.q_rmap=q_rmap;
    tembi.sums[2]=temp_sum[q_rmap+kmer_size]-temp_sum[q_rmap];

    return(tembi);


}



inline int disth(const point_t &a, const point_t &b, int rad) {
    int distc = 0;
    for (size_t i = 0; i < a.size(); i++) {
        if(abs(a[i]-b[i])>rad){
            distc++;
            break;
        }
    }

    if(distc==0)
        return 1;
    else
        return 0;
}

vector<bilabel> extract_bilabels(int kmer_size, int dsize, int start_locus, vector<int>& temp, vector<int>& temp_sum, int orientation, int chromosome, string rmapname,int rmap_id){

    vector<bilabel> retvar;
    int locus=start_locus;

    for(int i=0;i<temp.size()-(2*kmer_size+3)+1;i++){

        bilabel tembi;

        tembi.rmap_from=rmapname;
        tembi.rmap_id=rmap_id;
        tembi.i_rmap=i;
        tembi.chromosome=chromosome;
        tembi.locus=locus;
        tembi.orientation=orientation;

        int k1_end=i+kmer_size-1;
        int p_rmap=k1_end+1;
        int q_rmap=p_rmap;

        if(i>0)
            tembi.sums[0]=temp_sum[k1_end]-temp_sum[i-1];
        else
            tembi.sums[0]=temp_sum[k1_end];


        while(temp_sum[q_rmap]-temp_sum[k1_end]<dsize){
            q_rmap++;
            if(q_rmap+kmer_size>=temp.size())
                return(retvar);
        }

        tembi.sums[1]=temp_sum[q_rmap]-temp_sum[k1_end];
        tembi.q_rmap=q_rmap;
        tembi.sums[2]=temp_sum[q_rmap+kmer_size]-temp_sum[q_rmap];
        locus+=(temp[i]);

        point_t thisbib=uncover_frags(kmer_size,tembi.i_rmap,tembi.q_rmap, temp);

        retvar.push_back(tembi);

    }

    return(retvar);


}




int min(int a,int b){

    if(a<=b)
        return a;
    else
        return b;
}

int max(int a,int b){

    if(a>=b)
        return a;
    else
        return b;
}


int main(int argc, char *argv[]){

    if(argc<4){
        cout<<"usage: ./buildgraph <rmap_file.val> <kmer_size> <dsize> <t_f>"<<endl;
        return(1);
	}

	int t_sum = 2000;
	int t_frag = atoi(argv[4]);

	int t_nmerge=1500;

	int cutoff_freq=4,supp=3;

    int minmerge=10,maxmerge=1000;

	int dsize=atoi(argv[3]);

    int kmer_size = atoi(argv[2]);

    ofstream bilabelout("bilabels_from_data.txt");
    ofstream prefix_suffix("prefix_suffix.txt");
    ofstream rmapsfile("rmaps.txt");

    ifstream infile(argv[1]);

    map <string,int> nameid;
    int lastname=0;

    string sttt;

    vector<int> starts;
    vector<int> ends;
    vector<int> orientation;
    vector<int> chromosome;

    starts.push_back(-1);
    ends.push_back(-1);
    orientation.push_back(-1);
    chromosome.push_back(-1);

    map<string,int> infomap; // maps Rmap name to Rmap id

    int rmapno=0,lastm=0;

    vector< vector < int > > rmap,rmap_sum;

    string str;

	vector < vector < vector <float> > > bilabel_container;

	vector < vector < int > > bilabel_freq;

	map<string,int> bilabel_pair_map;

	int rmap_no=0;
    int rmap_cnt=0;

	int lastpairmapped=0;

	vector<string> bilabelsums;

	int totnumoffrags=0;

    int selrmaps=0;

    int ik;

    int coverage=0,merged=0;

    int bicnt=0,pre_suf_cnt=0;

    vector< bilabel > all_bilabels,ps_bilabels;

    vector< bilabel > pre_suf;

    while (getline(infile,str)){ // Read Rmap data and extract all bi-labels

        istringstream ss2(str);

        string rmap_name;
        ss2>>rmap_name;

        getline(infile,str);

        istringstream ss1(str);

        string dummy;

        ss1>>dummy>>dummy;

        float input;

        vector<int> temp;
        vector<int> temp_sum;

        int sum=0;

        while(ss1>>input){
            totnumoffrags++;
            int intin=input*1000;
            if(intin<1000)
                continue;

            sum+=intin;
            temp.push_back(intin);
            temp_sum.push_back(sum);
        }

        //----------------------------uncomment for Rmapper2.0

        vector<int> temp_rev;
        vector<int> temp_sum_rev;
        sum=0;
        for(int i=temp.size()-1;i>=0;i--){
            sum+=temp[i];
            temp_rev.push_back(temp[i]);
            temp_sum_rev.push_back(sum);
        }



        vector<bilabel> bilabels_from_this_rmap = extract_bilabels(kmer_size, dsize, starts[0], temp, temp_sum, orientation[0], chromosome[0], rmap_name,rmap_cnt);

        for(int ii=0;ii<bilabels_from_this_rmap.size();ii++){

            bilabel tembi=extract_one_bilabel(bilabels_from_this_rmap[ii].i_rmap, kmer_size-1, dsize, starts[0], temp, temp_sum, orientation[0], chromosome[0], rmap_name,rmap_cnt);

            ps_bilabels.push_back(tembi);

            prefix_suffix<<"bi_"<<pre_suf_cnt++<<" "<<tembi.rmap_from<<" "<<tembi.i_rmap<<" "<<tembi.q_rmap<<" "<<tembi.locus<<" "<<tembi.chromosome<<" "<<tembi.orientation<<" ";
            for(int iii=0;iii<3;iii++){
                prefix_suffix<<" d"<<iii+1<<": "<<tembi.sums[iii]<<" ";
            }

            prefix_suffix<<" km1: ";

            for(int iii=tembi.i_rmap;iii<tembi.i_rmap+kmer_size-1;iii++){
                prefix_suffix<<temp[iii]<<" ";
            }

            prefix_suffix<<" km2: ";

            for(int iii=tembi.q_rmap+1;iii<tembi.q_rmap+kmer_size-1+1;iii++){
                prefix_suffix<<temp[iii]<<" ";
            }

            prefix_suffix<<endl;


            tembi=extract_one_bilabel(bilabels_from_this_rmap[ii].i_rmap+1, kmer_size-1, dsize, starts[0], temp, temp_sum, orientation[0], chromosome[0], rmap_name,rmap_cnt);
            ps_bilabels.push_back(tembi);

            prefix_suffix<<"bi_"<<pre_suf_cnt++<<" "<<tembi.rmap_from<<" "<<tembi.i_rmap<<" "<<tembi.q_rmap<<" "<<tembi.locus<<" "<<tembi.chromosome<<" "<<tembi.orientation<<" ";
            for(int iii=0;iii<3;iii++){
                prefix_suffix<<" d"<<iii+1<<": "<<tembi.sums[iii]<<" ";
            }

            prefix_suffix<<" km1: ";

            for(int iii=tembi.i_rmap;iii<tembi.i_rmap+kmer_size-1;iii++){
                prefix_suffix<<temp[iii]<<" ";
            }

            prefix_suffix<<" km2: ";

            for(int iii=tembi.q_rmap+1;iii<tembi.q_rmap+kmer_size-1+1;iii++){
                prefix_suffix<<temp[iii]<<" ";
            }

            prefix_suffix<<endl;


            bilabelout<<"bi_"<<bicnt++<<" "<<bilabels_from_this_rmap[ii].rmap_from<<" "<<bilabels_from_this_rmap[ii].i_rmap<<" "<<bilabels_from_this_rmap[ii].q_rmap<<" "<<bilabels_from_this_rmap[ii].locus<<" "<<bilabels_from_this_rmap[ii].chromosome<<" "<<bilabels_from_this_rmap[ii].orientation<<" ";
            for(int iii=0;iii<3;iii++){
                bilabelout<<" d"<<iii+1<<": "<<bilabels_from_this_rmap[ii].sums[iii]<<" ";
            }

            tembi=bilabels_from_this_rmap[ii];

            bilabelout<<" km1: ";

            for(int iii=tembi.i_rmap;iii<tembi.i_rmap+kmer_size;iii++){
                bilabelout<<temp[iii]<<" ";
            }

            bilabelout<<" km2: ";

            for(int iii=tembi.q_rmap+1;iii<tembi.q_rmap+kmer_size+1;iii++){
                bilabelout<<temp[iii]<<" ";
            }

            bilabelout<<endl;

        }


        rmapsfile<<rmap_name<<endl;
        for(int ii=0;ii<temp.size();ii++){
            rmapsfile<<temp[ii]<<" ";
        }
        rmapsfile<<endl;
        for(int ii=0;ii<temp_sum.size();ii++){
            rmapsfile<<temp_sum[ii]<<" ";
        }

        rmapsfile<<endl<<endl;


        all_bilabels.insert(all_bilabels.end(),bilabels_from_this_rmap.begin(),bilabels_from_this_rmap.end());

        rmap.push_back(temp);
        rmap_sum.push_back(temp_sum);

        rmap_cnt++;


 //---------------------------------------------------------------------------------------
        rmap_name=rmap_name+"_r";
        temp=temp_rev;
        temp_sum=temp_sum_rev;

        bilabels_from_this_rmap = extract_bilabels(kmer_size, dsize, starts[0], temp, temp_sum, orientation[0], chromosome[0], rmap_name,rmap_cnt);

        for(int ii=0;ii<bilabels_from_this_rmap.size();ii++){

            bilabel tembi=extract_one_bilabel(bilabels_from_this_rmap[ii].i_rmap, kmer_size-1, dsize, starts[0], temp, temp_sum, orientation[0], chromosome[0], rmap_name,rmap_cnt);
            ps_bilabels.push_back(tembi);

            prefix_suffix<<"bi_"<<pre_suf_cnt++<<" "<<tembi.rmap_from<<" "<<tembi.i_rmap<<" "<<tembi.q_rmap<<" "<<tembi.locus<<" "<<tembi.chromosome<<" "<<tembi.orientation<<" ";
            for(int iii=0;iii<3;iii++){
                prefix_suffix<<" d"<<iii+1<<": "<<tembi.sums[iii]<<" ";
            }

            prefix_suffix<<" km1: ";

            for(int iii=tembi.i_rmap;iii<tembi.i_rmap+kmer_size-1;iii++){
                prefix_suffix<<temp[iii]<<" ";
            }

            prefix_suffix<<" km2: ";

            for(int iii=tembi.q_rmap+1;iii<tembi.q_rmap+kmer_size-1+1;iii++){
                prefix_suffix<<temp[iii]<<" ";
            }

            prefix_suffix<<endl;


            tembi=extract_one_bilabel(bilabels_from_this_rmap[ii].i_rmap+1, kmer_size-1, dsize, starts[0], temp, temp_sum, orientation[0], chromosome[0], rmap_name,rmap_cnt);
            ps_bilabels.push_back(tembi);


            prefix_suffix<<"bi_"<<pre_suf_cnt++<<" "<<tembi.rmap_from<<" "<<tembi.i_rmap<<" "<<tembi.q_rmap<<" "<<tembi.locus<<" "<<tembi.chromosome<<" "<<tembi.orientation<<" ";
            for(int iii=0;iii<3;iii++){
                prefix_suffix<<" d"<<iii+1<<": "<<tembi.sums[iii]<<" ";
            }

            prefix_suffix<<" km1: ";

            for(int iii=tembi.i_rmap;iii<tembi.i_rmap+kmer_size-1;iii++){
                prefix_suffix<<temp[iii]<<" ";
            }

            prefix_suffix<<" km2: ";

            for(int iii=tembi.q_rmap+1;iii<tembi.q_rmap+kmer_size-1+1;iii++){
                prefix_suffix<<temp[iii]<<" ";
            }

            prefix_suffix<<endl;


            bilabelout<<"bi_"<<bicnt++<<" "<<bilabels_from_this_rmap[ii].rmap_from<<" "<<bilabels_from_this_rmap[ii].i_rmap<<" "<<bilabels_from_this_rmap[ii].q_rmap<<" "<<bilabels_from_this_rmap[ii].locus<<" "<<bilabels_from_this_rmap[ii].chromosome<<" "<<bilabels_from_this_rmap[ii].orientation<<" ";
            for(int iii=0;iii<3;iii++){
                bilabelout<<" d"<<iii+1<<": "<<bilabels_from_this_rmap[ii].sums[iii]<<" ";
            }

            tembi=bilabels_from_this_rmap[ii];

            bilabelout<<" km1: ";

            for(int iii=tembi.i_rmap;iii<tembi.i_rmap+kmer_size;iii++){
                bilabelout<<temp[iii]<<" ";
            }

            bilabelout<<" km2: ";

            for(int iii=tembi.q_rmap+1;iii<tembi.q_rmap+kmer_size+1;iii++){
                bilabelout<<temp[iii]<<" ";
            }

            bilabelout<<endl;

        }


        rmapsfile<<rmap_name<<endl;
        for(int ii=0;ii<temp.size();ii++){
            rmapsfile<<temp[ii]<<" ";
        }
        rmapsfile<<endl;
        for(int ii=0;ii<temp_sum.size();ii++){
            rmapsfile<<temp_sum[ii]<<" ";
        }

        rmapsfile<<endl<<endl;


        all_bilabels.insert(all_bilabels.end(),bilabels_from_this_rmap.begin(),bilabels_from_this_rmap.end());

        rmap.push_back(temp);
        rmap_sum.push_back(temp_sum);

        rmap_cnt++;


 //---------------------------------------------------------------------------------------



        getline(infile,str);

    }


    vector< vector < int > > binned_bilabel;
    vector<string> bin_id;
    map<string,int> bin_map;
    int lastbin=0;
    int bilabel_counter=0;

    vector< vector < int > > storage_details;

    vector <int> smalv1(200,-1);
    vector <vector <int> > smalv2(200,smalv1);
    vector< vector < vector < int> > > container3d(200,smalv2);



    for(int i=0;i<all_bilabels.size();i++){

        bilabel tembi=all_bilabels[i];

        float d1f=(float)tembi.sums[0]/t_sum;
        float d2f=(float)tembi.sums[1]/t_sum;
        float d3f=(float)tembi.sums[2]/t_sum;

        int a = round(d1f), b = round(d2f), c = round(d3f);

        if(a<200 && b<200 && c<200){

            if(container3d[a][b][c]==-1){
                vector<int> temv;
                temv.push_back(i);
                container3d[a][b][c]=binned_bilabel.size();
                binned_bilabel.push_back(temv);

            }
            else{
                binned_bilabel[container3d[a][b][c]].push_back(i);
            }

            continue;
        }


        int dee1 [] = {a,a+1};
        int dee2 [] = {b,b+1};
        int dee3 [] = {c,c+1};


        for(int h=0;h<1;h++){

            for(int hh=0;hh<1;hh++){

                for(int hhh=0;hhh<1;hhh++){

                    stringstream sss;
                    sss<<dee1[h]<<"_"<<dee2[hh]<<"_"<<dee3[hhh];
                    if(bin_map.find(sss.str())==bin_map.end()){
                        vector<int> temv;
                        temv.push_back(i);

                        bin_map[sss.str()]=binned_bilabel.size();
                        binned_bilabel.push_back(temv);
                        bin_id.push_back(sss.str());
                    }
                    else{
                        binned_bilabel[bin_map.find(sss.str())->second].push_back(i);
                    }

                }

            }

        }



    }

    double avgcnt=0;

//    ofstream recordmerges("recordmerges.txt");

//    ofstream details("storage_details");

//    ofstream kmerdet("merge_kmers_distribution.txt");


    for(int i=0;i<binned_bilabel.size();i++){
//        details<<i<<" "<<binned_bilabel[i].size()<<endl;
        avgcnt+=binned_bilabel[i].size();
    }

    cout<<"Number of bi-labels: "<<all_bilabels.size()<<endl;


    cout<<"Average density: "<<avgcnt/binned_bilabel.size()<<endl;


    //build k-d tree for all containers with atleast 10 bi-labels

    int trees_built=0;

    vector<KDTree> alltrees;
    vector<int> tree_loc(binned_bilabel.size(),-1);

//    ofstream logfile("fullog");

    for(int i=0;i<binned_bilabel.size();i++){

        if(binned_bilabel[i].size()<5){
//            treeptrs[i]=NULL;

        }
        else{

            pointVec points(binned_bilabel[i].size());

            for(int ii=0;ii<binned_bilabel[i].size();ii++){


                bilabel tembi=all_bilabels[binned_bilabel[i][ii]];

                points[ii] = uncover_frags(kmer_size, tembi.i_rmap,tembi.q_rmap,rmap[tembi.rmap_id]);

            }

            tree_loc[i]=trees_built;

            KDTree tree(points);
            alltrees.push_back(tree);
            trees_built++;


        }

//        if(i%10000==0)
//            cout<<i<<" "<<trees_built<<endl;

    }

    cout<<"Trees built: "<<trees_built<<endl;

    ofstream bilmerge("bilabel_merges.txt");

    vector< vector <int> > bilabel_merges(all_bilabels.size());

    vector< int > merge_status(all_bilabels.size(),-1);

    double tot_time=0;

    set< pair<int,int> >::reverse_iterator rit;

    int totcor=0,totmer=0;
    int numgenbi=0;


    for(int i=0;i<all_bilabels.size();i++){

        time_t tstart=clock();

        bilabel tembi=all_bilabels[i];

        point_t pt = uncover_frags(kmer_size, tembi.i_rmap,tembi.q_rmap,rmap[tembi.rmap_id]);


        float d1f=(float)tembi.sums[0]/t_sum;
        float d2f=(float)tembi.sums[1]/t_sum;
        float d3f=(float)tembi.sums[2]/t_sum;

        int a = round(d1f), b = round(d2f), c = round(d3f);


        int dee1 [] = {a-1,a,a+1};
        int dee2 [] = {b-1,b,b+1};
        int dee3 [] = {c-1,c,c+1};

//        cout<< "bi_"<<i<<" ";

        for(int h=0;h<3;h++){

            for(int hh=0;hh<3;hh++){

                for(int hhh=0;hhh<3;hhh++){

                    if(dee1[h]<200 && dee2[hh]<200 && dee3[hhh]<200){

                        if(container3d[dee1[h]][dee2[hh]][dee3[hhh]]==-1)
                            continue;

                        if(tree_loc[container3d[dee1[h]][dee2[hh]][dee3[hhh]]]<0)
                            continue;

                        auto res = alltrees[tree_loc[container3d[dee1[h]][dee2[hh]][dee3[hhh]]]].neighborhood_indices(pt, t_frag);

                        for (int a : res) {

                                if(a>=binned_bilabel[container3d[dee1[h]][dee2[hh]][dee3[hhh]]].size()){
                                    continue;
                                }
                                int merged_bil=binned_bilabel[container3d[dee1[h]][dee2[hh]][dee3[hhh]]][a];
                                if(merged_bil!=i && merge_status[merged_bil]!=i){
                                    bilabel_merges[i].push_back(merged_bil);
                                    merge_status[merged_bil]=i;
                                }
                        }

                        continue;
                    }

                    stringstream sss;
                    sss<<dee1[h]<<"_"<<dee2[hh]<<"_"<<dee3[hhh];

                    map <string,int>::iterator it= bin_map.find(sss.str());

                    if(it==bin_map.end())
                        continue;

                    if(tree_loc[it->second]<0){
                        continue;
                    }
                    else{



                        auto res = alltrees[tree_loc[it->second]].neighborhood_indices(pt, t_frag);

                        for (int a : res) {

                                if(a>=binned_bilabel[it->second].size()){
                                    continue;
                                }

                                int merged_bil=binned_bilabel[it->second][a];
                                if(merged_bil!=i && merge_status[merged_bil]!=i){
                                    bilabel_merges[i].push_back(merged_bil);
                                    merge_status[merged_bil]=i;
                                }
                        }

                    }


                }

            }

        }

//        cout<<endl;

        bilmerge<<"bi_"<<i<<" ";

        for(int ii=0;ii<bilabel_merges[i].size();ii++){
            bilmerge<<bilabel_merges[i][ii]<<" ";

        }
        bilmerge<<endl;



        bilmerge<<endl;

        time_t tend=clock();

        tot_time+=double( tend - tstart) / CLOCKS_PER_SEC;

        if(i>0 && i%10000==0)
            cout<<i<<" "<< " Average_time: "<<tot_time/(i+1)<<" approx_time_remaining: "<<((all_bilabels.size()-i)*tot_time)/(i+1)<<" "<<numgenbi<<endl;


    }



    //---------------------------------------------------------------------Now build bilabelled de bruijn graph

    minmerge=4;maxmerge=2000;

    t_frag=1000;
    t_nmerge=1000;
    cutoff_freq=4,supp=4;


    vector<KDTree>().swap(alltrees);

    infile.close();
    infile.open("bilabel_merges.txt");

    vector<int> selected;

    set< pair<int,int> > selset;

    vector<int> bilstatus(all_bilabels.size(),-2);

    vector<int> bilabel_m_correctness(all_bilabels.size(),-2);

    int badl=0;

    vector < vector<int > > bimerge;

    bimerge.resize(all_bilabels.size());

    int bilcnt=0,good=0;

    while (getline(infile,str)){

        istringstream ss2(str);

        string biname;
        ss2>>biname;
//
//        if(bilmap.find(biname)==bilmap.end())
//            continue;

        int temin;
        vector<int> temvec;

        while(ss2>>temin){

                temvec.push_back(temin);
        }

        bimerge[bilcnt]=temvec;


        if(temvec.size()>=minmerge && temvec.size()<=maxmerge){
            good++;
            bilstatus[bilcnt]++;
            selset.insert(make_pair(temvec.size(),bilcnt));

        }

        bilcnt++;
        getline(infile,str);

    }

    cout<<"Selected bils: "<<selset.size()<<endl;
    cout<<"Percent good: "<<(float)good*100/bimerge.size()<<endl;


    vector< vector <int> > edges_to_bil;

    int edgecnt=0;

//    ofstream edgebil("edge_to_bil.txt");
//    ofstream check_correctedness("check_corr.txt");
//    ofstream edge_range("edge_ranges.txt");


    for(rit=selset.rbegin();rit!=selset.rend();rit++){

        selected.push_back(rit->second);
    }


    vector<an_edge> all_edges;

    vector<a_node> all_nodes;

    vector< vector<int> > graph_fwd,graph_rev;

    vector<int> node_status;

    pointVec node_pts;
    int nodecnt=0;

    pointVec nodepts;

//    ofstream initgraph("initial_graph.txt");

//    ofstream outprepts("investigate_prepts.txt");

    vector<int> prefixmap;
    int lastpre=0;
    pointVec allprefixes;

    int edges_not_considered=0;


    for(int i=0;i<selected.size();i++){

            if(bilstatus[selected[i]]>=0)
                continue;

            vector<int> thisedgebils;
            thisedgebils.push_back(selected[i]);

            bilstatus[selected[i]]=edgecnt;

            for(int ii=0;ii<bimerge[selected[i]].size();ii++){

                if(bimerge[selected[i]][ii]==selected[i])
                    continue;

                bilstatus[bimerge[selected[i]][ii]]=edgecnt;
                thisedgebils.push_back(bimerge[selected[i]][ii]);

            }

//            edgebil<<"edge_"<<edgecnt<<" ";
//            edge_range<<"edge_"<<edgecnt<<" ";

            an_edge newedge;
            newedge.bilabels=thisedgebils;
            newedge.rep_bil=selected[i];
//            newedge.frags=kmers[selected[i]];

            a_node suf;

            vector<a_node> prefixes;

            pointVec sufpts;

            pointVec vecprepts;

            vector<pointVec> storeallpre;

            vector<int> freq;

            int a=0,b=0,c=0,ia=-1,ib=-1,ic=-1;

            point_t precon,sufcon;

            for(int otherbils : thisedgebils ){

                bilabel tembp = ps_bilabels[otherbils*2];
                bilabel tembs = ps_bilabels[otherbils*2+1];

                point_t prefixhere = uncover_frags(kmer_size-1, tembp.i_rmap,tembp.q_rmap,rmap[tembp.rmap_id]);
                point_t suffixhere = uncover_frags(kmer_size-1, tembs.i_rmap,tembs.q_rmap,rmap[tembs.rmap_id]);

                int flag=0;
                for(int ii=0;ii<vecprepts.size();ii++){

                    if(disth(vecprepts[ii],prefixhere,t_nmerge)){
                        storeallpre[ii].push_back(prefixhere);
                        freq[ii]++;
                        flag++;
                        break;
                    }
                }

                if(flag==0){
                    freq.push_back(1);
                    vecprepts.push_back(prefixhere);
                    pointVec ptemp;
                    ptemp.push_back(prefixhere);
                    storeallpre.push_back(ptemp);
                }


//                prepts.push_back(ps_kmers[otherbils*2]);
                sufpts.push_back(suffixhere);

            }

            for(int ii=0;ii<freq.size();ii++){

                if(freq[ii]>=a){
                    a=freq[ii];
                    ia=ii;
                }
                else if(freq[ii]<a && freq[ii]>=b){
                    b=freq[ii];
                    ib=ii;
                }
                else if(freq[ii]<b && freq[ii]>=c){
                    c=freq[ii];
                    ic=ii;
                }
            }



            pointVec consensusprefix;

            for(int jjj=0;jjj<storeallpre.size();jjj++){

                    point_t thiscon;

                    for(int jj=0; jj<storeallpre[jjj][0].size(); jj++){

                        int total=0;
                        for(int jjjj=0;jjjj<storeallpre[jjj].size();jjjj++){

                            total+=storeallpre[jjj][jjjj][jj];

                        }
                        thiscon.push_back(total/storeallpre[jjj].size());
                    }

                    consensusprefix.push_back(thiscon);

            }


//            outprepts<<"edge_"<<edgecnt<<endl;

            pointVec finalprepts;

            for(int ii=0;ii<vecprepts.size();ii++){

                    if(freq[ii]<thisedgebils.size()/3 || freq[ii]<c)
                        continue;

                    a_node candidatenode;
                    candidatenode.frags=consensusprefix[ii];
                    candidatenode.edge_id=edgecnt;
                    candidatenode.isprefix=1;
                    finalprepts.push_back(consensusprefix[ii]);

                    prefixes.push_back(candidatenode);

//                    for(int jj=0; jj<consensusprefix[ii].size();jj++){
//                            outprepts<<consensusprefix[ii][jj]<<" ";

//                    }
//                    outprepts<<" -> "<<freq[ii]<<endl;
            }

//            outprepts<<" ------------------------"<<endl;

            if(finalprepts.size()==0){
                edges_not_considered++;
                continue;
            }



            int arr[sufpts.size()];
            int n=sufpts.size();

            for(int jj=0; jj<sufpts[0].size(); jj++){

                for(int jjj=0;jjj<sufpts.size();jjj++){

                    arr[jjj]=sufpts[jjj][jj];
                }

                sufcon.push_back(kthSmallest(arr, 0, n - 1, n/2));

            }



            suf.frags = sufcon;
            suf.isprefix=0;

//            nodepts.push_back(precon);

            nodepts.insert(nodepts.end(),finalprepts.begin(),finalprepts.end());
            nodepts.push_back(sufcon);

            suf.edge_id=edgecnt;

            newedge.prefixes=prefixes;
            newedge.suffix=suf;

            vector<int> temf,temr,teme;

            all_edges.push_back(newedge);

            for(int ii=0;ii<prefixes.size();ii++){
//                initgraph<<"node_"<<nodecnt<<" edge: "<<edgecnt<<" frags: ";
                prefixmap.push_back(nodecnt);
                allprefixes.push_back(prefixes[ii].frags);


//                for(int ff : prefixes[ii].frags)
//                    initgraph<<ff<<" ";

                temr.push_back(nodecnt++);
//                initgraph<<endl;

            }



//            initgraph<<"node_"<<nodecnt<<" edge: "<<edgecnt<<" frags: ";

//            for(int ff:suf.frags)
//                    initgraph<<ff<<" ";

//            initgraph<<endl;

            temf.push_back(nodecnt++);


            for(int ii=0;ii<prefixes.size();ii++){
                graph_fwd.push_back(temf);
                graph_rev.push_back(teme);
            }


            graph_fwd.push_back(teme);
            graph_rev.push_back(temr);

            all_nodes.insert(all_nodes.end(),prefixes.begin(),prefixes.end());

            all_nodes.push_back(suf);


            edges_to_bil.push_back(thisedgebils);

            edgecnt++;

    }

//    return 1;

    cout<<"Edges_not_con: "<<edges_not_considered<<endl;

//    ofstream graphrevnow("grev_now.txt");

//    for(int index=0;index< graph_rev.size();index++){
//
//        graphrevnow<<"node_"<<index;
//
//        for(int ii:graph_rev[index]){
//            graphrevnow<<" "<<ii;
//        }
//        graphrevnow<<endl;
//    }

    cout<<"Number of edges: "<<edgecnt<<endl;
    cout<<"Number of nodes: "<<nodecnt<<endl;
    cout<<"Number of prefixes: "<<allprefixes.size()<<endl;

    KDTree tree(allprefixes);

//    ofstream explore("explorenodes.txt");


    node_status.resize(all_nodes.size(),-1);

    vector <vector <int> > node_merges;

    vector<int> final_nodes;

    ofstream graph_nodes("graph.txt");

//    ofstream nodemerges("nodemerges.txt");

    int allmerge=0,correctmerge=0;

    vector<int> matchrmap(rmap.size(),-1);

//    ofstream seesupport("seesupport.txt");

    int understood=0,toofew=0,toomany=0;

    for(int i=0;i<all_nodes.size();i++){

        time_t tstart=clock();

        if(all_nodes[i].isprefix==1){
            continue;
        }
//        vector<int> thismerges;

        final_nodes.push_back(i);

        node_status[i]=i;
//        thismerges.push_back(i);

        point_t thisnodept=all_nodes[i].frags;

        auto res = tree.neighborhood_indices(thisnodept, t_frag);

//        nodemerges<<"node_"<<i<<" ";
        int edge1=all_nodes[i].edge_id;

        for (int prenode : res) {

                int a=prefixmap[prenode];
                int suffix_of_a = graph_fwd[a][0];

                if(graph_rev[suffix_of_a].back()==i){
                    node_status[a]=i;
                    vector<int>::iterator it=find(graph_rev[suffix_of_a].begin(),graph_rev[suffix_of_a].end(),a);
                    if(it!=graph_rev[suffix_of_a].end())
                        graph_rev[suffix_of_a].erase(it);

                    continue;
                }

                int edge2=all_nodes[a].edge_id;

                for(int ebilabel : edges_to_bil[edge1]){
                    matchrmap[all_bilabels[ebilabel].rmap_id]=i;
                }
                int support=0;
                for(int ebilabel : edges_to_bil[edge2]){
                    if(matchrmap[all_bilabels[ebilabel].rmap_id]==i)
                        support++;
                }

//                seesupport<<"node_"<<i<<" "<<a<<" "<<support<<endl;

                allmerge++;

                if(support<supp)
                    continue;

                correctmerge++;

//                nodemerges<<a<<" ";
                node_status[a]=i;
//                thismerges.push_back(a);

                graph_fwd[i].push_back(suffix_of_a);
                vector<int>::iterator it=find(graph_rev[suffix_of_a].begin(),graph_rev[suffix_of_a].end(),a);
                if(it!=graph_rev[suffix_of_a].end())
                    graph_rev[suffix_of_a].erase(it);

                graph_rev[suffix_of_a].push_back(i);

        }

//        nodemerges<<endl;


    }



    vector<int> deegdist(1000,0);

    cout<<"Total node merges: "<<allmerge<<" Correct: "<<correctmerge<<" "<<(float)correctmerge*100/allmerge<<endl;


    int prefixesmerged=0;
    int edge_del=0,trim_del=0,prefix_trim=0;

    set <int> delnodes;

    for(int i=0;i<node_status.size();i++){

        if(all_nodes[i].isprefix==1){  // if node is a prefix

            if(node_status[i]>=0){ // if prefix already merged to a suffix then ignore
                prefixesmerged++;
            }
            else{ // if prefix node not merged to suffix then check if the suffix has an outgoing edge. if not then check if suffix has any other incoming edge. if not then delete both prefix and suffix.

                int nextnode=graph_fwd[i][0];

                if(graph_fwd[nextnode].empty() && graph_rev[nextnode].size()==1){
                    delnodes.insert(nextnode);
                    edge_del++;
                    node_status[nextnode]=-1;
//                    cout<<"edge_del node_"<<nextnode<<endl;
                    continue;
                }

                if(graph_rev[nextnode].size()>1){
                    prefix_trim++;

                    vector<int>::iterator it=find(graph_rev[nextnode].begin(),graph_rev[nextnode].end(),i);
                    if(it!=graph_rev[nextnode].end())
                        graph_rev[nextnode].erase(it);

                    continue;
                }

                final_nodes.push_back(i);
            }


        }
        else{ // if node is a suffix

            if(graph_fwd[i].size()==0){ // if suffix was not merged to any prefix i.e no outgoing edge

                int flag_del_node=0;

                for(int ii=0;ii<graph_rev[i].size();ii++){

                    int prenode = graph_rev[i][ii]; // a prefix node of the suffix

                    if(graph_fwd[prenode].size()>1){ // if prefix node has other outgoing edges then trim i

                        delnodes.insert(i);
                        trim_del++;
                        flag_del_node++;
                        break;
                    }


                }

                if(flag_del_node>0){

//                    cout<<"to trim node_"<<i<<endl;

                    for(int ii=0;ii<graph_rev[i].size();ii++){

                        int prenode = graph_rev[i][ii]; // a prefix node of the suffix

//                        cout<<"node_"<<i<<" trimmed prenode: node_"<<prenode<<endl;

                        vector<int>::iterator it=find(graph_fwd[prenode].begin(),graph_fwd[prenode].end(),i);
                        if(it!=graph_fwd[prenode].end())
                            graph_fwd[prenode].erase(it);

                    }

                }

            }

        }


    }

    cout<<"Num prefixes merged: "<<prefixesmerged<<endl;
    cout<<"Num edges deleted: "<<edge_del<<endl;
    cout<<"Num nodes trimmed: "<<trim_del<<endl;
    cout<<"Num prefixes trimmed: "<<prefix_trim<<endl;

    ofstream graphrev("graph_rev.txt");

    ofstream indees("indegrees.txt");

    vector<int> indegrees(1000,0);


    for(int i=0;i<final_nodes.size();i++){

            if(delnodes.find(final_nodes[i])==delnodes.end()){


                graph_nodes<<"node_"<<final_nodes[i]<<" ";

                point_t thisnodept=all_nodes[final_nodes[i]].frags;

                for(int nfrag: thisnodept)
                    graph_nodes<<nfrag<<" ";

                graph_nodes<<" fwd_conns: ";
                for(int fwdn: graph_fwd[final_nodes[i]])
                    graph_nodes<<"node_"<<fwdn<<" ";

                graph_nodes<<endl;

                graphrev<<"node_"<<final_nodes[i]<<" ";

                for(int nfrag: thisnodept)
                    graphrev<<nfrag<<" ";

                graphrev<<" rev_conns: ";
                for(int revn: graph_rev[final_nodes[i]])
                    graphrev<<"node_"<<revn<<" ";

                graphrev<<endl;

                if(graph_rev[final_nodes[i]].size()<1000)
                    indegrees[graph_rev[final_nodes[i]].size()]++;


                if(graph_fwd[final_nodes[i]].size()<1000)
                    deegdist[graph_fwd[final_nodes[i]].size()]++;
            }

    }


    int ij=0;
    for(int i : deegdist){
        ij++;
        if(i==0)
            continue;
        cout<<"degree_"<<ij-1<<": "<<i<< " " <<(float)i*100/(final_nodes.size()-delnodes.size())<<endl;
    }
    ij=0;
    for(int i : indegrees){
        ij++;
        if(i==0)
            continue;
        indees<<"degree_"<<ij-1<<": "<<i<< " " <<(float)i*100/(final_nodes.size()-delnodes.size())<<endl;
    }


}
