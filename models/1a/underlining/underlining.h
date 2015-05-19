#ifndef __TEST_H__
#define __TEST_H__
#include<stdio.h>
#include<assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include<string.h>
#include<dirent.h>
#include<sys/stat.h>
#include<sys/param.h>
#include<sys/types.h>
#include<stdint.h>
#include <iostream>
#include <fstream>
#include<string>
#include<complex>
#include<limits>
#include <vector>
#include <map>
#include <algorithm>
#include <float.h>
using namespace std;

const double tol=6;
const int loop_max=15 ;
int updated_doc=0 ;
double LOGZERO= log(0);
double UNINITIALIZED=log(-1);

vector <string> file_names;

const char* document_file="document_name.txt";
const char* labels_file="labels.bin";

 double *F;
 int *back_ptr ;
 int GAMMA=200;
double inline getValue_double_2D_matrix(double *ptr, int row,int col,int row_length);
void inline setValue_double_2D_matrix(double *ptr, int row,int col,double value,int row_length);
int inline getValue_int_3D_matrix(int *ptr, int x,int y,int z,int x_size,int y_size);
void inline setValue_int_3D_matrix(int *ptr, int x,int y,int z,int value,int x_size,int y_size);
int inline getValue_int_2D_matrix(int *ptr, int row,int col,int row_length);
void inline setValue_int_2D_matrix(int *ptr, int row,int col,int value,int row_length);

double inline get_d_2D_m(double* ptr,int row,int col,int row_length);

const char *TOPIC_FILE="topics.bin";

#define topic_count 5  // 1 to 30 !
#define total_words_count 10000

//ALWAYS   SETS  THIS COREECTLY FOR  lOGGING PURPOSE !
#define m_diag 1
double *M ;
double *THETA ;
double *PHI ;

int loop=0;
int* m_ij_count_word ; 


typedef struct {
  bool* word_caps; // keeps tag if word at this pos was  
  bool* comma_present_before; 
  bool* period_present_before; 
  vector<string> word_list ; 
  uint32_t *int_word_list ; // Its an array to words
  uint32_t  *topic_usage;
  int real_word_count ; //Number of words in document
  int true_label ;
  int sampled_label;
  uint8_t *topic_list ; //Topic associated with each word 
  char* document_name ;
  bool* priority_flag ;
  bool  hasPriorityWord;
  int mygamma;
}document ;

document *d_ptr ;
const char* M_file="../M.bin";
const char* THETA_file="../THETA.bin";
const char* PHI_file="../PHI.bin";
void load_model_params();
//init specific document
void init() ;
void init_document_params();
void init_topic_assignment(document* current_doc);
double  get_current_sum(document*  current_doc_ptr);
int linear_scan(double*cdf , int size,double data);
void  update_topic_usage(document* current_doc_ptr, int  new_topic,int old_topic  ); //update topic ratio !
//evaluate model
void test_model();
int inline  getValue_int_3D_matrix(int *ptr, int x,int y,int z,int x_size,int y_size);
void inline setValue_int_3D_matrix(int *ptr, int x,int y,int z,int value,int x_size,int y_size);
void        update_m_i_j_count_word(int from_topic_i,int  from_topic_j,int word);
void        print_topic_corpus();
void        print_model_params();
void        print_system_config();
double      get_d_3D_m(const double *ptr ,uint32_t x,uint32_t y,uint32_t z,uint32_t y_size,uint32_t z_size);
void        set_d_3D_m(double *ptr,uint32_t x,uint32_t y,uint32_t z,double  value,uint32_t y_size,uint32_t z_size);
double     *log_initial_cost ;
void        read_topic();
void        make_labels();
void        read_document_names();
void        initialize_memory_for_DP(document*  doc);
void initialize_base_case_for_DP(document * doc);

//Takes sum of two log numbers 
double elsum(double ln_x , double ln_y);
double under_lining(document* doc,int i, int j , int gamma);
void print_topic_corpus(document *doc);

int* path;
double *factor_matrix_log;
void make_factor_matrix_log(document* doc);
double *single_word_factor ;
void make_single_word_factor_matrix_log(double * single_word_factor);

void print_path(int i , int j , int gamma);
void print_filled_table(document* doc);
bool isBaseCase(int word_index , int j , int gamma);

void  reconstruct_doc(document *doc);
void make_int_words(document* doc,vector<string>& word_list );
//void make_int_words(document* doc);

//These things cam because of new data interpretation
//Now theyi are in strings so,I will use these data structures for 
//string interpretation.
void make_word_id(string word);
std::map<string , int > word_id_map;
std::map<int , string > id_word_map;
string  get_word_from_id(const int id);

int  get_id_from_word(string word);
//while we read document we will silently populate  this map for furtehr uses !
void  create_word_id_id_word_map(string word);

//helper  function for day to word 
bool isWordStartsWithCaps(string word);
bool isCommaBeforeMe(string word) ;
bool isPeriodBeforeMe(string word) ;

double alpha_comma_Y=1;
double alpha_comma_N=1;
double alpha_period_Y=1;
double alpha_period_N=1;
double alpha_upper_Y=1;
double alpha_upper_N=1;
// shape  and scale of \alpha

//make string lower case
string make_lowercase_string(string word);
bool doesWordStartsUpper(string data);

const string highlight_tag_start="<:t:>";
const string highlight_tag_end="<:/t:>";
const string comma=",";
const string period=".";

double get_punctuation_multipler(document* doc, int pos );
#define  total_configs_punctuation 8
double multiplier_values[total_configs_punctuation];
void fill_multiplier_values();

double get_fjwTheta(int j, int w );
double get_punctuation_multipler_num(document* doc, int pos,int t_last, int t_next );
double get_punctuation_multipler_den(document* doc, int pos,int t_last, int t_next );

string src1="/data1/rahul/projects/text_mining/branch/punctuation/chris_data/CF";
//string src1="/data1/rahul/projects/text_mining/branch/punctuation/dec10_run/1a/CF";
const char* dict_file="/data1/rahul/projects/text_mining/branch/punctuation/chris_data/word_list.txt" ;

string write_dir="test_results" ;

void  make_Dictionary();
void run_test(document* useme);
void document_fee_mem(document* d_curr);
void  test_document(string directory_name); // creates document  structure !
#endif 
