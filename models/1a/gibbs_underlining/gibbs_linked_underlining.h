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
const int loop_max=200 ;
int updated_doc=0 ;
#define topic_count 20  // 1 to 30 !
#define total_words_count 10000
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
int linear_scan(double*cdf , int size,double data);
//void  update_topic_usage(document* current_doc_ptr, int  new_topic,int old_topic  ); //update topic ratio !

  double get_current_sum(document* doc,const uint32_t N_d, int sz);
inline void update_topic_usage(document* doc,uint8_t new_topic,const int old_topic);
//evaluate model
void test_model();
int inline  getValue_int_3D_matrix(int *ptr, int x,int y,int z,int x_size,int y_size);
void inline setValue_int_3D_matrix(int *ptr, int x,int y,int z,int value,int x_size,int y_size);
void  print_topic_corpus();
void  print_model_params();
void  print_system_config();
void  read_topic();
void  make_labels();

//Takes sum of two log numbers 
double elsum(double ln_x , double ln_y);
double under_lining(document* doc,int i, int j , int gamma);
void print_topic_corpus(document *doc);
void make_factor_matrix_log(document* doc);
void  reconstruct_doc(document *doc);
void make_int_words(document* doc,vector<string>& word_list );
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

//Set/ read complex arrays 
void inline set_d_2D_m(double *ptr,int row,int col,double value,int row_length);
double inline get_d_2D_m(double* ptr,int row,int col,int row_length);
double  get_f_3D_m(const double *ptr ,int x,int y,int z,int y_size,int z_size);
void  set_f_3D_m(double *ptr,int x,int y,int z,double  value,int y_size,int z_size);
char  get_c_3D_m(const char *ptr ,int x,int y,int z,int y_size,int z_size);
void  set_c_3D_m(char *ptr,int x,int y,int z,char  value,int y_size,int z_size);

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
double get_punctuation_multipler_num(document* doc, int pos,int t_last, int t_next );
double get_punctuation_multipler_den(document* doc, int pos,int t_last, int t_next );

const int document_count=54;
//const string data="/data1/rahul/projects/text_mining/branch/punctuation/chris_data/BF";
const char* data="/home/rkumar2/final_test/data/BF" ;
const char* dict_file="/home/rkumar2/final_test/data/word_list.txt" ;

string write_dir="test_results" ;
void  make_Dictionary();
void run_test(document* useme);
void document_fee_mem(document* d_curr);
void  test_document(string directory_name); // creates document  structure !
string make_uppercase_string(string data);
inline void  scale_pdf( double * pdf ,const int size, double cdf_sum  );
inline void make_cdf_from_pdf(double* pdf ,   double * cdf ,const int size);
inline double rand_0_1();
#endif 
