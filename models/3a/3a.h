#ifndef __DUMMY_H__
#define __DUMMY_H__
#include<stdio.h>
//#define NDEBUG
#include<stdint.h>
#include<omp.h>
#include<float.h>
#include<assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include<string.h>
#include<dirent.h>
#include<sys/stat.h>
#include<sys/param.h>
#include<sys/types.h>
#include<signal.h>
//GSL library
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>
#include<gsl/gsl_matrix.h>
#include <utility>
#include <string>
#include <iostream>
#include<iterator>
#include<map>
#include<vector>
#include <cctype>
#include<algorithm>
#include<climits>
using namespace std ;

const uint16_t loop_max=500;
uint16_t loop=0;
const double  epsilon= 0.05; //error associated with user variable !
const double  epsilon_m=    0.0000001; //error associated with user variable !
const double epsilon_phi=   0.0000001; //error associated with user variable !
const double epsilon_theta= 0.0000001; //error associated with user variable !
const double lower_limit_theta= 0.0000000001;
const double lower_limit_m =    0.0000000001;
const double lower_limit_phi =    0.000000001;
clock_t  start;
clock_t  end;
gsl_rng *rand_obj ;

double maxima_phi_i_possible=0.00000000000000000000000000000000000001;
double log_Fx_maxima_phi_i_possible=-DBL_MAX;

double * diversity;
uint32_t *freq_single_den; //  freqency of theta_w in denominator

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

#define topic_count 5   // 1 to 30 !
#define document_count 790
#define total_words_count 10000

document *d_ptr ;
//documet
//Number of topics defined in begining 
int updated_doc=0 ;

//GIBBS  samples count
//
int gibbs_counter=0;
const int gibbs_sample_count =25 ;
double *M_gibbs  ; 
double *THETA_gibbs  ;
double *PHI_gibbs ;
void init_gibbs_params();
void write_model_gibbs();
void do_gibbs_expectation();
double alpha_comma_Y_gibbs=0;
double alpha_comma_N_gibbs=0;
double alpha_period_Y_gibbs=0;
double alpha_period_N_gibbs=0;
double alpha_upper_Y_gibbs=0;
double alpha_upper_N_gibbs=0;

double alpha_scale=1;
double alpha_shape=1;


const double  log_ratio=log(1000);


double *M ; 
double *THETA ;
double *PHI;

//File names
const char* M_file="M.bin";
const char* THETA_file="THETA.bin";
const char* PHI_file="PHI.bin";
const char* TOPIC_file="topics_list.bin";

//For keeping last maxima , this will  speed up  search
double *last_maxima_theta_iw;
double *last_maxima_m_ij;
double *last_maxima_phi ;

//Thread info 
const  int m_thread_cnt=8;
const  int theta_thread_cnt=8;

//Priors
const double m_diag=500;
const double m_nondiag=1;
const double m_beta=1;

const double theta_k=1;
const double theta_beta=1;
/*
double phi_sigma_square=10000;
double phi_sigma_square=1;
double phi_mu=0;
*/
//FUNCTIONS
//INIT  functions
void init(); //wrapper   for  all the inits 
void init_model_parameters() ; //Inits model params 
void init_document_params(); //INIT  document structure
void init_variables();
void init_topic_assignment(document* current_doc);
void test_model();
void free_things();

//Topic update for each document
inline void  update_topic_usage(document* current_doc_ptr, const uint8_t  new_topic,const int old_topic  );
void topic_update(document* current_doc_ptr); 
void update_TOPIC();
double get_factor12(const int word_curr,const int word_next,const int t_prime,const int t_last,const int t_next,document* doc,const uint32_t mypos);
//double get_factor12(const uint32_t word_curr,const int word_next,const uint8_t t_prime,const int t_last,const int t_next,document* doc,const uint32_t mypos);
inline double get_factor2_t_i(document* current_doc_ptr, double sum);
double get_current_sum(document* doc,const uint32_t N_d, int sz);
inline void scale_pdf(double* pdf ,const int size,double sum );
inline void make_cdf_from_pdf(double* pdf ,double* cdf , const int size );

//update m_i_j
void  update_M();
void  update_m_ij(const int row,const int col,double * sum_den,uint32_t* freq_den,uint32_t*  freq_num,uint32_t* ref);
void get_max_f_m_ij(const uint8_t i, const uint8_t j,double*log_F_max,double* maxima,double*sum_den,uint32_t*freq_den,uint32_t* freq_num,uint32_t* ref,double * maxima_possible,double * val_maxima);
double evaluate_m_ij_ll(double x0,const uint8_t i,const uint8_t j,double* sum_den,uint32_t* freq_den,uint32_t* freq_num,uint32_t* ref, double *maxima_m_ij_possible,double *log_Fx_maxima_m_ij_possible);
void m_compute_double_term_word_stat(const uint8_t i,double* sum_den ,uint32_t* freq_den,uint32_t* freq_num_i,uint32_t* w_topic_count_map);


//UPDATE THETA matrix
void update_THETA();
void  update_theta_iw(const uint8_t i,const uint32_t w, double* den_single_term_header,double * sum_den,uint32_t* freq_den,const uint32_t freq_single_num,uint32_t *  freq_num,const uint32_t freq_single_den_w,uint32_t * ref);
void get_max_f_theta_iw(const uint8_t i, const uint32_t w,double*log_F_max,double* maxima, double* den_single_term_header,double * sum_den,uint32_t* freq_den,const uint32_t freq_single_num,uint32_t *  freq_num,const uint32_t freq_single_den_w, uint32_t* ref, double * maxima_possible,double *value_possible);
double evaluate_theta_iw_ll(double x0,const uint8_t i,const uint32_t w,const double den_single_term_header,double * sum_den,uint32_t * freq_den,const uint32_t freq_single_num,uint32_t* freq_num,const  uint32_t freq_single_den_w,uint32_t* ref , double *maxima_theta_iw_possible,double *log_Fx_maxima_theta_iw_possible);
double get_den_theta_jw(const uint8_t j , const uint32_t w);
double compute_theta_single_term_denominator(const uint32_t  w );
void theta_compute_double_term_word_stat(const uint32_t w,double* sum_den,uint32_t* freq_den,uint32_t* freq_single_num_w,uint32_t* freq_double_num,uint32_t* ref );

//Common codes
double inline  get_number_between_two_numbers(double low, double high);
inline double rand_0_1();
void ctrl_c_handler(int single_num);
void  write_model();
double timetaken(clock_t start,clock_t  end , const char* mesaage);
double gamma_approx_sampler(double r,double log_F_max,double log_f_r,double maxima);
void   guess_lebel(document* current_doc_ptr);
double m_diagonistics(const int i,const int j,double* sum_den,uint32_t* freq_den,uint32_t* freq_num , uint32_t* ref ,double *maxima_possible,double*value_possible,const char*msg , double maxima);


double theta_diagonistics(const int i, const int w,const double den_single_term_header,double * sum_den,uint32_t* freq_den,const uint32_t freq_single_num,uint32_t *  freq_num,const  uint32_t freq_single_den_w,uint32_t* ref, double *possible_maxima,double * val_at_max,const char* msg,double maxima);
int linear_scan(double*cdf , const int size,double data);
double inline get_d_2D_m(double *ptr, const int row,const int col,const int row_length);
void inline set_d_2D_m(double *ptr, const int row,const int col,double value,const int row_length);
uint32_t  get_i_3D_m(const uint32_t *ptr ,int x,int y,int z,int y_size,int z_size);
void  set_i_3D_m(uint32_t *ptr,int x,int y,int z,int value,int x_size,int y_size);
uint32_t inline get_i_2D_m(uint32_t *ptr, int row,int col,int row_length);
void inline set_i_2D_m(uint32_t *ptr, int row,int col,int value,int row_length);
void load_model_params();
void load_topics();

//update PHI 
double  evaluate_phi_i(const double x0,const uint8_t i );
void update_PHI();
void   get_max_f_phi_i(const uint8_t i, double* log_Fx,double* maxima ) ; // calculate maximum of function of log !

void print_topic_corpus();
void print_model_params();
void print_system_config();


void reset_topic_usage();

uint32_t  get_key_3D_m(int x,int y,int z,int x_size,int y_size);
void print_time();

//Makes word  of document 
//void make_int_words(document* doc);

void make_int_words(document* doc,vector<string>& word_list );

//These things cam because of new data interpretation
//Now theyi are in strings so,I will use these data structures for 
//string interpretation.
void make_word_id(string word);
std::map<string , int > word_id_map;
std::map<int , string > id_word_map;
string  get_word_from_id(int id);

int  get_id_from_word(string word);
//while we read document we will silently populate  this map for furtehr uses !
void  create_word_id_id_word_map(string word);

//helper  function for day to word 
bool isWordStartsWithCaps(string word);
bool isCommaBeforeMe(string word) ;
bool isPeriodBeforeMe(string word) ;


//New variables  for  \alpha_c_Y and so on !
double alpha_comma_Y=0;
double alpha_comma_N=0;
double alpha_period_Y=0;
double alpha_period_N=0;
double alpha_upper_Y=0;
double alpha_upper_N=0;
// shape  and scale of \alpha


//make string lower case
string make_lowercase_string(string word);
bool doesWordStartsUpper(string data);

const string highlight_tag_start="<:t:>";
const string highlight_tag_end="<:/t:>";
const string comma=",";
const string period=".";

double get_punctuation_multipler_num(document* doc, int pos,int t_last,int t_next );
double get_punctuation_multipler_den(document* doc, int pos,int t_last, int t_next );

#define  total_configs_punctuation 8
double multiplier_values[total_configs_punctuation];
void fill_multiplier_values();

bool  didReturnDBL_MAX=false;

//alpha_comma_Y
void sample_alpha_comma_Y();
void  make_alpha_comma_Y_statistics(int &num_count,double* sum_den,map<int,vector<int> >& ref);
double evaluate_alpha_comma_Y_log_posterior(double x0,const int &num_count, double* sum_den,const map<int , vector<int> > & ref,double *maxima_alpha_comma_Y_possible,double *log_Fx_maxima_alpha_comma_Y_possible);
void get_max_f_alpha_comma_Y(double* log_F_max,double *maxima ,int &num_count,double*sum_den,map<int,vector<int> >& ref,double * maxima_possible,double * val_maxima); 
double alpha_comma_Y_diagonistics(int &num_count,double*sum_den,map<int,vector<int> >& ref,double * possible_maxima,double * val_at_max ,const char* msg,double maxima);

//alpha_comma_N
void sample_alpha_comma_N();
void  make_alpha_comma_N_statistics(int &num_count,double* sum_den,map<int,vector<int> >& ref);
double evaluate_alpha_comma_N_log_posterior(double x0,const int &num_count, double* sum_den,const map<int , vector<int> > & ref,double *maxima_alpha_comma_N_possible,double *log_Fx_maxima_alpha_comma_N_possible);
void get_max_f_alpha_comma_N(double* log_F_max,double *maxima ,int &num_count,double*sum_den,map<int,vector<int> >& ref,double * maxima_possible,double * val_maxima); 
double alpha_comma_N_diagonistics(int &num_count,double*sum_den,map<int,vector<int> >& ref,double * possible_maxima,double * val_at_max ,const char* msg,double maxima);

//alpha_upper_Y
void sample_alpha_upper_Y();
void  make_alpha_upper_Y_statistics(int &num_count,double* sum_den,map<int,vector<int> >& ref);
double evaluate_alpha_upper_Y_log_posterior(double x0,const int &num_count, double* sum_den,const map<int , vector<int> > & ref,double *maxima_alpha_upper_Y_possible,double *log_Fx_maxima_alpha_upper_Y_possible);
void get_max_f_alpha_upper_Y(double* log_F_max,double *maxima ,int &num_count,double*sum_den,map<int,vector<int> >& ref,double * maxima_possible,double * val_maxima); 
double alpha_upper_Y_diagonistics(int &num_count,double*sum_den,map<int,vector<int> >& ref,double * possible_maxima,double * val_at_max ,const char* msg,double maxima);

//alpha_upper_N
void sample_alpha_upper_N();
void  make_alpha_upper_N_statistics(int &num_count,double* sum_den,map<int,vector<int> >& ref);
double evaluate_alpha_upper_N_log_posterior(double x0,const int &num_count, double* sum_den,const map<int , vector<int> > & ref,double *maxima_alpha_upper_N_possible,double *log_Fx_maxima_alpha_upper_N_possible);
void get_max_f_alpha_upper_N(double* log_F_max,double *maxima ,int &num_count,double*sum_den,map<int,vector<int> >& ref,double * maxima_possible,double * val_maxima); 
double alpha_upper_N_diagonistics(int &num_count,double*sum_den,map<int,vector<int> >& ref,double * possible_maxima,double * val_at_max ,const char* msg,double maxima);

//alpha_upper_Y
void sample_alpha_period_Y();
void  make_alpha_period_Y_statistics(int &num_count,double* sum_den,map<int,vector<int> >& ref);
double evaluate_alpha_period_Y_log_posterior(double x0,const int &num_count, double* sum_den,const map<int , vector<int> > & ref,double *maxima_alpha_period_Y_possible,double *log_Fx_maxima_alpha_period_Y_possible);
void get_max_f_alpha_period_Y(double* log_F_max,double *maxima ,int &num_count,double*sum_den,map<int,vector<int> >& ref,double * maxima_possible,double * val_maxima); 
double alpha_period_Y_diagonistics(int &num_count,double*sum_den,map<int,vector<int> >& ref,double * possible_maxima,double * val_at_max ,const char* msg,double maxima);

//alpha_period_N
void sample_alpha_period_N();
void  make_alpha_period_N_statistics(int &num_count,double* sum_den,map<int,vector<int> >& ref);
double evaluate_alpha_period_N_log_posterior(double x0,const int &num_count, double* sum_den,const map<int , vector<int> > & ref,double *maxima_alpha_period_N_possible,double *log_Fx_maxima_alpha_period_N_possible);
void get_max_f_alpha_period_N(double* log_F_max,double *maxima ,int &num_count,double*sum_den,map<int,vector<int> >& ref,double * maxima_possible,double * val_maxima); 
double alpha_period_N_diagonistics(int &num_count,double*sum_den,map<int,vector<int> >& ref,double * possible_maxima,double * val_at_max ,const char* msg,double maxima);

///home/rkumar2/final_test/data
///data1/rahul/projects/text_mining/branch/punctuation/chris_data
const char* data="/data1/rahul/projects/text_mining/branch/punctuation/chris_data/AB" ;
const char* dict_file="/data1/rahul/projects/text_mining/branch/punctuation/chris_data/word_list.txt" ;
#endif 
