#include "2a.h"
#include "M_IJ.h"
#include "THETA_IW.h"
#include  "PHI_I.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnPrint.h"
#include <algorithm>
#include <vector>
#include<sstream>
using namespace ROOT::Minuit2;
using namespace std;

void  make_Dictionary(){
  char buffer[1000];
  FILE* fp;
  fp= fopen(dict_file,"r");
  assert(fp!=NULL);
  memset(buffer,'\0',1000);
  char* pch ;
  //get words and remove \n from there end !
  while(fgets(buffer,1000,fp)!=NULL){
    pch=strtok(buffer," = \n"); // these are the tokens from file!
    int i=1;
    int  count=-1;
    string word ;
    while (pch != NULL)
    {
      if(i==1){
        word=pch;
      }else{
        count=atoi(pch);
      }
      pch = strtok (NULL, "  = \n");
      i++;
      // cout << "word: " <<  doc << " " << " count:" <<  count << endl ;

      create_word_id_id_word_map(word);
    }
    //    dictionary.push_back(doc);
    //doc_offsetmap[doc]=count;
    memset(buffer,0,1000);
  }//while ends here !      
  fclose(fp);  
}

void init(){
  make_Dictionary();
  init_variables();
  init_model_parameters();
  init_document_params();
  for(uint32_t w=0;w<total_words_count;w++){
    for(int q2=0;q2<document_count;q2++){
      document* d=d_ptr+q2;
      if(d->int_word_list[0]==w){
        freq_single_den[w]++;
      }
    }//frequncy how many times this word appeared at first pos and theta_iw was used !
  }
  //GIBBS initialisation !
  //init_gibbs_params();
}

void init_variables(){
  freq_single_den= (uint32_t*)calloc(total_words_count ,sizeof(uint32_t));
  assert(freq_single_den!=NULL);

  last_maxima_phi=(double*)malloc(sizeof(double)* topic_count);
  assert(last_maxima_phi!=NULL);
  memset(last_maxima_phi,0,sizeof(double)*topic_count);
  last_maxima_theta_iw=(double*) calloc(topic_count*total_words_count,sizeof(double));
  assert(last_maxima_theta_iw!=NULL);
  last_maxima_m_ij=(double*) calloc(topic_count*topic_count,sizeof(double));
  assert(last_maxima_m_ij!=NULL);

  rand_obj = gsl_rng_alloc(gsl_rng_mt19937);
  M=(double*)malloc(sizeof(double)*topic_count*topic_count);
  assert(M!=NULL);
  memset(M,0,sizeof(double)*topic_count*topic_count);

  d_ptr=(document *) calloc(document_count,sizeof(document));
  assert(d_ptr!=NULL);

  THETA=(double*)malloc(sizeof(double)*topic_count*total_words_count);
  assert(THETA != NULL);
  memset(THETA,0,sizeof(double)*topic_count*total_words_count);

  PHI=(double*)malloc(sizeof(double)*topic_count);
  assert(PHI!=NULL);
  memset(PHI,0,sizeof(double)*topic_count);

  //Laplace spread parameter !
  diversity=(double*)calloc(topic_count,sizeof(double));
  assert(diversity!=NULL);

}

void init_gibbs_params(){
  M_gibbs=(double*)malloc(sizeof(double)*topic_count*topic_count);
  assert(M_gibbs!=NULL);
  memset(M_gibbs,0,sizeof(double)*topic_count*topic_count);

  THETA_gibbs=(double*)malloc(sizeof(double)*topic_count*total_words_count);
  assert(THETA_gibbs != NULL);
  memset(THETA_gibbs,0,sizeof(double)*topic_count*total_words_count);

  PHI_gibbs=(double*)malloc(sizeof(double)*topic_count);
  assert(PHI_gibbs!=NULL);
  memset(PHI_gibbs,0,sizeof(double)*topic_count);
}

void do_gibbs_expectation(){
  bool  hasTimeCome = loop >loop_max-gibbs_sample_count ? true :false ;
  if(hasTimeCome){

    gibbs_counter++ ;
    for(int i= 0; i<topic_count;i++){
      PHI_gibbs[i]+=PHI[i];
      for(int j= 0; j<topic_count;j++){
        double addme=get_d_2D_m(M,i,j,topic_count);
        double before_addme=get_d_2D_m(M_gibbs,i,j,topic_count);
        set_d_2D_m(M_gibbs,i,j,addme+before_addme,topic_count);
      }

    }
    for(int i= 0; i<topic_count;i++){
      for(int j= 0; j<total_words_count;j++){
        double addme=get_d_2D_m(THETA,i,j,total_words_count);
        double before_addme=get_d_2D_m(THETA_gibbs,i,j,total_words_count);
        set_d_2D_m(THETA_gibbs,i,j,addme+before_addme,total_words_count);
      }

    }
    //Now average  for expectation
    if(loop==loop_max-1){
      for(int i= 0; i<topic_count;i++){
        PHI_gibbs[i]=PHI_gibbs[i]/gibbs_counter;
        for(int j= 0; j<topic_count;j++){
          double before=get_d_2D_m(M_gibbs,i,j,topic_count);
          set_d_2D_m(M_gibbs,i,j,before/gibbs_counter,topic_count);
        }
      }
      for(int i= 0; i<topic_count;i++){
        for(int j= 0; j<total_words_count;j++){
          double before=get_d_2D_m(THETA_gibbs,i,j,total_words_count);
          set_d_2D_m(THETA_gibbs,i,j,before/gibbs_counter,total_words_count);
        }
      }
    }

  }

}

void init_model_parameters(){
  for(int i=0;i<topic_count;i++){
    for(int j=0;j<topic_count;j++)
    {
      double m_k= (i==j)? m_diag : m_nondiag ;
      double value=gsl_ran_gamma(rand_obj,m_k,m_beta);
      set_d_2D_m(M,i,j,value,topic_count);
      // M[i][j]=value;
      set_d_2D_m(last_maxima_m_ij,i,j,value,topic_count);			
    }
  }

  //Init topic word affinity matrix matrix
  for(uint8_t i=0;i<topic_count;i++){
    for(uint32_t j=0;j<total_words_count;j++){
      double value=gsl_ran_gamma(rand_obj,theta_k,theta_beta);
      //THETA[i][j]=value;
      set_d_2D_m(THETA,i,j,value,total_words_count);
      set_d_2D_m(last_maxima_theta_iw,i,j,value,total_words_count);			
    }
  }
  //Initialise prior on Laplace
  //Intialise by IG(1,1) , Chris suggested so 
  //a=1(shape),b=1(scale) in  terms of Gaussian vocab !
  for(int  i=0;i<topic_count;i++){
    diversity[i]=1.0/gsl_ran_gamma(rand_obj,1,1.0/1);
  } 

  for(int i=0;i<topic_count;i++){
    double temp=gsl_ran_laplace(rand_obj,diversity[i]);
    PHI[i]=temp;
  }

} //init_model ends here !

void init_document_params(){
  DIR *dir;
  struct dirent *ent ;
  dir=opendir(data);
  while((ent=readdir(dir))!=NULL){
    if(strchr(ent->d_name,'.')==NULL){
      FILE*  fp;
      char filename[1000];
      char *pch;
      document * current_doc_ptr=d_ptr+ (updated_doc++) ;
      //Extract name
      strcpy(filename,data);
      strcat(filename,"/");
      strcat(filename,ent->d_name);
      strcat(filename,"\0");
      current_doc_ptr->document_name=(char*)malloc(sizeof(char)*1000);
      assert(current_doc_ptr->document_name!=NULL);
      strcpy(current_doc_ptr->document_name,ent->d_name); 
      if(current_doc_ptr->document_name[0]=='p'){
        current_doc_ptr->true_label= 1;
      }else if (current_doc_ptr->document_name[0]=='n'){
        current_doc_ptr->true_label=-1;
      }else{
        printf("some  problem in file reading !\n");
        printf("%s\n",current_doc_ptr->document_name);
        exit(0);
      }
      current_doc_ptr->hasPriorityWord=false;
      //for topic  usage
      current_doc_ptr->topic_usage=(uint32_t*)malloc(sizeof(uint32_t)*(topic_count));
      assert(current_doc_ptr->topic_usage!=NULL);
      memset(current_doc_ptr->topic_usage,0,sizeof(uint32_t)*(topic_count));
      fp= fopen(filename,"r");
      assert(fp!=NULL);
      memset(filename,'\0',1000);
      vector<string> w_list ;
      current_doc_ptr->word_list=w_list;
      while(fgets(filename,1000,fp)!=NULL){
        pch=strtok(filename,"\n"); // these are the tokens from file!
        string  word=pch; //these are the words 
        w_list.push_back(word);
      }//while ends here !      
      fclose(fp); 
      assert(current_doc_ptr->document_name!=NULL);
      vector<string>& ref= w_list ;
      assert(current_doc_ptr->document_name!=NULL);
      make_int_words(current_doc_ptr,ref);
      w_list.clear();
      init_topic_assignment(current_doc_ptr);
    }//if ends here nn   ONE DOCUMENT
  } //while ends here 
  closedir(dir);
}

void make_int_words(document* doc,vector<string>& word_list ){
  assert(doc->document_name!=NULL);
  bool isUnderliningON=false;
  doc->mygamma=0;
  vector<int> word_index;
  vector<bool>word_caps_vec;
  vector<bool>priority_flag_vec;
  vector <string> delimiter_vec;
  string space=" ";
  string delimiter=space;
  string word="";
  doc->real_word_count=0;
  for(vector<string>::const_iterator it=word_list.begin(); it!=word_list.end(); ++it ){
    word=*it;
    if(word.compare(highlight_tag_start)==0){
      isUnderliningON=true;
      doc->hasPriorityWord=true;
      continue ;
    }
    if(word.compare(highlight_tag_end)==0){
      isUnderliningON=false;
      continue ;
    }
    if(!(word.compare(comma)==0 || word.compare(period)==0 || word.compare(highlight_tag_start)==0 || word.compare(highlight_tag_end)==0)){
      int id =  get_id_from_word(word); 
      if(id==-1){

        //  cout << "Word: " << word << endl ;
        //  cout << "Parsing error  in document: " << doc->document_name   << endl ;
        continue;
      }
      doc->real_word_count++ ;
      word_index.push_back(id);
      //priority flag
      if(isUnderliningON){
        priority_flag_vec.push_back(true);
        doc->mygamma++ ;
      }else{
        priority_flag_vec.push_back(false);
      }
      if(doesWordStartsUpper(word)){
        word_caps_vec.push_back(true);
      }else{
        word_caps_vec.push_back(false);
      }
      delimiter_vec.push_back(delimiter);
      delimiter=space;
    }else{
      //set the delimiter
      if(!(word.compare(highlight_tag_start)==0 || word.compare(highlight_tag_end)==0)){
        if(word.compare(comma)==0){
          delimiter=comma;
        }else if(word.compare(period)==0){
          delimiter=period;
        }else{
          delimiter=space ;
        }      
      }else{
        delimiter=space;
      }
    }
    word.clear();
  }

  int sz=doc->real_word_count;
  assert(word_index.size()==sz);
  assert(delimiter_vec.size()==sz);
  assert(priority_flag_vec.size()==sz);
  assert(word_caps_vec.size()==sz);
  doc->int_word_list=(uint32_t*)calloc(sz,sizeof(uint32_t)) ;
  assert(doc->int_word_list!=NULL);

  doc->priority_flag=(bool*)calloc(sz,sizeof(bool));
  assert(doc->priority_flag!=NULL);
  for(int i=0;i<sz;i++){
    //integer word list
    doc->int_word_list[i]=word_index[i];
    //priority_flag 
    doc->priority_flag[i]=priority_flag_vec[i];

    if(delimiter_vec[i].compare(comma)==0){

    }else if(delimiter_vec[i].compare(period)==0){

    }else{
      // means delimiter was white space 
    } 

  }


  //clear every thing 
  //
  word_index.clear();
  word_caps_vec.clear();
  priority_flag_vec.clear();
  delimiter_vec.clear();
  space.clear();
  delimiter.clear();
  word.clear();

}

void inline set_d_2D_m(double *ptr,int row,int col,double value,int row_length){
  ptr[row_length*row+col]=value;
}
double inline get_d_2D_m(double* ptr,int row,int col,int row_length){
  assert(row>= 0 && col >=0 && row_length >=0 );
  return *(ptr+(row_length)*row+col);
}
uint32_t  get_i_3D_m(const uint32_t *ptr ,int x,int y,int z,int y_size,int z_size){
  assert(x>= 0 && y >=0 && z >=0 && y_size >=0 &&  z_size >=0 );
  return *(ptr+ x*y_size*z_size+y*z_size+z ) ;
}
void  set_i_3D_m(uint32_t *ptr,int x,int y,int z,int  value,int y_size,int z_size){
  assert(x>= 0 && y >=0 && z >=0 && z_size >=0 &&  y_size >=0 );
  assert(value >=0 );
  *(ptr+ x*y_size*z_size+y*z_size+z )=value ;
}

void inline set_i_2D_m(uint32_t *ptr,int row,int col,int value,int row_length){
  *(ptr+(row_length)*row+col)=value;
}
uint32_t inline get_i_2D_m(uint32_t* ptr,int row,int col,int row_length){
  return *(ptr+(row_length)*row+col);
}

void init_topic_assignment(document * doc){
  time_t seconds ;
  time(&seconds);
  uint8_t  topic=0 ;
  uint32_t N_d=doc->real_word_count;
  doc->topic_list= (uint8_t*)malloc(N_d*sizeof(uint8_t));
  memset(doc->topic_list,0,N_d*sizeof(uint8_t));
  assert(doc->topic_list!=NULL);
  bool hasPriorityWord=doc->hasPriorityWord;
  for(uint32_t i=0;i<N_d;i++){
    if(doc->priority_flag[i]==true){
      doc->topic_list[i]=topic_count-1;
      update_topic_usage(doc,topic_count-1,-1);
      continue;
    }
    const uint32_t word = doc->int_word_list[i];
    double den=0;
    for(int ii=0;ii<topic_count;ii++){
      den+=get_d_2D_m(THETA,ii,word,total_words_count);
    }
    double  *pdf = (double*)calloc(topic_count,sizeof(double));
    double  *cdf = (double*)calloc(topic_count,sizeof(double));
    double  cdf_sum=0;
    for(int ii=0;ii<topic_count;ii++){
      pdf[ii]=get_d_2D_m(THETA,ii,word,total_words_count)/den ;
      //probablity of fixed topic !
      if(ii==topic_count-1 && hasPriorityWord){
        pdf[ii]=0;
      }
      cdf_sum+=pdf[ii];
    }
    scale_pdf(pdf,topic_count,cdf_sum);
    make_cdf_from_pdf(pdf,cdf,topic_count);
    double cdf_rand=rand_0_1();
    topic=linear_scan(cdf,topic_count,cdf_rand) ;
    free(pdf);
    free(cdf);
    //        assert(topic>=0 && topic<topic_count-1);
    doc->topic_list[i] = topic;
    const int old_topic=-1;
    update_topic_usage(doc,topic,old_topic);
  }
}

inline void update_topic_usage(document* doc,uint8_t new_topic,const int old_topic){
  uint32_t   temp=0 ;
  if(old_topic==-1){
    //means first word of document !and we are still initing hence 
    temp= doc->topic_usage[new_topic];
    doc->topic_usage[new_topic]=temp+1;
  }else{
    //old topic is being replaced by new topic 
    temp= doc->topic_usage[new_topic];
    doc->topic_usage[new_topic]=temp+1;
    temp= doc->topic_usage[old_topic];
    assert(temp>0);
    doc->topic_usage[old_topic]=temp-1;
  }
}

double get_factor12(const uint32_t word_curr,const int word_next,const uint8_t t_prime,const int t_last,const int t_next,document* doc,const uint32_t mypos){
  bool isFirstWord=t_last==-1 ? true :false ;
  bool isLastWord=word_next==-1 ? true :false ;
  double result =0;
  if(!isLastWord && !isFirstWord){
    double num=get_d_2D_m(M,t_last,t_prime,topic_count) *
      get_d_2D_m(M,t_prime,t_next,topic_count)*
      get_d_2D_m(THETA,t_prime,word_curr,total_words_count);
    double den=0;
    //for(uint8_t i=0;i<topic_count-1;i++){
    for(uint8_t i=0;i<topic_count;i++){
      den+=get_d_2D_m(M,t_prime,i,topic_count)
        * get_d_2D_m(THETA,i,word_next,total_words_count);
    }
    assert(num>0);
    result=num/den;
  }else if (isFirstWord){
    double num=get_d_2D_m(THETA,t_prime,word_curr,total_words_count)* get_d_2D_m(M,t_prime,t_next,topic_count);
    double den=0;
    //for(uint8_t i=0;i<topic_count-1;i++){
    for(uint8_t i=0;i<topic_count-1;i++){
      den+=get_d_2D_m(M,t_prime,i,topic_count) *get_d_2D_m(THETA,i,word_next,total_words_count);
    }
    assert(num>0);
    result=num/den;
  }else if(isLastWord){
    double num=get_d_2D_m(THETA,t_prime,word_curr,total_words_count)* get_d_2D_m(M,t_last,t_prime,topic_count);

    assert(num>0);
    result=num;
  }else{
    printf("Error state !\n");
    exit(0);
  }

  return result;
  }

  inline double  get_factor2_t_i(document* doc,double sum){
    double res=0;
    //  double x=0;
    const int label= doc->true_label;
    res=1/(1+exp(-1*label*sum));
    assert(res>0);
    assert(isinf(res)==0);
    assert(isnan(res)==0);
    return res;
  }

  double get_current_sum(document* doc,const uint32_t N_d, int sz){
    double sum=0;
    for(uint32_t i=0;i<sz;i++){
      sum+=PHI[i] * (double)(doc->topic_usage[i])/N_d ;
    }

    //    printf("sum=%f \n",sum);
    assert(isinf(sum)==0);
    assert(isnan(sum)==0);
    return sum ;
  }

  void topic_update(document*  doc){
    const uint32_t N_d= doc->real_word_count;
    double sum=get_current_sum(doc,N_d,topic_count);//This is same
    double *t_prime_pdf=(double*)malloc(sizeof(double)*(topic_count));
    assert(t_prime_pdf!=NULL);
    memset(t_prime_pdf,0,sizeof(double)*(topic_count));   
    double *t_prime_cdf=(double*)malloc(sizeof(double)*(topic_count));
    assert(t_prime_cdf!=NULL);
    memset(t_prime_cdf,0,sizeof(double)*(topic_count)); 
    bool  isDocumentUnderLined= doc->hasPriorityWord ;
    for(uint32_t i=0;i<N_d;i++)//update all word-topic association !
    {
      const uint32_t word_curr=doc->int_word_list[i];
      const uint8_t before_update= doc->topic_list[i];
      const int  t_last= ((i!=0 )? doc->topic_list[i-1]:-1);
      const int t_next=((i!=(N_d-1)) ? doc->topic_list[i+1]:-1);//set next to one if we are on last word 
      const int word_next=(i==(N_d-1))?-1 :doc->int_word_list[i+1];
      double delta_minus = PHI[before_update]/(double)N_d;
      bool isUnderLinedWord=doc->priority_flag[i];
      //-1 because we dont want to assign last  topic to anybody!       //Not  true, Entirely !
      int new_topic=-1;
      if(isUnderLinedWord){
        new_topic=topic_count-1;
        assert(isDocumentUnderLined==true);
        //continue;
      }else{
        double sum_cdf=0.0;
        for(uint8_t j=0;j<topic_count;j++){ 
          const uint8_t t_prime=j;
          double fact12_t_prime= get_factor12(word_curr,word_next,t_prime,t_last,t_next,doc,i);
          assert(fact12_t_prime>0);
          double delta_plus=PHI[t_prime]/(double)N_d;//j is t_i_prime !
          double sum_prime= sum - delta_minus +delta_plus ; 
          double fact2_t_prime=get_factor2_t_i(doc,sum_prime);
          //        sum_cdf+=(fact12_t_prime*fact2_t_prime);
          t_prime_pdf[j]=fact12_t_prime*fact2_t_prime;
          if(j==topic_count-1 && isDocumentUnderLined){
            t_prime_pdf[j]=0;
          }
          sum_cdf+=t_prime_pdf[j];
        }
        scale_pdf(t_prime_pdf,topic_count,sum_cdf);
        make_cdf_from_pdf(t_prime_pdf,t_prime_cdf,topic_count);
        double cdf_rand=rand_0_1();
        new_topic=linear_scan(t_prime_cdf,topic_count,cdf_rand) ;
        if(isDocumentUnderLined){
          assert(new_topic!=topic_count-1);
        }
      }
      update_topic_usage(doc,new_topic,before_update);
      doc->topic_list[i]=new_topic ;
      memset(t_prime_pdf,0,sizeof(double)*(topic_count));
      memset(t_prime_cdf,0,sizeof(double)*(topic_count));
      //I must update sum ALSO !
      double delta_plus_new=PHI[new_topic]/(double)N_d;//j is t_i_prime !
      sum+=delta_plus_new- delta_minus;
    }//done over the document !
    free(t_prime_pdf);
    free(t_prime_cdf);
    double b_p= 1/(1+exp(-1*sum));
    assert(b_p >0);
    assert(isnan(b_p)==0);
    assert(isinf(b_p)==0);
    double  rand_1= rand_0_1();
    doc->sampled_label= ( b_p>rand_1) ? 1 : -1;
  }//function ends here !

  inline void  scale_pdf( double * pdf ,const int size, double cdf_sum  ){
    double scaling_factor= 1/cdf_sum ;
    int i=0;
    for(i=0;i<size;i++){
      pdf[i]=pdf[i]*scaling_factor;
    }
  }

  inline void make_cdf_from_pdf(double* pdf ,   double * cdf ,const int size){
    int i=0;
    for(i=0;i<size;i++){
      if(i!=0){
        cdf[i]=pdf[i]+cdf[i-1];
      }else{
        cdf[i]=pdf[i];
      }
    }
  }

  int linear_scan(double* cdf ,const int size, double data){
    int index=-1;
    int i=0;
    for(i=0 ; i<size;i++){
      if (cdf[i]>data){
        index= i;
        break ;
      }
    }
    if (index==-1){
      index=size-1;
    }
    return index;
  }

  inline double rand_0_1(){
    return (double)rand()/(double)RAND_MAX;
  }

  void  update_m_ij(const int row,const int col,double * sum_den,uint32_t* freq_den,uint32_t*  freq_num){
    double maxima_possible=0.0000000000000000000000000001;
    double log_Fx_possible=-DBL_MAX;
    double  maxima=0;
    double  log_F_max=0;
    double log_f_r=0;
    // double log_f_l=0;
    //double m_old=M[row][col];
    double m_old=get_d_2D_m(M,row,col,topic_count);
    get_max_f_m_ij(row,col,&log_F_max,&maxima,sum_den,freq_den,freq_num,&maxima_possible,&log_Fx_possible); 
    double r=maxima ;
    bool done =false ;
    double step_size_m=0.01*maxima;
    while(!done ){
      r+=(step_size_m) ;
      log_f_r=evaluate_m_ij_ll(r,row,col,sum_den,freq_den,freq_num,&maxima_possible,&log_Fx_possible);
      step_size_m=2*step_size_m;
      if( log_F_max-log_f_r > log_ratio+epsilon ){
        done=true;
      }
    }
    double sampled_number= gamma_approx_sampler(r,log_F_max,log_f_r,maxima);
    //UPDATE SHARED VARIABLES 
    for(uint32_t s1=0;s1<total_words_count;s1++){
      //sum_den[s1]=sum_den[s1]+THETA[col][s1]*(sampled_number-m_old);
      sum_den[s1]=sum_den[s1]+get_d_2D_m(THETA,col,s1,total_words_count)*(sampled_number-m_old);
    }
    //reset  the possible
    maxima_possible= 0.0000000000000000000000000001;
    log_Fx_possible=-DBL_MAX;
    //update the values
    set_d_2D_m(M,row,col,sampled_number,topic_count);
  }

  void get_max_f_m_ij(const uint8_t i, const uint8_t j,double*log_F_max,double* maxima,double*sum_den,uint32_t*freq_den,uint32_t* freq_num,double * maxima_possible,double * val_maxima){
    M_IJ fcn(i,j,sum_den,freq_den,freq_num,maxima_possible,val_maxima);
    //Call to minimiser !
    MnUserParameters upar;
    MnUserParameterState state(upar);
    double m_ij_start=get_d_2D_m(last_maxima_m_ij,i,j,topic_count);
    upar.Add("m_ij", m_ij_start, epsilon_m);
    upar.SetLowerLimit("m_ij",lower_limit_m);
    MnMigrad migrad(fcn, upar);
    FunctionMinimum min = migrad();
    if(min.IsValid()){
      *maxima=min.UserState().Value("m_ij");
      *log_F_max=-1*min.Fval(); // log(f(x)) is returned !
    }else{
      printf("\nMinimization failed in M[%d][%d]\n",i,j);
      *maxima=*maxima_possible;
      *log_F_max=*val_maxima; // log(f(x)) is returned !
      // std::cout << min  << std::endl << std::endl;
      printf("possible maxima=%.20f, val=%.20f \n" ,*maxima,*log_F_max );
      *maxima=m_diagonistics(i,j,sum_den,freq_den,freq_num,maxima_possible,val_maxima,"failed in get max_f_m",*maxima);
      *log_F_max=evaluate_m_ij_ll(*maxima,i,j,sum_den,freq_den,freq_num,maxima_possible,val_maxima);
      printf("Improved maxima= %.20f , and fn val=%.20f\n",*maxima,*log_F_max);
    }
    set_d_2D_m(last_maxima_m_ij,i,j,*maxima,topic_count);			
  }

  void get_max_f_theta_iw(const uint8_t i, const uint32_t w,double*log_F_max,double* maxima, double* den_single_term_header,double * sum_den,uint32_t* freq_den,const uint32_t freq_single_num,uint32_t *  freq_num,const uint32_t freq_single_den_w,double * maxima_possible,double *value_possible){
    THETA_IW fcn(i,w,den_single_term_header,sum_den,freq_den,freq_single_num,freq_num,freq_single_den_w,maxima_possible,value_possible);
    //Call to minimiser !
    MnUserParameters upar;
    MnUserParameterState state(upar);
    double theta_iw_start=get_d_2D_m(last_maxima_theta_iw,i,w,total_words_count);
    upar.Add("theta_iw", theta_iw_start,epsilon_theta);
    upar.SetLowerLimit("theta_iw",lower_limit_theta);
    MnMigrad migrad(fcn, upar);
    FunctionMinimum min = migrad();
    if(min.IsValid()){
      *maxima=min.UserState().Value("theta_iw");
      *log_F_max=-1*min.Fval(); // - log(f(x)) is returned !a
    }else{
      printf("\nMinimization failed in Theta[%d][%d]\n",i,w);
      // std::cout << min  << std::endl << std::endl;
      *maxima=*maxima_possible;
      *log_F_max=*value_possible; // log(f(x)) is returned !
      printf("possible maxima=%.25f, val=%.20f \n" ,*maxima,*log_F_max );
      *maxima=theta_diagonistics(i,w,*den_single_term_header,sum_den,freq_den,freq_single_num,freq_num,freq_single_den_w,maxima_possible,value_possible,"couldn't find maxima",*maxima);
      *log_F_max=evaluate_theta_iw_ll(*maxima, i, w, *den_single_term_header, sum_den, freq_den, freq_single_num, freq_num,freq_single_den_w,maxima_possible,value_possible);
      printf("Improved maxima=%.25f ,fn_val=%.20f\n",*maxima,*log_F_max);
      // exit(0);
    }
    set_d_2D_m(last_maxima_theta_iw,i,w,*maxima,total_words_count);	
  }

  void get_max_f_phi_i(const uint8_t i,double* log_F_max,double* maxima){

    PHI_I fcn(i);
    MnUserParameters upar;
    MnUserParameterState state(upar);
    upar.Add("phi_i", last_maxima_phi[i],epsilon_phi);
    MnMigrad migrad(fcn, upar);
    FunctionMinimum min = migrad();

    if(min.IsValid()){
      *maxima=min.UserState().Value("phi_i");
      *log_F_max=-1*min.Fval(); //  log(f(x)) is returned !
    }else{
      printf("Minimization failed of PHI[%d] \n",i);
      *maxima=maxima_phi_i_possible;
      *log_F_max=log_Fx_maxima_phi_i_possible; // log(f(x)) is returned !
      printf("possible maxima=%.20f\n" ,*maxima );
      //  std::cout << min  << std::endl << std::endl;
    }
    last_maxima_phi[i]=*maxima;
  }

  double evaluate_m_ij_ll(double x0,const uint8_t i,const uint8_t j,double* sum_den,uint32_t* freq_den,uint32_t* freq_num,double *maxima_m_ij_possible,double *log_Fx_maxima_m_ij_possible){
    double value=0;
    const double m_k=( i==j ? m_diag:m_nondiag);
    if(x0<= 0){
      printf("Called  value was =%.20f ",x0);
      exit(0);
    }
    assert(isnan(x0)==0);
    assert(x0>0);
    value=(m_k-1)*log(x0)+(-1*x0/m_beta); // prior done 
    double x0_old=get_d_2D_m(M,i,j,topic_count);//save it for later updating it 
    value+= (freq_num[0]*log(x0));
    for(uint32_t s2=0;s2<total_words_count;s2++){
      //value+= ((freq_num[s2]*log(x0))-(freq_den[s2]*log(sum_den[s2]+THETA[j][s2]*(x0-x0_old))));
      value+= -(freq_den[s2]*log(sum_den[s2]+get_d_2D_m(THETA,j,s2,total_words_count)*(x0-x0_old)));
    }
    //update den as per new x0
    if(value >*log_Fx_maxima_m_ij_possible){
      *log_Fx_maxima_m_ij_possible=value;
      *maxima_m_ij_possible=x0;
    }
    assert(isnan(value)==0);
    assert(isinf(value)==0);
    return value;
  }

  void update_M(){

    omp_set_num_threads(m_thread_cnt);
#pragma omp parallel for

    for(uint8_t i=0; i<topic_count;i++){
      double*  sum_den=(double*)calloc(total_words_count,sizeof(double));
      assert(sum_den!=NULL);
      uint32_t *  freq_den=(uint32_t*)calloc(total_words_count,sizeof(uint32_t));
      assert(freq_den!=NULL);
      //uint32_t *freq_num_i=(uint32_t*)calloc(topic_count,sizeof(uint32_t)* total_words_count);
      //uint32_t *freq_num_i=(uint32_t*)calloc(topic_count*total_words_count,sizeof(uint32_t));
      uint32_t *freq_num_i=(uint32_t*)calloc(topic_count,sizeof(uint32_t));
      assert(freq_num_i!=NULL);
      m_compute_double_term_word_stat(i,sum_den,freq_den,freq_num_i);
      for(uint8_t j=0;j<topic_count;j++){
        uint32_t temp= freq_num_i[j];
        //uint32_t* freq_num=freq_num_i+j*topic_count;
        uint32_t* freq_num= &temp;
        update_m_ij(i,j,sum_den,freq_den,freq_num);
      }
      free(freq_num_i);
      free(sum_den);
      free(freq_den);
    }
  }



  double evaluate_theta_iw_ll(double x0,const uint8_t i,const uint32_t w,const double den_single_term_header,double * sum_den,uint32_t * freq_den,const uint32_t freq_single_num,uint32_t* freq_num,const  uint32_t freq_single_den_w,double *maxima_theta_iw_possible,double *log_Fx_maxima_theta_iw_possible){
    double value=0;
    value=(theta_k-1)*log(x0)+(-1*x0/theta_beta); // prior done 
    if(x0<= 0){
      printf("Called  value was =%.20f ",x0);
      exit(0);
    }
    assert(x0>0);
    assert(isnan(x0)==0);
    double  x0_old= get_d_2D_m(THETA,i,w,total_words_count);
    //have all denominator for f(t,w,theta,M);
    double den_single_term=den_single_term_header+ x0-x0_old;
    value+=(freq_single_num*log(x0) - freq_single_den_w*log(den_single_term));
    double updated_sum[topic_count];
    memcpy(updated_sum,sum_den,sizeof(double)*topic_count);
    for(uint32_t s1=0;s1<topic_count;s1++){
      //updated_sum[s1]+=(M[s1][i]*(x0-x0_old));
      updated_sum[s1]+=(get_d_2D_m(M,s1,i,topic_count)*(x0-x0_old));
    }
    //NUMERATOR
    value+= (freq_num[0]*(log(x0)));
    for(uint32_t s2=0;s2<topic_count;s2++){
      value-= (freq_den[s2]*log(updated_sum[s2]));// - because in denominator
    }
    //update den as per new x0
    if(*log_Fx_maxima_theta_iw_possible<value){
      //set  maxima & value
      *log_Fx_maxima_theta_iw_possible=value;
      *maxima_theta_iw_possible=x0;
    }
    assert(isnan(value)==0) ;
    assert(isinf(value)==0) ;
    return value;
  }

  void  update_theta_iw(const uint8_t i,const uint32_t w, double* den_single_term_header,double * sum_den,uint32_t* freq_den,const uint32_t freq_single_num,uint32_t *  freq_num,const uint32_t freq_single_den_w){
    double log_Fx_possible=-DBL_MAX;
    double maxima_possible=0.0000000000000000000000001;
    //double theta_old= THETA[i][w] ;
    double theta_old= get_d_2D_m(THETA,i,w,total_words_count) ;
    double log_f_r=0;
    double  maxima=0;
    double  log_F_max=0;
    get_max_f_theta_iw(i,w,&log_F_max,&maxima,den_single_term_header,sum_den,freq_den,freq_single_num,freq_num,freq_single_den_w,&maxima_possible,&log_Fx_possible);
    double r=maxima ;
    bool    done =false ;
    //double step_size_theta=0.0000001;
    double step_size_theta=0.01*maxima;
    //double step_size_theta=0.000000000000000001;
    //right limit
    while(!done ){
      //r+=(step_size_theta*maxima) ;
      r+=step_size_theta ;
      log_f_r=evaluate_theta_iw_ll(r,i,w,*den_single_term_header,sum_den,freq_den,freq_single_num,freq_num,freq_single_den_w,&maxima_possible,&log_Fx_possible);
      step_size_theta=1.3*step_size_theta;
      if( log_F_max-log_f_r > log_ratio+epsilon){
        done=true;
      }
    }
    double sampled_number= gamma_approx_sampler( r,log_F_max,log_f_r,maxima);

    //UPDATE SHARED VALUES
    (*den_single_term_header)+=(sampled_number-theta_old );
    //update sum vector becuase we have sampled a new number !
    for(uint32_t s1=0;s1<topic_count;s1++){
      //sum_den[s1]+=M[s1][i]*(sampled_number-theta_old);
      sum_den[s1]+=get_d_2D_m(M,s1,i,topic_count)*(sampled_number-theta_old);
    }

    set_d_2D_m(THETA,i,w,sampled_number,total_words_count);
    //  printf("THETA[%d][%d]=%f ",i,w,sampled_number);
  }

  void update_THETA(){

    omp_set_num_threads(theta_thread_cnt);
#pragma omp parallel for
    for(int w=0;w<total_words_count;w++){
      uint32_t* freq_double_num=(uint32_t*)malloc(sizeof(uint32_t)* topic_count*total_words_count);
      assert(freq_double_num!=NULL);
      memset(freq_double_num,0,sizeof(uint32_t)*topic_count*total_words_count);
      uint32_t * freq_single_num_w=(uint32_t*)malloc(sizeof(uint32_t)*topic_count);
      assert(freq_single_num_w!=NULL);
      memset(freq_single_num_w,0,sizeof(uint32_t)*topic_count);
      double *sum_den=(double*)malloc(topic_count*sizeof(double));  
      assert(sum_den!=NULL);
      memset(sum_den,0,sizeof(double)*topic_count);
      uint32_t *freq_den= (uint32_t*)malloc(topic_count*sizeof(uint32_t));
      assert(freq_den!=NULL);
      memset(freq_den,0,sizeof(uint32_t)*topic_count);
      const uint32_t freq_single_den_w=freq_single_den[w];
      theta_compute_double_term_word_stat(w,sum_den,freq_den,freq_single_num_w,freq_double_num);
      double  den_single_term_header=compute_theta_single_term_denominator(w);
      for(uint8_t i=0; i<topic_count;i++){
        //  uint32_t* freq_double_num_iw=freq_double_num+topic_count*i;
        uint32_t   freq_num_wi=freq_double_num[i];
        uint32_t* freq_double_num_iw=&freq_num_wi;

        update_theta_iw(i,w,&den_single_term_header,sum_den,freq_den,freq_single_num_w[i],freq_double_num_iw,freq_single_den_w);
      } 
      free(freq_double_num);
      free(freq_single_num_w);
      free(sum_den);
      free(freq_den);
    }
  }

  //Clamping always happen at 0,be it left tructaion or right side tructaion!
  double evaluate_phi_i(double x,const uint8_t i){
    double x_old=PHI[i];
    assert(isnan(x)==0);

    PHI[i]=x; // updating  , but need to restore later to older value !
    double value=0;
    //    value= -0.5 *pow((x-phi_mu)/sqrt(phi_sigma_square),2) ;
    value= (x>0) ? -1*x/diversity[i] : x/diversity[i] ;
    assert(value<=0);
    double v=0; 
    //comes from likelihood  or data !
    for(uint32_t d=0;d<document_count;d++){
      document* doc=d_ptr+d ;
      double den=0;
      const int label= doc->true_label;
      const uint32_t N_d=doc->real_word_count;
      for(uint32_t s=0;s<topic_count;s++){
        den+= (PHI[s]*((double)(doc->topic_usage[s])/N_d));
      }
      v=exp(-1*label*den);
      //TAKE PRECAUTIONARY MEASURES
      //assert(isinf(v)==0);
      assert(isnan(v)==0);
      if(isinf(v)!=0){
        if(isinf(v)==1){
          //value+=(-1*log(1+very_big));
          value+=(-1*(-1*label*den));
        }else if(isinf(v)==-1){
          //value+=(-1*log(1+very_small));
          value+=(-1*log(1));
        }else{
          printf("some thing wrong here in eval_phi\n");
          exit(0);
        }
      }else{
        value+=(-1*log(1+exp(-1*label*den)));
      }
    }
    if(log_Fx_maxima_phi_i_possible<value){
      //set  maxima & value
      log_Fx_maxima_phi_i_possible=value;
      maxima_phi_i_possible=x;
    }
    PHI[i]=x_old;
    assert(isnan(value)==0);
    assert(isinf(value)==0);
    //    cout <<"i = " <<  i  << " x= " <<  x << "  value =  " << value << "\n";
    return value;
  }

  void  update_phi_i(const uint8_t i){
    double l=0;
    double r=0;
    double log_f_r=0;
    double log_f_l=0;
    bool done =false;
    bool  isLeftSideClamped= i==topic_count-1 ? true : false;
    bool  isRightSideClamped= !isLeftSideClamped;
    double    step_size_phi=0.0000000000001;
    double  maxima=0;
    double  log_F_max=0;

    log_Fx_maxima_phi_i_possible=-DBL_MAX ;
    get_max_f_phi_i(i,&log_F_max,&maxima);

    r=maxima;
    while(!done){
      r+=(step_size_phi) ;
      log_f_r=evaluate_phi_i(r,i);
      if( log_F_max-log_f_r > log_ratio+epsilon ){
        done=true;
      }
      step_size_phi= 2*step_size_phi;
    }

    l=maxima;
    done =false;
    step_size_phi=0.0000000000001;
    while(!done){
      l-=step_size_phi ;
      log_f_l=evaluate_phi_i(l,i);
      if( log_F_max-log_f_l > log_ratio+epsilon ){
        done=true;
      }
      step_size_phi=2*step_size_phi;
    }

    if(isLeftSideClamped){
      l=l<0 ? 0: l;
      log_f_l=evaluate_phi_i(l,i);
      //    assert(r>0);
      //World is not perfect
      if(r<0 || r<=l  ){
        //find new R
        r=l;
        done =false;
        step_size_phi=0.0000000000001;

        while(!done){
          r+=(step_size_phi) ;
          log_f_r=evaluate_phi_i(r,i);
          if( fabs(log_f_l-log_f_r) > 2*log_ratio+epsilon ){
            done=true;
          }
          step_size_phi= 2*step_size_phi;
        } // have a new right boundary
      }
      //reset maxima and log_F_max
      maxima=log_f_l > log_f_r ? l:r ;
      log_F_max= evaluate_phi_i(maxima,i);
    }

    if(isRightSideClamped){
      r=r >0 ? 0: r;
      log_f_r=evaluate_phi_i(r,i);
      // assert(l<0);
      if(l >=0 || l>=r  ){
        //find new R
        l=r;
        done =false;
        step_size_phi=0.0000000000001;

        while(!done){
          l-=(step_size_phi) ;
          log_f_l=evaluate_phi_i(l,i);
          if( fabs(log_f_l-log_f_r) > 2*log_ratio+epsilon ){
            done=true;
          }
          step_size_phi= 2*step_size_phi;
        } // have a new right boundary
      }
      //reset maxima and log_F_max
      maxima=log_f_l > log_f_r ? l:r ;
      log_F_max= evaluate_phi_i(maxima,i);

    }

    /*
       assert(l<=maxima);
       assert(r > l && r >= maxima);
       assert(log_F_max>=log_f_l && log_F_max >= log_f_r);
       */
    //Rejection sampling !
    bool sampled= false ;
    double sampled_number=0;
    double p_z0= 0;
    int reject_cnt=0;
    while(!sampled){
      double z0=get_number_between_two_numbers(l,r) ;
      double log_u0=  log(rand_0_1());// last 2 terms are kq(z0)
      p_z0= evaluate_phi_i(z0,i);
      if( (log_u0 < p_z0- log_F_max) ){
        //accept
        sampled_number=z0;
        sampled=true;
        //      printf("accepeted\n");
      }else{
        //    printf("rejected\n");
        reject_cnt++;
        if(reject_cnt>10000){
          printf("sample  was selected after max count reached to %d\n",reject_cnt);
          printf("Here  maxima=%.25f,log_F_max=%.25f \n",maxima,log_F_max);
          printf("value at lower=%.25f,log_f_l=%.25f\n",l,log_f_l) ;          
          printf("value at right=%.25f,log_f_r=%.25f\n",r,log_f_r) ;          
          sampled_number=maxima ;
          sampled=true;
        }
      }
    }
    //set to the new value !
    PHI[i]=sampled_number;
    //   printf("Reject count=%d\n ",reject_cnt) ;
    //    printf("PHI[%d]=%f  ",i,sampled_number);
  }

  void update_PHI(){
    for(uint8_t s=0;s<topic_count;s++){
      update_phi_i(s);
    }
    //  printf("\n");
  }

  double inline  get_number_between_two_numbers(double low,double  high){
    double   result=0;
    result=low+(high-low)*rand_0_1();
    return result;
  }

  void write_topic_state(){
    document* d_ptr1 = d_ptr ;
    FILE* fp;
    fp=fopen(TOPIC_file,"w");
    assert(fp!=NULL);

    printf("writting topic to disk!\n");
    for(uint32_t i=0;i<document_count;i++){
      d_ptr1=d_ptr+i;
      uint32_t sz=d_ptr1->real_word_count;
      uint8_t * buf=(uint8_t*)calloc(sz,sizeof(uint8_t));
      assert(buf!=NULL);
      memcpy(buf,d_ptr1->topic_list,sizeof(uint8_t)*sz);
      fwrite(buf,sizeof(uint8_t),sz,fp);
      free(buf);
    }
    fclose(fp);
  }

  void write_model(){
    FILE* fp;
    fp=fopen(THETA_file,"w");
    assert(fp!=NULL);
    fwrite(THETA,sizeof(double),topic_count*total_words_count,fp);
    fclose(fp);
    fp=fopen(M_file,"w");
    assert(fp!=NULL);
    fwrite(M,sizeof(double),topic_count*topic_count,fp);
    fclose(fp);
    fp=fopen(PHI_file,"w");
    assert(fp!=NULL);
    fwrite(PHI,sizeof(double),topic_count,fp);
    fclose(fp);
  }
  void write_model_gibbs(){
    
    cout <<  "writting model : LOOP_COUNT=  "  << loop  << "\n" ;
    cout << "GIBBS SAMPLER EXPRECTATION" << "\n";
    FILE* fp;
    fp=fopen("THETA_gibbs.bin","w");
    assert(fp!=NULL);
    fwrite(THETA_gibbs,sizeof(double),topic_count*total_words_count,fp);
    fclose(fp);
    fp=fopen("M_gibbs.bin","w");
    assert(fp!=NULL);
    fwrite(M_gibbs,sizeof(double),topic_count*topic_count,fp);
    fclose(fp);
    fp=fopen("PHI_gibbs.bin","w");
    assert(fp!=NULL);
    fwrite(PHI_gibbs,sizeof(double),topic_count,fp);
    fclose(fp);
  }


  double timetaken(clock_t start,clock_t end,const char* message){
    double  elapsed= ((long double)(end-start))/CLOCKS_PER_SEC;
    printf("Time taken is %.10f  %s\n", elapsed,message);
    return elapsed;
  }

  void m_compute_double_term_word_stat(const uint8_t i,double* sum_den ,uint32_t* freq_den,uint32_t* freq_num_i){
    for(uint32_t w=0;w<total_words_count;w++){
      double sum=0;
      for(uint8_t ti_prime=0;ti_prime<topic_count;ti_prime++){
        sum+=get_d_2D_m(M,i,ti_prime,topic_count) *get_d_2D_m(THETA,ti_prime,w,total_words_count);
      }
      assert(sum >=0);
      sum_den[w]=sum;
    }//

    for(int d=0;d<document_count;d++){
      document *doc=d_ptr+d;
      int real_word_count=doc->real_word_count;
      for(int j=1;j<real_word_count-1 ;j++){
        uint8_t t_last=doc->topic_list[j-1];
        uint8_t t_curr=doc->topic_list[j];
        if(t_last!=i  ){
          continue;
        }else{
          int word = doc->int_word_list[j];
          freq_den[word]++ ;
          freq_num_i[t_curr]++ ;
        }
      }
    }
  }

  void theta_compute_double_term_word_stat(const uint32_t w,double* sum_den,uint32_t* freq_den,uint32_t* freq_single_num_w,uint32_t* freq_double_num){
    for(uint8_t t_last=0;t_last<topic_count;t_last++){
      double sum=0;
      for(uint8_t  ti_prime=0;ti_prime<topic_count;ti_prime++){
        sum+=get_d_2D_m(M,t_last,ti_prime,topic_count)*get_d_2D_m(THETA,ti_prime,w,total_words_count);
      }
      assert(sum>=0);
      sum_den[t_last]=sum;
    }
    //Make numerator frequency in single term 
    for(uint32_t q2=0;q2<document_count;q2++){
      document* d=d_ptr+q2;
      if((d->int_word_list[0]==w)){
        uint32_t  t= d->topic_list[0];
        ( freq_single_num_w[t])++ ;
      }
    }

    for(int d=0;d<document_count;d++){
      document *doc=d_ptr+d;
      int real_word_count=doc->real_word_count;
      for(int j=1;j<real_word_count ;j++){
        int word = doc->int_word_list[j];
        if(word!=w){
          continue;
        }else{
          uint8_t t_last=doc->topic_list[j-1];
          uint8_t t_curr=doc->topic_list[j];
          freq_den[t_last]++;
          freq_double_num[t_curr]++;

        }
      }
    }
  }

  void update_TOPIC(){
    for(uint32_t i=0;i < document_count;i++){
      topic_update(d_ptr+i);
    }
  }

  void test_model(){
    int p2p=0;
    int p2n=0;
    int n2p=0;
    int n2n=0;
    int i=0;
    document* curr_doc= d_ptr;
    for(i=0;i<document_count;i++){
      curr_doc= d_ptr+i;
      int original_label= curr_doc->true_label;
      int  predicted_label= curr_doc-> sampled_label;
      if(original_label==1 &&  predicted_label==1){
        p2p++ ;
      }else if(original_label==1 &&  predicted_label==-1){
        p2n++;
      }else if(original_label==-1 &&  predicted_label==1){
        n2p++;
      }else if(original_label==-1 &&  predicted_label==-1){
        n2n++;
      }else{
        printf("some thing wrong has happened !\n ");
        exit(0);
      }

    }
    double precision=p2p/(double)(p2p+p2n);
    double recall=p2p/(double)(p2p+n2p);
    double  accuracy=(p2p+n2n)/(double)(p2p+p2n+n2n+n2p);

    double F1_score=2*precision*recall/(precision+recall);
    printf("p2p=%d ,p2n=%d,n2p=%d, n2n=%d \n",p2p,p2n,n2p,n2n);
    printf("precision=%f",precision);
    // printf(" recall=%f  ",recall);
    //  printf(" F1_score=%f  ",F1_score);
    printf(" accuracy=%f \n ",accuracy);
  }

  void print_model_params(){
    printf("printing M matrix  \n ");
    for(int i=0;i<topic_count;i++){
      for(int j=0;j<topic_count;j++){
        if(i!=j){
          //printf("M[%d][%d]= %.20f ",i,j,M[i][j]);
          printf("M[%d][%d]= %.20f ",i,j,get_d_2D_m(M,i,j,topic_count));
        }else{
          printf("**M[%d][%d]= %.20f** ",i,j,get_d_2D_m(M,i,j,topic_count));
        }
      }
      printf("\n");
    }
    printf("printing M matrix  \n ");

    printf("printing THETA matrix  \n ");
    for(int i=0;i<total_words_count;i++){
      for(int j=0;j<topic_count;j++){
        //printf("THETA[%d]->[%d]= %.20f ",i,j,THETA[j][i]);
        printf("THETA[%d]->[%d]= %.20f ",i,j,get_d_2D_m(THETA,j,i,total_words_count));
      }
      printf("\n");
    }
    printf("PHI vector\n");
    for(int  i=0;i<topic_count;i++)
      printf("PHI[%d]= %.20f ",i,PHI[i]);
    printf("\n");
  }


  void print_topic_corpus(){
    for(int i=0;i<document_count;i++){
      document* d=d_ptr+i ;
      int sz=d->real_word_count;
      printf("printing topics and words of document: %s\n ",d->document_name);
      printf("Word/Topic mapping  \n");
      for(int j=0;j<sz;j++){
        printf("word %d-> topic: %d \n",d->int_word_list[j],d->topic_list[j]) ;
      }
    }
  }

  void  free_things(){
    for(int  i=0;i<document_count;i++){
      free((d_ptr+i)-> topic_list);
      free((d_ptr+i)-> int_word_list);
      free((d_ptr+i)->priority_flag);
      free((d_ptr+i)->topic_usage);
    }
    free(d_ptr);
    gsl_rng_free(rand_obj);
    free(last_maxima_theta_iw);
    free(last_maxima_m_ij);
    free(last_maxima_phi);
    free(freq_single_den);
  }

  double gamma_approx_sampler(double r,double log_F_max,double log_f_r,double maxima){
    double shape=0;
    double scale=0;
    shape=1+(log_F_max-log_f_r)/((r-maxima)/maxima +(log(maxima)-log(r)));
    scale=maxima/(shape-1);
    assert(shape>0 && scale>0);
    double sampled_number=gsl_ran_gamma(rand_obj,shape,scale);
    return sampled_number;

  }

  void   guess_lebel(document* doc){
    const uint32_t N_d= doc->real_word_count;
    double sum=get_current_sum(doc,N_d,topic_count);
    double b_p= 1/(1+exp(-sum));
    double  rand_1= rand_0_1();
    doc->sampled_label= (rand_1 < b_p) ? 1 : -1;
  }

  double m_diagonistics(const int i,const int j,double* sum_den,uint32_t* freq_den,uint32_t* freq_num ,double *maxima_possible,double*value_possible,const char*msg , double maxima){
    printf("%s \n",msg);
    maxima=(maxima!=-1)?maxima :0.1 ;
    double lower=maxima/30;
    double upper=maxima*20;
    double step1=(upper-lower)/600;
    double uu=lower ;
    double found_maxima=maxima;
    double  vaL_maxima=evaluate_m_ij_ll(found_maxima,i,j,sum_den,freq_den,freq_num,maxima_possible,value_possible);
    double  found_val=vaL_maxima;
    printf("In m_diagnostics \n");
    printf("Maxima=%.15f,Value calculate at maxima=%.20f\n",maxima,vaL_maxima);

    for(int xx=0;uu<upper;xx++){
      uu=lower+xx*step1;
      double temp=evaluate_m_ij_ll(uu,i,j,sum_den,freq_den,freq_num,maxima_possible,value_possible);
      if(temp > found_val){
        found_val=temp;
        found_maxima=uu;
      }else{
      }
      printf("xx= %.20f ,val= %.20f \n",uu,temp);
    }
    return  found_maxima;
  }

  double theta_diagonistics(const int i, const int w,const double den_single_term_header,double * sum_den,uint32_t* freq_den,const uint32_t freq_single_num,uint32_t *  freq_num,const  uint32_t freq_single_den_w,double *possible_maxima,double * val_at_max,const char* msg,double maxima){
    printf("%s \n",msg);
    maxima=(maxima!=-1)?maxima :0.1 ;
    double lower=maxima/25;
    double upper=maxima*25;
    double step1=(upper-lower)/1000;
    double uu=lower ;
    double found_maxima=maxima;
    double  vaL_maxima=evaluate_theta_iw_ll(found_maxima, i, w, den_single_term_header, sum_den, freq_den, freq_single_num, freq_num,freq_single_den_w,possible_maxima,val_at_max);
    double  found_val=vaL_maxima;
    printf("-----------------------In THETA  DIAGONISTICS\n-----------------------");
    for(int xx=0;uu<upper;xx++){
      uu=lower+xx*step1;
      double temp=evaluate_theta_iw_ll(uu, i, w, den_single_term_header, sum_den, freq_den, freq_single_num, freq_num,freq_single_den_w,possible_maxima,val_at_max);
      if(temp > found_val){
        found_val=temp;
        found_maxima=uu;
      }else{
      }
      // printf("xx= %.20f ,val= %.20f \n",uu,temp);
    }
    return found_maxima;
  }

  double compute_theta_single_term_denominator(uint32_t w){
    double res=0;
    for(uint8_t g=0;g<topic_count;g++){
      res+=get_d_2D_m(THETA,g,w,total_words_count);
    }
    return res ;
  } 

  void print_system_config(){
    printf("---------SYSTEM CONFIG-----------\n");
    printf("2a setup\n");
    cout<< data << endl ;
    printf("TOPIC COUNT =%d,M_DIAG=%f\n",topic_count,m_diag);

  }

  void print_time(){
    time_t tim=time(NULL);
    char *s=ctime(&tim);
    s[strlen(s)-1]=0;        // remove \n
    printf("Time= %s \n", s);

  }

  void  create_word_id_id_word_map(string word1){
    string word=make_lowercase_string(word1);
    static int new_id=0;

    if(word.compare(comma)==0 || word.compare(period)==0 || word.compare(highlight_tag_start)==0 || word.compare(highlight_tag_end)==0){
      //No id for  , and . 
      return ;
    }
    map<string,int> ::const_iterator it = word_id_map.find(word);
    if(it==word_id_map.end()){
      word_id_map[word]=new_id;
      id_word_map[new_id]=word ;
      //cout << "word: " <<word << "  ID: " << new_id <<  "\n";
      new_id++ ;
    }

  }

  int  get_id_from_word(string word1){
    int id=-1;
    string word=make_lowercase_string(word1);
    map<string,int >::const_iterator it=word_id_map.find(word);
    //should always pass because ideally wh have created them as soon as we saw them for first time 
    if(it!= word_id_map.end()){
      id=it->second;
    }
    return id;
  }

  string get_word_from_id(int id){
    map<int,string >::const_iterator it=id_word_map.find(id);
    //should always pass because ideally wh have created them as soon as we saw them for first time 
    assert(it != id_word_map.end());
    return it->second;
  }

  string make_lowercase_string(string data){
    //cout <<  "word size: " << data.size() << "\n";
    char pch[1000] ; 
    memset(pch,'\0',1000);
    strcpy(pch,data.c_str());
    string res;
    for(int i=0;i<data.length();i++){
      pch[i]=tolower(pch[i]);
    }
    res=pch;
    return  res;
  }

  bool doesWordStartsUpper(string data){
    //cout <<  "word size: " << data.size() << "\n";
    bool isFirstCharUpper=false;
    char pch[1000] ; 
    memset(pch,'\0',1000);
    strcpy(pch,data.c_str());
    char firstChar=pch[0];
    if(isalpha(firstChar) && isupper(firstChar)){
      isFirstCharUpper=true;
    }
    return  isFirstCharUpper;
  }


  int main(){
    // exit(0);
    print_system_config(); 
    init();
    printf("Init  done \n ");
    for( loop=0;loop<loop_max;loop++){
      printf("Iteration=%d \n",loop+1);
      print_time();
      printf("TOPIC update\n");
      update_TOPIC();
      print_time();

      printf("THETA update\n");
      update_THETA();
      print_time();

      printf("M update\n");
      update_M();
      print_time();
      printf("PHI UPDATE\n");
      update_PHI();
      print_time();

      test_model();
      fflush(stdout);

      //GIBBS  expectation
      //do_gibbs_expectation();

      if(loop==loop_max-1){
        write_model();
        print_model_params();
        print_topic_corpus();
        //           write_model_gibbs();
      }
    }
    free_things();
    return 0;
  }
