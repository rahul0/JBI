#include "gibbs_linked_underlining.h"

void  make_Dictionary(){
  char buffer[2000];
  FILE* fp;
  fp= fopen(dict_file,"r");
  assert(fp!=NULL);
  memset(buffer,'\0',2000);
  char* pch ;
  while(fgets(buffer,2000,fp)!=NULL){
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
        create_word_id_id_word_map(word);
    }
    memset(buffer,0,2000);
  }//while ends here !      
  fclose(fp);  
}

void init(){
  d_ptr=(document *) calloc(document_count,sizeof(document));
  assert(d_ptr!=NULL);

  M=(double*)malloc(sizeof(double)*topic_count*topic_count);
  assert(M!=NULL);
  memset(M,0,sizeof(double)*topic_count*topic_count);
  THETA=(double*)malloc(sizeof(double)*topic_count*total_words_count);
  assert(THETA != NULL);
  memset(THETA,0,sizeof(double)*topic_count*total_words_count);
  PHI=(double*)malloc(sizeof(double)*topic_count);
  assert(PHI!=NULL);
  memset(PHI,0,sizeof(double)*topic_count);
  load_model_params(); // gets theta, m, phi params
  make_Dictionary();
  init_document_params();
  fill_multiplier_values();
}

void load_model_params(){
  FILE  *fp ;
  cout<< "Loaded model params" << endl ;
  fp=fopen(M_file,"r");
  assert(fp!=NULL);
  int res=  fread(M,sizeof(double),topic_count*topic_count,fp);
  assert(res==topic_count*topic_count);

  fp=fopen(THETA_file,"r");
  assert(fp!=NULL);
  res=  fread(THETA,sizeof(double),topic_count*total_words_count,fp);
  assert(res==topic_count*total_words_count);

  fp=fopen(PHI_file,"r");
  assert(fp!=NULL);
  res=  fread(PHI,sizeof(double),topic_count,fp);
  assert(res==topic_count);
  fclose(fp);
}

void init_document_params(){
  DIR *dir;
  struct dirent *ent ;
  dir=opendir(data);
  while((ent=readdir(dir))!=NULL){
    if(strchr(ent->d_name,'.')==NULL){
      FILE*  fp;
      char filename[2000];
      char *pch;
      document * current_doc_ptr=d_ptr+ (updated_doc++) ;
      //Extract name
      strcpy(filename,data);
      strcat(filename,"/");
      strcat(filename,ent->d_name);
      strcat(filename,"\0");
      current_doc_ptr->document_name=(char*)malloc(sizeof(char)*2000);
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
      memset(filename,'\0',2000);
      vector<string> w_list ;
      current_doc_ptr->word_list=w_list;
      while(fgets(filename,2000,fp)!=NULL){
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
      //make_link(current_doc_ptr);
    }//if ends here nn   ONE DOCUMENT
  } //while ends here 
  closedir(dir);
}

//This function makes word  and also sets flag that word belongs to high priority topic!
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
//        cout << "Parsing error" << endl ;
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
  doc->int_word_list=(uint32_t*)malloc(sizeof(uint32_t)*sz) ;
  assert(doc->int_word_list!=NULL);
  memset(doc->int_word_list,0,sizeof(uint32_t)*sz);

  doc->priority_flag=(bool*)malloc(sizeof(bool)*sz);
  assert(doc->priority_flag!=NULL);
  memset(doc->priority_flag,false,sizeof(bool)*sz);

  doc->comma_present_before=(bool*)malloc(sizeof(bool)*sz);
  assert(doc->comma_present_before!=NULL);
  memset(doc->comma_present_before,false,sizeof(bool)*sz);

  doc->period_present_before=(bool*)malloc(sizeof(bool)*sz);
  assert(doc->period_present_before!=NULL);
  memset(doc->period_present_before,false,sizeof(bool)*sz);

  doc->word_caps=(bool*)malloc(sizeof(bool)*sz);
  assert(doc->word_caps!=NULL);
  memset(doc->word_caps,false,sizeof(bool)*sz);

  for(int i=0;i<sz;i++){
    //integer word list
    doc->int_word_list[i]=word_index[i];
    //priority_flag 
    doc->priority_flag[i]=priority_flag_vec[i];

    if(delimiter_vec[i].compare(comma)==0){
      doc->comma_present_before[i]=true;
    }else if(delimiter_vec[i].compare(period)==0){
      doc->period_present_before[i]=true;
    }else{
      // means delimiter was white space 
    } 
    doc->word_caps[i]=word_caps_vec[i];
  }
  //clear every thing 
  //
  word_index.clear();
  word_caps_vec.clear();
  priority_flag_vec.clear();
  delimiter_vec.clear();
//  doc->word_list.clear();
  space.clear();
 delimiter.clear();
 word.clear();

}


void init_topic_assignment(document * doc){
    time_t seconds ;
  time(&seconds);
  uint8_t  topic=0 ;
  uint32_t N_d=doc->real_word_count;
  doc->topic_list= (uint8_t*)malloc(N_d*sizeof(uint8_t));
  memset(doc->topic_list,0,N_d*sizeof(uint8_t));
  assert(doc->topic_list!=NULL);
  for(uint32_t w=0;w<N_d;w++){
    const uint32_t word = doc->int_word_list[w];
    double  *pdf = (double*)calloc(topic_count,sizeof(double));
    double  *cdf = (double*)calloc(topic_count,sizeof(double));
    double  cdf_sum=0;
    for(int ii=0;ii<topic_count;ii++){
      if(w==0){
        pdf[ii]=get_d_2D_m(THETA,ii,word,total_words_count) ;
      }else{
        uint8_t t_last=doc->topic_list[w-1];
        //different numerator
        if(ii==t_last){
          double  punctuation_multiple1_num=get_punctuation_multipler_num(doc,w,t_last,ii);
          pdf[ii]=punctuation_multiple1_num*get_d_2D_m(THETA,ii,word,total_words_count)*get_d_2D_m(M,t_last ,ii,topic_count);
        }else{
          pdf[ii]=get_d_2D_m(THETA,ii,word,total_words_count)*get_d_2D_m(M,t_last ,ii,topic_count);
        }
      }
      cdf_sum+=pdf[ii];
    }
    scale_pdf(pdf,topic_count,cdf_sum);
    make_cdf_from_pdf(pdf,cdf,topic_count);
    double cdf_rand=rand_0_1();
    topic=linear_scan(cdf,topic_count,cdf_rand) ;
    free(pdf);
    free(cdf);
    doc->topic_list[w] = topic;
    const int old_topic=-1;
    update_topic_usage(doc,topic,old_topic);
  }
//  cout << endl ;
    double sum=get_current_sum(doc,N_d,topic_count);//This is same
    double b_p= 1/(1+exp(-1*sum));
    assert(b_p >0);
    assert(isnan(b_p)==0);
    assert(isinf(b_p)==0);
    double  rand_1= rand_0_1();
    doc->sampled_label= ( b_p>rand_1) ? 1 : -1;
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


double get_factor12(const int word_curr,const int word_next,const int t_prime,const int t_last,const int t_next,document* doc,const uint32_t mypos){
  bool isFirstWord=mypos==0 ? true :false ;
  bool isLastWord=word_next==-1 ? true :false ;
  double result =0;
  if(!isLastWord && !isFirstWord){
    double num=get_d_2D_m(M,t_last,t_prime,topic_count) *  get_d_2D_m(M,t_prime,t_next,topic_count)* get_d_2D_m(THETA,t_prime,word_curr,total_words_count);
    assert(num>0);
    double  punctuation_multiple1_num= get_punctuation_multipler_num(doc,mypos,t_last,t_prime);

    if(t_last==t_prime){
      num=num*punctuation_multiple1_num;
    }
    //effect of punctuation from next words topic
    double  punctuation_multiple2_num= get_punctuation_multipler_num(doc,mypos+1,t_prime,t_next);
    if(t_prime==t_next){
      num=num*punctuation_multiple2_num;
    }
    double den=0;
    for(uint8_t i=0;i<topic_count;i++){
      if(i==t_prime){
        double  punctuation_multiple2_den= get_punctuation_multipler_den(doc,mypos+1,t_prime,t_next);
        den+=punctuation_multiple2_den*get_d_2D_m(M,t_prime,i,topic_count) * get_d_2D_m(THETA,i,word_next,total_words_count);
      }else{
        den+=get_d_2D_m(M,t_prime,i,topic_count) * get_d_2D_m(THETA,i,word_next,total_words_count);
      }
    }
    result=num/den;
  assert(result>0);
    }else if (isFirstWord){
      double num=get_d_2D_m(THETA,t_prime,word_curr,total_words_count)* get_d_2D_m(M,t_prime,t_next,topic_count);
    assert(num>0);
      if(t_prime==t_next){
        double  punctuation_multiple2_num= get_punctuation_multipler_num(doc,mypos+1,t_prime,t_next);
        num=num*punctuation_multiple2_num;
      }
      double den=0;
      for(uint8_t i=0;i<topic_count;i++){
        if(i==t_prime){
          double  punctuation_multiple2_den= get_punctuation_multipler_den(doc,mypos+1,t_prime,t_next);
          den+=punctuation_multiple2_den*get_d_2D_m(M,t_prime,i,topic_count) * get_d_2D_m(THETA,i,word_next,total_words_count);
        }else{
          den+=get_d_2D_m(M,t_prime,i,topic_count) * get_d_2D_m(THETA,i,word_next,total_words_count);
        }
      }
      result=num/den;
  assert(result>0);
    }else if(isLastWord){
      double num=get_d_2D_m(THETA,t_prime,word_curr,total_words_count)* get_d_2D_m(M,t_last,t_prime,topic_count);
    assert(num>0);
      if(t_last==t_prime){
        double  punctuation_multiple1_num= get_punctuation_multipler_num(doc,mypos,t_last,t_prime);
        num=num*punctuation_multiple1_num;
      }
      result=num;
  assert(result>0);
    }else{
      printf("Error state !\n");
      exit(0);
    }
    return result;
  }

inline double  get_factor2_t_i(document* doc,double sum){
    double res=0;
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
    for(uint32_t i=0;i<N_d;i++)//update all word-topic association !
    {
      const uint32_t word_curr=doc->int_word_list[i];
      const uint8_t before_update= doc->topic_list[i];
      const int  t_last= ((i!=0 )? doc->topic_list[i-1]:-1);
      const int t_next=((i!=(N_d-1)) ? doc->topic_list[i+1]:-1);//set next to one if we are on last word 
      const int word_next=(i==(N_d-1))?-1 :doc->int_word_list[i+1];
      double delta_minus = PHI[before_update]/(double)N_d;
      int  new_topic=-1;
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
          sum_cdf+=t_prime_pdf[j];
        }
        scale_pdf(t_prime_pdf,topic_count,sum_cdf);
        make_cdf_from_pdf(t_prime_pdf,t_prime_cdf,topic_count);
        double cdf_rand=rand_0_1();
        new_topic=linear_scan(t_prime_cdf,topic_count,cdf_rand) ;
      
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


void inline set_d_2D_m(double *ptr,int row,int col,double value,int row_length){
  ptr[row_length*row+col]=value;
}

double inline get_d_2D_m(double* ptr,int row,int col,int row_length){
  assert(row>= 0 && col >=0 && row_length >0 );
  return *(ptr+(row_length)*row+col);
}

//double  get_d_3D_m(const double *ptr ,int x,int y,int z,int y_size,int z_size){
double  get_f_3D_m(const double *ptr ,int x,int y,int z,int y_size,int z_size){
  assert(x>= 0 && y >=0 && z >=0 && y_size >0 &&  z_size >0 );
  return *(ptr+ x*y_size*z_size+y*z_size+z ) ;
}

//void  set_d_3D_m(double *ptr,int x,int y,int z,double  value,int y_size,int z_size){
void  set_f_3D_m(double *ptr,int x,int y,int z,double  value,int y_size,int z_size){
  assert(x>= 0 && y >=0 && z >=0 && z_size >0 &&  y_size >0 );
  // assert(value >=0 );
  *(ptr+ x*y_size*z_size+y*z_size+z )=value ;
}

char get_c_3D_m(const char *ptr ,int x,int y,int z,int y_size,int z_size){
  assert(x>= 0 && y >=0 && z >=0 && y_size >0 &&  z_size >0 );
  return *(ptr+ x*y_size*z_size+y*z_size+z ) ;
}

void  set_c_3D_m( char* ptr,int x,int y,int z,char  value,int y_size,int z_size){
  assert(x>= 0 && y >=0 && z >=0 && z_size >0 &&  y_size >0 );
  *(ptr+ (x*y_size*z_size+y*z_size+z) )=value ;
}

void print_system_config(){
  printf("---------SYSTEM CONFIG-----------\n");
  printf("topic_count=%d,m_diag=%d\n",topic_count,m_diag);
  cout<< "Make sure you  have set ALPHA params correctly" <<  endl;
}


//exponential sum 
double elsum(double ln_x , double ln_y){
  double res=0;
  assert(isnan(ln_y)==0 && isnan(ln_x)==0);
  if(ln_x==0 || ln_y==0 ){
    res= 0;
  }else{
    res=  ln_x + ln_y ;
  }
  assert(isnan(res)==0 && isinf(res)==0);
  return res ;
}

/*
  void  reconstruct_doc(document *doc){
    ofstream myfile;
    string original_name= doc->document_name ;
    string file_name=write_dir+"/"+original_name+".html" ;
    myfile.open(file_name.c_str());
    string html_header="<!DOCTYPE html> <html> <body>";
    myfile << html_header << endl ;
    string html_heading="<h1>" +original_name+ "</h1>";
    myfile << html_heading << endl ;
    int N_d=doc->real_word_count;
    string start=" <font color=\"red\"><b> " ;  
    string end=" </font></b> ";  
    double match=0;
    double special_topic_occur=0;
    for(int i=0;i<N_d;i++){
      int id= doc->int_word_list[i] ;
      map<int,string>::iterator it =id_word_map.find(id);
      if(doc->period_present_before[i]){
        myfile << " . "  ;
      }
      if(doc->comma_present_before[i]){
        myfile << " , "  ;
      }
      string word; 
      if(it!=id_word_map.end()){
        word= it->second;
      }else{
        cout << "Error in  finding map for  id= " <<id <<  " \n" ;
        exit(0);
      } //got word !
      //If word has flag set for upper case then show as upper case
      if(doc->word_caps[i]){
       word=make_uppercase_string(word);
      }
      //decide how to append or prepend 
      if( doc->topic_list[i]==topic_count-1){ 
        ++special_topic_occur;
        if(doc->priority_flag[i]==true ){
        match++ ;
        }
        word=  start + " "+ word + " "+ end  ;
      }
      myfile << word << " " ;
    word.clear();
    }
    myfile << "\n" ;
    string html_end_tag="</p> </body> </html>";
    myfile << html_end_tag ;
    myfile.close();
//    cout <<  doc->document_name << " = " << (double)match/(doc->mygamma)*100 <<  " GAMMA_ESTIMATED: " << special_topic_occur <<  "  GAMMA_TRUE= " << doc->mygamma <<  "     ratio=  " << (double)special_topic_occur/(doc->mygamma)*100    << " precision underlining="  << (double)match/special_topic_occur << endl ;
    cout <<  doc->document_name << " = " << (double)match/(doc->mygamma)*100 <<  "  ratio(predicted/TRUE)GAMMA=  " << (double)special_topic_occur/(doc->mygamma)*100  << " precision underlining="  << (double)match/special_topic_occur*100 << endl ;
   original_name.clear();
   file_name.clear();
   html_end_tag.clear();
   html_heading.clear();
   start.clear();
   end.clear();
   html_end_tag.clear();
  
  }
*/


  void  reconstruct_doc(document *doc){
    int tp=0;
    int fn=0;
    int fp=0;
    int tn=0;

    ofstream myfile;
    string original_name= doc->document_name ;
    string file_name=write_dir+"/"+original_name+".html" ;
    myfile.open(file_name.c_str());
    string html_header="<!DOCTYPE html> <html> <body>";
    myfile << html_header << endl ;
    string html_heading="<h1>" +original_name+ "</h1>";
    myfile << html_heading << endl ;
    int N_d=doc->real_word_count;
    string start=" <font color=\"red\"><b> " ;  
    string end=" </font></b> ";  
    double match=0;
    double special_topic_occur=0;
    for(int i=0;i<N_d;i++){
      int id= doc->int_word_list[i] ;
      map<int,string>::iterator it =id_word_map.find(id);
      if(doc->period_present_before[i]){
        myfile << " . "  ;
      }
      if(doc->comma_present_before[i]){
        myfile << " , "  ;
      }
      string word; 
      if(it!=id_word_map.end()){
        word= it->second;
      }else{
        cout << "Error in  finding map for  id= " <<id <<  " \n" ;
        exit(0);
      } //got word !
      //If word has flag set for upper case then show as upper case
      if(doc->word_caps[i]){
       word=make_uppercase_string(word);
      }
      //decide how to append or prepend 
      bool isWordUnderLined_Really= doc->priority_flag[i];
      bool isWordUnderLined_Prediction=doc->topic_list[i]==topic_count-1 ? true: false;
      if( isWordUnderLined_Prediction ==true &&  isWordUnderLined_Really==true ){
        ++tp;
        word=  start + " "+ word + " "+ end  ;
        }else if(isWordUnderLined_Prediction ==true &&  isWordUnderLined_Really==false ){
        ++fp;

        word=  start + " "+ word + " "+ end  ;
        }else if(isWordUnderLined_Prediction ==false &&  isWordUnderLined_Really==false)
        {
       ++ tn ;
        }else{
       ++ fn;
        }


      myfile << word << " " ;
      word.clear();
    }
    myfile << "\n" ;
    string html_end_tag="</p> </body> </html>";
    myfile << html_end_tag ;
    myfile.close();
    //cout <<  doc->document_name << " = " << (double)match/(doc->mygamma)*100 <<  "  ratio(predicted/TRUE)GAMMA=  " << (double)special_topic_occur/(doc->mygamma)*100  << " precision underlining="  << (double)match/special_topic_occur*100 << endl ;
    double accuracy= (double)(tp+tn)/(tp+tn+fp+fn)*100 ;
    double precision= (double)tp/(tp+fp);
    double recall=(double)tp/(tp+fn) ;
    double F= 2*(recall*precision)/(precision+recall);

    cout <<  doc->document_name << " = " << accuracy <<   "          precision= " << precision  << " recall: "  << recall<< " F: "<< F << 
     "tp+fn " << tp+fn  << " GAMMA =" << doc->mygamma << endl ;
   original_name.clear();
   file_name.clear();
   html_end_tag.clear();
   html_heading.clear();
   start.clear();
   end.clear();
   html_end_tag.clear();
  
  }
bool doesWordStartsUpper(string data){
    bool isFirstCharUpper=false;
    char pch[2000] ; 
    memset(pch,'\0',2000);
    strcpy(pch,data.c_str());
    char firstChar=pch[0];
    if(isalpha(firstChar) && isupper(firstChar)){
      isFirstCharUpper=true;
    }
    return  isFirstCharUpper;
  }

double get_punctuation_multipler_num(document* doc, int pos,int t_last, int t_next ){
    double factor=1;
    if(t_last==t_next){
      factor=doc->comma_present_before[pos] ? factor*alpha_comma_Y : factor*alpha_comma_N ;
      factor=doc->period_present_before[pos] ? factor*alpha_period_Y : factor*alpha_period_N ;
      factor=doc->word_caps[pos] ?  factor*alpha_upper_Y : factor*alpha_upper_N ;
    }else{

    }
    return factor;
  }

double get_punctuation_multipler_den(document* doc, int pos,int t_last, int t_next ){
    double factor=1;
    factor=doc->comma_present_before[pos] ? factor*alpha_comma_Y : factor*alpha_comma_N ;
    factor=doc->period_present_before[pos] ? factor*alpha_period_Y : factor*alpha_period_N ;
    factor=doc->word_caps[pos] ?  factor*alpha_upper_Y : factor*alpha_upper_N ;
    return factor;
  }

void fill_multiplier_values(){
    for(int i=0 ;i<total_configs_punctuation;i++ ){
      double val=-1;
      switch(i){
        case 0 : val=alpha_comma_N*alpha_upper_N*alpha_period_N;
                 break;
        case 1 : val=alpha_comma_N*alpha_upper_N*alpha_period_Y;
                 break;
        case 2 : val=alpha_comma_N*alpha_upper_Y*alpha_period_N;
                 break;
        case 3 : val=alpha_comma_N*alpha_upper_Y*alpha_period_Y;
                 break;
        case 4 : val=alpha_comma_Y*alpha_upper_N*alpha_period_N;
                 break;
        case 5 : val=alpha_comma_Y*alpha_upper_N*alpha_period_Y;
                 break;
        case 6 : val=alpha_comma_Y*alpha_upper_Y*alpha_period_N;
                 break;
        case 7 : val=alpha_comma_Y*alpha_upper_Y*alpha_period_Y;
                 break;
        default :  break;  
      }
      multiplier_values[i]=val;
      assert(multiplier_values[i]!=-1);
    }
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
    new_id++ ;
  }
  //should never fail
  assert(new_id <= total_words_count);
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

string get_word_from_id(const int id){
    map<int,string >::const_iterator it=id_word_map.find(id);
    //should always pass because ideally wh have created them as soon as we saw them for first time 
    assert(it != id_word_map.end());
    return it->second;
  }

string make_lowercase_string(string data){
   char pch[2000] ; 
    memset(pch,'\0',2000);
    strcpy(pch,data.c_str());
    string res;
    for(int i=0;i<data.length();i++){
      pch[i]=tolower(pch[i]);
    }
    res=pch;
    return  res;
  }

string make_uppercase_string(string data){
    char pch[2000] ; 
    memset(pch,'\0',2000);
    strcpy(pch,data.c_str());
    string res;
    for(int i=0;i<data.length();i++){
      pch[i]=toupper(pch[i]);
    }
    res=pch;
    return  res;
  }
  void test_model(){
    int tp=0;
    int fn=0;
    int fp=0;
    int tn=0;
    int i=0;
    document* curr_doc= d_ptr;
    for(i=0;i<document_count;i++){
      curr_doc= d_ptr+i;
      int original_label= curr_doc->true_label;
      int  predicted_label= curr_doc-> sampled_label;
      if(original_label==1 &&  predicted_label==1){
        tp++ ;
      }else if(original_label==1 &&  predicted_label==-1){
        fn++;
      }else if(original_label==-1 &&  predicted_label==1){
        fp++;
      }else if(original_label==-1 &&  predicted_label==-1){
        tn++;
      }else{
        printf("some thing wrong has happened !\n ");
        exit(0);
      }

    }
    double  accuracy=(tp+tn)/(double)(tp+fp+tn+fn);
    double precision=tp/(double)(tp+fp);
    double recall=tp/(double)(tp+fn);
    /*
   if(loop>50 && accuracy > 0.75){
   for(int i=0;i < document_count;i++){
     // topic_update(d_ptr+i);
      reconstruct_doc(d_ptr+i);
    }
   exit(0);

   }
*/
    double F1_score=2*precision*recall/(precision+recall);
    printf("tp=%d ,fn=%d,fp=%d, tn=%d \n",tp,fn,fp,tn);
    printf("precision=%f",precision);
    printf(" recall=%f  ",recall);
    printf(" F1_score=%f  ",F1_score);
    printf(" accuracy=%f \n ",accuracy);
  }


  void document_fee_mem(document* d_curr){
    free(d_curr->word_caps);
    free(d_curr->comma_present_before);
    free(d_curr->period_present_before);
    free(d_curr->int_word_list);
    free(d_curr->priority_flag);
    free(d_curr->topic_usage);
    free(d_curr->document_name);
    free(d_curr);
  }


void print_time(){
    time_t tim=time(NULL);
    char *s=ctime(&tim);
    s[strlen(s)-1]=0;        // remove \n
    printf("Time= %s \n", s);

  }

void update_TOPIC(){
    for(int i=0;i < document_count;i++){
      topic_update(d_ptr+i);
    }
  }

inline double rand_0_1(){
    return (double)rand()/(double)RAND_MAX;
  }

int main(){
  //  print_system_config();
alpha_upper_N=1;
alpha_upper_Y=1;
alpha_comma_N=1;
alpha_comma_Y=1;
alpha_period_Y=1;
alpha_period_N=1;

  init();
    for( loop=0;loop<loop_max;loop++){
      update_TOPIC();
    }

    for(int i=0;i < document_count;i++){
     // topic_update(d_ptr+i);
      reconstruct_doc(d_ptr+i);
    }

}
