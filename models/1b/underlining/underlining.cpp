#include "underlining.h"

void  make_Dictionary(){
  char buffer[1000];
  FILE* fp;
  fp= fopen(dict_file,"r");
  assert(fp!=NULL);
  memset(buffer,'\0',1000);
  char* pch ;
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
        create_word_id_id_word_map(word);
    }
    memset(buffer,0,1000);
  }//while ends here !      
  fclose(fp);  
}

void init(){
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
}

void load_model_params(){
  FILE  *fp ;

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

void test_document(string directory_name){
//  static bool magic=false;
  cout << "=TESTING FOR DATASET: " << directory_name << endl ; 
  DIR *dir;
  struct dirent *ent ;
  const char*data=directory_name.c_str();
  dir=opendir(data);
  while((ent=readdir(dir))!=NULL){
    if(strchr(ent->d_name,'.')==NULL){
      FILE*  fp;
      char filename[1000];
      char *pch;
     // document * current_doc_ptr=d_ptr+ (updated_doc++) ;
      document * current_doc_ptr=(document*)malloc(sizeof(document)) ;
      //Extract name
      strcpy(filename,data);
      strcat(filename,"/");
      strcat(filename,ent->d_name);
      strcat(filename,"\0");
      current_doc_ptr->document_name=(char*)malloc(sizeof(char)*1000);
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

      assert(current_doc_ptr->document_name!=NULL);
      current_doc_ptr->hasPriorityWord=false;
      current_doc_ptr->topic_usage=(uint32_t*)malloc(sizeof(uint32_t)*(topic_count));
      assert(current_doc_ptr->topic_usage!=NULL);
      memset(current_doc_ptr->topic_usage,0,sizeof(uint32_t)*(topic_count));
 
      /*
      if( magic ||strcmp("/home/rahul/projects/text_mining/branch/punctuation/chris_data/F/p_71",filename)!=0){
        magic=true;
      continue ;
      }
      */
      fp= fopen(filename,"r");
      assert(fp!=NULL);
      memset(filename,'\0',sizeof(char)*1000);
      vector<string> w_list ;
      while(fgets(filename,1000,fp)!=NULL){
      assert(current_doc_ptr->document_name!=NULL);
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
      run_test(current_doc_ptr);
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
      //  cout << "Parsing error" << endl ;
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
  memset(doc->int_word_list,0,sz);

  doc->priority_flag=(bool*)malloc(sizeof(bool)*sz);
  assert(doc->priority_flag!=NULL);
  memset(doc->priority_flag,false,sz);

  doc->comma_present_before=(bool*)malloc(sizeof(bool)*sz);
  assert(doc->comma_present_before!=NULL);
  memset(doc->comma_present_before,false,sz);

  doc->period_present_before=(bool*)malloc(sizeof(bool)*sz);
  assert(doc->period_present_before!=NULL);
  memset(doc->period_present_before,false,sz);

  doc->word_caps=(bool*)malloc(sizeof(bool)*sz);
  assert(doc->word_caps!=NULL);
  memset(doc->word_caps,false,sz);

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

void inline set_d_2D_m(double *ptr,int row,int col,double value,int row_length){
  ptr[row_length*row+col]=value;
}

double inline get_d_2D_m(double* ptr,int row,int col,int row_length){
  assert(row>= 0 && col >=0 && row_length >0 );
  return *(ptr+(row_length)*row+col);
}

double  get_d_3D_m(const double *ptr ,int x,int y,int z,int y_size,int z_size){
  assert(x>= 0 && y >=0 && z >=0 && y_size >0 &&  z_size >0 );
  return *(ptr+ x*y_size*z_size+y*z_size+z ) ;
}

void  set_d_3D_m(double *ptr,int x,int y,int z,double  value,int y_size,int z_size){
  assert(x>= 0 && y >=0 && z >=0 && z_size >0 &&  y_size >0 );
  // assert(value >=0 );
  *(ptr+ x*y_size*z_size+y*z_size+z )=value ;
}
double  get_i_3D_m(const int *ptr ,int x,int y,int z,int y_size,int z_size){
  assert(x>= 0 && y >=0 && z >=0 && y_size >0 &&  z_size >0 );
  return *(ptr+ x*y_size*z_size+y*z_size+z ) ;
}

void  set_i_3D_m(int *ptr,int x,int y,int z,double  value,int y_size,int z_size){
  assert(x>= 0 && y >=0 && z >=0 && z_size >0 &&  y_size >0 );
  // assert(value >=0 );
  *(ptr+ x*y_size*z_size+y*z_size+z )=value ;
}

void print_system_config(){
  printf("---------SYSTEM CONFIG-----------\n");
  printf("----MODEL TEXTING\n---");
  printf("topic_count=%d,m_diag=%d\n",topic_count,m_diag);
  cout<< "Make sure you  have set ALPHA params correctly" <<  endl;
}

void initialize_memory_for_DP(document*  doc){
  int N_d = doc->real_word_count;
  assert(N_d > 0);
  int sz=N_d*topic_count*GAMMA;
  F= (double*)malloc(sizeof(double)*sz);
  assert(F!=NULL);
  back_ptr= (int*)malloc(sizeof(int)*sz);
  assert(back_ptr!=NULL);
  for(int i=0 ; i<sz;i++){
    //In the begining everthing is uninitialised !
    F[i]=UNINITIALIZED;
    back_ptr[i]=-1;
  }
}

void initialize_base_case_for_DP(document *doc){
  int N_d=doc->real_word_count;
  for(int word_index=1;word_index<=N_d; word_index++){
    for(int j=0;j<topic_count;j++ ){
      bool isTopic_in_U=(j==topic_count-1) ? true :false;
      for(int gamma=0;gamma<GAMMA ; gamma++){
        double log_val=0;
        int i =word_index-1;
        if(isBaseCase(word_index,j,gamma )){
        if(word_index==1){  // this translated to  i=1  in paper
          switch(gamma){
            case 1  : log_val= isTopic_in_U ?  get_fjwTheta(j,doc->int_word_list[0]):0;
                      break;
            case 0  : log_val= !isTopic_in_U ?  get_fjwTheta(j,doc->int_word_list[0]) :0 ;
                      break;
            default: log_val= 0;
                     break;
          }
        }// only when you are first word 
        assert(isinf(log_val)==0 &&  isnan(log_val)==0);
        set_d_3D_m(F,i,j,gamma,log_val,topic_count,GAMMA);
        }
        //       cout << "F[" <<word_index << "," <<  j << "," << gamma << "]= " << log_val << "\n" ;
      }
    }
  }
}

double get_fjwTheta(int j, int w ){
  double res=0;
  assert( j>=0 &&j<topic_count);
  assert(w<total_words_count);
  double val=get_d_2D_m(THETA,j,w,total_words_count);
  double sum=0;
  for(int j_prime=0; j_prime< topic_count;j_prime++){
    sum+=get_d_2D_m(THETA,j_prime,w,total_words_count);
  }
res=val/sum ;
assert(res>0);
  return log(res);

}


double under_lining(document* doc,int i, int j , int gamma){
  int N_d=doc->real_word_count;
  //cout <<  "calling underlining(" << i <<"," << j  <<"," << gamma << ")\n" ;
//  int w_i=doc->int_word_list[i];
  int word_index=i+1;
  //  if(i==0 || gamma > i){
  if( isBaseCase(word_index,j,gamma) )
     {
    //BASE CASE
    double log_val=get_d_3D_m(F,i,j,gamma,topic_count,GAMMA) ;
    assert(isinf(log_val)==0 &&  isnan(log_val)==0);
    //cout <<  "log_value of log(F[ " <<word_index <<"," << j <<","<< gamma << "])=" << log_val << "\n";
    return log_val;
  }else{
    bool isTopic_in_U=(j==topic_count-1) ? true :false;
    if(isTopic_in_U){
      double max_log_val=-DBL_MAX; 
      double max_index= -1;
      for(int k=0;k<topic_count;k++){
        double temp=get_d_3D_m(F,i-1,k,gamma-1,topic_count,GAMMA);
        bool isInitialised= isnan(temp)==0 ? true:false;
        if(!isInitialised){
          temp=under_lining(doc,i-1,k,gamma-1); 
          assert(isnan(temp)==0 && isinf(temp)==0);
          set_d_3D_m(F,i-1,k,gamma-1,temp,topic_count,GAMMA);
          //  cout <<  "log_value of log(F[" <<i <<"," << k <<","<< gamma-1 << "])=" << temp << "\n";
        }
        //double res=elsum(temp,get_d_3D_m(factor_matrix_log,k,j,w_i,topic_count,N_d)); 
        double res=elsum(temp,get_d_3D_m(factor_matrix_log,k,j,i,topic_count,N_d)); 
        if(res==0){ // why , because everyone is one log scale , hence not a fair comparision !
          continue ;
        }else{
          max_index=max_log_val>=res ? max_index:k;
          max_log_val=max_log_val>=res ? max_log_val:res;
          assert(isinf(res)==0 && isnan(res)==0);
          assert(isinf(max_log_val)==0 && isnan(max_log_val)==0);
        }
      }
      // cout <<  "log_value of log(F[" <<i+1 <<"," << j <<","<< gamma << "])=" << max_log_val << "\n";
      assert(max_log_val!=0);
      set_d_3D_m(F,i,j,gamma,max_log_val,topic_count,GAMMA);
      assert(max_log_val!=0);
      if(max_index==-1){
        cout << "Debug me !" <<  "\n";
      }
      assert(max_index!=-1);
      set_i_3D_m(back_ptr,i,j,gamma,max_index,topic_count,GAMMA);
      //      cout << "max_index " << max_index  << "\n" ;
      assert(isinf(max_log_val)==0 &&  isnan(max_log_val)==0);
      return max_log_val;
    }else{//NOT IN U
      double max_log_val=-DBL_MAX; 
      double max_index= -1;
      for(int k=0;k<topic_count;k++){
        double temp=get_d_3D_m(F,i-1,k,gamma,topic_count,GAMMA);
        bool isInitialised= isnan(temp)==0 ? true:false;
        if(!isInitialised){
          temp=under_lining(doc,i-1,k,gamma); 
          assert(isinf(temp)==0 &&  isnan(temp)==0);
          set_d_3D_m(F,i-1,k,gamma,temp,topic_count,GAMMA);
          // cout <<  "log_value of log(F[" <<i <<"," << k <<","<< gamma << "])=" << temp << "\n";
        }
        //double  res=elsum(temp,get_d_3D_m(factor_matrix_log,k,j,w_i,topic_count,N_d)); 
        double  res=elsum(temp,get_d_3D_m(factor_matrix_log,k,j,i,topic_count,N_d)); 
        if(res==0){ // why , because everyone is one log scale , hence not a fair comparision !
          continue ;
        }else{
          max_index=res > max_log_val ? k : max_index ;
          max_log_val=res > max_log_val ? res:max_log_val;
          assert(isinf(res)==0 && isnan(res)==0);
          assert(isinf(max_log_val)==0 && isnan(max_log_val)==0);
        }
      }
      // cout <<  "log_value of log(F[" <<i+1 <<"," << j <<","<< gamma << "])=" << max_log_val << "\n";
       assert(max_log_val!=0);
      set_d_3D_m(F,i,j,gamma,max_log_val,topic_count,GAMMA);
      if(max_index==-1){
        cout << "Debug me !" <<  "\n" ;
      }
      assert(max_index!=-1);
      set_i_3D_m(back_ptr,i,j,gamma,max_index,topic_count,GAMMA);
      //    cout << "max_index " << max_index << "\n" ;
      assert(isinf(max_log_val)==0 &&  isnan(max_log_val)==0);
      return max_log_val;
    }
  }//else ends here 
}//function ends here \n!

//exponential sum 
double elsum(double ln_x , double ln_y){
  double res=0;
  assert(isnan(ln_y)==0 && isnan(ln_x)==0);
  if(ln_x==0 || ln_y==0 ){
    res= 0;
  }else{
    //assert( isnan(ln_x + ln_y)==0 && isinf(ln_x + ln_y)==0);
    res=  ln_x + ln_y ;
  }
  assert(isnan(res)==0 && isinf(res)==0);
  return res ;
}


void make_factor_matrix_log(document *doc){
  int N_d=doc->real_word_count;
  for(int i=0;i<topic_count;i++){
    for(int j=0;j<topic_count;j++){
      for(int w=0;w<N_d;w++){
        int word=doc->int_word_list[w];
        double den=0;
        int t_last=i;
        int t_prime=j;
        double num=get_d_2D_m(M,i,j,topic_count)*get_d_2D_m(THETA,j,word,total_words_count);
        double punctuation_multiple_num=get_punctuation_multipler_num(doc,w,t_last,t_prime);
        num=num*punctuation_multiple_num;
        for(int x=0;x<topic_count;x++){
          if(i==x){
            double punctuation_multiple_den=get_punctuation_multipler_den(doc,w,i,x);
            den+=(punctuation_multiple_den*get_d_2D_m(M,i,x,topic_count)*get_d_2D_m(THETA,x,word,total_words_count));
          }else{
            den+=get_d_2D_m(M,i,x,topic_count)*get_d_2D_m(THETA,x,word,total_words_count);
          }
        }
        double res=  num/den;
        assert(res>0);
        set_d_3D_m(factor_matrix_log,i,j,w,log(res),topic_count,N_d);
//        cout << "FACTOR[" <<i+1 << "," << j << "," << w <<"]=" << log(res) << "\n" ;
      }
    }
  }
}

  void print_path(int i , int j , int gamma){
    //BASE CASE:
    if(i==0){
      return ;
    }
    else{
      bool isTopic_in_U=j==topic_count-1  ? true : false;
      int k=get_i_3D_m(back_ptr,i,j,gamma,topic_count,GAMMA);
      path[i]=k;
      j=path[i];
      i=i-1;
      gamma=isTopic_in_U ? gamma-1 :gamma ;
      print_path(i,j,gamma);
    }
  }

  void print_filled_table(document* doc){
    cout << "printing table \n" ;
    int word_index=doc->real_word_count;
    for(int i=0;i<word_index;i++)
      for(int j=0;j<topic_count;j++){
        for (int gamma=0;gamma<GAMMA ;gamma ++){
        cout << "F[" <<i+1 << "," << j << "," << gamma << "]= " << get_d_3D_m(F,i,j,gamma,topic_count,GAMMA) << "\n" ;
        }
      }
  }

  bool isBaseCase(int word_index , int j , int gamma){

    bool isTopic_in_U= j==topic_count-1 ? true : false ;
    bool res=false ;
    //  if(word_index==1 ||(gamma>word_index)|| (!isTopic_in_U && gamma ==word_index ) || (isTopic_in_U && word_index >gamma+1)){
    //if(word_index==1 ||(gamma>word_index)|| (!isTopic_in_U && gamma ==word_index ) || (isTopic_in_U && word_index ==gamma+1)){
    if(word_index==1){
      res=true;
    }
    else{
      if(word_index >1 && gamma>word_index){
        res=true;
      }
      else if ( word_index >1 &&!isTopic_in_U && gamma ==word_index ){
        res=true;
      }
      else if (word_index>1 &&isTopic_in_U && gamma==0){
        res=true;
      }
      else{
      }
    }
    return res;
  }

  
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
      //decide how to append or prepend 
      if( path[i]==topic_count-1){ 
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
//    cout << "file_name:" << doc->document_name << " match%= " << v(double)match/(doc->mygamma)*100 <<  " word_count= " << doc->document_name <<" " << " GAMMA= " << doc->mygamma << endl ;
   cout << doc->document_name << "    " << (double)match/(doc->mygamma)*100 << endl; 
   // cout<< "closing after reconstructing test!" << endl ;
   // exit(0);
   original_name.clear();
   file_name.clear();
   html_end_tag.clear();
   html_heading.clear();
   start.clear();
   end.clear();
   html_end_tag.clear();
  
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
   // assert(it != word_id_map.end());
    //  cout << "for word: " <<word << " id was: " << it->second << " fetching" "\n";
//    assert(it->second < total_words_count  && it->second >=0);
//    return it->second;
  return id;
  }

  string get_word_from_id(const int id){
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



  //C(fixed) and calculate statistics over  UP

  void run_test(document* d_curr){
    if(d_curr->mygamma==0){
    cout << "GAMMA=0 "  << d_curr->document_name << endl ;
      document_fee_mem(d_curr);
      return ;
    }
    GAMMA=d_curr->mygamma ;
    int N_d=d_curr->real_word_count ;
    int sz= topic_count*topic_count*(d_curr->real_word_count);
    factor_matrix_log=(double*)malloc(sz*sizeof(double));
    assert(factor_matrix_log!=NULL);
    memset(factor_matrix_log,0,sz*sizeof(double));
    make_factor_matrix_log(d_curr);
    path= (int*)malloc(sizeof(int)*(d_curr->real_word_count));
    assert(path!=NULL);
    memset(path,-1,sizeof(int)*(d_curr->real_word_count));
    initialize_memory_for_DP(d_curr);
    initialize_base_case_for_DP(d_curr) ;
    for(int j=0;j<topic_count ;j++){
      for(int gamma=GAMMA-1;gamma>=0;gamma--){
        double log_val=under_lining(d_curr,N_d-1,j,gamma);
        set_d_3D_m(F,N_d-1,j,gamma,log_val,topic_count,GAMMA);
      }
    }
    print_path(N_d-1,topic_count-1,GAMMA-1);
    reconstruct_doc(d_curr);
    free(factor_matrix_log);
    free(F);
    free(path);
    free(back_ptr);
    document_fee_mem(d_curr);
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

int main(){
  //  print_system_config();
  alpha_upper_N=1;
  alpha_upper_Y=1;
  alpha_comma_N=1;
  alpha_comma_Y=1;
  alpha_period_Y=1;
  alpha_period_N=1;

  init();

  test_document(src1); // creates document  structure !

}
