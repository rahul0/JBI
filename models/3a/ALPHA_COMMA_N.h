#include "3a.h"
#include "Minuit2/FCNBase.h"

namespace ROOT {
  namespace Minuit2 {
    class ALPHA_COMMA_N : public FCNBase {
      public:
        //constructor !
        ALPHA_COMMA_N(const int &num_count1, double* sum_den1,const map<int , vector<int> > & ref1, double * possible_maxima1,double *val_maxima1)
        {
          num_count=num_count1;
          sum_den=sum_den1;
          ref=ref1;
          possible_maxima=possible_maxima1;
          val_maxima=val_maxima1;
        }         
        //destructor
        ~ALPHA_COMMA_N(){  }  

        double operator()(const std::vector<double>& alpha_comma_N_vec)const{
          double val ;
          double  x0 =alpha_comma_N_vec[0];
          assert(x0>0);
          assert(isinf(x0)==0);
          assert(isnan(x0)==0);

          val=-1*evaluate_alpha_comma_N_log_posterior(x0,num_count,sum_den,ref,possible_maxima,val_maxima);
          return val; //becuase we want to  maximise but these guys minimise !
        } // operator ends here !

        double Up() const {
          return  1;
        }//why > becz manual says so !

      private:
        int num_count;
        double *sum_den;
        map<int,vector<int> > ref;   
        double*possible_maxima;
        double* val_maxima;
 
    };
  }  // namespace Minuit2
} 
