//#include "linearmodel.h"
#include "1b.h"
#include "Minuit2/FCNBase.h"

namespace ROOT {
  namespace Minuit2 {
    class THETA_IW : public FCNBase {
      public:
        //constructor !
        THETA_IW(int i1,int w1,double* den_single_term_header,double* sum_den,uint32_t*freq_den,uint32_t freq_single_num,uint32_t* freq_num,uint32_t freq_single_den_w1,double * possible_maxima1,double *val_maxima1)
          //: i(i1),w(w1)
        {
          i=i1;
          w=w1;
          denominator=den_single_term_header;
          sum_den1=sum_den;
          freq_den1=freq_den;
          freq_single_num1=freq_single_num;
          freq_num1=freq_num;
          freq_single_den_w=freq_single_den_w1;
          possible_maxima=possible_maxima1;
          val_maxima=val_maxima1;
        }         
        //destructor
        ~THETA_IW(){  }  

        double operator()(const std::vector<double>& theta_i_w_vec)const{
          double val ;
          double  x0 =theta_i_w_vec[0];
          assert(x0>0);
          assert(isinf(x0)==0);
          assert(isnan(x0)==0);
          val=-1*evaluate_theta_iw_ll(x0,i,w,*denominator,sum_den1,freq_den1,freq_single_num1,freq_num1,freq_single_den_w,possible_maxima,val_maxima);
          return val; //becuase we want to  maximise but these guys minimise !
        } // operator ends here !

        double Up() const {
          return  1;
        }//why > becz manual says so !

      private:
        int i;
        int w;
        double*  denominator ;
        double *sum_den1;
        uint32_t * freq_den1;
        uint32_t freq_single_num1;
       uint32_t*  freq_num1 ;
       uint32_t freq_single_den_w;
        double*possible_maxima;
        double* val_maxima;
    };
  }  // namespace Minuit2
}  // namespace ROOT
