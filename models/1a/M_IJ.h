#include "1a.h"
#include "Minuit2/FCNBase.h"

namespace ROOT {
  namespace Minuit2 {
    class M_IJ : public FCNBase {

      public:
        M_IJ(int i1,int j1,double* sum_den1,uint32_t* freq_den1,uint32_t* freq_num1,double * max_pos1,double* log_F_max_possible1 )
        //  : i(i1),j(j1)
      {
      i=i1;
      j=j1;
      sum_den=sum_den1;
      freq_den=freq_den1;
      freq_num=freq_num1;
      max_pos=max_pos1;
      log_F_max_possible=log_F_max_possible1;
      } //constructor !
        ~M_IJ(){
        }  //destructor
        double operator()(const std::vector<double>& m_i_j_vec)const{
          double  x0 =m_i_j_vec[0];
         assert(x0>0);
         assert(isinf(x0)==0);
         assert(isnan(x0)==0);
         double  val=-1*evaluate_m_ij_ll(x0,i,j,sum_den,freq_den,freq_num,max_pos,log_F_max_possible);
          return val;
        } // operator ends here !

        double Up() const {
          return 1;
        }//why > becz manual says so !
      private:
        int i;
        int j;
        double* sum_den;
        uint32_t* freq_den;
        uint32_t* freq_num;
        double * max_pos;
        double * log_F_max_possible;
    };
  }  // namespace Minuit2
}  // namespace ROOT
