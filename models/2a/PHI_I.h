//#include "linearmodel.h"
#include "2a.h"
#include "Minuit2/FCNBase.h"

namespace ROOT {
  namespace Minuit2 {
    class PHI_I : public FCNBase {
    public:
        PHI_I(int i1):
          i(i1)
      { } //constructor !
        ~PHI_I(){
        } 
          //destructor

        double operator()(const std::vector<double>& phi_i_vec)const{
          double val ;
          double x0=phi_i_vec[0];
          assert(isnan(x0)==0);
          assert(isinf(x0)==0);
 
          val=-1*evaluate_phi_i(x0,i); // for minima 
          return val;         
        } // operator ends here !

        double Up() const {
          return 1 ;
        }//why > becz manual says so !
      private:
        int i;
          };
  }  // namespace Minuit2
}  // namespace ROOT
