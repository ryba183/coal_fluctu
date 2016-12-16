#include <iostream>
#include <libcloudph++/lgrngn/factory.hpp>
#include <boost/assign/ptr_map_inserter.hpp>
#include <stdio.h>
#include <libcloudph++/common/hydrostatic.hpp>
#include <libcloudph++/common/theta_std.hpp>
#include <libcloudph++/common/theta_dry.hpp>
#include <libcloudph++/common/lognormal.hpp>
#include <libcloudph++/common/unary_function.hpp>
#include <time.h>
#include <libcloudph++/common/earth.hpp>


using namespace std;
using namespace libcloudphxx::lgrngn;

  namespace hydrostatic = libcloudphxx::common::hydrostatic;
  namespace theta_std = libcloudphxx::common::theta_std;
  namespace theta_dry = libcloudphxx::common::theta_dry;
  namespace lognormal = libcloudphxx::common::lognormal;

  //aerosol bimodal lognormal dist. 
  const quantity<si::length, float>
    mean_rd1 = float(17e-6) * si::metres,
    mean_rd2 = float(.15e-6 / 2) * si::metres;
  const quantity<si::dimensionless, float>
    sdev_rd1 = float(1.4),
    sdev_rd2 = float(1.6);
  const quantity<power_typeof_helper<si::length, static_rational<-3>>::type, float>
    n1_stp = float(10e6) / si::cubic_metres,
    n2_stp = float(40e6) / si::cubic_metres;



// lognormal aerosol distribution
template <typename T>
struct log_dry_radii : public libcloudphxx::common::unary_function<T>
{
  T funval(const T lnrd) const
  {   
    return T(( 
        lognormal::n_e(mean_rd1, sdev_rd1, n1_stp, quantity<si::dimensionless, float>(lnrd))
      // +  lognormal::n_e(mean_rd2, sdev_rd2, n2_stp, quantity<si::dimensionless, float>(lnrd)) 
      ) * si::cubic_metres
    );  
  }   

  log_dry_radii *do_clone() const 
  { return new log_dry_radii( *this ); }
};  


void two_step(particles_proto_t<float> *prtcls, 
             arrinfo_t<float> th,
             arrinfo_t<float> rhod,
             arrinfo_t<float> rv,
             opts_t<float> opts)
{
    prtcls->step_sync(opts,th,rv,rhod);
    cout << prtcls->step_async(opts) << endl;
}


int main(){
  opts_init_t<float> opts_init;

  int sim_time=500;//2500; // 2500 steps

  opts_init.dt=sim_time;
  opts_init.sstp_coal = sim_time; 
  opts_init.kernel = kernel_t::hall;
  opts_init.terminal_velocity = vt_t::beard77fast;
  opts_init.dx = 1e-2;
  opts_init.dy = 1e-2;
  opts_init.dz = 1e-2; 

//  opts_init.dx = 1;
//  opts_init.dy = 1;
//  opts_init.dz = 1; 

  opts_init.nx = 1; 
  opts_init.ny = 1; 
  opts_init.nz = 1; 
  opts_init.x1 = opts_init.nx * opts_init.dx;
  opts_init.y1 = opts_init.ny * opts_init.dy;
  opts_init.z1 = opts_init.nz * opts_init.dz;
  opts_init.rng_seed = time(NULL);

//  opts_init.sd_conc = 100;
  opts_init.sd_const_multi = 1;
  opts_init.n_sd_max = 60e6;

  boost::assign::ptr_map_insert<
    log_dry_radii<float> // value type
  >(  
    opts_init.dry_distros // map
  )(  
    0. // key
  ); 

  particles_proto_t<float> *prtcls;
     prtcls = factory<float>(
        (backend_t)OpenMP, 
        opts_init
      );

  using libcloudphxx::common::earth::rho_stp;
  const float rho_stp_f = (rho_stp<float>() / si::kilograms * si::cubic_metres);
  std::cout << "rho stp f = " << rho_stp_f << std::endl;

  float pth[] = {300.};
  float prhod[] = {rho_stp_f};
  float prv[] = {.01};
  long int strides[] = {sizeof(float)};

  arrinfo_t<float> th(pth, strides);
  arrinfo_t<float> rhod(prhod, strides);
  arrinfo_t<float> rv(prv, strides);

  prtcls->init(th,rv,rhod);

  opts_t<float> opts;
  opts.adve = 0;
  opts.sedi = 0;
  opts.cond = 0;
  opts.coal = 1;
  opts.rcyc = 0;

  prtcls->diag_sd_conc();
  std::cout << "sd conc: " << prtcls->outbuf()[0] << std::endl;

  prtcls->diag_all();
  prtcls->diag_dry_mom(0);
  auto drop_no = prtcls->outbuf()[0];
  std::cout << "number of droplets: " << drop_no << std::endl;

  prtcls->diag_all();
  prtcls->diag_wet_mom(1);
  std::cout << "mean wet radius: " << prtcls->outbuf()[0] / drop_no << std::endl;

  prtcls->diag_all();
  prtcls->diag_dry_mom(1);
  std::cout << "mean dry radius: " << prtcls->outbuf()[0] / drop_no << std::endl;

  prtcls->step_sync(opts,th,rv,rhod);
  cout << prtcls->step_async(opts) << endl;

  std::cout << "po symulacji:" << std::endl;

  prtcls->diag_sd_conc();
  std::cout << "sd conc: " << prtcls->outbuf()[0] << std::endl;

  prtcls->diag_all();
  prtcls->diag_dry_mom(0);
  drop_no = prtcls->outbuf()[0];
  std::cout << "number of droplets: " << drop_no << std::endl;

  prtcls->diag_all();
  prtcls->diag_wet_mom(1);
  std::cout << "mean wet radius: " << prtcls->outbuf()[0] / drop_no << std::endl;

  prtcls->diag_all();
  prtcls->diag_dry_mom(1);
  std::cout << "mean dry radius: " << prtcls->outbuf()[0] / drop_no << std::endl;

//  debug::print(prtcls->impl->n);

//  for(int i=0;i<100;++i)
//  {
//    two_step(prtcls,th,rhod,rv,opts);
//  }

}
