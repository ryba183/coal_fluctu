/*
 * calculate coalescence fluctuations and statistics
 * ensamble of nx same cells * n_rep repetitions
*/

#include <iostream>
#include <fstream>
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
  mean_rd1 = float(15e-6) * si::metres;
//  mean_rd2 = float(.15e-6 / 2) * si::metres;
const quantity<si::dimensionless, float>
  sdev_rd1 = float(1.4);
//  sdev_rd2 = float(1.6);
const quantity<power_typeof_helper<si::length, static_rational<-3>>::type, float>
  n1_stp = float(142e6) / si::cubic_metres;
//  n2_stp = float(40e6) / si::cubic_metres;


//globals
std::array<float, 1201> rad_bins;
int n_cell;
float rho_stp_f;
const int n_rep = 5;
const int sim_time=1;//500;//2500; // 2500 steps
const int nx = 1;
const int ny = 1;
const int nz = 1;
const float dt = 1.;


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

// aerosol distribution exponential in droplet volume as a function of ln(r)
template <typename T>
struct exp_dry_radii : public libcloudphxx::common::unary_function<T>
{
  T funval(const T lnrd) const
  {   
    T r = exp(lnrd);
    return (n1_stp * si::cubic_metres) * 3. * pow(r,3) / pow(mean_rd1 / si::metres, 3) * exp( - pow(r/(mean_rd1 / si::metres), 3));
  }   

  exp_dry_radii *do_clone() const 
  { return new exp_dry_radii( *this ); }
};  

void two_step(particles_proto_t<float> *prtcls, 
             arrinfo_t<float> th,
     //        arrinfo_t<float> rhod,
             arrinfo_t<float> rv,
             opts_t<float> opts)
{
    prtcls->step_sync(opts,th,rv);//,rhod);
    prtcls->step_async(opts);
}


void diag(particles_proto_t<float> *prtcls, std::array<float, 1201> &res_bins)
{
  prtcls->diag_sd_conc();
  std::cout << "sd conc: " << prtcls->outbuf()[0] << std::endl;

/*
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

*/
  prtcls->diag_all();
  prtcls->diag_wet_mom(3);
  float sum = 0;
  auto out = prtcls->outbuf();
  for(int c=0; c < n_cell; ++c)
    sum += out[c];
  std::cout << "3rd wet mom mean: " << sum / n_cell << std::endl;

  // get mass density function
  for (int i=0; i <rad_bins.size() -1; ++i)
  {
//    prtcls->diag_all();
  //  prtcls->diag_wet_size_spectr( rad_bins[i], 0.62 ); //sigma0 = 0.62 like in Shima (2009)
    prtcls->diag_wet_rng(rad_bins[i], rad_bins[i+1]);
    prtcls->diag_wet_mom(0);
    float rad = (rad_bins[i] + rad_bins[i+1]) / 2.;
    float mean = 0;
    auto buf = prtcls->outbuf();
    for(int c=0; c < n_cell; ++c)
    {
//      std::cout << buf[c] << " ";
      mean += buf[c];
    }
    mean /= n_cell;
    
    res_bins[i]= mean * rho_stp_f / 1e6 // now its number per cm^3
                     * 3.14 *4. / 3. *rad * rad * rad * 1e3 * 1e3;  // to get mass in grams
                     // / (rad_bins[i+1] - rad_bins[i]); // to get density function
//    of_size_spectr << rad * 1e6 << " " << res_bins_pre[i] << " " << res_bins_post[i] << std::endl; 
  }
    std::cout << "res_bins sum (LWC?): " << std::accumulate(res_bins.begin(), res_bins.end(), 0.) << std::endl;
}



int main(){
  std::ofstream of_size_spectr("size_spectr.dat");
  std::ofstream of_max_rw_std_dev("max_rw_std_dev.dat");

  std::array<float, nx> init_cloud_mass;
  std::array<float, nx> init_rain_mass;
  std::array<float, nx * n_rep> t10;
  std::array<std::array<float, nx * n_rep>, sim_time> tau; // ratio of rain mass to LWC
  std::array<std::array<float, nx * n_rep>, sim_time> max_rw; // max wet radius in cell
  t10.fill(0.);

  std::array<std::array<float, 1201>, n_rep> res_bins_pre;
  std::array<std::array<float, 1201>, n_rep> res_bins_post;
  std::iota(rad_bins.begin(), rad_bins.end(), 0);
  for (auto &rad_bin : rad_bins)
  {
    rad_bin = rad_bin * 1e-6 + 10e-6; 
  }
  // repetitions loop
  for(int rep = 0; rep < n_rep; ++rep)
  {
    opts_init_t<float> opts_init;
  
    opts_init.dt=dt;
    opts_init.sstp_coal = 1; 
    opts_init.sstp_cond = 1; 
    opts_init.kernel = kernel_t::hall_pinsky_1000mb_grav;
    opts_init.terminal_velocity = vt_t::beard76;
    opts_init.dx = 1.06e4 /  (n1_stp * si::cubic_metres);
    opts_init.dy = 1;
    opts_init.dz = 1; 
  
    opts_init.sedi_switch=0;
    opts_init.src_switch=0;
    opts_init.chem_switch=0;
  
    opts_init.nx = nx; 
    opts_init.ny = ny; 
    opts_init.nz = 1; 
    opts_init.x1 = opts_init.nx * opts_init.dx;
    opts_init.y1 = opts_init.ny * opts_init.dy;
    opts_init.z1 = opts_init.nz * opts_init.dz;
    opts_init.rng_seed = time(NULL);
  
    n_cell = opts_init.nx * opts_init.ny * opts_init.nz;
  
  //  opts_init.sd_conc = int(128);
    opts_init.sd_const_multi = 1;
  //  opts_init.n_sd_max = 20e6 * opts_init.x1 * opts_init.y1 * opts_init.z1 + 1;
    opts_init.n_sd_max = 3e7;// 20e6 * opts_init.x1 * opts_init.y1 * opts_init.z1 + 1;
  std::cout << "opts_init.n_sd_max: " << opts_init.n_sd_max << std::endl; 
  
  //  for (auto rad_bin : rad_bins)
    //  std::cout << rad_bin;
  /*
    boost::assign::ptr_map_insert<
      log_dry_radii<float> // value type
    >(  
      opts_init.dry_distros // map
    )(  
      0. // key
    ); 
  */
    boost::assign::ptr_map_insert<
      exp_dry_radii<float> // value type
    >(  
      opts_init.dry_distros // map
    )(  
      0. // key
    ); 
  //  opts_init.dry_sizes[0.] = {{17e-6, 20e6}, {21.4e-6, 10e6}};
  
    std::unique_ptr<particles_proto_t<float>> prtcls(
      factory<float>(
        (backend_t)CUDA, 
        opts_init
      )
    );
  
    using libcloudphxx::common::earth::rho_stp;
    rho_stp_f = (rho_stp<float>() / si::kilograms * si::cubic_metres);
    std::cout << "rho stp f = " << rho_stp_f << std::endl;
  
    std::array<float, 1> pth;
    std::array<float, 1> prhod;
    std::array<float, 1> prv;
  
    pth.fill(300.);
    prhod.fill(rho_stp_f);
    prv.fill(.01);
  
    long int strides[] = {/*sizeof(float)*/ 0, 0, 0};
  
    arrinfo_t<float> th(pth.data(), strides);
    arrinfo_t<float> rhod(prhod.data(), strides);
    arrinfo_t<float> rv(prv.data(), strides);
  
    prtcls->init(th,rv,rhod);
  
    opts_t<float> opts;
    opts.adve = 0;
    opts.sedi = 0;
    opts.cond = 0;
    opts.coal = 1;
    opts.rcyc = 0;
  
    std::fill(res_bins_pre[rep].begin(), res_bins_pre[rep].end(), 0.);
    std::fill(res_bins_post[rep].begin(), res_bins_post[rep].end(), 0.);
    diag(prtcls.get(), res_bins_pre[rep]);
  
  //  prtcls->step_sync(opts,th,rv);//,rhod);
  //  cout << prtcls->step_async(opts) << endl;
  
    prtcls->diag_wet_rng(0, 40e-6); // cloud water (like in Onishi)
    prtcls->diag_wet_mom(3);
    auto arr = prtcls->outbuf();
    for(int j=0; j<n_cell; ++j)
    {
      init_cloud_mass[j] = arr[j];
    }
  
    prtcls->diag_wet_rng(40e-6, 1); // rain water (like in Onishi)
    prtcls->diag_wet_mom(3);
    arr = prtcls->outbuf();
    for(int j=0; j<n_cell; ++j)
    {
      init_rain_mass[j] = arr[j];
    }
  
  
    float glob_max_rw = 0.;
    for(int i=0; i<sim_time; ++i)
    {
      two_step(prtcls.get(),th,rv,opts);
  
      // get std dev of max rw
      prtcls->diag_max_rw();
      arr = prtcls->outbuf();
      for(int j=0; j<n_cell; ++j)
      {
        if(arr[j] > glob_max_rw) glob_max_rw = arr[j];
        max_rw[i][j + rep * nx] = arr[j];
      }
      prtcls->diag_sd_conc();
      arr = prtcls->outbuf();
      float mean_sd_conc = 0;
      for(int j=0; j<n_cell; ++j)
      {
        mean_sd_conc += arr[j]; 
      }
      mean_sd_conc /= float(n_cell);
  
  printf("\rrep no: %3d progress: %3d%%: rw_max %lf mean_sd_conc %lf", rep, int(float(i) / sim_time * 100), glob_max_rw, mean_sd_conc);
  std::cout << std::flush;
  
  
      // get t10 (time to conver 10% of cloud water into rain water)
      prtcls->diag_wet_rng(40e-6, 1); // rain water (like in Onishi)
      prtcls->diag_wet_mom(3);
      arr = prtcls->outbuf();
      for(int j=0; j<n_cell; ++j)
      {
        if(t10[j + rep * nx] == 0. && arr[j] >= init_cloud_mass[j] * .1)
          t10[j + rep * nx] = i * opts_init.dt;
        tau[i][j + rep * nx] = arr[j] / init_cloud_mass[j];
      }
    }
  
    std::cout << std::endl << "po symulacji, max_rw: " << glob_max_rw << std::endl;
  
    diag(prtcls.get(), res_bins_post[rep]);
    std::cout << std::endl;
//    auto raw_ptr = prtcls.release();
//    delete raw_ptr;
  }

  // calc and print max rw (mass) std dev
  for(int i=0; i<sim_time; ++i)
  {
    float mean_max_rw = 0.;
    float std_dev_max_rw = 0.;
    for(int j=0; j < nx * n_rep; ++j)
      mean_max_rw += std::pow(max_rw[i][j], 3); // mass std dev...

    mean_max_rw /= float(nx * n_rep);
    for(int j=0; j < nx * n_rep; ++j)
      std_dev_max_rw += std::pow(std::pow(max_rw[i][j], 3) / mean_max_rw - 1, 2);
    std_dev_max_rw = std::sqrt(std_dev_max_rw / (nx * n_rep-1));
    
    of_max_rw_std_dev << i * dt << " " << std_dev_max_rw << std::endl;
  }

  /*
  // calc and print out mean t10 and t10 std_dev
  float mean_t10 = 0.;
  for(int j=0; j<n_cell; ++j)
  {
    if(t10[j] == 0.)
      std::cerr << "t10 == 0. !! too short simulation" << std::endl;
    mean_t10 += t10[j];
  }
  mean_t10/=n_cell;
  float std_dev=0.;
  for(int j=0; j<n_cell; ++j)
    std_dev += pow(t10[j] - mean_t10, 2);
  std_dev /= (n_cell-1);
  std_dev = sqrt(std_dev);
  std::cout << "mean(t10%) = " << mean_t10 << std::endl;
  std::cout << "realtive std_dev(t10%) = " << std_dev / mean_t10 << std::endl;

  // calc and print out mean tau10 and tau10 std_dev
  const int mean_t10_idx = mean_t10 / opts_init.dt + 0.5;
  float mean = 0.;
  for(int j=0; j<n_cell; ++j)
  {
    mean += tau[mean_t10_idx][j];
  }
  mean/=n_cell;
  std_dev=0.;
  for(int j=0; j<n_cell; ++j)
    std_dev += pow(tau[mean_t10_idx][j] - mean, 2);
  std_dev /= (n_cell-1);
  std_dev = sqrt(std_dev);
  std::cout << "mean(tau10%) = " << mean << std::endl;
  std::cout << "realtive std_dev(tau10%) = " << std_dev / mean << std::endl;

// output
  for (int i=0; i <rad_bins.size() -1; ++i)
  {
    float rad = (rad_bins[i] + rad_bins[i+1]) / 2.;
    of_size_spectr << rad * 1e6 << " " << res_bins_pre[i] << " " << res_bins_post[i] << std::endl; 
  }

*/
//  debug::print(prtcls->impl->n);

//  for(int i=0;i<100;++i)
//  {
//    two_step(prtcls,th,rhod,rv,opts);
//  }

}
