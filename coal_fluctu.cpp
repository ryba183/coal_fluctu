/*
 * calculate coalescence fluctuations and statistics
 * ensamble of nx cells * n_rep repetitions
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

#define Onishi


using namespace std;
using namespace libcloudphxx::lgrngn;

namespace hydrostatic = libcloudphxx::common::hydrostatic;
namespace theta_std = libcloudphxx::common::theta_std;
namespace theta_dry = libcloudphxx::common::theta_dry;
namespace lognormal = libcloudphxx::common::lognormal;

//aerosol bimodal lognormal dist. 
const quantity<si::length, float>
  mean_rd1 = float(15e-6) * si::metres;  // Onishi
//  mean_rd1 = float(9.3e-6) * si::metres;  // WANG 2007 (and Unterstrasser 2017)
const quantity<si::dimensionless, float>
  sdev_rd1 = float(1.4);
const quantity<power_typeof_helper<si::length, static_rational<-3>>::type, float>
  n1_stp = float(142e6) / si::cubic_metres;


//globals
std::array<float, 1201> rad_bins;
int n_cell;
float rho_stp_f;
const int n_rep = 1;
const int sim_time=2500;//500;//2500; // 2500 steps
const int nx = 1e3;
const int ny = 1;
const int nz = 1;
const float dt = 1.;
const float Np = 4e4;
#ifdef Onishi
  const float dx = Np /  (n1_stp * si::cubic_metres); // for Onishi comparison
#else
  const float dx = 1e-6; // for bi-disperse (alfonso) comparison
#endif
const int dev_id=1;
const int sstp_coal = 10;
float cell_vol;

//const int sd_const_multi = 1; const float sd_conc = 0; const bool tail = 0;

const int sd_const_multi = 0; const float sd_conc = 1024; const bool tail = true;



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

  std::cout << "max possible rad (based on mean 3rd wet mom): " << pow(sum/n_cell * rho_stp_f * cell_vol, 1./3.) << std::endl;

  // get spectrum
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
    mean = mean * rho_stp_f / n_cell; // mean number of droplets of radius rad [1/m^3]
    
    // to get mass in bins in [g/cm^3]
    
    res_bins[i]= mean / 1e6 // now its number per cm^3
                     * 3.14 *4. / 3. *rad * rad * rad * 1e3 * 1e3;  // to get mass in grams
    

    // to get mass density function (not through libcloudphxx estimator)
    /*
    res_bins[i] = mean / (rad_bins[i+1] - rad_bins[i]) // number density 
                    * 4./3.*3.14*pow(rad,4)*1e3        // * vol * rad * density
                    * 1e3;                             // to get grams
    */
  }
    std::cout << "res_bins sum (LWC?): " << std::accumulate(res_bins.begin(), res_bins.end(), 0.) << std::endl;
}



int main(){
  std::ofstream of_size_spectr("size_spectr.dat");
  std::ofstream of_max_drop_vol("max_drop_vol.dat");

  std::array<float, nx> init_cloud_mass;
  std::array<float, nx> init_rain_mass;
  std::array<float, nx * n_rep> t10;
  //std::array<std::array<float, nx * n_rep>, sim_time> float tau; // ratio of rain mass to LWC
  auto tau = new float[sim_time+1][nx*n_rep]; // ratio of rain mass to LWC
  //std::array<std::array<float, nx * n_rep>, sim_time> max_rw; // max wet radius in cell
  auto max_rw = new float[sim_time+1][nx*n_rep]; // ratio of rain mass to LWC
  t10.fill(0.);

  std::vector<std::array<float, 1201>> res_bins_pre(n_rep);
//  auto res_bins_pre = new float[n_rep][1201];
  std::vector<std::array<float, 1201>> res_bins_post(n_rep);
//  auto res_bins_post = new float[n_rep][1201];
  std::iota(rad_bins.begin(), rad_bins.end(), 0);
  for (auto &rad_bin : rad_bins)
  {
    rad_bin = rad_bin * 1e-6;// + 10e-6; 
  }

#ifdef Onishi
  std::cout << "Onishi (expvolume) run!" << std::endl;
  std::cout << "Np = " << Np << std::endl;
#else
  std::cout << "Alfonso (bi-disperse) run!" << std::endl;
  std::cout << "dx = " << dx << std::endl;
#endif

  std::cout << "n_rep = " << n_rep 
            << " nx = " << nx
            << " sim_time = " << sim_time
            << " dt = " << dt
            << " sstp_coal = " << sstp_coal
            << " dev_id = " << dev_id
            << " const_multi = " << sd_const_multi
            << " sd_conc = " << sd_conc
            << " tail = " << tail
            << std::endl;

  // repetitions loop
  for(int rep = 0; rep < n_rep; ++rep)
  {
    opts_init_t<float> opts_init;
  
    opts_init.dev_id=dev_id;
    opts_init.dt=dt;
    opts_init.sstp_coal = sstp_coal; 
    opts_init.sstp_cond = 1; 
//    opts_init.kernel = kernel_t::hall_pinsky_1000mb_grav;
    opts_init.kernel = kernel_t::hall;
//    opts_init.kernel = kernel_t::Long;

    opts_init.terminal_velocity = vt_t::beard76;
  //  opts_init.terminal_velocity = vt_t::beard77fast;
    opts_init.dx = dx;
    opts_init.dy = 1;
    opts_init.dz = 1; 
    cell_vol = opts_init.dx;
  
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

  
    opts_init.sd_conc = sd_conc;//int(1024);
    opts_init.sd_conc_large_tail = tail;
    opts_init.sd_const_multi = sd_const_multi;
  //  opts_init.n_sd_max = 20e6 * opts_init.x1 * opts_init.y1 * opts_init.z1 + 1;
    opts_init.n_sd_max = 1e8;// 20e6 * opts_init.x1 * opts_init.y1 * opts_init.z1 + 1;
//  std::cout << "opts_init.n_sd_max: " << opts_init.n_sd_max << std::endl; 
  
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

#ifdef Onishi
    boost::assign::ptr_map_insert<
      exp_dry_radii<float> // value type
    >(  
      opts_init.dry_distros // map
    )(  
      0. // key
    ); 
#else
    opts_init.dry_sizes[0.] = {{17e-6, 20e6}, {21.4e-6, 10e6}};
#endif
  
    std::unique_ptr<particles_proto_t<float>> prtcls(
      factory<float>(
        (backend_t)OpenMP,//CUDA, 
        opts_init
      )
    );
  
    using libcloudphxx::common::earth::rho_stp;
    rho_stp_f = (rho_stp<float>() / si::kilograms * si::cubic_metres);
//    std::cout << "rho stp f = " << rho_stp_f << std::endl;
  
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
    // get max rw
    prtcls->diag_max_rw();
    arr = prtcls->outbuf();
    for(int j=0; j<n_cell; ++j)
    {
      if(arr[j] > glob_max_rw) glob_max_rw = arr[j];
      max_rw[0][j + rep * nx] = arr[j];
    }
  
    for(int i=1; i<=sim_time; ++i)
    {
      two_step(prtcls.get(),th,rv,opts);
  
      // get max rw
      prtcls->diag_max_rw();
      arr = prtcls->outbuf();
      for(int j=0; j<n_cell; ++j)
      {
        if(arr[j] > glob_max_rw) glob_max_rw = arr[j];
        max_rw[i][j + rep * nx] = arr[j];
      }

      // get mean_sd_conc
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
  for(int i=0; i<=sim_time; ++i)
  {
    float mean_max_vol = 0.;
    float std_dev_max_vol = 0.;
    for(int j=0; j < nx * n_rep; ++j)
      mean_max_vol += std::pow(max_rw[i][j], 3); // mass std dev...

    mean_max_vol /= float(nx * n_rep);
    for(int j=0; j < nx * n_rep; ++j)
      std_dev_max_vol += std::pow(std::pow(max_rw[i][j], 3) / mean_max_vol - 1, 2);
    std_dev_max_vol = std::sqrt(std_dev_max_vol / (nx * n_rep-1));
    
    of_max_drop_vol << i * dt << " " << mean_max_vol << " " << std_dev_max_vol << std::endl; // save mean and relative std_dev
  }

  // calc and print out mean t10 and t10 std_dev
  float mean_t10 = 0.;
  int t10_sims = 0;
  for(int j=0; j<n_cell * n_rep; ++j)
  {
    if(t10[j] == 0.)
      std::cerr << "t10 == 0. !! too short simulation" << std::endl;
    else
    {
      ++t10_sims;
      mean_t10 += t10[j];
    }
  }
  mean_t10/=t10_sims;//n_cell*n_rep;
  float std_dev=0.;
  for(int j=0; j<n_cell*n_rep; ++j)
    std_dev += pow(t10[j] - mean_t10, 2);
  std_dev /= (t10_sims-1);
  std_dev = sqrt(std_dev);
  std::cout << "mean(t10%) = " << mean_t10 << std::endl;
  std::cout << "realtive std_dev(t10%) = " << std_dev / mean_t10 << std::endl;

  // calc and print out mean tau10 and tau10 std_dev
  const int mean_t10_idx = mean_t10 / dt + 0.5;
  float mean = 0.;
  if(mean_t10 > 0.)
  {
    for(int j=0; j<n_cell*n_rep; ++j)
    {
      mean += tau[mean_t10_idx][j];
    }
  mean/=n_cell*n_rep;
  std_dev=0.;
  for(int j=0; j<n_cell*n_rep; ++j)
    std_dev += pow(tau[mean_t10_idx][j] - mean, 2);
  std_dev /= (n_cell*n_rep-1);
  std_dev = sqrt(std_dev);
  }
  else 
  {
    mean = NAN;
    std_dev = NAN;
  }
  std::cout << "mean(tau10%) = " << mean << std::endl;
  std::cout << "realtive std_dev(tau10%) = " << std_dev / mean << std::endl;

// output
  for (int i=0; i <rad_bins.size() -1; ++i)
  {
    float rad = (rad_bins[i] + rad_bins[i+1]) / 2.;
    float pre = 0, post = 0;
    for(int j=0; j< n_rep; ++j)
    {
      pre += res_bins_pre[j][i];
      post += res_bins_post[j][i];
    }
    pre /= n_rep;
    post /= n_rep;
    of_size_spectr << rad * 1e6 << " " << pre << " " << post << std::endl; 
  }

//  debug::print(prtcls->impl->n);

//  for(int i=0;i<100;++i)
//  {
//    two_step(prtcls,th,rhod,rv,opts);
//  }

}
