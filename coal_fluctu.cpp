/*
 * calculate coalescence fluctuations and statistics
 * ensamble of nx cells * n_rep repetitions
 * statistics of max_rw can be calcul;ated using larger cells, conataining more than one cimulation cell
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
#include <numeric>


 #define Onishi
 #define cutoff 40e-6
// #define HallDavis

#define HIST_BINS 5001
#define BACKEND multi_CUDA
#define N_SD_MAX 1e8 //1e8
#define NXNYNZ 720 // number of cells in each direction
#define SEDI 1
#define RCYC 0
#define N_REP 1e0
#define SIMTIME 4000 // number of steps 
#define NP 1e0 // init number of droplets per cell
#define DT 0.1 // [s]
#define DISS_RATE 1 // [cm^2 / s^3]

using namespace std;
using namespace libcloudphxx::lgrngn;

using real_t = float;

namespace hydrostatic = libcloudphxx::common::hydrostatic;
namespace theta_std = libcloudphxx::common::theta_std;
namespace theta_dry = libcloudphxx::common::theta_dry;
namespace lognormal = libcloudphxx::common::lognormal;

const quantity<si::length, real_t>
  mean_rd1 = real_t(15e-6) * si::metres;  // Onishi
//  mean_rd1 = real_t(0.02e-6) * si::metres;  // api_lgrngn
//  mean_rd1 = real_t(9.3e-6) * si::metres;  // WANG 2007 (and Unterstrasser 2017)
//  mean_rd1 = real_t(10.177e-6) * si::metres;  // Shima small
const quantity<si::dimensionless, real_t>
  sdev_rd1 = real_t(1.4);
const quantity<power_typeof_helper<si::length, static_rational<-3>>::type, real_t>
  n1_stp = real_t(142e6) / si::cubic_metres; // Onishi
//  n1_stp = real_t(297e6) / si::cubic_metres; // WANG 2007 (and Unter 2017)
//  n1_stp = real_t(60e6) / si::cubic_metres; // api_lgrngn
//  n1_stp = real_t(226.49e6) / si::cubic_metres; // Shima small


//globals
std::array<real_t, HIST_BINS> rad_bins;
real_t rho_stp_f;
const int n_rep = N_REP; // number of repetitions of simulation
const int sim_time=SIMTIME; //2500;//500;//2500; // 2500 steps
const int nx = NXNYNZ; // total number of collision cells
const int ny = NXNYNZ;
const int nz = NXNYNZ;
constexpr int n_cell = NXNYNZ * NXNYNZ * NXNYNZ;

const real_t dt = DT;
const real_t Np = NP; // number of droplets per simulation (collision cell)
const real_t Np_in_avg_r_max_cell = Np; // number of droplets per large cells in which we look for r_max
//#ifdef Onishi
  const int n_cells_per_avg_r_max_cell = Np_in_avg_r_max_cell / Np;
  const real_t cell_vol = Np /  (n1_stp * si::cubic_metres); // for Onishi comparison
  const real_t dx = pow(cell_vol, real_t(1./3.));
  const real_t dy = pow(cell_vol, real_t(1./3.));
  const real_t dz = pow(cell_vol, real_t(1./3.));
//#else
//  const int n_cells_per_avg_r_max_cell = 1; // r_max in each small cell separately
//  const real_t dx = 1e-6; // for bi-disperse (alfonso) comparison
//  const real_t dx = 1e6; // for Shima comparison
//#endif
const int n_large_cells = (nx * ny * nz) / n_cells_per_avg_r_max_cell;
const int sstp_coal = 1;


 const int sd_const_multi = 1; const real_t sd_conc = 0; const bool tail = 0;

//  const int sd_const_multi = 0; const real_t sd_conc = 1e3; const bool tail = 1;

// lognormal aerosol distribution
template <typename T>
struct log_dry_radii : public libcloudphxx::common::unary_function<T>
{
  T funval(const T lnrd) const
  {   
    return T(( 
        lognormal::n_e(mean_rd1, sdev_rd1, n1_stp, quantity<si::dimensionless, real_t>(lnrd))
      // +  lognormal::n_e(mean_rd2, sdev_rd2, n2_stp, quantity<si::dimensionless, real_t>(lnrd)) 
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
#ifdef cutoff
    if(r>= cutoff) return 0.; else 
#endif
return (n1_stp * si::cubic_metres) * 3. * pow(r,3) / pow(mean_rd1 / si::metres, 3) * exp( - pow(r/(mean_rd1 / si::metres), 3));
  }   

  exp_dry_radii *do_clone() const 
  { return new exp_dry_radii( *this ); }
};  

void two_step(particles_proto_t<real_t> *prtcls, 
             arrinfo_t<real_t> th,
     //        arrinfo_t<real_t> rhod,
             arrinfo_t<real_t> rv,
             arrinfo_t<real_t> diss_rate,
             opts_t<real_t> opts)
{
    prtcls->step_sync(opts,th,rv,arrinfo_t<real_t>(),arrinfo_t<real_t>(),arrinfo_t<real_t>(),arrinfo_t<real_t>(), diss_rate);
    prtcls->step_async(opts);
}


void diag(particles_proto_t<real_t> *prtcls, std::array<real_t, HIST_BINS> &res_bins, std::array<real_t, HIST_BINS> &res_stddev_bins)
{
  prtcls->diag_all();
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
  real_t sum = 0;
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
    real_t rad = (rad_bins[i] + rad_bins[i+1]) / 2.;
    real_t mean = 0;
    auto buf = prtcls->outbuf();
    for(int c=0; c < n_cell; ++c)
    {
//      std::cout << buf[c] << " ";
      mean += buf[c];
    }
    mean /= n_cell;

    real_t std_dev = 0;
    for(int c=0; c < n_cell; ++c)
    {
      std_dev += pow(buf[c] - mean, 2.);
    }
    std_dev = sqrt(std_dev / n_cell);

    mean = mean * rho_stp_f; // mean number of droplets of radius rad [1/m^3]
    std_dev *= rho_stp_f;
    
    // to get mass in bins in [g/cm^3]
    
    res_bins[i]= mean / 1e6 // now its number per cm^3
                     * 3.14 *4. / 3. *rad * rad * rad * 1e3 * 1e3;  // to get mass in grams
    res_stddev_bins[i]= std_dev / 1e6 // now its number per cm^3
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
  std::cerr << "main start" << std::endl;

  // sanity check
#ifdef Onishi
  if(n_cells_per_avg_r_max_cell * Np != Np_in_avg_r_max_cell)
    throw std::runtime_error("Np_in_avg_r_max_cell nie jest wilokrotnoscia Np");
  if(n_large_cells * n_cells_per_avg_r_max_cell != n_cell)
    throw std::runtime_error("n_cell nie jest wilokrotnoscia n_cells_per_avg_r_max_cell");
#endif

#ifdef cutoff
  std::cout << "init distr cutoff at " << cutoff << " microns!" << std::endl;
#endif

#ifdef Onishi
  std::cout << "Onishi (expvolume) run!" << std::endl;
  std::cout << "Np = " << Np << std::endl;
  std::cout << "Np per avg cell = " << Np_in_avg_r_max_cell << std::endl;
#else
  std::cout << "Alfonso (bi-disperse) run!" << std::endl;
#endif
  std::cout << "dx = " << dx * 1e2  << "cm (cell vol = " << cell_vol * 1e6 << " cm^3)"<< std::endl;
  std::cout << "x1 = " << dx * nx * 1e2  << "cm (domain vol = "<< dx * nx * dy * ny * dz * nz  << " m^3)" << std::endl;

  std::cout << "n_rep = " << n_rep 
            << " n_large_cells = " << n_large_cells
            << " n_cell = " << n_cell
            << " sim_time = " << sim_time
            << " dt = " << dt
            << " sstp_coal = " << sstp_coal
            << " const_multi = " << sd_const_multi
            << " sd_conc = " << sd_conc
            << " tail = " << tail
            << " mean_rd1 = " << mean_rd1
            << " n1_stp = " << n1_stp
            << " sedi = " << SEDI
            << " rcyc = " << RCYC
            << " backend = " << BACKEND
            << std::endl;

  std::cout << std::flush;

  std::ofstream of_size_spectr("size_spectr.dat");
  std::ofstream of_max_drop_vol("max_drop_vol.dat");

  std::vector<real_t> init_cloud_mass(n_cell);
  std::vector<real_t> init_rain_mass(n_cell);
  real_t init_tot_cloud_mass;
  real_t init_tot_rain_mass;

//  std::vector<real_t> t10(n_cell * n_rep, 0);
  std::vector<real_t> t10_tot(n_rep, 0);
//  std::vector<real_t> t_max_40(n_cell * n_rep, 0);// = new real_t[n_cell * n_rep];
//  std::vector<std::vector<real_t>> tau           (sim_time+1, std::vector<real_t>(n_cell*n_rep)); // ratio of rain mass to LWC
//  std::vector<std::vector<real_t>> nrain         (sim_time+1, std::vector<real_t>(n_cell*n_rep)); // number of rain drops
//  std::vector<std::vector<real_t>> max_rw        (sim_time+1, std::vector<real_t>(n_rep * n_large_cells)); // max rw per large (averaging) cell
//  std::vector<std::vector<real_t>> max_rw_small  (sim_time+1, std::vector<real_t>(n_rep * n_cell)); // max rw^3 per small cells (to compare with Alfonso)

  std::vector<std::array<real_t, HIST_BINS>> res_bins_pre(n_rep);
  std::vector<std::array<real_t, HIST_BINS>> res_stddev_bins_pre(n_rep);
//  auto res_bins_pre = new real_t[n_rep][HIST_BINS];
  std::vector<std::array<real_t, HIST_BINS>> res_bins_post(n_rep);
  std::vector<std::array<real_t, HIST_BINS>> res_stddev_bins_post(n_rep);
//  auto res_bins_post = new real_t[n_rep][HIST_BINS];
  std::iota(rad_bins.begin(), rad_bins.end(), 0);
  for (auto &rad_bin : rad_bins)
  {
    rad_bin = rad_bin * 1e-6;// + 10e-6; 
  }

  // repetitions loop
  for(int rep = 0; rep < n_rep; ++rep)
  {
    opts_init_t<real_t> opts_init;
  
    opts_init.dt=dt;
    opts_init.sstp_coal = sstp_coal; 
    opts_init.sstp_cond = 1; 
//    opts_init.kernel = kernel_t::hall_pinsky_1000mb_grav;
#ifdef HallDavis
    opts_init.kernel = kernel_t::hall_davis_no_waals;
#else
    opts_init.kernel = kernel_t::hall;
#endif
//    opts_init.kernel = kernel_t::Long;
  //  opts_init.kernel = kernel_t::geometric;

    opts_init.terminal_velocity = vt_t::beard76;
  //  opts_init.terminal_velocity = vt_t::beard77fast;
    opts_init.dx = dx;
    opts_init.dy = dy;
    opts_init.dz = dz; 
  //  cell_vol = opts_init.dx * opts_init.dy * opts_init.dz;
  
    opts_init.sedi_switch=1;
    opts_init.src_switch=0;
    opts_init.chem_switch=0;
    opts_init.turb_adve_switch=1;

    opts_init.periodic_topbot_walls = 1;
  
    opts_init.nx = nx; 
    opts_init.ny = ny; 
    opts_init.nz = nz; 
    opts_init.x0 = 0;
    opts_init.y0 = 0;
    opts_init.z0 = 0;
    opts_init.x1 = opts_init.nx * opts_init.dx;
    opts_init.y1 = opts_init.ny * opts_init.dy;
    opts_init.z1 = opts_init.nz * opts_init.dz;
    opts_init.rng_seed = time(NULL);
  
  
    opts_init.sd_conc = sd_conc;//int(1024);
    opts_init.sd_conc_large_tail = tail;
    opts_init.sd_const_multi = sd_const_multi;
  //  opts_init.n_sd_max = 20e6 * opts_init.x1 * opts_init.y1 * opts_init.z1 + 1;
    opts_init.n_sd_max = N_SD_MAX;// 20e6 * opts_init.x1 * opts_init.y1 * opts_init.z1 + 1;
//  std::cout << "opts_init.n_sd_max: " << opts_init.n_sd_max << std::endl; 
  
  //  for (auto rad_bin : rad_bins)
    //  std::cout << rad_bin;
//    boost::assign::ptr_map_insert<
//      log_dry_radii<real_t> // value type
//    >(  
//      opts_init.dry_distros // map
//    )(  
//      0. // key
//    ); 

#ifdef Onishi
    opts_init.dry_distros.emplace(
      0, //0 // key (kappa)
      std::make_shared<exp_dry_radii<real_t>> () // value
    );

#else
    opts_init.dry_sizes.emplace(
      0, //0 // key (kappa)
      std::map<real_t, std::pair<real_t, int> > {
        {17e-6  , {20e6, 20e6 * cell_vol}}, // radius, STP concentration, number of SD
        {21.4e-6, {10e6, 10e6 * cell_vol}}, // radius, STP concentration, number of SD
      }
    );
#endif

    opts_init.SGS_mix_len = std::vector<real_t>(nz, opts_init.z1); // z1 because the whole domain is like a LES cell in which flow is not resolved
  
    std::unique_ptr<particles_proto_t<real_t>> prtcls(
      factory<real_t>(
        (backend_t)BACKEND, 
        opts_init
      )
    );
  
    using libcloudphxx::common::earth::rho_stp;
    rho_stp_f = (rho_stp<real_t>() / si::kilograms * si::cubic_metres);
  
    std::vector<real_t> pth(n_cell, 300);
    std::vector<real_t> prhod(n_cell, rho_stp_f);
    std::vector<real_t> prv(n_cell, .01);
    std::vector<real_t> pdiss_rate(n_cell, DISS_RATE * 1e-4); // 1e-4 to turn cm^2/s^3 to m^2/s^3
  
    //long int strides[] = {sizeof(real_t) * NXNYNZ * NXNYNZ, sizeof(real_t) * NXNYNZ, sizeof(real_t)};
    long int strides[] = {1 * NXNYNZ * NXNYNZ, 1 * NXNYNZ, 1};
  
    arrinfo_t<real_t> th(pth.data(), strides);
    arrinfo_t<real_t> rhod(prhod.data(), strides);
    arrinfo_t<real_t> rv(prv.data(), strides);
    arrinfo_t<real_t> diss_rate(pdiss_rate.data(), strides);

    prtcls->init(th,rv,rhod);
  
    opts_t<real_t> opts;
    opts.adve = 0;
    opts.sedi = SEDI;
    opts.cond = 0;
    opts.coal = 1;
    opts.rcyc = RCYC;
    opts.turb_adve = 1;
  
    std::fill(res_bins_pre[rep].begin(), res_bins_pre[rep].end(), 0.);
    std::fill(res_bins_post[rep].begin(), res_bins_post[rep].end(), 0.);
    std::fill(res_stddev_bins_pre[rep].begin(), res_stddev_bins_pre[rep].end(), 0.);
    std::fill(res_stddev_bins_post[rep].begin(), res_stddev_bins_post[rep].end(), 0.);
    diag(prtcls.get(), res_bins_pre[rep], res_stddev_bins_pre[rep]);
  
  //  prtcls->step_sync(opts,th,rv);//,rhod);
  //  cout << prtcls->step_async(opts) << endl;

    init_tot_cloud_mass = 0;
    init_tot_rain_mass = 0;
  
    prtcls->diag_wet_rng(0, 40e-6); // cloud water (like in Onishi)
    prtcls->diag_wet_mom(3);
    auto arr = prtcls->outbuf();
    for(int j=0; j<n_cell; ++j)
    {
      init_cloud_mass[j] = arr[j];
      init_tot_cloud_mass += arr[j];
    }
  
    prtcls->diag_wet_rng(40e-6, 1); // rain water (like in Onishi)
    prtcls->diag_wet_mom(3);
    arr = prtcls->outbuf();
    for(int j=0; j<n_cell; ++j)
    {
      init_rain_mass[j] = arr[j];
      init_tot_rain_mass += arr[j];
    }
  
    real_t rep_max_rw = 0.;
    // get max rw
    prtcls->diag_max_rw();
    arr = prtcls->outbuf();
    int large_cell_idx = -1;
//    for(int j=0; j<n_cell; ++j)
//    {
////      real_t large_cell_max_rw;
//      if(j % n_cells_per_avg_r_max_cell == 0)
//      {
//  //      large_cell_max_rw = 0.;
//        large_cell_idx++;
//        max_rw[0][large_cell_idx + rep * n_large_cells] = 0.;//large_cell_max_rw;
//      }
//      if(arr[j] > rep_max_rw) rep_max_rw = arr[j];
//      if(arr[j] > max_rw[0][large_cell_idx + rep * n_large_cells])
//        max_rw[0][large_cell_idx + rep * n_large_cells] = arr[j];
//      max_rw_small[0][j + rep * n_cell] = arr[j];
//    }
  
    for(int i=1; i<=sim_time; ++i)
    {
      two_step(prtcls.get(),th,rv,diss_rate,opts);
      // get max rw
//      prtcls->diag_max_rw();
//      arr = prtcls->outbuf();
//      int large_cell_idx = -1;
//      for(int j=0; j<n_cell; ++j)
//      {
//  //      real_t large_cell_max_rw;
//        if(j % n_cells_per_avg_r_max_cell == 0)
//        {
//    //      large_cell_max_rw = 0.;
//          large_cell_idx++;
//          max_rw[i][large_cell_idx + rep * n_large_cells] = 0.;//large_cell_max_rw;
//        }
//        if(arr[j] > rep_max_rw) rep_max_rw = arr[j];
//        if(arr[j] > max_rw[i][large_cell_idx + rep * n_large_cells])
//          max_rw[i][large_cell_idx + rep * n_large_cells] = arr[j];
//        max_rw_small[i][j + rep * n_cell] = arr[j];
//
//        // get time for max_rw to reach 40um
//        if(t_max_40[j + rep * n_cell] == 0. && arr[j] >= 40e-6)
//          t_max_40[j + rep * n_cell] = i * opts_init.dt; 
//      }

      // get mean_sd_conc
      prtcls->diag_all();
      prtcls->diag_sd_conc();
      arr = prtcls->outbuf();
      real_t mean_sd_conc = 0;
      for(int j=0; j<n_cell; ++j)
      {
        mean_sd_conc += arr[j]; 
      }
      mean_sd_conc /= real_t(n_cell);
  
  printf("\rrep no: %3d progress: %3d%%: rw_max %lf mean_sd_conc %lf t10_tot %lf", rep, int(real_t(i) / sim_time * 100), rep_max_rw, mean_sd_conc, t10_tot[rep]);
  std::cout << std::flush;
  
      // get t10 (time to conver 10% of cloud water into rain water)
      prtcls->diag_wet_rng(40e-6, 1); // rain water (like in Onishi)
      prtcls->diag_wet_mom(3);
      arr = prtcls->outbuf();
      real_t cloud_mass_tot = 0;
      for(int j=0; j<n_cell; ++j)
      {
//        if(t10[j + rep * n_cell] == 0. && arr[j] >= init_cloud_mass[j] * .1)
//          t10[j + rep * n_cell] = i * opts_init.dt;
//        tau[i][j + rep * n_cell] = arr[j];// / init_cloud_mass[j];
        cloud_mass_tot += arr[j];
      }
      if(t10_tot[rep] == 0. && cloud_mass_tot >= init_tot_cloud_mass * .1)
        t10_tot[rep] =  i * opts_init.dt;

//      prtcls->diag_wet_mom(0);
//      arr = prtcls->outbuf();
//      for(int j=0; j<n_cell; ++j)
//      {
//        nrain[i][j + rep * n_cell] = arr[j];
//      }
//
//      prtcls->diag_wet_rng(0, 1); 
//      prtcls->diag_wet_mom(3);
//      arr = prtcls->outbuf();
//      for(int j=0; j<n_cell; ++j)
//        if(arr[j] > 0) tau[i][j + rep * n_cell] /= arr[j]; // to avoid (small) variability in LWC?
    }
  
    std::cout << std::endl << "po symulacji, max_rw: " << rep_max_rw << std::endl;
  
    diag(prtcls.get(), res_bins_post[rep], res_stddev_bins_post[rep]);
    std::cout << std::endl;
//    auto raw_ptr = prtcls.release();
//    delete raw_ptr;
  }


  // find history in which max_rw grew the most 
//    real_t max_growth=0.;
//    int max_growth_idx;
//    for(int j=0; j < n_cell * n_rep; ++j)
//    {
//      real_t max_rw_local_growth = max_rw_small[sim_time][j] - max_rw_small[0][j];
//      if(max_rw_local_growth > max_growth)
//      {
//        max_growth = max_rw_local_growth;
//        max_growth_idx = j;
//      }
//    }

  real_t mean_max_rad_small = 0.;
  // calc and print max rw (and mass) stats
  for(int i=0; i<=sim_time; ++i)
  {
    real_t glob_max_rad = 0.;
    real_t mean_max_rad = 0.;
    real_t mean_max_vol_small = 0.;
    real_t mean_tau = 0.;
    real_t mean_nrain = 0.;
    real_t std_dev_max_vol_small = 0.;
    real_t std_dev_max_rad = 0.;
    real_t std_dev_max_rad_small = 0.;
    real_t std_dev_tau = 0.;
    real_t skew_max_rad_small = 0.;
    real_t kurt_max_rad_small = 0.;
    real_t skew_max_rad_large = 0.;
    real_t kurt_max_rad_large = 0.;

//    for(int j=0; j < n_large_cells * n_rep; ++j)
//    {
//      mean_max_rad += max_rw[i][j];
//      if(max_rw[i][j] > glob_max_rad) glob_max_rad = max_rw[i][j];
//    }
//    for(int j=0; j < n_cell * n_rep; ++j)
//    {
//      mean_max_vol_small += pow(max_rw_small[i][j], 3.);
//      mean_max_rad_small += max_rw_small[i][j];
//      mean_tau += tau[i][j];
//      mean_nrain += nrain[i][j];
//    }
//
//
//    mean_max_rad /= real_t(n_large_cells * n_rep);
//    mean_max_rad_small /= real_t(n_cell * n_rep);
//    mean_max_vol_small /= real_t(n_cell * n_rep);
//    mean_tau /= real_t(n_cell * n_rep);
//    mean_nrain /= real_t(n_cell * n_rep);
//
//
//    for(int j=0; j < n_large_cells * n_rep; ++j)
//    {
//      std_dev_max_rad += std::pow(max_rw[i][j] - mean_max_rad, 2);
//      skew_max_rad_large += std::pow(max_rw[i][j] - mean_max_rad, 3); 
//      kurt_max_rad_large += std::pow(max_rw[i][j] - mean_max_rad, 4); 
//    }
//    for(int j=0; j < n_cell * n_rep; ++j)
//    {
//      std_dev_max_vol_small += std::pow(pow(max_rw_small[i][j], 3.) / mean_max_vol_small - 1, 2); // relative std dev of mass of largest one
//      std_dev_max_rad_small += std::pow(max_rw_small[i][j] - mean_max_rad_small, 2); 
//      std_dev_tau += std::pow(tau[i][j] - mean_tau, 2); 
//      skew_max_rad_small += std::pow(max_rw_small[i][j] - mean_max_rad_small, 3); 
//      kurt_max_rad_small += std::pow(max_rw_small[i][j] - mean_max_rad_small, 4); 
//    }
//
//    kurt_max_rad_small = kurt_max_rad_small / (n_cell * n_rep) / pow(std_dev_max_rad_small / (n_cell * n_rep), 2.);
//    std_dev_max_vol_small = std::sqrt(std_dev_max_vol_small / (n_cell * n_rep-1));
//    std_dev_max_rad_small = std::sqrt(std_dev_max_rad_small / (n_cell * n_rep-1));
//    std_dev_tau = std::sqrt(std_dev_tau / (n_cell * n_rep-1));
//    skew_max_rad_small = skew_max_rad_small / (n_cell * n_rep) / pow(std_dev_max_rad_small, 3.);
//
//    kurt_max_rad_large = kurt_max_rad_large / (n_large_cells * n_rep) / pow(std_dev_max_rad / (n_large_cells * n_rep), 2.);
//    std_dev_max_rad = std::sqrt(std_dev_max_rad / (n_large_cells * n_rep-1));
//    skew_max_rad_large = skew_max_rad_large / (n_large_cells * n_rep) / pow(std_dev_max_rad, 3.);
//    
//    // 1 - time 2,3 - mean_max_vol_small 4,5 - mean_max_rad(large) 6 - glob max rad 7 - max growth rw
//    // 8 - mean max rw small 9 - std dev max rw small 10 - skew max rw small 11 - kurt max rw small
//    // 12 - skew max rw large 13 - kurt max rw large 14 - mean tau(rain_mass/tot_mass) 
//    // 15 - std_dev tau 16 - mean_nrain
//    of_max_drop_vol << i * dt << " " << mean_max_vol_small << " " << std_dev_max_vol_small << " " << mean_max_rad << " " << std_dev_max_rad << " " << glob_max_rad << " " << max_rw_small[i][max_growth_idx] << " " << mean_max_rad_small << " " << std_dev_max_rad_small << " " << skew_max_rad_small << " " << kurt_max_rad_small << " " << skew_max_rad_large << " " << kurt_max_rad_large << " " << mean_tau << " " << std_dev_tau << " " << mean_nrain << std::endl; 
  } // end of the simulation loop

//
//  // cailc how much the radius of the lucky fraction of small cells increased!!
////all
//  {
//    real_t ensf = n_cell * n_rep;
//    int ens = int(ensf);
//    // cailc how quickkly the lucky fraction of small cells reached r_max=40um
//    std::sort(std::begin(t_max_40), std::end(t_max_40), std::less<real_t>());
//    auto first_nonzero = std::find_if( begin(t_max_40), end(t_max_40), [](real_t x) { return x != 0; });
//    cout << "first nonzero - begin: " << first_nonzero - begin(t_max_40) << endl;
//    if(first_nonzero + ens - begin(t_max_40) > (end(t_max_40)-begin(t_max_40))) 
//      cout << "too short simulation, too small mean sample!" << endl;
//    //else
//    ens = end(t_max_40) - first_nonzero;
//    {
//      real_t lucky_mean_t_max_40 = std::accumulate(first_nonzero, first_nonzero + ens, 0.) / real_t(ens);
//      real_t sq_sum = std::inner_product(first_nonzero, first_nonzero + ens, first_nonzero, 0.0);
//      real_t stdev = std::sqrt(sq_sum / ens - lucky_mean_t_max_40 * lucky_mean_t_max_40);
//      cout << "mean time to reach r=40um: " << lucky_mean_t_max_40 << "  std_dev: " << stdev << endl; 
//    }
//  }
////1e-1
//  {
//    real_t ensf = nx * n_rep / 1e1;
//    int ens = int(ensf);
//    std::sort(std::begin(max_rw_small[sim_time]), std::end(max_rw_small[sim_time]), std::greater<real_t>());
//    real_t lucky_mean_rw = std::accumulate(std::begin(max_rw_small[sim_time]), std::begin(max_rw_small[sim_time]) + ens, 0.) / real_t(ens);
//    cout << "ratoi of lucky 1e-1 fraction final radius to mean final radius: " << lucky_mean_rw / mean_max_rad_small << endl; 
//
//    // cailc how quickkly the lucky fraction of small cells reached r_max=40um
//    auto first_nonzero = std::find_if( begin(t_max_40), end(t_max_40), [](real_t x) { return x != 0; });
//    if(first_nonzero + ens - begin(t_max_40) > (end(t_max_40)-begin(t_max_40))) cout << "too short simulation, too small 1e-1 sample!" << endl;
//    else
//    {
//      real_t lucky_mean_t_max_40 = std::accumulate(first_nonzero, first_nonzero + ens, 0.) / real_t(ens);
//      real_t sq_sum = std::inner_product(first_nonzero, first_nonzero + ens, first_nonzero, 0.0);
//      real_t stdev = std::sqrt(sq_sum / ens - lucky_mean_t_max_40 * lucky_mean_t_max_40);
//      cout << "time for the luckiest 1e-1 to reach r=40um: " << lucky_mean_t_max_40 << "  std_dev: " << stdev << endl; 
//    }
//  }
////1e-2
//  {
//    real_t ensf = nx * n_rep / 1e2;
//    int ens = int(ensf);
//    real_t lucky_mean_rw = std::accumulate(std::begin(max_rw_small[sim_time]), std::begin(max_rw_small[sim_time]) + ens, 0.) / real_t(ens);
//    cout << "ratoi of lucky 1e-2 fraction final radius to mean final radius: " << lucky_mean_rw / mean_max_rad_small << endl; 
//
//    // cailc how quickkly the lucky fraction of small cells reached r_max=40um
//    auto first_nonzero = std::find_if( begin(t_max_40), end(t_max_40), [](real_t x) { return x != 0; });
//    if(first_nonzero + ens - begin(t_max_40) > (end(t_max_40)-begin(t_max_40))) cout << "too short simulation, too small 1e-2 sample!" << endl;
//    else {
//      real_t lucky_mean_t_max_40 = std::accumulate(first_nonzero, first_nonzero + ens, 0.) / real_t(ens);
//      real_t sq_sum = std::inner_product(first_nonzero, first_nonzero + ens, first_nonzero, 0.0);
//      real_t stdev = std::sqrt(sq_sum / ens - lucky_mean_t_max_40 * lucky_mean_t_max_40);
//      cout << "time for the luckiest 1e-2 to reach r=40um: " << lucky_mean_t_max_40 << "  std_dev: " << stdev << endl; 
//    }
//  }
////1e-3
//  {
//    // cailc how much the radius of the lucky fraction of small cells increased!!
//    real_t ensf = nx * n_rep / 1e3;
//    int ens = int(ensf);
//    real_t lucky_mean_rw = std::accumulate(std::begin(max_rw_small[sim_time]), std::begin(max_rw_small[sim_time]) + ens, 0.) / real_t(ens);
//    cout << "ratoi of lucky 1e-3 fraction final radius to mean final radius: " << lucky_mean_rw / mean_max_rad_small << endl; 
//
//    // cailc how quickkly the lucky fraction of small cells reached r_max=40um
//    auto first_nonzero = std::find_if( begin(t_max_40), end(t_max_40), [](real_t x) { return x != 0; });
//    if(first_nonzero + ens - begin(t_max_40) > (end(t_max_40)-begin(t_max_40))) cout << "too short simulation, too small 1e-3 sample!" << endl;
//    else{
//      real_t lucky_mean_t_max_40 = std::accumulate(first_nonzero, first_nonzero + ens, 0.) / real_t(ens);
//      real_t sq_sum = std::inner_product(first_nonzero, first_nonzero + ens, first_nonzero, 0.0);
//      real_t stdev = std::sqrt(sq_sum / ens - lucky_mean_t_max_40 * lucky_mean_t_max_40);
//      cout << "time for the luckiest 1e-3 to reach r=40um: " << lucky_mean_t_max_40 << "  std_dev: " << stdev << endl; 
//    }
//  }
////1e-4
//  {
//    // cailc how much the radius of the lucky fraction of small cells increased!!
//    real_t ensf = nx * n_rep / 1e4;
//    int ens = int(ensf);
//    real_t lucky_mean_rw = std::accumulate(std::begin(max_rw_small[sim_time]), std::begin(max_rw_small[sim_time]) + ens, 0.) / real_t(ens);
//    cout << "ratoi of lucky 1e-4 fraction final radius to mean final radius: " << lucky_mean_rw / mean_max_rad_small << endl; 
//
//    // cailc how quickkly the lucky fraction of small cells reached r_max=40um
//    auto first_nonzero = std::find_if( begin(t_max_40), end(t_max_40), [](real_t x) { return x != 0; });
//    if(first_nonzero + ens - begin(t_max_40) > (end(t_max_40)-begin(t_max_40))) cout << "too short simulation, too small 1e-4 sample!" << endl;
//    else {
//      real_t lucky_mean_t_max_40 = std::accumulate(first_nonzero, first_nonzero + ens, 0.) / real_t(ens);
//      real_t sq_sum = std::inner_product(first_nonzero, first_nonzero + ens, first_nonzero, 0.0);
//      real_t stdev = std::sqrt(sq_sum / ens - lucky_mean_t_max_40 * lucky_mean_t_max_40);
//      cout << "time for the luckiest 1e-4 to reach r=40um: " << lucky_mean_t_max_40 << "  std_dev: " << stdev << endl; 
//    }
//  }
//  // calc and print out mean t10 and t10 std_dev
//  real_t mean_t10 = 0.;
//  int t10_sims = 0;
//  for(int j=0; j<n_cell * n_rep; ++j)
//  {
//    if(t10[j] != 0.)
//    {
//      ++t10_sims;
//      mean_t10 += t10[j];
//    }
//  }
//  std::cout << "no. of cells that didn't reach t10%: " << n_cell * n_rep - t10_sims << std::endl;
//  mean_t10/=t10_sims;//n_cell*n_rep;
//  real_t std_dev=0.;
//  for(int j=0; j<n_cell*n_rep; ++j)
//    std_dev += pow(t10[j] - mean_t10, 2);
//  std_dev /= (t10_sims-1);
//  std_dev = sqrt(std_dev);
//  std::cout << "mean(t10%) = " << mean_t10 << std::endl;
//  std::cout << "realtive std_dev(t10%) = " << std_dev / mean_t10 << std::endl;
//
//  real_t mean_t10_tot;
//  for(int j=0; j<n_rep; ++j)
//    mean_t10_tot += t10_tot[j];
//  mean_t10_tot /= n_rep;
//  std::cout << "mean(t10% in the domain) = " << mean_t10_tot << std::endl;
//
//  // calc and print out mean tau10 and tau10 std_dev
//  const int mean_t10_idx = mean_t10 / dt + 0.5;
//  real_t mean = 0.;
//  if(mean_t10 > 0.)
//  {
//    for(int j=0; j<n_cell*n_rep; ++j)
//    {
//      mean += tau[mean_t10_idx][j];
//    }
//  mean/=n_cell*n_rep;
//  std_dev=0.;
//  for(int j=0; j<n_cell*n_rep; ++j)
//    std_dev += pow(tau[mean_t10_idx][j] - mean, 2);
//  std_dev /= (n_cell*n_rep-1);
//  std_dev = sqrt(std_dev);
//  }
//  else 
//  {
//    mean = NAN;
//    std_dev = NAN;
//  }
//  std::cout << "mean(tau10%) = " << mean << std::endl;
//  std::cout << "realtive std_dev(tau10%) = " << std_dev / mean << std::endl;

// output
  for (int i=0; i <rad_bins.size() -1; ++i)
  {
    real_t rad = (rad_bins[i] + rad_bins[i+1]) / 2.;
    real_t pre = 0, post = 0, stddev_post = 0;
    for(int j=0; j< n_rep; ++j)
    {
      pre += res_bins_pre[j][i];
      post += res_bins_post[j][i];
      stddev_post += res_stddev_bins_post[j][i];
    }
    pre /= n_rep;
    post /= n_rep;
    stddev_post /= n_rep;
    of_size_spectr << rad * 1e6 << " " << pre << " " << post << " " << stddev_post << std::endl; 
  }

//  debug::print(prtcls->impl->n);

//  for(int i=0;i<100;++i)
//  {
//    two_step(prtcls,th,rhod,rv,opts);
//  }

}
