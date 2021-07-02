/**
 * Simple simulation to show model effects.  
 *
 * @author Francesca Capel
 * @date November 2018
 */

functions {
  
  /**
  * Energy spectrum sampling.
  * 
  * @author Francesca Capel
  * @date June 2018
  */

  /**
  * Shape of the energy spectrum: a power law.
  */
  real dNdE_pl(real E, real alpha) {
    
    real spec = pow(E, -alpha);
    return spec; 
    
  }

  /**
  * Sample an energy from a power law spectrum defined by alpha.
  * Uses rejection sampling.
  * Sampled energy is in units of EeV.
  * @param alpha the spectral index
  * @param Emin the minimum energy
  */  
  real spectrum_rng(real alpha, real Emin) {
    
    real E;
    real d;
    real d_upp_lim = dNdE_pl(Emin, alpha);

    int accept = 0;
    
    while(accept != 1) {
      
      E = uniform_rng(Emin, 1e4);
      d = uniform_rng(0, d_upp_lim);
      
      if (d < dNdE_pl(E, alpha)) {
        accept = 1;
      }
    }

    return E;
  }

  /**
  * Log pdf for the power law distribution, bounded
  * below by Emin.
  * @param E energy in EeV
  * @param alpha spectral index
  * @param Emin minimum energy in EeV
  */
  real spectrum_lpdf(real E, real alpha, real Emin) {
    
    real norm = log(alpha - 1) - log(Emin);
    real lprob = -alpha * log(E / Emin);
    
    return lprob + norm;
  }

  /**
  * Continuous energy loss approximation for the 
  * propagation of UHE protons.
  * Based on the implementation in De Domenico & Insolia (2012).
  *
  * @author Francesca Capel
  * @date July 2018
  */

  /**
  * Analytic approx of GZK losses due to photomeson production.
  * Originally parametrised by Anchordorqui et al. (1997)
  * Based on Berezinsky and Grigor'eva (1988) suggestion of Aexp(-B/E) form.
  * @param z the redshift
  * @apram E the energy in eV
  */
  real beta_pi(real z, real E) {
    
    real output;
    real check = 6.86 * exp(-0.807 * z) * 1e20;
    real p[3];
    p[1] = 3.66e-8;
    p[2] = 2.87e20;
    p[3] = 2.42e-8;
    
    if (E <= check) {
      output = p[1] * pow(1 + z, 3) * exp(-p[2] / ((1 + z) * E));
    }
    
    if (E > check) {
      output = p[3] * pow(1 + z, 3);
    }
    
    return output;
  }

  /**
  * Losses due to adiabatic expansion. 
  * Makes use of the hubble constant converted into units of [yr^-1].
  * @param z the redshift
  */
  real beta_adi(real z) {
    
    real lCDM[2];
    real a;
    real b;
    real H0_si = 70; // s^-1 km Mpc^-1
    real H0 = ((H0_si / 3.086e22) * 1.0e3) * 3.154e7; // yr^-1
    lCDM[1] = 0.3;
    lCDM[2] = 0.7;
    
    a = lCDM[1] * pow(1 + z, 3);
    b = 1 - sum(lCDM) * pow(1 + z, 2);
    
    return H0 * pow(a + lCDM[2] + b, 0.5);
  }

  /**
  * Approximation of the get_phi(xi) function as xi->inf.
  * Used in calculation of beta_bh.
  * Described in Chodorowski et al. (1992).
  */
  real phi_inf(real xi) {
    
    real d[4]; 
    real sum_term = 0;
    d[1] = -86.07;
    d[2] = 50.96;
    d[3] = -14.45;
    d[4] = 8.0 / 3.0;
    
    for (i in 1:num_elements(d)) {
      sum_term += d[i] * pow(log(xi), i - 1);
    }
    
    return xi * sum_term;
  }

  /**
  * Approximation of the get_phi(xi) function for different regimes
  * of xi.
  * Used in calculation of beta_bh.
  * Described in Chodorowski et al. (1992). 
  */ 
  real get_phi(real xi) {
    
    real output;
    real sum_term;
    real phi_inf_term;
    real c[4];
    real f[3];
    
    if (xi == 2) { 
      output = (pi() / 12); // * pow(xi - 2, 4);
    }
    
    else if (xi < 25) {
      
      c[1] = 0.8048;
      c[2] = 0.1459;
      c[3] = 1.137e-3;
      c[4] = -3.879e-6;
      sum_term = 0;
      
      for (i in 1:num_elements(c)) {
        sum_term += c[i] * pow(xi - 2, i); 
      }
      
      output = (pi() / 12) * pow(xi - 2, 4) / (1 + sum_term);
    }
    
    else if (xi >= 25) {
      
      f[1] = 2.910;
      f[2] = 78.35;
      f[3] = 1837;
      
      phi_inf_term = phi_inf(xi);
      sum_term = 0;
      for (i in 1:num_elements(f)) {
        sum_term += f[i] * pow(xi, -i);
      }
      
      output = phi_inf_term / (1 - sum_term);
    }
    
    return output;
  }

  /**
  * Integrand for integral over get_phi(xi) in beta_bh calculation. 
  */
  real phi_integrand(real xi, real E, real z) {
    
    real num = get_phi(xi);
    real B = 1.02;
    real denom = exp( (B * xi) / ((1 + z) * E) ) - 1;
    
    return num / denom;
  }

  /**
  * Integrand as a functional of get_phi(xi) used in calcultion of beta_bh.
  * Described in de Domenico and Insolia (2013).
  */
  real[] integrand(real xi, real[] state, real[] params, real[] x_r, int[] x_i) { 
    
    real dstatedxi[1];
    
    real E = params[1];
    real z = params[2];
    
    dstatedxi[1] = phi_integrand(xi, E, z);
    return dstatedxi;
  }

  /**
  * Losses due to Bethe-Heitler pair production.
  * Described in de Domenico amnd Insolia (2013).    
  */
  real beta_bh(real z, real E, real[] xiout, real[] x_r, int[] x_i) {
    
    real params[2];
    real integration_result[1,1];
    real state0[1];
    real integ;
    real A = 3.44e-18;  
      
    state0[1] = phi_integrand(2.0, E, z);
    
    params[1] = E;
    params[2] = z;
    
    integration_result = integrate_ode_rk45(integrand, state0, 2.0, xiout, params, x_r, x_i);
    integ = integration_result[1,1];
    
    return (A  / pow(E, 3)) * integ;
  }

  /**
  * Losses due to Bethe-Heitler pair production.
  * Approximation of inner integral to allow Stan to fit.
  */
  real beta_bh_approx(real z, real E) {
    
    real A = 3.44e-18;
    real B = (3.10075045 * z) - 4.68161976;
    real C = (2.08663474 * z) + 3.8072154;
    real integ = pow(E, 2.4211487 * exp(B / E)) * exp(C);
    
    return (A  / pow(E, 3)) * integ;
    
  }

  /**
  * Total energy losses as a function of distance.
  * @param E the energy in eV (!)
  */
  real dEdr(real r, real E, real[] xiout, real[] x_r, int[] x_i) {
    
    real c = 3.066e-7; // Mpc yr^-1
    real DH = 4285.714; // Mpc
    real z = r / DH;
    real beta_pp;
    real Ltot;
    
    beta_pp = 3.154e7 * beta_bh(z, E / 1.0e18, xiout, x_r, x_i);
    Ltot = c / (beta_adi(z) + beta_pi(z, E) + beta_pp);
    
    return - E / Ltot;
    
  }

  /**
  * ODE system to be solved for arrival energy.
  * NB: E in eV (!)
  */
  real[] E_ode_sim(real r, real[] state, real[] params, real[] x_r, int[] x_i) {
    
    real dstatedr[1];
    
    real E = state[1];
    
    dstatedr[1] = dEdr(r, E, x_r, x_r, x_i);
    return dstatedr;
      
  }

  /**
  * ODE system for be solved for equivalent source energy.
  * Nb: E in eV (!)
  */
  real[] E_ode_rev_sim(real r, real[] state, real[] params, real[] x_r, int[] x_i) {
    
    real dstatedr[1];
    real D = params[1];
    real r_rev = D - r;
    real E = state[1];
    
    dstatedr[1] = - dEdr(r_rev, E, x_r, x_r, x_i);
    return dstatedr;
      
  }

  /**
  * Calculate equivalent energy at the source for an arrival energy of Eth.
  * Solves the ODE system  dE/dr = E / Lloss.
  */
  real get_source_threshold_energy_sim(real Eth, data real[] D, data real[] x_r, int[] x_i) {

    real Eth_src;
    real params[1];
    real integration_result[1, 1];
    real Eth_in[1];
    
    Eth_in[1] = Eth * 1.0e18; // eV
    params[1] = D[1];
    
    integration_result = integrate_ode_rk45(E_ode_rev_sim, Eth_in, 0.0, D, params, x_r, x_i);
    Eth_src = integration_result[1, 1] / 1.0e18; // EeV   
    
    return Eth_src;
  }

  /**
  * Get the vector of source threshold energies for all sources, 
  * including the background component.
  */
  vector get_Eth_src_sim(real Eth, data real[,] D, data real[] x_r, int[] x_i) {

    int N = num_elements(D);
    vector[N] Eth_src;
    
    for (k in 1:N) {
      Eth_src[k] = get_source_threshold_energy_sim(Eth, D[k], x_r, x_i);
    }    
    
    return Eth_src;
    
  }

  /**
  * Calculate the arrival energy taking into account all propagation losses
  * for a given intial energy E.
  * Solves the ODE system dE/dr = - E / Lloss
  */
  real get_arrival_energy_sim(real E, data real[] D, data real[] x_r, int[] x_i) {
    
    real Earr;
    real params[0];
    real integration_result[1, 1];
    real E_in[1];
    
    E_in[1] = E * 1.0e18; // eV
    
    integration_result = integrate_ode_rk45(E_ode_sim, E_in, 0.0, D, params, x_r, x_i);
    Earr = integration_result[1, 1] / 1.0e18; // EeV   
    
    return Earr;   
  }

  /**
  * Total energy losses as a function of distance.
  * Uses beta_bh_approx to allow stan to fit.
  * @param E the energy in eV (!)
  */
  real dEdr_approx(real r, real E) {
    
    real c = 3.066e-7; // Mpc/yr
    real DH = 4285.714; // Mpc
    real z = r / DH;
    real beta_pp;
    real Ltot;
    
    beta_pp = 3.154e7 * beta_bh_approx(z, E / 1.0e18);
    Ltot = c / (beta_adi(z) + beta_pi(z, E) + beta_pp);
    
    return - E / Ltot;
    
  }

  /**
  * ODE system to be solved for arrival energy.
  * NB: E in eV (!)
  */
  real[] E_ode(real r, real[] state, real[] params, real[] x_r, int[] x_i) {
    
    real dstatedr[1];
    
    real E = state[1];
    
    dstatedr[1] = dEdr_approx(r, E);
    return dstatedr;
    
  }

  /**
  * ODE system for be solved for equivalent source energy.
  * Nb: E in eV (!)
  */
  real[] E_ode_rev(real r, real[] state, real[] params, real[] x_r, int[] x_i) {
    
    real dstatedr[1];
    real D = params[1];
    real r_rev = D - r;
    real E = state[1];
    
    dstatedr[1] = - dEdr_approx(r_rev, E);
    return dstatedr;
    
  }

  /**
  * Calculate equivalent energy at the source for an arrival energy of Eth.
  * Solves the ODE system  dE/dr = E / Lloss.
  */
  real get_source_threshold_energy(real Eth, data real[] D, data real[] x_r, int[] x_i) {
    
    real Eth_src;
    real params[1];
    real integration_result[1, 1];
    real Eth_in[1];
    
    Eth_in[1] = Eth * 1.0e18; // eV
    params[1] = D[1];
    
    integration_result = integrate_ode_rk45(E_ode_rev, Eth_in, 0.0, D, params, x_r, x_i);
    Eth_src = integration_result[1, 1] / 1.0e18; // EeV
    
    return Eth_src;
  }

  /**
  * Get the vector of source threshold energies for all sources,
  * including the background component.
  */
  vector get_Eth_src(real Eth, data real[,] D, data real[] x_r, int[] x_i) {
    
    int N = num_elements(D);
    vector[N] Eth_src;
    
    for (k in 1:N) {
      Eth_src[k] = get_source_threshold_energy(Eth, D[k], x_r, x_i);
    }
    
    return Eth_src;
    
  }

  /**
  * Calculate the arrival energy taking into account all propagation losses
  * for a given intial energy E.
  * Solves the ODE system dE/dr = - E / Lloss
  */
  real get_arrival_energy(real E, data real[] D, data real[] x_r, int[] x_i) {
    
    real Earr;
    real params[0];
    real integration_result[1, 1];
    real E_in[1];
    
    E_in[1] = E * 1.0e18; // eV
    
    integration_result = integrate_ode_rk45(E_ode, E_in, 0.0, D, params, x_r, x_i);
    Earr = integration_result[1, 1] / 1.0e18; // EeV
    
    return Earr;
  }

  /**
  * Functions for sampling and fitting the vMF distribution in Stan.
  *
  * @author Francesca Capel
  * @date May 2018
  */

  /**
  * Compute the absolute value of a vector. 
  */
  real abs_val(vector input_vector) {
    
    real av;
    int n = num_elements(input_vector);
    
    real sum_squares = 0;
    for (i in 1:n) {
      sum_squares += (input_vector[i] * input_vector[i]);
    }
    av = sqrt(sum_squares);
    return av;
    
  }

  /**
  * Sample point on sphere orthogonal to mu.
  */
  vector sample_orthonormal_to_rng(vector mu) {
    
    int dim = num_elements(mu);
    vector[dim] v;
    vector[dim] proj_mu_v;
    vector[dim] orthto;
    
    for (i in 1:dim) {
      v[i] = normal_rng(0, 1);
    }
    
    proj_mu_v = mu * dot_product(mu, v) / abs_val(mu);
    orthto = v - proj_mu_v;
    
    return (orthto / abs_val(orthto));
    
  }

  /**
  * Rejection sampling scheme for sampling distance from center on
  * surface of the sphere.
  */
  real sample_weight_rng(real kappa, int dim) {
    
    real sdim = dim - 1; /* as S^{n-1} */
    real b = sdim / (sqrt(4. * pow(kappa, 2) + pow(sdim, 2)) + 2 * kappa);
    real x = (1 - b) / (1 + b);
    real c = kappa * x + sdim * log(1 - pow(x, 2));
    
    int i = 0;
    real z;
    real w;
    real u;

    while (i == 0) {

      z = beta_rng(sdim / 2, sdim / 2);
      w = (1 - (1 + b) * z) / (1 - (1 - b) * z);
      u = uniform_rng(0, 1);

      if (kappa * w + sdim * log(1 - x * w) - c >= log(u)) {
        i = 1;
      }

    }
    
    return w;
  }

  /**
  * Generate an N-dimensional sample from the von Mises - Fisher
  * distribution around center mu in R^N with concentration kappa.
  */
  vector vMF_rng(vector mu, real kappa) {
    
    int dim = num_elements(mu);
    vector[dim] result;
    
    real w = sample_weight_rng(kappa, dim);
    vector[dim] v = sample_orthonormal_to_rng(mu);
    
    result = ( v * sqrt(1 - pow(w, 2)) ) + (w * mu);
    return result;
    
  }

  /**
  * Sample a point uniformly from the surface of a sphere of 
  * a certain radius.
  */
  vector sphere_rng(real radius) {
    
    vector[3] result;
    real u = uniform_rng(0, 1);
    real v = uniform_rng(0, 1);
    real theta = 2 * pi() * u;
    real phi = acos( (2 * v) - 1 );
    
    result[1] = radius * cos(theta) * sin(phi); 
    result[2] = radius * sin(theta) * sin(phi); 
    result[3] = radius * cos(phi);
    
    return result;
    
  }

  /**
  * Functions for modelling the exposure of a UHECR observatory.
  *
  * @author Francesca Capel
  * @date May 2018
  */


  /**
  * Calculate xi part of exposure.
  * @param theta from 0 to pi.
  * @param p observatory dependent parameters.
  */
  real xi_exp(real theta, real[] p) { 
    return (p[3] - (p[2] * cos(theta))) / (p[1] * sin(theta));
  }

  /**
  * Calculate alpha_m part of exposure.
  * @param theta from 0 to pi.
  * @param p observatory dependent parameters.
  */
  real alpha_m(real theta, real[] p) {
    
    real am;
    
    real xi_val = xi_exp(theta, p);
    if (xi_val > 1) {
      am = 0;
    }
    else if (xi_val < -1) {
      am = pi();
    }
    else {
      am = acos(xi_val);
    }
    
    return am;
  }

  /**
  * Calculate the exposure factor for a given position on the sky. 
  * @param theta from 0 to pi.
  * @param p observatory dependent parameters.
  */
  real m(real theta, real[] p) {
    return (p[1] * sin(theta) * sin(alpha_m(theta, p)) 
      + alpha_m(theta, p) * p[2] * cos(theta));
  }

  /**
  * Utils for use with Stan simulations and models.
  *
  * @author Francesca Capel
  * @date October 2018
  */
    
  /**
  * Calculate weights from source distances.
  */
  vector get_source_weights(real[] L, real[] D) {
    
    int N = num_elements(D);
    vector[N] weights;
    
    real normalisation = 0;
    
    for (k in 1:N) {
      normalisation += (L[k] / pow(D[k], 2));
    }
    for (k in 1:N) {
      weights[k] = (L[k] / pow(D[k], 2)) / normalisation;
    }
    
    return weights;
  }

  /**
  * Calculate weights for each source accounting for exposure 
  * and propagation effects.
  */
  vector get_exposure_weights(vector F, vector eps, real alpha_T, vector Eth_src, real Eth, real alpha) {
    
    int N = num_elements(F);
    vector[N] weights;
    
    real normalisation = 0;
    
    for (k in 1:N-1) {
      normalisation += F[k] * eps[k] * pow(Eth_src[k] / Eth, 1 - alpha);
    }
    normalisation += F[N] * (alpha_T / (4 * pi()));
    
    for (k in 1:N-1) {
      weights[k] = (F[k] * eps[k] * pow(Eth_src[k] / Eth, 1 - alpha)) / normalisation;
    }
    weights[N] = (F[N] * (alpha_T / (4 * pi()))) / normalisation;
    
    return weights;
  }

  /**
  * Convert from unit vector omega to theta of spherical coordinate system.
  * @param omega a 3D unit vector.
  */
  real omega_to_theta(vector omega) {
    
    real theta;
    
    int N = num_elements(omega);
    
    if (N != 3) {
      print("Error: input vector omega must be of 3 dimensions");
    }
    
    theta = acos(omega[3]);
    
    return theta;
  }

  /**
  * Calculate the expected value of N for the generative model.
  */
  real get_Nex_sim(vector F, vector eps, real alpha_T, vector Eth_src, real Eth, real alpha) {
    
    int N = num_elements(F);
    real Nex = 0;
    
    for (k in 1:N-1) {
      Nex += F[k] * eps[k] * pow(Eth_src[k] / Eth, 1 - alpha);
    }
    Nex += F[N] * (alpha_T / (4 * pi()));

    return Nex;
  }
    
  /**
  * Calculate the total source flux.
  * @param L the luminosity in s^-1
  * @param D the distance in Mpc
  */
  real get_Fs(real[] L, real[] D) {
    
    int N = num_elements(D);
    real Fs = 0;

    for (k in 1:N) {
      Fs += L[k] / (4 * pi() * pow(D[k], 2));
    }
    
    return Fs;
  }

  /**
  * Sample from the vMF centred on varpi with spread kappa, 
  * accounting for detector exposure effects.
  * Uses rejection sampling.
  */
  vector exposure_limited_vMF_rng(vector varpi, real kappa, real a0, real theta_m) {
    
    real params[3];
    real m_max;
    real accept;
    int count;
    vector[3] omega;
    real theta;
    real pdet;
    vector[2] p;
    
    /* exposure */
    params[1] = cos(a0);
    params[2] = sin(a0);
    params[3] = cos(theta_m);
    m_max = m(pi(), params);
    
    accept = 0;
    count = 0;

    while (accept != 1) {	  

      omega = vMF_rng(varpi, kappa);
      theta = omega_to_theta(omega);
      pdet = m(theta, params) / m_max;
      p[1] = pdet;
      p[2] = 1 - pdet;
      accept = categorical_rng(p);
      count += 1;
      if (count > 1.0e7) {
        
        print("Was stuck in exposure_limited_rng");
        accept = 1;

      }
      
    } 
    
      return omega;
    }

  /**
  * Sample uniformly from the unit sphere, 
  * accounting for detector exposure effects.
  * Uses rejection sampling.
  */
  vector exposure_limited_sphere_rng(real a0, real theta_m) {
    
    real params[3];
    real m_max;
    real accept;
    vector[3] omega;
    real theta;
    real pdet;
    vector[2] p;
    
    params[1] = cos(a0);
    params[2] = sin(a0);
    params[3] = cos(theta_m);
    m_max = m(pi(), params);
    
    accept = 0;
    while (accept != 1) {	  

      omega = sphere_rng(1);
      theta = omega_to_theta(omega);
      pdet = m(theta, params) / m_max;
      p[1] = pdet;
      p[2] = 1 - pdet;
      accept = categorical_rng(p);

    } 

      return omega;
    }


  /**
  * Interpolate x from a given set of x and y values.
  */
  real interpolate(vector x_values, vector y_values, real x) {
    real x_left;
    real y_left;
    real x_right;
    real y_right;
    real dydx;
    
    int Nx = num_elements(x_values);
    real xmin = x_values[1];
    real xmax = x_values[Nx];
    int i = 1;
    
    if (x > xmax || x < xmin) {
      
      if(x > xmax) {
        return y_values[Nx];
      }
      else if (x < xmin) {
        return y_values[1];
      }
    }
    
    if( x >= x_values[Nx - 1] ) {
      i = Nx - 1;
    }
    else {
      while( x > x_values[i + 1] ) { i += 1; }
    }
    
    x_left = x_values[i];
    y_left = y_values[i];
    x_right = x_values[i + 1];
    y_right = y_values[i + 1];
    
    dydx = (y_right - y_left) / (x_right - x_left);
    
    return y_left + dydx * (x - x_left);
  }

  /**
  * Calculate the Nex for a given kappa by
  * interpolating over a vector of eps values
  * for each source.
  */
  real get_Nex(vector F, vector[] eps, vector kappa_grid, vector kappa, real alpha_T, vector Eth_src, real Eth, real alpha) {
    
    int Ns = num_elements(F);
    vector[Ns] N;
    real eps_from_kappa;
    
    for (k in 1:Ns-1) {
      eps_from_kappa = interpolate(kappa_grid, eps[k], kappa[k]);
      N[k] = F[k] * eps_from_kappa * pow(Eth_src[k] / Eth, 1 - alpha); 
    }
    N[Ns] = F[Ns] * (alpha_T / (4 * pi()));
    
    return sum(N);
  }

  /**
  * Calculate the Nex for a given kappa by
  * interpolating over a vector of eps values
  * for each source.
  * This version ignores energy information and is 
  * to be used in the arrival direction only models.
  */
  real get_Nex_arr(vector F, vector[] eps, vector kappa_grid, real kappa, real alpha_T) {
    
    int Ns = num_elements(F);
    vector[Ns] N;
    real eps_from_kappa;
    
    for (k in 1:Ns-1) {
      eps_from_kappa = interpolate(kappa_grid, eps[k], kappa);
      N[k] = F[k] * eps_from_kappa; 
    }
    N[Ns] = F[Ns] * (alpha_T / (4 * pi()));
    
    return sum(N);
  }


  /**
  * Define the fik PDF.
  * NB: Cannot be vectorised.
  * Uses sinh(kappa) ~ exp(kappa)/2 
  * approximation for kappa > 100.
  */
  real fik_lpdf(vector v, vector mu, real kappa, real kappa_d) {
    
    real lprob;
    real inner = abs_val((kappa_d * v) + (kappa * mu));
    
    if (kappa > 100 || kappa_d > 100) {
      lprob = log(kappa * kappa_d) - log(4 * pi() * inner) + inner - (kappa + kappa_d) + log(2);
    }
    else {   
      lprob = log(kappa * kappa_d) - log(4 * pi() * sinh(kappa) * sinh(kappa_d)) + log(sinh(inner)) - log(inner);
    }
    
    return lprob;   
  }

  /**
  * Calculate the deflection parameter of the vMF distribution.
  * @param E energy in EeV
  * @param B rms magnetic field strength in nG
  * @param D distance in Mpc / 10
  */
  real get_kappa(real E, real B, real D) {
    
    return 2.3 * inv_square( 0.0401 * inv(E / 50) * B * sqrt(D) );
  }

  /**
  * Calculate the vector of kappa values for different sources.
  * @param E energy in EeV
  * @param B rms magnetic field strength in nG
  * @param D distance in Mpc / 10
  */
  vector get_kappa_ex(vector E, real B, vector D) {
    
    int Ns = num_elements(E);
    vector[Ns] kappa_ex = 2.3 * inv_square( 0.0401 * inv(E / 50) * B .* sqrt(D) );
    return kappa_ex;
  }

  /**
  * Calculate the vector of expected E values for different sources.
  */
  vector get_Eex(real alpha, vector Eth_src) {
    
    int N = num_elements(Eth_src);
    vector[N] Eex = pow(2, 1 / (alpha - 1)) * Eth_src[1:N];
    
    return Eex;
  }

  vector get_eps_from_kappa(vector kappa_grid, vector[] eps, vector kappa_ex) {
    
    int N = num_elements(kappa_ex);
    vector[N] eps_from_kappa;
    
    for (k in 1:N) {
      eps_from_kappa[k] = interpolate(kappa_grid, eps[k], kappa_ex[k]);
    }
    
    return eps_from_kappa;
  } 
  
}


data {

  /* source position */
  unit_vector[3] varpi;
  
  /* source spectrum */
  real alpha;
  real<lower=0> Eth; // EeV
  real<lower=0> Eerr; // EeV
  
  /* flux */
  int<lower=0> N;
 
  /* deflection */
  real<lower=0> B; // nG
  real<lower=0> kappa_d;

}

transformed data {
  
  real x_r[1];
  int x_i[0];

  /* approximate infinity */
  x_r[1] = 1.0e4;
 
}

generated quantities {

  unit_vector[3] omega[N];
  unit_vector[3] omega_det[N];
  
  real E[N];
  real kappa[N];
  real Earr[N];
  real Edet[N];
  real Eth_src[N];
  real D[N];
  real D_in[N, 1];
  
  for (i in 1:N) {

    /* Uniformly distributed sources */
    D[i] = uniform_rng(0, 400);
    D_in[i, 1] = D[i];
    Eth_src[i] = get_source_threshold_energy_sim(Eth, D_in[i], x_r, x_i);

    /* Energy sampling and propagation */
    E[i] = spectrum_rng(alpha, Eth_src[i]);
    kappa[i] = get_kappa(E[i], B, D[i] / 10);
    omega[i] = vMF_rng(varpi, kappa[i]);      
    Earr[i] = get_arrival_energy_sim(E[i], D_in[i], x_r, x_i);
      
    /* Detection effects */
    omega_det[i] = vMF_rng(omega[i], kappa_d);  	  
    Edet[i] = normal_rng(Earr[i], Eerr * Earr[i]);

  }
  
}

