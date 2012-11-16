<!DOCTYPE HTML PUBLIC "-//IETF//DTD HTML//EN">
<html><head>
<meta http-equiv="Content-Type" CONTENT="text/html; charset=iso-8859-1">
<meta name="Author" CONTENT="Leon Kos">
<meta name="Keywords" CONTENT="IAPWS IF-74, properties of water and steam, javascript">

<title>IAPWS IF-97</title>

<script language=javascript>
// Copyright 2001, Leon Kos
// All properties in SI units

var R = 0.461526e3;      // J kg^-1 K^-1

function Equation(T_reducing, p_reducing, n, J, I)
{
  this.T_reducing = T_reducing;
  this.p_reducing = p_reducing; // pressure or density
  this.n = n;          // coefficients
  this.J = J || null;  // exponents J
  this.I = I || null;  // exponents I
  this.tau = function (temperature)
  { return this.T_reducing/temperature; }
  this.pi =  function (pressure)
  {return pressure/this.p_reducing; }
  this.delta = function (density)
  {return density/this.p_reducing; }
}


// EQUATION DATA

var equation_3_1 = new Equation( 1386, 16.53e6, // Region 1, basic equation 3.1
  [ 0.14632971213167,     -0.84548187169114,    -0.37563603672040e1,
    0.33855169168385e1,   -0.95791963387872,     0.15772038513228,
   -0.16616417199501e-1,   0.81214629983568e-3,  0.28319080123804e-3,
   -0.60706301565874e-3,  -0.18990068218419e-1, -0.32529748770505e-1,
   -0.21841717175414e-1,  -0.52838357969930e-4, -0.47184321073267e-3,
   -0.30001780793026e-3,   0.47661393906987e-4, -0.44141845330846e-5,
   -0.72694996297594e-15, -0.31679644845054e-4, -0.28270797985312e-5,
   -0.85205128120103e-9,  -0.22425281908000e-5, -0.65171222895601e-6,
   -0.14341729937924e-12, -0.40516996860117e-6, -0.12734301741641e-8,
   -0.17424871230634e-9,  -0.68762131295531e-18, 0.14478307828521e-19,
    0.26335781662795e-22, -0.11947622640071e-22, 0.18228094581404e-23,
   -0.93537087292458e-25 ],
  [-2, -1, 0, 1, 2, 3, 4, 5, -9, -7, -1, 0, 1, 3, -3, 0, 1, 3, 17,
  -4, 0, 6, -5, -2, 10, -8, -11, -6, -29, -31, -38, -39, -40, -41],
  [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3,
   4, 4, 4, 5, 8, 8, 21, 23, 29, 30, 31, 32]
  );

// Region 2

var equation_3_9 = new Equation ( 540, 1e6, // Region 2, equation 3.9
 [ -0.96927686500217E+01,  0.10086655968018E+02, -0.56087911283020E-02,
    0.71452738081455E-01, -0.40710498223928E+00,  0.14240819171444E+01, 
   -0.43839511319450E+01, -0.28408632460772E+00,  0.21268463753307E-01],
 [0, 1, -5, -4, -3, -2, -1, 2,  3]);

var equation_3_10 = new Equation ( 540, 1e6, // Region 2, eq 3.10, residual
 [ -0.17731742473213E-02,  -0.17834862292358E-01,  -0.45996013696365E-01, 
   -0.57581259083432E-01,  -0.50325278727930E-01,  -0.33032641670203E-04, 
   -0.18948987516315E-03,  -0.39392777243355E-02,  -0.43797295650573E-01, 
   -0.26674547914087E-04,   0.20481737692309E-07,   0.43870667284435E-06, 
   -0.32277677238570E-04,  -0.15033924542148E-02,  -0.40668253562649E-01, 
   -0.78847309559367E-09,   0.12790717852285E-07,   0.48225372718507E-06, 
    0.22922076337661E-05,  -0.16714766451061E-10,  -0.21171472321355E-02, 
   -0.23895741934104E+02,  -0.59059564324270E-17,  -0.12621808899101E-05, 
   -0.38946842435739E-01,   0.11256211360459E-10,  -0.82311340897998E+01, 
    0.19809712802088E-07,   0.10406965210174E-18,  -0.10234747095929E-12, 
   -0.10018179379511E-08,  -0.80882908646985E-10,   0.10693031879409E+00, 
   -0.33662250574171E+00,   0.89185845355421E-24,   0.30629316876232E-12, 
   -0.42002467698208E-05,  -0.59056029685639E-25,   0.37826947613457E-05, 
   -0.12768608934681E-14,   0.73087610595061E-28,   0.55414715350778E-16, 
   -0.94369707241210E-06 ],
 [0, 1, 2, 3, 6, 1, 2, 4, 7, 36, 0, 1, 3, 6, 35, 1, 2, 3, 7, 3, 16, 35, 0, 11,
 25, 8, 36, 13, 4, 10, 14, 29, 50, 57, 20, 35, 48, 21, 53, 39, 26, 40, 58 ],
 [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 5, 6, 6, 6, 7, 7, 7,
   8, 8, 9, 10, 10, 10, 16, 16, 18, 20, 20, 20, 21, 22, 23, 24, 24, 24 ]
);

// Region 3, basic equation 3.21 , note that second parameter is density
var equation_3_21 = new Equation (647.096, 322, 
 [ 0.10658070028513E+01, -0.15732845290239E+02,  0.20944396974307E+02,
  -0.76867707878716E+01,  0.26185947787954E+01, -0.28080781148620E+01,
   0.12053369696517E+01, -0.84566812812502E-02, -0.12654315477714E+01,
  -0.11524407806681E+01,  0.88521043984318E+00, -0.64207765181607E+00,
   0.38493460186671E+00, -0.85214708824206E+00,  0.48972281541877E+01,
  -0.30502617256965E+01,  0.39420536879154E-01,  0.12558408424308E+00,
  -0.27999329698710E+00,  0.13899799569460E+01, -0.20189915023570E+01,
  -0.82147637173963E-02, -0.47596035734923E+00,  0.43984074473500E-01,
  -0.44476435428739E+00,  0.90572070719733E+00,  0.70522450087967E+00,
   0.10770512626332E+00, -0.32913623258954E+00, -0.50871062041158E+00,
  -0.22175400873096E-01,  0.94260751665092E-01,  0.16436278447961E+00,
  -0.13503372241348E-01, -0.14834345352472E-01,  0.57922953628084E-03,
   0.32308904703711E-02,  0.80964802996215E-04, -0.16557679795037E-03,
  -0.44923899061815E-04 ],
  [0, 0, 1, 2, 7, 10, 12, 23, 2, 6, 15, 17, 0, 2, 6, 7, 22, 26, 0, 2, 4,
  16, 26, 0, 2, 4, 26, 1, 3, 26, 0, 2, 26, 2, 26, 2, 26, 0, 1, 26],
  [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3,
   4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 8, 9, 9, 10, 10, 11]
);


// REGION 1
// specific volume
function v_1(T, p)
{
  var dg_pi = 0.0;
  
  with (equation_3_1)
    {
      for (var i = 0; i < 34; i++)
        dg_pi += -n[i] * I[i] * Math.pow(7.1 - pi(p), I[i] - 1)
          * Math.pow(tau(T) - 1.222, J[i]);
      
      return pi(p) * dg_pi * R * T / p;
    }
}


// specific internal energy
function u_1(T, p)
{
  var dg_pi = 0.0;
  var dg_tau = 0.0;
  with ( equation_3_1 )
    {
      for (var i = 0; i < 34; i++)
        {
          dg_pi += -n[i] * I[i] * Math.pow(7.1 - pi(p), I[i] - 1)
            * Math.pow(tau(T) - 1.222, J[i]);
          dg_tau += n[i] * Math.pow( 7.1 - pi(p), I[i])
            * J[i] * Math.pow(tau(T) - 1.222, J[i] - 1);
        }
      return R * T * (tau(T) * dg_tau - pi(p) * dg_pi);
    }
}



// specific entropy
function s_1(T, p)
{
  var dg = 0.0;
  var dg_tau = 0.0;
  with ( equation_3_1 )
    {
      for ( var i = 0; i < 34; i++)
        {
          dg += n[i] * Math.pow( 7.1 - pi(p), I[i] )
            * Math.pow(tau(T) - 1.222, J[i]);
          dg_tau += n[i] * Math.pow( 7.1 - pi(p), I[i])
            * J[i] * Math.pow(tau(T) - 1.222, J[i] - 1);
        }
      return R * (tau(T) * dg_tau - dg);
    }
}


// specific_enthalpy
function h_1(T, p)
{
  var dg_tau = 0.0;
  with ( equation_3_1 )
    {
      for (var i = 0; i < 34; i++)
        dg_tau += n[i] * Math.pow( 7.1 - pi(p), I[i])
          * J[i] * Math.pow(tau(T) - 1.222, J[i] - 1);
      return R * T * tau(T) * dg_tau;
    }
}

// specific_isobaric_heat_capacity
function cp_1(T, p)
{
  var dg_tt = 0.0;
  with ( equation_3_1 )
    {
      for ( var i = 0; i < 34; i++)
        dg_tt += n[i] * Math.pow( 7.1 - pi(p), I[i])
          * J[i] * ( J[i] - 1) * Math.pow(tau(T) - 1.222, J[i] - 2);

      return -R * tau(T) * tau(T) * dg_tt;
    }
}

// specific isochoric heat capacity
function cv_1(T, p)
{
  var dg_tt = 0.0;
  var dg_pi = 0.0;
  var dg_pt = 0.0;
  var dg_pp = 0.0;
  function sqr(x) { return x * x;}
  
  with ( equation_3_1)
    {
      for ( var i = 0; i < 34; i++)
        {
          dg_tt += n[i] * Math.pow( 7.1 - pi(p), I[i])
            * J[i] * ( J[i] - 1) * Math.pow(tau(T) - 1.222, J[i] - 2);
          dg_pi += -n[i] * I[i] * Math.pow(7.1 - pi(p), I[i] - 1)
            * Math.pow(tau(T) - 1.222, J[i]);
          dg_pt += -n[i] * I[i] * Math.pow(7.1 - pi(p), I[i] - 1)
            * J[i] * Math.pow( pi(p) - 1.222, J[i] - 1);
          dg_pp += n[i] * I[i] * (I[i] - 1) * Math.pow(7.1 - pi(p), I[i] - 2)
            * Math.pow(tau(T) - 1.222, J[i]);
        }
      return R * (-sqr(tau(T)) * dg_tt + sqr(dg_pi - tau(T)*dg_pt)/dg_pp);
    }
}


//REGION 2
// Specific volume
function v_2(T, p)
{
  var g0_pi;
  var gr_pi = 0.0;
  
  with (equation_3_9)
    {
      g0_pi = 1/pi(p);
    }
  with (equation_3_10)
    {
      for (var i = 0; i < 43; i++)
        gr_pi += n[i] * I[i] * Math.pow(pi(p), I[i] - 1)
          *  Math.pow(tau(T) - 0.5, J[i]);
      return R * T * pi(p) * (g0_pi + gr_pi) / p;
    }
}

//Specific internal energy
function u_2(T, p)
{
  var g0_tau = 0.0;
  var gr_tau = 0.0;
  var g0_pi;
  var gr_pi = 0.0;
  with (equation_3_9)
    {
      g0_pi = 1/pi(p);
      for (var i = 0; i < 9; i++)
        g0_tau += n[i] * J[i] * Math.pow(tau(T), J[i] - 1);
    }
  with (equation_3_10)
    {
      for (var i = 0; i < 43; i++)
        {
          gr_pi += n[i] * I[i] * Math.pow(pi(p), I[i] - 1)
            *  Math.pow(tau(T) - 0.5, J[i]);
          gr_tau += n[i] * Math.pow(pi(p), I[i])
            * J[i] * Math.pow(tau(T) - 0.5, J[i] - 1);
        }
      return R * T * (tau(T)*(g0_tau+gr_tau) - pi(p)*(g0_pi+gr_pi)) ;
    }     
}

function s_2(T, p)
{
  var g0 = Math.log(equation_3_9.pi(p));
  var gr = 0.0;
  var g0_tau = 0.0;
  var gr_tau = 0.0;
  with (equation_3_9)
    {
      for (var i = 0; i < 9; i++)
        {
          g0 += n[i] * Math.pow(tau(T), J[i]);
          g0_tau += n[i] * J[i] * Math.pow(tau(T), J[i] - 1);
        }
    }
  with (equation_3_10)
    {
      for (var i = 0; i < 43; i++)
        {
          gr += n[i] * Math.pow(pi(p), I[i])
            *  Math.pow(tau(T) - 0.5, J[i]);
          gr_tau += n[i] * Math.pow(pi(p), I[i])
            * J[i] * Math.pow(tau(T) - 0.5, J[i] - 1);
        }
      return R * (tau(T)*(g0_tau+gr_tau) - (g0 +gr)) ;
    }     
}

// specific enthalpy
function h_2(T, p)
{
  var g0_tau = 0.0;
  var gr_tau = 0.0;
  with (equation_3_9)
    {
      for (var i = 0; i < 9; i++)
        g0_tau += n[i] * J[i] * Math.pow(tau(T), J[i] - 1);
    }
  with (equation_3_10)
    {
      for (var i = 0; i < 43; i++)
          gr_tau += n[i] * Math.pow(pi(p), I[i])
            * J[i] * Math.pow(tau(T) - 0.5, J[i] - 1);
      return R * T * tau(T) * ( g0_tau + gr_tau ) ;
    }       
}

// specific isobaric heat capacity
function cp_2(T, p)
{
  var g0_tt = 0.0;
  var gr_tt = 0.0;
  with (equation_3_9)
    {
      for (var i = 0; i < 9; i++)
        g0_tt += n[i] * J[i] * (J[i] - 1) * Math.pow(tau(T), J[i] - 2);
    }
  with (equation_3_10)
    {
      for (var i = 0; i < 43; i++)
        gr_tt += n[i] * Math.pow(pi(p), I[i]) * J[i] * (J[i] - 1)
          * Math.pow(tau(T) - 0.5, J[i] - 2);
      return -R * tau(T) * tau(T) * ( g0_tt + gr_tt ) ;
    }       
}

// REGION 2,3

// Auxiliary equation between Regions 2 and 3, Equ. 3.6 and Equ. 3.7
var n23 = [ 0.34805185628969E+03, -0.11671859879975E+01,
          0.10192970039326E-02,  0.57254459862746E+03, 0.13918839778870E+02];

function p_B(T) // equation 3.6
{
  var n = n23;
  return 1e6 * (n[0] + n[1] * T + n[2] * T * T );
}

function T_B(p) // equation 3.7
{
  var n = n23;
  return  n[3] + Math.sqrt((p/1e6 - n[4])/n[2]);
}  


// REGION 3
// pressure
function p_3(T, rho)
{
  var dh_d;
  with (equation_3_21)
    {
      dh_d = n[0]/delta(rho);
      for (var i = 1; i < 40; i++)
        dh_d += n[i] * I[i] * Math.pow(delta(rho), I[i] - 1)
          * Math.pow(tau(T), J[i]);
      return rho * R * T * delta(rho) * dh_d;
    }
}

// density must be calculated iteratively with Newton method
function rho_3(T, p)
{
  var densold, diffdens;
  var dc = 322; // critical density in kg m^-33
  var Tc = 647.096; //critical temperature in K
  densold = 750; // a good downhill starting point
  for (var j = 0; j < 1000; j++)
    {
      var dh_d;
      var dh_dd;
      with (equation_3_21)
        {
          dh_d = n[0] / delta(densold);
          dh_dd = -n[0] * dc / densold / densold;
          for ( var i = 1; i < 40; i ++)
            {
              dh_d += n[i] * I[i] * Math.pow(delta(densold), I[i] - 1)
                * Math.pow(tau(T), J[i]);
              dh_dd += n[i] * I[i] * (I[i]-1) /dc
                * Math.pow(delta(densold), I[i] - 2) * Math.pow(tau(T), J[i]);
            }
        }
      var derivprho = R * T / dc * (2 * densold * dh_d
                                    + densold * densold * dh_dd);
      var densnew = densold + (p - R * T * densold * densold / dc * dh_d )
        / derivprho;
      if (densnew < 0) // sanity check
        densnew = 700;
      var diffdens = Math.abs(densnew-densold)/densnew;
      if ( diffdens < 0.0000000005)
        return densnew;
      densold = densnew;
    }
  alert("Accuracy problem: " + diffdens );
  return densold;
  return -1 ; // fault
}

// specific internal energy
function u_3(T, rho)
{
  var dh_tau = 0.0;
  with (equation_3_21)
    {
      for (var i = 1; i < 40; i++)
        dh_tau += n[i] * Math.pow(delta(rho), I[i])
          * J[i] * Math.pow(tau(T), J[i] - 1);
      return R * T * tau(T) * dh_tau;
    }
}

// specific entropy
function s_3(T, rho)
{
  var dh;
  var dh_tau = 0.0;
  with (equation_3_21)
    {
      dh = n[0] * Math.log(delta(rho));
      for (var i = 1; i < 40; i++)
        {
          dh += n[i] * Math.pow(delta(rho), I[i])
            * Math.pow(tau(T), J[i]);
          dh_tau += n[i] * Math.pow(delta(rho), I[i])
            * J[i] * Math.pow(tau(T), J[i] - 1);
        }
      return R * (tau(T) * dh_tau - dh);
    }
}

// specific enthalpy
function h_3(T, rho)
{
  var dh_d;
  var dh_tau = 0.0;
  with (equation_3_21)
    {
      dh_d = n[0]/delta(rho);
      for (var i = 1; i < 40; i++)
        {
          dh_d += n[i] * I[i] * Math.pow(delta(rho), I[i] - 1)
            * Math.pow(tau(T), J[i]);
          dh_tau += n[i] * Math.pow(delta(rho), I[i])
            * J[i] * Math.pow(tau(T), J[i] - 1);
        }
      return R * T * (tau(T) * dh_tau + delta(rho) * dh_d);
    }
}

// specific isochoric heat capacity
function cv_3(T, rho)
{
  var dh_tt = 0.0;
  with (equation_3_21)
    {
      for (var i = 1; i < 40; i++)
        dh_tt += n[i] *  Math.pow(delta(rho), I[i] )
          * J[i] * (J[i] - 1 ) * Math.pow(tau(T), J[i] - 2);
      return -R * tau(T) * tau(T) * dh_tt;
    }
}


// REGION 4
// coefficients of equations 3.22 to 3.24
var r4_coeff =  [ 0.11670521452767E+04,-0.72421316703206E+06,
                -0.17073846940092E+02,  0.12020824702470E+05,
                -0.32325550322333E+07,  0.14915108613530E+02,
                -0.48232657361591E+04,  0.40511340542057E+06,
                -0.23855557567849E+00,  0.65017534844798E+03];

// saturation pressure
function p_s(T)
{
  var T_reducing = 1;
  var p_reducing = 1e6;
  var n = r4_coeff;
  var theta = T/T_reducing +  n[8]/((T/T_reducing)-n[9]);
  var A = theta*theta + n[0] * theta + n[1];
  var B = n[2] * theta*theta + n[3] * theta + n[4];
  var C = n[5]*theta*theta + n[6] * theta + n[7];
  function quad(x) { return x*x*x*x; }
  return p_reducing*quad(2*C/(-B+Math.sqrt(B*B-4*A*C)));
}


// for region 3 parameter p_rho means density elsewere it means pressure
function Region(v, rho, u, s, h, cp, cv, w)
{
  this.v =  v || function (T, p) { return null; } ;
  this.rho = rho || function (T, p) { return null; } ;
  this.u = u;
  this.s = s;
  this.h = h;
  this.cp = cp || function (T, p_rho) { return null };
  this.cv = cv || function (T, p_rho) { return null };
  this.w = w || function (T, p_rho) { return null };
  this.toHTML = function(T, p)
  {
    var p_rho = p;
    var r = new String();
    with (this)
    {
      if ( v(T, p) == null) // region3
        {
          p_rho = rho(T, p);
          r = "rho = " + p_rho + " kg m<sup>-3</sup><br>\n";
        }
      else
        r =  "v = " + v(T, p_rho) + " m<sup>3</sup> kg<sup>-1</sup><br>\n";
      r += "h = " + h(T, p_rho)/1000 + " kJ kg<sup>-1</sup><br>\n"
        + "u = " + u(T, p_rho)/1000 + " kJ kg<sup>-1</sup><br>\n"
        + "s = " + s(T, p_rho)/1000
        + " kJ kg<sup>-1</sup> K<sup>-1</sup><br>\n";
      
      if ( cp(T,p_rho) != null)
        r += "c<sub>p</sub> = " + cp(T, p_rho)/1000
          + " kJ kg<sup>-1</sup> K<sup>-1</sup><br>\n";
      if ( cv(T,p_rho) != null )
        r += "c<sub>v</sub> = " + cv(T, p_rho)/1000
          + " kJ kg<sup>-1</sup> K<sup>-1</sup><br>\n";
    }
    return r;
  }
}

var region = new Array();
region[1] = new Region(v_1, null, u_1, s_1, h_1, cp_1);
region[2] = new Region(v_2, null, u_2, s_2, h_2, cp_2);
region[3] = new Region(null, rho_3, u_3, s_3, h_3, null, cv_3);

// returns region number for a given temperature and pressure
function which_region(T, p)
{
 if (273.15 <= T && T <= 623.15 && p_s(T) <= p && p <= 100e6)
   return 1;
 else
   if ( (273.15 <= T &&  T <= 623.15 &&  0 < p && p <= p_s(T))
        || (623.15 < T  && T <= 863.15 && 0 < p && p <= p_B(T))
        || (863.15 < T && T <= 1073.15 && 0 < p && p <= 100e6) )
     return 2;
   else
     if ( 623.5 <= T && T <= T_B(p) && p_B(T) <= p && p <= 100e6)
       return 3;


 return 0; // out of range  
}


//User interface functions

// New window display
function display_message(title, html_text)
{
  var w = window.open("", null,  "width=300,height=300,"
                      + "menubar=yes,scrollbars=yes,resizable=yes");
  var d = w.document.open("text/html");
  d.writeln("<html><head>\n<title>" + title + "</title>\n</head>\n<body>");
  d.writeln("<h1>" + title + "</h1>");
  d.writeln(html_text);
  d.writeln("</body>\n</html>");
  d.close();
}

function show_properties(T, p)
{
  var p_rho;
  var html_text;
  var r = which_region(T, p);
  
  if (r)
    {
      html_text = region[r].toHTML(T, p);
      display_message("Properties for Region " + r, html_text);
    }
  else
    alert("Out of range!\n" + "T = " + T + " K")
}

var email= { username: "leon.kos",  domain : "lecad.uni-lj.si" }


</script>
</head>

<body>
<h1>IAPWS Industrial Formulation 1997</h1>

For background information about the IAPWS Industrial
Formulation 1997 for the Thermodynamic Properties of
Water and Steam (IAPWS-IF97) and references see
<a href=http://www.iapws.org/relguide/IF97.pdf>IAPWS-IF97 release </a>. 

<p>
The IAPWS-IF97 divides the thermodynamic surface
into five regions (see figure below): 

<ul>
<li> region 1 for the liquid state from low to high  pressures, 
<li> region 2 for the vapor and ideal gas state, 
<li> region 3 for the thermodynamic state around the critical point, 
<li> region 4 for the saturation curve (vapor-liquid equilibrium), 
<li> region 5 for high temperatures above 1073.15 K (800 °C) and
                pressures up to 10 MPa (100 bar). 
</ul>

<img src=regions.png alt="IAPWS 97 regions">

<p>
<h2> Properties of Water and Steam calculator</h3>
<form>
    Temperature : 
    <input type=text name=T size=10 value=300>
    <select name=T_unit>
      <option value=0> K
      <option value=273.15> &#176;C
    </select><br>
    Pressure : 
    <input type=text name=p size=10 value=3> 
    <select name=p_unit>
      <option value=1e6> MPa
      <option value=1e5> bar
      <option value=1> Pa
    </select>

  <input type=button value="Calculate properties"
     onClick="show_properties((T.value - 0) + (T_unit.options[T_unit.selectedIndex].value - 0), p.value * p_unit.options[p_unit.selectedIndex].value);">
    </form>

    <h2> Program description</h2>

Program in Javascript is designed to provide a fremework for further
development. Only properties from basic equations are covered for regions
1, 2 and 3. Backward equations are not implemented. For Region 3 density
must be calculated from pressure with Newton iteration. All regions give
correct verification results.

    <p>It seems that not all users realise the kind of software which is
described here. When you calculate properties the program runs on your
machine! Please use <kbd>view Page source</kbd> to review the source code
or click <a href=listing.html>this link</a> for fancy listing. This is
NOT external program calculator! That's why it is licensed by the following
licence. </p>
    
<h3> Licence </h3>
Copyright  2001 Leon Kos, University of Ljubljana
<p>
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License , or
    (at your option) any later version.
<p>
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
<p>
    The software is released under GNU General Public License.

<address>
    <script language=javascript> 
      document.writeln("<a href=mailto:" 
	+ email.username + "@" + email.domain + ">"
	+ " Leon Kos </a> <br> \n");
    </script>  
    University of Ljubljana, Slovenia 
</address>

<hr>
<!-- hhmts start -->Last modified: Mon Mar 14 14:26:22 CET 2005 <!-- hhmts end -->
</body> </html>

syntax highlighted by Code2HTML, v. 0.9