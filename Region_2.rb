# Ruby Document
# Constants
R = 0.461526  #kJ kg-1 K-1
T_c = 647.096 #K
p_c = 22.064  #MPa
rho_c = 322   #kg m-3

#T input
#p input
#region 2

No = [ -0.96927686500217E+01,  0.10086655968018E+02, -0.56087911283020E-02,
   0.71452738081455E-01, -0.40710498223928E+00,  0.14240819171444E+01, 
  -0.43839511319450E+01, -0.28408632460772E+00,  0.21268463753307E-01]
Jo = [0, 1, -5, -4, -3, -2, -1, 2,  3]

N = [ -0.17731742473213E-02,  -0.17834862292358E-01,  -0.45996013696365E-01, 
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
  -0.94369707241210E-06 ]
J = [0, 1, 2, 3, 6, 1, 2, 4, 7, 36, 0, 1, 3, 6, 35, 1, 2, 3, 7, 3, 16, 35, 0, 11,
25, 8, 36, 13, 4, 10, 14, 29, 50, 57, 20, 35, 48, 21, 53, 39, 26, 40, 58 ]
I = [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 5, 6, 6, 6, 7, 7, 7,
  8, 8, 9, 10, 10, 10, 16, 16, 18, 20, 20, 20, 21, 22, 23, 24, 24, 24 ]
  
def Region_2(temp,press)
  temp, press = temp.to_f, press.to_f
  pi = press/1.0
  tau = 540 / temp
  gamma_o = Math.log(pi)
  gamma_o_pi = 1.0/pi
  gamma_o_pi_pi = -1.0/pi**2.0
  gamma_o_tau=0
  gamma_o_tau_tau=0
  gamma_o_pi_tau = 0
  for i in 0..8
    gamma_o += No[i] * tau ** Jo[i]
    gamma_o_tau += No[i] * Jo[i] * tau**(Jo[i]-1)
    gamma_o_tau_tau += No[i] * Jo[i] * (Jo[i]-1) * tau**(Jo[i]-2)
  end
  gamma_r,gamma_r_pi, gamma_r_pi_pi,gamma_r_tau,gamma_r_tau_tau,gamma_r_pi_tau = 0,0,0,0,0,0
  for i in 0..42
    gamma_r += N[i] * pi ** I[i] * (tau-0.5)**J[i]
    gamma_r_pi += N[i] * I[i] * pi**(I[i]-1) * (tau-0.5)**J[i]
    gamma_r_pi_pi += N[i] * I[i] * (I[i]-1) * pi**(I[i]-2)*(tau-0.5)**J[i]
    gamma_r_tau += N[i] * pi ** I[i] * J[i] * (tau-0.5)**(J[i]-1)
    gamma_r_tau_tau += N[i] * pi**I[i]*J[i]*(J[i]-1) * (tau-0.5)**(J[i]-2)
    gamma_r_pi_tau += N[i] * I[i] * pi**(I[i]-1) * J[i] * (tau-0.5)**(J[i]-1)
  end

  
  v = pi*(gamma_o_pi+gamma_r_pi)*R*temp/press
  u = R*temp*( tau* (gamma_o_tau-gamma_r_tau) - pi*(gamma_o_pi+gamma_r_pi))
  s = R*tau*(gamma_o_tau + gamma_r_tau) - R*(gamma_o +gamma_r)
  h = tau*(gamma_o_tau + gamma_r_tau)*R*temp
  cp= -tau**2*(gamma_o_tau_tau + gamma_r_tau_tau)*R
  cv= -tau**2*(gamma_o_tau_tau + gamma_r_tau_tau)*R - R*((1+pi*gamma_r_pi-tau*pi*gamma_r_pi_tau)**2/(1-pi**2*gamma_r_pi_pi))
  w_num = (1+2*pi*gamma_r_pi+pi**2*gamma_r_pi**2)
  w_den = (1-pi**2*gamma_r_pi_pi)+((1+pi*gamma_r_pi-tau*pi*gamma_r_pi_tau)**2)/(tau**2*(gamma_o_tau_tau + gamma_r_tau_tau))
  w = Math.sqrt(R*temp*w_num/w_den)
  # w calculation is still wrong for some reason. 
  # Cv is probably wrong too
  region_2 = [v,u,s,h,cp,cv,w]
  names = ["v","u","s","h","cp","cv","w"] 
  for j in 0..region_2.length-1
    puts names[j] + " = " + region_2[j].to_s
  end
  
end

Region_2(300,0.0035)
