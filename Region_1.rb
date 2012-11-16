# Ruby Document
# Constants
R = 0.461526  #kJ kg-1 K-1
T_c = 647.096 #K
p_c = 22.064  #MPa
rho_c = 322   #kg m-3

#T input
#p input
#region 1
I = 
[
    0,0,0,0,0,0,0,0,
    1,1,1,1,1,1,2,2,
    2,2,2,3,3,3,4,4,
    4,5,8,8,21,23,29,30,
    31,32
]
J = [-2, -1, 0, 1, 2, 3, 4, 5, -9, -7, -1, 0, 1, 3, -3, 0, 1, 3, 17,
  -4, 0, 6, -5, -2, 10, -8, -11, -6, -29, -31, -38, -39, -40, -41]
N = [
     0.14632971213167e1  ,-0.84548187169114e1  ,-0.37563603672040e1  ,
     0.33855169168385e1  ,-0.95791963387872e1  , 0.15772038513228e1  ,
    -0.16616417199501e-1 , 0.81214629983568e-3 , 0.28319080123804e-3 ,
    -0.60706301565874e-3 ,-0.18990068218419e-1 ,-0.32529748770505e-1 ,
    -0.21841717175414e-1 ,-0.52838357969930e-4 ,-0.47184321073267e-3 ,
    -0.30001780793026e-3 , 0.47661393906987e-4 ,-0.44141845330846e-5 ,
    -0.72694996297594e-15,-0.31679644845054e-4 ,-0.28270797985312e-5 ,
    -0.85205128120103e-9 ,-0.22425281908000e-5 ,-0.65171222895601e-6 ,
    -0.14341729937924e-12,-0.40516996860117e-6 ,-0.12734301741641e-8 ,
    -0.17424871230634e-9 ,-0.68762131295531e-18, 0.14478307828521e-19,
     0.26335781662795e-22,-0.11947622640071e-22, 0.18228094581404e-23,
    -0.93537087292458e-25
    ]

$Eq_7_const = [I,J,n] #Call constant by Eq_7_const[i][1] to reference item within nested array

def Eq_7 (temp,press)
    gamma = 0
    i = 0
    pi = press/16.53
    tau = temp/1386
    until i == 34 do
        gamma = N[i]*(7.1-pi)**I[i] * (tau-1.222)**J[i]
        gamma_pi = -N[i]*I[i]*((7.1-pi)**(I[i]-1))*(tau-1.222)**J[i]
        gamma_pi_pi = N[i]*I[i]*(I[i]-pi)**(I[i]-2)*(tau-1.222)**J[i]
        gamma_tau = N[i]*(7.1-pi)**I[i]*J[i]*(tau-1.222)**(J[i]-1)
        gamma_tau_tau = N[i]*(7.1-pi)**I[i]*J[i]*(J[i]-1)*(tau-1.222)**(J[i]-1)
        gamma_pi_tau = -N[i]*I[i]*(7.1-pi)**(I[i]-1)*J[i]*(tau-1.222)**(J[i]-1)
        i += 1
    end
    puts gamma_pi
    v = pi*gamma_pi*R*temp
    u = R*temp*(tau*gamma_tau - pi*gamma_pi)
    s = R*(tau*gamma_tau - gamma)
    h = R*temp*tau*gamma_tau
    cp = -R*(tau**2)*gamma_tau_tau
    cv = -R*(tau**2)*gamma_tau_tau + R*(gamma_tau - tau*gamma_pi_tau)**2/gamma_pi_pi
    w = Math.sqrt(gamma_pi**2/ ((gamma_pi-tau*gamma_pi_tau)**2/(tau**2*gamma_pi_pi)-gamma_pi_pi))
    
    return region_1 = [v,u,s,h,cp,cv,w]
end

Eq_7 (300,3)