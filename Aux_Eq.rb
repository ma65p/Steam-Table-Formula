# Ruby Document
# Constants
R = 0.461526  #kJ kg-1 K-1
T_c = 647.096 #K
p_c = 22.064  #MPa
rho_c = 322   #kg m-3

#T input
#p input

#Equation 5


def Eq_5_pi (temp)
n = [0.34805185628969e3, -0.11671859879975e1, 0.10192970039326e-2, 0.57254459862746e3, 0.13918839778870e2]
    theta = temp
	pi = n[0] + n[1]*theta + n[2]*theta**2
end

def Eq_5_theta (press)
n = [0.34805185628969e3, -0.11671859879975e1, 0.10192970039326e-2, 0.57254459862746e3, 0.13918839778870e2]
    pi = press
	theta = n[3] + ((pi-n[4])/n[2])**0.5
end

test_temp = Eq_5_theta(0.165291643e2)
test_press = Eq_5_pi (0.623150000e3)

