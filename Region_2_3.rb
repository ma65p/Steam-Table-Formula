# Ruby Document
# Constants
R = 0.461526  #kJ kg-1 K-1
T_c = 647.096 #K
p_c = 22.064  #MPa
rho_c = 322   #kg m-3

#T input
#p input
#Auxillary region between 2 and 3

N = [ 0.34805185628969E+03, -0.11671859879975E+01,
      0.10192970039326E-02,  0.57254459862746E+03, 0.13918839778870E+02]
def aux_temp(temp)
  temp = temp.to_f
  press_B = 1e6 * (N[0] + N[1] * temp + N[2] * temp * temp)
end 

def aux_press(press)
  press = press.to_f
  temp_B = N[3] + Math.sqrt((press/1e6 - N[4])/N[2])
end

puts aux_temp(623.150)