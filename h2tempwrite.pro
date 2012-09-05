; Write the calculated warm H2 gas temps (err) and masses (err) to a sav file


ohmh2_temp = float([353, 0, 0,309,290,340,0,0,340,356,459,655,0,340,270,465,289,340,367,288,290,0])
ohmh2_temp_err = float([13,0,0,0,749,37,0,0,37,0,8,0,0,37,5,18,14,37,15,10,16,0])

ohmh2_mass = [3.41,0,0,0,1.42,1.42,0,0,1.42,0,0.89,0,0,1.02,2.74,0.38,1.68,2.76,2.61,5.40,2.32,0]*1d7
ohmh2_mass_err = [0.73,0,0,0,28.77,0.86,0,0,1.04,0,0.06,0,0,0.59,0.44,0.06,0.49,1.62,0.48,1.22,0.72,0]*1d7

save,ohmh2_temp,ohmh2_temp_err,ohmh2_mass,ohmh2_mass_err,filename='~/Astronomy/Research/Spitzer/OHM/ohmh2temp.sav'


conh2_temp = float([419,297,322,0,409,430,349,0,373,300,296,299])
conh2_temp_err = float([14,12,36,0,22,10,58,0,9,17,6,12])

conh2_mass = [0.81,3.47,1.36,0,0.87,0.36,1.02,0,0.69,3.95,2.40,1.03]
conh2_mass_err = [0.20,1.01,0.73,0.0,0.17,0.05,1.27,0,0.11,1.62,0.35,0.25]

save,conh2_temp,conh2_temp_err,conh2_mass,conh2_mass_err,filename='~/Astronomy/Research/Spitzer/control/conh2temp.sav'

end
