############################################################################################################################
#############################################################################################################################
lw        = 5#pre_mask_lw
LINES_STK = ([],[],[],[],[],[],[],[],[],[],[],[])
###############CENTER#####################WIDTH-FIT##################WIDTH-PLT##############CTR-FIT-BNDS###################CTR-OFFSET#############LNE_AMP_BNDS############LNE_NME_HDR_CMT################
LINES_STK[0].append( 912.0),LINES_STK[1].append(lw) ,LINES_STK[2].append(15),LINES_STK[7].append(0.002),LINES_STK[8].append(00.),LINES_STK[10].append(0.001),LINES_STK[3].append('LyL')     #Continuum Break
LINES_STK[0].append( 972.0),LINES_STK[1].append(3.0),LINES_STK[2].append(15),LINES_STK[7].append(0.002),LINES_STK[8].append(00.),LINES_STK[10].append(0.001),LINES_STK[3].append('Lyg')     #HI absoprtion
LINES_STK[0].append(1025.2),LINES_STK[1].append(lw) ,LINES_STK[2].append(15),LINES_STK[7].append(0.002),LINES_STK[8].append(00.),LINES_STK[10].append(0.001),LINES_STK[3].append('LyB')     #HI absoprtion
LINES_STK[0].append(1192.0),LINES_STK[1].append(lw) ,LINES_STK[2].append(15),LINES_STK[7].append(0.002),LINES_STK[8].append(00.),LINES_STK[10].append(0.001),LINES_STK[3].append('Si2-1')   #ISM, blend 1190+1193 
LINES_STK[0].append(1215.7),LINES_STK[1].append(5.0),LINES_STK[2].append(20),LINES_STK[7].append(0.020),LINES_STK[8].append(00.),LINES_STK[10].append(0.001),LINES_STK[3].append('Lya')    #HI emission/absorption
LINES_STK[0].append(1238.0),LINES_STK[1].append(6.0),LINES_STK[2].append(14),LINES_STK[7].append(0.020),LINES_STK[8].append(00.),LINES_STK[10].append(0.001),LINES_STK[3].append('N5-1')   #
LINES_STK[0].append(1242.0),LINES_STK[1].append(6.0),LINES_STK[2].append(14),LINES_STK[7].append(0.020),LINES_STK[8].append(-2.),LINES_STK[10].append(0.001),LINES_STK[3].append('N5-2')   #
LINES_STK[0].append(1260.4),LINES_STK[1].append(4.0),LINES_STK[2].append(18),LINES_STK[7].append(0.020),LINES_STK[8].append(00.),LINES_STK[10].append(0.001),LINES_STK[3].append('Si2-1')   #ISM Si II 1263 triplet 
LINES_STK[0].append(1264.0),LINES_STK[1].append(4.0),LINES_STK[2].append(18),LINES_STK[7].append(0.020),LINES_STK[8].append(00.),LINES_STK[10].append(0.001),LINES_STK[3].append('Si2-*')   #SiII** Nebular Shapley
LINES_STK[0].append(1303.2),LINES_STK[1].append(5.0),LINES_STK[2].append(19),LINES_STK[7].append(0.020),LINES_STK[8].append(03.),LINES_STK[10].append(0.001),LINES_STK[3].append('OSD')    #ISM, blend 1303 OI-SiIII
LINES_STK[0].append(1309.0),LINES_STK[1].append(4.0),LINES_STK[2].append(18),LINES_STK[7].append(0.020),LINES_STK[8].append(00.),LINES_STK[10].append(0.001),LINES_STK[3].append('Si2-*')   #SiII** Nebular Shapley
LINES_STK[0].append(1334.5),LINES_STK[1].append(3.0),LINES_STK[2].append(19),LINES_STK[7].append(0.002),LINES_STK[8].append(04.),LINES_STK[10].append(0.001),LINES_STK[3].append('C2-1')   #ISM Triplet
LINES_STK[0].append(1343.0),LINES_STK[1].append(2.0),LINES_STK[2].append(19),LINES_STK[7].append(0.022),LINES_STK[8].append(04.),LINES_STK[10].append(0.001),LINES_STK[3].append('O4')    #Wind Line http://www.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php ****O stars
LINES_STK[0].append(1371.0),LINES_STK[1].append(4.0),LINES_STK[2].append(19),LINES_STK[7].append(0.022),LINES_STK[8].append(04.),LINES_STK[10].append(0.001),LINES_STK[3].append('O5')    #New-Update Wind Line #1371 OV OTELO DATABASE http://research.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
LINES_STK[0].append(1393.8),LINES_STK[1].append(3.0),LINES_STK[2].append(14),LINES_STK[7].append(0.001),LINES_STK[8].append(04.),LINES_STK[10].append(0.001),LINES_STK[3].append('Si4-1')  #ISM
LINES_STK[0].append(1402.8),LINES_STK[1].append(3.0),LINES_STK[2].append(16),LINES_STK[7].append(0.002),LINES_STK[8].append(-2,),LINES_STK[10].append(0.001),LINES_STK[3].append('Si4-2')  #ISM
LINES_STK[0].append(1417.2),LINES_STK[1].append(3.0),LINES_STK[2].append(16),LINES_STK[7].append(0.020),LINES_STK[8].append(00.),LINES_STK[10].append(0.001),LINES_STK[3].append('Si3-3')  #Halliday+08GMASS photospheric abs lin ERROR IN ID AS Si2 sholuld be Si3   ****O stars
LINES_STK[0].append(1453.0),LINES_STK[1].append(3.0),LINES_STK[2].append(16),LINES_STK[7].append(0.020),LINES_STK[8].append(05.),LINES_STK[10].append(0.001),LINES_STK[3].append('Fe5')    #http://research.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
LINES_STK[0].append(1485.4),LINES_STK[1].append(4.0),LINES_STK[2].append(18),LINES_STK[7].append(0.020),LINES_STK[8].append(00.),LINES_STK[10].append(0.001),LINES_STK[3].append('Si2-4')  #Halliday+08GMASS photospheric abs lin
LINES_STK[0].append(1486.0),LINES_STK[1].append(4.0),LINES_STK[2].append(16),LINES_STK[7].append(0.020),LINES_STK[8].append(00.),LINES_STK[10].append(0.001),LINES_STK[3].append('N4-1')   #http://research.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
LINES_STK[0].append(1501.8),LINES_STK[1].append(4.0),LINES_STK[2].append(16),LINES_STK[7].append(0.002),LINES_STK[8].append(00.),LINES_STK[10].append(0.001),LINES_STK[3].append('S5')     #Halliday+08GMASS photospheric abs lin ****O stars
LINES_STK[0].append(1526.7),LINES_STK[1].append(5.0),LINES_STK[2].append(19),LINES_STK[7].append(0.002),LINES_STK[8].append(00.),LINES_STK[10].append(0.001),LINES_STK[3].append('Si2-5')  #ISM photospheric abs lin
LINES_STK[0].append(1548.2),LINES_STK[1].append(4.0),LINES_STK[2].append(18),LINES_STK[7].append(0.002),LINES_STK[8].append(00.),LINES_STK[10].append(0.001),LINES_STK[3].append('C4-1')   #ISM Blend 1548.2 + 1550.8 1
LINES_STK[0].append(1550.7),LINES_STK[1].append(4.0),LINES_STK[2].append(18),LINES_STK[7].append(0.002),LINES_STK[8].append(00.),LINES_STK[10].append(0.001),LINES_STK[3].append('C4-2')   #ISM Blend 1548.2 + 1550.8 2
LINES_STK[0].append(1608.5),LINES_STK[1].append(4.0),LINES_STK[2].append(20),LINES_STK[7].append(0.002),LINES_STK[8].append(00.),LINES_STK[10].append(0.001),LINES_STK[3].append('Fe2-1')  #ISM
LINES_STK[0].append(1640.0),LINES_STK[1].append(4.0),LINES_STK[2].append(18),LINES_STK[7].append(0.020),LINES_STK[8].append(-5.),LINES_STK[10].append(0.001),LINES_STK[3].append('He2')    #Nebular
LINES_STK[0].append(1656.0),LINES_STK[1].append(3.0),LINES_STK[2].append(18),LINES_STK[7].append(0.020),LINES_STK[8].append(-5.),LINES_STK[10].append(0.001),LINES_STK[3].append('C1')     #New-Update Fanelli #1656 CI OTELO DATABASE http://research.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
LINES_STK[0].append(1660.0),LINES_STK[1].append(3.0),LINES_STK[2].append(18),LINES_STK[7].append(0.020),LINES_STK[8].append(-5.),LINES_STK[10].append(0.001),LINES_STK[3].append('O3]-1')   #New-Update Fanelli #1660 OII-Doublet OTELO DATABASE http://research.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
LINES_STK[0].append(1666.0),LINES_STK[1].append(3.0),LINES_STK[2].append(18),LINES_STK[7].append(0.020),LINES_STK[8].append(-5.),LINES_STK[10].append(0.001),LINES_STK[3].append('O3]-2')   #New-Update Fanelli #1666 OII-Doublet OTELO DATABASE http://research.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
LINES_STK[0].append(1670.8),LINES_STK[1].append(4.0),LINES_STK[2].append(18),LINES_STK[7].append(0.020),LINES_STK[8].append(00.),LINES_STK[10].append(0.001),LINES_STK[3].append('Al2')    #ISM
LINES_STK[0].append(1709.6),LINES_STK[1].append(3.0),LINES_STK[2].append(18),LINES_STK[7].append(0.002),LINES_STK[8].append(00.),LINES_STK[10].append(0.001),LINES_STK[3].append('Ni2-1')  #Halliday+08GMASS 
LINES_STK[0].append(1718.5),LINES_STK[1].append(2.5),LINES_STK[2].append(18),LINES_STK[7].append(0.020),LINES_STK[8].append(-6.),LINES_STK[10].append(0.001),LINES_STK[3].append('N4-2')   #Halliday+08GMASS photospheric abs lin
LINES_STK[0].append(1722.0),LINES_STK[1].append(3.0),LINES_STK[2].append(18),LINES_STK[7].append(0.020),LINES_STK[8].append(+5.),LINES_STK[10].append(0.001),LINES_STK[3].append('Si4-3')  #Fanelli BL 1719 http://www.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
LINES_STK[0].append(1741.5),LINES_STK[1].append(4.0),LINES_STK[2].append(18),LINES_STK[7].append(0.020),LINES_STK[8].append(-6.),LINES_STK[10].append(0.001),LINES_STK[3].append('Ni2-2')  #Halliday+08GMASS 
LINES_STK[0].append(1751.9),LINES_STK[1].append(2.0),LINES_STK[2].append(18),LINES_STK[7].append(0.020),LINES_STK[8].append(-2.),LINES_STK[10].append(0.001),LINES_STK[3].append('Ni2-3')  #Halliday+08GMASS 
LINES_STK[0].append(1808.0),LINES_STK[1].append(3.0),LINES_STK[2].append(18),LINES_STK[7].append(0.020),LINES_STK[8].append(+2.),LINES_STK[10].append(0.001),LINES_STK[3].append('Si2-6')  #Halliday+08GMASS 
LINES_STK[0].append(1854.7),LINES_STK[1].append(4.0),LINES_STK[2].append(17),LINES_STK[7].append(0.020),LINES_STK[8].append(-2.),LINES_STK[10].append(0.001),LINES_STK[3].append('Al3-1')  #ISM H+? Fe2-2 Old
LINES_STK[0].append(1862.8),LINES_STK[1].append(4.0),LINES_STK[2].append(17),LINES_STK[7].append(0.020),LINES_STK[8].append(+2.),LINES_STK[10].append(0.001),LINES_STK[3].append('Al3-2')  #ISM H+? Fe2-3 Old
LINES_STK[0].append(1889.0),LINES_STK[1].append(4.0),LINES_STK[2].append(17),LINES_STK[7].append(0.020),LINES_STK[8].append(+2.),LINES_STK[10].append(0.001),LINES_STK[3].append('Si3]')  #New-Update ISM H+? #1889 SiIII OTELO DATABASE http://research.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
LINES_STK[0].append(1907.0),LINES_STK[1].append(4.0),LINES_STK[2].append(15),LINES_STK[7].append(0.020),LINES_STK[8].append(-10),LINES_STK[10].append(0.001),LINES_STK[3].append('C3]-1')   #Nebular Blend 1907 + 1909
LINES_STK[0].append(1909.0),LINES_STK[1].append(4.0),LINES_STK[2].append(15),LINES_STK[7].append(0.020),LINES_STK[8].append(+15),LINES_STK[10].append(0.001),LINES_STK[3].append('C3]-2')   #Nebular Blend 1907 + 1909
LINES_STK[0].append(2028.4),LINES_STK[1].append(4.0),LINES_STK[2].append(15),LINES_STK[7].append(0.022),LINES_STK[8].append(00.),LINES_STK[10].append(0.001),LINES_STK[3].append('Mg1-1')  #Halliday+08GMASS 
LINES_STK[0].append(2062.6),LINES_STK[1].append(3.0),LINES_STK[2].append(15),LINES_STK[7].append(0.020),LINES_STK[8].append(+2.),LINES_STK[10].append(0.001),LINES_STK[3].append('Zn2')    #Halliday+08GMASS 
LINES_STK[0].append(2325.5),LINES_STK[1].append(8.0),LINES_STK[2].append(15),LINES_STK[7].append(0.020),LINES_STK[8].append(00.),LINES_STK[10].append(0.001),LINES_STK[3].append('C2]')   #http://www.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
LINES_STK[0].append(2343.5),LINES_STK[1].append(3.0),LINES_STK[2].append(15),LINES_STK[7].append(0.020),LINES_STK[8].append(00.),LINES_STK[10].append(0.001),LINES_STK[3].append('Fe2-4')  #ISM
LINES_STK[0].append(2370.5),LINES_STK[1].append(4.0),LINES_STK[2].append(15),LINES_STK[7].append(0.020),LINES_STK[8].append(00.),LINES_STK[10].append(0.001),LINES_STK[3].append('Fe2-*')  #ISM
LINES_STK[0].append(2402.6),LINES_STK[1].append(3.0),LINES_STK[2].append(15),LINES_STK[7].append(0.020),LINES_STK[8].append(00.),LINES_STK[10].append(0.001),LINES_STK[3].append('Fe2-6')  #ISM
LINES_STK[0].append(2608.0),LINES_STK[1].append(5.0),LINES_STK[2].append(10),LINES_STK[7].append(0.020),LINES_STK[8].append(-2.),LINES_STK[10].append(0.001),LINES_STK[3].append('Fe2-7')  #ISM previouus 2593.7 2608CORRECTED OTHELO
LINES_STK[0].append(2796.0),LINES_STK[1].append(6.0),LINES_STK[2].append(15),LINES_STK[7].append(0.200),LINES_STK[8].append(+4.),LINES_STK[10].append(0.001),LINES_STK[3].append('Mg2-1')  #ISM
LINES_STK[0].append(2803.0),LINES_STK[1].append(4.0),LINES_STK[2].append(15),LINES_STK[7].append(0.200),LINES_STK[8].append(-4.),LINES_STK[10].append(0.001),LINES_STK[3].append('Mg2-2')  #ISM
LINES_STK[0].append(2852.0),LINES_STK[1].append(3.0),LINES_STK[2].append(15),LINES_STK[7].append(0.020),LINES_STK[8].append(+2.),LINES_STK[10].append(0.001),LINES_STK[3].append('Mg1-2')  #ISM
###############CENTER#####################WIDTH-FIT##################WIDTH-PLT##############CTR-FIT-BNDS###################CTR-OFFSET#############LNE_AMP_BNDS############LNE_NME_HDR_CMT################

###########MARKER##################LNE-PLT-TXT######################LNE-HDR-NME###################LNE-LBL-NME#########
LINES_STK[6].append('X'),LINES_STK[4].append('Ly-L')    ,LINES_STK[5].append('L01'),LINES_STK[9].append('Ly$_{limit}$')     ##Continuum Break
LINES_STK[6].append('X'),LINES_STK[4].append('Ly-g')    ,LINES_STK[5].append('L02'),LINES_STK[9].append('Ly$_{\\gamma}$')     #HI absoprtion
LINES_STK[6].append('o'),LINES_STK[4].append('Ly-B')    ,LINES_STK[5].append('L03'),LINES_STK[9].append('Ly$_{\\beta}$')     #HI absoprtion
LINES_STK[6].append('v'),LINES_STK[4].append('SiII-1')  ,LINES_STK[5].append('S04'),LINES_STK[9].append('SiII-1')   #ISM, blend 1190+1193 
LINES_STK[6].append('^'),LINES_STK[4].append('Ly-a')    ,LINES_STK[5].append('L05'),LINES_STK[9].append('Ly$_{\\alpha}$')      #HI emission/absorption
LINES_STK[6].append('<'),LINES_STK[4].append('NV-1')    ,LINES_STK[5].append('N06'),LINES_STK[9].append('NV')        #
LINES_STK[6].append('<'),LINES_STK[4].append('NV-2')    ,LINES_STK[5].append('N07'),LINES_STK[9].append('NV')        #
LINES_STK[6].append('>'),LINES_STK[4].append('SiII')    ,LINES_STK[5].append('S09'),LINES_STK[9].append('SiII')       #ISM
LINES_STK[6].append('>'),LINES_STK[4].append('SiII*')   ,LINES_STK[5].append('S10'),LINES_STK[9].append('SiII*')      #SiII** Nebular Shapley
LINES_STK[6].append('8'),LINES_STK[4].append('OI+SiII') ,LINES_STK[5].append('O12'),LINES_STK[9].append('OI+SiII')   #ISM, blend 1303
LINES_STK[6].append('>'),LINES_STK[4].append('SiII*')   ,LINES_STK[5].append('S11'),LINES_STK[9].append('SiII*')      #SiII** Nebular Shapley
LINES_STK[6].append('s'),LINES_STK[4].append('CII-1')   ,LINES_STK[5].append('C13'),LINES_STK[9].append('CII')       #ISM
LINES_STK[6].append('s'),LINES_STK[4].append('OIV')     ,LINES_STK[5].append('O14'),LINES_STK[9].append('OIV')       #Wind Line http://www.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
LINES_STK[6].append('^'),LINES_STK[4].append('OV')      ,LINES_STK[5].append('O90'),LINES_STK[9].append('OV')       #New-Update Wind Line #1371 OV OTELO DATABASE http://research.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
LINES_STK[6].append('*'),LINES_STK[4].append('SiIV-1')  ,LINES_STK[5].append('S15'),LINES_STK[9].append('SiIV')      #ISM
LINES_STK[6].append('p'),LINES_STK[4].append('SiIV-2')  ,LINES_STK[5].append('S16'),LINES_STK[9].append('SiIV')      #ISM
LINES_STK[6].append('s'),LINES_STK[4].append('SiIII-3') ,LINES_STK[5].append('S18'),LINES_STK[9].append('SiIII')      #Halliday+08GMASS photospheric abs lin
LINES_STK[6].append('p'),LINES_STK[4].append('FeV')     ,LINES_STK[5].append('F19'),LINES_STK[9].append('FeV')       #http://research.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
LINES_STK[6].append('P'),LINES_STK[4].append('SiII-4')  ,LINES_STK[5].append('S20'),LINES_STK[9].append('SiII')      #Halliday+08GMASS photospheric abs lin
LINES_STK[6].append('*'),LINES_STK[4].append('NIV-1')   ,LINES_STK[5].append('N21'),LINES_STK[9].append('NIV')       #http://research.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
LINES_STK[6].append('h'),LINES_STK[4].append('SV')      ,LINES_STK[5].append('S22'),LINES_STK[9].append('SV')        #Halliday+08GMASS photospheric abs lin
LINES_STK[6].append('H'),LINES_STK[4].append('SiII-5')  ,LINES_STK[5].append('S23'),LINES_STK[9].append('SiII')      #ISM photospheric abs lin
LINES_STK[6].append('X'),LINES_STK[4].append('CIV-1')   ,LINES_STK[5].append('C24'),LINES_STK[9].append('CIV')       #ISM Blend 1548.2 + 1550.8 1
LINES_STK[6].append('P'),LINES_STK[4].append('CIV-2')   ,LINES_STK[5].append('C25'),LINES_STK[9].append('CIV')       #ISM Blend 1548.2 + 1550.8 2
LINES_STK[6].append('D'),LINES_STK[4].append('FeII-1')  ,LINES_STK[5].append('F27'),LINES_STK[9].append('FeII')      #ISM
LINES_STK[6].append('d'),LINES_STK[4].append('HeII')    ,LINES_STK[5].append('H28'),LINES_STK[9].append('HeII')      #Nebular
LINES_STK[6].append('H'),LINES_STK[4].append('CI')      ,LINES_STK[5].append('C91'),LINES_STK[9].append('CI')              #New-Update Fanelli #1656 CI OTELO DATABASE http://research.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
LINES_STK[6].append('X'),LINES_STK[4].append('OII]-1')  ,LINES_STK[5].append('O92'),LINES_STK[9].append('OIII]')          #New-Update Fanelli #1660 OII-Doublet OTELO DATABASE http://research.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
LINES_STK[6].append('P'),LINES_STK[4].append('OII]-2')  ,LINES_STK[5].append('O93'),LINES_STK[9].append('OIII]')          #New-Update Fanelli #1666 OII-Doublet OTELO DATABASE http://research.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
LINES_STK[6].append('o'),LINES_STK[4].append('AlII')    ,LINES_STK[5].append('A29'),LINES_STK[9].append('AlII')      #ISM
LINES_STK[6].append('X'),LINES_STK[4].append('NiII-1')  ,LINES_STK[5].append('N30'),LINES_STK[9].append('NiII')      #Halliday+08GMASS 
LINES_STK[6].append('o'),LINES_STK[4].append('NIV-2')   ,LINES_STK[5].append('N31'),LINES_STK[9].append('NIV')       #Halliday+08GMASS photospheric abs lin
LINES_STK[6].append('o'),LINES_STK[4].append('SiIV-3')  ,LINES_STK[5].append('S32'),LINES_STK[9].append('SIV')       #Fanelli BL 1719 http://www.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
LINES_STK[6].append('v'),LINES_STK[4].append('NiII-2')  ,LINES_STK[5].append('N33'),LINES_STK[9].append('NiII')      #Halliday+08GMASS 
LINES_STK[6].append('^'),LINES_STK[4].append('NiII-3')  ,LINES_STK[5].append('N34'),LINES_STK[9].append('NiII')      #Halliday+08GMASS 
LINES_STK[6].append('>'),LINES_STK[4].append('SiII-6')  ,LINES_STK[5].append('S36'),LINES_STK[9].append('SiII')      #Halliday+08GMASS 
LINES_STK[6].append('8'),LINES_STK[4].append('AlIII-1') ,LINES_STK[5].append('A37'),LINES_STK[9].append('AlIII')      #ISM H+? Old FeII
LINES_STK[6].append('s'),LINES_STK[4].append('AlIII-2') ,LINES_STK[5].append('A38'),LINES_STK[9].append('AlIII')      #ISM H+? Old FeII
LINES_STK[6].append('s'),LINES_STK[4].append('SiIII]')  ,LINES_STK[5].append('F94'),LINES_STK[9].append('SiIII]')     #New-Update ISM H+? #1889 SiIII OTELO DATABASE http://research.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
LINES_STK[6].append('*'),LINES_STK[4].append('CIII]-1') ,LINES_STK[5].append('C39'),LINES_STK[9].append('CIII]')      #Nebular Blend 1907 + 1909
LINES_STK[6].append('*'),LINES_STK[4].append('CIII]-2') ,LINES_STK[5].append('C40'),LINES_STK[9].append('CIII]')      #Nebular Blend 1907 + 1909
LINES_STK[6].append('p'),LINES_STK[4].append('MgI')     ,LINES_STK[5].append('M41'),LINES_STK[9].append('MgI')       #Halliday+08GMASS 
LINES_STK[6].append('8'),LINES_STK[4].append('ZnII')    ,LINES_STK[5].append('Z42'),LINES_STK[9].append('ZnII')      #Halliday+08GMASS 
LINES_STK[6].append('s'),LINES_STK[4].append('CII]')    ,LINES_STK[5].append('C43'),LINES_STK[9].append('CII]')       #http://www.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
LINES_STK[6].append('p'),LINES_STK[4].append('FeII-4')  ,LINES_STK[5].append('F44'),LINES_STK[9].append('FeII')      #ISM
LINES_STK[6].append('p'),LINES_STK[4].append('FeII-5')  ,LINES_STK[5].append('F45'),LINES_STK[9].append('FeII*')     #ISM
LINES_STK[6].append('P'),LINES_STK[4].append('FeII-6')  ,LINES_STK[5].append('F46'),LINES_STK[9].append('FeII')      #ISM
LINES_STK[6].append('*'),LINES_STK[4].append('FeII-7')  ,LINES_STK[5].append('F47'),LINES_STK[9].append('FeII')      #ISM
LINES_STK[6].append('h'),LINES_STK[4].append('MgII-1')  ,LINES_STK[5].append('M48'),LINES_STK[9].append('MgII')      #ISM
LINES_STK[6].append('h'),LINES_STK[4].append('MgII-2')  ,LINES_STK[5].append('M49'),LINES_STK[9].append('MgII')      #ISM
LINES_STK[6].append('h'),LINES_STK[4].append('MgI-2')   ,LINES_STK[5].append('M50'),LINES_STK[9].append('MgI')       #ISM
###########MARKER##################LNE-PLT-TXT######################LNE-HDR-NME###################LNE-LBL-NME#########
#############################################################################################################################
#############################################################################################################################

#########################################################################################################################
#############################################################################################################################
lw        = 5#pre_mask_lw
LINES_PLT_FG = ([],[],[],[],[],[],[],[],[],[],[],[])
###############CENTER#####################WIDTH-FIT##################WIDTH-PLT##############CTR-FIT-BNDS###################CTR-OFFSET#############LNE_AMP_BNDS############LNE_NME_HDR_CMT################
#LINES_PLT_FG[0].append( 912.0),LINES_PLT_FG[1].append(lw),LINES_PLT_FG[2].append(15),LINES_PLT_FG[7].append(0.002),LINES_PLT_FG[8].append(00),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('LyL')      #Continuum Break
#LINES_PLT_FG[0].append( 972.0),LINES_PLT_FG[1].append(3.0),LINES_PLT_FG[2].append(15),LINES_PLT_FG[7].append(0.002),LINES_PLT_FG[8].append(00),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('Lyg')     #HI absoprtion
#LINES_PLT_FG[0].append(1025.2),LINES_PLT_FG[1].append(lw),LINES_PLT_FG[2].append(15),LINES_PLT_FG[7].append(0.002),LINES_PLT_FG[8].append(00.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('LyB')     #HI absoprtion
LINES_PLT_FG[0].append(1192.0),LINES_PLT_FG[1].append(lw),LINES_PLT_FG[2].append(15),LINES_PLT_FG[7].append(0.002),LINES_PLT_FG[8].append(00.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('Si2-1')   #ISM, blend 1190+1193 
LINES_PLT_FG[0].append(1215.7),LINES_PLT_FG[1].append(5.0),LINES_PLT_FG[2].append(20),LINES_PLT_FG[7].append(0.020),LINES_PLT_FG[8].append(00.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('Lya')    #HI emission/absorption
LINES_PLT_FG[0].append(1238.0),LINES_PLT_FG[1].append(6.0),LINES_PLT_FG[2].append(18),LINES_PLT_FG[7].append(0.020),LINES_PLT_FG[8].append(00.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('N5-1')   #
LINES_PLT_FG[0].append(1242.0),LINES_PLT_FG[1].append(6.0),LINES_PLT_FG[2].append(14),LINES_PLT_FG[7].append(0.020),LINES_PLT_FG[8].append(-2.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('N5-2')   #
##LINES_PLT_FG[0].append(1240.0),LINES_PLT_FG[1].append(6.0),LINES_PLT_FG[2].append(14),LINES_PLT_FG[7].append(0.020),LINES_PLT_FG[8].append(00.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('N5Dbl')  #
LINES_PLT_FG[0].append(1260.4),LINES_PLT_FG[1].append(4.0),LINES_PLT_FG[2].append(18),LINES_PLT_FG[7].append(0.020),LINES_PLT_FG[8].append(00.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('Si2-1')   #ISM Si II 1263 triplet 
LINES_PLT_FG[0].append(1264.0),LINES_PLT_FG[1].append(4.0),LINES_PLT_FG[2].append(18),LINES_PLT_FG[7].append(0.020),LINES_PLT_FG[8].append(00.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('Si2-*')   #SiII** Nebular Shapley
LINES_PLT_FG[0].append(1303.2),LINES_PLT_FG[1].append(5.0),LINES_PLT_FG[2].append(19),LINES_PLT_FG[7].append(0.020),LINES_PLT_FG[8].append(03.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('OSD')    #ISM, blend 1303 OI-SiIII
LINES_PLT_FG[0].append(1309.0),LINES_PLT_FG[1].append(4.0),LINES_PLT_FG[2].append(18),LINES_PLT_FG[7].append(0.020),LINES_PLT_FG[8].append(00.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('Si2-*')   #SiII** Nebular Shapley
LINES_PLT_FG[0].append(1334.5),LINES_PLT_FG[1].append(3.0),LINES_PLT_FG[2].append(19),LINES_PLT_FG[7].append(0.002),LINES_PLT_FG[8].append(04.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('C2-1')   #ISM Triplet
LINES_PLT_FG[0].append(1343.0),LINES_PLT_FG[1].append(4.0),LINES_PLT_FG[2].append(22),LINES_PLT_FG[7].append(0.022),LINES_PLT_FG[8].append(00.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('O4')    #Wind Line http://www.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php ****O stars
#LINES_PLT_FG[0].append(1371.0),LINES_PLT_FG[1].append(4.0),LINES_PLT_FG[2].append(19),LINES_PLT_FG[7].append(0.022),LINES_PLT_FG[8].append(04.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('O5')    #New-Update Wind Line #1371 OV OTELO DATABASE http://research.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
LINES_PLT_FG[0].append(1393.8),LINES_PLT_FG[1].append(2.0),LINES_PLT_FG[2].append(18),LINES_PLT_FG[7].append(0.001),LINES_PLT_FG[8].append(00.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('Si4-1')  #ISM
LINES_PLT_FG[0].append(1402.8),LINES_PLT_FG[1].append(2.0),LINES_PLT_FG[2].append(16),LINES_PLT_FG[7].append(0.062),LINES_PLT_FG[8].append(02,),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('Si4-2')  #ISM
#LINES_PLT_FG[0].append(1404.0),LINES_PLT_FG[1].append(3.0),LINES_PLT_FG[2].append(18),LINES_PLT_FG[7].append(0.002),LINES_PLT_FG[8].append(04.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('Si4Dbl') #ISM
LINES_PLT_FG[0].append(1417.2),LINES_PLT_FG[1].append(3.0),LINES_PLT_FG[2].append(16),LINES_PLT_FG[7].append(0.020),LINES_PLT_FG[8].append(00.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('Si3-3')  #Halliday+08GMASS photospheric abs lin ERROR IN ID AS Si2 sholuld be Si3   ****O stars
#LINES_PLT_FG[0].append(1453.0),LINES_PLT_FG[1].append(3.0),LINES_PLT_FG[2].append(16),LINES_PLT_FG[7].append(0.020),LINES_PLT_FG[8].append(05.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('Fe5')    #http://research.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
#LINES_PLT_FG[0].append(1485.4),LINES_PLT_FG[1].append(4.0),LINES_PLT_FG[2].append(18),LINES_PLT_FG[7].append(0.020),LINES_PLT_FG[8].append(00.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('Si2-4')  #Halliday+08GMASS photospheric abs lin
#LINES_PLT_FG[0].append(1486.0),LINES_PLT_FG[1].append(4.0),LINES_PLT_FG[2].append(16),LINES_PLT_FG[7].append(0.020),LINES_PLT_FG[8].append(00.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('N4-1')   #http://research.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
LINES_PLT_FG[0].append(1501.8),LINES_PLT_FG[1].append(4.0),LINES_PLT_FG[2].append(16),LINES_PLT_FG[7].append(0.002),LINES_PLT_FG[8].append(00.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('S5')     #Halliday+08GMASS photospheric abs lin ****O stars
LINES_PLT_FG[0].append(1526.7),LINES_PLT_FG[1].append(5.0),LINES_PLT_FG[2].append(19),LINES_PLT_FG[7].append(0.002),LINES_PLT_FG[8].append(00.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('Si2-5')  #ISM photospheric abs lin
LINES_PLT_FG[0].append(1548.2),LINES_PLT_FG[1].append(4.0),LINES_PLT_FG[2].append(18),LINES_PLT_FG[7].append(0.002),LINES_PLT_FG[8].append(00.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('C4-1')   #ISM Blend 1548.2 + 1550.8 1
LINES_PLT_FG[0].append(1550.7),LINES_PLT_FG[1].append(4.0),LINES_PLT_FG[2].append(18),LINES_PLT_FG[7].append(0.002),LINES_PLT_FG[8].append(00.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('C4-2')   #ISM Blend 1548.2 + 1550.8 2
##LINES_PLT_FG[0].append(1549.1),LINES_PLT_FG[1].append(4.0),LINES_PLT_FG[2].append(18),LINES_PLT_FG[7].append(0.002),LINES_PLT_FG[8].append(00.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('CDbl')   #ISM Blend 1548.2 + 1550.8
LINES_PLT_FG[0].append(1608.5),LINES_PLT_FG[1].append(4.0),LINES_PLT_FG[2].append(20),LINES_PLT_FG[7].append(0.002),LINES_PLT_FG[8].append(00.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('Fe2-1')  #ISM
LINES_PLT_FG[0].append(1640.0),LINES_PLT_FG[1].append(5.0),LINES_PLT_FG[2].append(20),LINES_PLT_FG[7].append(0.020),LINES_PLT_FG[8].append(00.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('He2')    #Nebular
#LINES_PLT_FG[0].append(1656.0),LINES_PLT_FG[1].append(2.0),LINES_PLT_FG[2].append(18),LINES_PLT_FG[7].append(0.020),LINES_PLT_FG[8].append(00.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('C1')     #New-Update Fanelli #1656 CI OTELO DATABASE http://research.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
LINES_PLT_FG[0].append(1660.0),LINES_PLT_FG[1].append(2.0),LINES_PLT_FG[2].append(18),LINES_PLT_FG[7].append(0.020),LINES_PLT_FG[8].append(00.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('O3]-1')   #New-Update Fanelli #1660 OII-Doublet OTELO DATABASE http://research.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
LINES_PLT_FG[0].append(1666.0),LINES_PLT_FG[1].append(4.0),LINES_PLT_FG[2].append(18),LINES_PLT_FG[7].append(0.020),LINES_PLT_FG[8].append(02.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('O3]-2')   #New-Update Fanelli #1666 OII-Doublet OTELO DATABASE http://research.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
LINES_PLT_FG[0].append(1670.8),LINES_PLT_FG[1].append(4.0),LINES_PLT_FG[2].append(18),LINES_PLT_FG[7].append(0.020),LINES_PLT_FG[8].append(00.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('Al2')    #ISM
LINES_PLT_FG[0].append(1709.6),LINES_PLT_FG[1].append(2.5),LINES_PLT_FG[2].append(18),LINES_PLT_FG[7].append(0.002),LINES_PLT_FG[8].append(00.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('Ni2-1')  #Halliday+08GMASS 
LINES_PLT_FG[0].append(1718.5),LINES_PLT_FG[1].append(2.5),LINES_PLT_FG[2].append(18),LINES_PLT_FG[7].append(0.020),LINES_PLT_FG[8].append(00.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('N4-2')   #Halliday+08GMASS photospheric abs lin
#LINES_PLT_FG[0].append(1722.0),LINES_PLT_FG[1].append(3.0),LINES_PLT_FG[2].append(18),LINES_PLT_FG[7].append(0.020),LINES_PLT_FG[8].append(+5.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('Si4-3')  #Fanelli BL 1719 http://www.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
LINES_PLT_FG[0].append(1741.5),LINES_PLT_FG[1].append(4.0),LINES_PLT_FG[2].append(18),LINES_PLT_FG[7].append(0.020),LINES_PLT_FG[8].append(00.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('Ni2-2')  #Halliday+08GMASS 
LINES_PLT_FG[0].append(1751.9),LINES_PLT_FG[1].append(4.0),LINES_PLT_FG[2].append(18),LINES_PLT_FG[7].append(0.020),LINES_PLT_FG[8].append(00.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('Ni2-3')  #Halliday+08GMASS 
##LINES_PLT_FG[0].append(1746.9),LINES_PLT_FG[1].append(3.0),LINES_PLT_FG[2].append(18),LINES_PLT_FG[7].append(0.002),LINES_PLT_FG[8].append(00.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('Ni2Dbl') #Halliday+08GMASS Doublet
LINES_PLT_FG[0].append(1808.0),LINES_PLT_FG[1].append(3.0),LINES_PLT_FG[2].append(18),LINES_PLT_FG[7].append(0.020),LINES_PLT_FG[8].append(+2.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('Si2-6')  #Halliday+08GMASS 
LINES_PLT_FG[0].append(1854.7),LINES_PLT_FG[1].append(2.0),LINES_PLT_FG[2].append(20),LINES_PLT_FG[7].append(0.020),LINES_PLT_FG[8].append(-2.),LINES_PLT_FG[10].append(0.091),LINES_PLT_FG[3].append('Al3-1')  #ISM H+? Fe2-2 Old
LINES_PLT_FG[0].append(1862.8),LINES_PLT_FG[1].append(4.0),LINES_PLT_FG[2].append(17),LINES_PLT_FG[7].append(0.020),LINES_PLT_FG[8].append(+2.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('Al3-2')  #ISM H+? Fe2-3 Old
LINES_PLT_FG[0].append(1889.0),LINES_PLT_FG[1].append(4.0),LINES_PLT_FG[2].append(17),LINES_PLT_FG[7].append(0.020),LINES_PLT_FG[8].append(+2.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('Si3]')  #New-Update ISM H+? #1889 SiIII OTELO DATABASE http://research.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
LINES_PLT_FG[0].append(1907.0),LINES_PLT_FG[1].append(4.0),LINES_PLT_FG[2].append(15),LINES_PLT_FG[7].append(0.020),LINES_PLT_FG[8].append(00.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('C3]-1')   #Nebular Blend 1907 + 1909
LINES_PLT_FG[0].append(1909.0),LINES_PLT_FG[1].append(4.0),LINES_PLT_FG[2].append(15),LINES_PLT_FG[7].append(0.020),LINES_PLT_FG[8].append(+15),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('C3]-2')   #Nebular Blend 1907 + 1909
LINES_PLT_FG[0].append(2028.4),LINES_PLT_FG[1].append(4.0),LINES_PLT_FG[2].append(15),LINES_PLT_FG[7].append(0.022),LINES_PLT_FG[8].append(00.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('Mg1-1')  #Halliday+08GMASS 
LINES_PLT_FG[0].append(2062.6),LINES_PLT_FG[1].append(3.0),LINES_PLT_FG[2].append(15),LINES_PLT_FG[7].append(0.020),LINES_PLT_FG[8].append(+2.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('Zn2')    #Halliday+08GMASS 
LINES_PLT_FG[0].append(2325.5),LINES_PLT_FG[1].append(8.0),LINES_PLT_FG[2].append(15),LINES_PLT_FG[7].append(0.020),LINES_PLT_FG[8].append(00.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('C2]')   #http://www.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
LINES_PLT_FG[0].append(2343.5),LINES_PLT_FG[1].append(3.0),LINES_PLT_FG[2].append(15),LINES_PLT_FG[7].append(0.020),LINES_PLT_FG[8].append(00.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('Fe2-4')  #ISM
LINES_PLT_FG[0].append(2370.5),LINES_PLT_FG[1].append(5.0),LINES_PLT_FG[2].append(20),LINES_PLT_FG[7].append(0.020),LINES_PLT_FG[8].append(02.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('Fe2-*')  #ISM
#LINES_PLT_FG[0].append(2402.6),LINES_PLT_FG[1].append(3.0),LINES_PLT_FG[2].append(15),LINES_PLT_FG[7].append(0.020),LINES_PLT_FG[8].append(00.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('Fe2-6')  #ISM
LINES_PLT_FG[0].append(2608.0),LINES_PLT_FG[1].append(5.0),LINES_PLT_FG[2].append(10),LINES_PLT_FG[7].append(0.020),LINES_PLT_FG[8].append(-2.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('Fe2-7')  #ISM previouus 2593.7 2608CORRECTED OTHELO
#LINES_PLT_FG[0].append(2796.0),LINES_PLT_FG[1].append(6.0),LINES_PLT_FG[2].append(15),LINES_PLT_FG[7].append(0.200),LINES_PLT_FG[8].append(+4.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('Mg2-1')  #ISM
#LINES_PLT_FG[0].append(2803.0),LINES_PLT_FG[1].append(4.0),LINES_PLT_FG[2].append(20),LINES_PLT_FG[7].append(0.200),LINES_PLT_FG[8].append(00.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('Mg2-2')  #ISM
LINES_PLT_FG[0].append(2852.0),LINES_PLT_FG[1].append(3.0),LINES_PLT_FG[2].append(15),LINES_PLT_FG[7].append(0.020),LINES_PLT_FG[8].append(+2.),LINES_PLT_FG[10].append(0.001),LINES_PLT_FG[3].append('Mg1-2')  #ISM
###############CENTER#####################WIDTH-FIT##################WIDTH-PLT##############CTR-FIT-BNDS###################CTR-OFFSET#############LNE_AMP_BNDS############LNE_NME_HDR_CMT################

###########MARKER##################LNE-PLT-TXT######################LNE-HDR-NME###################LNE-LBL-NME#########
##LINES_PLT_FG[6].append('.'),LINES_PLT_FG[4].append('Ly-l')    ,LINES_PLT_FG[5].append('L01'),LINES_PLT_FG[11].append(7),LINES_PLT_FG[9].append('Ly$_{limit}$')      #Continuum Break
##LINES_PLT_FG[6].append('X'),LINES_PLT_FG[4].append('Ly-g')    ,LINES_PLT_FG[5].append('L02') ,LINES_PLT_FG[11].append(7),LINES_PLT_FG[9].append('Ly$_{\\gamma}$')     #HI absoprtion
##LINES_PLT_FG[6].append('o'),LINES_PLT_FG[4].append('Ly-B')    ,LINES_PLT_FG[5].append('L03') ,LINES_PLT_FG[11].append(7),LINES_PLT_FG[9].append('Ly$_{\\beta}$')     #HI absoprtion
LINES_PLT_FG[6].append('v'),LINES_PLT_FG[4].append('SiII-1')  ,LINES_PLT_FG[5].append('S04') ,LINES_PLT_FG[11].append(1),LINES_PLT_FG[9].append('SiII-1')   #ISM, blend 1190+1193 
LINES_PLT_FG[6].append('^'),LINES_PLT_FG[4].append('Ly-a')    ,LINES_PLT_FG[5].append('L05'),LINES_PLT_FG[11].append(7),LINES_PLT_FG[9].append('Ly$_{\\alpha}$')      #HI emission/absorption
LINES_PLT_FG[6].append('<'),LINES_PLT_FG[4].append('NV')      ,LINES_PLT_FG[5].append('N06'),LINES_PLT_FG[11].append(6),LINES_PLT_FG[9].append('NV')        #
LINES_PLT_FG[6].append('<'),LINES_PLT_FG[4].append('  ')      ,LINES_PLT_FG[5].append('N07'),LINES_PLT_FG[11].append(6),LINES_PLT_FG[9].append('NV')        #
##LINES_PLT_FG[6].append('<'),LINES_PLT_FG[4].append('NV-D')    ,LINES_PLT_FG[5].append('N08'),LINES_PLT_FG[11].append(0),LINES_PLT_FG[9].append('NV')        #
LINES_PLT_FG[6].append('>'),LINES_PLT_FG[4].append('SiII')    ,LINES_PLT_FG[5].append('S09'),LINES_PLT_FG[11].append(1),LINES_PLT_FG[9].append('SiII')     #ISM
LINES_PLT_FG[6].append('>'),LINES_PLT_FG[4].append('SiII*')   ,LINES_PLT_FG[5].append('S10'),LINES_PLT_FG[11].append(3),LINES_PLT_FG[9].append('SiII*')     #SiII** Nebular Shapley
LINES_PLT_FG[6].append('8'),LINES_PLT_FG[4].append('OI+SiII') ,LINES_PLT_FG[5].append('O12'),LINES_PLT_FG[11].append(1),LINES_PLT_FG[9].append('OI+SiII')   #ISM, blend 1303
LINES_PLT_FG[6].append('>'),LINES_PLT_FG[4].append('SiII*')   ,LINES_PLT_FG[5].append('S11'),LINES_PLT_FG[11].append(3),LINES_PLT_FG[9].append('SiII*')     #SiII** Nebular Shapley
LINES_PLT_FG[6].append('s'),LINES_PLT_FG[4].append('CII')     ,LINES_PLT_FG[5].append('C13'),LINES_PLT_FG[11].append(1),LINES_PLT_FG[9].append('CII')       #ISM
LINES_PLT_FG[6].append('s'),LINES_PLT_FG[4].append('OIV')     ,LINES_PLT_FG[5].append('O14'),LINES_PLT_FG[11].append(6),LINES_PLT_FG[9].append('OIV')       #Wind Line http://www.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
#LINES_PLT_FG[6].append('^'),LINES_PLT_FG[4].append('OV')      ,LINES_PLT_FG[5].append('O90'),LINES_PLT_FG[11].append(3),LINES_PLT_FG[9].append('OV')       #New-Update Wind Line #1371 OV OTELO DATABASE http://research.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
LINES_PLT_FG[6].append('*'),LINES_PLT_FG[4].append('SiIV')    ,LINES_PLT_FG[5].append('S15'),LINES_PLT_FG[11].append(2),LINES_PLT_FG[9].append('SiIV')      #ISM
LINES_PLT_FG[6].append('p'),LINES_PLT_FG[4].append('SiIV')  ,LINES_PLT_FG[5].append('S16'),LINES_PLT_FG[11].append(2),LINES_PLT_FG[9].append('SiIV')      #ISM
#LINES_PLT_FG[6].append('8'),LINES_PLT_FG[4].append('SiIV-D')  ,LINES_PLT_FG[5].append('S17'),LINES_PLT_FG[11].append(2),LINES_PLT_FG[9].append('SiIV')      #ISM
LINES_PLT_FG[6].append('s'),LINES_PLT_FG[4].append('SiIII')   ,LINES_PLT_FG[5].append('S18'),LINES_PLT_FG[11].append(6),LINES_PLT_FG[9].append('SiIII')      #Halliday+08GMASS photospheric abs lin
#LINES_PLT_FG[6].append('p'),LINES_PLT_FG[4].append('FeV')     ,LINES_PLT_FG[5].append('F19'),LINES_PLT_FG[11].append(0),LINES_PLT_FG[9].append('FeV')       #http://research.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
#LINES_PLT_FG[6].append('P'),LINES_PLT_FG[4].append('SiII')    ,LINES_PLT_FG[5].append('S20'),LINES_PLT_FG[11].append(3),LINES_PLT_FG[9].append('SiII')      #Halliday+08GMASS photospheric abs lin
#LINES_PLT_FG[6].append('*'),LINES_PLT_FG[4].append('NIV')     ,LINES_PLT_FG[5].append('N21'),LINES_PLT_FG[11].append(0),LINES_PLT_FG[9].append('NIV')       #http://research.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
LINES_PLT_FG[6].append('h'),LINES_PLT_FG[4].append('SV')      ,LINES_PLT_FG[5].append('S22'),LINES_PLT_FG[11].append(6),LINES_PLT_FG[9].append('SV')        #Halliday+08GMASS photospheric abs lin
LINES_PLT_FG[6].append('H'),LINES_PLT_FG[4].append('SiII')    ,LINES_PLT_FG[5].append('S23'),LINES_PLT_FG[11].append(1),LINES_PLT_FG[9].append('SiII')      #ISM photospheric abs lin
LINES_PLT_FG[6].append('X'),LINES_PLT_FG[4].append('CIV')     ,LINES_PLT_FG[5].append('C24'),LINES_PLT_FG[11].append(2),LINES_PLT_FG[9].append('CIV')       #ISM Blend 1548.2 + 1550.8 1
LINES_PLT_FG[6].append('P'),LINES_PLT_FG[4].append('     ')   ,LINES_PLT_FG[5].append('C25'),LINES_PLT_FG[11].append(2),LINES_PLT_FG[9].append('CIV')       #ISM Blend 1548.2 + 1550.8 2
##LINES_PLT_FG[6].append('X'),LINES_PLT_FG[4].append('CIV-D')   ,LINES_PLT_FG[5].append('C26'),LINES_PLT_FG[11].append(2),LINES_PLT_FG[9].append('CIV')       #ISM Blend 1548.2 + 1550.8
LINES_PLT_FG[6].append('D'),LINES_PLT_FG[4].append('FeII')    ,LINES_PLT_FG[5].append('F27'),LINES_PLT_FG[11].append(1),LINES_PLT_FG[9].append('FeII')      #ISM
LINES_PLT_FG[6].append('d'),LINES_PLT_FG[4].append('HeII')    ,LINES_PLT_FG[5].append('H28'),LINES_PLT_FG[11].append(5),LINES_PLT_FG[9].append('HeII')      #Nebular
#LINES_PLT_FG[6].append('H'),LINES_PLT_FG[4].append('CI')      ,LINES_PLT_FG[5].append('C91'),LINES_PLT_FG[11].append(0),LINES_PLT_FG[9].append('CI')              #New-Update Fanelli #1656 CI OTELO DATABASE http://research.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
LINES_PLT_FG[6].append('X'),LINES_PLT_FG[4].append('OIII]')   ,LINES_PLT_FG[5].append('O92'),LINES_PLT_FG[11].append(4),LINES_PLT_FG[9].append('OIII]')     #New-Update Fanelli #1660 OII-Doublet OTELO DATABASE http://research.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
LINES_PLT_FG[6].append('P'),LINES_PLT_FG[4].append('OIII]')   ,LINES_PLT_FG[5].append('O93'),LINES_PLT_FG[11].append(4),LINES_PLT_FG[9].append('OIII]')     #New-Update Fanelli #1666 OII-Doublet OTELO DATABASE http://research.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
LINES_PLT_FG[6].append('o'),LINES_PLT_FG[4].append('AlII')    ,LINES_PLT_FG[5].append('A29'),LINES_PLT_FG[11].append(1),LINES_PLT_FG[9].append('AlII')      #ISM
LINES_PLT_FG[6].append('X'),LINES_PLT_FG[4].append('NiII')    ,LINES_PLT_FG[5].append('N30'),LINES_PLT_FG[11].append(1),LINES_PLT_FG[9].append('NiII')      #Halliday+08GMASS 
LINES_PLT_FG[6].append('o'),LINES_PLT_FG[4].append('NIV')     ,LINES_PLT_FG[5].append('N31'),LINES_PLT_FG[11].append(5),LINES_PLT_FG[9].append('NIV')       #Halliday+08GMASS photospheric abs lin
#LINES_PLT_FG[6].append('o'),LINES_PLT_FG[4].append('SiIV')    ,LINES_PLT_FG[5].append('S32'),LINES_PLT_FG[11].append(0),LINES_PLT_FG[9].append('SIV')       #Fanelli BL 1719 http://www.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
LINES_PLT_FG[6].append('v'),LINES_PLT_FG[4].append('NiII')    ,LINES_PLT_FG[5].append('N33'),LINES_PLT_FG[11].append(1),LINES_PLT_FG[9].append('NiII')      #Halliday+08GMASS 
LINES_PLT_FG[6].append('^'),LINES_PLT_FG[4].append('NiII')    ,LINES_PLT_FG[5].append('N34'),LINES_PLT_FG[11].append(1),LINES_PLT_FG[9].append('NiII')      #Halliday+08GMASS 
##LINES_PLT_FG[6].append('<'),LINES_PLT_FG[4].append('NiII')    ,LINES_PLT_FG[5].append('N35'),LINES_PLT_FG[11].append(0),LINES_PLT_FG[9].append('NiII')      #Halliday+08GMASS Doublet
LINES_PLT_FG[6].append('>'),LINES_PLT_FG[4].append('SiII')    ,LINES_PLT_FG[5].append('S36'),LINES_PLT_FG[11].append(1),LINES_PLT_FG[9].append('SiII')      #Halliday+08GMASS 
LINES_PLT_FG[6].append('8'),LINES_PLT_FG[4].append('AlIII')   ,LINES_PLT_FG[5].append('A37'),LINES_PLT_FG[11].append(2),LINES_PLT_FG[9].append('AlIII')      #ISM H+? Old FeII
LINES_PLT_FG[6].append('s'),LINES_PLT_FG[4].append('AlIII')   ,LINES_PLT_FG[5].append('A38'),LINES_PLT_FG[11].append(2),LINES_PLT_FG[9].append('AlIII')      #ISM H+? Old FeII
LINES_PLT_FG[6].append('s'),LINES_PLT_FG[4].append('SiIII]')  ,LINES_PLT_FG[5].append('F94'),LINES_PLT_FG[11].append(4),LINES_PLT_FG[9].append('SiIII]')    #New-Update ISM H+? #1889 SiIII OTELO DATABASE http://research.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
LINES_PLT_FG[6].append('*'),LINES_PLT_FG[4].append('CIII]')   ,LINES_PLT_FG[5].append('C39'),LINES_PLT_FG[11].append(4),LINES_PLT_FG[9].append('CIII]')      #Nebular Blend 1907 + 1909
LINES_PLT_FG[6].append('*'),LINES_PLT_FG[4].append('CIII]')   ,LINES_PLT_FG[5].append('C40'),LINES_PLT_FG[11].append(4),LINES_PLT_FG[9].append('CIII]')      #Nebular Blend 1907 + 1909
LINES_PLT_FG[6].append('p'),LINES_PLT_FG[4].append('MgI')     ,LINES_PLT_FG[5].append('M41'),LINES_PLT_FG[11].append(1),LINES_PLT_FG[9].append('MgI')       #Halliday+08GMASS 
LINES_PLT_FG[6].append('8'),LINES_PLT_FG[4].append('ZnII')    ,LINES_PLT_FG[5].append('Z42'),LINES_PLT_FG[11].append(1),LINES_PLT_FG[9].append('ZnII')      #Halliday+08GMASS 
LINES_PLT_FG[6].append('s'),LINES_PLT_FG[4].append('CII]')    ,LINES_PLT_FG[5].append('C43'),LINES_PLT_FG[11].append(4),LINES_PLT_FG[9].append('CII]')       #http://www.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
LINES_PLT_FG[6].append('p'),LINES_PLT_FG[4].append('FeII')    ,LINES_PLT_FG[5].append('F44'),LINES_PLT_FG[11].append(1),LINES_PLT_FG[9].append('FeII')      #ISM
LINES_PLT_FG[6].append('p'),LINES_PLT_FG[4].append('FeII*')   ,LINES_PLT_FG[5].append('F45'),LINES_PLT_FG[11].append(3),LINES_PLT_FG[9].append('FeII*')     #ISM
#LINES_PLT_FG[6].append('P'),LINES_PLT_FG[4].append('FeII')    ,LINES_PLT_FG[5].append('F46'),LINES_PLT_FG[11].append(1),LINES_PLT_FG[9].append('FeII')      #ISM
LINES_PLT_FG[6].append('*'),LINES_PLT_FG[4].append('FeII')    ,LINES_PLT_FG[5].append('F47'),LINES_PLT_FG[11].append(1),LINES_PLT_FG[9].append('FeII')      #ISM
#LINES_PLT_FG[6].append('h'),LINES_PLT_FG[4].append('MgII')   ,LINES_PLT_FG[5].append('M48'),LINES_PLT_FG[11].append(1),LINES_PLT_FG[9].append('MgII')      #ISM
#LINES_PLT_FG[6].append('h'),LINES_PLT_FG[4].append('     ')   ,LINES_PLT_FG[5].append('M49'),LINES_PLT_FG[11].append(1),LINES_PLT_FG[9].append('MgII')      #ISM
LINES_PLT_FG[6].append('h'),LINES_PLT_FG[4].append('MgI')    ,LINES_PLT_FG[5].append('M50'),LINES_PLT_FG[11].append(0),LINES_PLT_FG[9].append('MgI')       #ISM
###########MARKER##################LNE-PLT-TXT######################LNE-HDR-NME###################LNE-LBL-NME#########
#############################################################################################################################
#############################################################################################################################

#############################################################################################################################
#############################################################################################################################

lw        = 5#pre_mask_lw
LINES_PLT_BG = ([],[],[],[],[],[],[],[],[],[],[],[])
###############CENTER#####################WIDTH-FIT##################WIDTH-PLT##############CTR-FIT-BNDS###################CTR-OFFSET#############LNE_AMP_BNDS############LNE_NME_HDR_CMT################
#LINES_PLT_BG[0].append( 912.0),LINES_PLT_BG[1].append(lw),LINES_PLT_BG[2].append(15),LINES_PLT_BG[7].append(0.002),LINES_PLT_BG[8].append(00),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('LyL')      #Continuum Break
#LINES_PLT_BG[0].append( 972.0),LINES_PLT_BG[1].append(3.0),LINES_PLT_BG[2].append(15),LINES_PLT_BG[7].append(0.002),LINES_PLT_BG[8].append(00),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('Lyg')     #HI absoprtion
#LINES_PLT_BG[0].append(1025.2),LINES_PLT_BG[1].append(lw),LINES_PLT_BG[2].append(15),LINES_PLT_BG[7].append(0.002),LINES_PLT_BG[8].append(00.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('LyB')     #HI absoprtion
#LINES_PLT_BG[0].append(1192.0),LINES_PLT_BG[1].append(lw),LINES_PLT_BG[2].append(15),LINES_PLT_BG[7].append(0.002),LINES_PLT_BG[8].append(00.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('Si2-1')   #ISM, blend 1190+1193 
LINES_PLT_BG[0].append(1215.7),LINES_PLT_BG[1].append(5.0),LINES_PLT_BG[2].append(20),LINES_PLT_BG[7].append(0.020),LINES_PLT_BG[8].append(00.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('Lya')    #HI emission/absorption
#LINES_PLT_BG[0].append(1238.0),LINES_PLT_BG[1].append(6.0),LINES_PLT_BG[2].append(18),LINES_PLT_BG[7].append(0.020),LINES_PLT_BG[8].append(00.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('N5-1')   #
#LINES_PLT_BG[0].append(1242.0),LINES_PLT_BG[1].append(6.0),LINES_PLT_BG[2].append(14),LINES_PLT_BG[7].append(0.020),LINES_PLT_BG[8].append(-2.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('N5-2')   #
##LINES_PLT_BG[0].append(1240.0),LINES_PLT_BG[1].append(6.0),LINES_PLT_BG[2].append(14),LINES_PLT_BG[7].append(0.020),LINES_PLT_BG[8].append(00.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('N5Dbl')  #
#LINES_PLT_BG[0].append(1260.4),LINES_PLT_BG[1].append(4.0),LINES_PLT_BG[2].append(18),LINES_PLT_BG[7].append(0.020),LINES_PLT_BG[8].append(00.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('Si2-1')   #ISM Si II 1263 triplet 
#LINES_PLT_BG[0].append(1264.0),LINES_PLT_BG[1].append(4.0),LINES_PLT_BG[2].append(18),LINES_PLT_BG[7].append(0.020),LINES_PLT_BG[8].append(00.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('Si2-*')   #SiII** Nebular Shapley
LINES_PLT_BG[0].append(1303.2),LINES_PLT_BG[1].append(3.0),LINES_PLT_BG[2].append(19),LINES_PLT_BG[7].append(0.020),LINES_PLT_BG[8].append(00.),LINES_PLT_BG[10].append(0.601),LINES_PLT_BG[3].append('OSD')    #ISM, blend 1303 OI-SiIII
#LINES_PLT_BG[0].append(1309.0),LINES_PLT_BG[1].append(4.0),LINES_PLT_BG[2].append(18),LINES_PLT_BG[7].append(0.020),LINES_PLT_BG[8].append(00.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('Si2-*')   #SiII** Nebular Shapley
LINES_PLT_BG[0].append(1334.5),LINES_PLT_BG[1].append(3.0),LINES_PLT_BG[2].append(19),LINES_PLT_BG[7].append(0.002),LINES_PLT_BG[8].append(04.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('C2-1')   #ISM Triplet
#LINES_PLT_BG[0].append(1343.0),LINES_PLT_BG[1].append(4.0),LINES_PLT_BG[2].append(22),LINES_PLT_BG[7].append(0.022),LINES_PLT_BG[8].append(00.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('O4')    #Wind Line http://www.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php ****O stars
#LINES_PLT_BG[0].append(1371.0),LINES_PLT_BG[1].append(4.0),LINES_PLT_BG[2].append(19),LINES_PLT_BG[7].append(0.022),LINES_PLT_BG[8].append(04.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('O5')    #New-Update Wind Line #1371 OV OTELO DATABASE http://research.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
LINES_PLT_BG[0].append(1393.8),LINES_PLT_BG[1].append(2.0),LINES_PLT_BG[2].append(18),LINES_PLT_BG[7].append(0.061),LINES_PLT_BG[8].append(00.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('Si4-1')  #ISM
#LINES_PLT_BG[0].append(1402.8),LINES_PLT_BG[1].append(2.0),LINES_PLT_BG[2].append(16),LINES_PLT_BG[7].append(0.062),LINES_PLT_BG[8].append(00.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('Si4-2')  #ISM
#LINES_PLT_BG[0].append(1404.0),LINES_PLT_BG[1].append(5.0),LINES_PLT_BG[2].append(18),LINES_PLT_BG[7].append(0.002),LINES_PLT_BG[8].append(04.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('Si4Dbl') #ISM
#LINES_PLT_BG[0].append(1417.2),LINES_PLT_BG[1].append(3.0),LINES_PLT_BG[2].append(16),LINES_PLT_BG[7].append(0.020),LINES_PLT_BG[8].append(00.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('Si3-3')  #Halliday+08GMASS photospheric abs lin ERROR IN ID AS Si2 sholuld be Si3   ****O stars
#LINES_PLT_BG[0].append(1453.0),LINES_PLT_BG[1].append(3.0),LINES_PLT_BG[2].append(16),LINES_PLT_BG[7].append(0.020),LINES_PLT_BG[8].append(05.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('Fe5')    #http://research.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
#LINES_PLT_BG[0].append(1485.4),LINES_PLT_BG[1].append(4.0),LINES_PLT_BG[2].append(18),LINES_PLT_BG[7].append(0.020),LINES_PLT_BG[8].append(00.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('Si2-4')  #Halliday+08GMASS photospheric abs lin
#LINES_PLT_BG[0].append(1486.0),LINES_PLT_BG[1].append(4.0),LINES_PLT_BG[2].append(16),LINES_PLT_BG[7].append(0.020),LINES_PLT_BG[8].append(00.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('N4-1')   #http://research.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
#LINES_PLT_BG[0].append(1501.8),LINES_PLT_BG[1].append(4.0),LINES_PLT_BG[2].append(16),LINES_PLT_BG[7].append(0.002),LINES_PLT_BG[8].append(00.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('S5')     #Halliday+08GMASS photospheric abs lin ****O stars
LINES_PLT_BG[0].append(1526.7),LINES_PLT_BG[1].append(5.0),LINES_PLT_BG[2].append(19),LINES_PLT_BG[7].append(0.002),LINES_PLT_BG[8].append(00.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('Si2-5')  #ISM photospheric abs lin
LINES_PLT_BG[0].append(1548.2),LINES_PLT_BG[1].append(4.0),LINES_PLT_BG[2].append(18),LINES_PLT_BG[7].append(0.002),LINES_PLT_BG[8].append(00.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('C4-1')   #ISM Blend 1548.2 + 1550.8 1
#LINES_PLT_BG[0].append(1550.7),LINES_PLT_BG[1].append(4.0),LINES_PLT_BG[2].append(18),LINES_PLT_BG[7].append(0.002),LINES_PLT_BG[8].append(00.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('C4-2')   #ISM Blend 1548.2 + 1550.8 2
##LINES_PLT_BG[0].append(1549.1),LINES_PLT_BG[1].append(4.0),LINES_PLT_BG[2].append(18),LINES_PLT_BG[7].append(0.002),LINES_PLT_BG[8].append(00.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('CDbl')   #ISM Blend 1548.2 + 1550.8
#LINES_PLT_BG[0].append(1608.5),LINES_PLT_BG[1].append(5.5),LINES_PLT_BG[2].append(20),LINES_PLT_BG[7].append(0.002),LINES_PLT_BG[8].append(00.),LINES_PLT_BG[10].append(0.901),LINES_PLT_BG[3].append('Fe2-1')  #ISM
#LINES_PLT_BG[0].append(1640.0),LINES_PLT_BG[1].append(5.0),LINES_PLT_BG[2].append(20),LINES_PLT_BG[7].append(0.020),LINES_PLT_BG[8].append(00.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('He2')    #Nebular
#LINES_PLT_BG[0].append(1656.0),LINES_PLT_BG[1].append(2.0),LINES_PLT_BG[2].append(18),LINES_PLT_BG[7].append(0.020),LINES_PLT_BG[8].append(00.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('C1')     #New-Update Fanelli #1656 CI OTELO DATABASE http://research.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
#LINES_PLT_BG[0].append(1660.0),LINES_PLT_BG[1].append(2.0),LINES_PLT_BG[2].append(18),LINES_PLT_BG[7].append(0.020),LINES_PLT_BG[8].append(00.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('O3]-1')   #New-Update Fanelli #1660 OII-Doublet OTELO DATABASE http://research.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
#LINES_PLT_BG[0].append(1666.0),LINES_PLT_BG[1].append(4.0),LINES_PLT_BG[2].append(18),LINES_PLT_BG[7].append(0.020),LINES_PLT_BG[8].append(02.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('O3]-2')   #New-Update Fanelli #1666 OII-Doublet OTELO DATABASE http://research.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
LINES_PLT_BG[0].append(1670.8),LINES_PLT_BG[1].append(4.0),LINES_PLT_BG[2].append(18),LINES_PLT_BG[7].append(0.020),LINES_PLT_BG[8].append(00.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('Al2')    #ISM
#LINES_PLT_BG[0].append(1709.6),LINES_PLT_BG[1].append(2.5),LINES_PLT_BG[2].append(18),LINES_PLT_BG[7].append(0.002),LINES_PLT_BG[8].append(00.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('Ni2-1')  #Halliday+08GMASS 
#LINES_PLT_BG[0].append(1718.5),LINES_PLT_BG[1].append(2.5),LINES_PLT_BG[2].append(18),LINES_PLT_BG[7].append(0.020),LINES_PLT_BG[8].append(00.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('N4-2')   #Halliday+08GMASS photospheric abs lin
#LINES_PLT_BG[0].append(1722.0),LINES_PLT_BG[1].append(3.0),LINES_PLT_BG[2].append(18),LINES_PLT_BG[7].append(0.020),LINES_PLT_BG[8].append(+5.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('Si4-3')  #Fanelli BL 1719 http://www.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
#LINES_PLT_BG[0].append(1741.5),LINES_PLT_BG[1].append(4.0),LINES_PLT_BG[2].append(18),LINES_PLT_BG[7].append(0.020),LINES_PLT_BG[8].append(00.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('Ni2-2')  #Halliday+08GMASS 
#LINES_PLT_BG[0].append(1751.9),LINES_PLT_BG[1].append(4.0),LINES_PLT_BG[2].append(18),LINES_PLT_BG[7].append(0.020),LINES_PLT_BG[8].append(00.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('Ni2-3')  #Halliday+08GMASS 
##LINES_PLT_BG[0].append(1746.9),LINES_PLT_BG[1].append(3.0),LINES_PLT_BG[2].append(18),LINES_PLT_BG[7].append(0.002),LINES_PLT_BG[8].append(00.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('Ni2Dbl') #Halliday+08GMASS Doublet
#LINES_PLT_BG[0].append(1808.0),LINES_PLT_BG[1].append(3.0),LINES_PLT_BG[2].append(18),LINES_PLT_BG[7].append(0.020),LINES_PLT_BG[8].append(+2.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('Si2-6')  #Halliday+08GMASS 
#LINES_PLT_BG[0].append(1854.7),LINES_PLT_BG[1].append(4.0),LINES_PLT_BG[2].append(20),LINES_PLT_BG[7].append(0.020),LINES_PLT_BG[8].append(-2.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('Al3-1')  #ISM H+? Fe2-2 Old
LINES_PLT_BG[0].append(1862.8),LINES_PLT_BG[1].append(4.0),LINES_PLT_BG[2].append(17),LINES_PLT_BG[7].append(0.020),LINES_PLT_BG[8].append(+2.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('Al3-2')  #ISM H+? Fe2-3 Old
#LINES_PLT_BG[0].append(1889.0),LINES_PLT_BG[1].append(4.0),LINES_PLT_BG[2].append(17),LINES_PLT_BG[7].append(0.020),LINES_PLT_BG[8].append(+2.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('Si3]')  #New-Update ISM H+? #1889 SiIII OTELO DATABASE http://research.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
#LINES_PLT_BG[0].append(1907.0),LINES_PLT_BG[1].append(4.0),LINES_PLT_BG[2].append(15),LINES_PLT_BG[7].append(0.020),LINES_PLT_BG[8].append(00.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('C3]-1')   #Nebular Blend 1907 + 1909
#LINES_PLT_BG[0].append(1909.0),LINES_PLT_BG[1].append(4.0),LINES_PLT_BG[2].append(15),LINES_PLT_BG[7].append(0.020),LINES_PLT_BG[8].append(+15),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('C3]-2')   #Nebular Blend 1907 + 1909
LINES_PLT_BG[0].append(2028.4),LINES_PLT_BG[1].append(4.0),LINES_PLT_BG[2].append(15),LINES_PLT_BG[7].append(0.022),LINES_PLT_BG[8].append(00.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('Mg1-1')  #Halliday+08GMASS 
LINES_PLT_BG[0].append(2062.6),LINES_PLT_BG[1].append(3.0),LINES_PLT_BG[2].append(15),LINES_PLT_BG[7].append(0.020),LINES_PLT_BG[8].append(+2.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('Zn2')    #Halliday+08GMASS 
LINES_PLT_BG[0].append(2325.5),LINES_PLT_BG[1].append(8.0),LINES_PLT_BG[2].append(15),LINES_PLT_BG[7].append(0.020),LINES_PLT_BG[8].append(00.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('C2]')   #http://www.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
LINES_PLT_BG[0].append(2343.5),LINES_PLT_BG[1].append(3.0),LINES_PLT_BG[2].append(15),LINES_PLT_BG[7].append(0.020),LINES_PLT_BG[8].append(00.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('Fe2-4')  #ISM
LINES_PLT_BG[0].append(2370.5),LINES_PLT_BG[1].append(5.0),LINES_PLT_BG[2].append(20),LINES_PLT_BG[7].append(0.020),LINES_PLT_BG[8].append(02.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('Fe2-*')  #ISM
#LINES_PLT_BG[0].append(2402.6),LINES_PLT_BG[1].append(3.0),LINES_PLT_BG[2].append(15),LINES_PLT_BG[7].append(0.020),LINES_PLT_BG[8].append(00.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('Fe2-6')  #ISM
LINES_PLT_BG[0].append(2608.0),LINES_PLT_BG[1].append(5.0),LINES_PLT_BG[2].append(10),LINES_PLT_BG[7].append(0.020),LINES_PLT_BG[8].append(-2.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('Fe2-7')  #ISM previouus 2593.7 2608CORRECTED OTHELO
#LINES_PLT_BG[0].append(2796.0),LINES_PLT_BG[1].append(6.0),LINES_PLT_BG[2].append(15),LINES_PLT_BG[7].append(0.200),LINES_PLT_BG[8].append(+4.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('Mg2-1')  #ISM
#LINES_PLT_BG[0].append(2803.0),LINES_PLT_BG[1].append(4.0),LINES_PLT_BG[2].append(20),LINES_PLT_BG[7].append(0.200),LINES_PLT_BG[8].append(00.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('Mg2-2')  #ISM
LINES_PLT_BG[0].append(2852.0),LINES_PLT_BG[1].append(3.0),LINES_PLT_BG[2].append(15),LINES_PLT_BG[7].append(0.020),LINES_PLT_BG[8].append(+2.),LINES_PLT_BG[10].append(0.001),LINES_PLT_BG[3].append('Mg1-2')  #ISM
###############CENTER#####################WIDTH-FIT##################WIDTH-PLT##############CTR-FIT-BNDS###################CTR-OFFSET#############LNE_AMP_BNDS############LNE_NME_HDR_CMT################

###########MARKER##################LNE-PLT-TXT######################LNE-HDR-NME###################LNE-LBL-NME#########
##LINES_PLT_BG[6].append('.'),LINES_PLT_BG[4].append('Ly-l')    ,LINES_PLT_BG[5].append('L01'),LINES_PLT_BG[11].append(7),LINES_PLT_BG[9].append('Ly$_{limit}$')      #Continuum Break
##LINES_PLT_BG[6].append('X'),LINES_PLT_BG[4].append('Ly-g')    ,LINES_PLT_BG[5].append('L02') ,LINES_PLT_BG[11].append(7),LINES_PLT_BG[9].append('Ly$_{\\gamma}$')     #HI absoprtion
##LINES_PLT_BG[6].append('o'),LINES_PLT_BG[4].append('Ly-B')    ,LINES_PLT_BG[5].append('L03') ,LINES_PLT_BG[11].append(7),LINES_PLT_BG[9].append('Ly$_{\\beta}$')     #HI absoprtion
#LINES_PLT_BG[6].append('v'),LINES_PLT_BG[4].append('SiII-1')  ,LINES_PLT_BG[5].append('S04') ,LINES_PLT_BG[11].append(1),LINES_PLT_BG[9].append('SiII-1')   #ISM, blend 1190+1193 
LINES_PLT_BG[6].append('^'),LINES_PLT_BG[4].append('Ly-a')    ,LINES_PLT_BG[5].append('L05'),LINES_PLT_BG[11].append(7),LINES_PLT_BG[9].append('Ly$_{\\alpha}$')      #HI emission/absorption
#LINES_PLT_BG[6].append('<'),LINES_PLT_BG[4].append('NV')      ,LINES_PLT_BG[5].append('N06'),LINES_PLT_BG[11].append(6),LINES_PLT_BG[9].append('NV')        #
#LINES_PLT_BG[6].append('<'),LINES_PLT_BG[4].append('  ')      ,LINES_PLT_BG[5].append('N07'),LINES_PLT_BG[11].append(6),LINES_PLT_BG[9].append('NV')        #
##LINES_PLT_BG[6].append('<'),LINES_PLT_BG[4].append('NV-D')    ,LINES_PLT_BG[5].append('N08'),LINES_PLT_BG[11].append(0),LINES_PLT_BG[9].append('NV')        #
#LINES_PLT_BG[6].append('>'),LINES_PLT_BG[4].append('SiII')    ,LINES_PLT_BG[5].append('S09'),LINES_PLT_BG[11].append(1),LINES_PLT_BG[9].append('SiII')     #ISM
#LINES_PLT_BG[6].append('>'),LINES_PLT_BG[4].append('SiII*')   ,LINES_PLT_BG[5].append('S10'),LINES_PLT_BG[11].append(3),LINES_PLT_BG[9].append('SiII*')     #SiII** Nebular Shapley
LINES_PLT_BG[6].append('8'),LINES_PLT_BG[4].append('OI+SiII') ,LINES_PLT_BG[5].append('O12'),LINES_PLT_BG[11].append(1),LINES_PLT_BG[9].append('OI+SiII')   #ISM, blend 1303
#LINES_PLT_BG[6].append('>'),LINES_PLT_BG[4].append('SiII*')   ,LINES_PLT_BG[5].append('S11'),LINES_PLT_BG[11].append(3),LINES_PLT_BG[9].append('SiII*')     #SiII** Nebular Shapley
LINES_PLT_BG[6].append('s'),LINES_PLT_BG[4].append('CII')     ,LINES_PLT_BG[5].append('C13'),LINES_PLT_BG[11].append(1),LINES_PLT_BG[9].append('CII')       #ISM
#LINES_PLT_BG[6].append('s'),LINES_PLT_BG[4].append('OIV')     ,LINES_PLT_BG[5].append('O14'),LINES_PLT_BG[11].append(6),LINES_PLT_BG[9].append('OIV')       #Wind Line http://www.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
#LINES_PLT_BG[6].append('^'),LINES_PLT_BG[4].append('OV')      ,LINES_PLT_BG[5].append('O90'),LINES_PLT_BG[11].append(3),LINES_PLT_BG[9].append('OV')       #New-Update Wind Line #1371 OV OTELO DATABASE http://research.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
LINES_PLT_BG[6].append('*'),LINES_PLT_BG[4].append('SiIV')    ,LINES_PLT_BG[5].append('S15'),LINES_PLT_BG[11].append(2),LINES_PLT_BG[9].append('SiIV')      #ISM
#LINES_PLT_BG[6].append('p'),LINES_PLT_BG[4].append('SiIV')    ,LINES_PLT_BG[5].append('S16'),LINES_PLT_BG[11].append(2),LINES_PLT_BG[9].append('SiIV')      #ISM
#LINES_PLT_BG[6].append('8'),LINES_PLT_BG[4].append('SiIV-D')  ,LINES_PLT_BG[5].append('S17'),LINES_PLT_BG[11].append(2),LINES_PLT_BG[9].append('SiIV')      #ISM
#LINES_PLT_BG[6].append('s'),LINES_PLT_BG[4].append('SiIII')   ,LINES_PLT_BG[5].append('S18'),LINES_PLT_BG[11].append(6),LINES_PLT_BG[9].append('SiIII')      #Halliday+08GMASS photospheric abs lin
#LINES_PLT_BG[6].append('p'),LINES_PLT_BG[4].append('FeV')     ,LINES_PLT_BG[5].append('F19'),LINES_PLT_BG[11].append(0),LINES_PLT_BG[9].append('FeV')       #http://research.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
#LINES_PLT_BG[6].append('P'),LINES_PLT_BG[4].append('SiII')    ,LINES_PLT_BG[5].append('S20'),LINES_PLT_BG[11].append(3),LINES_PLT_BG[9].append('SiII')      #Halliday+08GMASS photospheric abs lin
#LINES_PLT_BG[6].append('*'),LINES_PLT_BG[4].append('NIV')     ,LINES_PLT_BG[5].append('N21'),LINES_PLT_BG[11].append(0),LINES_PLT_BG[9].append('NIV')       #http://research.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
#LINES_PLT_BG[6].append('h'),LINES_PLT_BG[4].append('SV')      ,LINES_PLT_BG[5].append('S22'),LINES_PLT_BG[11].append(6),LINES_PLT_BG[9].append('SV')        #Halliday+08GMASS photospheric abs lin
LINES_PLT_BG[6].append('H'),LINES_PLT_BG[4].append('SiII')    ,LINES_PLT_BG[5].append('S23'),LINES_PLT_BG[11].append(1),LINES_PLT_BG[9].append('SiII')      #ISM photospheric abs lin
LINES_PLT_BG[6].append('X'),LINES_PLT_BG[4].append('CIV')     ,LINES_PLT_BG[5].append('C24'),LINES_PLT_BG[11].append(2),LINES_PLT_BG[9].append('CIV')       #ISM Blend 1548.2 + 1550.8 1
#LINES_PLT_BG[6].append('P'),LINES_PLT_BG[4].append('CIV')     ,LINES_PLT_BG[5].append('C25'),LINES_PLT_BG[11].append(2),LINES_PLT_BG[9].append('CIV')       #ISM Blend 1548.2 + 1550.8 2
##LINES_PLT_BG[6].append('X'),LINES_PLT_BG[4].append('CIV-D')   ,LINES_PLT_BG[5].append('C26'),LINES_PLT_BG[11].append(2),LINES_PLT_BG[9].append('CIV')       #ISM Blend 1548.2 + 1550.8
#LINES_PLT_BG[6].append('D'),LINES_PLT_BG[4].append('FeII')    ,LINES_PLT_BG[5].append('F27'),LINES_PLT_BG[11].append(1),LINES_PLT_BG[9].append('FeII')      #ISM
#LINES_PLT_BG[6].append('d'),LINES_PLT_BG[4].append('HeII')    ,LINES_PLT_BG[5].append('H28'),LINES_PLT_BG[11].append(5),LINES_PLT_BG[9].append('HeII')      #Nebular
#LINES_PLT_BG[6].append('H'),LINES_PLT_BG[4].append('CI')      ,LINES_PLT_BG[5].append('C91'),LINES_PLT_BG[11].append(0),LINES_PLT_BG[9].append('CI')              #New-Update Fanelli #1656 CI OTELO DATABASE http://research.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
#LINES_PLT_BG[6].append('X'),LINES_PLT_BG[4].append('OIII]')   ,LINES_PLT_BG[5].append('O92'),LINES_PLT_BG[11].append(4),LINES_PLT_BG[9].append('OIII]')     #New-Update Fanelli #1660 OII-Doublet OTELO DATABASE http://research.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
#LINES_PLT_BG[6].append('P'),LINES_PLT_BG[4].append('OIII]')   ,LINES_PLT_BG[5].append('O93'),LINES_PLT_BG[11].append(4),LINES_PLT_BG[9].append('OIII]')     #New-Update Fanelli #1666 OII-Doublet OTELO DATABASE http://research.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
LINES_PLT_BG[6].append('o'),LINES_PLT_BG[4].append('AlII')    ,LINES_PLT_BG[5].append('A29'),LINES_PLT_BG[11].append(1),LINES_PLT_BG[9].append('AlII')      #ISM
#LINES_PLT_BG[6].append('X'),LINES_PLT_BG[4].append('NiII')    ,LINES_PLT_BG[5].append('N30'),LINES_PLT_BG[11].append(1),LINES_PLT_BG[9].append('NiII')      #Halliday+08GMASS 
#LINES_PLT_BG[6].append('o'),LINES_PLT_BG[4].append('NIV')     ,LINES_PLT_BG[5].append('N31'),LINES_PLT_BG[11].append(5),LINES_PLT_BG[9].append('NIV')       #Halliday+08GMASS photospheric abs lin
#LINES_PLT_BG[6].append('o'),LINES_PLT_BG[4].append('SiIV')    ,LINES_PLT_BG[5].append('S32'),LINES_PLT_BG[11].append(0),LINES_PLT_BG[9].append('SIV')       #Fanelli BL 1719 http://www.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
#LINES_PLT_BG[6].append('v'),LINES_PLT_BG[4].append('NiII')    ,LINES_PLT_BG[5].append('N33'),LINES_PLT_BG[11].append(1),LINES_PLT_BG[9].append('NiII')      #Halliday+08GMASS 
#LINES_PLT_BG[6].append('^'),LINES_PLT_BG[4].append('NiII')    ,LINES_PLT_BG[5].append('N34'),LINES_PLT_BG[11].append(1),LINES_PLT_BG[9].append('NiII')      #Halliday+08GMASS 
##LINES_PLT_BG[6].append('<'),LINES_PLT_BG[4].append('NiII')    ,LINES_PLT_BG[5].append('N35'),LINES_PLT_BG[11].append(0),LINES_PLT_BG[9].append('NiII')      #Halliday+08GMASS Doublet
#LINES_PLT_BG[6].append('>'),LINES_PLT_BG[4].append('SiII')    ,LINES_PLT_BG[5].append('S36'),LINES_PLT_BG[11].append(1),LINES_PLT_BG[9].append('SiII')      #Halliday+08GMASS 
#LINES_PLT_BG[6].append('8'),LINES_PLT_BG[4].append('AlIII')   ,LINES_PLT_BG[5].append('A37'),LINES_PLT_BG[11].append(2),LINES_PLT_BG[9].append('AlIII')      #ISM H+? Old FeII
LINES_PLT_BG[6].append('s'),LINES_PLT_BG[4].append('AlIII')   ,LINES_PLT_BG[5].append('A38'),LINES_PLT_BG[11].append(2),LINES_PLT_BG[9].append('AlIII')      #ISM H+? Old FeII
#LINES_PLT_BG[6].append('s'),LINES_PLT_BG[4].append('SiIII]')  ,LINES_PLT_BG[5].append('F94'),LINES_PLT_BG[11].append(4),LINES_PLT_BG[9].append('SiIII]')    #New-Update ISM H+? #1889 SiIII OTELO DATABASE http://research.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
#LINES_PLT_BG[6].append('*'),LINES_PLT_BG[4].append('CIII]')   ,LINES_PLT_BG[5].append('C39'),LINES_PLT_BG[11].append(4),LINES_PLT_BG[9].append('CIII]')      #Nebular Blend 1907 + 1909
#LINES_PLT_BG[6].append('*'),LINES_PLT_BG[4].append('CIII]')   ,LINES_PLT_BG[5].append('C40'),LINES_PLT_BG[11].append(4),LINES_PLT_BG[9].append('CIII]')      #Nebular Blend 1907 + 1909
LINES_PLT_BG[6].append('p'),LINES_PLT_BG[4].append('MgI')     ,LINES_PLT_BG[5].append('M41'),LINES_PLT_BG[11].append(1),LINES_PLT_BG[9].append('MgI')       #Halliday+08GMASS 
LINES_PLT_BG[6].append('8'),LINES_PLT_BG[4].append('ZnII')    ,LINES_PLT_BG[5].append('Z42'),LINES_PLT_BG[11].append(1),LINES_PLT_BG[9].append('ZnII')      #Halliday+08GMASS 
LINES_PLT_BG[6].append('s'),LINES_PLT_BG[4].append('CII]')    ,LINES_PLT_BG[5].append('C43'),LINES_PLT_BG[11].append(4),LINES_PLT_BG[9].append('CII]')       #http://www.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php
LINES_PLT_BG[6].append('p'),LINES_PLT_BG[4].append('FeII')    ,LINES_PLT_BG[5].append('F44'),LINES_PLT_BG[11].append(1),LINES_PLT_BG[9].append('FeII')      #ISM
LINES_PLT_BG[6].append('p'),LINES_PLT_BG[4].append('FeII*')   ,LINES_PLT_BG[5].append('F45'),LINES_PLT_BG[11].append(3),LINES_PLT_BG[9].append('FeII*')     #ISM
#LINES_PLT_BG[6].append('P'),LINES_PLT_BG[4].append('FeII')    ,LINES_PLT_BG[5].append('F46'),LINES_PLT_BG[11].append(1),LINES_PLT_BG[9].append('FeII')      #ISM
LINES_PLT_BG[6].append('*'),LINES_PLT_BG[4].append('FeII')    ,LINES_PLT_BG[5].append('F47'),LINES_PLT_BG[11].append(1),LINES_PLT_BG[9].append('FeII')      #ISM
#LINES_PLT_BG[6].append('h'),LINES_PLT_BG[4].append('MgII')   ,LINES_PLT_BG[5].append('M48'),LINES_PLT_BG[11].append(1),LINES_PLT_BG[9].append('MgII')      #ISM
#LINES_PLT_BG[6].append('h'),LINES_PLT_BG[4].append('     ')   ,LINES_PLT_BG[5].append('M49'),LINES_PLT_BG[11].append(1),LINES_PLT_BG[9].append('MgII')      #ISM
LINES_PLT_BG[6].append('h'),LINES_PLT_BG[4].append('MgI')    ,LINES_PLT_BG[5].append('M50'),LINES_PLT_BG[11].append(0),LINES_PLT_BG[9].append('MgI')       #ISM
###########MARKER##################LNE-PLT-TXT######################LNE-HDR-NME###################LNE-LBL-NME#########
#############################################################################################################################
#############################################################################################################################

