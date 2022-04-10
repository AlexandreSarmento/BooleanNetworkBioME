function solutions = systemOfOdes(t,Y0,p)
  
  HIF1ac = Y0(1);
  HIF1aac = Y0(2);
  HIF1aan = Y0(3);
  HIF1an = Y0(4);
  HIF1anAC = Y0(5);
  FIHa = Y0(6);
  PHD_tot = Y0(7);
  PHDa = Y0(8);
  PNUTS = Y0(9);
  ATRa = Y0(10);
  p53c = Y0(11);
  p53n = Y0(12);
  p53nP = Y0(13);
  p53nAC = Y0(14);
  Mdm2c = Y0(15);
  Mdm2n = Y0(16);
  SHIP = Y0(17);
  PIP3 = Y0(18);
  AKTa = Y0(19);
  PFKL = Y0(20);
  VEGF = Y0(21);
  Bax = Y0(22);
  
  p.O2 = partialPressureO2(t);
  
  dHIF1acdt = p.ksHIF1ac - p.kinHIF1ac*HIF1ac + p.koutHIF1an*HIF1an - p.dHIF1ac*HIF1ac - ...
              p.kdeHIF1ac*FIHa*HIF1ac/(HIF1ac + p.jdeHIF1ac) - p.kd1HIF1ac*PHDa*HIF1ac/(HIF1ac + p.jd1HIF1ac) -...
              p.kd2HIF1ac*Mdm2c*HIF1ac/p.jd2HIF1ac/(1+HIF1ac/p.jd2HIF1ac+HIF1aac/p.jd2HIF1aac + p53c/p.jdp53c);

  dHIF1aacdt = p.kdeHIF1ac*FIHa*HIF1ac/(HIF1ac + p.jdeHIF1ac) - p.dHIF1aac*HIF1aac - ...
               p.kd1HIF1aac*PHDa*HIF1aac/(HIF1aac + p.jd1HIF1aac) -...
               p.kd2HIF1aac*Mdm2c*HIF1aac/p.jd2HIF1aac/(1+HIF1ac/p.jd2HIF1ac + HIF1aac/p.jd2HIF1aac + p53c/p.jdp53c) -...
               p.kinHIF1aac*HIF1aac + p.koutHIF1aan*HIF1aan;

  dHIF1aandt = p.kinHIF1aac*HIF1aac - p.koutHIF1aan*HIF1aan - p.dHIF1aan*HIF1aan - ...
               p.kd2HIF1aan*Mdm2n*HIF1aan/p.jd2HIF1aan/(1+HIF1aan/p.jd2HIF1aan + HIF1an/p.jd2HIF1an + p53n/p.jdp53n + 
               p53nP/p.jdp53nP);

  dHIF1andt = p.kinHIF1ac*HIF1ac - p.koutHIF1an*HIF1an - p.dHIF1an*HIF1an - ...
              p.kd2HIF1an*Mdm2n*HIF1an/p.jd2HIF1an/(1 + HIF1aan/p.jd2HIF1aan + ...
              HIF1an/p.jd2HIF1an + p53n/p.jdp53n + p53nP/p.jdp53nP)- ...
              p.kacHIF1an*p.p300*HIF1an/p.jacHIF1an/(1+HIF1an/p.jacHIF1an + ...
              p53nP/p.jacp53nP) + p.kdeHIF1anAC*HIF1anAC/(HIF1anAC + p.jdeHIF1anAC);

  dHIF1anACdt = p.kacHIF1an*p.p300*HIF1an/p.jacHIF1an/(1 + HIF1an/p.jacHIF1an + p53nP/p.jacp53nP) - ...
                p.kdeHIF1anAC*HIF1anAC/(HIF1anAC + p.jdeHIF1anAC) - p.dHIF1anAC*HIF1anAC;

  %HIF1a_tot=HIF1ac+HIF1aac+HIF1aan+HIF1an+HIF1anAC
  
  FIH = p.FIH_tot - FIHa;

  dFIHadt = p.kacFIH*p.O2*p.ss/(p.O2*p.ss + p.jO2FIH)*FIH/(FIH + p.jacFIH) - p.kdeFIH*FIHa/(FIHa + p.jdeFIH);

  dPHD_totdt = p.ksPHD_tot0 + p.ksPHD_tot1*HIF1aan^4/(HIF1aan^4 + p.jksPHD_tot1^4) + p.ksPHD_tot2*HIF1anAC^4/(HIF1anAC^4 + p.jksPHD_tot2^4) + ...
               p.ksPHD_tot3*HIF1an^4/(HIF1an^4 + p.jksPHD_tot3^4) - p.dPHD_tot*PHD_tot;

  PHD = PHD_tot - PHDa;
  
  dPHDadt = p.kacPHD*p.O2*p.ss/(p.O2*p.ss + p.jO2PHD)*PHD - p.kdePHD*PHDa;

  dPNUTSdt = p.ksPNUTS0 + p.ksPNUTS*HIF1anAC^4/(HIF1anAC^4 + p.jsPNUTS^4) - p.dPNUTS*PNUTS;
  
  ATR = p.ATR_tot - ATRa;
  
  dATRadt=(p.kacATR0+p.kacATR*p.KO2/(p.O2*p.ss+p.KO2)*ATRa)*ATR/(ATR+p.jacATR)-p.kdeATR*ATRa/(ATRa+p.jdeATR);

  dp53cdt = p.ksp53c - p.kinp53c*p53c + p.koutp53n*p53n - p.dp53c*p53c - ...
            p.kdp53c*Mdm2c*p53c/p.jdp53c/(1 + HIF1ac/p.jd2HIF1ac + HIF1aac/p.jd2HIF1aac + p53c/p.jdp53c);

  dp53ndt = p.kinp53c*p53c - p.koutp53n*p53n - p.dp53n*p53n - ...
           (p.kacp53n1*ATRa + p.kacp53n2*PNUTS)*p53n/(p53n + p.jacp53n) + ...
            p.kdep53nP*p53nP/(p53nP + p.jdep53nP)-....
            p.kdp53n*Mdm2n*p53n/p.jdp53n/(1 + HIF1aan/p.jd2HIF1aan + HIF1an/p.jd2HIF1an + p53n/p.jdp53n + p53nP/p.jdp53nP);

  dp53nPdt = (p.kacp53n1*ATRa + p.kacp53n2*PNUTS)*p53n/(p53n + p.jacp53n) - p.kdep53nP*p53nP/(p53nP + p.jdep53nP) - ...
              p.kdp53nP*Mdm2n*p53nP/p.jdp53nP/(1 + HIF1aan/p.jd2HIF1aan + HIF1an/p.jd2HIF1an + p53n/p.jdp53n + p53nP/p.jdp53nP) - ...
              p.kacp53nP*p.p300*p53nP/p.jacp53nP/(1 + HIF1an/p.jacHIF1an + p53nP/p.jacp53nP) + p.kdep53nAC*p53nAC/(p53nAC + p.jdep53nAC) - ...
              p.dp53nP*p53nP;

  dp53nACdt = p.kacp53nP*p.p300*p53nP/p.jacp53nP/(1 + HIF1an/p.jacHIF1an + p53nP/p.jacp53nP)- ...
              p.kdep53nAC*p53nAC/(p53nAC + p.jdep53nAC) - p.dp53nAC*p53nAC;

  %p53tot = p53c + p53n + p53nP + p53nAC;

  dMdm2cdt = p.ksMdm2c0 + p.ksMdm2c*p53nAC^4/(p53nAC^4 + p.jsMdm2c^4) - p.kdMdm2c*Mdm2c - ...
             p.kinMdm2c*AKTa*Mdm2c/(Mdm2c + p.jinMdm2c) + p.koutMdm2n*Mdm2n/(Mdm2n + p.joutMdm2n);

  dMdm2ndt = p.kinMdm2c*AKTa*Mdm2c/(Mdm2c + p.jinMdm2c) - p.koutMdm2n*Mdm2n/(Mdm2n + p.joutMdm2n) - p.dMdm2n*Mdm2n;

  %Mdm2_tot = Mdm2c + Mdm2n;

  dSHIPdt = p.ksSHIP0 + p.ksSHIP*p53nAC^4/(p53nAC^4 + p.jsSHIP^4) - p.dSHIP*SHIP;
  
  PIP2 = p.PIP_tot - PIP3;

  dPIP3dt = p.kp2*PIP2/(PIP2 + p.jp2) - p.kp3*SHIP*PIP3/(PIP3 + p.jp3);
  
  AKT = p.AKT_tot - AKTa;

  dAKTadt = p.kacAKT*PIP3*AKT/(AKT + p.jacAKT) - p.kdeAKT*AKTa/(AKTa + p.jdeAKT);

  dPFKLdt = p.ksPFKL0 + p.ksPFKL1*HIF1an^4/(HIF1an^4 + p.jsPFKL1^4) + p.ksPFKL2*HIF1aan^4/(HIF1aan^4 + p.jsPFKL2^4) - p.dPFKL*PFKL;

  dVEGFdt = p.ksVEGF0 + p.ksVEGF*HIF1anAC^4/(HIF1anAC^4 + p.jsVEGF^4) - p.dVEGF*VEGF;

  dBaxdt = p.ksBax0 + p.ksBax*p53nAC^4/(p53nAC^4 + p.jsBax^4) - p.dBax*Bax;
  
  solutions = [dHIF1acdt;dHIF1aacdt;dHIF1aandt;dHIF1andt;dHIF1anACdt;...
               dFIHadt;dPHD_totdt;dPHDadt;dPNUTSdt;dATRadt;...
               dp53cdt;dp53ndt;dp53nPdt;dp53nACdt;dMdm2cdt;...
               dMdm2ndt;dSHIPdt;dPIP3dt;dAKTadt;dPFKLdt;...
               dVEGFdt;dBaxdt];

end