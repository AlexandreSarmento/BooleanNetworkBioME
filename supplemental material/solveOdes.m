function solveOdes(p,initial_conditions,timeInterval,timePoints)
  
  [TIMES,OUTPUTS] = ode15s(@(time,output) systemOfOdes(time,output,p),timeInterval,initial_conditions);
  
  IndexOfTimePoints = timepointidx(timePoints,TIMES);
  HIF1a_tot = OUTPUTS(IndexOfTimePoints,1) + OUTPUTS(IndexOfTimePoints,2) + ...
              OUTPUTS(IndexOfTimePoints,3) + OUTPUTS(IndexOfTimePoints,4) + ...
              OUTPUTS(IndexOfTimePoints,5);
  TP53_tot = OUTPUTS(IndexOfTimePoints,11) + OUTPUTS(IndexOfTimePoints,12) + ...
             OUTPUTS(IndexOfTimePoints,13) + OUTPUTS(IndexOfTimePoints,14);

  MDM2_tot = OUTPUTS(IndexOfTimePoints,15) + OUTPUTS(IndexOfTimePoints,16);

  timeSerie2bin = horzcat(OUTPUTS(IndexOfTimePoints,:),HIF1a_tot,TP53_tot,MDM2_tot);
 
  inputParamsSet = strcat('norm_hyp2','_tf_',num2str(timePoints(end)/60),'hours_',...
                          '_ksp53c_',num2str(p.ksp53c),'_jsBax_',num2str(p.jsBax),'_jsPFKL1_',num2str(p.jsPFKL1),...
			                    '_jsVEGF_',num2str(p.jsVEGF),'_jsSHIP_',num2str(p.jsSHIP),'_ksPFKL1_',num2str(p.ksPFKL1),...
			                    '_ksVEGF_',num2str(p.ksVEGF),'_dVEGF_',num2str(p.dVEGF),'_dPFKL_',num2str(p.dPFKL));
  timeSerieFilePath = 'C:\Users\Alexandre Sarmento\Documents\octave\timeSerieWang\';
  %% GENERATE CSV FILE TO BINARIZATION
  cHeader = {'HIF1ac','HIF1aac','HIF1aan','HIF1an','HIF1anAC','FIH','PHD_tot','PHDa','PNUTS','ATR',...
             'p53c','p53n','p53nP','p53nAC','Mdm2c','Mdm2n','SHIP','PIP3','AKTa','PFKL','VEGF','Bax',...
             'HIF1a_{tot}','TP53_{tot}','MDM2_{tot}'};
  commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
  commaHeader = commaHeader(:)';
  commaHeader(end) = [];
  textHeader = cell2mat(commaHeader); %cHeader in text with commas
  %write header to file
  fid = fopen(strcat(timeSerieFilePath,inputParamsSet,'.csv'),'w'); 
  fprintf(fid,'%s\n',textHeader);
  fclose(fid);
  %write data to end of file
  dlmwrite(strcat(timeSerieFilePath,inputParamsSet,'.csv'),timeSerie2bin,'-append','precision',6);

  ## GENERATE PLOTS
  figsFilePath = strcat('C:\Users\Alexandre Sarmento\Documents\octave\figuresWang\',inputParamsSet);
  if not(isfolder(figsFilePath))
    mkdir(figsFilePath);
  end
  
  for loop = 1:length(cHeader)
    h = figure(loop);
    set(h, 'Visible', 'off')
    plot(TIMES(IndexOfTimePoints),timeSerie2bin(:,loop),'-ok','LineWidth',2);
    #axis([0 timePoints(end) 0 max(timeSerie2bin(:,loop))]);
    xlabel('time (min)','FontSize',14,'FontWeight','bold');
    ylabel('concentration (A.U)','FontSize',14,'FontWeight','bold');
    title(cHeader{loop});
    set(gca,'FontSize',14);
    fileName = cHeader{loop};
    saveas(gca,fullfile(figsFilePath,fileName),'png');  
  end
  close all;
  clear all; 
end