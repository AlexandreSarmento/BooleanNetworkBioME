function idxs = timepointidx(ListOfTimePoints,TIMES)
  
  [~,idx]=min(abs(TIMES-ListOfTimePoints));
  minVal=TIMES(idx);
  IdentifyIndex = ismember(TIMES,minVal);
  idxs = find(IdentifyIndex);
  
end