function O2 = partialPressureO2(time)
  
  if time <= 4*60
    O2 = 21;
  else
    O2 = 0.5;
  end
end