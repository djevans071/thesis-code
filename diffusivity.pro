

function diffusivity, RF1, RF2, RF3, fps, mpp

  rcut = poscut(RF1)
  RVE  = bayes(RF1,rcut,fps,mpp) & $
    print, 'RF1 complete'
  RVE2 = bayes(RF2,rcut,fps/2.,mpp) & $
    print, 'RF2 complete'
  RVE3 = bayes(RF3,rcut,fps/3.,mpp) & $
    print, 'RF3 complete'
  
  r = reform(rve[0,*])
  v = reform(rve[1,*])
  v_err = reform(rve[2,*])      ;error in v
  posdist = reform(rve[5,*])
  
  eta1 = reform(rve[3,*])
    eta1_err = reform(rve[4,*])
  eta2 = reform(rve2[3,*])
    eta2_err = reform(rve2[4,*])
  eta3 = reform(rve3[3,*])
    eta3_err = reform(rve3[4,*])
    
  Dif1 = eta1*fps/2.
  Dif2 = (eta2-eta1)*fps/2.     ;diffusion relations
  Dif3 = (eta3-eta2)*fps/2.
  Dif4 = (eta3-eta1)*fps/4.
  
  Dif_err = fps*sqrt(eta1_err^2. + eta2_err^2. + $
    eta3_err^2.)/4.
  ;Dif_err = eta1_err*fps

return, [transpose(r), $
           transpose(Dif4), transpose(Dif_err), $
           transpose(v), transpose(v_err), $
           transpose(posdist)]

end           