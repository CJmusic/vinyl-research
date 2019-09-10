HawfA = fdesign.audioweighting('WT,Class','A',1,44.1e3);
HawfITUR = fdesign.audioweighting('WT','ITUR4684',44.1e3);
Afilter = design(HawfA,'SystemObject',true);
ITURfilter = design(HawfITUR,'SystemObject',true);
hfvt = fvtool(Afilter,ITURfilter);
axis([0.1 22 -80 20]);
legend(hfvt,'A-weighting','ITU-R 468-4');