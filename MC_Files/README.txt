## README file to describe sequence to generate/reconstruct/analyze neutral kaon decays with the flexible MC

-----------------------------------------
---------- Generate MC files ------------
-----------------------------------------
1. Modify the macro file to generate the desired decay type in the decay region. Use the K2pipi.mac as an example.
2. Create a flexible MC geometry file in the  $NA62TOOLSSOURCE/Conditions/MC/BeamlineConfiguration/ directory (e.g. BeamlineDetectors_KS_KL_test.dat)
3. Once you have completed the two above steps run the following sequence of commands:
   - cd NA62MC
   - source scripts/env.sh
   - NA62MC -m macros/K2pipi.mac -n 10000 -r 9001 -s $seed -o prod_xxxdecay_xxxk.root --autoflush 70


-----------------------------------------
--------- Reconstruct MC files ----------
-----------------------------------------

1. Modify $NA62RECOSOURCE/config/NA62Reconstruction.MC.StdRecoSequence.conf to only reconstruct the needed detectors
2. Oncle you have made the necessary changes run the following sequence of commands:
   - cd $NA62RECOSOURCE
   - source scripts/env.sh
   - NA62Reco -c config/NA62Reconstruction.MC.NoOverlay.conf -i prod_xxxdecay_xxxk.root -e1 -o reco_xxxdecay_xxxk.root --autoflush 70
