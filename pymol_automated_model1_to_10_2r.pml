
load ./model_1.pdb
load ./model_2.pdb
load ./model_3.pdb
load ./model_4.pdb
load ./model_5.pdb
load ./model_6.pdb
load ./model_7.pdb
load ./model_8.pdb
load ./model_9.pdb
load ./model_10.pdb

set_view (\
     1.000000000,    0.000000000,    0.000000000,\
     0.000000000,    1.000000000,    0.000000000,\
     0.000000000,    0.000000000,    1.000000000,\
     0.000000790,   -0.000000045, -246.803405762,\
   -24.021657944,   -1.887244582,  -33.041385651,\
  -3452.230468750, 3945.836914062,  -20.000000000 )
 
set bg_rgb, white 
show surface
hide surface, lig_1.pdb
hide surface, lig_2.pdb
hide surface, lig_3.pdb
hide surface, lig_4.pdb
hide surface, lig_5.pdb
hide surface, lig_6.pdb
hide surface, lig_7.pdb
hide surface, lig_8.pdb
hide surface, lig_9.pdb
hide surface, lig_10.pdb

set cartoon_color, purple, lig_2.pdb
set cartoon_color, green, lig_3.pdb
set cartoon_color, yellow, lig_4.pdb
set cartoon_color, grey, lig_5.pdb
set cartoon_color, pink, lig_6.pdb
set cartoon_color, blue, lig_7.pdb
set cartoon_color, black, lig_8.pdb
set cartoon_color, teal, lig_9.pdb
set cartoon_color, red, lig_10.pdb

set transparency, 0.0
show cartoon, rec

###################
### Save a copy ###
###################
set antialias, 2
set hash_max, 220
set ray_shadows,0
png model1to10_r1.png, width=45cm, height=45cm, dpi=300, ray=1200,1200

set_view (\
    -0.699680567,   -0.353625149,    0.620799839,\
     0.036005698,    0.850358725,    0.524968684,\
    -0.713549316,    0.389661342,   -0.582249641,\
    -0.000023380,    0.000080887, -311.091522217,\
   -23.161323547,   -4.897098541,  -21.675985336,\
  -10184.497070312, 10806.663085938,  -20.000000000 )

png model1to10_r2.png, width=45cm, height=45cm, dpi=300, ray=1200,1200

quit

