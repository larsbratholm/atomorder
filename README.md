# atomorder
(Soon to be L'arsign(TM))
Maps atoms from reactants to products in a reaction and other secret cool stuff.

examples:

./main.py -r examples/reaction/R{2,1}.xyz -p examples/reaction/P.xyz --print-level 4 -m full --atomic-sybyl-weight 1 --bond-weight 50

./main.py -r examples/water_clusters/org_wat_ah_adf_agm_bbf_ava.xyz -p examples/water_clusters/org_wat_bn_ayk_bbf_adx_bvp.xyz --print-level 4 -m full --atomic-sybyl-weight 1 --bond-weight 50

./main.py -r examples/reaction2/R{1,2}.xyz -p examples/reaction2/P{1,2}.xyz --print-level 4 -m full --atomic-sybyl-weight 1 --bond-weight 0.01
