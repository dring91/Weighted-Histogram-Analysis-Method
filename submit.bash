#!/bin/bash
#PBS -N analyze
#PBS -l nodes=1:ppn=1,pmem=500mb,walltime=01:00:00
#PBS -d /home/daring/wavy/
#PBS -q p_rrig           
#PBS -j oe

. ~/py-3.6/bin/activate

#time python prog/calcCentersOfMass.py inputs/N100_S0_R3_A1_W2_P1_out.conf -g 1 -g 2 3 --format conf
#time python prog/undo_bondswap.py -i N100_S0_R3_A1_W8_P0_E0515_K0_D23.8_equil.conf -p 2249059.walnut-login1.internal/ -s .conf -l 100 --polymer 1 --unwrap
#time python prog/undo_bondswap.py -i conf.40000000 -p 2249060.walnut-login1.internal/ -s .40000000 -l 100 --polymer 1 --unwrap
#time python prog/undo_bondswap.py -i conf.70000000 -p 2249060.walnut-login1.internal/ -s .70000000 -l 100 --polymer 1 --unwrap
#time python prog/undo_bondswap.py -i conf.100000000 -p 2249060.walnut-login1.internal/ -s .100000000 -l 100 --polymer 1 --unwrap
#time python prog/thermo.py 2249059.walnut-login1.internal/ 2249060.walnut-login1.internal/ -l 26
#N=100; S=0; R=3; A=1; W=8; P=(0 1); E=0515; K=4000; D=(24.35 25.0); t=(40000000 70000000 100000000)
#dirs=(2250662 2250663 2250664 2250665 2250666 2250667)
#N=100; S=0; R=3; A=1; W=8; P=1; E=0515; K=4000; D=(24.95 25.0 25.05); t=(40000000 70000000 100000000)
#dirs=(2253606 2253607 2253608 2253609 2253610 2253611 2253612 2253613 2253614)
N=100; S=0; R=3; A=1; W=8; P=1; E=0515; K=4000; D=(25.05 25.0 24.95 24.9 24.85 24.8 24.75 24.7 24.65 24.6 24.55 24.5 24.45 24.4 24.35 24.3 24.25 24.2 24.15 24.1 24.05 24.0 23.95 23.9); t=(40000000 70000000 100000000)
dirs=(2253612 2253613 2253614 2253609 2253610 2253611 2253606 2253607 2253608 2254028 2254029 2254030 2254031 2254032 2254033 2254034 2254035 2254036 2254037 2254038 2254039 2254040 2254041 2254042 2254043 2254044 2254045 2254046 2254047 2254048 2254049 2254050 2254051 2254052 2254053 2254054 2254055 2254056 2254057 2254058 2254059 2254060 2254061 2254062 2254063 2254064 2254065 2254066 2254067 2254068 2254069 2254070 2254071 2254072 2254073 2254074 2254075 2254076 2254077 2254078 2254079 2254080 2254081 2254082 2254083 2254084 2254085 2254086 2254087 2254088 2254089 2254090)
#N=100; S=0; R=3; A=1; W=8; P=1; E=0515; K=4000; D=(25.05 25.0 24.9 24.85 24.8); t=(40000000 70000000 100000000)
#dirs=(2254028 2254029 2254030 2254031 2254032 2254033 2254034 2254035 2254036 2254037 2254038 2254039 2254040 2254041 2254042 2254043 2254044 2254045 2254046 2254047 2254048 2254049 2254050 2254051 2254052 2254053 2254054 2254055 2254056 2254057 2254058 2254059 2254060 2254061 2254062 2254063 2254064 2254065 2254066 2254067 2254068 2254069 2254070 2254071 2254072 2254073 2254074 2254075 2254076 2254077 2254078 2254079 2254080 2254081 2254082 2254083 2254084 2254085 2254086 2254087 2254088 2254089 2254090)
for i in {0..23}; do
  for j in 0 1 2; do
    flat=$i*3+$j
    dir[$flat]=${dirs[$flat]}.walnut-login1.internal
    file=N${N}_S${S}_R${R}_A${A}_W${W}_P${P}_E${E}_K${K}_D${D[$i]}_${t[$j]}.lammpstrj
    filenames[$flat]=${dir[$flat]}/coms.dat
    #time python prog/calcMeltHeight.py -i ${dir[$flat]}/$file -o coms.dat --polymer 1 --capillary 2 3 --dir ${dir[$flat]}/
  done
done
#time python prog/wham.py -i ${dir[*]} -o wham.dat -p ./ -k $K -d ${D[*]}
#time python prog/plot-ComSeparation.py ${filenames[*]} -l ${D[*]}

#time python prog/undo_bondswap.py -i N100_S0_R3_A1_W8_P1_out.conf -p inputs/ -s _N100_S0_R3_A1_W8_P1.conf -l 100 --polymer 1 --unwrap

dirs=(2271339 2271340 2271341 2271342 2271343); labels=(0 4 8 12 16)
dirs=(2273282 2273283 2273284 2273285 2273286); labels=(0 4 8 12 16)
dirs=(2274806 2274807 2274808 2274809 2274810); labels=(20 24 28 32 36)
dirs=(2275473 2275474 2275475 2275477); labels=(40 44 48 56)
dirs=(2273282 2273283 2273284 2273285 2273286 2274806 2274807 2274808 2274809 2274810 2275473 2275474 2275475 2275476 2275477); labels=(0 4 8 12 16 20 24 28 32 36 40 44 48 52 56); shifts=(0 -10 -24 -40 -59 -78 -99 -122 -148 -175 -175 -175 -175 -175 -175)
dirs=(2318723 2318724 2318725 2318726); N=100; S=(0 0 5 5); R=(3 3 6 6); A=(1.5 1 3 2); W=(4 16 4 16); P=1; E=0400
for i in {0..3}; do
  #files[$i]="${dirs[$i]}.walnut-login1.internal/plumed.out"
  #files[$i]="${dirs[$i]}.walnut-login1.internal/heights.dat"
  #files[$i]="${dirs[$i]}.walnut-login1.internal/N100_S0_R3_A1_W8_P1_E0515_M${labels[$i]}_biased.lammpstrj"
  files[$i]="${dirs[$i]}.walnut-login1.internal/N${N}_S${S[$i]}_R${R[$i]}_A${A[$i]}_W${W[$i]}_P${P}_E${E}_equil.lammpstrj"
  #echo "${files[$i]}"
  #time python scripts/calcMeltHeight.py -i ${files[$i]} -o heights.dat --dir ${dirs[$i]}.walnut-login1.internal/ --polymer 1 --capillary 2 3
done

time python scripts/simple-pmf.py -i ${files[0]} -p equilibrate/

#time python scripts/wham.py -i ${files[*]} -o wham -p indus/ -k 0.5 -x ${labels[*]} -g ${shifts[*]}
#time python scripts/plot-indus.py -i ${files[*]:10:15} -l ${labels[*]:10:15} -k 0.5
#time python scripts/uwham.py -i ${files[*]} -K 0.5 -N ${labels[*]} -s ${shifts[*]}
#time python scripts/plot-height_vs_number.py ${files[*]} -l ${labels[*]}
#time python scripts/shift_config.py 2249060.walnut-login1.internal/unshuffled.conf -o inputs/indus/shifted

#time python scripts/bond_lengths.py inputs/cylinders/N5_S0_R3_out.conf
#time python scripts/undo_bondswap.py -i N100_S0_R3_out.conf -p inputs/cylinders/ -s _N100_S0_R3.conf -l 100 --polymer 1 --unwrap
#time python scripts/shorten_chains.py -i inputs/cylinders/unshuffled_N100_S0_R3 -o N5_S0_R3_out -l 5 --polymer 1
