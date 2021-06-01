#!/bin/sh

# conducts whole experiment putting density max-curve to stdout
# all experiment settings must be in environment variables (see errors
# for details)
# Rules: 5<->5 and 10<->16 with prob in PROB_MAIN

# TODO: check env args

# environment variables:
# PROB_MAIN
date
make purge 2>err
python system/gen_rules.py \
	5-5=$PROB_MAIN 10-16=$PROB_MAIN 16-10=$PROB_MAIN \
	5-10=$PROB_SEC 10-5=$PROB_SEC 16-5=$PROB_SEC \
	> rules.txt 2>err || exit 1

python system/gen_config_flatwave.py 2>err > config.c || exit 1

echo "EXPERIMENT: PROB_MAIN=$PROB_MAIN PROB_SEC=$PROB_SEC"
make run 2>err

cwd=$PWD

mkdir exp 2>err
cd density
python ../system/build_max.py *.xls > $cwd/exp/exp_${PROB_MAIN}_${PROB_SEC}.xls
