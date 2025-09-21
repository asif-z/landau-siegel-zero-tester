#ifndef PRESETS_H
#define PRESETS_H

// presets of optimized values of lambda and the corresponding phi and E values. Given some number of primes we want to sum, we have presets for the optimal choice of lambda

enum Preset {
    smallX1,
    mediumX1,
    largeX1,
    smallX2,
    mediumX2,
    largeX2,
    hugeX2
  };

void initializeLambda(enum Preset preset, arb_t lambda, arb_t phi, arb_t E, long prec);

#endif