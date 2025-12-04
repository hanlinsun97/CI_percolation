#ifndef RNG_H
#define RNG_H

#include <stdint.h>

// Thread-safe RNG state (PCG algorithm - fast and high quality)
typedef struct {
  uint64_t state;
  uint64_t inc;
} rng_state_t;

// Initialize RNG state
static inline void rng_init(rng_state_t *rng, uint64_t seed) {
  rng->state = 0;
  rng->inc = (seed << 1) | 1;
  // Advance state
  rng->state = rng->state * 6364136223846793005ULL + rng->inc;
  rng->state += seed;
  rng->state = rng->state * 6364136223846793005ULL + rng->inc;
}

// Generate random 32-bit integer
static inline uint32_t rng_next(rng_state_t *rng) {
  uint64_t oldstate = rng->state;
  rng->state = oldstate * 6364136223846793005ULL + rng->inc;
  uint32_t xorshifted = ((oldstate >> 18) ^ oldstate) >> 27;
  uint32_t rot = oldstate >> 59;
  return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

// Generate random integer in range [0, bound)
static inline uint32_t rng_range(rng_state_t *rng, uint32_t bound) {
  uint32_t threshold = -bound % bound;
  for (;;) {
    uint32_t r = rng_next(rng);
    if (r >= threshold)
      return r % bound;
  }
}

// Generate random double in [0, 1)
static inline double rng_double(rng_state_t *rng) {
  return (double)rng_next(rng) / (double)0x100000000ULL;
}

#endif // RNG_H
