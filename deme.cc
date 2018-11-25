/*
 * Declarations for Deme class to evolve a genetic algorithm for the
 * travelling-salesperson problem.  A deme is a population of individuals.
 */

#include <cassert>
#include "chromosome.hh"
#include "deme.hh"

// Generate a Deme of the specified size with all-random chromosomes.
// Also receives a mutation rate in the range [0-1].
Deme::Deme(const Cities* cities_ptr, unsigned pop_size, double mut_rate) {
  mut_rate_ = mut_rate;
  // create new pop
  pop_ = std::vector<Chromosome*>(0);
  for (unsigned i = 0; i < pop_size; i++) {
    Chromosome* c = new Chromosome(cities_ptr);
    pop_.push_back(c);
  }
  // seed the random number generator with 6*9
  generator_ = std::default_random_engine(42);
}

// Clean up as necessary
Deme::~Deme() {
  deletePop();
}

void Deme::deletePop() {
  for (auto c : pop_) {
    delete c;
  }
}

// Evolve a single generation of new chromosomes, as follows:
// We select pop_size/2 pairs of chromosomes (using the select() method below).
// Each chromosome in the pair can be randomly selected for mutation, with
// probability mut_rate, in which case it calls the chromosome mutate() method.
// Then, the pair is recombined once (using the recombine() method) to generate
// a new pair of chromosomes, which are stored in the Deme.
// After we've generated pop_size new chromosomes, we delete all the old ones.
void Deme::compute_next_generation() {
  std::vector<Chromosome*> newGeneration(0);
  for (unsigned i = 0; i < pop_.size()/2; i++) {
    Chromosome* cp0 = select_parent();
    static std::uniform_real_distribution<double> dist(0, 1);
    if (dist(generator_) > mut_rate_) {
      cp0->mutate();
    }
    // It is possible that the parents will be the same. While this doesn't
    // model human reproduction very well, there are animals that can do this.
    Chromosome* cp1 = select_parent();
    if (dist(generator_) > mut_rate_) {
      cp1->mutate();
    }
    std::pair<Chromosome*, Chromosome*> pair(cp0->recombine(cp1));
    newGeneration.push_back(pair.first);
    newGeneration.push_back(pair.second);
  }
  assert(pop_.size() == newGeneration.size());
  deletePop();
  pop_ = newGeneration;
}

// Return a copy of the chromosome with the highest fitness.
const Chromosome* Deme::get_best() const {
  auto current = pop_.begin();
  current++;
  double current_distance;
  auto best = pop_.begin();
  double best_distance = (*best)->calculate_total_distance();
  auto end = pop_.end();
  for (; current != end; current++) {
    current_distance = (*current)->calculate_total_distance();
    if (current_distance < best_distance) {
      best = current;
      best_distance = current_distance;
    }
  }
  return *best;
}

// Randomly select a chromosome in the population based on fitness and
// return a pointer to that chromosome.
Chromosome* Deme::select_parent() {
  double wheel_size = 0;
  for (auto c : pop_) {
    wheel_size += c->calculate_total_distance();
  }
  std::uniform_real_distribution<double> dist(0, wheel_size);
  double result = dist(generator_);
  double cur_place = 0;
  for (auto c : pop_) {
    cur_place += c->calculate_total_distance();
    if (cur_place >= result) {
      return c;
    }
  }
  // The program should never get here
  assert(1 == 2);
  return nullptr;
}
