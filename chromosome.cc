/*
 * Implementation for Chromosome class
 */

#include <algorithm>
#include <cassert>
#include <random>

#include "chromosome.hh"

//////////////////////////////////////////////////////////////////////////////
// Generate a completely random permutation from a list of cities
Chromosome::Chromosome(const Cities* cities_ptr) : cities_ptr_(cities_ptr),
                           order_(random_permutation(cities_ptr->size())) {
  assert(is_valid());
  // Seed the random number generator with 6*9
  generator_ = std::default_random_engine(42);
}

//////////////////////////////////////////////////////////////////////////////
// Clean up as necessary
Chromosome::~Chromosome() {
  assert(is_valid());
}

//////////////////////////////////////////////////////////////////////////////
// Perform a single mutation on this chromosome
void Chromosome::mutate() {
  // if the length of the permutation changed over the lifetime of the
  // chromosome, then the distribution could not be static. but the header file
  // specifies that this will not happen.
  static std::uniform_int_distribution<unsigned> dist0(0, order_.size()-1);
  // pick a second random number not equal to the first by adding a nonzero but
  // otherwise random number to the first number. mod by order_.size() so that
  // it is still a valid index.
  static std::uniform_int_distribution<unsigned> dist1(1, order_.size()-1);
  unsigned p0 = dist0(generator_);
  unsigned p1 = (p0 + dist1(generator_)) % order_.size();
  std::swap(order_[p0], order_[p1]);
  assert(is_valid());
}

//////////////////////////////////////////////////////////////////////////////
// Return a pair of offsprings by recombining with another chromosome
// Note: this method allocates memory for the new offsprings
std::pair<Chromosome*, Chromosome*>
Chromosome::recombine(const Chromosome* other) {
  assert(is_valid());
  assert(other->is_valid());
  // Assert that the other chromosome points to the same cities as this one
  // does, because if it doesn't bad stuff is going to happen.
  assert(other->compare_cities_ptr(cities_ptr_));
  // Ultimately I decided to stick with this copypaste. putting "generate a pair
  // of nonequal indices" into a method would make the return value an int
  // array, and that would be uglier than just putting it inline. And making
  // dist0 and dist1 for the whole class requires they be initialized after
  // order does, and getting that to work was not straightforward.
  static std::uniform_int_distribution<unsigned> dist0(0, order_.size()-1);
  static std::uniform_int_distribution<unsigned> dist1(1, order_.size()-1);
  unsigned p0 = dist0(generator_);
  unsigned p1 = (p0 + dist1(generator_)) % order_.size();
  Chromosome* c1 = create_crossover_child(this, other, p0, p1);
  Chromosome* c2 = create_crossover_child(this, other, p0, p1);
  std::pair<Chromosome*, Chromosome*> retval(c1, c2);
  return retval;
}

//////////////////////////////////////////////////////////////////////////////
// For an ordered set of parents, return a child using the ordered crossover.
// The child will have the same values as p1 in the range [b,e),
// and all the other values in the same order as in p2.
Chromosome* Chromosome::create_crossover_child(const Chromosome* p1, const
                           Chromosome* p2, unsigned b, unsigned e) const {
  Chromosome* child = p1->clone();

  // We iterate over both parents separately, copying from parent1 if the
  // value is within [b,e) and from parent2 otherwise
  unsigned i = 0, j = 0;
  // cpplint didn't like the following code. i'm leaving it in because it was
  // provided.
  for ( ; i < p1->order_.size() && j < p2->order_.size(); ++i) {
    if (i >= b and i < e) {
      child->order_[i] = p1->order_[i];
    }
    else { // Increment j as long as its value is in the [b,e) range of p1
      while (p1->is_in_range(p2->order_[j], b, e)) {
        ++j;
      }
      assert(j < p2->order_.size());
      child->order_[i] = p2->order_[j];
      j++;
    }
  }

  assert(child->is_valid());
  return child;
}

// Return a positive fitness value, with higher numbers representing
// fitter solutions (shorter total-city traversal path).
double Chromosome::get_fitness() const {
  return cities_ptr_->total_path_distance(order_);
}

// A chromsome is valid if it has no repeated values in its permutation,
// as well as no indices above the range (length) of the chromosome.
// We implement this check with a sort, which is a bit inefficient, but simple
bool Chromosome::is_valid() const {
  // first, sanity check
  if (cities_ptr_->size() != order_.size()) {
    return false;
  }
  // use counting sort. if during the counting sort you reach a previously
  // counted city, return false. otherwise return true.
  std::vector<bool> city_already_in_perm(cities_ptr_->size(), false);
  for (auto i : order_) {
    if (city_already_in_perm[i]) {
      return false;
    }
    city_already_in_perm[i] = true;
  }
  return true;
}

// Find whether a certain value appears in a given range of the chromosome.
// Returns true if value is within the specified the range specified
// [begin, end) and false otherwise.
bool Chromosome::is_in_range(unsigned value, unsigned begin, unsigned end)
                                                                   const {
  if (begin >= end) {
    return false;
  }
  // let's just treat this like an array so we don't have to get two new
  // iterators and increment them up.
  for (unsigned i = begin; i < end; i++) {
    if (order_[i] == value) {
      return true;
    }
  }
  return false;
}
