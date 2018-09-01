#ifndef __VECTOR_SORTING_H__
#define __VECTOR_SORTING_H__

#include <algorithm> // std::count
#include <assert.h>
#include <vector>

std::vector< std::vector<int>> all_permutations_00p( std::vector<int> n0)
{
  std::vector< std::vector<int>> np;
  std::vector<int> n1 = n0;
  n1[2] = -n1[2];
  std::sort( n1.begin(), n1.end());
  np.push_back(n0); np.push_back(n1);
  while (std::next_permutation( n0.begin(), n0.end())) { np.push_back(n0); };
  while (std::next_permutation( n1.begin(), n1.end())) { np.push_back(n1); };
  return np;
}

std::vector< std::vector<int>> all_permutations_0pp( std::vector<int> n0)
{
  std::vector< std::vector<int>> np;
  std::vector<int> n1 = n0;
  std::vector<int> n2 = n0;
  n1[2] = -n1[2]; n2[1] = -n2[1]; n2[2] = -n2[2];
  std::sort( n1.begin(), n1.end()); std::sort( n2.begin(), n2.end());
  np.push_back(n0); np.push_back(n1); np.push_back(n2);
  while (std::next_permutation( n0.begin(), n0.end())) { np.push_back(n0); };
  while (std::next_permutation( n1.begin(), n1.end())) { np.push_back(n1); };
  while (std::next_permutation( n2.begin(), n2.end())) { np.push_back(n2); };
  return np;
}

std::vector< std::vector<int>> all_permutations_ppp( std::vector<int> n0)
{
  std::vector< std::vector<int>> np0,np;
  std::vector<int> n1 = n0;
  n1[2] = -n1[2];
  std::sort( n1.begin(), n1.end());
  np0.push_back(n0); np0.push_back(n1);
  while (std::next_permutation( n1.begin(), n1.end())) { np0.push_back(n1); };
  for (auto nx = np0.begin(); nx != np0.end(); nx++) {
   n1 = *nx; n1[0] = -n1[0]; n1[1] = -n1[1]; n1[2] = -n1[2];
   np.push_back(*nx); np.push_back(n1); };
  return np;
}

std::vector< std::vector<int>> all_permutations_0pq( std::vector<int> n0)
{
  std::vector< std::vector<int>> np0,np;
  std::vector<int> n1 = n0;
  n1[2] = -n1[2];
  std::sort( n1.begin(), n1.end());
  np0.push_back(n0); np0.push_back(n1);
  while (std::next_permutation( n0.begin(), n0.end())) { np0.push_back(n0); };
  while (std::next_permutation( n1.begin(), n1.end())) { np0.push_back(n1); };
  for (auto nx = np0.begin(); nx != np0.end(); nx++) {
   n1 = *nx; n1[0] = -n1[0]; n1[1] = -n1[1]; n1[2] = -n1[2];
   np.push_back(*nx); np.push_back(n1); };
  return np;
}

std::vector< std::vector<int>> all_permutations_ppq( std::vector<int> n0)
{
  std::vector< std::vector<int>> np0,np;
  std::vector<int> n1 = n0;
  std::vector<int> n2 = n0;
  n1[1] = -n1[1]; n2[2] = -n2[2];
  std::sort( n1.begin(), n1.end()); std::sort( n2.begin(), n2.end());
  np0.push_back(n0); np0.push_back(n1); np0.push_back(n2);
  while (std::next_permutation( n0.begin(), n0.end())) { np0.push_back(n0); };
  while (std::next_permutation( n1.begin(), n1.end())) { np0.push_back(n1); };
  while (std::next_permutation( n2.begin(), n2.end())) { np0.push_back(n2); };
  for (auto nx = np0.begin(); nx != np0.end(); nx++) {
   n1 = *nx; n1[0] = -n1[0]; n1[1] = -n1[1]; n1[2] = -n1[2];
   np.push_back(*nx); np.push_back(n1); };
  return np;
}

std::vector< std::vector<int>> all_permutations_pqq( std::vector<int> n0)
{
  std::vector< std::vector<int>> np0,np;
  std::vector<int> n1 = n0;
  std::vector<int> n2 = n0;
  n1[0] = -n1[0]; n2[2] = -n2[2];
  std::sort( n1.begin(), n1.end()); std::sort( n2.begin(), n2.end());
  np0.push_back(n0); np0.push_back(n1); np0.push_back(n2);
  while (std::next_permutation( n0.begin(), n0.end())) { np0.push_back(n0); };
  while (std::next_permutation( n1.begin(), n1.end())) { np0.push_back(n1); };
  while (std::next_permutation( n2.begin(), n2.end())) { np0.push_back(n2); };
  for (auto nx = np0.begin(); nx != np0.end(); nx++) {
   n1 = *nx; n1[0] = -n1[0]; n1[1] = -n1[1]; n1[2] = -n1[2];
   np.push_back(*nx); np.push_back(n1); };
  return np;
}

std::vector< std::vector<int>> all_permutations_pqr( std::vector<int> n0)
{
  std::vector< std::vector<int>> np0,np;
  std::vector<int> n1 = n0;
  std::vector<int> n2 = n0;
  std::vector<int> n3 = n0;
  n1[1] = -n1[1]; n2[2] = -n2[2]; n3[1] = -n3[1]; n3[2] = -n3[2];
  std::sort( n1.begin(), n1.end()); std::sort( n2.begin(), n2.end());
  std::sort( n3.begin(), n3.end());
  np0.push_back(n0); np0.push_back(n1); np0.push_back(n2); np0.push_back(n3);
  while (std::next_permutation( n0.begin(), n0.end())) { np0.push_back(n0); };
  while (std::next_permutation( n1.begin(), n1.end())) { np0.push_back(n1); };
  while (std::next_permutation( n2.begin(), n2.end())) { np0.push_back(n2); };
  while (std::next_permutation( n3.begin(), n3.end())) { np0.push_back(n3); };
  for (auto nx = np0.begin(); nx != np0.end(); nx++) {
   n1 = *nx; n1[0] = -n1[0]; n1[1] = -n1[1]; n1[2] = -n1[2];
   np.push_back(*nx); np.push_back(n1); };
  return np;
}

// from a general starting vector, add all permutations to a vector
std::vector< std::vector<int>> all_permutations( std::vector<int> n0)
{
  assert( n0.size() == 3); // only 3-vectors

  // count number of zeros
  int zeroCount = std::count( n0.begin(), n0.end(), 0);
  assert( zeroCount < 3); // must be nonzero vector

  // count number of unique, taking into account signs
  std::vector<int> n0a;
  for (auto nx = n0.begin(); nx != n0.end(); nx++) { n0a.push_back( abs(*nx)); }
  std::sort( n0a.begin(), n0a.end());
  int uniqueCount = 1;
  if (n0a[0] != n0a[1]) { uniqueCount++; }
  if (n0a[1] != n0a[2]) { uniqueCount++; }

  // check for ppq vs pqq case
  bool sameFront = ( n0a[0] == n0a[1] );

  // decide which class it belongs to
  if      (zeroCount == 2) { return all_permutations_00p( n0a); }
  else if (zeroCount == 1) {
    if      (uniqueCount == 2) { return all_permutations_0pp( n0a); }
    else if (uniqueCount == 3) { return all_permutations_0pq( n0a); }
  }
  else if (zeroCount == 0) {
    if      (uniqueCount == 3) { return all_permutations_pqr( n0a); }
    else if (uniqueCount == 2) {
      if (sameFront) { return all_permutations_ppq( n0a); }
      else           { return all_permutations_pqq( n0a); }
    }
    else if (uniqueCount == 1) { return all_permutations_ppp( n0a); }
  }
  assert(0); // failed to find a class for the momentum
}

int norm2( std::vector<int> vec) { return vec[0]*vec[0] +vec[1]*vec[1] +vec[2]*vec[2]; }

// figure out all of the permutations that have a given p^2
std::vector< std::vector<int>> all_combos( int p2)
{
  double p = sqrt(p2);
  std::vector< std::vector<int>> np;
  for (int nx = ceil( p); nx >= floor( p/sqrt(3.) ); nx--) {
    for (int ny = ceil( sqrt(p2 -nx*nx)); ny >= floor( sqrt(p2 -nx*nx)/sqrt(2.) ); ny--) {
      int nz = floor( sqrt(p2 -nx*nx -ny*ny));
      if (ny > nx) { continue; }
      if (nz > ny) { continue; }
      if ( norm2({nx,ny,nz}) == p2 ) { np.push_back({nx,ny,nz}); }
    }
  }
  return np;
}

#endif
