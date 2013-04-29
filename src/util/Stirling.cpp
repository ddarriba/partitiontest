/*
 * Stirling.cpp
 *
 *  Created on: Mar 20, 2013
 *      Author: diego
 */

#include "Stirling.h"
#include <stdexcept>

namespace partest {

Stirling::iterator::iterator(unsigned n, bool first) :
    kappa(n), M(n) {
  if ( n <= 0 ) throw std::invalid_argument("Stirling::iterator: n must be positive");

  if ( !first ) for ( unsigned i = 1; i < n; ++i ) {
    kappa[i] = i;
    M[i] = i;
  }

  // integrityCheck();
}

Stirling::iterator&
Stirling::iterator::operator++() {
  const unsigned n = size();

  for ( unsigned i = n - 1; i > 0; --i )
    if ( kappa[i] <= M[i - 1] ) {
      ++kappa[i];

      const unsigned new_max = std::max(M[i], kappa[i]);
      M[i] = new_max;
      for ( unsigned j = i + 1; j < n; ++j ) {
        kappa[j] = 0;
        M[j] = new_max;
      }

      // integrityCheck();
      return *this;
    }

  throw std::overflow_error("no following partitions");
}

Stirling::iterator&
Stirling::iterator::operator--() {
  const unsigned n = size();

  for ( unsigned i = n - 1; i > 0; --i )
    if ( kappa[i] > 0 ) {
      --kappa[i];

      const unsigned m_i = M[i - 1];
      M[i] = m_i;
      for ( unsigned j = i + 1; j < n; ++j ) {
          const unsigned new_max = m_i + j - i;
          kappa[j] = new_max;
          M[j] = new_max;
        }

        // integrityCheck();
        return *this;
      }

    throw std::underflow_error("no preceding partitions");
  }

  void Stirling::iterator::integrityCheck() {
    const unsigned n = size();
    unsigned max = kappa[0];

    for ( unsigned i = 0; i < n; ++i ) {
      max = std::max(max, kappa[i]);

      if ( max != M[i] ) throw std::domain_error("integrity check failed");
    }
  }

  std::ostream&
  operator<<(std::ostream& out, Stirling::iterator &it) {
    out << '(';

    if ( it->size() > 1 )
      std::copy(it->begin(), it->end() - 1, std::ostream_iterator<unsigned>(out, " "));

    out << *(it->end() - 1) << ')';

    return out;
  }

  Stirling::iterator_k::iterator_k(unsigned n, unsigned psize, bool first) :
      iterator(n, first), psize(psize) {
    if ( psize <= 0 || psize > n )
      throw std::invalid_argument("Stirling::iterator_k: psize must be in [1..n]");

    if ( first ) {
      const unsigned offset = n - psize;
      for ( unsigned i = offset + 1; i < n; ++i ) {
        kappa[i] = i - offset;
        M[i] = i - offset;
      }
    } else {
      std::fill(kappa.begin() + psize, kappa.end(), psize - 1);
      std::fill(M.begin() + psize, M.end(), psize - 1);
    }

    it_size = 0;
    for ( unsigned int i = 0; i <= psize; i++ ) {
      it_size += pow(-1, psize - i) * Utilities::combinatorial(psize, i) * pow(i, n);
    }
    it_size /= Utilities::factorial(psize);
  }

  Stirling::iterator_k&
  Stirling::iterator_k::operator++() {
    const unsigned n = size();

    for ( unsigned i = n - 1; i > 0; --i )
      if ( kappa[i] < psize - 1 && kappa[i] <= M[i - 1] ) {
        ++kappa[i];

        const unsigned new_max = std::max(M[i], kappa[i]);
        M[i] = new_max;

        for ( unsigned j = i + 1; j <= n - (psize - new_max); ++j ) {
          kappa[j] = 0;
          M[j] = new_max;
        }

        for ( unsigned j = n - (psize - new_max) + 1; j < n; ++j ) {
          const unsigned new_max = psize - (n - j);
          kappa[j] = new_max;
          M[j] = new_max;
        }

        // integrityCheck();
        return *this;
      }

    throw std::overflow_error("no following k-partitions");
  }

  Stirling::iterator_k&
  Stirling::iterator_k::operator--() {
    const unsigned n = size();

    for ( unsigned i = n - 1; i > 0; --i )
      if ( kappa[i] > 0 && psize - M[i - 1] <= n - i ) {
        --kappa[i];

        const unsigned m_i = M[i - 1];
        M[i] = m_i;

        for ( unsigned j = i + 1; j < i + (psize - m_i); ++j ) {
          const unsigned new_max = m_i + j - i;
          kappa[j] = new_max;
          M[j] = new_max;
        }

        for ( unsigned j = i + (psize - m_i); j < n; ++j ) {
          kappa[j] = psize - 1;
          M[j] = psize - 1;
        }
        // integrityCheck();
        return *this;
      }

    throw std::underflow_error("no preceding k-partitions");
  }

  void Stirling::iterator_k::integrityCheck() {
    iterator::integrityCheck();

    if ( subsets() != iterator::subsets() )
      throw std::domain_error("integrity check 2 failed");
  }

  Stirling::iterator_k * Stirling::get_partitions(int k, t_schemesVector* &schemeVector,
      t_partitioningScheme* &elements) {
	  Stirling::iterator_k * it;
    if ( !elements ) {
      elements = new t_partitioningScheme;
    } else {
      elements->clear();
    }
    if ( !schemeVector ) {
      schemeVector = new t_schemesVector;
    }
    if ( schemeVector->size() > 0 ) {
      //build t_partition_elements
      for ( unsigned int i = 0; i < schemeVector->size(); i++ ) {
        unsigned int p_index = 0;
        for ( unsigned int j = 0; j < (schemeVector->at(i)).size(); j++ ) {
          p_index += (schemeVector->at(i))[j];
        }
        elements->push_back(p_index);
      }
      it = new Stirling::iterator_k(schemeVector->size(), k);
    } else {
      unsigned int p_index = 1;
      for ( unsigned int i = 0; i < k; i++ ) {
        elements->push_back(p_index);
        p_index *= 2;
      }
      it = new Stirling::iterator_k(k, k);
    }
    return it;
  }

} /* namespace partest */
