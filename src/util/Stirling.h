/*  PartitionTest, fast selection of the best fit partitioning scheme for
 *  multi-gene data sets.
 *  Copyright May 2013 by Diego Darriba
 *
 *  This program is free software; you may redistribute it and/or modify its
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation; either version 3 of the License, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *  for more details.
 *
 *  For any other inquiries send an Email to Diego Darriba
 *  ddarriba@udc.es
 */

/**
 * @file Stirling.h
 */

#ifndef STIRLING_H_
#define STIRLING_H_

#include <iterator>
#include <vector>
#include <iostream>
#include <memory>
#include <algorithm>
#include <math.h>
#include "Utilities.h"

namespace partest
{

  class Stirling
  {
  public:
  public:
    class iterator : public std::iterator<std::bidirectional_iterator_tag,
        const std::vector<unsigned> >
    {
    public:
      iterator (unsigned n, bool first = true);

      reference operator* () const
      {
        return kappa;
      }
      pointer operator-> () const
      {
        return &kappa;
      }

      unsigned size () const
      {
        return kappa.size ();
      }
      unsigned subsets () const
      {
        return M[size () - 1] + 1;
      }

      iterator& operator++ ();
      iterator& operator-- ();

      template<typename Elem>
        std::vector<std::vector<Elem> > *
        operator[] (const std::vector<Elem> &v) const;
//        std::vector<uint> *
//        operator[](const std::vector<uint> &v) const;

    protected:
      std::vector<unsigned> kappa, M;

      void integrityCheck ();
    };

    class iterator_k : public iterator
    {
    public:
      iterator_k (unsigned n, unsigned psize, bool first = true);

      // optimized version
      unsigned subsets () const
      {
        return psize;
      }

      iterator_k& operator++ ();
      iterator_k& operator-- ();

      long get_it_size ()
      {
        return it_size;
      }

    private:
      const unsigned psize;
      long it_size;

      void integrityCheck ();
    };

    static iterator_k * get_partitions (int k, t_schemesVector* &schemeVector,
                                        t_partitioningScheme* &ptr_elements);
  };

  extern std::ostream& operator<< (std::ostream& out, Stirling::iterator &it);

  template<typename Elem>
    std::vector<std::vector<Elem> > * Stirling::iterator::operator[] (
        const std::vector<Elem> &v) const
    {
      std::vector<std::vector<Elem> > *part =
          new std::vector<std::vector<Elem> > (subsets ());

      for (unsigned i = 0; i < size (); ++i)
        (*part)[kappa[i]].push_back (v[i]);

      return part;
    }

  template<typename Elem>
    std::ostream&
    operator<< (std::ostream &out, std::vector<Elem> &s)
    {
      out << '(';

      for (unsigned i = 0; i < s.size () - 1; ++i)
      {
        out << s[i] << ' ';
      }
      out << *(s.end () - 1) << ')';

      return out;
    }

  template<typename Elem>
    std::ostream&
    operator<< (std::ostream &out, std::vector<std::vector<Elem> > &part)
    {
      out << '(';

      for (unsigned i = 0; i < part.size () - 1; ++i)
      {
        out << part[i] << ' ';
      }
      out << *(part.end () - 1) << ')';

      return out;
    }

} /* namespace partest */
#endif /* STIRLING_H_ */
