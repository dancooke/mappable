// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include <vector>
#include <deque>
#include <string>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <random>
#include <utility>
#include <cassert>
#include <iostream>

#include "mappable/mappable_fwd.hpp"

using namespace mappable;

// Some types used in this example

// Mappable provdes two fundamental region types which all Mappable types must use.
//
// ContigRegion is the simpler of the two; it is just a start and end co-ordinate. Use this
// when you know, or don't need to know the contig mapping.
// 
// GenomicRegion is a contig name + a ContigRegion, and thus defines a complete genomic mapping.
// It therefore uses a few more bytes compared to ContigRegion.
//
// Note both ContigRegion and GenomicRegion use half open intervals [begin, end).
//
// All Mappable methods, algorithms, and containers work the same with both ContigRegion and GenomicRegion,
// but you can't directly compare one with the other, which means you can't mix the two types in the same
// collection.

//
// Example using ContigRegion
//
struct Read : public Mappable<Read>
{
    ContigRegion region;
    unsigned quality;
    const ContigRegion& mapped_region() const noexcept { return region; } // required for Mappable
    Read(unsigned begin, unsigned end, unsigned quality) : region {begin, end}, quality {quality} {}
};

std::ostream& operator<<(std::ostream& os, const Read& read)
{
    os << "(" << read.region << " " << read.quality << ")";
    return os;
}

using MappingQualitySet = MappableFlatMultiSet<Read>;

//
// Example using GenomicRegion
//
struct Allele : public Mappable<Allele>, Comparable<Allele>
{
    GenomicRegion region;
    std::string sequence;
    const GenomicRegion& mapped_region() const noexcept { return region; } // required for Mappable
    template <typename R, typename S> Allele(R&& region, S&& sequence)
    : region {std::forward<R>(region)}, sequence {std::forward<S>(sequence)} {}
};

// We define operators == and < for Allele as we want to make sure alleles with
// different sequences are distinct; the default comparitors provided by Mappable
// only compare regions.
// Note is important the any user defined operator< comparitor for Mappable types
// satisfies the properties of Mappables default operator<. That is, any Mappable type must 
// first sort by the types region; only objects deemed equivilant by the default operator<
// can be futher sorted.
bool operator==(const Allele& lhs, const Allele& rhs) noexcept
{
    return lhs.region == rhs.region && lhs.sequence == rhs.sequence;
}
bool operator<(const Allele& lhs, const Allele& rhs) noexcept
{
    return lhs.region == rhs.region ? lhs.sequence < rhs.sequence : lhs.region < rhs.region;
}

std::ostream& operator<<(std::ostream& os, const Allele& allele)
{
    os << "{" << allele.region << " " << allele.sequence << "}";
    return os;
}

using AlleleSet = MappableFlatSet<Allele>;

struct Genotype : public Mappable<Genotype>
{
    Allele first, second;
    const GenomicRegion& mapped_region() const noexcept { return first.mapped_region(); } // required for Mappable
    template <typename A1, typename A2> Genotype(A1&& first, A2&& second)
    : first {std::forward<A1>(first)}, second {std::forward<A2>(second)} {}
    template <typename R, typename S1, typename S2>
    Genotype(R&& region, S1&& sequence1, S2&& sequence2)
    : first {region, std::forward<S1>(sequence1)}
    , second {std::forward<R>(region), std::forward<S2>(sequence2)}
    {} 
};

std::ostream& operator<<(std::ostream& os, const Genotype& genotype)
{
    os << "[" << genotype.first << " " << genotype.second << "]";
    return os;
}

using AlleleReference = MappableReferenceWrapper<Allele>;

std::ostream& operator<<(std::ostream& os, const AlleleReference& allele)
{
    os << allele.get();
    return os;
}

// Some utility print methods

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
{
    std::copy(std::cbegin(v), std::cend(v), std::ostream_iterator<T> {os, " "});
    return os;
}

template <typename It>
std::ostream& operator<<(std::ostream& os, const OverlapRange<It>& range)
{
    using T = typename std::iterator_traits<It>::value_type;
    std::copy(std::cbegin(range), std::cend(range), std::ostream_iterator<T> {os, " "});
    return os;
}

template <typename It>
std::ostream& operator<<(std::ostream& os, const ContainedRange<It>& range)
{
    using T = typename std::iterator_traits<It>::value_type;
    std::copy(std::cbegin(range), std::cend(range), std::ostream_iterator<T> {os, " "});
    return os;
}

//
// This example shows basic usage with standard C++ containers
//
void std_container_example()
{
    using std::cbegin; using std::cend;
    std::cout << "Running std_container_example..." << std::endl;
    
    const ContigRegion r1 {0, 2}, r2 {1, 2}, r3 {1, 3}, r4 {2, 5}, r5 {3, 4}, r6 {4, 5};
    
    // Mappable works with standard containers
    const std::vector<ContigRegion> regions_vector {r1, r2, r3, r4, r5, r6};
    const std::deque<ContigRegion> regions_deque {r1, r2, r3, r4, r5, r6};
    
    // Mappable algorithms require sorted ranges
    assert(std::is_sorted(cbegin(regions_vector), cend(regions_vector)));
    assert(std::is_sorted(cbegin(regions_deque), cend(regions_deque)));
    
    std::cout << "regions_vector: " << regions_vector << std::endl;
    std::cout << "regions_deque: " << regions_vector << std::endl;
    
    const ContigRegion test {2, 4};
    
    const auto overlapped_vector = overlap_range(regions_vector, test);
    std::cout << "There are " << count_overlapped(regions_vector, test)
              << " regions in regions_vector overlapped with " << test << ": "
              << overlapped_vector << std::endl;
    const auto overlapped_deque = overlap_range(regions_deque, test);
    std::cout << "There are " << count_overlapped(regions_deque, test)
              << " regions in regions_deque overlapped with " << test << ": "
              << overlapped_deque << std::endl;
    
    // Note: You wouldn't use count_overlapped like this if you already had
    // the result for overlap_range; just call size(overlapped_vector) instead.
}

void mappable_multiset_example()
{
    std::cout << "Running mappable_set_example..." << std::endl;
    
    constexpr std::size_t num_reads {100'000};
    constexpr unsigned read_size {150};
    constexpr ContigRegion::Size contig_size {2'000'000};
    
    static std::default_random_engine gen {};
    std::uniform_int_distribution<ContigRegion::Size> position_dist {0, contig_size};
    std::uniform_int_distribution<unsigned> quality_dist {0, 100};
    
    MappingQualitySet reads {};
    
    // First let's fill the container with some randomly generated data.
    // We don't need to worry about sorting as MappableFlatMultiSet will
    // handle this.
    std::cout << "Inserting " << num_reads << " randomly generated 'reads' into the MappableFlatMultiSet" << std::endl;
    reads.reserve(num_reads);
    std::generate_n(std::inserter(reads, std::begin(reads)), num_reads, [&] () -> Read { 
        auto begin = position_dist(gen);
        auto quality = quality_dist(gen);
        return {begin, std::min(begin + read_size, contig_size), quality};
    });
    // Note this is a very inefficient way to input data into a MappableFlatMultiSet!
    // A much better way is to insert pre-sorted data.
    
    // We can now query the MappableFlatMultiSet using any mappable algorithm.
    const ContigRegion test_region {contig_size / 2 - contig_size / 4, contig_size / 2 + contig_size / 4};
    std::cout << "There are " << count_overlapped(reads, test_region) << " reads overlapping with " << test_region << std::endl;
    std::cout << "There are " << count_contained(reads, test_region) << " reads contained within " << test_region << std::endl;
    
    // Let's do something a bit more interesting; calculate the average read quality
    // of reads contained in a region.
    const auto contained = contained_range(reads, test_region);
    const auto quality_sum = std::accumulate(std::cbegin(contained), std::cend(contained), 0.0,
                                             [] (auto curr, const auto& read) { return curr + read.quality; });
    const auto mean_quality = quality_sum / size(contained);
    std::cout << "The mean quality of reads contained in " << test_region << " is " << mean_quality << std::endl;
    
    // Which positions on are contig have read coverage?
    const auto covered_regions = extract_covered_regions(reads);
    const auto covered_length = sum_region_sizes(covered_regions);
    const auto covered_fraction = 100 * static_cast<double>(covered_length) / contig_size;
    std::cout << "The generated reads covered " << covered_fraction << "% of the contig: "
              << covered_regions << std::endl;
    
    // We can also easily get any 'intervening' regions between the covered regions.
    const auto intervening_regions = extract_intervening_regions(covered_regions, ContigRegion {0, contig_size});
    const auto intervening_length = sum_region_sizes(intervening_regions);
    std::cout << "covered_length + intervening_length = " << covered_length + intervening_length << std::endl;
    
    // We can also efficiently look at the coverage in a region
    const auto depths = calculate_positional_coverage(reads, reads[num_reads / 2]);
    std::cout << "The depth of each position in the region " << mapped_region(reads[num_reads / 2])
              << " is " << depths << std::endl;
    
    // It's easy to remove items from a MappableFlatMultiSet...
    reads.erase_overlapped(test_region);
    std::cout << "After removing reads overlapped with " << test_region 
              << " there are " << reads.size() << " reads remaining" << std::endl;
}

void mappable_set_example()
{
    std::cout << "Running mappable_multiset_example..." << std::endl;
    
    // MappableFlatSet is similar to MappableFlatMultiSet, but duplicates are not allowed.
    // There is also no reserve method for MappableFlatSet.
    
    AlleleSet alleles {};
    alleles.emplace(GenomicRegion {"X", 100, 101}, "A");
    alleles.emplace(GenomicRegion {"X", 101, 102}, "C");
    alleles.emplace(GenomicRegion {"X", 101, 102}, "AC");
    assert(alleles.size() == 3);
    alleles.emplace(GenomicRegion {"X", 100, 101}, "A");
    assert(alleles.size() == 3);
}

void complex_usage_example()
{
    std::cout << "Running complex_usage_example..." << std::endl;
    
    // One of the most powerful of the Mappable library is the ability to mix
    // different Mappable types (assuming same underlying region type - see intro).
    // For example, we can suppose we have some genotypes
    
    const std::vector<Genotype> genotypes {
        Genotype {GenomicRegion {"X", 0, 2}, "CC", "CC"},
        Genotype {GenomicRegion {"X", 2, 3}, "A", "C"},
        Genotype {GenomicRegion {"X", 3, 4}, "G", "T"},
        Genotype {GenomicRegion {"X", 4, 9}, "ACGT", ""},
        Genotype {GenomicRegion {"X", 9, 10}, "A", "A"}
    };
    
    assert(std::is_sorted(std::cbegin(genotypes), std::cend(genotypes)));
    
    // What alleles are present?
    AlleleSet alleles {};
    for (const auto& genotype : genotypes) {
        alleles.insert(genotype.first);
        alleles.insert(genotype.second);
    }
    
    // Which alleles are overlapped with each genotype?
    for (const auto& genotype : genotypes) {
        auto overlapped = overlap_range(alleles, genotype);
        std::cout << genotype << ": " << overlapped << std::endl;
    }
    // Note we didn't need to pass a region to overlap_range, we just passed
    // the genotype directly. Because both Allele and Genotype are mappable,
    // we can easily compare them in region space.
}

void mappable_reference_wrapper_example()
{
    std::cout << "Running complex_usage_example..." << std::endl;
    
    // Sometimes you will want to store a collection of references to Mappable objects.
    // As standard C++ solutions to this (e.g. pointers or std::reference_wrapper) are not
    // themselves Mappable, you wouldn't be able to use any mappable algorithms on this collection.
    // MappableReferenceWrapper provides a solution to this.
    
    std::vector<Allele> alleles {};
    alleles.emplace_back(GenomicRegion {"X", 100, 101}, "A");
    alleles.emplace_back(GenomicRegion {"X", 101, 104}, "GTC");
    alleles.emplace_back(GenomicRegion {"X", 102, 103}, "T");
    
    std::cout << "alleles: " << alleles << std::endl;
    
    std::vector<AlleleReference> allele_refs {};
    allele_refs.push_back(alleles[0]);
    allele_refs.push_back(alleles[2]);
    
    std::cout << "allele_refs: " << alleles << std::endl;
    
    assert(std::is_sorted(std::cbegin(allele_refs), std::cend(allele_refs)));
    
    std::cout << "There is " << count_overlapped(allele_refs, alleles[1])
              << " allele ref overlapped with " << alleles[1] << ": ";
    const auto overlapped = overlap_range(allele_refs, alleles[1]);
    std::cout << overlapped << std::endl;
}

int main()
{
    std_container_example();
    mappable_multiset_example();
    mappable_set_example();
    complex_usage_example();
    mappable_reference_wrapper_example();
}
