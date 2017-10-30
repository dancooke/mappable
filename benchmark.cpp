// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include <vector>
#include <string>
#include <algorithm>
#include <iterator>
#include <utility>
#include <random>
#include <chrono>
#include <iostream>

#include "mappable/mappable_fwd.hpp"
#include "IntervalTree.h" // https://github.com/ekg/intervaltree

using namespace mappable;

struct Read : public Mappable<Read>, Comparable<Read>
{
    ContigRegion region;
    std::string sequence;
    const auto& mapped_region() const noexcept { return region; } // required for Mappable
    template <typename R, typename S> Read(R&& region, S&& sequence)
    : region {std::forward<R>(region)}, sequence {std::forward<S>(sequence)} {}
};

bool operator==(const Read& lhs, const Read& rhs) noexcept
{
    return lhs.region == rhs.region && lhs.sequence == rhs.sequence;
}
bool operator<(const Read& lhs, const Read& rhs) noexcept
{
    return lhs.region == rhs.region ? lhs.sequence < rhs.sequence : lhs.region < rhs.region;
}

using ReadSet = MappableFlatSet<Read>;

using IntervalTreeRead = Interval<std::string>;

static constexpr std::array<char, 4> dnaBases {'A', 'C', 'G', 'T'};

std::string generate_sequence(const unsigned size)
{
    static std::default_random_engine gen {};
    std::uniform_int_distribution<char> dist {0, 4};
    std::string result(size, '.');
    std::generate_n(std::begin(result), size, [&] () { return dnaBases[dist(gen)]; });
    return result;
}

auto generate_mappable_data(const std::size_t num_reads, const unsigned read_size, const unsigned contig_size)
{
    static std::default_random_engine gen {};
    std::uniform_int_distribution<ContigRegion::Size> position_dist {0, contig_size};
    std::vector<Read> result {};
    result.reserve(num_reads);
    std::generate_n(std::back_inserter(result), num_reads, [&] () -> Read { 
        auto begin = position_dist(gen);
        return {ContigRegion {begin, std::min(begin + read_size, contig_size)}, generate_sequence(read_size)};
    });
    return result;
}

auto generate_interval_tree_data(const std::size_t num_reads, const unsigned read_size, const unsigned contig_size)
{
    static std::default_random_engine gen {};
    std::uniform_int_distribution<ContigRegion::Size> position_dist {0, contig_size};
    std::vector<IntervalTreeRead> result {};
    result.reserve(num_reads);
    std::generate_n(std::back_inserter(result), num_reads, [&] () -> IntervalTreeRead { 
        auto begin = position_dist(gen);
        return {begin, std::min(begin + read_size, contig_size), generate_sequence(read_size)};
    });
    return result;
}

auto generate_interval_tree_data(const std::vector<Read>& reads)
{
    std::vector<IntervalTreeRead> result {};
    result.reserve(reads.size());
    std::transform(std::cbegin(reads), std::cend(reads), std::back_inserter(result), 
                   [] (const auto& read) -> IntervalTreeRead {
                       return {mapped_begin(read), mapped_end(read), read.sequence};
                   });
    return result;
}

int main()
{
    constexpr std::size_t num_reads {10'000'000};
    constexpr unsigned read_size {150};
    constexpr ContigRegion::Size contig_size {50'000'000};
    
    std::cout << "Generating test data..." << std::endl;
    auto mappable_data = generate_mappable_data(num_reads, read_size, contig_size);
    auto interval_tree_data = generate_interval_tree_data(mappable_data);
    
    const ContigRegion small_test_region {contig_size / 2 - 200, contig_size / 2 + 200};
    const ContigRegion big_test_region {contig_size / 2 - contig_size / 4, contig_size / 2 + contig_size / 4};
    
    using D = std::chrono::microseconds;
    auto start = std::chrono::system_clock::now();
    auto end = start;
    auto duration = std::chrono::duration_cast<D>(end - start);
    
    std::cout << "Starting benchmarks..." << std::endl;
    
    //
    // IntervalTree tests
    //
    {
        start = std::chrono::system_clock::now();
        IntervalTree<std::string> tree {interval_tree_data}; // passes by ref
        end = std::chrono::system_clock::now();
        duration = std::chrono::duration_cast<D>(end - start);
        std::cout << "IntervalTree<std::string> constructed in " << duration.count() << "ms" << std::endl;
        
        start = std::chrono::system_clock::now();
        std::vector<IntervalTreeRead> interval_tree_overlapping_small;
        tree.findOverlapping(small_test_region.begin(), small_test_region.end(), interval_tree_overlapping_small);
        end = std::chrono::system_clock::now();
        duration = std::chrono::duration_cast<D>(end - start);
        std::cout << "interval_tree_overlapping_small.size() = " << interval_tree_overlapping_small.size() 
                  << ". Calculated in " << duration.count() << "ms" << std::endl;
        interval_tree_overlapping_small.clear();
        interval_tree_overlapping_small.shrink_to_fit();
        
        start = std::chrono::system_clock::now();
        std::vector<IntervalTreeRead> interval_tree_contained_small;
        tree.findContained(small_test_region.begin(), small_test_region.end(), interval_tree_contained_small);
        end = std::chrono::system_clock::now();
        duration = std::chrono::duration_cast<D>(end - start);
        std::cout << "interval_tree_contained_small.size() = " << interval_tree_contained_small.size() 
                  << ". Calculated in " << duration.count() << "ms" << std::endl;
        interval_tree_contained_small.clear();
        interval_tree_contained_small.shrink_to_fit();
        
        start = std::chrono::system_clock::now();
        std::vector<IntervalTreeRead> interval_tree_overlapping_big;
        tree.findOverlapping(big_test_region.begin(), big_test_region.end(), interval_tree_overlapping_big);
        end = std::chrono::system_clock::now();
        duration = std::chrono::duration_cast<D>(end - start);
        std::cout << "interval_tree_overlapping_big.size() = " << interval_tree_overlapping_big.size() 
                  << ". Calculated in " << duration.count() << "ms" << std::endl;
        interval_tree_overlapping_big.clear();
        interval_tree_overlapping_big.shrink_to_fit();
        
        start = std::chrono::system_clock::now();
        std::vector<IntervalTreeRead> interval_tree_contained_big;
        tree.findContained(big_test_region.begin(), big_test_region.end(), interval_tree_contained_big);
        end = std::chrono::system_clock::now();
        duration = std::chrono::duration_cast<D>(end - start);
        std::cout << "interval_tree_contained_big.size() = " << interval_tree_contained_big.size() 
                  << ". Calculated in " << duration.count() << "ms" << std::endl;
        interval_tree_contained_big.clear();
        interval_tree_contained_big.shrink_to_fit();
    }
    
    //
    // Mappable tests
    //
    {
        start = std::chrono::system_clock::now();
        ReadSet set {std::make_move_iterator(std::begin(mappable_data)), 
                     std::make_move_iterator(std::end(mappable_data))};
        end = std::chrono::system_clock::now();
        duration = std::chrono::duration_cast<D>(end - start);
        std::cout << "ReadSet constructed in " << duration.count() << "ms" << std::endl;
        
        start = std::chrono::system_clock::now();
        const auto mappable_overlapping_small = overlap_range(set, small_test_region);
        end = std::chrono::system_clock::now();
        duration = std::chrono::duration_cast<D>(end - start);
        std::cout << "size(mappable_overlapping_small) = " << size(mappable_overlapping_small)
                  << ". Calculated in " << duration.count() << "ms" << std::endl;
        
        start = std::chrono::system_clock::now();
        const auto mappable_contained_small = contained_range(set, small_test_region);
        end = std::chrono::system_clock::now();
        duration = std::chrono::duration_cast<D>(end - start);
        std::cout << "size(mappable_contained_small) = " << size(mappable_contained_small)
                  << ". Calculated in " << duration.count() << "ms" << std::endl;
        
        start = std::chrono::system_clock::now();
        const auto mappable_overlapping_big = overlap_range(set, big_test_region);
        end = std::chrono::system_clock::now();
        duration = std::chrono::duration_cast<D>(end - start);
        std::cout << "size(mappable_overlapping_big) = " << size(mappable_overlapping_big)
                  << ". Calculated in " << duration.count() << "ms" << std::endl;
        
        start = std::chrono::system_clock::now();
        const auto mappable_contained_big = contained_range(set, big_test_region);
        end = std::chrono::system_clock::now();
        duration = std::chrono::duration_cast<D>(end - start);
        std::cout << "size(mappable_contained_big) = " << size(mappable_contained_big)
                  << ". Calculated in " << duration.count() << "ms" << std::endl;
    }
}
