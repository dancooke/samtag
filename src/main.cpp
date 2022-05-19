// Copyright (c) 2022 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include <fstream>
#include <unordered_map>
#include <vector>
#include <string>
#include <string_view>
#include <algorithm>
#include <cstddef>
#include <iterator>
#include <tuple>
#include <iostream>
#include <memory>
#include <optional>
#include <filesystem>
#include <variant>
#include <type_traits>
#include <array>
#include <regex>
#include <functional>
#include <cassert>

#include <htslib/hts.h>
#include <htslib/sam.h>

#include "version.hpp"

namespace fs = std::filesystem;

const static std::string program {"samtag"};

void print_version()
{
    std::cout << program << " " << VERSION_MAJOR << "." << VERSION_MINOR;
}

void 
print_usage()
{
    std::cout << "Usage: " << program << " <command> [options] " << '\n';
    std::cout << '\n';
    std::cout << "commands:" << '\n';
    std::cout << "  tag    add tags by read name" << '\n';
    std::cout << "  stats  generate stats by tag" << '\n';
    std::cout << std::endl;
}

void 
print_usage(const std::string& subcommand, 
            const std::string& required)
{
    std::cout << "Usage: " << program << " " << subcommand << " [options] " << required << std::endl;
}

struct HtsFileDeleter
{
    void operator()(htsFile* file) const { hts_close(file); }
};
struct HtsHeaderDeleter
{
    void operator()(bam_hdr_t* header) const { bam_hdr_destroy(header); }
};
struct HtsIndexDeleter
{
    void operator()(hts_idx_t* index) const { hts_idx_destroy(index); }
};
struct HtsIteratorDeleter
{
    void operator()(hts_itr_t* iterator) const { sam_itr_destroy(iterator); }
};
struct HtsBam1Deleter
{
    void operator()(bam1_t* b) const { bam_destroy1(b); }
};

//
// tag
//

using ReadNameMap = std::unordered_map<std::string, std::string>;

auto estimate_line_count(const fs::path& filename)
{
    const auto bytes = fs::file_size(filename);
    std::ifstream file {filename};
    std::string line {};
    std::getline(file, line);
    return bytes / (line.length() + 1);
}

auto load_reads(const fs::path& qnames_tsv_path, const bool verbose = false)
{
    const std::size_t log_tick {10'000'000};
    const auto estimated_lines = estimate_line_count(qnames_tsv_path);
    std::ifstream qname_tsv {qnames_tsv_path};
    std::vector<std::pair<std::string, std::string>> reads {};
    reads.reserve(estimated_lines);
    std::string line {};
    for (std::size_t i {0}; std::getline(qname_tsv, line); ++i) {
        if (verbose && i > 0 and i % log_tick == 0) {
            std::clog << "Loaded " << i << " reads" << std::endl;
        }
        if (!line.empty()) {
            if (line.back() == '\r') line.pop_back();
            const auto name_end_pos = line.find('\t');
            reads.emplace_back(line.substr(0, name_end_pos), name_end_pos != std::string::npos ? line.substr(name_end_pos + 1) : "");
        }
    }
    return ReadNameMap {
        std::make_move_iterator(std::begin(reads)),
        std::make_move_iterator(std::end(reads))
    };
}

struct Tag
{
    using Name = std::array<char, 2>;
    using Value = std::variant<std::string_view, long, float>;
    Name id;
    Value value;
};

Tag::Value get_tag_value(const std::string_view value)
{
    try {
        if (value.find('.') == std::string::npos) {
            return std::stol(std::string{value});
        } else {
            return std::stof(std::string{value});
        }
    } catch (const std::invalid_argument&) {
        return value;
    }
}

auto get_tag(const std::string_view& tag)
{
    if (tag.size() < 2 || tag.size() == 3 || (tag.size() > 3 && tag[2] != ':')) {
        std::clog << "Invalid tag " << tag << " (required TAG:VALUE)" << std::endl;
        exit(1);
    }
    Tag result {};
    std::copy_n(std::cbegin(tag), 2, std::begin(result.id));
    if (tag.size() > 3) {
        result.value = get_tag_value(tag.substr(3));
    }
    return result;
}

template<class> inline constexpr bool always_false_v = false;

void add_tag(const Tag& tag, bam1_t* rec)
{
    std::visit([&] (auto&& value) {
        using T = std::remove_cv_t<std::decay_t<decltype(value)>>;
        if constexpr (std::is_same_v<T, long>) {
            bam_aux_update_int(rec, tag.id.data(), value);
        } else if constexpr (std::is_same_v<T, float>) {
            bam_aux_update_float(rec, tag.id.data(), value);
        } else if constexpr (std::is_same_v<T, std::string_view>) {
            if (!value.empty()) {
                bam_aux_update_str(rec, tag.id.data(), value.size() + 1, value.data());
            }
        } else { 
            static_assert(always_false_v<T>, "non-exhaustive visitor!");
        }
    }, tag.value);
}

void 
tag_reads(const fs::path& src_bam_path,
          const ReadNameMap& read_names,
          const std::optional<fs::path>& dst_bam_path,
          const std::optional<Tag>& tag = std::nullopt,
          const std::optional<std::uint16_t> flag = std::nullopt,
          const bool verbose = false)
{
    const std::size_t log_tick {10'000'000};
    std::unique_ptr<htsFile, HtsFileDeleter> src_bam {sam_open(src_bam_path.c_str(), "r"), HtsFileDeleter {}};
    std::unique_ptr<bam_hdr_t, HtsHeaderDeleter> header {sam_hdr_read(src_bam.get()), HtsHeaderDeleter {}};
    std::unique_ptr<htsFile, HtsFileDeleter> dst_bam {
        dst_bam_path ? sam_open(dst_bam_path->c_str(), "[w]b") : sam_open("-", "w"), 
        HtsFileDeleter {}};
    if (sam_hdr_write(dst_bam.get(), header.get()) < 0) {
        std::clog << "Error writing BAM" << std::endl;
        exit(1);
    }
    std::unique_ptr<bam1_t, HtsBam1Deleter> rec {bam_init1(), HtsBam1Deleter {}};
    int r;
    std::string name {};
    std::size_t i {0}, marked {0};
    std::vector<Tag> tags {};
    bool set_default_tag_value {};
    if (tag) {
        tags.push_back(*tag);
        try {
            // If no value is provided for the default tag then
            // the input file contains values
            set_default_tag_value = std::get<std::string_view>(tag->value).empty();
        } catch (const std::bad_variant_access&) {}
    }
    while ((r = sam_read1(src_bam.get(), header.get(), rec.get())) >= 0) {
        if (verbose && i > 0 and i % log_tick == 0) {
            std::clog << "Processed " << i << " reads -- marked " << marked
                      << " (~" << static_cast<int>(100 * marked / i) << "%)"
                      << std::endl;
        }
        name = bam_get_qname(rec.get());
        const auto read_itr = read_names.find(name);
        if (read_itr != std::cend(read_names)) {
            const std::string& edits {read_itr->second};
            auto read_flag = flag;
            if (!edits.empty()) {
                const auto tags_end = edits.find('\t');
                std::string_view reads_tags {edits};
                if (tags_end != std::string::npos) {
                    auto f = static_cast<std::uint16_t>(std::stoi(edits.substr(tags_end + 1)));
                    if (read_flag) {
                        *read_flag |= f;
                    } else {
                        read_flag = f;
                    }
                    reads_tags = edits.substr(0, tags_end);
                }
                // TODO - accept multiple tags seperated with ;
                if (set_default_tag_value) {
                    tags.front().value = get_tag_value(reads_tags);
                } else {
                    tags.push_back(get_tag(reads_tags));
                }
            }
            if (read_flag) {
                rec->core.flag |= *read_flag;
            }
            if (tags.empty() && !read_flag) {
                std::clog << "WARN: no tags or flags for read " << name << std::endl;
            }
            for (const auto& t : tags) {
                add_tag(t, rec.get());
            }
            if (tag) {
                tags.resize(1);
            } else {
                tags.clear();
            }
            ++marked;
        }
        if (sam_write1(dst_bam.get(), header.get(), rec.get()) < 0) {
            std::clog << "Error writing BAM" << std::endl;
            exit(1);
        }
        ++i;
    }
    if (verbose) {
        std::clog << "Processed " << i << " reads -- marked " << marked
                  << " (~" << static_cast<int>(100 * marked / i) << "%)"
                  << std::endl;
    }
}

struct TagOptions
{
    fs::path src_bam_path, qname_tsv_path;
    std::optional<fs::path> output;
    std::optional<Tag> tag;
    std::optional<std::uint16_t> flag;
    bool build_index;
    int verbose = 0;
};

void print_tag_help()
{
    std::cout << "Options:" << std::endl;
    std::cout << "  --help              Print help information" << std::endl;
    std::cout << "  -o, --output FILE   Output bam/cram FILE" << std::endl;
    std::cout << "  -t, --tag STR1:STR2 Add tag STR1 with value STR2 to selected reads" << std::endl;
    std::cout << "  -f, --flag FLAG     Add FLAG to selected reads" << std::endl;
    std::cout << "  -i, --index         Build the bam/cram index for the output" << std::endl;
    std::cout << "  --verbosity INT     Print verbose logging info [0]" << std::endl;
    std::cout << "  --version           Print version information";
}

auto 
pars_tag_args(const int argc, char** argv)
{
    assert(argc > 1);
    const std::string command {"tag"};
    const std::string usage_required {"<in.bam> <qnames.tsv>"};
    constexpr int num_positionals {2};
    TagOptions result {};
    std::vector<std::string_view> args(argv + 2, argv + argc), positionals {};
    std::optional<std::string_view> option {};
    for (int consumed {2}; const auto& arg : args) {
        if (arg == "-h" || arg == "--help") {
            print_usage(command, usage_required);
            std::cout << std::endl;
            print_tag_help();
            std::cout << std::endl;
            exit(0);
        } else if (arg == "--version") {
            print_version();
            std::cout << std::endl;
            exit(0);
        } else if (argc - consumed <= num_positionals - static_cast<int>(positionals.size())) {
            positionals.push_back(arg);
        } else if (arg == "-i" || arg == "--index") {
            result.build_index = true;
        } else if (arg.starts_with("-")) {
            option = arg;
        } else if (option == "-o" || option == "--output") {
            result.output = arg;
            option = std::nullopt;
        } else if (option == "-t" || option == "--tag") {
            result.tag = get_tag(arg);
            option = std::nullopt;
        } else if (option == "-f" || option == "--flag") {
            result.flag = static_cast<std::uint16_t>(std::stoi(std::string(arg)));
            option = std::nullopt;
        } else if (option == "--verbosity") {
            result.verbose = std::stoi(std::string(arg));
            option = std::nullopt;

        } else {
            positionals.push_back(arg);
        }
        ++consumed;
    }
    if (positionals.size() != num_positionals) {
        print_usage(command, usage_required);
        exit(1);
    }
    result.src_bam_path = positionals[0];
    result.qname_tsv_path = positionals[1];
    return result;
}

void check_input(const TagOptions& options)
{
    if (options.src_bam_path != "-" && !fs::exists(options.src_bam_path)) {
        std::cerr << "ERROR: input file " << options.src_bam_path << " does not exist." << std::endl;
        exit(1);
    }
    if (!fs::exists(options.qname_tsv_path)) {
        std::cerr << "ERROR: input file " << options.qname_tsv_path << " does not exist." << std::endl;
        exit(1);
    }
    if (options.build_index && !options.output) {
        std::clog << "Warn: cannot build bam index without --output!" << std::endl;
    }
}

void samtag_tag(const int argc, char** argv)
{
    const auto options = pars_tag_args(argc, argv);
    check_input(options);
    auto read_names = load_reads(options.qname_tsv_path, options.verbose);
    if (options.verbose > 0) std::clog << "Loaded " << read_names.size() << " read names" << std::endl;
    tag_reads(options.src_bam_path, read_names, options.output, options.tag, options.flag, options.verbose);
    if (options.build_index && options.output && sam_index_build(options.output->c_str(), 0) < 0) {
        std::clog << "Failed to build bam index" << std::endl;
    }
}

//
// stats
//

bool operator==(const hts_pair_pos_t& lhs, const hts_pair_pos_t& rhs) noexcept
{
    return lhs.beg == rhs.beg && lhs.end == rhs.end;
}
bool operator<(const hts_pair_pos_t& lhs, const hts_pair_pos_t& rhs) noexcept
{
    return lhs.beg < rhs.beg || (lhs.beg == rhs.beg && lhs.end < rhs.end);
}

auto read_bed_regions_by_contig(const fs::path& bed_path)
{
    std::ifstream bed {bed_path};
    std::unordered_map<std::string, std::vector<hts_pair_pos_t>> result {};
    std::string line {}, contig {}, start_str {}, stop_str {};
    std::size_t i {1};
    hts_pair_pos_t region {};
    while (std::getline(bed, line)) {
        auto contig_itr = std::find(std::cbegin(line), std::cend(line), '\t');
        contig.assign(std::cbegin(line), contig_itr);
        if (contig_itr == std::cend(line)) {
            std::cerr << "ERROR: malformed bed file " << bed_path << " at line " << i << std::endl;
            exit(1);
        }
        ++contig_itr;
        auto start_itr = std::find(contig_itr, std::cend(line), '\t');
        start_str.assign(contig_itr, start_itr);
        region.beg = std::stoul(start_str);
        if (start_itr == std::cend(line)) {
            std::cerr << "ERROR: malformed bed file " << bed_path << " at line " << i << std::endl;
            exit(1);
        }
        ++start_itr;
        auto stop_itr = std::find(start_itr, std::cend(line), '\t');
        stop_str.assign(start_itr, stop_itr);
        region.end = std::stoul(stop_str);
        result[contig].push_back(region);
        ++i;
    }
    return result;
}

auto merge(const std::vector<hts_pair_pos_t>& intervals)
{
    std::vector<hts_pair_pos_t> result {};
    if (intervals.empty()) return result;
    result.reserve(intervals.size());
    auto current = std::cbegin(intervals);
    auto overlapped = current;
    auto rightmost = current;
    for (const auto last = std::cend(intervals); current != last; ++current) {
        if (current->beg > rightmost->end) {
            if (result.empty() || result.back().end != rightmost->end) {
                result.push_back({overlapped->beg, rightmost->end});
            }
            rightmost = current;
            overlapped = current;
        } else if (current->end >= rightmost->end) {
            rightmost = current;
        }
    }
    result.push_back({overlapped->beg, rightmost->end});
    return result;
}

auto get_nonoverlapping_regions_by_contig(const fs::path& bed_path)
{
    auto regions = read_bed_regions_by_contig(bed_path);
    for (auto&& [contig, intervals] : regions) {
        std::sort(std::begin(intervals), std::end(intervals));
        intervals = merge(intervals);
    }
    return regions;
}

struct IntervalStats
{
    std::size_t num_contigs, num_targets, num_bases;
};

auto length(const hts_pair_pos_t& interval) noexcept
{
    return interval.end - interval.beg;
}

auto calculate_interval_stats(const std::unordered_map<std::string, std::vector<hts_pair_pos_t>>& regions)
{
    IntervalStats result {};
    for (const auto& [contig, intervals] : regions) {
        ++result.num_contigs;
        for (const auto& interval : intervals) {
            ++result.num_targets;
            result.num_bases += length(interval);
        }
    }
    return result;
}

auto get_tid(const std::string& contig, const bam_hdr_t& header)
{
    auto itr = std::find(header.target_name, header.target_name + header.n_targets, contig.c_str());
    return static_cast<int>(std::distance(header.target_name, itr));
}

auto make_hts_region_list(std::unordered_map<std::string, std::vector<hts_pair_pos_t>>& regions, const bam_hdr_t& header)
{
    std::vector<hts_reglist_t> result {};
    result.reserve(regions.size());
    for (auto& [contig, intervals] : regions) {
        assert(!intervals.empty());
        hts_reglist_t r {};
        r.reg = contig.c_str();
        r.tid = get_tid(contig, header);
        r.count = intervals.size();
        r.intervals = intervals.data();
        r.min_beg = intervals.front().beg;
        r.max_end = intervals.back().end;
        result.push_back(r);
    }
    return result;
}

struct SearchTag
{
    Tag::Name id;
    std::optional<std::string> value;
    std::optional<std::regex> pattern;
};

bool operator==(const SearchTag& lhs, const SearchTag& rhs) noexcept
{
    return lhs.id == rhs.id && lhs.value == rhs.value;
}

template <class T>
inline void hash_combine(std::size_t& seed, const T& v)
{
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
}

template <typename T, size_t N> 
struct std::hash<std::array<T, N>>
{
    std::size_t operator()(const std::array<T,N>& arr) const noexcept
    {
        std::size_t result {};
        for (const auto& v : arr) {
            hash_combine(result, v);
        }
        return result;
    }
};

struct SearchTagHash
{
    auto operator()(const SearchTag& tag) const noexcept
    {
        auto result = std::hash<decltype(tag.id)>{}(tag.id);
        hash_combine(result, tag.value);
        return result;
    }
};

using TagCountMap = std::unordered_map<SearchTag, std::size_t, SearchTagHash>;

struct TagStats
{
    TagCountMap counts;
    std::optional<TagCountMap> value_counts;
    std::size_t total_reads;
};

void add_read(const bam1_t& read, TagStats& stats)
{
    for (auto& [tag, count] : stats.counts) {
        const auto p = bam_aux_get(&read, tag.id.data());
        if (p) {
            if (tag.pattern) {
                const auto value = bam_aux2Z(p);
                // TODO - what if type is not string?
                if (value && std::regex_search(value, *tag.pattern)) {
                    ++count;
                    if (stats.value_counts) {
                        auto valued_tag = tag;
                        valued_tag.value = value;
                        ++((*stats.value_counts)[std::move(valued_tag)]);
                    }
                }
            } else {
                ++count;
                if (stats.value_counts) {
                    std::string value {};
                    switch (*p) {
                    case 'Z': value = bam_aux2Z(p); break;
                    case 'i': value = std::to_string(bam_aux2i(p)); break;
                    case 'f': value = std::to_string(bam_aux2f(p)); break;
                    }
                    auto valued_tag = tag;
                    valued_tag.value = std::move(value);
                    ++((*stats.value_counts)[std::move(valued_tag)]);
                }
            }
        }
    }
    ++stats.total_reads;
}

template <typename Range>
void write(const Range& values, std::ostream& os, const char delimiter='\t')
{
    using T = typename std::iterator_traits<typename Range::const_iterator>::value_type;
    std::copy(std::cbegin(values), std::prev(std::cend(values)), std::ostream_iterator<T> {os});
    os << values.back();
}

void write(const TagStats& stats, std::ostream& os, const bool sorted = false, const char delimiter='\t')
{
    os << "tag" << delimiter << "value" << delimiter << "count" << delimiter << "fraction" << '\n';
    os << '*' << delimiter << '*' << delimiter << stats.total_reads << delimiter << '1' << '\n';
    const auto write_row = [&] (const auto& tag, const auto count) {
        std::copy(std::cbegin(tag.id), std::cend(tag.id), std::ostreambuf_iterator<char> {os});
        os << delimiter;
        if (tag.value) {
            os << *tag.value;
        } else {
            os << '*';
        }
        os << delimiter;
        os << count << delimiter;
        if (stats.total_reads > 0) {
            os << static_cast<float>(count) / stats.total_reads;
        } else {
            os << '0';
        }
        os << std::endl;
    };
    if (sorted) {
        std::vector<std::pair<SearchTag, std::size_t>> counts {
            std::begin(stats.counts), std::end(stats.counts)
        };
        if (stats.value_counts) {
            counts.insert(std::end(counts), std::begin(*stats.value_counts), std::end(*stats.value_counts));
        }
        const static auto count_greater = [] (const auto& lhs, const auto& rhs) noexcept {
            return lhs.second > rhs.second;
        };
        std::sort(std::begin(counts), std::end(counts), count_greater);
        for (const auto& [tag, count] : counts) {
            write_row(tag, count);
        }
    } else {
        for (const auto& [tag, count] : stats.counts) {
            write_row(tag, count);
        }
        if (stats.value_counts) {
            for (const auto& [tag, count] : *stats.value_counts) {
                write_row(tag, count);
            }
        }
    }
}

struct StatsOptions
{
    fs::path bam_path;
    std::optional<fs::path> tag_path, bed_path, output;
    std::vector<SearchTag> tags;
    std::optional<std::uint16_t> require_flags, exclude_flag;
    std::optional<int> min_mapping_quality;
    bool split = false, sort = false;
    int verbose = 0;
};

void print_stats_help()
{
    std::cout << "Main options:" << std::endl;
    std::cout << "  -t, --tag STR1[:STR2]    Compute stats for tag STR1 with pattern STR2" << std::endl;
    std::cout << "  --tag-file FILE          Compute stats for tags listed in FILE" << std::endl;
    std::cout << "  --split                  Split tag values for matching patterns" << std::endl;
    std::cout << std::endl;
    std::cout << "Filtering options (Only include reads that...):" << std::endl;
    std::cout << "  -L, --target-regions     ...overlap (BED) regions in FILE" << std::endl;
    std::cout << "  -f, --require-flag FLAG  ...have all FLAGs present" << std::endl;
    std::cout << "  -F, --exclude-flag FLAG  ...have none of FLAGs present" << std::endl;
    std::cout << "  -q, --min-mapq INT       ...have mapping quality >= INT" << std::endl;
    std::cout << std::endl;
    std::cout << "Output options:" << std::endl;
    std::cout << "  -o, --output FILE        Output to FILE" << std::endl;
    std::cout << "  --sort                   Sort output by count" << std::endl;
    std::cout << std::endl;
    std::cout << "General options:" << std::endl;
    std::cout << "  --help                   Print help information" << std::endl;
    std::cout << "  --verbosity INT          Print verbose logging info [0]" << std::endl;
    std::cout << "  --version                Print version information";
}

auto parse_search_tag(const std::string_view tag)
{
    if (tag.size() < 2 || tag.size() == 3 || (tag.size() > 3 && tag[2] != ':')) {
        std::clog << "Invalid tag " << tag << " (required TAG[:VALUE])" << std::endl;
        exit(1);
    }
    SearchTag result {};
    std::copy_n(std::cbegin(tag), 2, std::begin(result.id));
    if (tag.size() > 3) {
        result.value = tag.substr(3);
        result.pattern = std::regex(*result.value);
    }
    return result;
}

auto parse_stats_args(const int argc, char** argv)
{
    assert(argc > 1);
    const std::string subcommand {"tag"};
    const std::string usage_required {"<in.bam>"};
    constexpr int num_positionals {1};
    StatsOptions result {};
    std::vector<std::string_view> args(argv + 2, argv + argc), positionals {};
    std::optional<std::string_view> option {};
    for (int consumed {2}; const auto& arg : args) {
        if (arg == "-h" || arg == "--help") {
            print_usage(subcommand, usage_required);
            std::cout << std::endl;
            print_stats_help();
            std::cout << std::endl;
            exit(0);
        } else if (arg == "--version") {
            print_version();
            std::cout << std::endl;
            exit(0);
        } else if (argc - consumed <= num_positionals - static_cast<int>(positionals.size())) {
            positionals.push_back(arg);
        } else if (arg == "--split") {
            result.split = true;
            option = std::nullopt;
        } else if (arg == "--sort") {
            result.sort = true;
            option = std::nullopt;
        } else if (arg.starts_with("-")) {
            option = arg;
        } else if (option == "-o" || option == "--output") {
            result.output = arg;
            option = std::nullopt;
        } else if (option == "-L" || option == "--target-regions") {
            result.bed_path = arg;
            option = std::nullopt;
        } else if (option == "-t" || option == "--tag") {
            result.tags.push_back(parse_search_tag(arg));
        } else if (option == "--tag-file") {
            result.tag_path = arg;
            option = std::nullopt;
        } else if (option == "-f" || option == "--require-flag") {
            result.require_flags = static_cast<std::uint16_t>(std::stoi(std::string(arg)));
            option = std::nullopt;
        } else if (option == "-F" || option == "--exclude-flag") {
            result.exclude_flag = static_cast<std::uint16_t>(std::stoi(std::string(arg)));
            option = std::nullopt;
        } else if (option == "-q" || option == "--min-mapq") {
            result.min_mapping_quality = std::stoi(std::string(arg));
            option = std::nullopt;
        } else if (option == "--verbosity") {
            result.verbose = std::stoi(std::string(arg));
            option = std::nullopt;
        } else {
            positionals.push_back(arg);
        }
        ++consumed;
    }
    if (positionals.size() != num_positionals) {
        print_usage(subcommand, usage_required);
        exit(1);
    }
    result.bam_path = positionals[0];
    return result;
}

void check_input(const StatsOptions& options)
{
    if (options.bam_path != "-" && !fs::exists(options.bam_path)) {
        std::cerr << "ERROR: input file " << options.bam_path << " does not exist." << std::endl;
        exit(1);
    }
    if (options.bed_path && !fs::exists(*options.bed_path)) {
        std::cerr << "ERROR: input file " << *options.bed_path << " does not exist." << std::endl;
        exit(1);
    }
    if (options.tags.empty() && !options.tag_path) {
        std::cerr << "ERROR: one of --tag or --tag-file is required." << std::endl;
        exit(1);
    }
}

auto init_stats(const StatsOptions& options)
{
    TagStats result {};
    result.counts.reserve(options.tags.size());
    for (const auto& tag : options.tags) {
        result.counts.emplace(tag, 0);
    }
    if (options.split) {
        result.value_counts.emplace(options.tags.size());
    }
    return result;
}

std::function<bool(const bam1_t&)> make_read_filter(const StatsOptions& options)
{
    std::vector<std::function<bool(const bam1_t&)>> conditions {};
    conditions.reserve(3);
    if (options.require_flags) {
        conditions.push_back([flag = *options.require_flags] (const bam1_t& read) { 
            return (read.core.flag & flag) == flag; });
    }
    if (options.exclude_flag) {
        conditions.push_back([flag = *options.exclude_flag] (const bam1_t& read) { 
            return (read.core.flag & flag) == 0; });
    }
    if (options.min_mapping_quality) {
        conditions.push_back([min_mq = options.min_mapping_quality] (const bam1_t& read) { 
            return read.core.qual >= min_mq; });
    }
    return [conditions = std::move(conditions)] (const bam1_t& read) -> bool { 
        return std::all_of(std::cbegin(conditions), std::cend(conditions), 
            [&] (const auto& condition) { return condition(read); }); 
        };
}

void log_progress(const TagStats& stats, const std::size_t total_reads, const int verbosity, const bool force = false)
{
    constexpr static std::size_t log_tick {10'000'000};
    if (verbosity > 0 && ((total_reads > 0 && (total_reads % log_tick == 0)) || force)) {
         std::clog << "Processed " << total_reads << " reads (used " << stats.total_reads << ")" << std::endl;
    }
}

void 
samtag_stats(const int argc, char** argv)
{
    const auto options = parse_stats_args(argc, argv);
    check_input(options);
    const auto read_filter = make_read_filter(options);
    std::unique_ptr<htsFile, HtsFileDeleter> bam {sam_open(options.bam_path.c_str(), "r"), HtsFileDeleter {}};
    std::unique_ptr<bam_hdr_t, HtsHeaderDeleter> header {sam_hdr_read(bam.get()), HtsHeaderDeleter {}};
    std::unique_ptr<bam1_t, HtsBam1Deleter> read {bam_init1(), HtsBam1Deleter {}};
    auto stats = init_stats(options);
    if (options.verbose > 0) {
        std::clog << "Loaded " << stats.counts.size() << " tags" << std::endl;
    }
    int r;
    std::size_t tot_reads {0};
    if (options.bed_path) {
        auto intervals = get_nonoverlapping_regions_by_contig(*options.bed_path);
        auto hts_intervals = make_hts_region_list(intervals, *header);
        if (options.verbose > 0) {
            auto interval_stats = calculate_interval_stats(intervals);
            std::clog << "Loaded " << interval_stats.num_targets << " non-overlapping targets ("
                      << interval_stats.num_bases << " bp) in " << interval_stats.num_contigs << " contigs " << std::endl;
        }
        std::unique_ptr<hts_idx_t, HtsIndexDeleter> index {sam_index_load(bam.get(), options.bam_path.c_str()), HtsIndexDeleter {}};
        auto itr = sam_itr_regions(index.get(), header.get(), hts_intervals.data(), hts_intervals.size());
        while ((r = sam_itr_next(bam.get(), itr, read.get())) >= 0) {
            if (read_filter(*read)) {
                add_read(*read, stats);
            }
            log_progress(stats, ++tot_reads, options.verbose);
        }
    } else {
        while ((r = sam_read1(bam.get(), header.get(), read.get())) >= 0) {
            if (read_filter(*read)) {
                add_read(*read, stats);
            }
            log_progress(stats, ++tot_reads, options.verbose);
        }
    }
    if (r < -1) {
        std::cerr << "Error reading " << options.bam_path << std::endl;
        exit(1);
    }
    log_progress(stats, tot_reads++, options.verbose, true);
    if (options.output) {
        std::ofstream output {*options.output};
        write(stats, output, options.sort);
    } else {
        write(stats, std::cout, options.sort);
    }
}

int main(int argc, char** argv)
{
    if (argc < 2) {
        print_usage();
        exit(1);
    }
    const std::string command {argv[1]};
    if (command == "tag") {
        samtag_tag(argc, argv);
    } else if (command == "stats") {
        samtag_stats(argc, argv);
    } else {
        std::cerr << "Unknown command " << command << std::endl;
        std::cerr << std::endl;
        print_usage();
        exit(1);
    }
}
